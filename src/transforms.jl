"""
High-level pixel ↔ world transforms.

Public API:
- `pixel_to_world(wcs, pixel)` — single coordinate
- `pixel_to_world(wcs, pixels)` — batch (matrix or vector-of-vectors)
- `world_to_pixel(wcs, world)` — single coordinate
- `world_to_pixel(wcs, worlds)` — batch

## Pixel coordinate convention

All pixel coordinates use the **FITS 1-based** convention:
- Pixel `[1, 1, ...]` is the centre of the first array element.
- `CRPIX` values from FITS headers are used directly.

This matches Julia's 1-based array indexing, so no offset is needed when
working with Julia arrays.

## Celestial coordinate convention

For axes with CTYPE codes `RA---*` / `DEC--*` (or galactic, ecliptic
equivalents), the world coordinates are in **degrees**.  The output vector
element at `lon_axis` is the right ascension (or longitude) and at `lat_axis`
is the declination (or latitude).

For purely linear axes the world coordinates are in whatever units the FITS
header specifies via `CDELT` and `CRVAL`.
"""

const _D2R_t = π / 180.0
const _R2D_t = 180.0 / π

# ──────────────────────────────────────────────────────────────────────────────
# Single-coordinate transforms
# ──────────────────────────────────────────────────────────────────────────────

"""
    pixel_to_world(wcs, pixel) -> world

Convert a single FITS pixel coordinate to world coordinates.

`pixel` must be a length-`wcs.naxis` vector (or tuple) of pixel positions.
Returns a length-`wcs.naxis` `Vector{Float64}` of world coordinates.

For purely linear WCS (no celestial projection) the result is:
    world[i] = CRVAL[i] + Σⱼ CD[i,j] * (pixel[j] - CRPIX[j])

For celestial axes the full pipeline (linear → deprojection → spherical
rotation) is applied.
"""
function pixel_to_world(wcs::WCSTransform, pixel::AbstractVector)
    length(pixel) == wcs.naxis ||
        throw(DimensionMismatch("pixel has length $(length(pixel)), expected $(wcs.naxis)"))

    # Step 1: Linear transform → intermediate world coordinates (degrees)
    x = pixel_to_intermediate(wcs, pixel)  # length naxis

    # Step 2: Non-linear (celestial) part
    if wcs.projection !== nothing
        lon_idx = wcs.lon_axis
        lat_idx = wcs.lat_axis

        x_lon = x[lon_idx]   # degrees
        x_lat = x[lat_idx]   # degrees

        # Deproject: (x, y) → (φ, θ) in radians
        phi, theta = intermediate_to_native(wcs.projection, x_lon, x_lat)

        # Spherical rotation: (φ, θ) → (α, δ) in radians
        # alpha_p, delta_p are the celestial coords of the native north pole.
        alpha_p = wcs.alpha_p * _D2R_t
        delta_p = wcs.delta_p * _D2R_t
        phi_p   = wcs.lonpole * _D2R_t

        alpha, delta = native_to_celestial(phi, theta, alpha_p, delta_p, phi_p)

        # Build world vector: start from the intermediate coords, then overwrite
        # the celestial axes with the spherical-rotation results.
        world = wcs.crval .+ x
        world[lon_idx] = mod(alpha * _R2D_t, 360.0)
        world[lat_idx] = delta * _R2D_t
        return world
    else
        # Purely linear: world = CRVAL + CD*(pixel - CRPIX)
        return wcs.crval .+ x
    end
end

# Convenience: accept tuples and static-array-likes
pixel_to_world(wcs::WCSTransform, pixel::Tuple) =
    pixel_to_world(wcs, collect(Float64, pixel))

"""
    world_to_pixel(wcs, world) -> pixel

Convert world coordinates to FITS pixel coordinates.

`world` must be a length-`wcs.naxis` vector (or tuple) of world positions
in the same units as `CRVAL` (degrees for celestial axes).
Returns a length-`wcs.naxis` `Vector{Float64}` of 1-based pixel coordinates.
"""
function world_to_pixel(wcs::WCSTransform, world::AbstractVector)
    length(world) == wcs.naxis ||
        throw(DimensionMismatch("world has length $(length(world)), expected $(wcs.naxis)"))

    if wcs.projection !== nothing
        lon_idx = wcs.lon_axis
        lat_idx = wcs.lat_axis

        # The fiducial celestial point is defined to have zero intermediate
        # coordinates; handling it directly avoids pole cancellation.
        lon_delta = mod(world[lon_idx] - wcs.crval[lon_idx] + 180.0, 360.0) - 180.0
        if abs(lon_delta) <= 1e-10 && abs(world[lat_idx] - wcs.crval[lat_idx]) <= 1e-10
            x = world .- wcs.crval
            x[lon_idx] = 0.0
            x[lat_idx] = 0.0
            return intermediate_to_pixel(wcs, x)
        end

        alpha_p = wcs.alpha_p * _D2R_t
        delta_p = wcs.delta_p * _D2R_t
        phi_p   = wcs.lonpole * _D2R_t

        alpha = world[lon_idx] * _D2R_t
        delta = world[lat_idx] * _D2R_t

        # Spherical rotation: (α, δ) → (φ, θ)
        phi, theta = celestial_to_native(alpha, delta, alpha_p, delta_p, phi_p)

        # Re-project: (φ, θ) → (x, y) in degrees
        x_lon, x_lat = native_to_intermediate(wcs.projection, phi, theta)

        # Build intermediate coordinate vector
        # Non-celestial axes: x_i = world_i - crval_i (trivial linear axes)
        x = world .- wcs.crval
        x[lon_idx] = x_lon
        x[lat_idx] = x_lat

        return intermediate_to_pixel(wcs, x)
    else
        x = world .- wcs.crval
        return intermediate_to_pixel(wcs, x)
    end
end

world_to_pixel(wcs::WCSTransform, world::Tuple) =
    world_to_pixel(wcs, collect(Float64, world))

# ──────────────────────────────────────────────────────────────────────────────
# Batch transforms
# ──────────────────────────────────────────────────────────────────────────────

"""
    pixel_to_world(wcs, pixels::AbstractMatrix) -> Matrix{Float64}

Batch pixel-to-world transform.

`pixels` must be an `naxis × N` matrix where each column is one pixel
coordinate.  Returns an `naxis × N` matrix of world coordinates.
"""
function pixel_to_world(wcs::WCSTransform, pixels::AbstractMatrix)
    naxis, N = size(pixels)
    naxis == wcs.naxis ||
        throw(DimensionMismatch("pixels has $(naxis) rows, expected $(wcs.naxis)"))
    result = similar(pixels, Float64, naxis, N)
    for k in 1:N
        result[:, k] = pixel_to_world(wcs, view(pixels, :, k))
    end
    return result
end

function pixel_to_world(wcs::WCSTransform, pixels::AbstractVector{<:AbstractVector})
    # Treat each nested vector as one coordinate in a simple batch.
    return [pixel_to_world(wcs, pixel) for pixel in pixels]
end

"""
    world_to_pixel(wcs, worlds::AbstractMatrix) -> Matrix{Float64}

Batch world-to-pixel transform.

`worlds` must be an `naxis × N` matrix where each column is one world
coordinate.  Returns an `naxis × N` matrix of pixel coordinates.
"""
function world_to_pixel(wcs::WCSTransform, worlds::AbstractMatrix)
    naxis, N = size(worlds)
    naxis == wcs.naxis ||
        throw(DimensionMismatch("worlds has $(naxis) rows, expected $(wcs.naxis)"))
    result = similar(worlds, Float64, naxis, N)
    for k in 1:N
        result[:, k] = world_to_pixel(wcs, view(worlds, :, k))
    end
    return result
end

function world_to_pixel(wcs::WCSTransform, worlds::AbstractVector{<:AbstractVector})
    # Treat each nested vector as one coordinate in a simple batch.
    return [world_to_pixel(wcs, world) for world in worlds]
end
