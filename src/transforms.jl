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

# ──────────────────────────────────────────────────────────────────────────────
# Single-coordinate transforms
# ──────────────────────────────────────────────────────────────────────────────

@inline _coordinate_float_type(::AbstractVector{T}) where {T<:Real} = _float_type(T)

@inline _coordinate_float_type(coords::Tuple{Vararg{Real}}) = _promote_float_type(coords...)

function _coordinate_vector(coords::Tuple{Vararg{Real}})
    # Preserve the tuple's promoted floating-point precision when materializing it.
    T = _coordinate_float_type(coords)
    return collect(T, coords)
end

function _world_from_intermediate(wcs::WCSTransform, intermediate::AbstractVector, ::Type{T}) where {T<:AbstractFloat}
    world = Vector{T}(undef, wcs.naxis)

    # Add CRVAL axis-by-axis so WCS metadata is converted to the coordinate type.
    for i in 1:wcs.naxis
        world[i] = T(wcs.crval[i]) + T(intermediate[i])
    end

    return world
end

function _world_offsets(wcs::WCSTransform, world::AbstractVector, ::Type{T}) where {T<:AbstractFloat}
    intermediate = Vector{T}(undef, wcs.naxis)

    # Subtract CRVAL axis-by-axis to avoid broadcast promotion through WCS storage.
    for i in 1:wcs.naxis
        intermediate[i] = T(world[i]) - T(wcs.crval[i])
    end

    return intermediate
end

"""
    pixel_to_world(wcs, pixel) -> world

Convert a single FITS pixel coordinate to world coordinates.

`pixel` must be a length-`wcs.naxis` vector (or tuple) of pixel positions.
Returns a length-`wcs.naxis` floating-point vector of world coordinates.

For purely linear WCS (no celestial projection) the result is:
    world[i] = CRVAL[i] + Σⱼ CD[i,j] * (pixel[j] - CRPIX[j])

For celestial axes the full pipeline (linear → deprojection → spherical
rotation) is applied.
"""
function pixel_to_world(wcs::WCSTransform, pixel::AbstractVector)
    length(pixel) == wcs.naxis ||
        throw(DimensionMismatch("pixel has length $(length(pixel)), expected $(wcs.naxis)"))
    T = _coordinate_float_type(pixel)

    # Step 1: Linear transform → intermediate world coordinates (degrees)
    x = pixel_to_intermediate(wcs, pixel)  # length naxis

    # Step 2: Non-linear (celestial) part
    if wcs.projection !== nothing
        lon_idx = wcs.lon_axis
        lat_idx = wcs.lat_axis
        x_lon = T(x[lon_idx])   # degrees
        x_lat = T(x[lat_idx])   # degrees

        # Deproject: (x, y) → (φ, θ) in radians
        phi, theta = intermediate_to_native(wcs.projection, x_lon, x_lat)

        # Spherical rotation: (φ, θ) → (α, δ) in radians
        # alpha_p, delta_p are the celestial coords of the native north pole.
        d2r = deg2rad(one(T))
        r2d = inv(d2r)
        alpha_p = T(wcs.alpha_p) * d2r
        delta_p = T(wcs.delta_p) * d2r
        phi_p = T(wcs.lonpole) * d2r

        alpha, delta = native_to_celestial(phi, theta, alpha_p, delta_p, phi_p)

        # Build world vector: start from the intermediate coords, then overwrite
        # the celestial axes with the spherical-rotation results.
        world = _world_from_intermediate(wcs, x, T)
        world[lon_idx] = mod(alpha * r2d, T(360))
        world[lat_idx] = delta * r2d
        return world
    else
        # Purely linear: world = CRVAL + CD*(pixel - CRPIX)
        return _world_from_intermediate(wcs, x, T)
    end
end

# Convenience: accept tuples and static-array-likes
pixel_to_world(wcs::WCSTransform, pixel::Tuple{Vararg{Real}}) =
    pixel_to_world(wcs, _coordinate_vector(pixel))

"""
    world_to_pixel(wcs, world) -> pixel

Convert world coordinates to FITS pixel coordinates.

`world` must be a length-`wcs.naxis` vector (or tuple) of world positions
in the same units as `CRVAL` (degrees for celestial axes).
Returns a length-`wcs.naxis` floating-point vector of 1-based pixel coordinates.
"""
function world_to_pixel(wcs::WCSTransform, world::AbstractVector)
    length(world) == wcs.naxis ||
        throw(DimensionMismatch("world has length $(length(world)), expected $(wcs.naxis)"))
    T = _coordinate_float_type(world)

    if wcs.projection !== nothing
        lon_idx = wcs.lon_axis
        lat_idx = wcs.lat_axis

        # The fiducial celestial point is defined to have zero intermediate
        # coordinates; handling it directly avoids pole cancellation.
        lon_delta = mod(T(world[lon_idx]) - T(wcs.crval[lon_idx]) + T(180), T(360)) - T(180)
        if abs(lon_delta) <= T(1e-10) && abs(T(world[lat_idx]) - T(wcs.crval[lat_idx])) <= T(1e-10)
            x = _world_offsets(wcs, world, T)
            x[lon_idx] = zero(T)
            x[lat_idx] = zero(T)
            return intermediate_to_pixel(wcs, x)
        end

        d2r = deg2rad(one(T))
        alpha_p = T(wcs.alpha_p) * d2r
        delta_p = T(wcs.delta_p) * d2r
        phi_p = T(wcs.lonpole) * d2r

        alpha = T(world[lon_idx]) * d2r
        delta = T(world[lat_idx]) * d2r

        # Spherical rotation: (α, δ) → (φ, θ)
        phi, theta = celestial_to_native(alpha, delta, alpha_p, delta_p, phi_p)

        # Re-project: (φ, θ) → (x, y) in degrees
        x_lon, x_lat = native_to_intermediate(wcs.projection, phi, theta)

        # Build intermediate coordinate vector
        # Non-celestial axes: x_i = world_i - crval_i (trivial linear axes)
        x = _world_offsets(wcs, world, T)
        x[lon_idx] = T(x_lon)
        x[lat_idx] = T(x_lat)

        return intermediate_to_pixel(wcs, x)
    else
        x = _world_offsets(wcs, world, T)
        return intermediate_to_pixel(wcs, x)
    end
end

world_to_pixel(wcs::WCSTransform, world::Tuple{Vararg{Real}}) =
    world_to_pixel(wcs, _coordinate_vector(world))

# ──────────────────────────────────────────────────────────────────────────────
# Batch transforms
# ──────────────────────────────────────────────────────────────────────────────

"""
    pixel_to_world(wcs, pixels::AbstractMatrix) -> Matrix{Float64}

Batch pixel-to-world transform.

`pixels` must be an `naxis × N` matrix where each column is one pixel
coordinate.  Returns an `naxis × N` floating-point matrix of world coordinates.
"""
function pixel_to_world(wcs::WCSTransform, pixels::AbstractMatrix)
    naxis, N = size(pixels)
    naxis == wcs.naxis ||
        throw(DimensionMismatch("pixels has $(naxis) rows, expected $(wcs.naxis)"))
    result = similar(pixels, _float_type(eltype(pixels)), naxis, N)
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
coordinate.  Returns an `naxis × N` floating-point matrix of pixel coordinates.
"""
function world_to_pixel(wcs::WCSTransform, worlds::AbstractMatrix)
    naxis, N = size(worlds)
    naxis == wcs.naxis ||
        throw(DimensionMismatch("worlds has $(naxis) rows, expected $(wcs.naxis)"))
    result = similar(worlds, _float_type(eltype(worlds)), naxis, N)
    for k in 1:N
        result[:, k] = world_to_pixel(wcs, view(worlds, :, k))
    end
    return result
end

function world_to_pixel(wcs::WCSTransform, worlds::AbstractVector{<:AbstractVector})
    # Treat each nested vector as one coordinate in a simple batch.
    return [world_to_pixel(wcs, world) for world in worlds]
end
