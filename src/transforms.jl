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

function _coordinate_vector(coords::Tuple{Vararg{Real, N}}) where {N}
    # Preserve tuple precision and length when materializing scalar arguments.
    T = _coordinate_float_type(coords)
    return SVector{N,T}(coords)
end

@inline function _world_from_intermediate(wcs::WCSTransform{N}, intermediate::AbstractVector, ::Type{T}) where {N, T <: AbstractFloat}
    return SVector{N,T}(ntuple(i -> T(wcs.crval[i]) + T(intermediate[i]), N))
end

@inline function _world_offsets(wcs::WCSTransform{N}, world::AbstractVector, ::Type{T}) where {N,T<:AbstractFloat}
    return SVector{N,T}(ntuple(i -> T(world[i]) - T(wcs.crval[i]), N))
end

function _set_celestial_axes!(coords::AbstractVector, lon_idx::Int, lat_idx::Int, lon, lat)
    # Mutable coordinate storage can be updated in place.
    coords[lon_idx] = lon
    coords[lat_idx] = lat
    return coords
end

function _set_celestial_axes!(coords::StaticVector{N,T}, lon_idx::Int, lat_idx::Int, lon, lat) where {N,T}
    # Immutable static coordinates are rebuilt with the celestial axes replaced.
    return SVector{N,T}(ntuple(i -> i == lon_idx ? T(lon) : i == lat_idx ? T(lat) : coords[i], N))
end

"""
    pixel_to_world(wcs, pixel) -> world

Convert a single FITS pixel coordinate to world coordinates.

`pixel` must be a length-`wcs.naxis` vector (or tuple) of pixel positions.
Returns a length-`wcs.naxis` floating-point vector of world coordinates.
Tuple inputs and scalar varargs are materialized as static coordinates and return
`SVector` results.

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
        x_lon = x[lon_idx]   # degrees
        x_lat = x[lat_idx]   # degrees

        # Deproject: (x, y) → (φ, θ) in radians
        phi, theta = intermediate_to_native(wcs.projection, x_lon, x_lat)

        # Spherical rotation: (φ, θ) → (α, δ) in radians
        # alpha_p, delta_p are the celestial coords of the native north pole.
        alpha_p = deg2rad(T(wcs.alpha_p))
        delta_p = deg2rad(T(wcs.delta_p))
        phi_p = deg2rad(T(wcs.lonpole))

        alpha, delta = native_to_celestial(phi, theta, alpha_p, delta_p, phi_p)
        # Build world vector: start from the intermediate coords, then overwrite
        # the celestial axes with the spherical-rotation results.
        world = _world_from_intermediate(wcs, x, T)
        return _set_celestial_axes!(world, lon_idx, lat_idx, mod(rad2deg(alpha), 360), rad2deg(delta))
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
Tuple inputs and scalar varargs are materialized as static coordinates and return
`SVector` results.
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
        lon_delta = mod(world[lon_idx] - T(wcs.crval[lon_idx]) + 180, 360) - 180
        if abs(lon_delta) <= T(1e-10) && abs(world[lat_idx] - T(wcs.crval[lat_idx])) <= T(1e-10)
            x = _world_offsets(wcs, world, T)
            x = _set_celestial_axes!(x, lon_idx, lat_idx, zero(T), zero(T))
            return intermediate_to_pixel(wcs, x)
        end

        alpha_p = deg2rad(T(wcs.alpha_p))
        delta_p = deg2rad(T(wcs.delta_p))
        phi_p = deg2rad(T(wcs.lonpole))

        alpha = deg2rad(T(world[lon_idx]))
        delta = deg2rad(T(world[lat_idx]))

        # Spherical rotation: (α, δ) → (φ, θ)
        phi, theta = celestial_to_native(alpha, delta, alpha_p, delta_p, phi_p)

        # Re-project: (φ, θ) → (x, y) in degrees
        x_lon, x_lat = native_to_intermediate(wcs.projection, phi, theta)

        # Build intermediate coordinate vector
        # Non-celestial axes: x_i = world_i - crval_i (trivial linear axes)
        x = _world_offsets(wcs, world, T)
        x = _set_celestial_axes!(x, lon_idx, lat_idx, T(x_lon), T(x_lat))

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

@inline _typed_matrix(::Type{T}, values::AbstractMatrix{T}) where {T<:AbstractFloat} = values
@inline _typed_matrix(::Type{T}, values::AbstractMatrix) where {T<:AbstractFloat} = Matrix{T}(values)

function _pixel_to_intermediate_batch(wcs::WCSTransform{A}, pixels::AbstractMatrix, ::Type{T}) where {A,T<:AbstractFloat}
    _, ncoords = size(pixels)
    if wcs.sip === nothing
        # Apply CD to all pixel coordinates, then fold in the constant CRPIX offset.
        cd = Matrix{T}(wcs.cd)
        intermediate = cd * _typed_matrix(T, pixels)
        origin = wcs.cd * wcs.crpix
        for k in 1:ncoords, i in 1:A
            intermediate[i, k] -= T(origin[i])
        end
        return intermediate
    end

    offsets = Matrix{T}(undef, A, ncoords)

    # SIP must be applied column-wise before the batched CD transform.
    for k in 1:ncoords
        focal = sip_pixel_to_focal(wcs.sip, view(pixels, :, k))
        for i in 1:A
            offsets[i, k] = T(focal[i]) - T(wcs.crpix[i])
        end
    end
    return Matrix{T}(wcs.cd) * offsets
end

function _intermediate_to_pixel_batch(wcs::WCSTransform{A}, intermediate::AbstractMatrix, ::Type{T}) where {A,T<:AbstractFloat}
    # Solve the inverse linear transform for all coordinates at once.
    pixels = Matrix{T}(wcs.cd) \ Matrix{T}(intermediate)
    for k in axes(pixels, 2), i in 1:A
        pixels[i, k] += T(wcs.crpix[i])
    end

    # Apply inverse SIP column-wise because the polynomial solve is per coordinate.
    wcs.sip === nothing && return pixels
    result = similar(pixels, T, A, size(pixels, 2))
    for k in axes(pixels, 2)
        result[:, k] = sip_focal_to_pixel(wcs.sip, view(pixels, :, k))
    end
    return result
end

function _linear_world_to_pixel_batch(wcs::WCSTransform{A}, worlds::AbstractMatrix, ::Type{T}) where {A,T<:AbstractFloat}
    # Solve CD * p = world for all columns, then fold in the CRPIX/CRVAL offset.
    cd = Matrix{T}(wcs.cd)
    pixels = cd \ _typed_matrix(T, worlds)
    origin = convert(SVector{A,T}, wcs.crpix) - (cd \ convert(SVector{A,T}, wcs.crval))
    for k in axes(pixels, 2), i in 1:A
        pixels[i, k] += origin[i]
    end
    return pixels
end

function _add_crval_rows!(wcs::WCSTransform{A}, coords::AbstractMatrix, ::Type{T}) where {A,T<:AbstractFloat}
    # Add CRVAL by row so all non-celestial axes preserve their header offsets.
    for k in axes(coords, 2), i in 1:A
        coords[i, k] += T(wcs.crval[i])
    end
    return coords
end

function _world_offsets_batch(wcs::WCSTransform{A}, worlds::AbstractMatrix, ::Type{T}) where {A,T<:AbstractFloat}
    offsets = similar(worlds, T, A, size(worlds, 2))

    # Subtract CRVAL by row to form intermediate coordinates for all axes.
    for k in axes(worlds, 2), i in 1:A
        offsets[i, k] = T(worlds[i, k]) - T(wcs.crval[i])
    end
    return offsets
end

"""
    pixel_to_world(wcs, pixels::AbstractMatrix) -> AbstractMatrix

Batched pixel-to-world transform. Achieves better throughput than calling the single-coordinate transform repeatedly.

`pixels` must be an `naxis × N` matrix where each column is one pixel
coordinate.  Returns an `naxis × N` floating-point matrix of world coordinates.
"""
function pixel_to_world(wcs::WCSTransform, pixels::AbstractMatrix)
    naxis, N = size(pixels)
    naxis == wcs.naxis ||
        throw(DimensionMismatch("pixels has $(naxis) rows, expected $(wcs.naxis)"))
    T = _float_type(eltype(pixels))

    # Apply the linear pixel-to-intermediate transform to the whole batch.
    intermediate = _pixel_to_intermediate_batch(wcs, pixels, T)

    if isnothing(wcs.projection)
        # Purely linear WCS only needs CRVAL added to each output axis.
        return _add_crval_rows!(wcs, intermediate, T)
    end

    lon_idx = wcs.lon_axis
    lat_idx = wcs.lat_axis
    d2r = deg2rad(one(T))
    r2d = inv(d2r)
    alpha_p = T(wcs.alpha_p) * d2r
    delta_p = T(wcs.delta_p) * d2r
    phi_p = T(wcs.lonpole) * d2r

    # Preserve non-celestial axes, then overwrite celestial axes after projection.
    world = _add_crval_rows!(wcs, copy(intermediate), T)
    for k in 1:N
        phi, theta = intermediate_to_native(wcs.projection, intermediate[lon_idx, k], intermediate[lat_idx, k])
        alpha, delta = native_to_celestial(phi, theta, alpha_p, delta_p, phi_p)
        world[lon_idx, k] = mod(alpha * r2d, T(360))
        world[lat_idx, k] = delta * r2d
    end
    return world
end

"""
    pixel_to_world(wcs, pixels::AbstractVector{<:AbstractVector}) -> Vector

Convenience batch transform for a vector of individual pixel coordinates.
Simply executes an equivalent of `map(pixel -> pixel_to_world(wcs, pixel), pixels)`.
Total throughput is better for an `naxis × N` matrix of pixel coordinates, so use that form
when throughput for large batches matters.
"""
function pixel_to_world(wcs::WCSTransform, pixels::AbstractVector{<:AbstractVector})
    return pixel_to_world.(Ref(wcs), pixels)
end

"""
    world_to_pixel(wcs, worlds::AbstractMatrix) -> AbstractMatrix

Batched world-to-pixel transform. Achieves better throughput than calling the single-coordinate transform repeatedly.

`worlds` must be an `naxis × N` matrix where each column is one world
coordinate.  Returns an `naxis × N` floating-point matrix of pixel coordinates.
"""
function world_to_pixel(wcs::WCSTransform, worlds::AbstractMatrix)
    naxis, N = size(worlds)
    naxis == wcs.naxis ||
        throw(DimensionMismatch("worlds has $(naxis) rows, expected $(wcs.naxis)"))
    T = _float_type(eltype(worlds))

    if wcs.projection === nothing
        # Purely linear WCS can solve all world coordinates directly.
        return _linear_world_to_pixel_batch(wcs, worlds, T)
    end

    # Start from linear world offsets; celestial axes may be overwritten below.
    intermediate = _world_offsets_batch(wcs, worlds, T)

    lon_idx = wcs.lon_axis
    lat_idx = wcs.lat_axis
    d2r = deg2rad(one(T))
    alpha_p = T(wcs.alpha_p) * d2r
    delta_p = T(wcs.delta_p) * d2r
    phi_p = T(wcs.lonpole) * d2r

    # Convert each celestial coordinate back to projection-plane coordinates.
    for k in 1:N
        lon_delta = mod(T(worlds[lon_idx, k]) - T(wcs.crval[lon_idx]) + T(180), T(360)) - T(180)
        if abs(lon_delta) <= T(1e-10) && abs(T(worlds[lat_idx, k]) - T(wcs.crval[lat_idx])) <= T(1e-10)
            intermediate[lon_idx, k] = zero(T)
            intermediate[lat_idx, k] = zero(T)
        else
            alpha = T(worlds[lon_idx, k]) * d2r
            delta = T(worlds[lat_idx, k]) * d2r
            phi, theta = celestial_to_native(alpha, delta, alpha_p, delta_p, phi_p)
            x_lon, x_lat = native_to_intermediate(wcs.projection, phi, theta)
            intermediate[lon_idx, k] = T(x_lon)
            intermediate[lat_idx, k] = T(x_lat)
        end
    end

    # Solve the inverse linear transform for all coordinates in one operation.
    return _intermediate_to_pixel_batch(wcs, intermediate, T)
end

"""
    world_to_pixel(wcs, worlds::AbstractVector{<:AbstractVector}) -> Vector

Convenience batch transform for a vector of individual world coordinates.
Simply executes an equivalent of `map(world -> world_to_pixel(wcs, world), worlds)`.
Total throughput is better for an `naxis × N` matrix of world coordinates, so use that form
when throughput for large batches matters.
"""
function world_to_pixel(wcs::WCSTransform, worlds::AbstractVector{<:AbstractVector})
    return world_to_pixel.(Ref(wcs), worlds)
end
