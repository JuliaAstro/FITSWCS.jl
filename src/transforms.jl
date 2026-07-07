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

## World coordinate convention

Return values from `pixel_to_world` and expected input values to
`world_to_pixel` follow these rules:

- **Celestial axes** (RA/DEC, GLON/GLAT, etc.): **degrees**.
- **Spectral axes** (FREQ, WAVE, VELO, ENER, AWAV, etc.): **SI units**
  (Hz, m, m/s, J, rad/s, or dimensionless for ZOPT/BETA).
- **Other linear axes**: whatever units the CD matrix encodes (typically the
  header's original CUNIT).

When `preserve_units=true` was set at construction, the values are scaled
back to the original header CUNIT for celestial and spectral axes.  See
`WCS()` for details.

The output vector element at `lon_axis` is right ascension/longitude and at
`lat_axis` is declination/latitude.
"""

# ──────────────────────────────────────────────────────────────────────────────
# Single-coordinate transforms
# ──────────────────────────────────────────────────────────────────────────────

@inline _coordinate_float_type(::AbstractVector{T}) where {T <: Real} = _float_type(T)

@inline _coordinate_float_type(coords::Tuple{Vararg{Real}}) = _promote_float_type(coords...)

function _coordinate_vector(coords::Tuple{Vararg{Real, N}}) where {N}
    # Preserve tuple precision and length when materializing scalar arguments.
    T = _coordinate_float_type(coords)
    return SVector{N, T}(coords)
end

@inline _tabular_data(::NoAuxiliaryWCSData) = NoTabularWCSData()
@inline _tabular_data(aux::AuxiliaryWCSData) = aux.tabular
@inline _spectral_data(::NoAuxiliaryWCSData) = NoSpectralWCSData()
@inline _spectral_data(aux::AuxiliaryWCSData) = aux.spectral
@inline _grism_data(::NoAuxiliaryWCSData) = NoGrismWCSData()
@inline _grism_data(aux::AuxiliaryWCSData) = aux.grism

function intermediate_to_tabular_world(::NoTabularWCSData, wcs::WCSTransform{N}, intermediate::AbstractVector, ::Type{T}) where {N, T <: AbstractFloat}
    spectral_data = _spectral_data(wcs.aux)
    grism_data = _grism_data(wcs.aux)
    if spectral_data isa NoSpectralWCSData && grism_data isa NoGrismWCSData
        # Fast path: no TAB, no spectral, no grism -- just add CRVAL.
        return SVector{N, T}(ntuple(i -> T(wcs.crval[i]) + T(intermediate[i]), N))
    end
    # Build a mutable world vector, apply spectral/grism conversion, then freeze.
    world = MVector{N, T}(undef)
    @inbounds for i in 1:N
        world[i] = T(wcs.crval[i]) + T(intermediate[i])
    end
    _apply_spectral_forward!(spectral_data, world, intermediate, T)
    _apply_grism_forward!(grism_data, world, intermediate, T)
    return SVector{N, T}(world)
end

function intermediate_to_tabular_world(tabular::TabularWCSData, wcs::WCSTransform{N}, intermediate::AbstractVector, ::Type{T}) where {N, T <: AbstractFloat}
    world = MVector{N, T}(undef)

    # Start with ordinary linear axes, then replace TAB-controlled axes.
    @inbounds for i in 1:N
        world[i] = T(wcs.crval[i]) + T(intermediate[i])
    end
    # Apply spectral conversion before TAB so that TAB axes (which are
    # never also spectral) overwrite the linear baseline correctly.
    _apply_spectral_forward!(_spectral_data(wcs.aux), world, intermediate, T)
    _apply_grism_forward!(_grism_data(wcs.aux), world, intermediate, T)
    @inbounds for i in eachindex(tabular.tables)
        table = tabular.tables[i]
        values = _tabular_forward(table, intermediate, wcs.crval)
        for j in eachindex(table.axes)
            axis = table.axes[j]
            world[axis] = T(values[j])
        end
    end

    return SVector{N, T}(world)
end

function tabular_world_to_intermediate(::NoTabularWCSData, wcs::WCSTransform{N}, world::AbstractVector, ::Type{T}) where {N, T <: AbstractFloat}
    spectral_data = _spectral_data(wcs.aux)
    grism_data = _grism_data(wcs.aux)
    if spectral_data isa NoSpectralWCSData && grism_data isa NoGrismWCSData
        # Fast path: no TAB, no spectral, no grism -- just subtract CRVAL.
        return SVector{N, T}(ntuple(i -> T(world[i]) - T(wcs.crval[i]), N))
    end
    intermediate = MVector{N, T}(undef)
    @inbounds for i in 1:N
        intermediate[i] = T(world[i]) - T(wcs.crval[i])
    end
    _apply_spectral_inverse!(spectral_data, intermediate, world, T)
    _apply_grism_inverse!(grism_data, intermediate, world, T)
    return SVector{N, T}(intermediate)
end

function tabular_world_to_intermediate(tabular::TabularWCSData, wcs::WCSTransform{N}, world::AbstractVector, ::Type{T}) where {N, T <: AbstractFloat}
    intermediate = MVector{N, T}(undef)

    # Start with ordinary linear axes, then replace TAB-controlled axes.
    @inbounds for i in 1:N
        intermediate[i] = T(world[i]) - T(wcs.crval[i])
    end
    # Apply spectral/grism inverse before TAB.
    _apply_spectral_inverse!(_spectral_data(wcs.aux), intermediate, world, T)
    _apply_grism_inverse!(_grism_data(wcs.aux), intermediate, world, T)
    @inbounds for i in eachindex(tabular.tables)
        table = tabular.tables[i]
        values = _tabular_inverse(table, world, wcs.crval)
        for j in eachindex(table.axes)
            axis = table.axes[j]
            intermediate[axis] = T(values[j])
        end
    end

    return SVector{N, T}(intermediate)
end

# ── Spectral pipeline helpers ─────────────────────────────────────────────────

function _apply_spectral_forward!(::NoSpectralWCSData, world, intermediate, ::Type{T}) where {T}
    return world  # no spectral axes — zero overhead
end

function _apply_spectral_forward!(spec::SpectralWCSData, world, intermediate, ::Type{T}) where {T}
    @inbounds for s in spec.specs
        _is_linear(s) && continue
        world[s.axis] = T(_spectral_x_to_world(T(intermediate[s.axis]), s))
    end
    return world
end

function _apply_spectral_inverse!(::NoSpectralWCSData, intermediate, world, ::Type{T}) where {T}
    return intermediate  # no spectral axes — zero overhead
end

function _apply_spectral_inverse!(spec::SpectralWCSData, intermediate, world, ::Type{T}) where {T}
    @inbounds for s in spec.specs
        _is_linear(s) && continue
        intermediate[s.axis] = T(_spectral_world_to_x(T(world[s.axis]), s))
    end
    return intermediate
end

# ── Grism pipeline helpers ────────────────────────────────────────────────────

function _apply_grism_forward!(::NoGrismWCSData, world, intermediate, ::Type{T}) where {T}
    return world  # no grism axes -- zero overhead
end

function _apply_grism_forward!(grism::GrismWCSData, world, intermediate, ::Type{T}) where {T}
    @inbounds for g in grism.specs
        world[g.axis] = T(_grism_x_to_world(T(intermediate[g.axis]), g))
    end
    return world
end

function _apply_grism_inverse!(::NoGrismWCSData, intermediate, world, ::Type{T}) where {T}
    return intermediate  # no grism axes -- zero overhead
end

function _apply_grism_inverse!(grism::GrismWCSData, intermediate, world, ::Type{T}) where {T}
    @inbounds for g in grism.specs
        intermediate[g.axis] = T(_grism_world_to_x(T(world[g.axis]), g))
    end
    return intermediate
end

function _set_celestial_axes!(coords::AbstractVector, lon_idx::Int, lat_idx::Int, lon, lat)
    # Mutable coordinate storage can be updated in place.
    coords[lon_idx] = lon
    coords[lat_idx] = lat
    return coords
end

function _set_celestial_axes!(coords::StaticVector{N, T}, lon_idx::Int, lat_idx::Int, lon, lat) where {N, T}
    # Immutable static coordinates are rebuilt with the celestial axes replaced.
    return SVector{N, T}(ntuple(i -> i == lon_idx ? T(lon) : i == lat_idx ? T(lat) : coords[i], N))
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

        # Deproject: (x, y) -> (phi, theta) in radians
        phi, theta = intermediate_to_native(wcs.projection, x_lon, x_lat)

        # Spherical rotation: (phi, theta) -> (alpha, delta) in radians
        # alpha_p, delta_p are the celestial coords of the native north pole.
        alpha_p = deg2rad(T(wcs.alpha_p))
        delta_p = deg2rad(T(wcs.delta_p))
        phi_p = deg2rad(T(wcs.lonpole))

        alpha, delta = native_to_celestial(phi, theta, alpha_p, delta_p, phi_p)
        # Build world vector: start from the intermediate coords, then overwrite
        # the celestial axes with the spherical-rotation results.
        world = intermediate_to_tabular_world(_tabular_data(wcs.aux), wcs, x, T)
        world = _set_celestial_axes!(world, lon_idx, lat_idx, mod(rad2deg(alpha), 360), rad2deg(delta))
    else
        # Purely linear: world = CRVAL + CD*(pixel - CRPIX)
        world = intermediate_to_tabular_world(_tabular_data(wcs.aux), wcs, x, T)
    end
    # Scale world coordinates back to CUNIT if preserve_units is set.
    # Convert unit_scaling to the world's element type so the broadcast
    # preserves type stability (Float32 in, Float32 out).
    if wcs.preserve_units
        return world ./ T.(wcs.unit_scaling)
    end
    return world
end

# Convenience: accept tuples, static-array-likes, and scalar varargs.
pixel_to_world(wcs::WCSTransform, pixel::Tuple{Vararg{Real}}) =
    pixel_to_world(wcs, _coordinate_vector(pixel))
pixel_to_world(wcs::WCSTransform, coords::Real...) =
    pixel_to_world(wcs, _coordinate_vector(coords))

"""
    world_to_pixel(wcs, world) -> pixel

Convert world coordinates to FITS pixel coordinates.

`world` must be a length-`wcs.naxis` vector (or tuple) of world positions
in the same units as `CRVAL` (degrees for celestial axes).
Returns a length-`wcs.naxis` floating-point vector of 1-based pixel coordinates.
Tuple inputs and scalar varargs are materialized as static coordinates and return
`SVector` results.
"""
function world_to_pixel(wcs::WCSTransform{N}, world::AbstractVector) where {N}
    length(world) == wcs.naxis ||
        throw(DimensionMismatch("world has length $(length(world)), expected $(wcs.naxis)"))
    T = _coordinate_float_type(world)

    # Scale world from CUNIT to canonical units if preserve_units is set.
    world_canon = if wcs.preserve_units
        SVector{N, T}(world) .* T.(wcs.unit_scaling)
    else
        SVector{N,T}(world)
    end

    if wcs.projection !== nothing
        lon_idx = wcs.lon_axis
        lat_idx = wcs.lat_axis

        # The fiducial celestial point is defined to have zero intermediate
        # coordinates; handling it directly avoids pole cancellation.
        btol = _boundary_tol(T)
        lon_delta = mod(world_canon[lon_idx] - T(wcs.crval[lon_idx]) + 180, 360) - 180
        if abs(lon_delta) <= btol && abs(world_canon[lat_idx] - T(wcs.crval[lat_idx])) <= btol
            x = tabular_world_to_intermediate(_tabular_data(wcs.aux), wcs, world_canon, T)
            x = _set_celestial_axes!(x, lon_idx, lat_idx, zero(T), zero(T))
            return intermediate_to_pixel(wcs, x)
        end

        alpha_p = deg2rad(T(wcs.alpha_p))
        delta_p = deg2rad(T(wcs.delta_p))
        phi_p = deg2rad(T(wcs.lonpole))

        alpha = deg2rad(T(world_canon[lon_idx]))
        delta = deg2rad(T(world_canon[lat_idx]))

        # Spherical rotation: (alpha, delta) -> (phi, theta)
        phi, theta = celestial_to_native(alpha, delta, alpha_p, delta_p, phi_p)

        # Re-project: (phi, theta) -> (x, y) in degrees
        x_lon, x_lat = native_to_intermediate(wcs.projection, phi, theta)

        # Build intermediate coordinate vector
        # Non-celestial axes: x_i = world_i - crval_i (trivial linear axes)
        x = tabular_world_to_intermediate(_tabular_data(wcs.aux), wcs, world_canon, T)
        x = _set_celestial_axes!(x, lon_idx, lat_idx, T(x_lon), T(x_lat))

        return intermediate_to_pixel(wcs, x)
    else
        x = tabular_world_to_intermediate(_tabular_data(wcs.aux), wcs, world_canon, T)
        return intermediate_to_pixel(wcs, x)
    end
end

world_to_pixel(wcs::WCSTransform, world::Tuple{Vararg{Real}}) =
    world_to_pixel(wcs, _coordinate_vector(world))
world_to_pixel(wcs::WCSTransform, coords::Real...) =
    world_to_pixel(wcs, _coordinate_vector(coords))

# ──────────────────────────────────────────────────────────────────────────────
# Batch transforms
# ──────────────────────────────────────────────────────────────────────────────

@inline _typed_matrix(::Type{T}, values::AbstractMatrix{T}) where {T <: AbstractFloat} = values
@inline _typed_matrix(::Type{T}, values::AbstractMatrix) where {T <: AbstractFloat} = Matrix{T}(values)

function _pixel_to_intermediate_batch(wcs::WCSTransform{A}, pixels::AbstractMatrix, ::Type{T}) where {A, T <: AbstractFloat}
    _, ncoords = size(pixels)
    if !has_distortion(wcs.pipeline)
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

    # Distortions are applied column-wise before the batched CD transform.
    for k in 1:ncoords
        focal = pixel_to_focal(wcs.pipeline, view(pixels, :, k), Val(A))
        for i in 1:A
            offsets[i, k] = T(focal[i]) - T(wcs.crpix[i])
        end
    end
    return Matrix{T}(wcs.cd) * offsets
end

function _intermediate_to_pixel_batch(wcs::WCSTransform{A}, intermediate::AbstractMatrix, ::Type{T}) where {A, T <: AbstractFloat}
    # Solve the inverse linear transform for all coordinates at once.
    pixels = Matrix{T}(wcs.cd) \ Matrix{T}(intermediate)
    for k in axes(pixels, 2), i in 1:A
        pixels[i, k] += T(wcs.crpix[i])
    end

    # Apply inverse distortions column-wise because the solve is per coordinate.
    !has_distortion(wcs.pipeline) && return pixels
    result = similar(pixels, T, A, size(pixels, 2))
    for k in axes(pixels, 2)
        result[:, k] = focal_to_pixel(wcs.pipeline, view(pixels, :, k), Val(A))
    end
    return result
end

function _linear_world_to_pixel_batch(wcs::WCSTransform{A}, worlds::AbstractMatrix, ::Type{T}) where {A, T <: AbstractFloat}
    # Solve CD * p = world for all columns, then fold in the CRPIX/CRVAL offset.
    cd = Matrix{T}(wcs.cd)
    pixels = cd \ _typed_matrix(T, worlds)
    origin = convert(SVector{A, T}, wcs.crpix) - (cd \ convert(SVector{A, T}, wcs.crval))
    for k in axes(pixels, 2), i in 1:A
        pixels[i, k] += origin[i]
    end
    return pixels
end

function _add_crval_rows!(wcs::WCSTransform{A}, coords::AbstractMatrix, ::Type{T}) where {A, T <: AbstractFloat}
    # Add CRVAL by row so all non-celestial axes preserve their header offsets.
    for k in axes(coords, 2), i in 1:A
        coords[i, k] += T(wcs.crval[i])
    end
    return coords
end

function _world_offsets_batch(wcs::WCSTransform{A}, worlds::AbstractMatrix, ::Type{T}) where {A, T <: AbstractFloat}
    offsets = similar(worlds, T, A, size(worlds, 2))

    # Subtract CRVAL by row to form intermediate coordinates for all axes.
    for k in axes(worlds, 2), i in 1:A
        offsets[i, k] = T(worlds[i, k]) - T(wcs.crval[i])
    end
    return offsets
end

function _intermediate_to_world_batch(wcs::WCSTransform{A}, intermediate::AbstractMatrix, ::Type{T}) where {A, T <: AbstractFloat}
    tabular = _tabular_data(wcs.aux)
    if !has_tabular(tabular)
        return _add_crval_rows!(wcs, intermediate, T)
    end

    # TAB interpolation is per coordinate, so evaluate it column-wise.
    world = similar(intermediate, T, A, size(intermediate, 2))
    for k in axes(intermediate, 2)
        world[:, k] = intermediate_to_tabular_world(tabular, wcs, view(intermediate, :, k), T)
    end
    return world
end

function _world_to_intermediate_batch(wcs::WCSTransform{A}, worlds::AbstractMatrix, ::Type{T}) where {A, T <: AbstractFloat}
    tabular = _tabular_data(wcs.aux)
    if !has_tabular(tabular)
        return _world_offsets_batch(wcs, worlds, T)
    end

    # TAB inverse is per coordinate, so evaluate it column-wise.
    intermediate = similar(worlds, T, A, size(worlds, 2))
    for k in axes(worlds, 2)
        intermediate[:, k] = tabular_world_to_intermediate(tabular, wcs, view(worlds, :, k), T)
    end
    return intermediate
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
        world = _intermediate_to_world_batch(wcs, intermediate, T)
    else
        lon_idx = wcs.lon_axis
        lat_idx = wcs.lat_axis
        alpha_p = deg2rad(T(wcs.alpha_p))
        delta_p = deg2rad(T(wcs.delta_p))
        phi_p = deg2rad(T(wcs.lonpole))

        # Preserve non-celestial axes, then overwrite celestial axes after projection.
        world = _intermediate_to_world_batch(wcs, copy(intermediate), T)
        for k in 1:N
            phi, theta = intermediate_to_native(wcs.projection, intermediate[lon_idx, k], intermediate[lat_idx, k])
            alpha, delta = native_to_celestial(phi, theta, alpha_p, delta_p, phi_p)
            world[lon_idx, k] = mod(rad2deg(alpha), 360)
            world[lat_idx, k] = rad2deg(delta)
        end
    end
    # Scale world coordinates back to CUNIT if preserve_units is set.
    # Convert unit_scaling to the world's element type for type stability.
    if wcs.preserve_units
        for i in 1:wcs.naxis
            world[i, :] ./= T(wcs.unit_scaling[i])
        end
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

    # Scale world coordinates from CUNIT to canonical if preserve_units is set.
    worlds_canon = if wcs.preserve_units
        ws = similar(worlds, T, naxis, N)
        for i in 1:naxis
            ws[i, :] .= worlds[i, :] .* T(wcs.unit_scaling[i])
        end
        ws
    else
        worlds
    end

    if wcs.projection === nothing
        if has_distortion(wcs.pipeline)
            # Distorted linear WCS still needs the per-coordinate inverse pipeline.
            intermediate = _world_to_intermediate_batch(wcs, worlds_canon, T)
            return _intermediate_to_pixel_batch(wcs, intermediate, T)
        end

        # Purely linear, undistorted WCS can solve all world coordinates directly.
        has_tabular(_tabular_data(wcs.aux)) &&
            return _intermediate_to_pixel_batch(wcs, _world_to_intermediate_batch(wcs, worlds_canon, T), T)
        return _linear_world_to_pixel_batch(wcs, worlds_canon, T)
    end

    # Start from linear world offsets; celestial axes may be overwritten below.
    intermediate = _world_to_intermediate_batch(wcs, worlds_canon, T)

    lon_idx = wcs.lon_axis
    lat_idx = wcs.lat_axis
    btol = _boundary_tol(T)
    alpha_p = deg2rad(T(wcs.alpha_p))
    delta_p = deg2rad(T(wcs.delta_p))
    phi_p = deg2rad(T(wcs.lonpole))

    # Convert each celestial coordinate back to projection-plane coordinates.
    for k in 1:N
        lon_delta = mod(T(worlds_canon[lon_idx, k]) - T(wcs.crval[lon_idx]) + T(180), T(360)) - T(180)
        if abs(lon_delta) <= btol && abs(T(worlds_canon[lat_idx, k]) - T(wcs.crval[lat_idx])) <= btol
            intermediate[lon_idx, k] = zero(T)
            intermediate[lat_idx, k] = zero(T)
        else
            alpha = deg2rad(T(worlds_canon[lon_idx, k]))
            delta = deg2rad(T(worlds_canon[lat_idx, k]))
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
