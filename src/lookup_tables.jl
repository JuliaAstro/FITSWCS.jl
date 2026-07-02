"""
Lookup-table interpolation utilities for external WCS data.
"""

"""
    LookupTable2D

Backend-independent two-dimensional Paper IV lookup table with FITS-style
linear index metadata.
"""
struct LookupTable2D{T <: Real, A <: AbstractMatrix{T}}
    data::A
    crpix::SVector{2, T}
    crval::SVector{2, T}
    cdelt::SVector{2, T}
end

function LookupTable2D(
        data::A;
        crpix = (one(T), one(T)),
        crval = (zero(T), zero(T)),
        cdelt = (one(T), one(T)),
    ) where {T <: Real, A <: AbstractMatrix{T}}
    isempty(data) && throw(ArgumentError("lookup table data must not be empty"))

    # Store table WCS metadata in the same element type as the table data.
    return LookupTable2D{T, A}(
        data,
        SVector{2, T}(crpix),
        SVector{2, T}(crval),
        SVector{2, T}(cdelt),
    )
end

function _lookup_fractional_index(table::LookupTable2D, x::Real, y::Real)
    R = promote_type(float(eltype(table.data)), _promote_float_type(x, y))

    # Convert world/pixel table coordinates into FITS-style fractional indices.
    coord = SVector{2, R}(x, y)
    crpix = SVector{2, R}(table.crpix)
    crval = SVector{2, R}(table.crval)
    cdelt = SVector{2, R}(table.cdelt)
    return (coord .- crval) ./ cdelt .+ crpix
end

function _lookup_cell(index::Real, n::Int)
    R = _promote_float_type(index)
    n > 0 || throw(ArgumentError("lookup table axes must not be empty"))

    # Singleton axes have no interpolation span, so the only cell is exact.
    n == 1 && return 1, 1, zero(R)

    # Clamp to the table domain and choose the lower corner of the containing cell.
    clamped = clamp(R(index), one(R), R(n))
    lo = min(floor(Int, clamped), n - 1)
    hi = lo + 1
    weight = clamped - R(lo)
    return lo, hi, weight
end

"""
    interpolate_lookup_table(table, x, y)

Return the bilinearly interpolated value from `table` at table coordinates
`(x, y)`.  Coordinates are converted through `crpix`, `crval`, and `cdelt`,
then clamped to the valid table domain.
"""
function interpolate_lookup_table(table::LookupTable2D, x::Real, y::Real)
    index = _lookup_fractional_index(table, x, y)

    # Locate the interpolation cell in each matrix dimension.
    i0, i1, ti = _lookup_cell(index[1], size(table.data, 1))
    j0, j1, tj = _lookup_cell(index[2], size(table.data, 2))

    # Blend the four cell corners, degenerating cleanly for singleton axes.
    v00 = table.data[i0, j0]
    v10 = table.data[i1, j0]
    v01 = table.data[i0, j1]
    v11 = table.data[i1, j1]
    return (1 - ti) * (1 - tj) * v00 +
        ti * (1 - tj) * v10 +
        (1 - ti) * tj * v01 +
        ti * tj * v11
end

function interpolate_lookup_table(table::LookupTable2D, coord::StaticVector{2, <:Real})
    # Allow distortion code to pass fixed-size coordinate vectors directly.
    return interpolate_lookup_table(table, coord[1], coord[2])
end

function (table::LookupTable2D)(x::Real, y::Real)
    # Treat lookup tables as coordinate-to-value functions at call sites.
    return interpolate_lookup_table(table, x, y)
end

function (table::LookupTable2D)(coord::StaticVector{2, <:Real})
    # Preserve the same callable interface for fixed-size coordinate vectors.
    return interpolate_lookup_table(table, coord)
end
