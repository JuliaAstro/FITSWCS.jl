"""
Small public API compatibility helpers.

These names mirror common WCS.jl workflows, but keep FITSWCS.jl's FITS
1-based pixel coordinate convention.
"""

"""
    WCSTransform(naxis; kwds...) -> WCSTransform

Construct a transform from WCS.jl-style keyword vectors and matrices.

Supported keywords are `crpix`, `crval`, `cdelt`, `ctype`, `cunit`, `pc`, `cd`,
`crota`, `lonpole`, and `latpole`.  More specialized wcslib fields remain
unsupported in this pure-Julia constructor.
"""
function WCSTransform(naxis::Integer; kwds...)
    naxis >= 1 || throw(ArgumentError("naxis must be >= 1, got $naxis"))

    # Translate the supported property-style inputs into ordinary FITS keys.
    header = Dict{String,Any}("NAXIS" => Int(naxis), "WCSAXES" => Int(naxis))
    for (key, value) in kwds
        _constructor_keyword_to_header!(header, Int(naxis), key, value)
    end

    # Reuse the main parser so validation and projection setup stay centralized.
    return WCS(header)
end

function _constructor_keyword_to_header!(header::Dict{String,Any},
                                         naxis::Int,
                                         key::Symbol,
                                         value)
    # Vector-valued WCS properties map directly to indexed FITS keywords.
    if key === :crpix
        _put_axis_values!(header, "CRPIX", value, naxis)
    elseif key === :crval
        _put_axis_values!(header, "CRVAL", value, naxis)
    elseif key === :cdelt
        _put_axis_values!(header, "CDELT", value, naxis)
    elseif key === :ctype
        _put_axis_values!(header, "CTYPE", value, naxis)
    elseif key === :cunit
        _put_axis_values!(header, "CUNIT", value, naxis)
    elseif key === :crota
        _put_axis_values!(header, "CROTA", value, naxis)
    elseif key === :pc
        _put_matrix_values!(header, "PC", value, naxis)
    elseif key === :cd
        _put_matrix_values!(header, "CD", value, naxis)
    elseif key === :lonpole
        header["LONPOLE"] = value
    elseif key === :latpole
        header["LATPOLE"] = value
    else
        throw(ArgumentError("unsupported WCSTransform constructor keyword: $key"))
    end

    return header
end

function _put_axis_values!(header::Dict{String,Any},
                           prefix::AbstractString,
                           values,
                           naxis::Int)
    length(values) == naxis ||
        throw(DimensionMismatch("$prefix values have length $(length(values)), expected $naxis"))

    # Store each value under the matching FITS numbered keyword.
    for i in 1:naxis
        header["$(prefix)$(i)"] = values[i]
    end

    return header
end

function _put_matrix_values!(header::Dict{String,Any},
                             prefix::AbstractString,
                             values,
                             naxis::Int)
    size(values) == (naxis, naxis) ||
        throw(DimensionMismatch("$prefix matrix has size $(size(values)), expected ($naxis, $naxis)"))

    # Store matrices in FITS world-axis by pixel-axis order.
    for i in 1:naxis, j in 1:naxis
        header["$(prefix)$(i)_$(j)"] = values[i, j]
    end

    return header
end

"""
    pix_to_world(wcs, pixels)
    pix_to_world(wcs, x, y, ...)

Compatibility alias for `pixel_to_world`.
"""
function pix_to_world(wcs::WCSTransform, pixels::Union{AbstractVector, AbstractMatrix, Tuple})
    # Delegate to the canonical transform implementation.
    return pixel_to_world(wcs, pixels)
end

function pix_to_world(wcs::WCSTransform, coords::Real...)
    # Materialize scalar coordinate arguments as static coordinates.
    return pixel_to_world(wcs, _coordinate_vector(coords))
end

"""
    pix_to_world!(wcs, pixels, worlds)

Mutating compatibility alias for `pixel_to_world`.
"""
function pix_to_world!(wcs::WCSTransform,
                       pixels::Union{AbstractVector, AbstractMatrix},
                       worlds::Union{AbstractVector, AbstractMatrix})
    size(worlds) == size(pixels) ||
        throw(DimensionMismatch("worlds has size $(size(worlds)), expected $(size(pixels))"))

    # Compute through the canonical API, then copy into the caller's output.
    worlds .= pixel_to_world(wcs, pixels)
    return worlds
end

"""
    world_to_pix(wcs, worlds)
    world_to_pix(wcs, x, y, ...)

Compatibility alias for `world_to_pixel`.
"""
function world_to_pix(wcs::WCSTransform, worlds::Union{AbstractVector, AbstractMatrix, Tuple})
    # Delegate to the canonical transform implementation.
    return world_to_pixel(wcs, worlds)
end

function world_to_pix(wcs::WCSTransform, coords::Real...)
    # Materialize scalar coordinate arguments as static coordinates.
    return world_to_pixel(wcs, _coordinate_vector(coords))
end

"""
    world_to_pix!(wcs, worlds, pixels)

Mutating compatibility alias for `world_to_pixel`.
"""
function world_to_pix!(wcs::WCSTransform,
                       worlds::Union{AbstractVector, AbstractMatrix},
                       pixels::Union{AbstractVector, AbstractMatrix})
    size(pixels) == size(worlds) ||
        throw(DimensionMismatch("pixels has size $(size(pixels)), expected $(size(worlds))"))

    # Compute through the canonical API, then copy into the caller's output.
    pixels .= world_to_pixel(wcs, worlds)
    return pixels
end
