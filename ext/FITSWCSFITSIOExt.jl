module FITSWCSFITSIOExt

import FITSIO
import FITSWCS:
    WCS,
    NoAuxiliaryWCSData,
    _auxiliary_wcs_data,
    _header_references_external_wcs_data,
    _lookup_table_from_image,
    _paper_iv_auxiliary_data

"""
    _fitsio_parameter_card(value)
Parse a FITSIO repeated Paper IV parameter-card payload into a `(name, value)` pair."""
function _fitsio_parameter_card(value)
    value isa AbstractString || return nothing
    parts = split(String(value), ':'; limit = 2)
    length(parts) == 2 || return nothing

    # FITSIO stores HIERARCH Paper IV parameters as repeated string cards.
    name = uppercase(strip(parts[1]))
    raw = strip(parts[2])
    parsed = tryparse(Float64, raw)
    return name, isnothing(parsed) ? raw : parsed
end

"""
    _set_fitsio_header_card!(dict::Dict{String, Any}, key::String, value)
Store a FITSIO header card and expand supported repeated parameter-card values."""
function _set_fitsio_header_card!(dict::Dict{String, Any}, key::String, value)
    # Preserve the normal keyword and expand parameter-card payloads when present.
    dict[key] = value
    parsed = _fitsio_parameter_card(value)
    isnothing(parsed) && return dict

    name, parsed_value = parsed
    if name == "EXTVER" || name == "NAXES" || startswith(name, "AXIS.")
        dict["$(key).$(name)"] = parsed_value
    end
    return dict
end

function _fitsio_header_dict(header::FITSIO.FITSHeader)
    # Copy keyword/value pairs into the dictionary shape used by core parsing.
    dict = Dict{String, Any}()
    for (key, value) in zip(keys(header), values(header))
        _set_fitsio_header_card!(dict, String(key), value)
    end
    return dict
end

function _fitsio_auxiliary_wcs_data(header::AbstractDict, fobj; alt::Char = ' ', minerr::Real = 0.0)
    alt_str = alt == ' ' ? "" : string(alt)

    # Keep the common header-only path cheap even when a FITSIO object is supplied.
    _header_references_external_wcs_data(header, alt_str) || return NoAuxiliaryWCSData()

    # Load referenced Paper IV image extensions and copy them into backend-neutral tables.
    return _paper_iv_auxiliary_data(
        header, spec -> begin
            hdu = fobj[spec.extname, spec.extver]
            data = FITSIO.read(hdu)
            table_header = _fitsio_header_dict(FITSIO.read_header(hdu))
            _lookup_table_from_image(data, table_header, spec.transpose)
        end; alt = alt, minerr = minerr
    )
end

function _auxiliary_wcs_data(header::AbstractDict, fobj::FITSIO.FITS; alt::Char = ' ', minerr::Real = 0.0)
    # Route FITSIO file containers through the extension-owned resolver.
    return _fitsio_auxiliary_wcs_data(header, fobj; alt = alt, minerr = minerr)
end

function _auxiliary_wcs_data(header::AbstractDict, fobj::FITSIO.HDU; alt::Char = ' ', minerr::Real = 0.0)
    # Allow callers to pass a FITSIO HDU when no file-level container is available.
    return _fitsio_auxiliary_wcs_data(header, fobj; alt = alt, minerr = minerr)
end

function WCS(header::FITSIO.FITSHeader; fobj = nothing, alt::Char = ' ', minerr::Real = 0.0)
    # Delegate all WCS validation and interpretation to the core parser.
    return WCS(_fitsio_header_dict(header); fobj = fobj, alt = alt, minerr = minerr)
end

function WCS(hdu::FITSIO.HDU; fobj = nothing, alt::Char = ' ', minerr::Real = 0.0)
    # Read the HDU header through FITSIO before using the header adapter.
    return WCS(FITSIO.read_header(hdu); fobj = fobj, alt = alt, minerr = minerr)
end

end # module FITSWCSFITSIOExt
