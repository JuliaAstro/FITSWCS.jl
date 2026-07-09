module FITSWCSFITSFilesExt

import FITSFiles
import FITSWCS:
    WCS,
    WCS_all,
    NoAuxiliaryWCSData,
    _auxiliary_wcs_data,
    _external_auxiliary_data,
    _header_references_external_wcs_data,
    _lookup_table_from_image

"""
    _fitsfiles_parameter_card(value)
Parse a FITSFiles repeated Paper IV parameter-card payload into a `(name, value)` pair."""
function _fitsfiles_parameter_card(value)
    value isa AbstractString || return nothing
    parts = split(String(value), ':'; limit = 2)
    length(parts) == 2 || return nothing

    # Some readers expose HIERARCH Paper IV parameters as repeated string cards.
    name = uppercase(strip(parts[1]))
    raw = strip(parts[2])
    parsed = tryparse(Float64, raw)
    return name, isnothing(parsed) ? raw : parsed
end

"""
    _set_fitsfiles_header_card!(dict::Dict{String, Any}, key::String, value)
Store a FITSFiles header card and expand supported repeated parameter-card values."""
function _set_fitsfiles_header_card!(dict::Dict{String, Any}, key::String, value)
    # Preserve the normal keyword and expand parameter-card payloads when present.
    dict[key] = value
    parsed = _fitsfiles_parameter_card(value)
    isnothing(parsed) && return dict

    name, parsed_value = parsed
    if name == "EXTVER" || name == "NAXES" || startswith(name, "AXIS.")
        dict["$(key).$(name)"] = parsed_value
    end
    return dict
end

function _fitsfiles_cards_dict(cards::FITSFiles.Cards)
    # Copy card keyword/value pairs into the dictionary shape used by core parsing.
    dict = Dict{String, Any}()
    for card in cards
        _set_fitsfiles_header_card!(dict, String(card.key), card.value)
    end
    return dict
end

function _fitsfiles_auxiliary_wcs_data(header::AbstractDict, fobj; alt::Char = ' ', minerr::Real = 0.0)
    alt_str = alt == ' ' ? "" : string(alt)

    # Keep the common header-only path cheap even when a FITSFiles object is supplied.
    _header_references_external_wcs_data(header, alt_str) || return NoAuxiliaryWCSData()

    # Load referenced external arrays and copy them into backend-neutral tables.
    return _external_auxiliary_data(
        header,
        spec -> begin
            hdu = _fitsfiles_lookup_hdu(fobj, spec.extname, spec.extver)
            _lookup_table_from_image(hdu.data, _fitsfiles_cards_dict(hdu.cards), spec.transpose)
        end,
        (extname, extver, extlev, column) -> begin
            hdu = _fitsfiles_lookup_hdu(fobj, extname, extver)
            return _fitsfiles_table_column(hdu, column)
        end; alt = alt, minerr = minerr
    )
end

function _fitsfiles_lookup_hdu(hdus::AbstractVector{<:FITSFiles.HDU}, extname::AbstractString, extver::Integer)
    # FITSFiles indexes by EXTNAME only, so scan manually to respect EXTVER.
    for hdu in hdus
        cards = hdu.cards
        haskey(cards, "EXTNAME") || continue
        uppercase(rstrip(String(cards["EXTNAME"]))) == uppercase(String(extname)) || continue
        Int(get(cards, "EXTVER", 1)) == Int(extver) && return hdu
    end
    throw(KeyError((String(extname), Int(extver))))
end

function _fitsfiles_table_column(hdu::FITSFiles.HDU, column)
    data = hdu.data
    name = Symbol(column)

    # FITSFiles table HDUs commonly expose columns through named tuple data.
    if data isa NamedTuple && haskey(data, name)
        return data[name]
    elseif data isa AbstractDict && haskey(data, column)
        return data[column]
    elseif data isa AbstractDict && haskey(data, name)
        return data[name]
    end

    throw(KeyError(column))
end

function _auxiliary_wcs_data(header::AbstractDict, fobj::AbstractVector{<:FITSFiles.HDU}; alt::Char = ' ', minerr::Real = 0.0)
    # Route FITSFiles full-file HDU vectors through the extension-owned resolver.
    return _fitsfiles_auxiliary_wcs_data(header, fobj; alt = alt, minerr = minerr)
end

function WCS(cards::FITSFiles.Cards; fobj::Union{Nothing, AbstractVector{<:FITSFiles.HDU}} = nothing, alt::Char = ' ', minerr::Real = 0.0)
    # Delegate all WCS validation and interpretation to the core parser.
    return WCS(_fitsfiles_cards_dict(cards); fobj = fobj, alt = alt, minerr = minerr)
end

function WCS(hdu::FITSFiles.HDU; fobj::Union{Nothing, AbstractVector{<:FITSFiles.HDU}} = nothing, alt::Char = ' ', minerr::Real = 0.0)
    # FITSFiles stores parsed header cards directly on each HDU.
    return WCS(hdu.cards; fobj = fobj, alt = alt, minerr = minerr)
end

function WCS_all(cards::FITSFiles.Cards; fobj::Union{Nothing, AbstractVector{<:FITSFiles.HDU}} = nothing, minerr::Real = 0.0, preserve_units::Bool = false)
    # Convert FITSFiles cards to Dict and delegate to the core parser.
    return WCS_all(_fitsfiles_cards_dict(cards); fobj, minerr, preserve_units)
end

function WCS_all(hdu::FITSFiles.HDU; fobj::Union{Nothing, AbstractVector{<:FITSFiles.HDU}} = nothing, minerr::Real = 0.0, preserve_units::Bool = false)
    # FITSFiles stores parsed header cards directly on each HDU.
    return WCS_all(hdu.cards; fobj, minerr, preserve_units)
end

end # module FITSWCSFITSFilesExt
