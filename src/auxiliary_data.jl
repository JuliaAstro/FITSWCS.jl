"""
Backend-agnostic auxiliary WCS data resolution.

The main package owns the data model and fallback behavior.  FITS backend
extensions add methods for their HDU-list/container types.
"""

"""Abstract supertype for backend-independent auxiliary WCS data payloads."""
abstract type AbstractAuxiliaryWCSData end

"""Auxiliary-data payload for WCS transforms whose headers need no external data."""
struct NoAuxiliaryWCSData <: AbstractAuxiliaryWCSData end

"""
    AuxiliaryWCSData

Backend-independent external WCS data resolved at construction time.
"""
struct AuxiliaryWCSData{D, C, T} <: AbstractAuxiliaryWCSData
    det2im::D
    cpdis::C
    tabular::T
end

AuxiliaryWCSData(; det2im = (nothing, nothing), cpdis = (nothing, nothing), tabular = nothing) =
    AuxiliaryWCSData(det2im, cpdis, tabular)

has_auxiliary(::NoAuxiliaryWCSData) = false
has_auxiliary(::AbstractAuxiliaryWCSData) = true

function _header_references_tabular_axis(header::AbstractDict, alt_str::AbstractString)
    suffix = uppercase(String(alt_str))

    # Scan only CTYPE keywords for the selected WCS alternate.
    for (key, value) in header
        key isa AbstractString || continue
        occursin(Regex("^CTYPE[1-9][0-9]*$(suffix)\$"), uppercase(String(key))) || continue
        occursin("-TAB", uppercase(String(value))) && return true
    end

    return false
end

function _header_references_external_wcs_data(header::AbstractDict, alt_str::AbstractString)
    # Paper IV distortion keywords and TAB axes are the current external-data signals.
    for key in keys(header)
        key isa AbstractString || continue
        is_lookup_distortion_keyword(key, alt_str) && return true
    end

    return _header_references_tabular_axis(header, alt_str)
end

function _auxiliary_wcs_data(header::AbstractDict, ::Nothing; alt::Char = ' ', minerr::Real = 0.0)
    alt_str = alt == ' ' ? "" : string(alt)

    # Header-only construction is valid unless the selected WCS references external arrays.
    if _header_references_external_wcs_data(header, alt_str)
        throw(ArgumentError("external WCS data are required; pass WCS(header; fobj=...)"))
    end

    return NoAuxiliaryWCSData()
end

function _auxiliary_wcs_data(header::AbstractDict, fobj; alt::Char = ' ', minerr::Real = 0.0)
    alt_str = alt == ' ' ? "" : string(alt)

    # Unknown fobj types are harmless for ordinary header-only WCS metadata.
    _header_references_external_wcs_data(header, alt_str) || return NoAuxiliaryWCSData()

    throw(
        ArgumentError(
            "no auxiliary WCS data resolver is defined for fobj type $(typeof(fobj))"
        )
    )
end
