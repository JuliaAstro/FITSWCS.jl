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

AuxiliaryWCSData(; det2im = (nothing, nothing), cpdis = (nothing, nothing), tabular = NoTabularWCSData()) =
    AuxiliaryWCSData(det2im, cpdis, tabular)

struct PaperIVLookupSpec
    extname::String
    extver::Int
    transpose::Bool
end

has_auxiliary(::NoAuxiliaryWCSData) = false
has_auxiliary(::AbstractAuxiliaryWCSData) = true

function is_lookup_distortion_keyword(key::AbstractString, alt_str::AbstractString = "")
    ukey = uppercase(String(key))
    suffix = uppercase(String(alt_str))

    # Astropy exposes Paper IV lookup tables as CPDIS and detector-to-image D2IM.
    occursin(Regex("^CPDIS[1-9][0-9]*$(suffix)\$"), ukey) && return true
    occursin(Regex("^D2IMDIS[1-9][0-9]*$(suffix)\$"), ukey) && return true
    occursin(Regex("^D2IMERR[1-9][0-9]*$(suffix)\$"), ukey) && return true
    ukey == "AXISCORR$(suffix)" && return true

    # wcslib Paper IV distortion parameters use DPja/DQia keyword families.
    occursin(Regex("^D[QP][1-9][0-9]*$(suffix)\\."), ukey) && return true

    return false
end

function _paper_iv_lookup_spec(
        header::AbstractDict,
        axis::Int,
        dist_prefix::AbstractString,
        param_prefix::AbstractString,
        err_prefix::AbstractString,
        extname::AbstractString,
        alt_str::AbstractString,
        minerr::Real,
    )
    dist_key = "$(dist_prefix)$(axis)$(alt_str)"
    haskey(header, dist_key) || return nothing
    uppercase(String(header[dist_key])) == "LOOKUP" ||
        throw(ArgumentError("unsupported Paper IV distortion type $(header[dist_key]) in $dist_key"))

    # Skip tables whose declared error is below the requested threshold.
    err_key = "$(err_prefix)$(axis)$(alt_str)"
    Float64(get(header, err_key, 0.0)) < Float64(minerr) && return nothing

    param = "$(param_prefix)$(axis)$(alt_str)"
    axis_key = "$(param).AXIS.$(axis)"
    haskey(header, axis_key) ||
        throw(ArgumentError("Paper IV lookup distortion $dist_key requires $axis_key"))
    table_axis = Int(header[axis_key])
    table_axis in (1, 2) ||
        throw(ArgumentError("$axis_key must be 1 or 2, got $table_axis"))

    # Astropy transposes lookup image data when the table axis differs.
    extver = Int(get(header, "$(param).EXTVER", 1))
    return PaperIVLookupSpec(String(extname), extver, table_axis != axis)
end

function _paper_iv_lookup_specs(header::AbstractDict; alt::Char = ' ', minerr::Real = 0.0)
    alt_str = alt == ' ' ? "" : string(alt)

    # Collect the two Paper IV image lookup families used by Astropy.
    det2im = ntuple(i -> _paper_iv_lookup_spec(header, i, "D2IMDIS", "D2IM", "D2IMERR", "D2IMARR", alt_str, minerr), 2)
    cpdis = ntuple(i -> _paper_iv_lookup_spec(header, i, "CPDIS", "DP", "CPERR", "WCSDVARR", alt_str, minerr), 2)
    return det2im, cpdis
end

function _lookup_table_from_image(data::AbstractMatrix, header::AbstractDict, transpose::Bool)
    table_data = transpose ? permutedims(data) : copy(data)

    # Preserve the lookup image's own linear coordinate metadata.
    return LookupTable2D(
        table_data;
        crpix = (get(header, "CRPIX1", 0.0), get(header, "CRPIX2", 0.0)),
        crval = (get(header, "CRVAL1", 0.0), get(header, "CRVAL2", 0.0)),
        cdelt = (get(header, "CDELT1", 1.0), get(header, "CDELT2", 1.0)),
    )
end

function _paper_iv_auxiliary_data(header::AbstractDict, loader; alt::Char = ' ', minerr::Real = 0.0)
    det2im_specs, cpdis_specs = _paper_iv_lookup_specs(header; alt = alt, minerr = minerr)

    # Resolve all referenced image arrays through the backend-provided loader.
    det2im = ntuple(i -> isnothing(det2im_specs[i]) ? nothing : loader(det2im_specs[i]), 2)
    cpdis = ntuple(i -> isnothing(cpdis_specs[i]) ? nothing : loader(cpdis_specs[i]), 2)

    if all(isnothing, det2im) && all(isnothing, cpdis)
        return NoAuxiliaryWCSData()
    end
    return AuxiliaryWCSData(det2im = det2im, cpdis = cpdis)
end

function _external_auxiliary_data(
        header::AbstractDict,
        paper_iv_loader,
        tabular_loader;
        alt::Char = ' ',
        minerr::Real = 0.0,
    )
    paper_iv = _paper_iv_auxiliary_data(header, paper_iv_loader; alt = alt, minerr = minerr)
    tabular = _tabular_auxiliary_data(header, tabular_loader; alt = alt)

    # Preserve the cheap no-auxiliary payload when neither external family is present.
    paper_iv isa NoAuxiliaryWCSData && tabular isa NoTabularWCSData &&
        return NoAuxiliaryWCSData()

    det2im = paper_iv isa AuxiliaryWCSData ? paper_iv.det2im : (nothing, nothing)
    cpdis = paper_iv isa AuxiliaryWCSData ? paper_iv.cpdis : (nothing, nothing)
    return AuxiliaryWCSData(det2im = det2im, cpdis = cpdis, tabular = tabular)
end

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
