module FITSWCSFITSIOExt

import FITSIO
import FITSWCS:
    WCS,
    NoAuxiliaryWCSData,
    _auxiliary_wcs_data,
    _header_references_external_wcs_data

function _fitsio_header_dict(header::FITSIO.FITSHeader)
    # Copy keyword/value pairs into the dictionary shape used by core parsing.
    dict = Dict{String, Any}()
    for (key, value) in zip(keys(header), values(header))
        dict[String(key)] = value
    end
    return dict
end

function _fitsio_auxiliary_wcs_data(header::AbstractDict, fobj; alt::Char=' ', minerr::Real=0.0)
    alt_str = alt == ' ' ? "" : string(alt)

    # Keep the common header-only path cheap even when a FITSIO object is supplied.
    _header_references_external_wcs_data(header, alt_str) || return NoAuxiliaryWCSData()

    # Methods for specific fobj types below; if none match, throw an error.
    throw(ArgumentError(
        "FITSIO auxiliary WCS data resolution is not implemented yet for fobj type $(typeof(fobj))"
    ))
end

function _auxiliary_wcs_data(header::AbstractDict, fobj::FITSIO.FITS; alt::Char=' ', minerr::Real=0.0)
    # Route FITSIO file containers through the extension-owned resolver.
    return _fitsio_auxiliary_wcs_data(header, fobj; alt=alt, minerr=minerr)
end

function _auxiliary_wcs_data(header::AbstractDict, fobj::FITSIO.HDU; alt::Char=' ', minerr::Real=0.0)
    # Allow callers to pass a FITSIO HDU when no file-level container is available.
    return _fitsio_auxiliary_wcs_data(header, fobj; alt=alt, minerr=minerr)
end

function WCS(header::FITSIO.FITSHeader; fobj=nothing, alt::Char=' ', minerr::Real=0.0)
    # Delegate all WCS validation and interpretation to the core parser.
    return WCS(_fitsio_header_dict(header); fobj=fobj, alt=alt, minerr=minerr)
end

function WCS(hdu::FITSIO.HDU; fobj=nothing, alt::Char=' ', minerr::Real=0.0)
    # Read the HDU header through FITSIO before using the header adapter.
    return WCS(FITSIO.read_header(hdu); fobj=fobj, alt=alt, minerr=minerr)
end

end # module FITSWCSFITSIOExt
