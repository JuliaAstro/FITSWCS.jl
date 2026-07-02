module FITSWCSFITSFilesExt

import FITSFiles
import FITSWCS:
    WCS,
    NoAuxiliaryWCSData,
    _auxiliary_wcs_data,
    _header_references_external_wcs_data

function _fitsfiles_cards_dict(cards::FITSFiles.Cards)
    # Copy card keyword/value pairs into the dictionary shape used by core parsing.
    dict = Dict{String, Any}()
    for card in cards
        dict[String(card.key)] = card.value
    end
    return dict
end

function _fitsfiles_auxiliary_wcs_data(header::AbstractDict, fobj; alt::Char = ' ', minerr::Real = 0.0)
    alt_str = alt == ' ' ? "" : string(alt)

    # Keep the common header-only path cheap even when a FITSFiles object is supplied.
    _header_references_external_wcs_data(header, alt_str) || return NoAuxiliaryWCSData()

    # Methods for specific fobj types below; if none match, throw an error.
    throw(
        ArgumentError(
            "FITSFiles auxiliary WCS data resolution is not implemented yet for fobj type $(typeof(fobj))"
        )
    )
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

end # module FITSWCSFITSFilesExt
