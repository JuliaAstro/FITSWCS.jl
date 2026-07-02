module FITSWCSFITSFilesExt

import FITSFiles
import FITSWCS:
    WCS,
    NoAuxiliaryWCSData,
    _auxiliary_wcs_data,
    _header_references_external_wcs_data,
    _lookup_table_from_image,
    _paper_iv_auxiliary_data

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

    # Load referenced Paper IV image HDUs and copy them into backend-neutral tables.
    return _paper_iv_auxiliary_data(
        header, spec -> begin
            hdu = _fitsfiles_lookup_hdu(fobj, spec.extname, spec.extver)
            _lookup_table_from_image(hdu.data, _fitsfiles_cards_dict(hdu.cards), spec.transpose)
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
