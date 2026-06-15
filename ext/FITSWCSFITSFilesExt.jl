module FITSWCSFITSFilesExt

import FITSFiles
import FITSWCS: from_header

function _fitsfiles_cards_dict(cards::FITSFiles.Cards)
    # Copy card keyword/value pairs into the dictionary shape used by core parsing.
    dict = Dict{String, Any}()
    for card in cards
        dict[String(card.key)] = card.value
    end
    return dict
end

function from_header(cards::FITSFiles.Cards; alt::Char=' ')
    # Delegate all WCS validation and interpretation to the core parser.
    return from_header(_fitsfiles_cards_dict(cards); alt=alt)
end

function from_header(hdu::FITSFiles.HDU; alt::Char=' ')
    # FITSFiles stores parsed header cards directly on each HDU.
    return from_header(hdu.cards; alt=alt)
end

end # module FITSWCSFITSFilesExt
