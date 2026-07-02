module FITSWCSFITSFilesExt

import FITSFiles
import FITSWCS: WCS

function _fitsfiles_cards_dict(cards::FITSFiles.Cards)
    # Copy card keyword/value pairs into the dictionary shape used by core parsing.
    dict = Dict{String, Any}()
    for card in cards
        dict[String(card.key)] = card.value
    end
    return dict
end

function WCS(cards::FITSFiles.Cards; alt::Char=' ')
    # Delegate all WCS validation and interpretation to the core parser.
    return WCS(_fitsfiles_cards_dict(cards); alt=alt)
end

function WCS(hdu::FITSFiles.HDU; alt::Char=' ')
    # FITSFiles stores parsed header cards directly on each HDU.
    return WCS(hdu.cards; alt=alt)
end

end # module FITSWCSFITSFilesExt
