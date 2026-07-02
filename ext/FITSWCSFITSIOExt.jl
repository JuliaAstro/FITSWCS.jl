module FITSWCSFITSIOExt

import FITSIO
import FITSWCS: WCS

function _fitsio_header_dict(header::FITSIO.FITSHeader)
    # Copy keyword/value pairs into the dictionary shape used by core parsing.
    dict = Dict{String, Any}()
    for (key, value) in zip(keys(header), values(header))
        dict[String(key)] = value
    end
    return dict
end

function WCS(header::FITSIO.FITSHeader; alt::Char=' ')
    # Delegate all WCS validation and interpretation to the core parser.
    return WCS(_fitsio_header_dict(header); alt=alt)
end

function WCS(hdu::FITSIO.HDU; alt::Char=' ')
    # Read the HDU header through FITSIO before using the header adapter.
    return WCS(FITSIO.read_header(hdu); alt=alt)
end

end # module FITSWCSFITSIOExt
