"""
Unit conversion factors to degrees.  Returns the multiplier so that
`value_in_deg = value * unit_to_deg(unit)`.

Supports the unit strings used in FITS `CUNITi` keywords.
"""
function unit_to_deg(unit::AbstractString)::Float64
    u = lowercase(strip(unit))
    if u == "deg" || u == ""
        return 1.0
    elseif u == "arcmin"
        return 1.0 / 60.0
    elseif u == "arcsec"
        return 1.0 / 3600.0
    elseif u == "mas"          # milli-arcsecond
        return 1.0 / 3_600_000.0
    elseif u == "rad"
        return 180.0 / π
    elseif u == "hr" || u == "hour"
        return 15.0
    elseif u == "min"          # minutes of time
        return 15.0 / 60.0
    elseif u == "s" || u == "sec"  # seconds of time
        return 15.0 / 3600.0
    else
        # Return NaN to signal unknown unit; callers should decide what to do.
        return NaN
    end
end
