# Time axis CTYPE prefixes (FITS Paper III / WCS Paper VI).
# These identify axes whose world coordinate is a time value.
const _TIME_SYSTEMS = Set(["TIME", "UTC", "TAI", "TDB", "TT", "TCG", "TCB", "LOCAL"])


# ── Time axis unit conversion ──────────────────────────────────────────────────

const _TIME_UNIT_TO_SECOND = Dict{String, Float64}(
    "S"  => 1.0,      "SEC" => 1.0,    "SECOND" => 1.0,
    "MS" => 1e-3,     "US"  => 1e-6,   "NS"  => 1e-9,
    "MIN"=> 60.0,     "MINUTE"=>60.0,
    "H"  => 3600.0,   "HR"  => 3600.0, "HOUR"=> 3600.0,
    "D"  => 86400.0,  "DAY" => 86400.0,
    "YR" => 31557600.0, "YEAR" => 31557600.0, "A" => 31557600.0,
    ""   => 1.0,
)

"""Convert a FITS time unit string to seconds (the canonical time unit)."""
function _unit_to_second(unit_str::AbstractString)
    u = uppercase(strip(unit_str))
    return get(_TIME_UNIT_TO_SECOND, u, 1.0)
end

# ── Time axis types ────────────────────────────────────────────────────────────

"""Abstract supertype for time WCS payloads."""
abstract type AbstractTimeWCSData end

"""No-op time payload for WCS transforms with no time axes."""
struct NoTimeWCSData <: AbstractTimeWCSData end

"""
    TimeSpec

Parsed time-axis specification.  Time axes are purely linear (pixel to seconds)
with no algorithm variants.  The keywords are stored as pass-through metadata
for downstream frame-correction code.
"""
struct TimeSpec
    axis::Int           # WCS axis number (1-based)
    mjdref::Float64     # MJDREF -- reference Modified Julian Date (NaN if absent)
    timesys::String     # TIMESYS -- time system ("" if absent)
    trefpos::String     # TREFPOS -- reference position ("" if absent)
    trefdir::String     # TREFDIR -- reference direction ("" if absent)
    timeunit::String    # TIMEUNIT -- unit for time keywords ("" if absent)
end

"""Resolved collection of all time-axis specifications for one WCS."""
struct TimeWCSData{T <: Tuple} <: AbstractTimeWCSData
    specs::T
end
