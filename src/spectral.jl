"""
Spectral coordinate transforms (FITS WCS Paper III).

Implements the S/P/X-type architecture: every spectral CTYPE encodes
three spectral variables — the S-type (what the user sees), the P-type
(the parent basic type of S), and the X-type (the basic type in which the
axis is linear).  The transform chain is:

    pixel → [CD] → intermediate (X-type) → [X→P] → P-type → [P→S] → world (S-type)

For linear spectral axes (blank algorithm code) X ≡ P ≡ S and no
conversion is needed.  For logarithmic axes the intermediate coordinate
is ``ln(X/X_ref)`` and the world coordinate is ``X_ref × exp(intermediate)``.
"""

# ── Type hierarchy ───────────────────────────────────────────────────────────

"""Abstract supertype for spectral WCS payloads."""
abstract type AbstractSpectralWCSData end

"""No-op spectral payload for WCS transforms with no spectral axes."""
struct NoSpectralWCSData <: AbstractSpectralWCSData end

"""
    SpectralSpec{X, P, S, A}

Parsed spectral-axis specification carrying the S/P/X types and algorithm
code as `Val` type parameters for compile-time dispatch in the hot path.

Type parameters:
- `X`: `Val{:F}`, `Val{:W}`, `Val{:A}`, or `Val{:V}` — X-type basic type
- `P`: same set — P-type parent basic type
- `S`: `Val{:FREQ}`, `Val{:WAVE}`, `Val{:VELO}`, etc. — S-type coordinate type
- `A`: `Val{:LINEAR}`, `Val{:LOG}`, `Val{:F2W}`, etc. — algorithm code
"""
struct SpectralSpec{X, P, S, A}
    axis::Int           # WCS axis number (1-based)
    restfrq::Float64    # RESTFRQa — rest frequency in Hz (NaN if absent)
    restwav::Float64    # RESTWAVa — rest wavelength in m (NaN if absent)
    crval_si::Float64   # CRVAL in S-type SI units (world coordinate at ref)
    x_r::Float64        # X-type reference value in SI (for offset recovery)
    cdelt_scale::Float64 # dX/dS evaluated at reference point (for CDELT scaling)
    # Reference frame metadata (stored for downstream frame-correction code;
    # not used in the core pixel↔world transform):
    specsys::String     # SPECSYSa  — frame of spectral coords (e.g. "LSRK")
    ssysobs::String     # SSYSOBSa — frame held constant during observation
    velosys::Float64    # VELOSYSa — observer→rest-frame radial velocity (m/s)
    zsource::Float64    # ZSOURCEa — source redshift (unitless)
    ssyssrc::String     # SSYSSRCa — reference frame for SOURCE cases
end

"""Resolved collection of all spectral-axis specifications for one WCS."""
struct SpectralWCSData{T <: Tuple} <: AbstractSpectralWCSData
    specs::T
end

"""Observation-level metadata needed for spectral reference frame corrections."""
struct ObservationSpec
    mjd_avg::Float64    # MJD-AVG    — representative MJD of observation
    date_avg::String    # DATE-AVG   — representative ISO date/time
    date_obs::String    # DATE-OBS   — observation start (fallback)
    obsgeo_x::Float64   # OBSGEO-X   — geocentric X (m, NaN if absent)
    obsgeo_y::Float64   # OBSGEO-Y   — geocentric Y (m, NaN if absent)
    obsgeo_z::Float64   # OBSGEO-Z   — geocentric Z (m, NaN if absent)
    velangl::Float64    # VELANGL    — true velocity angle from line of sight
end

# ── S-type ↔ P-type tables ──────────────────────────────────────────────────

# Map S-type symbol to its parent basic-type symbol.
const _S_TO_P = Dict{Symbol, Symbol}(
    :FREQ => :F,  :AFRQ => :F,  :ENER => :F,  :WAVN => :F,  :VRAD => :F,
    :WAVE => :W,  :VOPT => :W,  :ZOPT => :W,
    :AWAV => :A,
    :VELO => :V,  :BETA => :V,
)

# Map S-type string (4-char, uppercase) to its symbol.
const _STYPE_STR = Dict{String, Symbol}(
    "FREQ" => :FREQ, "AFRQ" => :AFRQ, "ENER" => :ENER, "WAVN" => :WAVN,
    "VRAD" => :VRAD, "WAVE" => :WAVE, "VOPT" => :VOPT, "ZOPT" => :ZOPT,
    "AWAV" => :AWAV, "VELO" => :VELO, "BETA" => :BETA,
)

# Map algorithm-code string (3-char) to its (X-type, P-type, algorithm symbol).
# A blank code (or "---") means linear: X ≡ P ≡ S-type parent.
const _ALGORITHM_MAP = Dict{String, Tuple{Symbol, Symbol, Symbol}}(
    ""    => (:F, :F, :LINEAR),   # placeholder — resolved after S-type is known
    "LOG" => (:F, :F, :LOG),     # placeholder — X/P resolved from S-type
    "F2W" => (:F, :W, :F2W),  "W2F" => (:W, :F, :W2F),
    "F2V" => (:F, :V, :F2V),  "V2F" => (:V, :F, :V2F),
    "W2V" => (:W, :V, :W2V),  "V2W" => (:V, :W, :V2W),
    "F2A" => (:F, :A, :F2A),  "A2F" => (:A, :F, :A2F),
    "W2A" => (:W, :A, :W2A),  "A2W" => (:A, :W, :A2W),
    "V2A" => (:V, :A, :V2A),  "A2V" => (:A, :V, :A2V),
)

# ── CDELT scaling: dX/dS derivatives ─────────────────────────────────────────

"""
    _dpds(s_type, s_si, spec) -> Float64

Derivative of P-type (SI) wrt S-type (display) evaluated at reference point.
From Paper III Table 4 (property equations).  For basic types (FREQ→F, WAVE→W,
VELO→V, AWAV→A) this is 1.0; for derived types the linear P↔S relation gives a
constant derivative.
"""
_dpds(::Val{:FREQ}, s_si, spec) = 1.0
_dpds(::Val{:WAVE}, s_si, spec) = 1.0
_dpds(::Val{:VELO}, s_si, spec) = 1.0
_dpds(::Val{:AWAV}, s_si, spec) = 1.0
_dpds(::Val{:AFRQ}, s_si, spec) = 2π
_dpds(::Val{:ENER}, s_si, spec) = inv(_H_PLANCK)
_dpds(::Val{:WAVN}, s_si, spec) = 1.0
function _dpds(::Val{:VRAD}, s_si, spec)
    return -spec.restfrq / _C_LIGHT
end
function _dpds(::Val{:VOPT}, s_si, spec)
    return spec.restwav / _C_LIGHT
end
function _dpds(::Val{:ZOPT}, s_si, spec)
    return spec.restwav
end
function _dpds(::Val{:BETA}, s_si, spec)
    return inv(_C_LIGHT)
end

"""
    _dxdp(x_type, p_type, algo, p_si, spec) -> Float64

Derivative of X-type (SI) wrt P-type (SI) evaluated at reference point.
From Paper III Table 3 (basic spectral transformations) plus the IUGG/Ciddor
air-wavelength derivative for air-wavelength types.
"""
function _dxdp(::X, ::P, ::A, p_si, spec) where {X, P, A}
    if X === P
        return 1.0
    elseif A === Val{:F2W}  # dν/dλ = -c/λ²
        return -_C_LIGHT / (p_si * p_si)
    elseif A === Val{:W2F}  # dλ/dν = -c/ν²
        return -_C_LIGHT / (p_si * p_si)
    elseif A === Val{:F2V}  # dν/dv = -c·ν₀ / ((c+v)·√(c²-v²))
        v = p_si
        return -_C_LIGHT * spec.restfrq / ((_C_LIGHT + v) * sqrt(_C_LIGHT*_C_LIGHT - v*v))
    elseif A === Val{:V2F}  # dv/dν = -4c·ν·ν₀²/(ν²+ν₀²)²
        ν = p_si
        ν₀ = spec.restfrq
        return -4*_C_LIGHT*ν*ν₀*ν₀ / ((ν*ν + ν₀*ν₀)^2)
    elseif A === Val{:W2V}  # dλ/dv = c·λ₀ / ((c-v)·√(c²-v²))
        v = p_si
        return _C_LIGHT * spec.restwav / ((_C_LIGHT - v) * sqrt(_C_LIGHT*_C_LIGHT - v*v))
    elseif A === Val{:V2W}  # dv/dλ = 4c·λ·λ₀²/(λ²+λ₀²)²
        λ = p_si
        λ₀ = spec.restwav
        return 4*_C_LIGHT*λ*λ₀*λ₀ / ((λ*λ + λ₀*λ₀)^2)
    elseif A === Val{:F2A}  # dν/dλ_a = (dν/dλ)·(dλ/dλ_a)
        return _dxdp(Val{:F}, Val{:W}, Val{:F2W}, p_si, spec) * _dwave_dawav(p_si)
    elseif A === Val{:A2F}  # dλ_a/dν = (dλ_a/dλ)·(dλ/dν)
        ν = p_si
        λ = _C_LIGHT / ν
        return _dawav_dwave(λ) * (-_C_LIGHT / (ν*ν))
    elseif A === Val{:W2A}  # dλ/dλ_a
        return _dwave_dawav(p_si)
    elseif A === Val{:A2W}  # dλ_a/dλ
        return _dawav_dwave(p_si)
    elseif A === Val{:V2A}  # dv/dλ_a = (dv/dλ)·(dλ/dλ_a)
        return _dxdp(Val{:V}, Val{:W}, Val{:V2W}, p_si, spec) * _dwave_dawav(p_si)
    elseif A === Val{:A2V}  # dλ_a/dv = (dλ_a/dλ)·(dλ/dv)
        v = p_si
        λ = _C_LIGHT / _velo_to_freq(v, spec.restfrq)
        return _dawav_dwave(λ) * _dxdp(Val{:W}, Val{:V}, Val{:W2V}, v, spec)
    end
    return 1.0
end

"""Derivative of vacuum wavelength wrt air wavelength (IUGG/Ciddor)."""
function _dwave_dawav(λ_vac)
    λ_μm = λ_vac * 1e6
    λ² = λ_μm * λ_μm
    n = _refractive_index_air(λ_μm)
    dn = -2e-6 / λ_μm * (1.62887/λ² + 2*0.01360/(λ²*λ²))
    return n - λ_vac * 1e6 * dn  # dλ_vac/dλ_air = n + λ_air·dn/dλ_air
end

"""Derivative of air wavelength wrt vacuum wavelength."""
function _dawav_dwave(λ_vac)
    return inv(_dwave_dawav(λ_vac))
end

"""
    _dxds(spec) -> Float64

Compute dX/dS at the reference point: the factor to multiply CDELT by so the
CD step produces intermediate in X-type SI units.
"""
function _dxds(spec::SpectralSpec{X, P, S, A}) where {X, P, S, A}
    if A === Val{:LOG} || A === Val{:LINEAR}
        return 1.0    # X = S for LOG/LINEAR, so dX/dS = 1
    end
    # Cross-type: dX/dS = dX/dP · dP/dS
    p_si = _s_to_p(spec.crval_si, S(), spec)
    return _dxdp(X(), P(), A(), p_si, spec) * _dpds(S(), spec.crval_si, spec)
end

# ── Physical constants (SI) ──────────────────────────────────────────────────

const _C_LIGHT = 2.99792458e8   # speed of light, m/s
const _H_PLANCK = 6.62607015e-34  # Planck constant, J·s

# ── Unit conversion table ────────────────────────────────────────────────────

const _SPECTRAL_UNIT_TO_SI = Dict{String, Float64}(
    # Frequency
    "HZ"  => 1.0,    "KHZ" => 1e3,   "MHZ" => 1e6,
    "GHZ" => 1e9,    "THZ" => 1e12,
    # Wavelength
    "M"   => 1.0,    "CM"  => 1e-2,  "MM"  => 1e-3,
    "UM"  => 1e-6,   "NM"  => 1e-9,  "A"   => 1e-10,
    "ANGSTROM" => 1e-10,
    # Velocity
    "M/S" => 1.0,    "KM/S"=> 1e3,
    # Dimensionless
    ""    => 1.0,
)

"""Convert a FITS spectral unit string to an SI scale factor."""
function _unit_to_si(unit_str::AbstractString)
    u = uppercase(strip(unit_str))
    return get(_SPECTRAL_UNIT_TO_SI, u, 1.0)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Hot-path transform functions below.
# These preserve input float type.  Parse-time code above (derivatives,
# CDELT scaling, unit conversion) uses Float64 as expected by the structs
# that will contain their output.
# ═══════════════════════════════════════════════════════════════════════════════

# ── P-type ↔ X-type ↔ S-type dispatch tables ────────────────────────────────

# Maps (from_type, to_type) to conversion function for P↔X step.
# These are the core non-linear conversions.

# X → P (forward): intermediate → P-type
function _x_to_p(x, ::Val{:F}, ::Val{:F}, ::Val{:LINEAR}, spec::SpectralSpec)
    return x  # linear
end
function _x_to_p(x, ::Val{:W}, ::Val{:W}, ::Val{:LINEAR}, spec::SpectralSpec)
    return x
end
function _x_to_p(x, ::Val{:A}, ::Val{:A}, ::Val{:LINEAR}, spec::SpectralSpec)
    return x
end
function _x_to_p(x, ::Val{:V}, ::Val{:V}, ::Val{:LINEAR}, spec::SpectralSpec)
    return x
end
# LOG: intermediate w is in S-type SI units, S = S_r * exp(w / S_r)
function _x_to_p(x, ::Val{:F}, ::Val{:F}, ::Val{:LOG}, spec::SpectralSpec)
    crval_si = oftype(x, spec.crval_si)
    return crval_si * exp(x / crval_si)
end
function _x_to_p(x, ::Val{:W}, ::Val{:W}, ::Val{:LOG}, spec::SpectralSpec)
    crval_si = oftype(x, spec.crval_si)
    return crval_si * exp(x / crval_si)
end
# F2W / W2F: frequency ↔ vacuum wavelength
function _x_to_p(x, ::Val{:F}, ::Val{:W}, ::Val{:F2W}, spec::SpectralSpec)
    return oftype(x, _C_LIGHT) / x
end
function _x_to_p(x, ::Val{:W}, ::Val{:F}, ::Val{:W2F}, spec::SpectralSpec)
    return oftype(x, _C_LIGHT) / x
end
# F2V / V2F: frequency ↔ relativistic velocity (already type-preserving)
function _x_to_p(x, ::Val{:F}, ::Val{:V}, ::Val{:F2V}, spec::SpectralSpec)
    return _freq_to_velo(x, oftype(x, spec.restfrq))
end
function _x_to_p(x, ::Val{:V}, ::Val{:F}, ::Val{:V2F}, spec::SpectralSpec)
    return _velo_to_freq(x, oftype(x, spec.restfrq))
end
# W2V / V2W: vacuum wavelength ↔ relativistic velocity
function _x_to_p(x, ::Val{:W}, ::Val{:V}, ::Val{:W2V}, spec::SpectralSpec)
    return _freq_to_velo(oftype(x, _C_LIGHT) / x, oftype(x, spec.restfrq))
end
function _x_to_p(x, ::Val{:V}, ::Val{:W}, ::Val{:V2W}, spec::SpectralSpec)
    return oftype(x, _C_LIGHT) / _velo_to_freq(x, oftype(x, spec.restfrq))
end
# Air-wavelength cross-conversions: route through wave↔awav.
function _x_to_p(x, ::Val{:F}, ::Val{:A}, ::Val{:F2A}, spec::SpectralSpec)
    return _wave_to_awav(oftype(x, _C_LIGHT) / x)
end
function _x_to_p(x, ::Val{:W}, ::Val{:A}, ::Val{:W2A}, spec::SpectralSpec)
    return _wave_to_awav(x)
end
function _x_to_p(x, ::Val{:V}, ::Val{:A}, ::Val{:V2A}, spec::SpectralSpec)
    return _wave_to_awav(oftype(x, _C_LIGHT) / _velo_to_freq(x, oftype(x, spec.restfrq)))
end
function _x_to_p(x, ::Val{:A}, ::Val{:F}, ::Val{:A2F}, spec::SpectralSpec)
    return oftype(x, _C_LIGHT) / _awav_to_wave(x)
end
function _x_to_p(x, ::Val{:A}, ::Val{:W}, ::Val{:A2W}, spec::SpectralSpec)
    return _awav_to_wave(x)
end
function _x_to_p(x, ::Val{:A}, ::Val{:V}, ::Val{:A2V}, spec::SpectralSpec)
    return _freq_to_velo(oftype(x, _C_LIGHT) / _awav_to_wave(x), oftype(x, spec.restfrq))
end

# P → S (forward): P-type → world (S-type).  Linear conversion for derived S-types.
_p_to_s(p, ::Val{:FREQ}, spec) = p
_p_to_s(p, ::Val{:AFRQ}, spec) = oftype(p, 2π) * p
_p_to_s(p, ::Val{:ENER}, spec) = oftype(p, _H_PLANCK) * p
_p_to_s(p, ::Val{:WAVN}, spec) = p
function _p_to_s(p, ::Val{:VRAD}, spec)
    T = typeof(p)
    return T(_C_LIGHT) * (T(spec.restfrq) - p) / T(spec.restfrq)
end
_p_to_s(p, ::Val{:WAVE}, spec) = p
function _p_to_s(p, ::Val{:VOPT}, spec)
    T = typeof(p)
    return T(_C_LIGHT) * (p - T(spec.restwav)) / T(spec.restwav)
end
function _p_to_s(p, ::Val{:ZOPT}, spec)
    restwav = oftype(p, spec.restwav)
    return (p - restwav) / restwav
end
_p_to_s(p, ::Val{:AWAV}, spec) = p
_p_to_s(p, ::Val{:VELO}, spec) = p
_p_to_s(p, ::Val{:BETA}, spec) = p / oftype(p, _C_LIGHT)

# S → P (inverse): world (S-type) → P-type
_s_to_p(s, ::Val{:FREQ}, spec) = s
_s_to_p(s, ::Val{:AFRQ}, spec) = s / oftype(s, 2π)
_s_to_p(s, ::Val{:ENER}, spec) = s / oftype(s, _H_PLANCK)
_s_to_p(s, ::Val{:WAVN}, spec) = s
function _s_to_p(s, ::Val{:VRAD}, spec)
    return oftype(s, spec.restfrq) * (1 - s / oftype(s, _C_LIGHT))
end
_s_to_p(s, ::Val{:WAVE}, spec) = s
function _s_to_p(s, ::Val{:VOPT}, spec)
    return oftype(s, spec.restwav) * (1 + s / oftype(s, _C_LIGHT))
end
function _s_to_p(s, ::Val{:ZOPT}, spec)
    return oftype(s, spec.restwav) * (1 + s)
end
_s_to_p(s, ::Val{:AWAV}, spec) = s
_s_to_p(s, ::Val{:VELO}, spec) = s
_s_to_p(s, ::Val{:BETA}, spec) = s * oftype(s, _C_LIGHT)

# P → X (inverse): P-type → intermediate.  Reverse of X→P.

# LINEAR inverse (identity for all basic types)
_p_to_x(p, ::Val{:F}, ::Val{:F}, ::Val{:LINEAR}, spec) = p
_p_to_x(p, ::Val{:W}, ::Val{:W}, ::Val{:LINEAR}, spec) = p
_p_to_x(p, ::Val{:A}, ::Val{:A}, ::Val{:LINEAR}, spec) = p
_p_to_x(p, ::Val{:V}, ::Val{:V}, ::Val{:LINEAR}, spec) = p

# LOG inverse: w = S_r * ln(S / S_r)  (w in S-type SI units)
function _p_to_x(p, ::Val{:F}, ::Val{:F}, ::Val{:LOG}, spec::SpectralSpec)
    crval_si = oftype(p, spec.crval_si)
    return crval_si * log(p / crval_si)
end
function _p_to_x(p, ::Val{:W}, ::Val{:W}, ::Val{:LOG}, spec::SpectralSpec)
    crval_si = oftype(p, spec.crval_si)
    return crval_si * log(p / crval_si)
end
# F2W / W2F: c/x is its own inverse.
function _p_to_x(p, ::Val{:W}, ::Val{:F}, ::Val{:F2W}, spec::SpectralSpec)
    return oftype(p, _C_LIGHT) / p
end
function _p_to_x(p, ::Val{:F}, ::Val{:W}, ::Val{:W2F}, spec::SpectralSpec)
    return oftype(p, _C_LIGHT) / p
end
# F2V / V2F (already type-preserving via _velo_to_freq/_freq_to_velo)
function _p_to_x(p, ::Val{:V}, ::Val{:F}, ::Val{:F2V}, spec::SpectralSpec)
    return _velo_to_freq(p, spec.restfrq)
end
function _p_to_x(p, ::Val{:F}, ::Val{:V}, ::Val{:V2F}, spec::SpectralSpec)
    return _freq_to_velo(p, oftype(p, spec.restfrq))
end
# W2V / V2W
function _p_to_x(p, ::Val{:V}, ::Val{:W}, ::Val{:W2V}, spec::SpectralSpec)
    return oftype(p, _C_LIGHT) / _velo_to_freq(p, oftype(p, spec.restfrq))
end
function _p_to_x(p, ::Val{:W}, ::Val{:V}, ::Val{:V2W}, spec::SpectralSpec)
    return _freq_to_velo(oftype(p, _C_LIGHT) / p, oftype(p, spec.restfrq))
end

# Air-wavelength inverse cross-conversions.
function _p_to_x(p, ::Val{:A}, ::Val{:F}, ::Val{:F2A}, spec::SpectralSpec)
    return oftype(p, _C_LIGHT) / _awav_to_wave(p)
end
function _p_to_x(p, ::Val{:A}, ::Val{:W}, ::Val{:W2A}, spec::SpectralSpec)
    return _awav_to_wave(p)
end
function _p_to_x(p, ::Val{:A}, ::Val{:V}, ::Val{:V2A}, spec::SpectralSpec)
    return _velo_to_freq(oftype(p, _C_LIGHT) / _awav_to_wave(p), oftype(p, spec.restfrq))
end
function _p_to_x(p, ::Val{:F}, ::Val{:A}, ::Val{:A2F}, spec::SpectralSpec)
    return _wave_to_awav(oftype(p, _C_LIGHT) / p)
end
function _p_to_x(p, ::Val{:W}, ::Val{:A}, ::Val{:A2W}, spec::SpectralSpec)
    return _wave_to_awav(p)
end
function _p_to_x(p, ::Val{:V}, ::Val{:A}, ::Val{:A2V}, spec::SpectralSpec)
    return _wave_to_awav(oftype(p, _C_LIGHT) / _freq_to_velo(p, oftype(p, spec.restfrq)))
end

# ── Core spectral conversion functions ───────────────────────────────────────

"""Frequency to relativistic velocity (apparent radial velocity)."""
function _freq_to_velo(ν, ν₀)
    T = float(promote_type(typeof(ν), typeof(ν₀)))
    r = ν₀ / ν
    return T(_C_LIGHT) * (r*r - 1) / (r*r + 1)
end

"""Relativistic velocity to frequency."""
function _velo_to_freq(v, ν₀)
    T = float(promote_type(typeof(v), typeof(ν₀)))
    return ν₀ * sqrt((T(_C_LIGHT) - v) / (T(_C_LIGHT) + v))
end

# ── Air-wavelength conversions (IUGG/Ciddor 1999, Paper III §4) ─────────────

"""
    _refractive_index_air(λ_μm)

Refractive index of dry air at standard temperature and pressure
evaluated at the given wavelength in μm.  Uses the IUGG 1999 /
Ciddor (1996) relation as published in Paper III Eq. 4.

!!! note
    WCSLIB uses the older Cox/Edlén (IAU 1957) relation.  The two
    differ by a nearly constant ratio IUGG/Cox ≈ 1.000015 across
    the optical range.  This implementation follows the published
    Paper III standard.
"""
function _refractive_index_air(λ_μm::T) where {T}
    FT = float(T)
    λ² = λ_μm * λ_μm
    return 1 + (FT(287.6155) + FT(1.62887) / λ² + FT(0.01360) / (λ² * λ²)) / 1_000_000
end

"""
    _wave_to_awav(λ_vac)

Vacuum wavelength → air wavelength (both in meters).

Paper III §4 notes that because the refractive index differs from unity
by only ∼3×10⁻⁴, the vacuum wavelength may be used in ``n(λ)`` with
negligible error (1 part in 10⁹).
"""
function _wave_to_awav(λ_vac)
    λ_μm = λ_vac * 1_000_000 # m → μm for the IUGG formula
    n = _refractive_index_air(λ_μm)
    return λ_vac / n
end

"""
    _awav_to_wave(λ_air)

Air wavelength → vacuum wavelength (both in meters).

Two fixed-point iterations suffice: start with ``λ_vac ≈ λ_air``,
evaluate ``n(λ_vac)``, update ``λ_vac = λ_air × n``, repeat.
Converges to 10⁻¹² in two steps.
"""
function _awav_to_wave(λ_air)
    λ_vac = λ_air
    for _ in 1:2
        n = _refractive_index_air(λ_vac * 1_000_000)
        λ_vac = λ_air * n
    end
    return λ_vac
end

# ── Pipeline helpers ─────────────────────────────────────────────────────────

"""Compile-time check: is this spectral spec for a linear (non-converting) axis?"""
_is_linear(::SpectralSpec{X, P, S, A}) where {X, P, S, A} = (A === Val{:LINEAR})



"""
    _spectral_x_to_world(x, spec::SpectralSpec)

Convert an intermediate coordinate offset ``x = X - X_r`` (X-type SI)
to a world coordinate (S-type SI) for a single spectral axis.

Adds ``X_r`` to recover the absolute X value before the X→P→S chain.
For LOG axes ``x = S_r·ln(S/S_r)`` is used directly.
"""
function _spectral_x_to_world(x, spec::SpectralSpec{X, P, S, A}) where {X, P, S, A}
    # For cross-type: x is the offset X - X_r, so add X_r to get absolute X.
    # For LOG: x = S_r·ln(S/S_r), the LOG formula uses x directly.
    x_abs = (A === Val{:LOG}) ? x : oftype(x, spec.x_r) + x
    p_si = _x_to_p(x_abs, X(), P(), A(), spec)
    return _p_to_s(p_si, S(), spec)
end

"""
    _spectral_world_to_x(world, spec::SpectralSpec)

Convert a world coordinate (S-type SI) to an intermediate coordinate
offset ``x = X - X_r`` (X-type SI) for a single spectral axis.
For LOG axes the offset is ``S_r·ln(S/S_r)``.
"""
function _spectral_world_to_x(world, spec::SpectralSpec{X, P, S, A}) where {X, P, S, A}
    p_si = _s_to_p(world, S(), spec)
    x_abs = _p_to_x(p_si, P(), X(), A(), spec)
    return (A === Val{:LOG}) ? x_abs : x_abs - oftype(x_abs, spec.x_r)
end
