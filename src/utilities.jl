@inline _float_type(::Type{T}) where {T<:Real} = float(T)

@inline _promote_float_type(x::Real) = _float_type(typeof(x))

@inline _promote_float_type(x::Real, y::Real) =
    promote_type(_promote_float_type(x), _promote_float_type(y))

@inline _promote_float_type(x::Real, y::Real, z::Real...) =
    promote_type(_promote_float_type(x, y), _promote_float_type(z...))

@inline _halfpi(::Type{T}) where {T<:AbstractFloat} = T(π / 2)
@inline _pi(::Type{T}) where {T<:AbstractFloat} = T(π)

# ──────────────────────────────────────────────────────────────────────────────
# Type-aware numerical tolerances
# ──────────────────────────────────────────────────────────────────────────────

"""
    _boundary_tol(T)

Slack for domain-boundary checks.  Returns a small multiple of machine epsilon
so that points landing fractionally outside the valid domain due to
floating-point roundoff are accepted.  Matches WCSLIB conventions for `Float64`
(≈ 1e-12).
"""
_boundary_tol(::Type{Float64}) = 1e-12
_boundary_tol(::Type{Float32}) = 1f-4
_boundary_tol(T::Type) = T(min(1000 * eps(T), _boundary_tol(Float32)))

"""
    _convergence_tol(T)

Tolerance for iterative convergence (Newton, bisection) and near-pole /
near-zero guards.  Matches WCSLIB conventions for `Float64` (≈ 1e-14).
"""
_convergence_tol(::Type{Float64}) = 1e-14
_convergence_tol(::Type{Float32}) = 1f-5
_convergence_tol(T::Type) = T(min(100 * eps(T), _convergence_tol(Float32)))

# ──────────────────────────────────────────────────────────────────────────────
# Small-integer powers (used by TPV, SIP, and other polynomial distortion)
# ──────────────────────────────────────────────────────────────────────────────

# Julia's ^ uses power_by_squares (generic loop with branching); for the
# exponent range used by TPD and SIP polynomials (max 9) unrolled multiplications
# are 2–3× faster with no dependency cost.

@inline _smallpow2(x) = x * x
@inline _smallpow3(x) = x * x * x
@inline _smallpow4(x) = let x2 = x * x; x2 * x2; end
@inline _smallpow5(x) = let x2 = x * x; x2 * x2 * x; end
@inline _smallpow6(x) = let x2 = x * x; x2 * x2 * x2; end
@inline _smallpow7(x) = let x2 = x * x; x2 * x2 * x2 * x; end
@inline _smallpow8(x) = let x2 = x * x; x4 = x2 * x2; x4 * x4; end
@inline _smallpow9(x) = let x2 = x * x; x4 = x2 * x2; x8 = x4 * x4; x8 * x; end

"""
    _smallpow(x, n)

Return ``x^n`` using unrolled multiplications for ``n ≤ 9``, falling back to
Julia's `^` for larger exponents.  Used by TPV, TPD, and SIP polynomial
evaluation.
"""
@inline function _smallpow(x::Real, n::Int)
    n == 0 && return one(x)
    n == 1 && return x
    n == 2 && return _smallpow2(x)
    n == 3 && return _smallpow3(x)
    n == 4 && return _smallpow4(x)
    n == 5 && return _smallpow5(x)
    n == 6 && return _smallpow6(x)
    n == 7 && return _smallpow7(x)
    n == 8 && return _smallpow8(x)
    n == 9 && return _smallpow9(x)
    return x^n   # fallback for exponents beyond the unrolled range
end

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
