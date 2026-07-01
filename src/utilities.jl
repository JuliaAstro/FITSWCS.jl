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
