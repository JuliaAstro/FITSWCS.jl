"""
Spherical rotation between native and celestial coordinate systems.

Implements Paper II (Calabretta & Greisen 2002), Section 2.4 and Eq. 2–5.

All angular inputs and outputs are in **radians**.

Notation following Paper II:
- (φ, θ)        – native spherical longitude and latitude
- (α, δ)        – celestial coordinates (e.g., RA/Dec)
- (αₚ, δₚ)     – celestial coordinates of the native north pole
- φₚ            – native longitude of the celestial pole  (= LONPOLE in degrees)
"""

"""
    native_to_celestial(phi, theta, alpha_p, delta_p, phi_p) -> (alpha, delta)

Convert native spherical coordinates `(phi, theta)` to celestial coordinates
`(alpha, delta)` using the spherical rotation defined by the native north pole
position `(alpha_p, delta_p)` and the celestial-pole native longitude `phi_p`.

All arguments and return values are in radians.

Paper II, Eq. 2.
"""
function native_to_celestial(phi::Real, theta::Real,
                              alpha_p::Real, delta_p::Real, phi_p::Real)
    dphi   = phi - phi_p
    sin_th = sin(theta)
    cos_th = cos(theta)
    sin_dp = sin(dphi)
    cos_dp = cos(dphi)
    sin_delp = sin(delta_p)
    cos_delp = cos(delta_p)

    # Paper II, Eq. 2
    delta = asin(clamp(sin_th * sin_delp + cos_th * cos_delp * cos_dp, -1.0, 1.0))
    alpha = alpha_p + atan(-cos_th * sin_dp, sin_th * cos_delp - cos_th * sin_delp * cos_dp)
    return alpha, delta
end

"""
    celestial_to_native(alpha, delta, alpha_p, delta_p, phi_p) -> (phi, theta)

Inverse of `native_to_celestial`: convert celestial coordinates `(alpha, delta)`
to native spherical coordinates `(phi, theta)`.

All arguments and return values are in radians.

Paper II, Eq. 5.
"""
function celestial_to_native(alpha::Real, delta::Real,
                              alpha_p::Real, delta_p::Real, phi_p::Real)
    da     = alpha - alpha_p
    sin_de = sin(delta)
    cos_de = cos(delta)
    sin_da = sin(da)
    cos_da = cos(da)
    sin_dp = sin(delta_p)
    cos_dp = cos(delta_p)

    # Paper II, Eq. 5
    theta = asin(clamp(sin_de * sin_dp + cos_de * cos_dp * cos_da, -1.0, 1.0))
    phi   = phi_p + atan(-cos_de * sin_da, sin_de * cos_dp - cos_de * sin_dp * cos_da)
    return phi, theta
end

"""
    native_theta0(proj::AbstractProjection) -> Float64

Return the native latitude of the fiducial point (θ₀) in degrees.
For all zenithal projections this is 90°.
"""
native_theta0(::TAN) = 90.0
native_theta0(::SIN) = 90.0
native_theta0(::STG) = 90.0
native_theta0(::ARC) = 90.0
native_theta0(::ZEA) = 90.0
native_theta0(::CAR) = 0.0
native_theta0(::CEA) = 0.0
native_theta0(::AIT) = 0.0
native_theta0(::UnknownProjection) = 90.0

"""
    native_phi0(proj::AbstractProjection) -> Float64

Return the native longitude of the fiducial point (φ₀) in degrees.
For all standard projections this is 0°.
"""
native_phi0(::AbstractProjection) = 0.0

"""
    compute_native_pole(alpha0, delta0, phi0, theta0, phi_p, latpole) -> (alpha_p, delta_p)

Compute the celestial coordinates of the native north pole `(alpha_p, delta_p)`.

Arguments (all in radians):
- `alpha0, delta0`   – CRVAL: celestial coords of the native fiducial point
- `phi0, theta0`     – native coords of the fiducial point
- `phi_p`            – LONPOLE: native longitude of the celestial north pole
- `latpole`          – LATPOLE hint (default π/2) to disambiguate two solutions

Paper II, Section 2.4, Eq. 11.
"""
function compute_native_pole(alpha0::Real, delta0::Real,
                              phi0::Real, theta0::Real, phi_p::Real,
                              latpole::Real = π/2)
    # Zenithal (theta0 = 90°): unique solution delta_p = delta0, alpha_p = alpha0.
    if abs(theta0 - π/2) < 1e-10
        return alpha0, delta0
    end

    # General: from Paper II Eq. 2 at the fiducial point,
    # sin(delta0) = sin(theta0)*sin(delta_p) + cos(theta0)*cos(delta_p)*cos(phi0 - phi_p)
    # = A*sin(delta_p) + B*cos(delta_p)  where A=sin(theta0), B=cos(theta0)*cos(phi0-phi_p).
    # This gives: delta_p = asin(rhs) - psi   OR   pi - asin(rhs) - psi
    # where R = hypot(A, B), rhs = sin(delta0)/R, psi = atan(A, B).

    dphi = phi0 - phi_p
    A    = sin(theta0)
    B    = cos(theta0) * cos(dphi)
    R    = hypot(A, B)

    local delta_p::Float64

    if R < 1e-12
        delta_p = latpole
    else
        rhs = clamp(sin(delta0) / R, -1.0, 1.0)
        # psi is the offset such that A*sin(x) + B*cos(x) = R*sin(x + psi)
        # → psi = atan(B, A)    [atan(y, x) such that sin(psi)=B/R, cos(psi)=A/R]
        # Wait: A*sin(x)+B*cos(x) = R*sin(x + atan(B/A)) only if using atan2.
        # More precisely: A*sin(x)+B*cos(x) = R*sin(x+psi) where psi=atan(B,A).
        # Then x + psi = asin(rhs) or pi - asin(rhs).
        psi = atan(B, A)   # atan2(B, A) — angle such that cos(psi)=A/R, sin(psi)=B/R

        s1  = asin(rhs)         # principal asin value in [-π/2, π/2]
        dp1 = s1 - psi          # first candidate
        dp2 = π - s1 - psi      # second candidate (supplementary angle)

        # Wrap both to valid latitude range [-π/2, π/2].
        # A latitude delta_p = d is the same as d + 2π*k. Valid range: [-π/2, π/2].
        # Equivalently, we require: -1 <= sin(delta_p) <= 1 (always true).
        # But delta_p as a geographic latitude must be in [-π/2, π/2].
        # Reduce modulo 2π, then reflect if outside [-π/2, π/2].
        dp1 = _reduce_lat(dp1)
        dp2 = _reduce_lat(dp2)

        # Disambiguate: choose candidate closer to latpole.
        delta_p = (abs(dp1 - latpole) <= abs(dp2 - latpole)) ? dp1 : dp2
    end

    # Compute alpha_p from Paper II Eq. 2 (alpha part):
    numer = -cos(theta0) * sin(dphi)
    denom =  sin(theta0) * cos(delta_p) - cos(theta0) * sin(delta_p) * cos(dphi)

    alpha_p = if abs(numer) < 1e-10 && abs(denom) < 1e-10
        alpha0   # Degenerate: delta_p = ±90°; alpha_p is arbitrary.
    else
        alpha0 - atan(numer, denom)
    end

    return alpha_p, delta_p
end

"""
Reduce angle d (radians) to the valid latitude range [-π/2, π/2].
Uses modular reduction to find the equivalent latitude value.
"""
function _reduce_lat(d::Real)
    # Map to [0, 2π) first.
    d = mod(d, 2π)
    # Map (π/2, 3π/2) → reflected, (3π/2, 2π) → negative.
    if d <= π/2
        return d              # already in [0, π/2]
    elseif d <= π
        return π - d          # reflect: (π/2, π] → [0, π/2)  [incorrect sign!]
    elseif d <= 3π/2
        return -(d - π)       # (π, 3π/2] → [0, -π/2]
    else
        return d - 2π         # (3π/2, 2π) → (-π/2, 0)
    end
end

"""
    default_lonpole(delta0, theta0) -> Float64

Compute the default LONPOLE (degrees) when not given explicitly.

Paper II, Section 2.4: `φₚ = 0°` when `δ₀ > θ₀`; `φₚ = 180°` otherwise.
"""
function default_lonpole(delta0::Real, theta0::Real)::Float64
    return delta0 > theta0 ? 0.0 : 180.0
end
