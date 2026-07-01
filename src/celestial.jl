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
    T = _promote_float_type(phi, theta, alpha_p, delta_p, phi_p)
    dphi = phi - phi_p
    sin_th, cos_th = sincos(theta)
    sin_dp, cos_dp = sincos(dphi)
    sin_delp, cos_delp = sincos(delta_p)

    # Paper II, Eq. 2
    delta = asin(clamp(sin_th * sin_delp + cos_th * cos_delp * cos_dp, -one(T), one(T)))
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
    T = _promote_float_type(alpha, delta, alpha_p, delta_p, phi_p)
    da = alpha - alpha_p
    sin_de, cos_de = sincos(delta)
    sin_da, cos_da = sincos(da)
    sin_dp, cos_dp = sincos(delta_p)

    # Paper II, Eq. 5
    theta = asin(clamp(sin_de * sin_dp + cos_de * cos_dp * cos_da, -one(T), one(T)))
    phi = phi_p + atan(-cos_de * sin_da, sin_de * cos_dp - cos_de * sin_dp * cos_da)
    # Normalise phi to [-π, π] (WCSLIB convention).
    phi = mod(phi + _pi(T), 2 * _pi(T)) - _pi(T)
    return phi, theta
end

"""
    native_theta0(proj::AbstractProjection) -> Float64

Return the native latitude of the fiducial point (θ₀) in degrees.
For all zenithal projections this is 90°.
"""
native_theta0(::AZP) = 90.0
native_theta0(::SZP) = 90.0
native_theta0(::TAN) = 90.0
native_theta0(::TPV) = 90.0   # TPV is TAN + polynomial, same native frame
native_theta0(::SIN) = 90.0
native_theta0(::STG) = 90.0
native_theta0(::ARC) = 90.0
native_theta0(::ZEA) = 90.0
native_theta0(::CAR) = 0.0
native_theta0(::CEA) = 0.0
native_theta0(::CYP) = 0.0
native_theta0(::MER) = 0.0
native_theta0(::SFL) = 0.0
native_theta0(::PAR) = 0.0
native_theta0(::MOL) = 0.0
native_theta0(::PCO) = 0.0
native_theta0(::AIT) = 0.0
# Zenithal polynomial
native_theta0(::ZPN) = 90.0
native_theta0(::AIR) = 90.0
# Conic projections: native_theta0 = sigma (the standard parallel).
native_theta0(p::COP) = p.sigma
native_theta0(p::COD) = p.sigma
native_theta0(p::COE) = p.sigma
native_theta0(p::COO) = p.sigma
# BON: wcslib calls prjoff(0,0) so native fiducial is the equator (theta=0).
native_theta0(::BON) = 0.0
# Quadrilateralized spherical cube / HEALPix
native_theta0(::TSC) = 0.0
native_theta0(::CSC) = 0.0
native_theta0(::QSC) = 0.0
native_theta0(::HPX) = 0.0
native_theta0(::XPH) = 90.0
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
    T = _promote_float_type(alpha0, delta0, phi0, theta0, phi_p) # Exclude latpole to avoid promoting to Float64 unnecessarily.
    tol = _convergence_tol(T)

    # Zenithal (theta0 = 90°): unique solution delta_p = delta0, alpha_p = alpha0.
    if abs(theta0 - _halfpi(T)) < tol
        return T(alpha0), T(delta0)
    end

    # General: from Paper II Eq. 2 at the fiducial point,
    # sin(delta0) = sin(theta0)*sin(delta_p) + cos(theta0)*cos(delta_p)*cos(phi0 - phi_p)
    # = A*sin(delta_p) + B*cos(delta_p)  where A=sin(theta0), B=cos(theta0)*cos(phi0-phi_p).
    # This gives: delta_p = asin(rhs) - psi   OR   pi - asin(rhs) - psi
    # where R = hypot(A, B), rhs = sin(delta0)/R, psi = atan(A, B).

    dphi = phi0 - phi_p
    sin_th, cos_th = sincos(T(theta0))
    A = sin_th
    B = cos_th * cos(dphi)
    R = hypot(A, B)

    delta_p = if R < tol
        T(latpole)
    else
        rhs = clamp(sin(T(delta0)) / R, -one(T), one(T))
        psi = atan(B, A)

        s1  = asin(rhs)         # principal asin value in [-π/2, π/2]
        dp1 = s1 - psi          # first candidate
        dp2 = _pi(T) - s1 - psi   # second candidate (supplementary angle)

        dp1 = _reduce_lat(dp1)
        dp2 = _reduce_lat(dp2)

        # Disambiguate: choose candidate closer to latpole.
        (abs(dp1 - T(latpole)) <= abs(dp2 - T(latpole))) ? dp1 : dp2
    end

    # Compute alpha_p from Paper II Eq. 2 (alpha part):
    numer = -cos_th * sin(dphi)
    denom =  sin_th * cos(delta_p) - cos_th * sin(delta_p) * cos(dphi)

    alpha_p = if abs(numer) < tol && abs(denom) < tol
        T(alpha0)   # Degenerate: delta_p = ±90°; alpha_p is arbitrary.
    else
        T(alpha0) - atan(numer, denom)
    end

    return alpha_p, delta_p
end

"""
Reduce angle d (radians) to the valid latitude range [-π/2, π/2].
Uses modular reduction to find the equivalent latitude value.
"""
function _reduce_lat(d::Real)
    T = _float_type(typeof(d))
    pi_T = _pi(T)
    halfpi = _halfpi(T)
    d = mod(d, 2 * pi_T)
    if d <= halfpi
        return d              # already in [0, π/2]
    elseif d <= pi_T
        return pi_T - d       # reflect: (π/2, π] → [0, π/2)  [incorrect sign!]
    elseif d <= pi_T + halfpi
        return -(d - pi_T)    # (π, 3π/2] → [0, -π/2]
    else
        return d - 2 * pi_T   # (3π/2, 2π) → (-π/2, 0)
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
