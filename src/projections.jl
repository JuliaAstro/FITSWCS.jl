"""
Spherical projection functions.

Each projection pair converts between native spherical coordinates (φ, θ)
expressed in **radians** and intermediate world coordinates (x, y) expressed
in **degrees**.  This matches the FITS WCS Paper II convention.

## References

- Calabretta & Greisen (2002), "Representations of celestial coordinates in
  FITS", Astronomy & Astrophysics, 395, 1077–1122.  (Paper II)
"""

# ──────────────────────────────────────────────────────────────────────────────
# Shared zenithal utilities
# ──────────────────────────────────────────────────────────────────────────────

"""
Shared azimuth angle φ from intermediate coordinates (x, y) [degrees].
Common to all zenithal projections.  Paper II, Eq. 14.
"""
@inline function _phi_zenithal(x::Real, y::Real)
    # WCSLIB defines native longitude as zero at the projection center.
    iszero(x) && iszero(y) && return zero(_promote_float_type(x, y))
    return atan(x, -y)
end

"""
Convert native azimuth + R_θ [degrees] back to (x, y) [degrees].
Common to all zenithal projections.
"""
@inline function _zenithal_native_to_xy(Rtheta_deg::Real, phi_rad::Real)
    x = Rtheta_deg * sin(phi_rad)
    y = -Rtheta_deg * cos(phi_rad)
    return x, y
end

"""
Wrap native longitude to the local projection branch around φ₀ = 0.
"""
@inline _wrap_native_phi(phi::T) where {T <: Real} = phi - T(2π) * round(phi / T(2π))

# ──────────────────────────────────────────────────────────────────────────────
# AZP/SZP default perspective forms
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(::AZP, x, y) -> (phi, theta)

Inverse AZP projection for the default parameter form.
"""
function intermediate_to_native(::AZP, x::Real, y::Real)
    # Default AZP is equivalent to the central gnomonic perspective.
    return intermediate_to_native(TAN(), x, y)
end

"""
    native_to_intermediate(::AZP, phi, theta) -> (x, y)

Forward AZP projection for the default parameter form.
"""
function native_to_intermediate(::AZP, phi::Real, theta::Real)
    # Default AZP is equivalent to TAN; non-default PV parameters are rejected at parse time.
    return native_to_intermediate(TAN(), phi, theta)
end

"""
    intermediate_to_native(::SZP, x, y) -> (phi, theta)

Inverse SZP projection for the default parameter form.
"""
function intermediate_to_native(::SZP, x::Real, y::Real)
    # Default SZP reduces to the same central perspective as TAN.
    return intermediate_to_native(TAN(), x, y)
end

"""
    native_to_intermediate(::SZP, phi, theta) -> (x, y)

Forward SZP projection for the default parameter form.
"""
function native_to_intermediate(::SZP, phi::Real, theta::Real)
    # Default SZP is equivalent to TAN; slant PV parameters are rejected at parse time.
    return native_to_intermediate(TAN(), phi, theta)
end

# ──────────────────────────────────────────────────────────────────────────────
# TAN – Gnomonic projection   (Paper II, Eq. 54–55)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(::TAN, x, y) -> (phi, theta)

Convert intermediate world coordinates `(x, y)` [degrees] to native spherical
coordinates `(phi, theta)` [radians] using the TAN (gnomonic) projection.

Singularity: `theta = 0` (horizon); R_θ → ∞.
"""
function intermediate_to_native(::TAN, x::Real, y::Real)
    T = _promote_float_type(x, y)
    phi = _phi_zenithal(x, y)
    Rth = hypot(x, y) # degrees
    theta = atan(rad2deg(one(T)), Rth)  # = atan2(180/π, Rth) in radians
    return phi, theta
end

"""
    native_to_intermediate(::TAN, phi, theta) -> (x, y)

Convert native spherical coordinates `(phi, theta)` [radians] to intermediate
world coordinates `(x, y)` [degrees] using the TAN projection.

Singularity: `theta ≤ 0`.
"""
function native_to_intermediate(::TAN, phi::Real, theta::Real)
    sth = sin(theta)
    if sth <= 0
        error("TAN projection: theta must be > 0 (got θ = $(rad2deg(theta))°)")
    end
    Rth = rad2deg(cos(theta) / sth)     # degrees
    return _zenithal_native_to_xy(Rth, phi)
end

# ──────────────────────────────────────────────────────────────────────────────
# SIN – Slant orthographic projection   (Paper II, Eq. 48–49)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(::SIN, x, y) -> (phi, theta)

Inverse SIN projection.  For the standard (non-slant) form ξ = η = 0,
this is straightforward.  The slant form solves a quadratic.
"""
function intermediate_to_native(sin_proj::SIN, x::Real, y::Real)
    T = _promote_float_type(x, y)
    xi  = T(sin_proj.xi)
    eta = T(sin_proj.eta)
    oneT = one(T)

    # Convert x, y to radians (projection formulas use dimensionless coords)
    xr = deg2rad(x)
    yr = deg2rad(y)
    r = hypot(xr, yr)

    # WCSLIB fixes the undefined native longitude at the projection center.
    iszero(r) && return zero(T), _halfpi(T)

    if iszero(xi) && iszero(eta)
        # Standard SIN: R_θ = cos(θ)  (in unit sphere coords)
        if r > oneT
            error("SIN projection: point outside valid domain (R_θ = $(r) > 1)")
        end
        theta = acos(r)   # θ = acos(R_θ)
        phi = atan(xr, -yr)
    else
        # Slant SIN: solve quadratic (Paper II, Eq. 49)
        a = xi^2 + eta^2 + oneT
        b = xi*(xr - xi) + eta*(yr - eta)
        c = (xr - xi)^2 + (yr - eta)^2 - oneT
        disc = b^2 - a*c
        if disc < zero(T)
            error("SIN projection: point outside valid domain (discriminant < 0)")
        end
        sth1 = (-b + sqrt(disc)) / a
        sth2 = (-b - sqrt(disc)) / a
        # Choose the solution with θ ≥ θ_0 = 0° (i.e., sin(θ) ≥ 0 preferred)
        # For the standard convention, take the larger sin(θ) value.
        sth = sth1 >= sth2 ? sth1 : sth2
        if abs(sth) > oneT
            sth = clamp(sth, -oneT, oneT)
        end
        theta = asin(sth)
        offset = oneT - sth
        phi = atan(xr - xi*offset, -(yr - eta*offset))
    end
    return phi, theta
end

"""
    native_to_intermediate(::SIN, phi, theta) -> (x, y)

Forward SIN projection.
"""
function native_to_intermediate(sin_proj::SIN, phi::Real, theta::Real)
    xi = sin_proj.xi
    eta = sin_proj.eta
    cth = cos(theta)
    sth = sin(theta)
    # Result in radians, then convert to degrees
    xr = cth * sin(phi) + xi  * (1 - sth)
    yr = -cth * cos(phi) + eta * (1 - sth)
    return rad2deg(xr), rad2deg(yr)
end

# ──────────────────────────────────────────────────────────────────────────────
# STG – Stereographic projection   (Paper II, Eq. 50)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(::STG, x, y) -> (phi, theta)

Inverse STG projection.
"""
function intermediate_to_native(::STG, x::Real, y::Real)
    T = _promote_float_type(x, y)
    phi = _phi_zenithal(x, y)
    Rth = hypot(x, y) # degrees
    # R_θ = 2 * (180/π) * cos(θ) / (1 + sin(θ))
    # → 1 + sin(θ) = 2*(180/π)*cos(θ)/R_θ
    # Solve: let s = sin(θ). Rth/R2D = 2cos(θ)/(1+s) = 2*sqrt(1-s²)/(1+s)
    # = 2*sqrt((1-s)(1+s))/(1+s) = 2*sqrt((1-s)/(1+s))
    # → (Rth/(2*R2D))² = (1-s)/(1+s)
    # → r = Rth/(2*R2D): s = (1-r²)/(1+r²)
    r = deg2rad(Rth) / 2
    sth = (1 - r^2) / (1 + r^2)
    theta = asin(clamp(sth, -one(T), one(T)))
    return phi, theta
end

"""
    native_to_intermediate(::STG, phi, theta) -> (x, y)

Forward STG projection.
"""
function native_to_intermediate(::STG, phi::Real, theta::Real)
    denom = 1 + sin(theta)
    if iszero(denom)
        error("STG projection: singularity at theta = -90°")
    end
    Rth = rad2deg(2 * cos(theta) / denom)
    return _zenithal_native_to_xy(Rth, phi)
end

# ──────────────────────────────────────────────────────────────────────────────
# ARC – Zenithal equidistant projection   (Paper II, Eq. 46)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(::ARC, x, y) -> (phi, theta)

Inverse ARC projection.
"""
function intermediate_to_native(::ARC, x::Real, y::Real)
    T = _promote_float_type(x, y)
    phi = _phi_zenithal(x, y)
    Rth = hypot(x, y) # degrees
    # R_θ = (180/π)*(π/2 − θ)  →  θ = π/2 − R_θ*(π/180)
    theta = _halfpi(T) - deg2rad(Rth)
    return phi, theta
end

"""
    native_to_intermediate(::ARC, phi, theta) -> (x, y)

Forward ARC projection.
"""
function native_to_intermediate(::ARC, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    Rth = rad2deg(_halfpi(T) - theta) # degrees
    return _zenithal_native_to_xy(Rth, phi)
end

# ──────────────────────────────────────────────────────────────────────────────
# ZEA – Lambert zenithal equal-area projection   (Paper II, Eq. 52)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(::ZEA, x, y) -> (phi, theta)

Inverse ZEA projection.
"""
function intermediate_to_native(::ZEA, x::Real, y::Real)
    T = _promote_float_type(x, y)
    phi = _phi_zenithal(x, y)
    Rth = hypot(x, y) # degrees
    # R_θ = 2*(180/π)*sin((π/2 − θ)/2)  →  sin((π/2−θ)/2) = R_θ*π/(360)
    arg = deg2rad(Rth) / 2 # = R_θ * π/360
    if abs(arg) > one(T) + T(1e-12)
        error("ZEA projection: point outside valid domain (|arg| = $(abs(arg)) > 1)")
    end
    arg = clamp(arg, -one(T), one(T))
    theta = _halfpi(T) - 2 * asin(arg)
    return phi, theta
end

"""
    native_to_intermediate(::ZEA, phi, theta) -> (x, y)

Forward ZEA projection.
"""
function native_to_intermediate(::ZEA, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    Rth = rad2deg(2 * sin((_halfpi(T) - theta) / 2)) # degrees
    return _zenithal_native_to_xy(Rth, phi)
end

# ──────────────────────────────────────────────────────────────────────────────
# CAR – Plate carrée / equirectangular projection   (Paper II, Eq. 84)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(::CAR, x, y) -> (phi, theta)

Inverse CAR projection.  x and y are native longitude and latitude in degrees.
"""
function intermediate_to_native(::CAR, x::Real, y::Real)
    phi = deg2rad(x)
    theta = deg2rad(y)
    return phi, theta
end

"""
    native_to_intermediate(::CAR, phi, theta) -> (x, y)

Forward CAR projection.
"""
function native_to_intermediate(::CAR, phi::Real, theta::Real)
    # Use the local longitude branch so inverse transforms prefer nearby pixels.
    phi_w = _wrap_native_phi(phi)
    x = rad2deg(phi_w)
    y = rad2deg(theta)
    return x, y
end

# ──────────────────────────────────────────────────────────────────────────────
# CEA – Cylindrical equal-area projection   (Paper II / wcslib CEA)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(proj::CEA, x, y) -> (phi, theta)

Inverse CEA projection.  `x` and `y` are intermediate coordinates in degrees.
"""
function intermediate_to_native(proj::CEA, x::Real, y::Real)
    T = _promote_float_type(x, y)
    lambda = T(proj.lambda)

    # Convert the linear longitude coordinate directly to native longitude.
    phi = deg2rad(x)

    # Recover latitude from the equal-area ordinate and check the finite domain.
    arg = lambda * deg2rad(y)
    abs(arg) <= one(T) + T(1e-13) ||
        error("CEA projection: point outside valid domain (|lambda*y*pi/180| = $(abs(arg)) > 1)")
    arg = clamp(arg, -one(T), one(T))
    theta = asin(arg)
    return phi, theta
end

"""
    native_to_intermediate(proj::CEA, phi, theta) -> (x, y)

Forward CEA projection.
"""
function native_to_intermediate(proj::CEA, phi::Real, theta::Real)
    # Longitude is linear, using the local branch around the fiducial meridian.
    T = _promote_float_type(phi, theta)
    phi_w = _wrap_native_phi(phi)
    x = rad2deg(phi_w)

    # Latitude maps by the equal-area sine relation with lambda scaling.
    y = rad2deg(sin(theta) / T(proj.lambda))
    return x, y
end

# ──────────────────────────────────────────────────────────────────────────────
# CYP – Cylindrical perspective projection   (Paper II, Eq. 76–77)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(proj::CYP, x, y) -> (phi, theta)

Inverse CYP projection.  `lambda` and `mu` are projection parameters.
"""
function intermediate_to_native(proj::CYP, x::Real, y::Real)
    T = _promote_float_type(x, y)
    lam = T(proj.lambda)
    mu  = T(proj.mu)

    # WCSLIB 8.9 cypx2s: phi_rad = (x / mu) * D2R  (w[1] = 1/mu, x0=0).
    phi = deg2rad(x) / mu

    # WCSLIB: eta = y / (R2D*(lambda+mu)), then atan2d(eta,1) + asind(eta*lambda/√(eta²+1)).
    eta = deg2rad(y) / (lam + mu)
    theta = atan(eta, one(T)) + asin(clamp((eta * lam) / hypot(eta, one(T)), -one(T), one(T)))
    return phi, theta
end

"""
    native_to_intermediate(proj::CYP, phi, theta) -> (x, y)

Forward CYP projection.
"""
function native_to_intermediate(proj::CYP, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    lam = T(proj.lambda)
    mu  = T(proj.mu)

    # WCSLIB 8.9 cyps2x: x = mu * phi_deg  (w[0] = mu, x0=0).
    phi_w = _wrap_native_phi(phi)
    x = rad2deg(mu * phi_w)

    # WCSLIB: y = R2D*(lambda+mu)*sin(theta) / (lambda + cos(theta)).
    denom = lam + cos(theta)
    denom != 0.0 ||
        error("CYP projection: singularity where lambda + cos(theta) = 0")
    y = rad2deg((lam + mu) * sin(theta) / denom)
    return x, y
end

# ──────────────────────────────────────────────────────────────────────────────
# MER – Mercator projection   (Paper II, Eq. 86)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(::MER, x, y) -> (phi, theta)

Inverse MER projection.
"""
function intermediate_to_native(::MER, x::Real, y::Real)
    # Longitude is linear and latitude is the inverse Mercator ordinate.
    phi = deg2rad(x)
    theta = 2 * atan(exp(deg2rad(y))) - _halfpi(_promote_float_type(x, y))
    return phi, theta
end

"""
    native_to_intermediate(::MER, phi, theta) -> (x, y)

Forward MER projection.
"""
function native_to_intermediate(::MER, phi::Real, theta::Real)
    # Mercator is singular at the native poles.
    abs(abs(theta) - π/2) > 1e-14 ||
        error("MER projection: singularity at theta = ±90°")
    phi_w = _wrap_native_phi(phi)
    x = rad2deg(phi_w)
    y = rad2deg(log(tan(_pi(_promote_float_type(phi, theta))/4 + theta/2)))
    return x, y
end

# ──────────────────────────────────────────────────────────────────────────────
# SFL – Sanson-Flamsteed projection   (Paper II, Eq. 88)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(::SFL, x, y) -> (phi, theta)

Inverse SFL projection.
"""
function intermediate_to_native(::SFL, x::Real, y::Real)
    # Latitude is linear; longitude expands by sec(theta).
    theta = deg2rad(y)
    cth = cos(theta)
    abs(cth) > 1e-14 ||
        error("SFL projection: longitude is undefined at theta = ±90°")
    phi = deg2rad(x) / cth
    return phi, theta
end

"""
    native_to_intermediate(::SFL, phi, theta) -> (x, y)

Forward SFL projection.
"""
function native_to_intermediate(::SFL, phi::Real, theta::Real)
    # Longitude contracts by cos(theta), with latitude preserved.
    phi_w = _wrap_native_phi(phi)
    x = rad2deg(phi_w * cos(theta))
    y = rad2deg(theta)
    return x, y
end

# ──────────────────────────────────────────────────────────────────────────────
# PAR – Parabolic projection   (Paper II, Eq. 89)
# ──────────────────────────────────────────────────────────────────────────────

const _PAR_EDGE_TOL = 1e-13

"""
    intermediate_to_native(::PAR, x, y) -> (phi, theta)

Inverse PAR projection.
"""
function intermediate_to_native(::PAR, x::Real, y::Real)
    T = _promote_float_type(x, y)

    # Recover latitude from the parabolic sine ordinate.
    arg = deg2rad(y) / _pi(T)
    abs(arg) <= one(T) ||
        error("PAR projection: point outside valid domain (|y*pi/180/pi| = $(abs(arg)) > 1)")
    theta = 3 * asin(clamp(arg, -one(T), one(T)))

    # Longitude uses the latitude-dependent parabolic scale factor.
    scale = 2 * cos(2 * theta / 3) - 1
    if abs(scale) <= T(_PAR_EDGE_TOL)
        # At the projected pole, WCSLIB accepts only x≈0 and sets phi=0.
        abs(T(x)) <= T(_PAR_EDGE_TOL) && return zero(T), theta
        error("PAR projection: longitude scale is zero")
    end
    phi = deg2rad(x) / scale
    return phi, theta
end

"""
    native_to_intermediate(::PAR, phi, theta) -> (x, y)

Forward PAR projection.
"""
function native_to_intermediate(::PAR, phi::Real, theta::Real)
    # Apply the parabolic longitude scale and sine latitude ordinate.
    phi_w = _wrap_native_phi(phi)
    x = rad2deg(phi_w * (2 * cos(2 * theta / 3) - 1))
    y = rad2deg(π * sin(theta / 3))
    return x, y
end

# ──────────────────────────────────────────────────────────────────────────────
# MOL – Mollweide projection   (Paper II, Eq. 90)
# ──────────────────────────────────────────────────────────────────────────────

const _MOL_MAXITER = 30
const _MOL_TOL = 1e-14
const _MOL_EDGE_TOL = 1e-12

function _mollweide_gamma(theta::Real)
    T = _promote_float_type(theta)

    # Solve 2γ + sin(2γ) = π sin(theta) with Newton iteration.
    abs(abs(theta) - _halfpi(T)) <= T(1e-14) && return copysign(_halfpi(T), theta)
    gamma = theta
    target = _pi(T) * sin(theta)
    for _ in 1:_MOL_MAXITER
        f = 2 * gamma + sin(2 * gamma) - target
        fp = 2 + 2 * cos(2 * gamma)
        !iszero(fp) || break
        step = f / fp
        gamma -= step
        abs(step) <= T(_MOL_TOL) && return gamma
    end
    error("MOL projection: auxiliary angle solve failed to converge")
end

"""
    intermediate_to_native(::MOL, x, y) -> (phi, theta)

Inverse MOL projection.
"""
function intermediate_to_native(::MOL, x::Real, y::Real)
    T = _promote_float_type(x, y)
    sqrt2 = sqrt(T(2))

    # Recover the auxiliary angle gamma from the vertical coordinate.
    sin_gamma = deg2rad(y) / sqrt2
    abs(sin_gamma) <= one(T) + T(_MOL_EDGE_TOL) ||
        error("MOL projection: point outside valid domain (|sin_gamma| = $(abs(sin_gamma)) > 1)")
    gamma = asin(clamp(sin_gamma, -one(T), one(T)))

    # Convert gamma to native latitude and undo the longitude scale.
    theta = asin(clamp((2 * gamma + sin(2 * gamma)) / _pi(T), -one(T), one(T)))
    cos_gamma = cos(gamma)
    if abs(cos_gamma) <= T(1e-14)
        # WCSLIB treats the Mollweide pole as valid only for x≈0.
        abs(T(x)) <= T(_MOL_EDGE_TOL) && return zero(T), theta
        error("MOL projection: longitude is undefined at projected pole")
    end
    phi = deg2rad(x) * _pi(T) / (2 * sqrt2 * cos_gamma)
    return phi, theta
end

"""
    native_to_intermediate(::MOL, phi, theta) -> (x, y)

Forward MOL projection.
"""
function native_to_intermediate(::MOL, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    sqrt2 = sqrt(T(2))
    # Solve the implicit Mollweide latitude equation before projecting.
    phi_w = _wrap_native_phi(phi)
    gamma = _mollweide_gamma(theta)
    x = rad2deg((2 * sqrt2 / _pi(T)) * phi_w * cos(gamma))
    y = rad2deg(sqrt2 * sin(gamma))
    return x, y
end

# ──────────────────────────────────────────────────────────────────────────────
# PCO – Polyconic projection   (Paper II)
# ──────────────────────────────────────────────────────────────────────────────

const _PCO_MAXITER = 64
const _PCO_TOL = 1e-12
const _PCO_SMALL_Y = 1e-4

"""
    intermediate_to_native(::PCO, x, y) -> (phi, theta)

Inverse PCO projection.
"""
function intermediate_to_native(::PCO, x::Real, y::Real)
    T = _promote_float_type(x, y)
    xj = T(x)
    yj = T(y)
    w = abs(yj)
    tol = T(_PCO_TOL)

    # The equator and native poles have direct limiting inverses.
    if w <= tol
        return deg2rad(xj), zero(T)
    elseif abs(w - 90) <= tol
        return zero(T), copysign(_halfpi(T), yj)
    elseif abs(xj) <= tol && w < 90
        return zero(T), deg2rad(yj)
    end

    # Near the equator, use WCSLIB's small-angle approximation to avoid cot(theta).
    theta_d = zero(T)
    ymtheta = zero(T)
    tantheta = zero(T)
    if w < T(_PCO_SMALL_Y)
        w3 = deg2rad(one(T)) / (2 * rad2deg(one(T)))
        theta_d = yj / (one(T) + w3 * xj^2)
        ymtheta = yj - theta_d
        tantheta = tan(deg2rad(theta_d))
    else
        # Solve only for theta with WCSLIB's bounded weighted interval division.
        theta_pos = yj
        theta_neg = zero(T)
        xx = xj^2
        fpos = xx
        fneg = -xx
        for _ in 1:_PCO_MAXITER
            lambda = clamp(fpos / (fpos - fneg), T(0.1), T(0.9))
            theta_d = theta_pos - lambda * (theta_pos - theta_neg)
            ymtheta = yj - theta_d
            tantheta = tan(deg2rad(theta_d))
            f = xx + ymtheta * (ymtheta - 2 * rad2deg(one(T)) / tantheta)

            # Stop once the scalar residue or bracket width reaches WCSLIB tolerance.
            (abs(f) < tol || abs(theta_pos - theta_neg) < tol) && break
            if f > zero(T)
                theta_pos = theta_d
                fpos = f
            else
                theta_neg = theta_d
                fneg = f
            end
        end
    end

    # Reconstruct phi from the solved theta and the eliminated forward equations.
    x1 = rad2deg(one(T)) - ymtheta * tantheta
    y1 = xj * tantheta
    phi_d = (iszero(x1) && iszero(y1)) ? zero(T) : rad2deg(atan(y1, x1)) / sin(deg2rad(theta_d))
    return deg2rad(phi_d), deg2rad(theta_d)
end

"""
    native_to_intermediate(::PCO, phi, theta) -> (x, y)

Forward PCO projection.
"""
function native_to_intermediate(::PCO, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    phi, theta = T(phi), T(theta)
    # Use the equatorial limit to avoid cot(theta) cancellation.
    phi_w = _wrap_native_phi(phi)
    if abs(theta) <= T(1e-14)
        return rad2deg(phi_w), zero(T)
    end

    # Project using the polyconic longitude argument phi*sin(theta).
    s = sin(theta)
    cotθ = cos(theta) / s
    a = phi_w * s
    x = rad2deg(cotθ * sin(a))
    y = rad2deg(theta + cotθ * (1 - cos(a)))
    return x, y
end

# ──────────────────────────────────────────────────────────────────────────────
# AIT – Hammer-Aitoff projection   (Paper II, Eq. 75)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(::AIT, x, y) -> (phi, theta)

Inverse AIT projection.  Valid for all sky positions.

The inverse uses the auxiliary variable z = sqrt(1 - u²/8 - v²/2)
where u = x·(π/180), v = y·(π/180), then recovers sin(θ) and sin(φ/2)
from the forward equations.
"""
function intermediate_to_native(::AIT, x::Real, y::Real)
    T = _promote_float_type(x, y)
    u = deg2rad(T(x))
    v = deg2rad(T(y))
    s = one(T) - u^2 / 8 - v^2 / 2
    if s < -T(2e-13)
        error("AIT projection: point outside valid domain (discriminant s = $s < 0)")
    end
    s = max(s, zero(T))
    z = sqrt(s)
    # g = sqrt(2/(1+z^2)); sinT = v/g; cosT = sqrt(1-sinT^2)
    g = sqrt(2 / (1 + z^2))
    sinT = v / g
    if abs(sinT) > one(T)
        sinT = clamp(sinT, -one(T), one(T))
    end
    theta = asin(sinT)
    cosT = cos(theta)
    # sin(phi/2) = u/(2*g*cosT); cos(phi/2) = z^2/cosT
    if abs(cosT) < T(1e-12)
        # At theta = ±90°, phi is undefined; return phi = 0
        return zero(T), theta
    end
    sinP2 = u / (2 * g * cosT)
    cosP2 = z^2 / cosT
    phi   = 2 * atan(sinP2, cosP2)
    return phi, theta
end

"""
    native_to_intermediate(::AIT, phi, theta) -> (x, y)

Forward AIT projection.  `phi` must be in the range [-π, π]; values outside
this range are wrapped to it.
"""
function native_to_intermediate(::AIT, phi::Real, theta::Real)
    # Wrap phi to (-π, π] so that phi/2 ∈ (-π/2, π/2].
    phi_w = _wrap_native_phi(phi)
    denom = 1 + cos(theta) * cos(phi_w / 2)
    if denom <= 0
        error("AIT projection: degenerate point (1 + cosθ·cos(φ/2) = $denom)")
    end
    g = sqrt(2 / denom)
    x = rad2deg(2 * g * cos(theta) * sin(phi_w / 2))
    y = rad2deg(g * sin(theta))
    return x, y
end

# ──────────────────────────────────────────────────────────────────────────────
# ZPN – Zenithal polynomial projection   (Paper II, Eq. 55)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(proj::ZPN, x, y) -> (phi, theta)

Inverse ZPN projection.  The native radius (in degrees) is evaluated as a
polynomial in the colatitude zd = π/2 − θ.  The polynomial is inverted
numerically.

Paper II, Section 5.2.
"""
function intermediate_to_native(proj::ZPN, x::Real, y::Real)
    T = _promote_float_type(x, y)
    pv = proj.pv

    # Native azimuth from zenithal geometry.
    phi = _phi_zenithal(x, y)

    # Native radius (degrees) and scaled to polynomial variable (radians).
    r_deg = hypot(x, y)

    # Degree of polynomial (length(pv)-1).
    npv = length(pv)

    # Evaluate the polynomial and its derivative to find the inflection/root.
    # If the polynomial is linear (or constant), solve analytically.
    # Otherwise use bisection in [0, π] to find zd such that P(zd) = r_deg/R2D.
    # Where R2D = 180/π.

    target = deg2rad(r_deg)  # dimensionless units (radians)

    # Special case: degree-0 polynomial (constant).
    if npv == 1
        abs(target - pv[1]) < T(1e-12) || error("ZPN projection: point outside valid domain")
        theta = _halfpi(T)
        return phi, theta
    end

    # Degree-1: r = pv[1] + pv[2]*zd  →  zd = (target - pv[1]) / pv[2]
    if npv == 2
        !iszero(pv[2]) || error("ZPN projection: degenerate polynomial (pv[2]=0)")
        zd = (target - T(pv[1])) / T(pv[2])
        abs(zd) <= π || error("ZPN projection: point outside valid domain")
        theta = _halfpi(T) - zd
        return phi, theta
    end

    # Degree-2: solve analytically to avoid bisection issues at the root.
    if npv == 3
        c = T(pv[1]) - target
        b = T(pv[2])
        a = T(pv[3])
        if abs(a) < T(1e-15)
            # Degenerate to degree-1.
            abs(b) > T(1e-15) || error("ZPN projection: degenerate polynomial")
            zd = -c / b
        else
            disc = b^2 - 4*a*c
            disc >= -T(1e-14) || error("ZPN projection: point outside valid domain")
            disc = max(zero(T), disc)
            sqrt_disc = sqrt(disc)
            # Pick the smaller non-negative root (closest to the projection center).
            zd1 = (-b + sqrt_disc) / (2*a)
            zd2 = (-b - sqrt_disc) / (2*a)
            if zd1 >= 0 && (zd1 <= zd2 || zd2 < 0)
                zd = zd1
            else
                zd = zd2
            end
            zd >= 0 || error("ZPN projection: point outside valid domain")
        end
        abs(zd) <= π || error("ZPN projection: point outside valid domain")
        theta = _halfpi(T) - zd
        return phi, theta
    end

    # General case: bisection on [0, π] to find zd such that P(zd) = target.
    # First check the endpoints and the zd=0 case.
    zfn(z) = T(evalpoly(z, pv)) - target

    # If target = 0 and pv[1] = 0, zd = 0 is a root (reference point).
    p_at_0 = zfn(zero(T))
    if abs(p_at_0) < T(1e-14)
        zd = zero(T)
        theta = _halfpi(T) - zd
        return phi, theta
    end

    # Find the first maximum of P (where P' turns negative) as the
    # upper limit.  The closest root to zd=0 is in [0, zd_max].
    dpv = [T(k) * T(pv[k+1]) for k in 1:npv-1]
    _dpoly(z) = evalpoly(z, dpv)
    zd_max = _pi(T)

    if length(dpv) >= 1
        dp0 = _dpoly(zero(T))
        if dp0 < 0
            zd_max = _pi(T)  # polynomial decreases monotonically — use full range
        else
            a, b = zero(T), _pi(T)
            fa, fb = dp0, _dpoly(b)
            if fa * fb < 0
                for _ in 1:64
                    m = (a + b) / 2
                    fm = _dpoly(m)
                    if fm * fa < 0
                        b, fb = m, fm
                    else
                        a, fa = m, fm
                    end
                    (b - a) < T(1e-14) && break
                end
                zd_max = (a + b) / 2
            end
        end
    end

    # Bisect on [0, zd_max] to find the smallest positive root.
    a, b = zero(T), zd_max
    pa = zfn(a)
    pb = zfn(b)

    # If pa is zero we already returned above; check pb.
    if abs(pb) < T(1e-14)
        zd = b
        theta = _halfpi(T) - zd
        return phi, theta
    end

    pa * pb > 0 && error("ZPN projection: point outside the valid domain of the polynomial")

    for _ in 1:64
        m  = (a + b) / 2
        pm = zfn(m)
        if abs(pm) < T(1e-15)
            zd = m
            theta = _halfpi(T) - zd
            return phi, theta
        end
        if pm * pa < 0
            b, pb = m, pm
        else
            a, pa = m, pm
        end
        (b - a) < T(1e-15) && break
    end

    zd = (a + b) / 2
    theta = _halfpi(T) - zd
    return phi, theta
end

"""
    native_to_intermediate(proj::ZPN, phi, theta) -> (x, y)

Forward ZPN projection.
"""
function native_to_intermediate(proj::ZPN, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    pv = proj.pv

    # Colatitude in radians.
    zd = _halfpi(T) - theta

    # evalpoly uses Horner's method internally and is faster than a manual loop.
    r = T(evalpoly(zd, pv))

    # Convert from radians (polynomial variable) to degrees (projection plane).
    r_deg = rad2deg(r)

    return _zenithal_native_to_xy(r_deg, phi)
end

# ──────────────────────────────────────────────────────────────────────────────
# AIR – Airy projection   (Paper II, Section 5.5, Eq. 30–31)
# ──────────────────────────────────────────────────────────────────────────────

const _AIR_TOL = 1e-14

"""
    intermediate_to_native(proj::AIR, x, y) -> (phi, theta)

Inverse AIR projection.  Bisection is used to invert the Airy formula.

Paper II, Section 5.5.
"""
function intermediate_to_native(proj::AIR, x::Real, y::Real)
    T = _promote_float_type(x, y)
    R2D = rad2deg(one(T)) # = 180 / π
    D2R = deg2rad(one(T)) # = π/180

    phi = _phi_zenithal(x, y)
    r_deg = hypot(x, y)

    theta_b = T(proj.theta_b)

    # Pre-compute the w constants from wcslib airset.
    # xi_b = (90 - theta_b) / 2  in radians
    xi_b = D2R * (90 - theta_b) / 2

    # w[0] = 2*R2D (the overall scale factor)
    w0 = 2 * R2D

    if abs(xi_b) < T(_AIR_TOL)
        # theta_b = 90°: w[1] = -0.5 (limit of log(cos(xi))/... as xi→0)
        w1 = T(-0.5)
        w2 = one(T)  # ensures w3 is finite
    else
        cos_xib = cos(xi_b)
        sin_xib = sin(xi_b)
        # w[1] = log(cos_xib) * cos_xib^2 / sin_xib^2
        w1 = log(cos_xib) * cos_xib^2 / sin_xib^2
        # w[2] = 0.5 - w1
        w2 = T(0.5) - w1
    end

    # w[3] = w[0] * w[2]
    w3 = w0 * w2

    # The target: r_deg / w[0].
    # AIR inverse: find xi such that -(log(cos(xi))/tan(xi) + w1*tan(xi)) = r_target
    r_target = r_deg / w0

    # Reference point: r_target = 0 → xi = 0, theta = π/2.
    if abs(r_target) < T(_AIR_TOL)
        return phi, _halfpi(T)
    end

    # Bisect on [0, π/2) to find xi.
    # At xi → 0: the function → 0.
    # At xi → π/2: the function → ∞.
    a, b = zero(T), _halfpi(T) * (1 - T(1e-10))

    for _ in 1:64
        m  = (a + b) / 2
        cos_m = cos(m)
        sin_m = sin(m)
        if abs(sin_m) < T(1e-15) || abs(cos_m) < T(1e-15)
            b = m
            continue
        end
        # g(xi) = -(log(cos(xi))/tan(xi) + w1*tan(xi)) - r_target
        # g(0) = -r_target < 0, g increases monotonically with xi.
        gm = -(log(cos_m) / (sin_m / cos_m) + w1 * (sin_m / cos_m)) - r_target
        if gm < 0
            a = m    # root is to the right
        else
            b = m    # root is to the left
        end
        (b - a) < T(1e-15) && break
    end

    xi = (a + b) / 2
    theta = _halfpi(T) - 2 * xi

    return phi, theta
end

"""
    native_to_intermediate(proj::AIR, phi, theta) -> (x, y)

Forward AIR projection.
"""
function native_to_intermediate(proj::AIR, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    R2D = rad2deg(one(T)) # = T(180 / π)
    D2R = deg2rad(one(T)) # = T(π / 180)

    theta_b = T(proj.theta_b)
    xi_b = D2R * (90 - theta_b) / 2

    w0 = 2 * R2D

    if abs(xi_b) < T(_AIR_TOL)
        w1 = T(-0.5)
    else
        cos_xib = cos(xi_b)
        sin_xib = sin(xi_b)
        w1 = log(cos_xib) * cos_xib^2 / sin_xib^2
    end
    w2 = T(0.5) - w1
    w3 = w0 * w2

    # xi = (π/2 - theta) / 2
    xi = (_halfpi(T) - theta) / 2

    if abs(xi) < T(_AIR_TOL)
        r_deg = xi * w3
    else
        cos_xi = cos(xi)
        tan_xi = tan(xi)
        r_deg = -w0 * (log(cos_xi) / tan_xi + w1 * tan_xi)
    end

    return _zenithal_native_to_xy(r_deg, phi)
end

# ──────────────────────────────────────────────────────────────────────────────
# Conic projection utilities
# ──────────────────────────────────────────────────────────────────────────────

# Shared inverse helper for all four conic projections:
# given (x, y) and the Y0 offset and cone half-angle C, return
# (r_deg, alpha_rad), where r_deg is unsigned if sigma ≥ 0.
# Caller is responsible for sign convention on r and computing phi/theta.

@inline function _conic_xy_to_r_alpha(x::T, y::T, Y0_deg::T, sigma_rad::T) where {T <: Real}
    dy = Y0_deg - y
    r = hypot(x, dy)
    sigma_rad < 0 && (r = -r)
    alpha = atan(x, dy)
    return r, alpha
end

# ──────────────────────────────────────────────────────────────────────────────
# COP – Conic perspective projection   (Paper II, Section 6.1)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(proj::COP, x, y) -> (phi, theta)

Inverse COP projection.

Paper II, Section 6.1, Eq. 57–58.
"""
function intermediate_to_native(proj::COP, x::Real, y::Real)
    T = _promote_float_type(x, y)
    R2D = rad2deg(one(T)) # = T(180 / π)

    sigma_rad = T(deg2rad(proj.sigma))
    delta_rad = T(deg2rad(proj.delta))

    # Cone constant C = sin(sigma)
    C = sin(sigma_rad)
    C != 0 || error("COP: sigma = 0 is a degenerate case (use CAR or SFL instead)")

    # Y0 = R2D * cos(delta) * cot(sigma) = R2D * cos(delta) * cos(sigma) / sin(sigma)
    Y0_deg = R2D * cos(delta_rad) * cos(sigma_rad) / sin(sigma_rad)

    r_deg, alpha_rad = _conic_xy_to_r_alpha(T(x), T(y), Y0_deg, sigma_rad)

    phi = alpha_rad / C

    # theta = sigma + atan(cot(sigma) - r_deg/(R2D * cos(delta)))
    theta = sigma_rad + atan(cos(sigma_rad) / sin(sigma_rad) -
                              deg2rad(r_deg) / cos(delta_rad))
    return phi, theta
end

"""
    native_to_intermediate(proj::COP, phi, theta) -> (x, y)

Forward COP projection.
"""
function native_to_intermediate(proj::COP, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    R2D = rad2deg(one(T)) # = T(180 / π)

    sigma_rad = T(deg2rad(proj.sigma))
    delta_rad = T(deg2rad(proj.delta))

    C = sin(sigma_rad)
    C != 0 || error("COP: sigma = 0 is degenerate")

    Y0_deg = R2D * cos(delta_rad) * cos(sigma_rad) / sin(sigma_rad)

    alpha_rad = C * _wrap_native_phi(phi)
    t = theta - sigma_rad
    cos_t = cos(t)
    abs(cos_t) > T(1e-15) || error("COP: singularity at theta = sigma ± 90°")

    r_deg = Y0_deg - R2D * cos(delta_rad) * tan(t)

    x = r_deg * sin(alpha_rad)
    y = -r_deg * cos(alpha_rad) + Y0_deg
    return x, y
end

# ──────────────────────────────────────────────────────────────────────────────
# COD – Conic equidistant projection   (Paper II, Section 6.2)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(proj::COD, x, y) -> (phi, theta)

Inverse COD projection.

Paper II, Section 6.2, Eq. 59–60.
"""
function intermediate_to_native(proj::COD, x::Real, y::Real)
    T = _promote_float_type(x, y)
    R2D = rad2deg(one(T)) # = T(180 / π)

    sigma_rad = T(deg2rad(proj.sigma))
    delta_rad = T(deg2rad(proj.delta))

    # Cone constant C = sin(sigma) * sinc(delta) (dimensionless)
    C = if abs(delta_rad) < T(1e-12)
        sin(sigma_rad)   # limit of sin(sigma)*sin(delta)/delta as delta→0
    else
        sin(sigma_rad) * sin(delta_rad) / delta_rad
    end
    abs(C) > T(1e-15) || error("COD: degenerate case (C ≈ 0)")

    Y0_deg = R2D * cos(delta_rad) * cos(sigma_rad) / C

    r_deg, alpha_rad = _conic_xy_to_r_alpha(T(x), T(y), Y0_deg, sigma_rad)

    phi   = alpha_rad / C
    theta = sigma_rad + delta_rad - deg2rad(r_deg) / (R2D / R2D)
    # Simplify: theta_rad = sigma_rad + delta_rad * (...) is the COD form.
    # From wcslib: theta = w3 - r; w3 = Y0_d + sigma_d (both in degrees)
    # So theta_rad = deg2rad(Y0_deg + rad2deg(sigma_rad) - r_deg)
    theta = deg2rad(Y0_deg + rad2deg(sigma_rad) - r_deg)

    return phi, theta
end

"""
    native_to_intermediate(proj::COD, phi, theta) -> (x, y)

Forward COD projection.
"""
function native_to_intermediate(proj::COD, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    R2D = rad2deg(one(T)) # = T(180 / π)

    sigma_rad = T(deg2rad(proj.sigma))
    delta_rad = T(deg2rad(proj.delta))

    C = if abs(delta_rad) < T(1e-12)
        sin(sigma_rad)
    else
        sin(sigma_rad) * sin(delta_rad) / delta_rad
    end
    abs(C) > T(1e-15) || error("COD: degenerate case (C ≈ 0)")

    Y0_deg = R2D * cos(delta_rad) * cos(sigma_rad) / C

    alpha_rad = C * _wrap_native_phi(phi)
    # r = (Y0_d + sigma_d) - theta_d
    r_deg = Y0_deg + rad2deg(sigma_rad) - rad2deg(theta)

    x = r_deg * sin(alpha_rad)
    y = -r_deg * cos(alpha_rad) + Y0_deg
    return x, y
end

# ──────────────────────────────────────────────────────────────────────────────
# COE – Conic equal-area projection   (Paper II, Section 6.3)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(proj::COE, x, y) -> (phi, theta)

Inverse COE projection.

Paper II, Section 6.3, Eq. 61–62.
"""
function intermediate_to_native(proj::COE, x::Real, y::Real)
    T = _promote_float_type(x, y)
    R2D = rad2deg(one(T)) # = T(180 / π)

    sigma_rad = T(deg2rad(proj.sigma))
    delta_rad = T(deg2rad(proj.delta))
    theta1_rad = sigma_rad - delta_rad
    theta2_rad = sigma_rad + delta_rad

    C = (sin(theta1_rad) + sin(theta2_rad)) / 2
    abs(C) > T(1e-15) || error("COE: degenerate case (C ≈ 0)")

    chi = R2D / C
    psi = one(T) + sin(theta1_rad) * sin(theta2_rad)
    Y0_deg = chi * sqrt(max(zero(T), psi - 2 * C * sin(sigma_rad)))

    # w[6] = chi^2 * psi,  w[7] = C / (2 * R2D^2)
    w6 = chi^2 * psi
    w7 = C / (2 * R2D^2)

    r_deg, alpha_rad = _conic_xy_to_r_alpha(T(x), T(y), Y0_deg, sigma_rad)

    phi = alpha_rad / C

    arg = (w6 - r_deg^2) * w7
    abs(arg) <= one(T) + T(1e-12) || error("COE: point outside valid domain")
    theta = asin(clamp(arg, -one(T), one(T)))

    return phi, theta
end

"""
    native_to_intermediate(proj::COE, phi, theta) -> (x, y)

Forward COE projection.
"""
function native_to_intermediate(proj::COE, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    R2D = rad2deg(one(T)) # = T(180 / π)

    sigma_rad = T(deg2rad(proj.sigma))
    delta_rad = T(deg2rad(proj.delta))
    theta1_rad = sigma_rad - delta_rad
    theta2_rad = sigma_rad + delta_rad

    C = (sin(theta1_rad) + sin(theta2_rad)) / 2
    abs(C) > T(1e-15) || error("COE: degenerate case")

    chi = R2D / C
    psi = one(T) + sin(theta1_rad) * sin(theta2_rad)
    Y0_deg = chi * sqrt(max(zero(T), psi - 2 * C * sin(sigma_rad)))

    alpha_rad = C * _wrap_native_phi(phi)

    if theta == -_halfpi(T)
        # Southern pole: r = w[8] = chi * sqrt(psi + 2*C)
        r_deg = chi * sqrt(max(zero(T), psi + 2 * C))
    else
        r_deg = chi * sqrt(max(zero(T), psi - 2 * C * sin(theta)))
    end

    x = r_deg * sin(alpha_rad)
    y = -r_deg * cos(alpha_rad) + Y0_deg
    return x, y
end

# ──────────────────────────────────────────────────────────────────────────────
# COO – Conic orthomorphic projection   (Paper II, Section 6.4)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(proj::COO, x, y) -> (phi, theta)

Inverse COO projection.

Paper II, Section 6.4, Eq. 63–64.
"""
function intermediate_to_native(proj::COO, x::Real, y::Real)
    T = _promote_float_type(x, y)
    R2D = rad2deg(one(T)) # = T(180 / π)

    sigma_rad = T(deg2rad(proj.sigma))
    delta_rad = T(deg2rad(proj.delta))
    theta1_rad = sigma_rad - delta_rad
    theta2_rad = sigma_rad + delta_rad

    # Cone constant C.
    C = if abs(delta_rad) < T(1e-14)
        sin(theta1_rad)
    else
        # tau1 = (90 - theta1) / 2 in radians,  tan1 = tan(tau1)
        tau1 = (_halfpi(T) - theta1_rad) / 2
        tau2 = (_halfpi(T) - theta2_rad) / 2
        log(cos(theta2_rad) / cos(theta1_rad)) / log(tan(tau2) / tan(tau1))
    end

    tau1 = (_halfpi(T) - theta1_rad) / 2
    tan1 = tan(tau1)
    psi  = R2D * cos(theta1_rad) / (C * tan1^C)
    Y0_deg = psi * tan((_halfpi(T) - sigma_rad) / 2)^C

    r_deg, alpha_rad = _conic_xy_to_r_alpha(T(x), T(y), Y0_deg, sigma_rad)

    phi = alpha_rad / C

    if abs(r_deg) < T(1e-14)
        theta = C < 0 ? -_halfpi(T) : error("COO: singularity at origin")
    else
        theta = _halfpi(T) - 2 * atan((r_deg / psi)^(one(T) / C))
    end

    return phi, theta
end

"""
    native_to_intermediate(proj::COO, phi, theta) -> (x, y)

Forward COO projection.
"""
function native_to_intermediate(proj::COO, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    R2D = rad2deg(one(T)) # = T(180 / π)

    sigma_rad = T(deg2rad(proj.sigma))
    delta_rad = T(deg2rad(proj.delta))
    theta1_rad = sigma_rad - delta_rad
    theta2_rad = sigma_rad + delta_rad

    C = if abs(delta_rad) < T(1e-14)
        sin(theta1_rad)
    else
        tau1_l = (_halfpi(T) - theta1_rad) / 2
        tau2_l = (_halfpi(T) - theta2_rad) / 2
        log(cos(theta2_rad) / cos(theta1_rad)) / log(tan(tau2_l) / tan(tau1_l))
    end

    tau1 = (_halfpi(T) - theta1_rad) / 2
    tan1 = tan(tau1)
    psi  = R2D * cos(theta1_rad) / (C * tan1^C)
    Y0_deg = psi * tan((_halfpi(T) - sigma_rad) / 2)^C

    alpha_rad = C * _wrap_native_phi(phi)

    if theta == -_halfpi(T)
        C < 0 || error("COO: singularity at southern pole when C ≥ 0")
        r_deg = zero(T)
    else
        # r_deg = psi * tan((π/2 - theta)/2)^C
        r_deg = psi * tan((_halfpi(T) - theta) / 2)^C
    end

    x = r_deg * sin(alpha_rad)
    y = -r_deg * cos(alpha_rad) + Y0_deg
    return x, y
end

# ──────────────────────────────────────────────────────────────────────────────
# BON – Bonne's projection   (Paper II, Section 7.4, Eq. 70)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(proj::BON, x, y) -> (phi, theta)

Inverse BON (Bonne) projection.

Paper II, Section 7.4, Eq. 70.
"""
function intermediate_to_native(proj::BON, x::Real, y::Real)
    T = _promote_float_type(x, y)
    R2D = rad2deg(one(T)) # = T(180 / π)

    theta1 = T(proj.theta1)

    # Degenerate case: theta1 == 0 degenerates to SFL.
    if theta1 == 0
        return intermediate_to_native(SFL(), x, y)
    end

    theta1_rad = deg2rad(theta1)

    # Y0 = R2D * (cot(theta1) + theta1_rad)
    Y0_deg = R2D * (cos(theta1_rad) / sin(theta1_rad) + theta1_rad)

    # wcslib bonx2s: dy = w[2] - (y + prj->y0) with prj->y0=0 → dy = Y0_deg - y.
    dy = Y0_deg - y
    r  = hypot(x, dy)
    theta1 < 0 && (r = -r)

    alpha_rad = atan(x, dy)

    theta_deg = Y0_deg - r  # since w[1]=1, theta_deg = Y0_deg - r_deg
    theta = deg2rad(theta_deg)
    cos_theta = cos(theta)

    if abs(cos_theta) < T(1e-12)
        phi = zero(T)
    else
        # phi_rad = alpha_rad * r_rad / cos_theta  where r_rad = r_deg/R2D
        phi = alpha_rad * (r / R2D) / cos_theta
    end

    return phi, theta
end

"""
    native_to_intermediate(proj::BON, phi, theta) -> (x, y)

Forward BON (Bonne) projection.
"""
function native_to_intermediate(proj::BON, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    R2D = rad2deg(one(T)) # = T(180 / π)

    theta1 = T(proj.theta1)

    # Degenerate case: theta1 == 0 degenerates to SFL.
    if theta1 == 0
        return native_to_intermediate(SFL(), phi, theta)
    end

    theta1_rad = deg2rad(theta1)

    Y0_deg = R2D * (cos(theta1_rad) / sin(theta1_rad) + theta1_rad)
    r_deg  = Y0_deg - rad2deg(theta)

    # alpha_rad = R2D * phi_rad * cos(theta) / r_deg
    # Derived from wcslib bons2x: s = r0*phi_deg; alpha_deg = s*cos(theta)/r_deg;
    # alpha_rad = deg2rad(alpha_deg) = R2D * phi_rad * cos(theta) / r_deg.
    # wcslib calls prjoff(0, 0) so prj->y0 = 0, and the y-offset is
    # y0_eff = prj->y0 - w[2] = -Y0_deg, giving y = -r*cos(alpha) + Y0_deg.
    if abs(r_deg) < T(1e-14)
        alpha_rad = zero(T)
    else
        alpha_rad = phi * cos(theta) * R2D / r_deg
    end

    x = r_deg * sin(alpha_rad)
    y = -r_deg * cos(alpha_rad) + Y0_deg
    return x, y
end

# ──────────────────────────────────────────────────────────────────────────────
# Quadrilateralized spherical cube utilities (shared by TSC, CSC, QSC)
# ──────────────────────────────────────────────────────────────────────────────

# Unit vector from (phi, theta) in radians → (l, m, n) = (x, y, z) in wcslib.
@inline function _sphere_to_uvec(phi::T, theta::T) where {T <: Real}
    cos_theta = cos(theta)
    return cos_theta * cos(phi), cos_theta * sin(phi), sin(theta)
end

# wcslib face conventions (tscs2x, cscs2x, qscs2x):
#   Face 0: north pole  (+z face): l = +z-dominant, xf = m/l (=y/z), yf = -l_face (= -x/z wait..)
# Actually wcslib uses (l,m,n) = (z,x,y) in its own convention... let me use its exact convention:
#   l = cos(theta)*cos(phi),  m = cos(theta)*sin(phi),  n = sin(theta)
# Then wcslib face selection (from tscs2x):
#   face 0: n dominates, n>0  → xf = m/n,  yf = l/n
#   face 1: l dominates, l>0  → xf = m/l,  yf = n/l
#   face 2: m dominates, m>0  → xf = -l/m, yf = n/m
#   face 3: l dominates, l<0  → xf = m/l,  yf = -n/l
#   face 4: m dominates, m<0  → xf = -l/m, yf = -n/m
#   face 5: n dominates, n<0  → xf = m/n,  yf = -l/n

# Face scale in degrees per face-coordinate unit (w[0] in WCSLIB).
const _FACESCALE = 45.0

# Face unit offsets (multiply by _FACESCALE to get degrees on canvas).
# WCSLIB 8.9: x = 45*(xf + x0u[face]), y = 45*(yf + y0u[face]).
const _FACE_X0U = (0.0, 0.0, 2.0, 4.0, 6.0, 0.0)
const _FACE_Y0U = (2.0, 0.0, 0.0, 0.0, 0.0, -2.0)

# Determine face and (xf, yf) ∈ [-1, 1]² from unit vector (l, m, n).
# WCSLIB 8.9 forward face detection (tscs2x, cscs2x, qscs2x).
function _xyz_to_cube_face(l::T, m::T, n::T) where {T <: Real}
    face = 0
    zeta = n
    if l > zeta; face = 1; zeta = l; end
    if m > zeta; face = 2; zeta = m; end
    if -l > zeta; face = 3; zeta = -l; end
    if -m > zeta; face = 4; zeta = -m; end
    if -n > zeta; face = 5; zeta = -n; end

    if face == 1
        xf =  m / zeta
        yf =  n / zeta
    elseif face == 2
        xf = -l / zeta
        yf =  n / zeta
    elseif face == 3
        xf = -m / zeta
        yf =  n / zeta
    elseif face == 4
        xf =  l / zeta
        yf =  n / zeta
    elseif face == 5
        xf =  m / zeta
        yf =  l / zeta
    else  # face == 0
        xf =  m / zeta
        yf = -l / zeta
    end
    return face, xf, yf
end

# Recover unit vector (l, m, n) from face and (xf, yf).
# WCSLIB 8.9 inverse face-to-sphere (tscx2s, cscx2s, qscx2s).
function _face_to_uvec(face::Int, xf::T, yf::T) where {T <: Real}
    if face == 1
        l =  one(T); m =  xf;     n =  yf
    elseif face == 2
        l = -xf;     m =  one(T); n =  yf
    elseif face == 3
        l = -one(T); m = -xf;     n =  yf
    elseif face == 4
        l =  xf;     m = -one(T); n =  yf
    elseif face == 5
        l =  yf;     m = -xf;     n = -one(T)
    else  # face == 0
        l = -yf;     m =  xf;     n =  one(T)
    end
    r = hypot(l, m, n)
    return l / r, m / r, n / r
end

# Convert canvas (x, y) [degrees] to face and (xf, yf) ∈ [-1, 1]².
# WCSLIB 8.9 inverse face determination (tscx2s, cscx2s, qscx2s).
# scale = 45.0 (w[0]) — face coordinate scaling in degrees/unit.
function _cube_xy_to_face(x_deg::T, y_deg::T, scale::T) where {T <: Real}
    TOL = T(1e-13)
    xf = x_deg / scale
    yf = y_deg / scale

    # Handle negative face wrapping (xf < -1 → add 8 face units).
    if xf < -one(T); xf += T(8); end

    if yf > one(T) + TOL
        face = 0; yf -= T(2)
    elseif yf < -one(T) - TOL
        face = 5; yf += T(2)
    elseif xf > T(5) + TOL
        face = 4; xf -= T(6)
    elseif xf > T(3) + TOL
        face = 3; xf -= T(4)
    elseif xf > one(T) + TOL
        face = 2; xf -= T(2)
    else
        face = 1
    end
    return face, xf, yf
end

# ──────────────────────────────────────────────────────────────────────────────
# TSC – Tangential spherical cube   (Paper II, Section 8.1)
# ──────────────────────────────────────────────────────────────────────────────

"""
    native_to_intermediate(::TSC, phi, theta) -> (x, y)

Forward TSC projection.  Paper II, Section 8.1.
"""
function native_to_intermediate(::TSC, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    l, m, n = _sphere_to_uvec(T(phi), T(theta))
    face, xf, yf = _xyz_to_cube_face(l, m, n)
    s = T(_FACESCALE)
    x = s * (xf + T(_FACE_X0U[face + 1]))
    y = s * (yf + T(_FACE_Y0U[face + 1]))
    return x, y
end

"""
    intermediate_to_native(::TSC, x, y) -> (phi, theta)

Inverse TSC projection.  Paper II, Section 8.1.
"""
function intermediate_to_native(::TSC, x::Real, y::Real)
    T = _promote_float_type(x, y)
    face, xf, yf = _cube_xy_to_face(T(x), T(y), T(_FACESCALE))
    l, m, n = _face_to_uvec(face, xf, yf)
    phi = atan(m, l)
    theta = asin(clamp(n, -one(T), one(T)))
    return phi, theta
end

# ──────────────────────────────────────────────────────────────────────────────
# CSC – COBE quadrilateralized spherical cube   (Paper II, Section 8.2)
# ──────────────────────────────────────────────────────────────────────────────

# Inverse (x2s) polynomial: maps (chi, psi) face coords to (xf, yf).
# chi is the primary face coordinate, psi is the secondary.
# From wcslib cscx2s: p[i][j] where result = sum chi^i * psi^j
const _CSC_INV_P = (
    (-0.27292696, -0.07629969, -0.22797056,  0.54852384, -0.62930065,  0.25795794,  0.02584375),
    (-0.02819452, -0.01471565,  0.48051509, -1.74114454,  1.71547508, -0.53022337),
    ( 0.27058160, -0.56800938,  0.30803317,  0.98938102, -0.83180469),
    (-0.60441560,  1.50880086, -0.93678576,  0.08693841),
    ( 0.93412077, -1.41601920,  0.33887446),
    (-0.63915306,  0.52032238),
    ( 0.14381585,),
)

# Forward (s2x) polynomial constants: maps (xf, yf) face coords to (chi, psi).
# From wcslib cscs2x notation, the formula for chi is:
#   chi = xf * (gstar + xf^2*(mm*yf^2 - gamma) - yf^2*(d0 + yf^2*d1))
#         + correction polynomial
# and symmetrically for psi with xf ↔ yf.
# The correction polynomial uses p00, p10, p01, p11, p20, p02.
const _CSC_FWD = (
    gstar  =  1.37484847732,
    mm     =  0.004869491981,
    gamma  = -0.13161671474,
    omega1 = -0.159596235474,
    d0     =  0.0759196200467,
    d1     = -0.0217762490699,
    p00    =  0.141189631152,
    p10    =  0.0809701286525,
    p01    = -0.281528535557,
    p11    =  0.15384112876,
    p20    = -0.178251207466,
    p02    =  0.106959469314,
)

# Evaluate the CSC forward (xf,yf → chi) single-axis polynomial.
# Returns chi for axis 1 (xf is primary); swap arguments for psi.
function _csc_fwd_axis(chi::T, psi::T) where {T <: Real}
    # WCSLIB 8.9 cscs2x: maps face coords (chi=m/l, psi=n/l) → projected coords.
    chi2 = chi^2
    psi2 = psi^2
    chi2co = one(T) - chi2
    psi2co = one(T) - psi2
    chipsi = abs(chi*psi)
    chi2psi2 = chipsi > T(1e-16) ? chi2*psi2 : zero(T)
    chi4 = chi2 > T(1e-16) ? chi2^2 : zero(T)
    psi4 = psi2 > T(1e-16) ? psi2^2 : zero(T)

    gstar  = T(1.37484847732)
    mm     = T(0.004869491981)
    gamma  = T(-0.13161671474)
    omega1 = T(-0.159596235474)
    d0     = T(0.0759196200467)
    d1     = T(-0.0217762490699)
    c00 = T(0.141189631152)
    c10 = T(0.0809701286525)
    c01 = T(-0.281528535557)
    c11 = T(0.15384112876)
    c20 = T(-0.178251207466)
    c02 = T(0.106959469314)

    return chi*(chi2 + chi2co*(gstar + psi2*(gamma*chi2co + mm*chi2 +
           psi2co*(c00 + c10*chi2 + c01*psi2 + c11*chi2psi2 + c20*chi4 +
           c02*psi4)) + chi2*(omega1 - chi2co*(d0 + d1*chi2))))
end

# Evaluate the CSC inverse (chi,psi → xf) single-axis polynomial.
# p is the primary, q is the secondary; swap for yf.
function _csc_inv_axis(xf::T, psi::T) where {T <: Real}
    # WCSLIB 8.9 cscx2s: maps projected coords back to face coords.
    xx = xf^2
    yy = psi^2

    p00 = T(-0.27292696); p10 = T(-0.07629969); p20 = T(-0.22797056)
    p30 = T( 0.54852384); p40 = T(-0.62930065); p50 = T( 0.25795794)
    p60 = T( 0.02584375)
    p01 = T(-0.02819452); p11 = T(-0.01471565); p21 = T( 0.48051509)
    p31 = T(-1.74114454); p41 = T( 1.71547508); p51 = T(-0.53022337)
    p02 = T( 0.27058160); p12 = T(-0.56800938); p22 = T( 0.30803317)
    p32 = T( 0.98938102); p42 = T(-0.83180469)
    p03 = T(-0.60441560); p13 = T( 1.50880086); p23 = T(-0.93678576)
    p33 = T( 0.08693841)
    p04 = T( 0.93412077); p14 = T(-1.41601920); p24 = T( 0.33887446)
    p05 = T(-0.63915306); p15 = T( 0.52032238)
    p06 = T( 0.14381585)

    z0 = p00 + xx*(p10 + xx*(p20 + xx*(p30 + xx*(p40 + xx*(p50 + xx*p60)))))
    z1 = p01 + xx*(p11 + xx*(p21 + xx*(p31 + xx*(p41 + xx*p51))))
    z2 = p02 + xx*(p12 + xx*(p22 + xx*(p32 + xx*p42)))
    z3 = p03 + xx*(p13 + xx*(p23 + xx*p33))
    z4 = p04 + xx*(p14 + xx*p24)
    z5 = p05 + xx*p15
    z6 = p06

    chi_corr = z0 + yy*(z1 + yy*(z2 + yy*(z3 + yy*(z4 + yy*(z5 + yy*z6)))))
    return xf + xf*(one(T) - xx)*chi_corr
end

"""
    native_to_intermediate(::CSC, phi, theta) -> (x, y)

Forward CSC projection.  Paper II, Section 8.2.

The CSC polynomial coefficients are stored as `float` (32-bit) in WCSLIB,
so results computed here in Float64 differ from Astropy/WCSLIB by up to
~2 µdeg (≈ 7 mas).  Both implementations are correct within their respective
precisions.
"""
function native_to_intermediate(::CSC, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    l, m, n = _sphere_to_uvec(T(phi), T(theta))
    face, xf, yf = _xyz_to_cube_face(l, m, n)
    chi = _csc_fwd_axis(xf, yf)
    psi = _csc_fwd_axis(yf, xf)
    s = T(_FACESCALE)
    x = s * (chi + T(_FACE_X0U[face + 1]))
    y = s * (psi + T(_FACE_Y0U[face + 1]))
    return x, y
end

"""
    intermediate_to_native(::CSC, x, y) -> (phi, theta)

Inverse CSC projection.  Paper II, Section 8.2.
"""
function intermediate_to_native(::CSC, x::Real, y::Real)
    T = _promote_float_type(x, y)
    face, chi, psi = _cube_xy_to_face(T(x), T(y), T(_FACESCALE))
    xf = _csc_inv_axis(chi, psi)
    yf = _csc_inv_axis(psi, chi)
    l, m, n = _face_to_uvec(face, xf, yf)
    phi   = atan(m, l)
    theta = asin(clamp(n, -one(T), one(T)))
    return phi, theta
end

function _qsc_forward(xi::T, eta::T, zeco::T, zeta::T) where {T <: Real}
    # WCSLIB 8.9 qscs2x: raw face components + angular distance → (xf, yf).
    xf = zero(T)
    yf = zero(T)
    if xi != 0 || eta != 0
        if -xi > abs(eta)
            omega = eta / xi
            tau   = one(T) + omega^2
            xf    = -sqrt(zeco / (one(T) - one(T) / sqrt(one(T) + tau)))
            yf    = (xf / 15) * (atand(omega) - asind(omega / sqrt(2 * tau)))
        elseif xi > abs(eta)
            omega = eta / xi
            tau   = one(T) + omega^2
            xf    = sqrt(zeco / (one(T) - one(T) / sqrt(one(T) + tau)))
            yf    = (xf / 15) * (atand(omega) - asind(omega / sqrt(2 * tau)))
        elseif -eta >= abs(xi)
            omega = xi / eta
            tau   = one(T) + omega^2
            yf    = -sqrt(zeco / (one(T) - one(T) / sqrt(one(T) + tau)))
            xf    = (yf / 15) * (atand(omega) - asind(omega / sqrt(2 * tau)))
        elseif eta >= abs(xi)
            omega = xi / eta
            tau   = one(T) + omega^2
            yf    = sqrt(zeco / (one(T) - one(T) / sqrt(one(T) + tau)))
            xf    = (yf / 15) * (atand(omega) - asind(omega / sqrt(2 * tau)))
        end
    end
    return xf, yf
end

# QSC inverse: (chi, psi) ∈ projected coords → (xf, yf).
# From wcslib qscx2s:
#   if |chi| >= |psi|: xf = chi / (pi/4 + (1-pi/4)*(chi^2/psi^2 - 1 + psi^2/chi^2))
#   ... this is the inverse of chi = xf*(pi/4 + (1-pi/4)*yf^2/xf^2)
# Solving: chi = a*(pi/4 + s*b^2/a^2) = a*pi/4 + s*b^2/a
# When |a|>=|b|: chi ≈ a*(pi/4) approximately. Solve by Newton:
function _qsc_inverse(xf::T, yf::T) where {T <: Real}
    # WCSLIB 8.9 qscx2s: face coords → (zeta, w, omega, direct).
    if abs(xf) < T(1e-15) && abs(yf) < T(1e-15)
        return one(T), zero(T), zero(T), true, xf, yf
    end
    SQRT2INV = one(T) / sqrt(T(2))
    direct = abs(xf) > abs(yf)
    if direct
        w = 15 * yf / xf
        omega = sind(w) / (cosd(w) - SQRT2INV)
        tau = one(T) + omega^2
        zeco = xf^2 * (one(T) - one(T) / sqrt(one(T) + tau))
        zeta = one(T) - zeco
    else
        w = 15 * xf / yf
        omega = sind(w) / (cosd(w) - SQRT2INV)
        tau = one(T) + omega^2
        zeco = yf^2 * (one(T) - one(T) / sqrt(one(T) + tau))
        zeta = one(T) - zeco
    end
    return zeta, w, omega, direct, xf, yf
end

"""
    native_to_intermediate(::QSC, phi, theta) -> (x, y)

Forward QSC projection.  Paper II, Section 8.3.
"""
function native_to_intermediate(::QSC, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    l, m, n = _sphere_to_uvec(T(phi), T(theta))
    face, _, _ = _xyz_to_cube_face(l, m, n)
    zeta = max(abs(l), abs(m), abs(n))
    zeco = 1 - zeta
    # Raw component values per wcslib face switch.
    if face == 1
        xi, eta = m, n
    elseif face == 2
        xi, eta = -l, n
    elseif face == 3
        xi, eta = -m, n
    elseif face == 4
        xi, eta = l, n
    elseif face == 5
        xi, eta = m, l
    else
        xi, eta = m, -l
    end
    xf, yf = _qsc_forward(xi, eta, zeco, zeta)
    s = T(_FACESCALE)
    x0v = (0.0, 0.0, 2.0, 4.0, 6.0, 0.0)
    y0v = (2.0, 0.0, 0.0, 0.0, 0.0, -2.0)
    x = s * (xf + x0v[face + 1])
    y = s * (yf + y0v[face + 1])
    return T(x), T(y)
end

"""
    intermediate_to_native(::QSC, x, y) -> (phi, theta)

Inverse QSC projection.  Paper II, Section 8.3.
"""
function intermediate_to_native(::QSC, x::Real, y::Real)
    T = _promote_float_type(x, y)
    face, xf_face, yf_face = _cube_xy_to_face(T(x), T(y), T(_FACESCALE))
    zeta, w, omega, direct, _, _ = _qsc_inverse(xf_face, yf_face)

    if zeta < -one(T)
        zeta = -one(T)
        w = zero(T)
    else
        zeco = 1 - zeta
        w = sqrt(max(zero(T), zeco * (2 - zeco) / (1 + omega^2)))
    end

    if face == 1
        l = zeta
        if direct
            m = w; if xf_face < 0; m = -m; end
            n = m * omega
        else
            n = w; if yf_face < 0; n = -n; end
            m = n * omega
        end
    elseif face == 2
        m = zeta
        if direct
            l = w; if xf_face > 0; l = -l; end
            n = -l * omega
        else
            n = w; if yf_face < 0; n = -n; end
            l = -n * omega
        end
    elseif face == 3
        l = -zeta
        if direct
            m = w; if xf_face > 0; m = -m; end
            n = -m * omega
        else
            n = w; if yf_face < 0; n = -n; end
            m = -n * omega
        end
    elseif face == 4
        m = -zeta
        if direct
            l = w; if xf_face < 0; l = -l; end
            n = l * omega
        else
            n = w; if yf_face < 0; n = -n; end
            l = n * omega
        end
    elseif face == 5
        n = -zeta
        if direct
            m = w; if xf_face < 0; m = -m; end
            l = m * omega
        else
            l = w; if yf_face < 0; l = -l; end
            m = l * omega
        end
    else  # face == 0
        n = zeta
        if direct
            m = w; if xf_face < 0; m = -m; end
            l = -m * omega
        else
            l = w; if yf_face > 0; l = -l; end
            m = -l * omega
        end
    end

    r = hypot(l, m, n)
    l /= r; m /= r; n /= r
    phi = atan(m, l)
    theta = asin(clamp(n, -one(T), one(T)))
    return phi, theta
end

# ──────────────────────────────────────────────────────────────────────────────
# HPX – HEALPix projection   (Calabretta & Roukema 2007, FITS version)
# ──────────────────────────────────────────────────────────────────────────────
#
# WCSLIB 8.9 hpxset/hpxx2s/hpxs2x.  Parameters H = pv[1] (default 4),
# K = pv[2] (default 3).  Derived constants:
#   w[2] = (K-1)/K            sin(theta) equatorial boundary (= 2/3 for K=3)
#   w[3] = 90*K/H             scale: y = w[3]*sin(theta) in equatorial zone
#   w[4] = (K+1)/2            half the number of polar diamonds per cap
#   w[5] = 90*(K-1)/H         y boundary between equatorial and polar (45° for H=4,K=3)
#   w[6] = 180/H              longitude step (= 45° for H=4)
#   w[8] = w[3]*w[0]          = w[3] (since w[0]=1 for R2D)
#   w[9] = w[6]*w[0]          = w[6]
#
# Equatorial: |y| ≤ w[5],  |sin(theta)| ≤ w[2]
#   Forward:  x = phi,  y = w[8]*sin(theta)
#   Inverse:  theta = asind(y/w[3])
#
# Polar:
#   Forward:  sigma = sqrt(K*(1-|sin|)),  y = ±w[9]*(w[4]-sigma)
#             x_c = diamond centre,  x = x_c + (phi - x_c)*sigma
#   Inverse:  sigma = w[4] - |y|/w[6],  sin(theta) = ±(1-sigma²/K)
#             phi = x + (x - x_c)*(1/sigma - 1) with diamond alignment

"""
    native_to_intermediate(proj::HPX, phi, theta) -> (x, y)

Forward HPX (HEALPix) projection.  WCSLIB 8.9 hpxs2x.
"""
function native_to_intermediate(proj::HPX, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    H = T(proj.H)
    K = T(proj.K)

    # WCSLIB w[] array: w[0]=1, w[1]=1 for r0=R2D.
    w2 = (K - 1) / K   # sin(theta) boundary
    w3 = 90 * K / H    # equatorial scale
    w4 = (K + 1) / 2   # half polar diamonds
    w9 = 180 / H       # w[6]*w[0] = 180/H

    sinth = sin(theta)
    abs_sin = abs(sinth)

    if abs_sin <= w2
        # Equatorial regime.  x = phi (w[0]=1, x0=0), y = w[3]*sin(theta).
        return rad2deg(phi), w3 * sinth
    end

    # Polar regime.
    sigma = sqrt(K * (1 - abs_sin))
    m = proj.H % 2    # WCSLIB: ((int)(H+0.5))%2 — 0 for even H
    n = proj.K % 2    # WCSLIB: ((int)(K+0.5))%2 — 1 for odd K
    offset = (n != 0 || sinth > 0) ? 0 : 1

    phi_deg = rad2deg(phi)
    # Diamond centre: x_c = -180 + (2*floor((phi+180)/360*H) + 1) * 180/H
    w7 = H / 360
    phi_c = -180 + (2 * floor((phi_deg + 180) * w7) + 1) * w9
    # phi - phi_c
    dx = phi_deg - phi_c

    xi = sigma - 1                      # distortion factor
    x = phi_deg + dx * xi                    # x = phi + (phi-phi_c)*(sigma-1)
    y = w9 * (w4 - sigma)                    # base y
    if sinth < 0
        y = -y
    end

    # Southern half-facet offset for even K
    if offset != 0
        h = floor(Int, phi_deg / w9) + m
        y += isodd(h) ? -w9 : w9
    end

    # Put phi=180 meridian in expected place
    if x > 180
        x = 360 - x
    end

    return x, y
end

"""
    intermediate_to_native(proj::HPX, x, y) -> (phi, theta)

Inverse HPX (HEALPix) projection.  WCSLIB 8.9 hpxx2s.
"""
function intermediate_to_native(proj::HPX, x::Real, y::Real)
    T = _promote_float_type(x, y)
    H = T(proj.H)
    K = T(proj.K)

    # WCSLIB w[] array (w[0]=1, w[1]=1 for r0=R2D).
    w3 = 90 * K / H  # equatorial scale
    w4 = (K + 1) / 2 # half polar diamonds
    w5 = 90 * (K - 1) / H  # |y| equatorial boundary
    w6 = 180 / H           # longitude step
    w9 = w6                # w[6]*w[0]
    w7 = H / 360

    absy = abs(y)

    if absy <= w5
        # Equatorial regime.
        phi = deg2rad(x)
        theta = asin(clamp(y / w3, -one(T), one(T)))
        return phi, theta
    end

    # Polar regime.
    n = proj.K % 2    # WCSLIB: ((int)(K+0.5))%2 — 1 for odd K
    m = proj.H % 2    # WCSLIB: ((int)(H+0.5))%2 — 0 for even H
    offset = (n != 0 || y > 0) ? 0 : 1

    sigma = w4 - absy / w6
    y_positive = y >= 0

    if sigma == 0
        sin_theta = one(T)
    else
        sin_theta = 1 - sigma^2 / K
        sin_theta = clamp(sin_theta, -one(T), one(T))
    end
    theta = y_positive ? asin(sin_theta) : -asin(sin_theta)

    # Diamond centre: x_c = -180 + (2*floor((x+180)/360*H) + 1) * 180/H
    phi_c = -180 + (2 * floor((x + 180) * w7) + 1) * w6
    dx = x - phi_c                              # x - x_c

    # Southern half-facet offset
    if offset != 0
        h = floor(Int, x / w6) + m
        dx += isodd(h) ? w6 : -w6
    end

    if sigma > 0
        phi = deg2rad(x + dx * (1 / sigma - 1))
    else
        phi = zero(T)
    end

    return phi, theta
end

# ──────────────────────────────────────────────────────────────────────────────
# XPH – HEALPix polar cap projection   (Calabretta & Roukema 2007)
# ──────────────────────────────────────────────────────────────────────────────
#
# WCSLIB 8.9 xphset/xphx2s/xphs2x.  XPH presents the north HEALPix polar cap
# as a rotated diamond covering ±90° in (x, y).  It uses its own scaling
# (w[0] = R2D/√2), a fiducial offset at (0°, 90°), and distinct quadrant logic.
# Derived constants for standard r0 = R2D:
#   w[0] = R2D / √2  ≈ 40.52   w[1] = 1/w[0]
#   w[2] = 2/3                 w[4] = √(2/3) * R2D  ≈ 60.96
#   w[5] = 90 - 1e-4*w[4]     w[6] = √1.5 * D2R    ≈ 0.0215

"""
    native_to_intermediate(::XPH, phi, theta) -> (x, y)

Forward XPH projection.  WCSLIB 8.9 xphs2x.
"""
function native_to_intermediate(::XPH, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)

    w0 = 1 / sqrt(T(2))
    w2 = T(2) / T(3)
    w4 = rad2deg(sqrt(w2))
    w5 = 90 - T(1e-4) * w4
    w6 = deg2rad(sqrt(T(1.5)))

    sinthe = sin(theta)
    abssin = abs(sinthe)

    chi = rad2deg(phi)
    if abs(chi) >= 180
        chi = mod(chi, 360)
        if chi < -180
            chi += 360
        elseif chi >= 180
            chi -= 360
        end
    end
    chi += 180
    psi = mod(chi, 90)
    chi -= 180

    if abssin <= w2
        xi = psi
        eta = T(67.5) * sinthe
    else
        sigma = if theta < w5
            sqrt(3 * (1 - abssin))
        else
            (90 - theta) * w6
        end
        xi  = 45 + (psi - 45) * sigma
        eta = 45 * (2 - sigma)
        if sinthe < 0
            eta = -eta
        end
    end

    xi -= 45
    eta -= 90

    if chi < -90
        x = w0 * (-xi + eta);  y = w0 * (-xi - eta)
    elseif chi < 0
        x = w0 * ( xi + eta);  y = w0 * (-xi + eta)
    elseif chi < 90
        x = w0 * ( xi - eta);  y = w0 * ( xi + eta)
    else
        x = w0 * (-xi - eta);  y = w0 * ( xi - eta)
    end
    return x, y
end

"""
    intermediate_to_native(::XPH, x, y) -> (phi, theta)

Inverse XPH projection.  WCSLIB 8.9 xphx2s.
"""
function intermediate_to_native(::XPH, x::Real, y::Real)
    T = _promote_float_type(x, y)

    w0 = inv(sqrt(T(2)))
    w0i = w0    # w[1] = w[0] = 1/√2 in WCSLIB xphset
    w4 = rad2deg(sqrt(T(2) / T(3)))

    xr = x * w0i
    yr = y * w0i

    if xr <= 0 && yr > 0
        xi1  = -xr - yr;  eta1 =  xr - yr;  phi_base = -T(180)
    elseif xr < 0 && yr <= 0
        xi1  =  xr - yr;  eta1 =  xr + yr;  phi_base = -T(90)
    elseif xr >= 0 && yr < 0
        xi1  =  xr + yr;  eta1 = -xr + yr;  phi_base = T(0)
    else
        xi1  = -xr + yr;  eta1 = -xr - yr;  phi_base = T(90)
    end

    xi = xi1  + 45
    eta = eta1 + 90
    abseta = abs(eta)

    if abseta <= 45
        phi = deg2rad(phi_base + xi)
        theta = asin(clamp(eta / T(67.5), -one(T), one(T)))
    elseif abseta <= 90
        sigma = (90 - abseta) / 45
        if xr == 0
            phi = yr <= 0 ? zero(T) : _pi(T)
        elseif yr == 0
            phi = xr < 0 ? -_halfpi(T) : _halfpi(T)
        else
            phi = deg2rad(phi_base + 45 + xi1 / sigma)
        end
        theta_abs = if sigma < T(1e-4)
            _halfpi(T) - sigma * deg2rad(w4)
        else
            asin(1 - sigma^2 / 3)
        end
        theta = eta >= 0 ? theta_abs : -theta_abs
    else
        error("XPH: point outside valid domain")
    end
    return phi, theta
end

# ──────────────────────────────────────────────────────────────────────────────
# UnknownProjection – catch-all error
# ──────────────────────────────────────────────────────────────────────────────

function intermediate_to_native(proj::UnknownProjection, x::Real, y::Real)
    error("Projection \"$(proj.code)\" is not implemented in FITSWCS.jl")
end

function native_to_intermediate(proj::UnknownProjection, phi::Real, theta::Real)
    error("Projection \"$(proj.code)\" is not implemented in FITSWCS.jl")
end
