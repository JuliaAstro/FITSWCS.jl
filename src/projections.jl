"""
Spherical projection functions.

Each projection pair converts between native spherical coordinates (φ, θ)
expressed in **radians** and intermediate world coordinates (x, y) expressed
in **degrees**.  This matches the FITS WCS Paper II convention.

## References

- Calabretta & Greisen (2002), "Representations of celestial coordinates in
  FITS", Astronomy & Astrophysics, 395, 1077–1122.  (Paper II)
"""

# @inline _float_type(::Type{T}) where {T<:AbstractFloat} = T
# @inline _float_type(::Type{T}) where {T<:Real} = Float64
# @inline _promote_float_type(x::Real) = _float_type(typeof(x))
@inline _promote_float_type(x::Real) = float(typeof(x))
@inline _promote_float_type(x::Real, y::Real) =
    promote_type(_promote_float_type(x), _promote_float_type(y))
@inline _halfpi(::Type{T}) where {T<:AbstractFloat} = T(π / 2)
@inline _pi(::Type{T}) where {T<:AbstractFloat} = T(π)

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
    # Longitude scales linearly by the cylindrical perspective lambda.
    phi = deg2rad(x) / T(proj.lambda)

    # Solve the perspective ordinate relation in closed form.
    eta = deg2rad(y) / (T(proj.mu) + T(proj.lambda))
    theta = atan(eta, one(T)) + asin(clamp((eta * T(proj.mu)) / hypot(eta, one(T)), -one(T), one(T)))
    return phi, theta
end

"""
    native_to_intermediate(proj::CYP, phi, theta) -> (x, y)

Forward CYP projection.
"""
function native_to_intermediate(proj::CYP, phi::Real, theta::Real)
    T = _promote_float_type(phi, theta)
    # Wrap longitude locally before applying the linear cylindrical scale.
    phi_w = _wrap_native_phi(phi)
    x = rad2deg(T(proj.lambda) * phi_w)

    # Project latitude by the perspective cylinder relation.
    denom = T(proj.mu) + cos(theta)
    denom != 0.0 ||
        error("CYP projection: singularity where mu + cos(theta) = 0")
    y = rad2deg((T(proj.mu) + T(proj.lambda)) * sin(theta) / denom)
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
# UnknownProjection – catch-all error
# ──────────────────────────────────────────────────────────────────────────────

function intermediate_to_native(proj::UnknownProjection, x::Real, y::Real)
    error("Projection \"$(proj.code)\" is not implemented in FITSWCS.jl")
end

function native_to_intermediate(proj::UnknownProjection, phi::Real, theta::Real)
    error("Projection \"$(proj.code)\" is not implemented in FITSWCS.jl")
end
