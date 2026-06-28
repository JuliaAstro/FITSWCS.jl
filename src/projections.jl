"""
Spherical projection functions.

Each projection pair converts between native spherical coordinates (φ, θ)
expressed in **radians** and intermediate world coordinates (x, y) expressed
in **degrees**.  This matches the FITS WCS Paper II convention.

## References

- Calabretta & Greisen (2002), "Representations of celestial coordinates in
  FITS", Astronomy & Astrophysics, 395, 1077–1122.  (Paper II)
"""

const _R2D = 180.0 / π    # radians → degrees
const _D2R = π / 180.0    # degrees → radians

# ──────────────────────────────────────────────────────────────────────────────
# Shared zenithal utilities
# ──────────────────────────────────────────────────────────────────────────────

"""
Shared azimuth angle φ from intermediate coordinates (x, y) [degrees].
Common to all zenithal projections.  Paper II, Eq. 14.
"""
@inline _phi_zenithal(x::Real, y::Real) = atan(x, -y)  # radians

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
@inline _wrap_native_phi(phi::Real) = phi - 2π * round(phi / (2π))

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
    phi   = _phi_zenithal(x, y)
    Rth   = sqrt(x^2 + y^2)             # degrees
    theta = atan(_R2D, Rth)             # = atan2(180/π, Rth) in radians
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
    if sth <= 0.0
        error("TAN projection: theta must be > 0 (got θ = $(rad2deg(theta))°)")
    end
    Rth = _R2D * cos(theta) / sth       # degrees
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
    xi  = sin_proj.xi
    eta = sin_proj.eta
    # Convert x, y to radians (projection formulas use dimensionless coords)
    xr = x * _D2R
    yr = y * _D2R
    if xi == 0.0 && eta == 0.0
        # Standard SIN: R_θ = cos(θ)  (in unit sphere coords)
        r2 = xr^2 + yr^2
        if r2 > 1.0
            error("SIN projection: point outside valid domain (R_θ = $(sqrt(r2)) > 1)")
        end
        theta = acos(sqrt(r2))   # θ = acos(R_θ)
        phi   = atan(xr, -yr)
    else
        # Slant SIN: solve quadratic (Paper II, Eq. 49)
        a = xi^2 + eta^2 + 1.0
        b = xi*(xr - xi) + eta*(yr - eta)
        c = (xr - xi)^2 + (yr - eta)^2 - 1.0
        disc = b^2 - a*c
        if disc < 0.0
            error("SIN projection: point outside valid domain (discriminant < 0)")
        end
        sth1 = (-b + sqrt(disc)) / a
        sth2 = (-b - sqrt(disc)) / a
        # Choose the solution with θ ≥ θ_0 = 0° (i.e., sin(θ) ≥ 0 preferred)
        # For the standard convention, take the larger sin(θ) value.
        sth = sth1 >= sth2 ? sth1 : sth2
        if abs(sth) > 1.0
            sth = clamp(sth, -1.0, 1.0)
        end
        theta = asin(sth)
        offset = 1.0 - sth
        phi   = atan(xr - xi*offset, -(yr - eta*offset))
    end
    return phi, theta
end

"""
    native_to_intermediate(::SIN, phi, theta) -> (x, y)

Forward SIN projection.
"""
function native_to_intermediate(sin_proj::SIN, phi::Real, theta::Real)
    xi  = sin_proj.xi
    eta = sin_proj.eta
    cth = cos(theta)
    sth = sin(theta)
    # Result in radians, then convert to degrees
    xr = cth * sin(phi) + xi  * (1.0 - sth)
    yr = -cth * cos(phi) + eta * (1.0 - sth)
    return xr * _R2D, yr * _R2D
end

# ──────────────────────────────────────────────────────────────────────────────
# STG – Stereographic projection   (Paper II, Eq. 50)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(::STG, x, y) -> (phi, theta)

Inverse STG projection.
"""
function intermediate_to_native(::STG, x::Real, y::Real)
    phi   = _phi_zenithal(x, y)
    Rth   = sqrt(x^2 + y^2)             # degrees
    # R_θ = 2 * (180/π) * cos(θ) / (1 + sin(θ))
    # → 1 + sin(θ) = 2*(180/π)*cos(θ)/R_θ
    # Solve: let s = sin(θ). Rth/R2D = 2cos(θ)/(1+s) = 2*sqrt(1-s²)/(1+s)
    # = 2*sqrt((1-s)(1+s))/(1+s) = 2*sqrt((1-s)/(1+s))
    # → (Rth/(2*R2D))² = (1-s)/(1+s)
    # → r = Rth/(2*R2D): s = (1-r²)/(1+r²)
    r   = Rth / (2.0 * _R2D)
    sth = (1.0 - r^2) / (1.0 + r^2)
    theta = asin(clamp(sth, -1.0, 1.0))
    return phi, theta
end

"""
    native_to_intermediate(::STG, phi, theta) -> (x, y)

Forward STG projection.
"""
function native_to_intermediate(::STG, phi::Real, theta::Real)
    denom = 1.0 + sin(theta)
    if denom == 0.0
        error("STG projection: singularity at theta = -90°")
    end
    Rth = 2.0 * _R2D * cos(theta) / denom
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
    phi   = _phi_zenithal(x, y)
    Rth   = sqrt(x^2 + y^2)             # degrees
    # R_θ = (180/π)*(π/2 − θ)  →  θ = π/2 − R_θ*(π/180)
    theta = π/2 - Rth * _D2R
    return phi, theta
end

"""
    native_to_intermediate(::ARC, phi, theta) -> (x, y)

Forward ARC projection.
"""
function native_to_intermediate(::ARC, phi::Real, theta::Real)
    Rth = _R2D * (π/2 - theta)          # degrees
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
    phi   = _phi_zenithal(x, y)
    Rth   = sqrt(x^2 + y^2)             # degrees
    # R_θ = 2*(180/π)*sin((π/2 − θ)/2)  →  sin((π/2−θ)/2) = R_θ*π/(360)
    arg   = Rth * _D2R / 2.0            # = R_θ * π/360
    if abs(arg) > 1.0
        error("ZEA projection: point outside valid domain (|arg| = $(abs(arg)) > 1)")
    end
    theta = π/2 - 2.0 * asin(arg)
    return phi, theta
end

"""
    native_to_intermediate(::ZEA, phi, theta) -> (x, y)

Forward ZEA projection.
"""
function native_to_intermediate(::ZEA, phi::Real, theta::Real)
    Rth = 2.0 * _R2D * sin((π/2 - theta) / 2.0)    # degrees
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
    phi   = x * _D2R
    theta = y * _D2R
    return phi, theta
end

"""
    native_to_intermediate(::CAR, phi, theta) -> (x, y)

Forward CAR projection.
"""
function native_to_intermediate(::CAR, phi::Real, theta::Real)
    # Use the local longitude branch so inverse transforms prefer nearby pixels.
    phi_w = _wrap_native_phi(phi)
    x = phi_w * _R2D
    y = theta * _R2D
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
    # Convert the linear longitude coordinate directly to native longitude.
    phi = x * _D2R

    # Recover latitude from the equal-area ordinate and check the finite domain.
    arg = proj.lambda * y * _D2R
    abs(arg) <= 1.0 ||
        error("CEA projection: point outside valid domain (|lambda*y*pi/180| = $(abs(arg)) > 1)")
    theta = asin(clamp(arg, -1.0, 1.0))
    return phi, theta
end

"""
    native_to_intermediate(proj::CEA, phi, theta) -> (x, y)

Forward CEA projection.
"""
function native_to_intermediate(proj::CEA, phi::Real, theta::Real)
    # Longitude is linear, using the local branch around the fiducial meridian.
    phi_w = _wrap_native_phi(phi)
    x = phi_w * _R2D

    # Latitude maps by the equal-area sine relation with lambda scaling.
    y = _R2D * sin(theta) / proj.lambda
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
    # Longitude scales linearly by the cylindrical perspective lambda.
    phi = (x * _D2R) / proj.lambda

    # Solve the perspective ordinate relation in closed form.
    eta = y * _D2R / (proj.mu + proj.lambda)
    theta = atan(eta, 1.0) + asin(clamp((eta * proj.mu) / hypot(eta, 1.0), -1.0, 1.0))
    return phi, theta
end

"""
    native_to_intermediate(proj::CYP, phi, theta) -> (x, y)

Forward CYP projection.
"""
function native_to_intermediate(proj::CYP, phi::Real, theta::Real)
    # Wrap longitude locally before applying the linear cylindrical scale.
    phi_w = _wrap_native_phi(phi)
    x = _R2D * proj.lambda * phi_w

    # Project latitude by the perspective cylinder relation.
    denom = proj.mu + cos(theta)
    denom != 0.0 ||
        error("CYP projection: singularity where mu + cos(theta) = 0")
    y = _R2D * (proj.mu + proj.lambda) * sin(theta) / denom
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
    phi = x * _D2R
    theta = 2.0 * atan(exp(y * _D2R)) - π/2
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
    x = phi_w * _R2D
    y = _R2D * log(tan(π/4 + theta/2))
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
    theta = y * _D2R
    cth = cos(theta)
    abs(cth) > 1e-14 ||
        error("SFL projection: longitude is undefined at theta = ±90°")
    phi = (x * _D2R) / cth
    return phi, theta
end

"""
    native_to_intermediate(::SFL, phi, theta) -> (x, y)

Forward SFL projection.
"""
function native_to_intermediate(::SFL, phi::Real, theta::Real)
    # Longitude contracts by cos(theta), with latitude preserved.
    phi_w = _wrap_native_phi(phi)
    x = _R2D * phi_w * cos(theta)
    y = _R2D * theta
    return x, y
end

# ──────────────────────────────────────────────────────────────────────────────
# PAR – Parabolic projection   (Paper II, Eq. 89)
# ──────────────────────────────────────────────────────────────────────────────

"""
    intermediate_to_native(::PAR, x, y) -> (phi, theta)

Inverse PAR projection.
"""
function intermediate_to_native(::PAR, x::Real, y::Real)
    # Recover latitude from the parabolic sine ordinate.
    arg = (y * _D2R) / π
    abs(arg) <= 1.0 ||
        error("PAR projection: point outside valid domain (|y*pi/180/pi| = $(abs(arg)) > 1)")
    theta = 3.0 * asin(clamp(arg, -1.0, 1.0))

    # Longitude uses the latitude-dependent parabolic scale factor.
    scale = 2.0 * cos(2.0 * theta / 3.0) - 1.0
    scale != 0.0 ||
        error("PAR projection: longitude scale is zero")
    phi = (x * _D2R) / scale
    return phi, theta
end

"""
    native_to_intermediate(::PAR, phi, theta) -> (x, y)

Forward PAR projection.
"""
function native_to_intermediate(::PAR, phi::Real, theta::Real)
    # Apply the parabolic longitude scale and sine latitude ordinate.
    phi_w = _wrap_native_phi(phi)
    x = _R2D * phi_w * (2.0 * cos(2.0 * theta / 3.0) - 1.0)
    y = _R2D * π * sin(theta / 3.0)
    return x, y
end

# ──────────────────────────────────────────────────────────────────────────────
# MOL – Mollweide projection   (Paper II, Eq. 90)
# ──────────────────────────────────────────────────────────────────────────────

const _MOL_MAXITER = 30
const _MOL_TOL = 1e-14

function _mollweide_gamma(theta::Real)
    # Solve 2γ + sin(2γ) = π sin(theta) with Newton iteration.
    abs(abs(theta) - π/2) <= 1e-14 && return copysign(π/2, theta)
    gamma = theta
    target = π * sin(theta)
    for _ in 1:_MOL_MAXITER
        f = 2.0 * gamma + sin(2.0 * gamma) - target
        fp = 2.0 + 2.0 * cos(2.0 * gamma)
        fp != 0.0 || break
        step = f / fp
        gamma -= step
        abs(step) <= _MOL_TOL && return gamma
    end
    error("MOL projection: auxiliary angle solve failed to converge")
end

"""
    intermediate_to_native(::MOL, x, y) -> (phi, theta)

Inverse MOL projection.
"""
function intermediate_to_native(::MOL, x::Real, y::Real)
    # Recover the auxiliary angle gamma from the vertical coordinate.
    sin_gamma = (y * _D2R) / sqrt(2.0)
    abs(sin_gamma) <= 1.0 ||
        error("MOL projection: point outside valid domain (|sin_gamma| = $(abs(sin_gamma)) > 1)")
    gamma = asin(clamp(sin_gamma, -1.0, 1.0))

    # Convert gamma to native latitude and undo the longitude scale.
    theta = asin(clamp((2.0 * gamma + sin(2.0 * gamma)) / π, -1.0, 1.0))
    cos_gamma = cos(gamma)
    abs(cos_gamma) > 1e-14 || return 0.0, theta
    phi = (x * _D2R) * π / (2.0 * sqrt(2.0) * cos_gamma)
    return phi, theta
end

"""
    native_to_intermediate(::MOL, phi, theta) -> (x, y)

Forward MOL projection.
"""
function native_to_intermediate(::MOL, phi::Real, theta::Real)
    # Solve the implicit Mollweide latitude equation before projecting.
    phi_w = _wrap_native_phi(phi)
    gamma = _mollweide_gamma(theta)
    x = _R2D * (2.0 * sqrt(2.0) / π) * phi_w * cos(gamma)
    y = _R2D * sqrt(2.0) * sin(gamma)
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
    u  = x * _D2R
    v  = y * _D2R
    s  = 1.0 - u^2 / 8.0 - v^2 / 2.0
    if s < 0.0
        error("AIT projection: point outside valid domain (discriminant s = $s < 0)")
    end
    z     = sqrt(s)
    # g = sqrt(2/(1+z^2)); sinT = v/g; cosT = sqrt(1-sinT^2)
    g     = sqrt(2.0 / (1.0 + z^2))
    sinT  = v / g
    if abs(sinT) > 1.0
        sinT = clamp(sinT, -1.0, 1.0)
    end
    theta = asin(sinT)
    cosT  = cos(theta)
    # sin(phi/2) = u/(2*g*cosT); cos(phi/2) = z^2/cosT
    if abs(cosT) < 1e-12
        # At theta = ±90°, phi is undefined; return phi = 0
        return 0.0, theta
    end
    sinP2 = u / (2.0 * g * cosT)
    cosP2 = z^2 / cosT
    phi   = 2.0 * atan(sinP2, cosP2)
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
    denom = 1.0 + cos(theta) * cos(phi_w / 2.0)
    if denom <= 0.0
        error("AIT projection: degenerate point (1 + cosθ·cos(φ/2) = $denom)")
    end
    g  = sqrt(2.0 / denom)
    x  = 2.0 * _R2D * g * cos(theta) * sin(phi_w / 2.0)
    y  =       _R2D * g * sin(theta)
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
