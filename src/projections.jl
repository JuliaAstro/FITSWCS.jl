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
    x = phi   * _R2D
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
    # Longitude is linear for cylindrical projections.
    x = phi * _R2D

    # Latitude maps by the equal-area sine relation with lambda scaling.
    y = _R2D * sin(theta) / proj.lambda
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
    phi_w = phi - 2π * round(phi / (2π))
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
