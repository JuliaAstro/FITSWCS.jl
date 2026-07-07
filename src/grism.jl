"""
Grism coordinate transforms (FITS WCS Paper III, Section 5).

Implements the ``GRI`` (grism in vacuum) and ``GRA`` (grism in air)
algorithm codes.  The axis is linear in a "grism parameter" Gamma that
combines the interference and refraction equations for gratings, prisms,
and grisms into a single world coordinate function.

The forward transform chain (five steps):

    w -> Gamma -> gamma -> lambda -> P -> S

where w is the intermediate world coordinate (linear in Gamma), gamma is
the exit angle, lambda is the wavelength, P is the basic spectral type,
and S is the requested spectral type.

Reference: wcslib spc.c:532-565 (setup), 600-611 (forward), 720-745 (inverse).
"""

# ── Type definition ───────────────────────────────────────────────────────────

"""
    GrismSpec

Parsed grism-axis specification.  The axis is linear in the grism parameter
Gamma = tan(gamma - gamma_r - theta), where gamma is the exit angle and
gamma_r is the reference exit angle at the reference wavelength lambda_r.
"""
struct GrismSpec
    axis::Int               # WCS axis number (1-based)
    crval_si::Float64       # lambda_r in SI (meters)
    cdelt_scale::Float64    # dGamma/dS at the reference point
    # Grism parameters from PVi_0 through PVi_6 (defaults applied)
    G::Float64              # PVi_0: grating ruling density (m^-1, default 0)
    m::Float64              # PVi_1: interference order (default 0)
    alpha_deg::Float64      # PVi_2: angle of incidence (deg, default 0)
    n_r::Float64            # PVi_3: refractive index at lambda_r (default 1)
    dn_r::Float64           # PVi_4: dn/dlambda at lambda_r (/m, default 0)
    epsilon_deg::Float64    # PVi_5: grating tilt angle (deg, default 0)
    theta_deg::Float64      # PVi_6: detector tilt angle (deg, default 0)
    # Pass-through for derived S-type derivative calculation (NaN if absent)
    restfrq::Float64        # RESTFRQ -- rest frequency in Hz
    restwav::Float64        # RESTWAV -- rest wavelength in m
    s_type::Symbol          # S-type for lambda<->S conversion (e.g. :VELO, :FREQ)
    algorithm::Symbol       # :GRI (vacuum) or :GRA (air)
    # Precomputed intermediates (radians where applicable)
    Gamma_r::Float64        # -tan(theta) -- reference grism parameter
    dGammadw::Float64       # dGamma/dw (constant)
    beta_r::Float64         # gamma_r, reference exit angle (radians)
    const_term::Float64     # (n_r - dn_r * lambda_r) * sin(alpha) -- numerator offset
    inv_denom::Float64      # 1 / t = 1 / (G*m/cos(epsilon) - dn_r*sin(alpha))
end

# ── Reference wavelength helper ────────────────────────────────────────────────

"""
    _grism_reference_wavelength(crval_si, s_type, algorithm, restfrq, restwav) -> Float64

Return the reference wavelength lambda_r used by the grism equation, in the
native wavelength type of the grism (vacuum for GRI, air for GRA).  For
wavelength S-types, `crval_si` may need vacuum<->air conversion.  For
frequency/velocity S-types, the S-type value is converted to wavelength.
"""
function _grism_reference_wavelength(crval_si::Float64, s_type::Symbol,
                                     algorithm::Symbol, restfrq::Float64,
                                     restwav::Float64)
    # Compute vacuum wavelength from the S-type reference value.
    lambda_vac = if s_type == :FREQ || s_type == :AFRQ || s_type == :ENER
        nu_r = s_type == :FREQ ? crval_si :
               s_type == :AFRQ ? crval_si / (2π) : crval_si / _H_PLANCK
        _C_LIGHT / nu_r
    elseif s_type == :VELO || s_type == :VRAD || s_type == :BETA
        nu_r = if s_type == :VELO
            _velo_to_freq(crval_si, restfrq)
        elseif s_type == :VRAD
            restfrq * (1 - crval_si / _C_LIGHT)
        else
            restfrq * sqrt((_C_LIGHT - crval_si) / (_C_LIGHT + crval_si))
        end
        _C_LIGHT / nu_r
    elseif s_type == :WAVN
        1.0 / crval_si
    elseif s_type == :AWAV
        _awav_to_wave(crval_si)      # air -> vacuum
    else
        crval_si  # WAVE, VOPT, ZOPT: already vacuum wavelength
    end
    # For GRA, the grism equation needs air wavelength.
    if algorithm == :GRA
        return _wave_to_awav(lambda_vac)
    end
    return lambda_vac
end

# ── Derivative helper ─────────────────────────────────────────────────────────

"""
    _grism_dlambda_ds(lambda_r, s_type, algorithm, restfrq, restwav) -> Float64

Compute dlambda/dS at the reference point, where lambda is the wavelength
variable used by the grism equation (vacuum for GRI, air for GRA) and S is
the requested spectral type (e.g. FREQ, WAVE, AWAV).
"""
function _grism_dlambda_ds(lambda_r::Float64, s_type::Symbol, algorithm::Symbol,
                           restfrq::Float64, restwav::Float64)
    # dlambda_vac / dS at reference, for the common S-types found on grism axes.
    nu_r = _C_LIGHT / lambda_r
    dlam_dnu = -_C_LIGHT / (nu_r * nu_r)  # dlambda/dnu (reused for several types)
    dlam_ds = if s_type == :FREQ
        dlam_dnu                         # dlambda_vac / dnu
    elseif s_type == :WAVE
        1.0                              # identity
    elseif s_type == :AWAV
        _dwave_dawav(lambda_r)           # dlambda_vac / dlambda_air
    elseif s_type == :VELO
        # dlambda/dv = dlambda/dnu * dnu/dv (relativistic).
        # At reference, v=0 for grism axes, dnu/dv = -c*nu0 / (c * c) ...
        # Use the known reference-point derivative from the Paper III table.
        dlam_dnu * (-_C_LIGHT * restfrq / (_C_LIGHT * _C_LIGHT))
    elseif s_type == :VRAD
        dlam_dnu * (-restfrq / _C_LIGHT)
    elseif s_type == :VOPT
        restwav / _C_LIGHT               # dlambda/dVOPT
    elseif s_type == :ZOPT
        restwav                           # dlambda/dZOPT
    elseif s_type == :BETA
        dlam_dnu * inv(_C_LIGHT)
    elseif s_type == :AFRQ
        dlam_dnu * 2π                    # dlambda/domega
    elseif s_type == :ENER
        dlam_dnu * inv(_H_PLANCK)
    elseif s_type == :WAVN
        1.0
    else
        1.0
    end
    # For GRA, lambda_grism is air wavelength, convert.
    if algorithm == :GRA
        return dlam_ds / _dwave_dawav(lambda_r)
    end
    return dlam_ds
end

# ── Setup ─────────────────────────────────────────────────────────────────────

"""
    _grism_setup(axis, lambda_r, pvs, dlambda_ds, restfrq, restwav, s_type, algorithm) -> GrismSpec

Precompute grism intermediates from PV parameters and the dlambda/dS derivative
chain.  `pvs` is a 7-element collection of PV values indexed 0..6.  Missing
entries should be `NaN` to trigger the default value.
"""
function _grism_setup(axis::Int, lambda_r::Float64, pvs, dlambda_ds::Float64,
                      restfrq::Float64, restwav::Float64, s_type::Symbol,
                      algorithm::Symbol)
    # Apply defaults for missing PV parameters (NaN signals "not present").
    G       = isfinite(pvs[1]) ? pvs[1] : 0.0      # PVi_0
    m       = isfinite(pvs[2]) ? pvs[2] : 0.0      # PVi_1
    alpha   = isfinite(pvs[3]) ? pvs[3] : 0.0      # PVi_2 (deg)
    n_r     = isfinite(pvs[4]) ? pvs[4] : 1.0      # PVi_3
    dn_r    = isfinite(pvs[5]) ? pvs[5] : 0.0      # PVi_4 (/m)
    epsilon = isfinite(pvs[6]) ? pvs[6] : 0.0      # PVi_5 (deg)
    theta   = isfinite(pvs[7]) ? pvs[7] : 0.0      # PVi_6 (deg)

    # Convert angles to radians for computation.
    alpha_rad = deg2rad(alpha)
    eps_rad   = deg2rad(epsilon)
    theta_rad = deg2rad(theta)

    # Compute intermediate quantities (matching wcslib spcset).
    # t = G*m / cos(epsilon)
    t = G * m / cos(eps_rad)
    # gamma_r = beta_r = asin(t * lambda_r - n_r * sin(alpha))
    sin_beta_r = t * lambda_r - n_r * sin(alpha_rad)
    # Clamp to [-1, 1] to avoid numerical edge cases.
    sin_beta_r = clamp(sin_beta_r, -1.0, 1.0)
    beta_r = asin(sin_beta_r)

    # Denominator of the grism equation: t_denom = t - dn_r * sin(alpha)
    t_denom = t - dn_r * sin(alpha_rad)

    # Precomputed intermediates (matching wcslib w[] array).
    Gamma_r   = -tan(theta_rad)                       # wcslib w[1]
    cos_beta  = cos(beta_r)
    cos_theta_sq = cos(theta_rad) * cos(theta_rad)
    dGammadw  = t_denom / (cos_beta * cos_theta_sq) * dlambda_ds  # wcslib w[2]
    const_term = (n_r - dn_r * lambda_r) * sin(alpha_rad)         # wcslib w[4]
    inv_denom  = 1.0 / t_denom                                     # wcslib w[5]

    return GrismSpec(axis, lambda_r, dGammadw,
                     G, m, alpha, n_r, dn_r, epsilon, theta,
                     restfrq, restwav, s_type, algorithm,
                     Gamma_r, dGammadw, beta_r, const_term, inv_denom)
end

# ── Forward transform: intermediate -> world ──────────────────────────────────

"""
    _grism_x_to_world(w, spec::GrismSpec) -> Float64

Convert an intermediate coordinate offset w (grism parameter Gamma minus Gamma_r)
to a world coordinate (S-type SI) for a grism axis.  The five-step chain:

1. Gamma = Gamma_r + w * dGamma/dw
2. gamma = atan(Gamma) + beta_r + theta
3. lambda = (sin(gamma) + const_term) * inv_denom
4. P = lambda  (X-type = P-type = wavelength; vacuum/air handled by caller)
5. S = P       (the caller applies P->S conversion via the spectral tables)
"""
function _grism_x_to_world(w, spec::GrismSpec)
    T = _float_type(typeof(w))
    # Step 1: w is already the grism parameter offset (dGamma/dS is baked into
    # the CD matrix via cdelt_scale).  The absolute grism parameter is Gamma_r + w.
    Gamma = T(spec.Gamma_r) + w
    # Step 2: gamma from Gamma.
    gamma = atan(Gamma) + T(spec.beta_r) + deg2rad(T(spec.theta_deg))
    # Step 3: lambda from gamma (grism equation).  This is vacuum wavelength for
    # GRI, air wavelength for GRA.
    lambda = (sin(gamma) + T(spec.const_term)) * T(spec.inv_denom)
    # Step 4: convert lambda to the requested S-type (e.g. VELO, FREQ).
    return _grism_lambda_to_s(lambda, spec)
end

# ── Lambda <-> S-type conversion helpers ──────────────────────────────────────

"""
    _grism_lambda_to_s(lambda, spec::GrismSpec) -> Float64

Convert the grism equation output (vacuum wavelength for GRI, air wavelength
for GRA) to the requested S-type world coordinate.
"""
function _grism_lambda_to_s(lambda, spec::GrismSpec)
    T = _float_type(typeof(lambda))
    st = spec.s_type
    # lambda from the grism equation: vacuum for GRI, air for GRA.
    # Convert to vacuum wavelength for non-wavelength S-types.
    lambda_vac = if spec.algorithm == :GRA
        st == :AWAV ? lambda : _awav_to_wave(lambda)
    else
        st == :AWAV ? _wave_to_awav(lambda) : lambda
    end
    if st == :WAVE
        return lambda_vac  # vacuum wavelength
    elseif st == :AWAV
        return spec.algorithm == :GRA ? lambda : _wave_to_awav(lambda)
    elseif st == :FREQ || st == :AFRQ || st == :ENER
        nu = T(_C_LIGHT) / lambda_vac
        if st == :FREQ;     return nu
        elseif st == :AFRQ; return nu * T(2π)
        else;               return nu * T(_H_PLANCK)
        end
    elseif st == :VELO
        nu = T(_C_LIGHT) / lambda_vac
        return _freq_to_velo(nu, T(spec.restfrq))
    elseif st == :VRAD
        nu = T(_C_LIGHT) / lambda_vac
        return T(_C_LIGHT) * (T(spec.restfrq) - nu) / T(spec.restfrq)
    elseif st == :VOPT
        return T(_C_LIGHT) * (lambda_vac - T(spec.restwav)) / T(spec.restwav)
    elseif st == :ZOPT
        return (lambda_vac - T(spec.restwav)) / T(spec.restwav)
    elseif st == :BETA
        nu = T(_C_LIGHT) / lambda_vac
        return _freq_to_velo(nu, T(spec.restfrq)) / T(_C_LIGHT)
    elseif st == :WAVN
        return 1 / lambda_vac
    else
        return lambda  # fallback
    end
end

"""
    _grism_s_to_lambda(world_si, spec::GrismSpec) -> Float64

Convert an S-type world coordinate to the grism equation's native wavelength.
This is the inverse of `_grism_lambda_to_s`.
"""
function _grism_s_to_lambda(world_si, spec::GrismSpec)
    T = _float_type(typeof(world_si))
    st = spec.s_type
    # Convert S-type to vacuum wavelength for the grism equation.
    lambda_vac = if st == :WAVE
        world_si  # already vacuum wavelength
    elseif st == :AWAV
        # S-type is air wavelength; for GRI we need vacuum.
        spec.algorithm == :GRI ? _awav_to_wave(world_si) : world_si
    elseif st == :FREQ
        T(_C_LIGHT) / world_si
    elseif st == :AFRQ
        T(_C_LIGHT) / (world_si / T(2π))
    elseif st == :ENER
        T(_C_LIGHT) / (world_si / T(_H_PLANCK))
    elseif st == :VELO
        nu = _velo_to_freq(world_si, T(spec.restfrq))
        T(_C_LIGHT) / nu
    elseif st == :VRAD
        nu = T(spec.restfrq) * (1 - world_si / T(_C_LIGHT))
        T(_C_LIGHT) / nu
    elseif st == :VOPT
        T(spec.restwav) * (1 + world_si / T(_C_LIGHT))
    elseif st == :ZOPT
        T(spec.restwav) * (1 + world_si)
    elseif st == :BETA
        nu = _velo_to_freq(world_si * T(_C_LIGHT), T(spec.restfrq))
        T(_C_LIGHT) / nu
    elseif st == :WAVN
        1 / world_si
    else
        world_si  # fallback
    end
    # For GRA + AWAV: the wavelength IS the air wavelength (already returned above).
    # For GRA + non-AWAV: convert vacuum -> air for the grism equation.
    if spec.algorithm == :GRA && st != :AWAV
        return _wave_to_awav(lambda_vac)
    end
    return lambda_vac
end

# ── Inverse transform: world -> intermediate ──────────────────────────────────

"""
    _grism_world_to_x(world_si, spec::GrismSpec) -> Float64

Convert a world coordinate (S-type SI, i.e. wavelength in meters) to an
intermediate coordinate offset w for a grism axis.  The inverse chain:

1. lambda = world (S = P = wavelength for grism)
2. gamma = asin(lambda / inv_denom - const_term)
3. Gamma = tan(gamma - beta_r - theta)
4. w = (Gamma - Gamma_r) / dGamma/dw

Returns `NaN` if the world coordinate maps outside the valid grism domain
(|sin(gamma)| > 1).
"""
function _grism_world_to_x(world_si, spec::GrismSpec)
    T = _float_type(typeof(world_si))
    # Step 1: convert S-type to wavelength.
    lambda = _grism_s_to_lambda(T(world_si), spec)
    # Step 2: lambda -> sin(gamma) -> gamma.
    s = lambda / T(spec.inv_denom) - T(spec.const_term)
    abs(s) <= 1 || return T(NaN)
    gamma = asin(s)
    # Step 3: gamma -> Gamma.
    Gamma = tan(gamma - T(spec.beta_r) - deg2rad(T(spec.theta_deg)))
    # Step 4: Gamma -> w (w is the offset from Gamma_r; dGamma/dS is handled
    # by the CD matrix inverse, so the intermediate is just Gamma - Gamma_r).
    return Gamma - T(spec.Gamma_r)
end

# ── Grism WCS data wrapper ────────────────────────────────────────────────────

"""Abstract supertype for grism WCS payloads."""
abstract type AbstractGrismWCSData end

"""No-op grism payload for WCS transforms with no grism axes."""
struct NoGrismWCSData <: AbstractGrismWCSData end

"""Resolved collection of all grism-axis specifications for one WCS."""
struct GrismWCSData{T <: Tuple} <: AbstractGrismWCSData
    specs::T
end
