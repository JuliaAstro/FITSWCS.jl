"""
    AbstractProjection

Supertype for all FITS WCS spherical projections.  Each concrete subtype
represents one projection algorithm as defined in Calabretta & Greisen (2002),
Paper II.
"""
abstract type AbstractProjection end

# ── Zenithal (azimuthal) projections ──────────────────────────────────────────

"""    TAN

Gnomonic (tangent-plane) projection.  FITS projection code `TAN`.
Paper II, Eq. 54–55.
"""
struct TAN <: AbstractProjection end

"""    SIN

Slant orthographic projection.  FITS projection code `SIN`.
Paper II, Eq. 48.  The standard (non-slant) form has projection parameters
ξ = η = 0.
"""
struct SIN <: AbstractProjection
    xi::Float64    # PV2_1; default 0
    eta::Float64   # PV2_2; default 0
end
SIN() = SIN(0.0, 0.0)

"""    STG

Stereographic projection.  FITS projection code `STG`.
Paper II, Eq. 50.
"""
struct STG <: AbstractProjection end

"""    ARC

Zenithal equidistant projection.  FITS projection code `ARC`.
Paper II, Eq. 46.
"""
struct ARC <: AbstractProjection end

"""    ZEA

Lambert zenithal equal-area projection.  FITS projection code `ZEA`.
Paper II, Eq. 52.
"""
struct ZEA <: AbstractProjection end

# ── Cylindrical projections ───────────────────────────────────────────────────

"""    CAR

Plate carrée (equirectangular) projection.  FITS projection code `CAR`.
Paper II, Eq. 84.  Native longitude/latitude map linearly to x/y.
"""
struct CAR <: AbstractProjection end

"""    CEA

Cylindrical equal-area projection.  FITS projection code `CEA`.
The `lambda` parameter is read from `PV<lat>_1` and defaults to 1.
"""
struct CEA <: AbstractProjection
    lambda::Float64
end
CEA() = CEA(1.0)

"""    CYP

Cylindrical perspective projection.  FITS projection code `CYP`.
Parameters `lambda` and `mu` are read from `PV<lat>_1` and `PV<lat>_2`.
"""
struct CYP <: AbstractProjection
    lambda::Float64
    mu::Float64
end
CYP() = CYP(1.0, 1.0)

"""    MER

Mercator projection.  FITS projection code `MER`.
"""
struct MER <: AbstractProjection end

# ── Pseudo-cylindrical projections ────────────────────────────────────────────

"""    SFL

Sanson-Flamsteed projection.  FITS projection code `SFL`.
"""
struct SFL <: AbstractProjection end

"""    PAR

Parabolic projection.  FITS projection code `PAR`.
"""
struct PAR <: AbstractProjection end

"""    MOL

Mollweide projection.  FITS projection code `MOL`.
"""
struct MOL <: AbstractProjection end

"""    AIT

Hammer-Aitoff projection.  FITS projection code `AIT`.
Paper II, Eq. 75.
"""
struct AIT <: AbstractProjection end

# ── Unknown / deferred ────────────────────────────────────────────────────────

"""    UnknownProjection

Placeholder for projection codes that are parsed from a FITS header but not
yet implemented.  Attempting a coordinate transform with this projection raises
an informative error.
"""
struct UnknownProjection <: AbstractProjection
    code::String
end

# ──────────────────────────────────────────────────────────────────────────────

"""
    SIPDistortion

Simple Imaging Polynomial distortion model.

SIP applies only to the first two pixel axes.  The forward coefficients `a`
and `b` map detector pixel coordinates to focal/image-plane pixel coordinates.
The optional inverse coefficients `ap` and `bp` map focal/image-plane
coordinates back to detector pixel coordinates.
"""
struct SIPDistortion
    crpix::Vector{Float64}
    a::Matrix{Float64}
    b::Matrix{Float64}
    ap::Union{Nothing, Matrix{Float64}}
    bp::Union{Nothing, Matrix{Float64}}
end

# ──────────────────────────────────────────────────────────────────────────────

"""
    WCSTransform

Parsed, validated FITS World Coordinate System transform.

## Coordinate conventions

Pixel coordinates follow the **FITS 1-based** convention: pixel 1 is the
centre of the first array element.  This matches the FITS standard and the
values stored in `CRPIX` header keywords.  When working with Julia 1-based
array indices the values are numerically identical; no offset is required.

## Fields

- `naxis`      – number of WCS axes.
- `crpix`      – reference pixel position (FITS 1-based), length `naxis`.
- `crval`      – world coordinate value at the reference pixel, length `naxis`.
- `cd`         – combined CD matrix (naxis × naxis), where
                 `cd[i,j] = CDELT_i * PC_i_j` (or the explicit `CD_i_j` value).
                 Units match `cunit`.
- `ctype`      – FITS `CTYPEi` strings, length `naxis`.
- `cunit`      – FITS `CUNITi` strings (empty string means degrees for
                 celestial axes), length `naxis`.
- `lonpole`    – native longitude of the celestial pole (degrees), φₚ in Paper II.
- `latpole`    – native latitude of the celestial pole hint (degrees), θₚ in
                 Paper II; used to resolve the `delta_p` ambiguity.
- `alpha_p`    – celestial longitude of the native north pole (degrees), αₚ in
                 Paper II.  Precomputed during construction; not a direct FITS
                 keyword.
- `delta_p`    – celestial latitude of the native north pole (degrees), δₚ in
                 Paper II.  Precomputed during construction.
- `projection` – spherical projection for the celestial axes, or `nothing` for
                 purely linear WCS.
- `sip`        – optional SIP distortion applied before the linear transform.
- `lon_axis`   – 1-based index of the longitude axis; 0 if no celestial axes.
- `lat_axis`   – 1-based index of the latitude axis; 0 if no celestial axes.
"""
struct WCSTransform
    naxis::Int
    crpix::Vector{Float64}
    crval::Vector{Float64}
    cd::Matrix{Float64}       # naxis × naxis
    ctype::Vector{String}
    cunit::Vector{String}
    lonpole::Float64          # degrees; φₚ (Paper II)
    latpole::Float64          # degrees; used during construction
    alpha_p::Float64          # degrees; celestial lon of native N pole
    delta_p::Float64          # degrees; celestial lat of native N pole
    projection::Union{Nothing, AbstractProjection}
    sip::Union{Nothing, SIPDistortion}
    lon_axis::Int             # 1-based index; 0 = no celestial lon axis
    lat_axis::Int             # 1-based index; 0 = no celestial lat axis
end
