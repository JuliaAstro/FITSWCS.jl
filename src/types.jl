"""
    AbstractProjection

Supertype for all FITS WCS spherical projections.  Each concrete subtype
represents one projection algorithm as defined in Calabretta & Greisen (2002),
Paper II.
"""
abstract type AbstractProjection end

# ── Zenithal (azimuthal) projections ──────────────────────────────────────────

"""    AZP

Zenithal perspective projection.  FITS projection code `AZP`.

Parameters: `mu` (PV<lat>_1, distance from sphere center, default 0),
`gamma` (PV<lat>_2, tilt angle in degrees, default 0).
When both are zero the projection degenerates to TAN (gnomonic).

Paper II, Eq. 25–27.
"""
struct AZP <: AbstractProjection
    mu::Float64
    gamma::Float64
end
AZP() = AZP(0.0, 0.0)

"""    SZP

Slant zenithal perspective projection.  FITS projection code `SZP`.

Parameters: `mu` (PV<lat>_1, distance from sphere center, default 0),
`phi_c` (PV<lat>_2, native longitude of slant, default 0),
`theta_c` (PV<lat>_3, native latitude of slant, default 90).

Paper II, Eq. 31–33.
"""
struct SZP <: AbstractProjection
    mu::Float64
    phi_c::Float64
    theta_c::Float64
end
SZP() = SZP(0.0, 0.0, 90.0)

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

"""    PCO

Polyconic projection.  FITS projection code `PCO`.
"""
struct PCO <: AbstractProjection end

"""    AIT

Hammer-Aitoff projection.  FITS projection code `AIT`.
Paper II, Eq. 75.
"""
struct AIT <: AbstractProjection end

# ── Zenithal polynomial projections ───────────────────────────────────────────

"""    ZPN

Zenithal/azimuthal polynomial projection.  FITS projection code `ZPN`.

The native radius is a polynomial in the colatitude `zd = π/2 − θ` (radians):

    r = Σ pv[m] * zd^m   (m = 0, 1, ..., N)

The polynomial coefficients are supplied through the PV header keywords
`PVlat_0`, `PVlat_1`, ... on the latitude-like axis.  The default
(coefficients all zero) projects to the origin for all sky directions.

Paper II, Eq. 55.
"""
struct ZPN <: AbstractProjection
    pv::Vector{Float64}   # polynomial coefficients pv[1] ≡ PV_0, pv[2] ≡ PV_1, ...
end
ZPN() = ZPN([0.0, 1.0])  # linear r = zd default (same as ARC)

"""    AIR

Airy projection.  FITS projection code `AIR`.

`theta_b` is the break latitude (degrees) at which the projection is zero;
default is 90°.

Paper II, Section 5.5, Eq. 30–31.
"""
struct AIR <: AbstractProjection
    theta_b::Float64   # break latitude in degrees
end
AIR() = AIR(90.0)

# ── Conic projections ──────────────────────────────────────────────────────────

"""    COP

Conic perspective projection.  FITS projection code `COP`.

Parameters `sigma` (PVlat_1) and `delta` (PVlat_2) define the cone aperture.
Both are in degrees.

Paper II, Section 6.1.
"""
struct COP <: AbstractProjection
    sigma::Float64   # degrees; native standard parallel
    delta::Float64   # degrees; half-opening angle
end

"""    COD

Conic equidistant projection.  FITS projection code `COD`.

Parameters `sigma` (PVlat_1) and `delta` (PVlat_2) in degrees.

Paper II, Section 6.2.
"""
struct COD <: AbstractProjection
    sigma::Float64
    delta::Float64
end

"""    COE

Conic equal-area projection.  FITS projection code `COE`.

Parameters `sigma` (PVlat_1) and `delta` (PVlat_2) in degrees.

Paper II, Section 6.3.
"""
struct COE <: AbstractProjection
    sigma::Float64
    delta::Float64
end

"""    COO

Conic orthomorphic projection.  FITS projection code `COO`.

Parameters `sigma` (PVlat_1) and `delta` (PVlat_2) in degrees.

Paper II, Section 6.4.
"""
struct COO <: AbstractProjection
    sigma::Float64
    delta::Float64
end

"""    BON

Bonne's projection.  FITS projection code `BON`.

`theta1` (PVlat_1) is the standard parallel in degrees.
When `theta1 == 0`, the projection degenerates to SFL (Sanson-Flamsteed).

Paper II, Section 7.4, Eq. 70.
"""
struct BON <: AbstractProjection
    theta1::Float64   # degrees
end

# ── Quadrilateralized spherical cube projections ───────────────────────────────

"""    TSC

Tangential spherical cube projection.  FITS projection code `TSC`.

Paper II, Section 8.1.
"""
struct TSC <: AbstractProjection end

"""    CSC

COBE quadrilateralized spherical cube projection.  FITS projection code `CSC`.

Paper II, Section 8.2.
"""
struct CSC <: AbstractProjection end

"""    QSC

Quadrilateralized spherical cube projection.  FITS projection code `QSC`.

Paper II, Section 8.3.
"""
struct QSC <: AbstractProjection end

# ── HEALPix projections ────────────────────────────────────────────────────────

"""    HPX

HEALPix projection.  FITS projection code `HPX`.

Parameters `H` (PVlat_1, default 4) and `K` (PVlat_2, default 3) define the
partition of the sphere.

Calabretta & Roukema (2007).
"""
struct HPX <: AbstractProjection
    H::Int   # number of facets around the equatorial zone
    K::Int   # number of facets in each polar cap
end
HPX() = HPX(4, 3)

"""    XPH

HEALPix polar cap projection (rotated HPX).  FITS projection code `XPH`.

Calabretta & Roukema (2007).
"""
struct XPH <: AbstractProjection end

# ── TPV / TPD ──────────────────────────────────────────────────────────────────

"""    TPV

TAN projection with sequent polynomial distortion (SCAMP convention).

TPV is not a fundamental projection; it is a TAN (gnomonic) projection combined
with a sequent polynomial distortion whose coefficients are stored in FITS
`PVi_m` keywords.  The polynomial operates on intermediate world coordinates
(x, y) in degrees **after** the CD matrix, using the TPD coefficient indexing
convention.

Coefficients are stored in **direct** form (wcslib convention): the polynomial
returns the corrected coordinate directly.  The identity polynomial is
``xcoeff = [0, 1]`` (i.e. ``x' = x``) and ``ycoeff = [0, 0, 1]`` (``y' = y``).

The CTYPE projection codes `TPV` and `TPD` both map to this type.

# Fields
- `xcoeff` – ``PV_{lon\\_axis,m}`` coefficients (direct form)
- `ycoeff` – ``PV_{lat\\_axis,m}`` coefficients (direct form)
"""
struct TPV <: AbstractProjection
    xcoeff::Vector{Float64}
    ycoeff::Vector{Float64}
end

# Default is identity polynomial (= plain TAN).  x' = 1·x, y' = 1·y.
TPV() = TPV(Float64[0.0, 1.0], Float64[0.0, 0.0, 1.0])

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
    crpix::SVector{2,Float64}
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
struct WCSTransform{N,L,P<:Union{Nothing,AbstractProjection},S<:Union{Nothing,SIPDistortion}}
    naxis::Int
    crpix::SVector{N,Float64}
    crval::SVector{N,Float64}
    cd::SMatrix{N,N,Float64,L}       # naxis × naxis
    ctype::Vector{String}
    cunit::Vector{String}
    lonpole::Float64          # degrees; φₚ (Paper II)
    latpole::Float64          # degrees; used during construction
    alpha_p::Float64          # degrees; celestial lon of native N pole
    delta_p::Float64          # degrees; celestial lat of native N pole
    projection::P
    sip::S
    lon_axis::Int             # 1-based index; 0 = no celestial lon axis
    lat_axis::Int             # 1-based index; 0 = no celestial lat axis
end
