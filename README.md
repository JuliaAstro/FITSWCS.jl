# FITSWCS.jl

FITSWCS.jl is a pure-Julia implementation of core FITS World Coordinate System
transforms.  It is intentionally dictionary-first: the main parser works with
FITS-like header dictionaries, while FITSIO.jl and FITSFiles.jl inputs are
adapted through package extensions.

This package is still under active development.  It currently focuses on
correct, maintainable coverage of common image WCS cases rather than full
wcslib parity.

## Quick Start

```julia
using FITSWCS

header = Dict(
    "NAXIS"  => 2,
    "CTYPE1" => "RA---TAN",
    "CTYPE2" => "DEC--TAN",
    "CRPIX1" => 512.0, "CRPIX2" => 512.0,
    "CRVAL1" => 83.8221, "CRVAL2" => -5.3911,
    "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4,
)

wcs = WCS(header)
world = pixel_to_world(wcs, [512.0, 512.0])
pixel = world_to_pixel(wcs, world)
```

Pixel coordinates use the FITS WCS convention: pixel `1` is the center of the
first pixel.  Coordinate vectors are ordered by FITS WCS axis number.  Batch
inputs are `naxis x npoints` matrices, where each column is one coordinate, or
vectors of coordinate vectors.

## WCS.jl-Style Names

For migration experiments, FITSWCS.jl also exports WCS.jl-style helper names:

- `WCS(header)` as the main parser/constructor for FITS-like headers
- `WCSTransform(naxis; kwds...)` for WCS.jl-style programmatic construction
  from core vectors/matrices such as `crpix`, `crval`, `cdelt`, `ctype`,
  `cunit`, `pc`, `cd`, and `crota`
- `pix_to_world(wcs, ...)` as an alias for `pixel_to_world(wcs, ...)`
- `world_to_pix(wcs, ...)` as an alias for `world_to_pixel(wcs, ...)`
- `pix_to_world!(wcs, pixels, worlds)` and
  `world_to_pix!(wcs, worlds, pixels)` for caller-provided output arrays

These aliases keep FITSWCS.jl's FITS 1-based pixel convention.  They do not
implement WCS.jl's full status-returning transform API, intermediate work-array
keywords, or an Astropy-style `origin` argument / shim for zero-based pixel
coordinates.

## Supported Header Keywords

The parser currently supports these image-WCS keyword families:

- axis count: `NAXIS`, `WCSAXES`
- per-axis values: `CTYPEi`, `CUNITi`, `CRPIXi`, `CRVALi`, `CDELTi`
- linear transforms: `PCi_ja`, `CDi_ja`, and legacy `CROTA2`
- alternate WCS suffixes through `WCS(header; alt='A')`
- celestial pole keywords: `LONPOLE`, `LATPOLE`
- projection parameters used by implemented projections, including:
  `PV<lat>_1`, `PV<lat>_2`, `PV<lat>_3`, `PV<lat>_0..30` for `ZPN`, and
  `PV<axis>_0..59` coefficients for `TPV` / `TPD`
- SIP distortion: `A_ORDER`, `B_ORDER`, `A_i_j`, `B_i_j`,
  `AP_ORDER`, `BP_ORDER`, `AP_i_j`, `BP_i_j`
- pre-2012 SCAMP TPV compatibility: `-TAN` celestial CTYPEs with high-index
  `PV` coefficients are interpreted as `-TPV`, and `TPV` / `TPD` take
  precedence over SIP when both are present

Celestial units are normalized to degrees at the public API boundary.  Linear,
spectral, time, and Stokes axes currently remain in the units encoded by their
header linear transform.

## Supported Celestial Coordinate Systems

Celestial axes are identified by their CTYPE prefix:

| System | Longitude | Latitude | Description |
|---|---|---|---|
| Equatorial | `RA` | `DEC` | Right ascension / declination |
| Galactic | `GLON` | `GLAT` | Galactic longitude / latitude |
| Ecliptic | `ELON` | `ELAT` | Ecliptic longitude / latitude |
| Helioecliptic | `HLON` | `HLAT` | Helioecliptic longitude / latitude |
| Supergalactic | `SLON` | `SLAT` | Supergalactic longitude / latitude |
| Helioprojective | `HPLN` | `HPLT` | Solar helioprojective longitude / latitude |

All six systems share the same spherical projection and rotation machinery;
only the fiducial native-pole coordinates differ.

## Supported Projections

All 28 WCSLIB spherical projections are implemented and checked against stored
Astropy / WCSLIB fixtures to sub-microarcsecond precision in the tested
regions, except CSC as noted below.  FITSWCS.jl also implements `TPV` / `TPD`
as TAN plus sequent polynomial distortion.

Zenithal: `AZP`, `SZP`, `TAN`, `SIN` (including slant), `STG`, `ARC`, `ZPN`, `ZEA`, `AIR`
Cylindrical: `CAR`, `CEA`, `CYP`, `MER`
Pseudocylindrical / conventional: `SFL`, `PAR`, `MOL`, `PCO`, `AIT`
Conic: `COP`, `COD`, `COE`, `COO`
Polyconic: `BON`
Quadrilateralized spherical cube: `TSC`, `CSC`¹, `QSC`
HEALPix: `HPX`, `XPH`
Distorted tangent plane: `TPV`, `TPD`

¹ CSC matches to ~9 mas due to WCSLIB storing its polynomial coefficients
  as 32-bit `float` while our implementation computes in 64-bit.

Unknown projection codes throw an informative error at transform time.

## FITS Loader Extensions

When the corresponding package is loaded, `WCS` accepts:

- `FITSIO.FITSHeader`
- `FITSIO.HDU`
- `FITSFiles.Card` vectors
- `FITSFiles.HDU`

FITSFiles.jl is currently a regular package dependency, while FITSIO.jl is a
weak dependency; both loader-specific methods are still isolated in extensions.

## Known Limitations

- **Paper III spectral algorithms**: implemented for all algorithm codes
  (`LOG`, `F2W`, `W2F`, `F2V`, `V2F`, `W2V`, `V2W`, air-wavelength cross-conversions)
  and all S-types (`FREQ`, `AFRQ`, `ENER`, `WAVN`, `VRAD`, `WAVE`, `VOPT`,
  `ZOPT`, `AWAV`, `VELO`, `BETA`).  Rest frequency/wavelength (`RESTFRQ`,
  `RESTWAV`) and reference-frame keywords (`SPECSYS`, `SSYSOBS`, `VELOSYS`,
  `ZSOURCE`, `SSYSSRC`, `MJD-AVG`, `OBSGEO-X/Y/Z`) are parsed and stored for
  downstream frame-correction code, but no frame-correction math is performed
  in the default transforms — matching WCSLIB and astropy.
- **Grism** algorithm codes `GRI`/`GRA` are not yet implemented.
- **Spectral world coordinates for algorithm-coded axes are returned in SI**
  (e.g. Hz, m, m/s) regardless of `CUNIT`; display-unit conversion at the API
  boundary is not yet wired in.
- **Time and Stokes axes**: transform linearly but carry no physical
  interpretation (e.g., `MJDREF`, `DATE-OBS`, `TIMESYS`, polarization state).
- **Full WCS.jl API compatibility**: only a partial compatibility layer exists
  (`WCS`, `pix_to_world`, `world_to_pix`, mutating `!` variants).

Benchmarks for representative paths live under `benchmark/`.
