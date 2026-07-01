# FITSWCS.jl

FITSWCS.jl is a pure-Julia implementation of core FITS World Coordinate System
transforms.  It is intentionally dictionary-first and FITS-loader agnostic:
FITSIO.jl and FITSFiles.jl support is provided through package extensions.

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

wcs = from_header(header)
world = pixel_to_world(wcs, [512.0, 512.0])
pixel = world_to_pixel(wcs, world)
```

Pixel coordinates use the FITS WCS convention: pixel `1` is the center of the
first pixel.  Coordinate vectors are ordered by FITS WCS axis number.  Batch
inputs are `naxis x npoints` matrices, where each column is one coordinate, or
vectors of coordinate vectors.

## WCS.jl-Style Names

For migration experiments, FITSWCS.jl also exports:

- `WCS(header)` as an alias for `from_header(header)`
- `WCSTransform(naxis; kwds...)` for WCS.jl-style programmatic construction
  from core vectors/matrices such as `crpix`, `crval`, `cdelt`, `ctype`,
  `cunit`, `pc`, `cd`, and `crota`
- `pix_to_world(wcs, ...)` as an alias for `pixel_to_world(wcs, ...)`
- `world_to_pix(wcs, ...)` as an alias for `world_to_pixel(wcs, ...)`
- `pix_to_world!(wcs, pixels, worlds)` and
  `world_to_pix!(wcs, worlds, pixels)` for caller-provided output arrays

These aliases keep FITSWCS.jl's FITS 1-based pixel convention.  They do not
implement WCS.jl's full status-returning transform API, intermediate work-array
keywords, or 0-origin pixel convention.

## Supported Header Keywords

The parser currently supports these image-WCS keyword families:

- axis count: `NAXIS`, `WCSAXES`
- per-axis values: `CTYPEi`, `CUNITi`, `CRPIXi`, `CRVALi`, `CDELTi`
- linear transforms: `PCi_ja`, `CDi_ja`, and legacy `CROTA2`
- alternate WCS suffixes through `from_header(header; alt='A')`
- celestial pole keywords: `LONPOLE`, `LATPOLE`
- projection parameters used by implemented projections:
  `PV<lat>_1` and `PV<lat>_2`
- SIP distortion: `A_ORDER`, `B_ORDER`, `A_i_j`, `B_i_j`,
  `AP_ORDER`, `BP_ORDER`, `AP_i_j`, `BP_i_j`

Celestial units are normalized to degrees at the public API boundary.  Linear,
spectral, time, and Stokes axes currently remain in the units encoded by their
header linear transform.

## Supported Projections

All 28 WCSLIB spherical projections are implemented and verified against
Astropy / WCSLIB to sub-microarcsecond precision (except CSC, see below).

Zenithal: `AZP`, `SZP`, `TAN`, `SIN` (including slant), `STG`, `ARC`, `ZPN`, `ZEA`, `AIR`
Cylindrical: `CAR`, `CEA`, `CYP`, `MER`
Pseudocylindrical / conventional: `SFL`, `PAR`, `MOL`, `PCO`, `AIT`
Conic: `COP`, `COD`, `COE`, `COO`
Polyconic: `BON`
Quadrilateralized spherical cube: `TSC`, `CSC`¹, `QSC`
HEALPix: `HPX`, `XPH`

¹ CSC matches to ~9 mas due to WCSLIB storing its polynomial coefficients
  as 32-bit `float` while our implementation computes in 64-bit.

Unknown projection codes throw an informative error at transform time.

## FITS Loader Extensions

When the corresponding package is loaded, `from_header` accepts:

- `FITSIO.FITSHeader`
- `FITSIO.HDU`
- `FITSFiles.Card` vectors
- `FITSFiles.HDU`

The core package does not depend on either FITS loader at runtime.

## Known Limitations

- **AZP / SZP**: only the default (central perspective) parameter forms are
  implemented.  Non-default `PV` parameters are rejected at parse time.
- **Paper III spectral algorithms**: plain linear spectral axes work, but
  algorithm-coded axes such as `FREQ-LOG`, `WAVE-F2W`, etc. throw an explicit
  parse error.
- **`-TAB` table-lookup axes**: throw an explicit parse error.
- **Paper IV distortion lookup tables** (`CPDIS`, `D2IMDIS`, `DP`, `DQ`): throw
  an explicit parse error.
- **TPV / TPD polynomial distortion**: the `-TPV` and `-TPD` projection suffixes
  are not implemented.
- **Time and Stokes axes**: transform linearly but carry no physical
  interpretation (e.g., `MJDREF`, `DATE-OBS`, `TIMESYS`, polarization state).
- **Full WCS.jl API compatibility**: only a partial compatibility layer exists
  (`WCS`, `pix_to_world`, `world_to_pix`, mutating `!` variants).

Benchmarks for representative paths live under `benchmark/`.
