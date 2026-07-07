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

World coordinates are returned in canonical units by default (degrees for
celestial axes, SI for spectral axes).  Pass `preserve_units=true` to `WCS()`
to return values in the original header `CUNIT` instead:

```julia
wcs_asec = WCS(header; preserve_units=true)
world_asec = pixel_to_world(wcs_asec, [512.0, 512.0])  # arcsec, not degrees
```

## Programmatic Construction

In addition to `WCS(header)`, a keyword-based constructor accepts core WCS
vectors and matrices directly:

```julia
wcs = WCS(2; ctype=["RA---TAN", "DEC--TAN"], crpix=[512.0, 512.0],
          crval=[83.8221, -5.3911], cdelt=[-2.7778e-4, 2.7778e-4])
```

Supported keywords: `crpix`, `crval`, `cdelt`, `ctype`, `cunit`, `pc`, `cd`,
`crota`, `lonpole`, `latpole`, `preserve_units`.

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
- spectral: `RESTFRQ`, `RESTWAV`, `SPECSYS`, `SSYSOBS`, `VELOSYS`, `ZSOURCE`,
  `SSYSSRC`
- tabular (`-TAB`): `PSi_0`, `PSi_1`, etc. with binary table coordinate arrays
  via `fobj`
- observation metadata: `MJD-AVG`, `DATE-AVG`, `OBSGEO-X/Y/Z`
- celestial reference frame: `RADESYS`, `EQUINOX`
- WCS identification: `WCSNAME`

Celestial coordinates are normalized to **degrees** and spectral coordinates to
**SI units** (Hz, m, m/s) at parse time.  `pixel_to_world` returns these
canonical units by default; pass `preserve_units=true` to `WCS()` to return
values in the original header `CUNIT` instead.  Linear, time, and Stokes axes
remain in the units encoded by their header linear transform.

## Spectral Coordinates (Paper III)

All FITS Paper III spectral types and algorithms are supported:

| S-type | Description | Algorithms |
|---|---|---|
| `FREQ` | Frequency | linear, `-LOG`, `-F2W`, `-F2V`, `-F2A` |
| `WAVE` | Vacuum wavelength | linear, `-LOG`, `-W2F`, `-W2V`, `-W2A` |
| `VELO` | Relativistic velocity | linear, `-V2F`, `-V2W`, `-V2A` |
| `AWAV` | Air wavelength | linear, `-A2F`, `-A2W`, `-A2V` |
| `AFRQ`, `ENER`, `WAVN` | Angular frequency, energy, wavenumber | linear |
| `VRAD`, `VOPT`, `ZOPT`, `BETA` | Derived velocity/redshift types | linear |

Air-wavelength conversions use the IUGG 1999 / Ciddor (1996) refractive-index
relation (Paper III eq. 4).  Cross-type algorithms (e.g. `WAVE-F2W`) compute
the full X→P→S chain, including derivative scaling of the CD matrix entry.
Tabular spectral axes (`-TAB`) use binary-table coordinate arrays via `fobj`.

Rest-frequency/wavelength and reference-frame keywords (`RESTFRQ`, `RESTWAV`,
`SPECSYS`, etc.) are parsed and stored on the transform for downstream
frame-correction code, but no velocity-frame correction is performed by this
package.

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

- Velocity-frame correction math (barycentric/LSRK conversion using `SPECSYS`,
  `SSYSOBS`, `VELOSYS`, `ZSOURCE`, `SSYSSRC`, `MJD-AVG`, `OBSGEO-X/Y/Z`)
  is not yet implemented; the keywords are parsed and stored for downstream
  use or a future package-level correction step.
- **Grism** algorithm codes `GRI`/`GRA` are not yet implemented.
- **Time and Stokes axes**: transform linearly but carry no physical
  interpretation (e.g., `MJDREF`, `DATE-OBS`, `TIMESYS`, polarization state).

Benchmarks for representative paths live under `benchmark/`.
