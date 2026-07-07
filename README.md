# FITSWCS.jl

FITSWCS.jl is a pure-Julia implementation of core FITS World Coordinate System
transformations. The main public interface for constructing WCS objects is `WCS`
which accepts dictionary input and returns a `WCSTransform`.
Other types (e.g., those from FITSIO.jl and FITSFiles.jl)
are parsed to dictionaries and then passed to `WCSTransform` in package extensions.

This package focuses on the published FITS WCS standard, though we are interested in contributions to
support non-standard FITS WCS features. Our implementations seek to optimize
scalar-path performance (allocation free `pixel_to_world` and `world_to_pixel`).
Batched versions of these functions are automatically multi-threaded and
take in `naxis x npoints` matrices, where each column is one coordinate.
Vectors of vectors (e.g., `[[1.0, 2.0], [2.0, 3.0]]`) are also supported.

## Quick Start

Pixel coordinates use the FITS WCS convention: pixel `1` is the center of the
first pixel.  Coordinate vectors are ordered by FITS WCS axis number.

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

### FITS Loader Extensions

When the corresponding package is loaded, `WCS` accepts FITSIO.jl and FITSFiles.jl
native types directly.  Headers are converted to `Dict` and passed through
to the core parser.

**FITSIO.jl:**

```julia
using FITSIO, FITSWCS

# From an in-memory FITSHeader.
fits = FITS("cube.fits")
hdr = read_header(fits[1])
wcs = WCS(hdr)

# Directly from an HDU; header is read internally.
wcs = WCS(fits[1])

# If external table data is specified in the header
# (-TAB, D2IM, CPDIS), pass the file as fobj.
wcs = WCS(fits[1]; fobj = fits)
```

**FITSFiles.jl:**

```julia
using FITSFiles, FITSWCS

# From a vector of parsed FITS cards.
fits = FITSFiles.read("cube.fits")
wcs = WCS(fits[1].cards)

# Directly from an HDU (cards extracted internally).
wcs = WCS(fits[1])

# If external table data is specified in the header,
# pass the full HDU vector as fobj.
wcs = WCS(fits[1]; fobj = fits)
```

## Programmatic Construction

In addition to `WCS(header)`, a keyword-based constructor accepts core WCS
vectors and matrices directly:

```julia
wcs = WCS(2; ctype=["RA---TAN", "DEC--TAN"], crpix=[512.0, 512.0],
          crval=[83.8221, -5.3911], cdelt=[-2.7778e-4, 2.7778e-4])
```

Supported keywords: `crpix`, `crval`, `cdelt`, `ctype`, `cunit`, `pc`, `cd`,
`crota`, `lonpole`, `latpole`, `radesys`, `equinox`, `wcsname`, `preserve_units`.

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
- time axes (`TIME`, `UTC`, `TAI`, `TDB`, `TT`, `TCG`, `TCB`, `LOCAL`):
  normalized to seconds; `MJDREF`, `TIMESYS`, `TREFPOS`, `TREFDIR`, `TIMEUNIT`
- observation metadata: `MJD-AVG`, `DATE-AVG`, `OBSGEO-X/Y/Z`
- celestial reference frame: `RADESYS`, `EQUINOX`
- WCS identification: `WCSNAME`

## Unit Conventions

Celestial coordinates are normalized to **degrees**, spectral coordinates to
**SI units** (Hz, m, m/s), and time coordinates to **seconds** at parse time.
`pixel_to_world` returns these
canonical units, and `world_to_pixel` also expects these canonical units.

Passing `preserve_units=true` to `WCS()` results in `pixel_to_world`
returning values in the original header `CUNIT` instead. Similarly,
`world_to_pixel` will expect input world coordinates in the proper CUNIT
as well.  Linear and Stokes axes remain
in the units encoded by their header linear transform.

```julia
wcs_asec = WCS(header; preserve_units=true)
world_asec = pixel_to_world(wcs_asec, [512.0, 512.0])  # returns arcsec, not degrees
```

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

## Spectral Coordinates

All FITS Paper III spectral types and algorithms are supported:

| S-type | Description | Algorithms |
|---|---|---|
| `FREQ` | Frequency | linear, `-LOG`, `-F2W`, `-F2V`, `-F2A` |
| `WAVE` | Vacuum wavelength | linear, `-LOG`, `-W2F`, `-W2V`, `-W2A` |
| `VELO` | Relativistic velocity | linear, `-V2F`, `-V2W`, `-V2A` |
| `AWAV` | Air wavelength | linear, `-A2F`, `-A2W`, `-A2V` |
| `AFRQ`, `ENER`, `WAVN` | Angular frequency, energy, wavenumber | linear |
| `VRAD`, `VOPT`, `ZOPT`, `BETA` | Derived velocity/redshift types | linear |
| Any wavelength type | As above | `-GRI` (grism in vacuum), `-GRA` (grism in air) |

Grism coordinates use the ideal disperser equation (grating interference +
prism refraction) from Paper III Section 5.  Seven PV parameters
(`PVi_0`–`PVi_6`) specify the disperser properties.

Air-wavelength conversions use the IUGG 1999 / Ciddor (1996) refractive-index
relation (Paper III eq. 4). Note that the WCSLIB code diverges from the published
paper and uses the Cox/Edlén (IAU 1957) relation.  The two differ by a nearly
constant ratio IUGG/Cox ≈ 1.000015 across the optical range.

Cross-type algorithms (e.g. `WAVE-F2W`) compute
the full X→P→S chain, including derivative scaling of the CD matrix entry.
Tabular spectral axes (`-TAB`) use binary-table coordinate arrays via `fobj`.

Rest-frequency/wavelength and reference-frame keywords (`RESTFRQ`, `RESTWAV`,
`SPECSYS`, etc.) are parsed and stored on the transform for downstream
frame-correction code, but no velocity-frame correction is performed by this
package.

## Benchmarks
Output from `benchmark/benchmarks.jl` run with on an Intel 12600K CPU with Julia 1.12.6.
*batch-* entries use the batched interface and are run
multi-threaded with 8 threads. *batch-100* use 100 coordinates
and *batch-1M* use 10^6 coordinates.
These benchmarks figures are re-run rarely so
performance on `main` may diverge from results here in the future.

Benchmark suite: pixel_to_world
| Benchmark | Median Time | Memory | Allocs |
|-----------|------------:|-------:|--------:|
| `2D-coupled-TAB/scalar` | 49.050 ns | 0 bytes | 0 |
| `3D-cube-TAB/batch-100` | 6.267 μs | 7.71 KiB | 45 |
| `3D-cube-TAB/scalar` | 87.500 ns | 0 bytes | 0 |
| `3D-cube-spec/scalar` | 71.580 ns | 0 bytes | 0 |
| `3D-cube/scalar` | 69.270 ns | 0 bytes | 0 |
| `AIT/scalar` | 84.030 ns | 0 bytes | 0 |
| `TAN-SIP-PaperIV/scalar` | 157.830 ns | 0 bytes | 0 |
| `TAN-SIP/scalar` | 101.000 ns | 0 bytes | 0 |
| `TAN/batch-100/Float32` | 4.947 μs | 5.30 KiB | 44 |
| `TAN/batch-100/Float64` | 5.531 μs | 6.02 KiB | 44 |
| `TAN/batch-1M/Float32` | 9.378 ms | 7.63 MiB | 45 |
| `TAN/batch-1M/Float64` | 14.073 ms | 15.26 MiB | 45 |
| `TAN/scalar` | 65.080 ns | 0 bytes | 0 |
| `TAN/scalar/SVector Float32` | 49.320 ns | 0 bytes | 0 |
| `TAN/scalar/SVector Float64` | 64.710 ns | 0 bytes | 0 |
| `TAN/scalar/Tuple` | 65.450 ns | 0 bytes | 0 |
| `TAN/scalar/preserve_units` | 65.310 ns | 0 bytes | 0 |
| `grism/AWAV-GRA/scalar` | 16.330 ns | 0 bytes | 0 |

Benchmark suite: world_to_pixel
| Benchmark | Median Time | Memory | Allocs |
|-----------|------------:|-------:|--------:|
| `2D-coupled-TAB/scalar` | 50.110 ns | 0 bytes | 0 |
| `3D-cube-TAB/batch-100` | 6.754 μs | 7.71 KiB | 45 |
| `3D-cube-TAB/scalar` | 99.030 ns | 0 bytes | 0 |
| `3D-cube-spec/scalar` | 82.890 ns | 0 bytes | 0 |
| `3D-cube/scalar` | 71.060 ns | 0 bytes | 0 |
| `AIT/scalar` | 75.060 ns | 0 bytes | 0 |
| `TAN-SIP-PaperIV/scalar` | 528.620 ns | 0 bytes | 0 |
| `TAN-SIP/scalar` | 236.460 ns | 0 bytes | 0 |
| `TAN/batch-100/Float32` | 4.894 μs | 5.30 KiB | 44 |
| `TAN/batch-100/Float64` | 5.781 μs | 6.02 KiB | 44 |
| `TAN/batch-1M/Float32` | 10.716 ms | 7.63 MiB | 45 |
| `TAN/batch-1M/Float64` | 13.327 ms | 15.27 MiB | 45 |
| `TAN/scalar` | 62.560 ns | 0 bytes | 0 |
| `TAN/scalar/preserve_units` | 62.670 ns | 0 bytes | 0 |
| `grism/AWAV-GRA/scalar` | 15.780 ns | 0 bytes | 0 |

Benchmark suite: parsing
| Benchmark | Median Time | Memory | Allocs |
|-----------|------------:|-------:|--------:|
| `WCS/3D-cube` | 368.853 μs | 68.83 KiB | 1850 |
| `WCS/3D-cube-TAB` | 320.728 μs | 69.78 KiB | 1839 |
| `WCS/3D-cube-spec` | 720.092 μs | 78.62 KiB | 2115 |
| `WCS/AIT` | 246.903 μs | 44.05 KiB | 1192 |
| `WCS/TAN` | 328.439 μs | 44.09 KiB | 1193 |
| `WCS/TAN-SIP` | 450.694 μs | 66.75 KiB | 1737 |
| `WCS/grism/AWAV-GRA` | 294.696 μs | 37.89 KiB | 984 |

## References
[This NASA page](https://fits.gsfc.nasa.gov/fits_wcs.html) links to most of the
relevant documentation on the FITS standard.
[wcslib](https://www.atnf.csiro.au/computing/software/wcs/wcslib/index.html)
is the de-facto standard software implementation of the FITS WCS standard
written by Mark Calabretta.

## Known Limitations

- Velocity-frame correction math (barycentric/LSRK conversion using `SPECSYS`,
  `SSYSOBS`, `VELOSYS`, `ZSOURCE`, `SSYSSRC`, `MJD-AVG`, `OBSGEO-X/Y/Z`)
  is not implemented; the keywords are parsed and stored for downstream
  use or a future package-level correction step.