# FITS WCS Implementation Notes

## 2026-06-14: Core Parsing And Projection Corrections

References reviewed:

- FITS WCS Paper I, "Basic concepts": linear transform pipeline, PC/CD defaults, WCSAXES dimensionality, 1-based pixel coordinate convention, and keyword defaults.
- FITS WCS Paper I, "Specification of units": angular unit strings `deg`, `arcmin`, `arcsec`, and `mas`.
- Astropy/wcslib search results:
  - `astropy/wcs/docstrings.py`: PC/CD/CROTA precedence and `LONPOLE`/`LATPOLE` behavior notes.
  - `cextern/wcslib/C/wcs.c`: CD-to-PC conversion, CROTA-to-PC construction, and celestial setup flow.
  - `astropy/modeling/projections.py`: TAN, STG, and CAR projection formula references.
  - `astropy/modeling/rotations.py`: Paper II spherical rotation convention background.
  - `astropy/wcs/src/sip.c`, `astropy/wcs/wcs.py`, and `astropy/modeling/polynomial.py`: SIP parsing and pipeline behavior for future work.
- ADS result `2005ASPC..347..491S`: Shupe et al. 2005, "The SIP Convention for Representing Distortion in FITS Image Headers."

Decisions made:

- Keep `from_header` dictionary-first and front-end agnostic.
- Infer `WCSAXES` from indexed WCS keywords when `WCSAXES` is absent but `NAXIS` is present, matching Paper I.
- Preserve the existing explicit-error behavior when both `NAXIS` and `WCSAXES` are missing from the simple dictionary input.
- Use Paper I defaults for core WCS keywords, including `CRPIXi = 0.0`.
- Reject mixed `PCi_ja` and `CDi_ja` matrices; allow `CDi_ja` to coexist with `CDELTi` because new readers ignore `CDELTi` in CD form.
- Parse strict FITS 4-3 `CTYPE` values and retain the base projection code for suffix conventions such as `-SIP`.
- Convert celestial `CRVAL` and matrix rows to degrees for the public API, but leave non-celestial axes in their header units.
- Fix AIT inverse projection branch handling by using the wrapped native longitude consistently.
- Parse SIN slant orthographic parameters from `PV<lat>_1` and `PV<lat>_2`, following the FITS `PVi_ma` convention where `i` is the latitude-like axis.
- Add an isolated SIP distortion layer:
  - parse `A_ORDER`/`B_ORDER`, `A_i_j`/`B_i_j`, and optional `AP_ORDER`/`BP_ORDER`, `AP_i_j`/`BP_i_j`;
  - require explicit `CRPIX1` and `CRPIX2` when SIP keywords are present;
  - apply forward SIP before the linear matrix, matching Astropy's pixel-to-focal pipeline;
  - use inverse coefficients when present and fixed-point iteration otherwise.
- Add package extensions for FITSIO.jl and FITSFiles.jl:
  - keep FITSIO.jl and FITSFiles.jl as weak dependencies;
  - adapt `FITSIO.FITSHeader`, `FITSIO.HDU`, FITSFiles card vectors, and `FITSFiles.HDU` into the core dictionary parser;
  - keep FITS parsing and HDU I/O behavior owned by the upstream packages.

Deferred:

- Reference-comparison fixtures generated from Astropy or wcslib.
- Real SIP header fixtures and Astropy-generated expected values.
- Broader Paper II validation for non-zenithal native pole edge cases.

## 2026-06-14: Mixed-Axis And Compliance Tracking

References reviewed:

- Astropy `astropy/wcs/wcsapi/tests/test_fitswcs.py`: spectral-cube tests with
  celestial axes split around a spectral axis.
- Astropy `astropy/wcs/tests/test_utils.py`: axis-correlation examples for
  `RA---TAN`, `FREQ`, `DEC--TAN` WCS objects.
- Astropy `astropy/wcs/docstrings.py`: `naxis` derivation from `NAXIS`,
  `WCSAXESa`, and the highest valid parameterized WCS keyword.
- WCS.jl public transform tests and docs: matrix-shaped coordinates use axes
  in rows and coordinate points in columns.

Decisions made:

- Keep mixed celestial/non-celestial support in the existing full-vector API:
  each output coordinate remains in its FITS WCS axis position.
- Treat spectral, time, and Stokes physical interpretation as future work; for
  now non-celestial axes use the linear Paper I transform.
- Add the compliance matrix required by the implementation plan and mark missing
  reference fixtures, Paper III features, Paper IV lookup distortions, API
  compatibility, and benchmarks as deferred.
- Add basic tests for time and Stokes axes as linear Paper I axes, while keeping
  physical interpretation explicitly deferred.
- Add a README documenting supported keywords, projections, loader extensions,
  pixel/axis conventions, and known limitations.
- Add a narrow WCS.jl-style API layer (`WCS`, `pix_to_world`, `world_to_pix`)
  that delegates to the canonical FITSWCS.jl API while retaining the FITS
  1-based pixel convention.

## 2026-06-14: CEA Projection

References reviewed:

- Astropy `astropy/modeling/projections.py`: CEA equations for pixel-to-sky and
  sky-to-pixel projection models.
- Astropy vendored wcslib `ceaset` comments and setup logic: `PV<lat>_1`
  defaults to 1 and must be in `(0, 1]`.
- WCS.jl header tests: a real-header-style CEA example appears in the existing
  wrapper package test suite.

Decisions made:

- Add `CEA(lambda)` as a small cylindrical projection type, with `lambda`
  parsed from `PV<lat>_1`.
- Use the same default native fiducial point as other cylindrical projections:
  `phi0 = 0`, `theta0 = 0`.
- Enforce the finite inverse projection domain rather than silently clamping
  coordinates that are outside the valid CEA latitude range.

## 2026-06-15: Astropy Regression Cross-Check

References reviewed:

- Astropy `astropy/wcs/docstrings.py`: `all_pix2world` applies detector,
  SIP, distortion-table, and wcslib core transforms; the `origin` argument
  controls pixel-coordinate origin.
- Astropy `astropy/wcs/wcs.py`: `all_world2pix` iteratively inverts the total
  transform, including distortion corrections.
- Astropy `astropy/wcs/utils.py`: high-level helpers pass through
  `all_pix2world`/`all_world2pix` in `mode="all"`.

Decisions made:

- Add `test/regression_astropy.py` as an optional Python-side integration
  verifier instead of adding Astropy to Julia test dependencies.
- Use `origin=1` for both `all_pix2world` and `all_world2pix`, matching the
  FITS 1-based convention used by FITSWCS.jl and the wcslib fixtures.
- Compare longitude-like first axes modulo 360 degrees to avoid false
  mismatches from phase-equivalent RA/longitude values.

Result:

- Astropy 6.1.7 matched all checked wcslib regression values in both
  directions: 9 cases, 31 pixel/world points, 0 mismatches.

## 2026-06-15: ARC, ZEA, And CEA Reference Fixtures

References reviewed:

- FITS WCS Paper II equations already cited in the projection implementation:
  ARC Eq. 46 and ZEA Eq. 52.
- Astropy 6.1.7 `astropy.wcs.WCS` with `all_pix2world(..., origin=1)` and
  `all_world2pix(..., origin=1)` for end-to-end reference values.

Decisions made:

- Add `test/regression_astropy_values.jl` for projections that did not yet
  have stored numeric external fixtures.
- Keep these fixtures separate from `test/regression_wcslib.jl` because their
  immediate provenance is Astropy, not WCS.jl/wcslib.
- Use non-reference and off-axis pixels for ARC and ZEA, and use CEA with
  `PV2_1 = 0.75` so projection-parameter parsing is covered.

Result:

- Added 30 Julia assertions comparing ARC, ZEA, and CEA pixel/world values and
  inverse transforms against Astropy-generated fixtures.

## 2026-06-15: Celestial Unit Reference Fixtures

References reviewed:

- Astropy 6.1.7 `astropy.wcs.WCS` behavior for headers with celestial
  `CUNIT1`/`CUNIT2` set to `arcsec` and `rad`.
- FITS WCS Paper I unit handling as already implemented in `unit_to_deg`.

Decisions made:

- Keep FITSWCS celestial world coordinates degree-valued at the public API
  boundary, matching Astropy's normalization of celestial CUNIT values.
- Add stored Astropy fixtures for both arcsecond and radian celestial headers,
  covering CRVAL and CDELT conversion through full TAN transforms.

Result:

- Added 12 Julia assertions comparing celestial CUNIT arcsec/rad pixel/world
  values and inverse transforms against Astropy-generated fixtures.

## 2026-06-15: Mutating WCS.jl-Style API Aliases

References reviewed:

- WCS.jl `src/WCS.jl`: `pix_to_world!(wcs, pixcoords, worldcoords)` and
  `world_to_pix!(wcs, worldcoords, pixcoords)` fill caller-provided arrays.
- WCS.jl docs: matrix-shaped coordinates use axes in rows and coordinate
  points in columns, matching the FITSWCS batch API.

Decisions made:

- Add basic mutating aliases that delegate through the canonical FITSWCS API
  and copy into the provided output arrays.
- Keep WCS.jl's optional status and intermediate work-array keyword API
  deferred; FITSWCS does not yet expose equivalent status information.
- Retain FITSWCS's FITS 1-based pixel-coordinate convention.

Result:

- Added vector and matrix tests for `pix_to_world!` and `world_to_pix!`,
  including output-shape validation.

## 2026-06-15: WCS.jl-Style Keyword Constructor

References reviewed:

- WCS.jl `WCSTransform(naxis; kwds...)`: programmatic construction accepts
  property-style vectors and matrices before calling wcslib setup.
- WCS.jl public property list: core fields include `crpix`, `crval`, `cdelt`,
  `ctype`, `cunit`, `pc`, `cd`, and legacy `crota`.

Decisions made:

- Add `WCSTransform(naxis; kwds...)` as a compatibility constructor for core
  WCS vectors and matrices.
- Translate supported constructor keywords into ordinary FITS header keys and
  call `from_header`, so validation, projection setup, SIP rejection, and unit
  handling remain owned by the parser.
- Reject specialized wcslib fields such as `restfrq` for now instead of
  accepting metadata that FITSWCS would not yet use.

Result:

- Added constructor tests for diagonal and PC-matrix linear transforms plus
  malformed or unsupported keyword errors.

## 2026-06-15: Vector Batch Transform Inputs

References reviewed:

- FITSWCS transform docstrings already described batch inputs as matrices or
  vectors of coordinate vectors.
- Existing matrix batch tests for mixed celestial/spectral WCS.

Decisions made:

- Keep matrices as the preferred dense batch representation.
- Add vector-of-vectors methods as a convenience path that simply maps each
  coordinate through the scalar transform, preserving one output vector per
  input coordinate.

Result:

- Added mixed-axis tests proving vector batches agree with matrix batches in
  both pixel-to-world and world-to-pixel directions.

## 2026-06-15: Explicit -TAB Deferral

References reviewed:

- Astropy/wcslib tabular-coordinate code paths: `tabprm` objects represent
  lookup-table coordinates and validate table dimensionality, axis mappings,
  coordinate arrays, and associated WTB metadata.
- Local Astropy 6.1.7 behavior for a minimal `FREQ-TAB` header without table
  metadata errors during WCS construction rather than applying a linear axis.

Decisions made:

- Keep full Paper III `-TAB` lookup support deferred.
- Reject any parsed `CTYPEia` algorithm code of `TAB` during `from_header`, so
  unsupported lookup coordinates are not silently treated as ordinary linear
  axes.
- Apply the same check to alternate WCS descriptions selected with `alt`.

Result:

- Added parser tests for `FREQ-TAB` recognition and primary/alternate
  `from_header` rejection.

## 2026-06-15: Explicit Paper IV Lookup-Distortion Deferral

References reviewed:

- Astropy `astropy/wcs/wcs.py`: the full pixel pipeline applies detector to
  image-plane lookup correction, SIP, FITS WCS distortion paper lookup
  correction, then wcslib core WCS.
- Astropy `astropy/wcs/docstrings.py`: `cpdis1`/`cpdis2` are pre-linear
  distortion lookup tables and `det2im1`/`det2im2` are detector to image-plane
  lookup corrections.
- Astropy/wcslib distortion code: Paper IV distortion parameter keywords use
  DP/DQ families and validate distortion parameters before use.

Decisions made:

- Keep Paper IV lookup-table distortion support deferred.
- Reject recognized lookup-distortion metadata (`CPDIS*`, `D2IMDIS*`,
  `D2IMERR*`, `AXISCORR`, `DP*.*`, and `DQ*.*`) during header parsing so
  distorted headers are not silently transformed as undistorted WCS.
- Keep SIP distortion unaffected because it is already implemented.

Result:

- Added parser tests for each guarded lookup-distortion keyword family and a
  regression check that SIP still parses through the supported distortion path.

## 2026-06-15: Explicit Spectral Algorithm Deferral

References reviewed:

- Astropy/wcslib spectral code paths: wcslib has dedicated spectral
  translation and spectral transform logic for algorithm-coded CTYPE values.
- Local Astropy 6.1.7 behavior: plain `FREQ` and `WAVE` axes are linear, while
  `FREQ-LOG` applies a non-linear spectral transform.

Decisions made:

- Continue supporting plain linear spectral axes such as `FREQ` and `WAVE`.
- Keep Paper III physical/spectral algorithms deferred.
- Reject non-celestial strict 4-3 CTYPE algorithm codes such as `FREQ-LOG` and
  `WAVE-F2W`, so those axes are not silently linearized.

Result:

- Added parser tests for spectral algorithm recognition and primary/alternate
  `from_header` rejection.

## 2026-06-15: Split Celestial/Spectral Reference Fixtures

References reviewed:

- Astropy 6.1.7 `astropy.wcs.WCS` behavior for a 3D cube with axis order
  `RA---TAN`, `FREQ`, `DEC--TAN`.
- Existing FITSWCS higher-dimensional tests for split celestial axes.

Decisions made:

- Keep linear spectral axes in mixed cubes as plain Paper I linear axes.
- Add stored Astropy values for a split celestial/spectral axis order so axis
  placement is checked against an external implementation, not just round-trip
  consistency.

Result:

- Added 8 Julia assertions comparing split-axis pixel/world values and inverse
  transforms against Astropy-generated fixtures.

## 2026-06-15: Linear Time And Stokes Reference Fixtures

References reviewed:

- Astropy 6.1.7 `astropy.wcs.WCS` behavior for a 4D WCS with `RA---TAN`,
  `DEC--TAN`, `TIME`, and `STOKES` axes.
- Existing FITSWCS higher-dimensional tests for basic linear time and Stokes
  axes.

Decisions made:

- Keep physical time and Stokes interpretation deferred.
- Add stored Astropy values for the supported subset: TIME and STOKES as plain
  Paper I linear axes embedded in a mixed celestial WCS.

Result:

- Added 8 Julia assertions comparing mixed celestial/TIME/STOKES pixel/world
  values and inverse transforms against Astropy-generated fixtures.
