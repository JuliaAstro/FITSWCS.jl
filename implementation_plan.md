# Agentic Coding Objective: Pure-Julia FITS WCS Implementation

## Objective

Create a pure-Julia implementation of the FITS World Coordinate System standard, aiming for feature parity with the existing `WCS.jl` package where practical, but without wrapping `wcslib`. The implementation should minimize dependencies, maximize standards coverage, provide excellent correctness-oriented test coverage, and remain simple enough to be maintainable by the JuliaAstro ecosystem over the long term. The package should be agnostic with respect to the front-end FITS loader; it should use package extensions for both FITSIO.jl and FITSFiles.jl to support loading WCS information from the header types exposed by both of these packages.

The implementation should be built incrementally. The agent should repeatedly:

1. Collect context from standards and reference implementations.
2. Implement a coherent feature slice.
3. Add tests that validate both interface behavior and numerical correctness.
4. Refactor for simplicity and maintainability.
5. Benchmark performance-sensitive paths when appropriate.
6. Continue to the next feature slice.

Do not attempt to implement the entire standard in one pass.

---

## Primary References

Use the following as truth/reference sources:

1. FITS WCS standard documents:

   * https://fits.gsfc.nasa.gov/fits_wcs.html
   * FITS WCS Papers I–IV where relevant.
   * Pay particular attention to Paper I and Paper II before implementing core linear WCS and celestial projections.

2. Existing Julia package:

   * `WCS.jl`, which currently wraps the C `wcslib`.
   * Use this to understand the existing Julia-facing API and expected user workflows.
   * Do not copy C-wrapper internals; use it mainly for interface and behavioral comparison.

3. Reference implementation:

   * Astropy WCS implementation.
   * Use `semble-search` on the Astropy repository:

     * https://github.com/astropy/astropy
   * Search specifically for:

     * FITS WCS header parsing.
     * `CTYPE`, `CRPIX`, `CRVAL`, `CDELT`, `CD`, `PC`, `CROTA`.
     * Celestial projection implementations.
     * SIP distortion handling.
     * WCSAxes / high-level API expectations only as background, not as implementation target.
     * Test files and known edge cases.
   * Astropy is not the formal standard, but it is useful as a compatibility and edge-case reference.

4. `wcslib` behavior:

   * Use as a behavioral reference where needed, but do not bind to it.
   * If available locally, compare outputs against `WCS.jl`/`wcslib` for generated and real FITS headers.

---

## Design Priorities

The package should prioritize, in order:

1. Correctness with respect to the FITS WCS standard.
2. Clear, maintainable, idiomatic Julia code.
3. Minimal dependency footprint.
4. Good test coverage with numerical correctness checks.
5. Good performance for common pixel/world transformations.
6. Extensibility toward broader WCS/GWCS-style functionality later.

Avoid premature abstraction. Prefer simple, explicit data structures and algorithms unless a more generic design is clearly justified. Avoid overly specific types in method definitions.

---

## Dependency Policy

Minimize dependencies. Before adding any dependency, evaluate:

* Is it necessary for correctness?
* Is it small and stable?
* Is it already common in the JuliaAstro ecosystem?
* Can the needed functionality be implemented cleanly in a few lines?
* Would vendoring or local implementation be clearer?

Expected acceptable dependencies may include:

* `LinearAlgebra` from the standard library.
* `Test` for testing.
* Possibly `StaticArrays.jl` only if there is a demonstrated performance or clarity benefit, but avoid it initially unless needed.

Avoid large dependencies for small utilities.

---

## Initial Repository Survey

You are beginning a new package in /home/cgarling/Development/julia/FITSWCS.jl. Edit and create files only in this directory.

---

## Research Phase

Use `semble-search` on Astropy before implementing each major feature area.

Suggested searches:

* `WCS CTYPE CRPIX CRVAL CDELT PC CD`
* `astropy wcs celestial projection TAN SIN AIT CAR`
* `astropy wcs sip distortion`
* `astropy wcs tests header`
* `astropy wcs wcslib comparison`
* `astropy coordinates pixel_to_world world_to_pixel`
* `astropy wcs FITS WCS Paper I`
* `astropy wcs distortion lookup table`
* `astropy wcs test_sip`
* `astropy wcs test_profiling` or performance-related tests

For each feature area, record brief notes in an implementation plan file, for example `docs/dev/implementation_plan.md` or `notes/wcs_implementation.md`, including:

* Relevant standard sections.
* Astropy files inspected.
* Important edge cases.
* Decisions made.
* Features deferred.

Do not let the notes become a substitute for tests or code.

---

## Core Architecture

Start with a small number of clear types.

Possible initial design:

```julia
struct WCS
    axes::Vector{WCSAxis}
    transform::WCSTransform
    metadata::Dict{String,Any}
end

struct WCSAxis
    ctype::String
    cunit::Union{Nothing,String}
    crpix::Float64
    crval::Float64
    cdelt::Float64
end

struct LinearTransform
    crpix::Vector{Float64}
    crval::Vector{Float64}
    matrix::Matrix{Float64}
end
```

The exact design may change, but keep the following conceptual layers separate:

1. Header parsing.
2. Internal WCS representation.
3. Pixel-to-intermediate coordinate transform.
4. Projection-specific transform.
5. World coordinate output.
6. Inverse transform.
7. Distortion corrections.

Avoid coupling FITS header parsing directly to numerical transform code.

---

## Implementation Milestones

### Milestone 1: Basic Header Parsing

Implement parsing for the core linear WCS keywords:

* `NAXIS`
* `WCSAXES`
* `CTYPEia`
* `CUNITia`
* `CRPIXia`
* `CRVALia`
* `CDELTia`
* `PCi_ja`
* `CDi_ja`
* `CROTAia`, if supported as a legacy compatibility path
* alternate WCS versions, if simple enough to support early

Requirements:

* Parse from a simple dictionary-like header representation first.
* Do not require FITS IO initially unless already present.
* Validate dimensions.
* Provide clear errors for malformed headers.
* Preserve unknown/non-core metadata where useful.

Tests:

* Minimal valid 1D, 2D, and 3D headers.
* Missing optional keywords.
* Invalid matrix dimensions.
* CD vs PC/CDELT precedence.
* Legacy CROTA behavior if implemented.
* Alternate WCS suffixes if implemented.

Correctness tests:

* Hand-computed linear transformations.
* Round-trip pixel → world → pixel for simple linear cases.
* Comparison to `WCS.jl`/`wcslib` where available.

---

### Milestone 2: Linear Pixel/World Transform

Implement basic linear WCS:

```julia
world = pixel_to_world(wcs, pixel)
pixel = world_to_pixel(wcs, world)
```

Also consider batch APIs:

```julia
pixel_to_world(wcs, pixels::AbstractMatrix)
world_to_pixel(wcs, worlds::AbstractMatrix)
```

Define the axis convention explicitly. FITS/WCS conventions can be confusing because FITS pixel axes, Julia array axes, and user-facing coordinates may differ.

Requirements:

* Document whether pixel coordinates are 1-indexed FITS-style coordinates or Julia array indices.
* Prefer exposing FITS/WCS pixel coordinates directly, because FITS WCS standards use 1-based pixel coordinates.
* If helper methods are added for Julia array indices, name them explicitly.
* Avoid hidden axis reversals.

Tests:

* Identity transform.
* Translation via `CRPIX`.
* Scaling via `CDELT`.
* Rotation/skew via `PC` or `CD`.
* Non-square axis counts where standard permits.
* Batch and scalar APIs produce identical results.

Correctness tests:

* Hand-computed transformations.
* Random linear WCS round-trip tests.
* Comparison against Astropy or `wcslib` for generated headers.

---

### Milestone 3: Celestial WCS Foundation

Implement parsing and internal representation for celestial WCS axes:

* Longitude axis.
* Latitude axis.
* `CTYPE` parsing such as `RA---TAN`, `DEC--TAN`.
* Native longitude/latitude concepts from the FITS WCS standard.
* Unit handling for degrees initially.

Start with `TAN` projection first.

Requirements:

* Keep projection logic separate from header parsing.
* Define a projection interface such as:

```julia
project(proj, lon, lat)
deproject(proj, x, y)
```

or equivalent.

* Be explicit about radians vs degrees internally.
* Prefer radians internally for trig-heavy code, with degree conversion at API boundaries.

Tests:

* Parse celestial `CTYPE`.
* Reject incompatible celestial axis pairs.
* Pixel → sky and sky → pixel for simple TAN headers.
* Round-trip tests near reference point.
* Behavior near projection singularities.

Correctness tests:

* Compare against examples from FITS WCS papers.
* Compare against Astropy and/or `wcslib` for representative headers.
* Include tolerances appropriate for floating-point projection math.

---

### Milestone 4: Common Celestial Projections

After `TAN`, add projections incrementally. Suggested order:

1. `TAN`
2. `SIN`
3. `CAR`
4. `AIT`
5. `STG`
6. `ARC`
7. `ZEA`
8. Other projections as needed by real HST/ground-based headers.

For each projection:

* Read the relevant FITS WCS paper section.
* Inspect Astropy tests and implementation behavior using `semble-search`.
* Implement forward transform.
* Implement inverse transform.
* Add standard-derived tests.
* Add comparison tests against Astropy/`wcslib`.

Do not add many projections without tests. One well-tested projection is better than five unvalidated ones.

---

### Milestone 5: Real FITS Header Compatibility

Collect a small suite of real-world FITS headers, preferably without large image data. Include examples from:

* HST ACS/WFC
* HST WFC3/UVIS
* HST WFC3/IR
* simple ground-based TAN images
* images with CD matrices
* images with PC/CDELT matrices
* headers using SIP distortion, if available

Tests should use header snippets or small fixture files.

Correctness tests:

* Compare pixel/world coordinates for selected points against Astropy and/or `WCS.jl`.
* Include center, corners, and random interior points.
* Include round-trip tests.
* Store expected values with clear provenance.

Avoid tests that merely check that code runs.

---

### Milestone 6: SIP Distortion

Implement SIP only after the core linear/celestial path is stable.

Research first:

* FITS SIP convention.
* Astropy SIP parser and tests.
* Real headers using SIP.

Implement:

* SIP coefficient parsing:

  * `A_ORDER`, `B_ORDER`
  * `A_i_j`, `B_i_j`
  * inverse coefficients if available:

    * `AP_ORDER`, `BP_ORDER`
    * `AP_i_j`, `BP_i_j`
* Forward SIP correction.
* Inverse SIP handling:

  * Use inverse polynomial if present.
  * Otherwise use iterative inversion with clear convergence criteria.

Tests:

* Polynomial evaluation correctness.
* Known SIP header comparison against Astropy.
* Round-trip tests at center, edges, and corners.
* Convergence/failure tests for inverse iteration.
* No-SIP path remains unaffected.

Keep SIP implementation isolated from the core projection code.

---

### Milestone 7: Higher-Dimensional WCS

Support non-celestial and mixed axes incrementally:

* Spectral axes.
* Time axes.
* Stokes axes.
* Mixed celestial + spectral cubes.

Start with parsing and linear transforms. Add physical interpretation only where standards coverage is clear.

Tests:

* 3D cube with RA/DEC/FREQ.
* Pixel/world round-trip.
* Axis order handling.
* Subsetting, if implemented.

Do not overbuild high-level coordinate object integration early.

---

### Milestone 8: API Compatibility Layer

Once core functionality is solid, compare against `WCS.jl` public API.

Decide which APIs should be supported directly, which should be adapted, and which should be deprecated or omitted.

Possible API goals:

```julia
wcs = WCS(header)
world = pix_to_world(wcs, x, y)
pixel = world_to_pix(wcs, ra, dec)
```

Also consider more Julian forms:

```julia
pixel_to_world(wcs, (x, y))
world_to_pixel(wcs, (ra, dec))
```

Requirements:

* Keep public names clear and stable.
* Avoid exposing internal representation too early.
* Document axis order carefully.
* Add migration notes if replacing existing `WCS.jl` internals.

Tests:

* API behavior.
* Type stability where practical.
* Error behavior.
* Compatibility with existing examples.

---

## Testing Strategy

Tests must validate correctness, not only interface behavior.

Use several categories of tests:

### 1. Unit Tests

For small pure functions:

* Header keyword parsing.
* Matrix construction.
* Projection math.
* Polynomial evaluation.
* Unit conversion.
* Axis classification.

### 2. Hand-Computed Tests

For simple cases where expected answers are obvious:

* Identity WCS.
* Pure shift.
* Pure scale.
* Simple rotation.
* Simple TAN near reference point.

These should not depend on external libraries.

### 3. Standard Example Tests

Use examples from FITS WCS papers where available.

Each test should cite the source section in a comment.

### 4. Reference Comparison Tests

Compare against one or both of:

* Astropy
* `WCS.jl`/`wcslib`

These tests can be generated offline and checked into fixtures as expected numerical values, avoiding Python/runtime dependency in normal CI.

If Python/Astropy comparison tests are added, mark them as optional or integration tests.

### 5. Property Tests

Use randomized but well-conditioned WCS objects:

* Pixel → world → pixel round-trip.
* World → pixel → world round-trip.
* Matrix inverse consistency.
* Projection inverse consistency away from singularities.

Avoid random tests that are nondeterministic; use fixed RNG seeds.

### 6. Real Header Regression Tests

Use real-world headers from astronomical instruments.

For each header, test selected pixel positions:

* Reference pixel.
* Image center.
* Four corners.
* Several interior points.

Expected values should be generated from a trusted implementation and stored with provenance.

### 7. Error Tests

Test malformed inputs:

* Missing required keywords.
* Invalid axis counts.
* Singular matrices.
* Unknown projection codes.
* Incompatible celestial axis pairs.
* Non-convergent inverse distortion solve.

Errors should be informative.

---

## Performance Strategy

Do not sacrifice clarity prematurely, but keep performance in mind.

Performance-sensitive paths:

* Scalar pixel → world transformation.
* Batch pixel → world transformation.
* Matrix application.
* Projection trig functions.
* SIP polynomial evaluation.
* Iterative inverse distortion correction.

Guidelines:

* Start with clear code.
* Benchmark before optimizing.
* Use allocation tests for hot scalar paths.
* Avoid unnecessary heap allocation in inner loops.
* Consider specialized small-dimensional methods for common 2D celestial WCS.
* Use `@inbounds` only after tests are strong.
* Avoid overly clever generated functions unless clearly justified.

Add benchmarks for:

* Single coordinate transform.
* Large batch transform.
* TAN projection.
* SIP transform.
* Comparison against current `WCS.jl` if available.

Benchmarks should not be required for normal tests. Benchmarks should follow the proper style needed to be run on CI with AirSpeedVelocity.jl.

---

## Documentation Requirements

For each major feature, add documentation explaining:

* Supported FITS WCS keywords.
* Axis order conventions.
* Pixel coordinate conventions.
* Supported projections.
* Known limitations.
* Examples using real FITS-like headers.
* Differences from `WCS.jl`/`wcslib`, if any.

Documentation should include simple examples first.

Every mathematically nontrivial projection should include references to the FITS WCS papers.

Maintain a docs/dev/wcs_compliance_matrix.md file as you work, listing every FITS WCS feature, its status, reference section, Astropy/wcslib comparison status, and test coverage status.

---

## Astropy Gap Backlog

These tasks come from comparing FITSWCS.jl with Astropy 6.1.7 / WCSLIB 8.3
using `semble search` and focused `python3.10` probes.

### Task A: Expand Projection Coverage To WCSLIB `PRJ_CODES`

Astropy exposes the WCSLIB projection set `AZP`, `SZP`, `TAN`, `STG`, `SIN`,
`ARC`, `ZPN`, `ZEA`, `AIR`, `CYP`, `CEA`, `CAR`, `MER`, `SFL`, `PAR`, `MOL`,
`AIT`, `COP`, `COE`, `COD`, `COO`, `BON`, `PCO`, `TSC`, `CSC`, `QSC`, `HPX`,
and `XPH`.

Implementation notes:

* Add each projection as an explicit Julia type and parser mapping.
* Parse projection parameters from `PV<lat>_m` consistently with WCSLIB.
* Prefer closed-form inverses where the standard gives them; use bounded Newton
  solves with clear convergence failures for polynomial/implicit projections.
* Add Astropy/WCSLIB-generated absolute-value regression tests, not just
  round-trip tests.
* Preserve numeric precision through projection methods: use promotion from the
  numeric arguments and typed constants so `Float32` inputs can produce
  `Float32` outputs without doing the internal work in `Float64`.
* Keep high-risk projections such as `ZPN`, `AIR`, conics, quadcube, and
  HEALPix in small independently reviewed slices.

Progress notes:

* **Complete.**  All 28 WCSLIB projection codes (`AZP`, `SZP`, `TAN`, `STG`,
  `SIN`, `ARC`, `ZPN`, `ZEA`, `AIR`, `CYP`, `CEA`, `CAR`, `MER`, `SFL`,
  `PAR`, `MOL`, `AIT`, `COP`, `COE`, `COD`, `COO`, `BON`, `PCO`, `TSC`,
  `CSC`, `QSC`, `HPX`, `XPH`) are implemented and verified against Astropy
  regression tests to sub-microarcsecond precision.
  - `AZP` and `SZP` support only the default (central perspective) parameter
    forms; non-default PV parameters are rejected at parse time.
  - `CSC` agrees with Astropy to ~9 mas because WCSLIB stores the CSC
    polynomial coefficients as 32-bit `float`; our implementation is in 64-bit
    and is actually closer to the mathematical ideal.

### Task B: Implement Paper IV Distortion Lookup Tables

Astropy's full pipeline applies detector-to-image lookup corrections, SIP,
Paper IV lookup corrections, then the core WCS transform.

Implementation notes:

* Represent `CPDIS1/2`, `D2IMDIS1/2`, `DP*`, `DQ*`, and related lookup-table
  metadata explicitly instead of rejecting them.
* Add image/table HDU readers through existing FITSIO.jl and FITSFiles.jl
  extensions.
* Extend the public parsing API beyond header-only inputs.  Astropy requires
  `WCS(header, fobj=hdulist)` for Paper IV lookup distortions because the
  header keywords only identify auxiliary image extensions such as `WCSDVARR`
  and `D2IMARR`; the array payloads must be read from the FITS HDUList.  A
  FITSWCS API needs an equivalent way to pass the owning HDU list or a lookup
  resolver into `from_header`.
* Preserve the distinction Astropy exposes between the full pipeline and the
  core wcslib transform: `all_pix2world` / `all_world2pix` apply detector
  lookup tables, SIP, Paper IV lookup tables, and then the core WCS, while
  `wcs_pix2world` / `wcs_world2pix` operate on the core WCS only.  In Astropy,
  loaded lookup distortions are also visible as `cpdis1`, `cpdis2`, `det2im1`,
  and `det2im2` fields on the `WCS` object.
* Do not assume `WCS.jl` can supply these fixtures as a reference path.  Local
  probes against WCS.jl 0.6.3 / wcslib 7.7.0 failed on Astropy's
  `dist_lookup.fits.gz` SCI header with `Invalid parameter value`; WCS.jl's
  public `from_header(header::String)` API has no place to pass the HDUList
  needed to populate the lookup arrays.
* Implement bilinear interpolation for lookup-table offsets.
* Invert the complete distortion pipeline with a vectorized iterative solver
  and convergence diagnostics.

### Task C: Implement Paper III Spectral And Tabular Coordinates

FITSWCS currently supports plain linear spectral axes but rejects `-TAB` and
non-celestial algorithm-coded axes such as `FREQ-LOG`.

This is a substantial body of physics code orthogonal to the celestial
projection work that is the package's primary focus.  A realistic
implementation would need to cover:

**1. Spectral CTYPE parsing.**  The CTYPE field encodes the spectral
coordinate system and algorithm in the form `TTTT-AAA`, where `TTTT` is one
of `FREQ`, `ENER`, `WAVN`, `VRAD`, `VOPT`, `ZOPT`, `AWAV`, `VELO`, `BETA`,
and `AAA` is the algorithm code.  Algorithm codes defined in Paper III:

- ```` (empty): linear in the coordinate as given
- `LOG`: logarithmic (CRVAL is the reference value, CDELT is a scale factor
  rather than a linear increment)
- `TAB`: tabular lookup (see below)
- `F2W`, `W2F`: frequency–wavelength conversion via c = ν·λ

Each algorithm has a well-defined conversion from pixel to intermediate
world coordinate that would need to be implemented.

**2. `-TAB` table-lookup axes.**  A `-TAB` axis references an extension HDU
containing a coordinate array and an indexing vector.  The implementation
would need to:

* Parse `PS<axis>_0`, `PS<axis>_1`, `PV<axis>_0`, `PV<axis>_1` keywords
  that identify the table extension by number.
* Read the binary table extension (coordinate array column + optional
  indexing column) from the FITS file through the FITSIO.jl or FITSFiles.jl
  extensions.
* Add an API for supplying the owning HDU list or a table resolver to the WCS
  parser.  Astropy's API is `WCS(header, fobj=hdulist)`; header-only
  construction raises `ValueError` for `-TAB` because the WCS-TABLE extension
  contains the coordinate/index arrays.  The resulting Astropy `WCS` object
  stores table descriptors in `w.wcs.wtb` and the public transforms operate
  through `all_pix2world` / `all_world2pix`.
* Treat WCS.jl as an incomplete comparison source for this feature.  Although
  WCS.jl exposes `from_header(header; table=true)` and mirrors wcslib `tabprm`
  pointers internally, local probes with Astropy's `example_4d_tab.fits` and
  `tab-time-last-axis.fits` failed at setup with `Invalid parameter value`.
  The missing piece is the same as for Astropy header-only construction: no
  public API supplies the HDUList table contents.
* Implement linear interpolation between table entries for pixel→world and
  an iterative inversion (or use an inverse index array) for world→pixel.
* This couples the WCS parser to the file loader in a way the current
  dictionary-first design intentionally avoids.  Supporting `-TAB` would
  require either a breaking API change to `from_header` (to pass a file
  handle or HDU list) or a deferred resolution pattern where the table is
  resolved later.

**3. Velocity conversions.**  The `VRAD`, `VOPT`, `ZOPT`, `VELO`, `BETA`
coordinate types all represent recessional velocity / redshift, differing in
the relativistic convention used.  Converting between them and
frequency/wavelength requires:

* A rest frequency or rest wavelength (`RESTFRQ` or `RESTWAV` keyword).
* The velocity convention: radio (`VRAD`), optical (`VOPT`), or
  relativistic (`ZOPT`, `VELO`).
* The observer/target frame specification (`SPECSYS` keyword).
* Unit conversion between m/s, km/s, and fractional redshift.

The formulas are standard (Paper III §3–4) but require careful attention to
the frame convention — whether the velocity is defined in the observer's
rest frame, the source's rest frame, or a barycentric/LSR frame.

**4. Keyword set.**  Beyond the core linear WCS keywords, Paper III
introduces:

* `RESTFRQ` / `RESTWAV` — rest frequency/wavelength of the spectral line
* `SPECSYS` — spectral reference frame (e.g., `TOPOCENT`, `GEOCENTR`,
  `BARYCENT`, `HELIOCEN`, `LSRK`, `LSRD`, `GALACTOC`, `LOCALGRP`)
* `SSYSOBS` / `SSYSSRC` — observer/source spectral system (only for
  `VELO`-type axes)
* `VELOSYS` — systemic velocity offset
* `ZSOURCE` — source redshift

**5. Scope assessment.**  This is a significant feature (~500+ lines of new
code, mostly physics rather than WCS plumbing) that would roughly double the
non-projection complexity of the package.  It serves a different audience
than the celestial focus and would be better approached as a separate
funding/contributor-driven effort rather than bundled into the current
projection work.  The honest parse-time rejection of algorithm-coded axes
(already implemented) is the correct near-term posture.

### Task D: Add TPV/TPD And SCAMP Compatibility

Astropy handles `-TPV`, removes conflicting SIP when TPV is explicit, and
contains compatibility handling for older SCAMP headers that encode TPV terms
with `-TAN`.

Implementation notes:

* Parse TPV/TPD polynomial coefficients from PV keywords.
* Define precedence rules for TPV versus SIP.
* Add a compatibility mode for pre-2012 SCAMP `-TAN` + PV headers, gated by an
  explicit parser option if needed.

### Task E: Add Time And Stokes Semantics

Linear `TIME` and `STOKES` axes currently transform numerically but do not carry
physical interpretation.

Implementation notes:

* Parse FITS time metadata such as `MJDREF`, `DATEREF`, `DATE-OBS`, `TIMESYS`,
  `TIMEUNIT`, `TREFPOS`, and observatory position keywords.
* Represent Stokes axes as quantized polarization states rather than ordinary
  continuous linear coordinates.
* Add metadata accessors while keeping numeric transform APIs stable.

### Task F: Add APE 14-Style Axis Metadata

Astropy exposes world-axis physical types, names, units, high-level object
classes, pixel shapes/bounds, and axis-correlation matrices.

Implementation notes:

* Add axis classification for celestial, spectral, temporal, Stokes, tabular,
  and generic linear axes.
* Add `world_axis_physical_types`, `world_axis_units`, `world_axis_names`, and
  `axis_correlation_matrix` accessors.
* Compute correlations from the CD/PC matrix and celestial axis coupling; mark
  all axes correlated when full-image distortion is present.

### Task G: Add Header Serialization, Validation, And Fixups

Astropy/WCSLIB can serialize WCS objects, find all alternate WCS descriptions,
validate files, and apply common fixups.

Implementation notes:

* Implement `to_header`/`to_fits` equivalents for the supported keyword subset.
* Add `find_all_wcs` over primary and alternate WCS suffixes.
* Add validation and optional fixups analogous to `cdfix`, `unitfix`,
  `datfix`, `spcfix`, and `cylfix`, keeping automatic mutation opt-in.

### Task H: Add WCS Subsetting, Slicing, And Comparison

Astropy can extract separable sub-WCS objects, slice WCSes, compare WCS
metadata, and transform pixel coordinates between WCS objects.

Implementation notes:

* Implement separability checks based on the linear matrix and distortion
  presence.
* Add subsetting helpers for celestial, spectral, temporal, and Stokes axes.
* Add WCS comparison options for ancillary keywords, CRPIX shifts, and tiling.
* Build pixel-to-pixel utilities on top of `pixel_to_world` and
  `world_to_pixel`.

---

## Development Loop

Work in small cycles.

For each cycle:

1. Choose one feature slice.
2. Read the relevant standard section.
3. Inspect Astropy implementation and tests using `semble-search`.
4. Inspect `WCS.jl` behavior where relevant.
5. Write or update a short implementation note.
6. Implement the minimal feature cleanly.
7. Add correctness tests.
8. Run the full test suite.
9. Add focused benchmarks if the feature is performance-sensitive.
10. Refactor for simplicity.
11. Commit with a clear message.

Do not proceed to the next feature until the current feature has meaningful correctness tests.

---

## Definition of Done for a Feature

A feature is not done unless:

* It is implemented.
* It has unit tests.
* It has at least one correctness test.
* Its behavior is documented or clearly covered by existing docs.
* It handles malformed inputs reasonably.
* Its axis and unit conventions are clear.
* It does not introduce unnecessary dependencies.
* It does not significantly complicate unrelated code.
* The full test suite passes.

---

## Guardrails

Do not:

* Copy large blocks of Astropy code.
* Treat Astropy as the standard when it disagrees with FITS WCS papers.
* Add broad abstractions before multiple concrete use cases exist.
* Implement many projections without correctness tests.
* Rely only on round-trip tests; round-trip tests can hide symmetric mistakes.
* Hide FITS-vs-Julia indexing conventions.
* Add large dependencies for small utilities.
* Optimize before measuring.
* Make public API decisions casually.

Prefer:

* Small, explicit types.
* Clear math.
* Comments citing standard sections.
* Tests with known expected values.
* Feature slices that can be reviewed independently.
* Internal APIs that can evolve before being exposed publicly.

---

## Suggested First Task

Begin with a research and design pass.

1. Inspect `WCS.jl` public API and current test suite.
2. Use `semble-search` on Astropy to locate:

   * Header parsing code.
   * Core WCS object construction.
   * Linear transform tests.
   * TAN projection tests.
3. Read FITS WCS Paper I sections on linear transformations and coordinate keywords.
4. Write a short implementation note proposing:

   * Core internal types.
   * Pixel coordinate convention.
   * Header representation.
   * First milestone scope.
5. Implement Milestone 1 and Milestone 2 only:

   * Basic header parsing.
   * Linear pixel/world transforms.
   * Correctness tests.

Do not implement celestial projections until the linear WCS layer is tested and stable.
