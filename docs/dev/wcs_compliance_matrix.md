# FITS WCS Compliance Matrix

This matrix tracks implementation status against the plan and the FITS WCS
paper families.  "Reference comparison" means checked against Astropy,
WCS.jl/wcslib, or another trusted implementation; hand-derived tests are noted
separately.

| Feature area | Status | Reference | Test coverage | Reference comparison |
| --- | --- | --- | --- | --- |
| Core axis count (`NAXIS`, `WCSAXES`) | Implemented | Paper I, Section 2.1 | Unit and error tests | Astropy/wcslib behavior reviewed |
| Core per-axis keywords (`CTYPE`, `CUNIT`, `CRPIX`, `CRVAL`, `CDELT`) | Implemented | Paper I, Section 2 | Unit and hand-computed tests | Astropy/wcslib behavior reviewed |
| Alternate WCS suffixes | Basic support | Paper I alternate coordinate descriptions | Primary/alternate parsing test | Astropy docstring behavior reviewed |
| PC/CD/CROTA linear matrices | Implemented | Paper I, Appendix A | Hand-computed and round-trip tests | Astropy/wcslib behavior reviewed |
| Celestial axis pairing and `CTYPE` parsing | Implemented | Paper II, Section 2 | Unit and error tests | Astropy/wcslib behavior reviewed |
| Celestial units | Basic support | Paper I units, Paper II degrees convention | Unit conversion, transform, and Astropy reference tests | **Astropy arcsec/rad reference values stored in test/regression_astropy_values.jl** |
| AZP and SZP default projections | Partial | Paper II perspective projections | Projection round-trip, parser error tests for unsupported PV parameters, randomized inverse, and Astropy reference tests | Default central-perspective forms match Astropy; non-default PV parameters are deferred with explicit errors |
| TAN projection | Implemented | Paper II, Eq. 54-55 | Projection and full WCS tests | **wcslib reference values stored in test/regression_wcslib.jl** |
| SIN projection including slant parameters | Implemented | Paper II, Eq. 48-49 | Projection and full WCS tests | **wcslib reference values stored in test/regression_wcslib.jl** |
| STG projection | Implemented | Paper II, Eq. 50 | Projection round-trip tests | **wcslib reference values stored in test/regression_wcslib.jl** |
| ARC projection | Implemented | Paper II, Eq. 46 | Projection round-trip and Astropy reference tests | **Astropy reference values stored in test/regression_astropy_values.jl** |
| ZEA projection | Implemented | Paper II, Eq. 52 | Projection round-trip and Astropy reference tests | **Astropy reference values stored in test/regression_astropy_values.jl** |
| CAR projection | Implemented | Paper II, Eq. 84 | Projection and full WCS tests | **wcslib reference values stored in test/regression_wcslib.jl** |
| CEA projection | Implemented | Paper II cylindrical equal-area projection | Projection, full WCS, and Astropy reference tests | **Astropy reference values stored in test/regression_astropy_values.jl** |
| CYP, MER, SFL, PAR, MOL, PCO projections | Implemented | Paper II cylindrical, pseudocylindrical, and polyconic projections | Projection round-trip, randomized inverse, and Astropy reference tests | **Astropy reference values stored in test/regression_astropy_values.jl** |
| AIT projection | Implemented | Paper II, Eq. 75 | Projection and full WCS tests | **wcslib reference values stored in test/regression_wcslib.jl** |
| Remaining WCSLIB projections (`ZPN`, `AIR`, `COP`, `COE`, `COD`, `COO`, `BON`, `TSC`, `CSC`, `QSC`, `HPX`, `XPH`) plus non-default AZP/SZP PV forms | Deferred | Paper II and WCSLIB projection set | Planned | Tracked in implementation_plan.md Astropy Gap Backlog |
| SIP distortion | Implemented | Shupe et al. 2005 SIP convention | Parser, polynomial, transform, and error tests | **wcslib reference values stored in test/regression_wcslib.jl** |
| FITSIO.jl extension | Implemented | Package API | Header and HDU tests | Local FITSIO API inspected |
| FITSFiles.jl extension | Implemented | Package API | Card vector and HDU tests | Local FITSFiles API inspected |
| Mixed celestial + spectral cubes | Basic linear/spherical support | Paper I linear axes; Paper II celestial pair | 3D split-axis and reference tests | **wcslib reference values stored in test/regression_wcslib.jl**; **Astropy split-axis values stored in test/regression_astropy_values.jl** |
| Spectral physical conversions | Deferred with explicit error for algorithm-coded axes | Paper III | Parser error tests plus linear spectral cube tests | Plain linear spectral axes work; non-linear algorithm CTYPEs such as `FREQ-LOG` are rejected instead of treated as linear |
| Time axes | Basic linear support | FITS time standard and Paper I linear axes | 4D mixed-axis and Astropy reference tests | Linear TIME axis checked against Astropy; no physical time conversion yet |
| Stokes axes | Basic linear support | FITS WCS conventions | 4D mixed-axis and Astropy reference tests | Linear STOKES axis checked against Astropy; no physical Stokes interpretation yet |
| Tabular lookup (`-TAB`) | Deferred with explicit error | Paper III | Parser error tests | Lookup-table transforms not implemented; `CTYPEia` with `TAB` algorithm is rejected instead of treated as linear |
| Distortion lookup tables | Deferred with explicit error | Paper IV | Parser error tests | CPDIS, D2IM, AXISCORR, DP, and DQ lookup-distortion keywords are rejected instead of ignored |
| Reference regression fixtures | **Implemented** | Instrument-style headers plus wcslib and Astropy | **106 tests in test/regression_wcslib.jl**; **122 tests in test/regression_astropy_values.jl**; optional `python3.10 test/regression_astropy.py` | wcslib (via WCS.jl from_header) for TAN-CDELT, TAN-PC-45deg, TAN-CD-HST, AIT, SIN, CAR, STG, 3D cube, TAN-SIP; Astropy 6.1.7 cross-check found 0 mismatches across 31 stored wcslib points; Astropy fixtures added for ARC, ZEA, CEA, AZP, SZP, CYP, MER, SFL, PAR, MOL, PCO, celestial CUNIT arcsec/rad, split celestial/spectral axes, and linear TIME/STOKES axes |
| API compatibility layer with WCS.jl names | Partial | WCS.jl public API | Header alias, keyword constructor, non-mutating alias, and mutating alias tests | FITS 1-based convention retained; constructor supports core vectors/matrices; full status/work-array API and specialized wcslib fields deferred |
| Type stability | **Implemented** | Julia `@inferred` macro | **13 `@inferred` tests in `test/runtests.jl` covering linear, TAN, SIP, and 3D transform paths** | `pixel_to_world` and `world_to_pixel` return `Vector{Float64}` (inferred), batch matrix returns `Matrix{Float64}` |
| Randomized property tests | **Implemented** | Fixed-seed MersenneTwister | **374 randomized tests in `test/runtests.jl`: linear WCS round-trips, PC-matrix CD consistency, 16 projection inverse round-trips, TAN celestial WCS round-trips** | Seed `0x5f3759df`; covers all implemented projections (AZP, SZP, TAN, SIN, STG, ARC, ZEA, CAR, CEA, CYP, MER, SFL, PAR, MOL, PCO, AIT) across 10 random configurations each |
| Benchmarks | **Implemented** | Plan performance section | **benchmark/benchmarks.jl (14 benchmarks, AirSpeedVelocity compatible)** | Baseline: TAN scalar ~245 ns / 560 B; AIT ~287 ns; SIP ~305 ns; parsing ~11 µs |
| User-facing documentation | Basic README | Plan documentation section | Manual review | Lists supported keywords, projections, conventions, and limitations |
