# Plan: SIMD Vectorization of Batch Transform Trigonometry

## Problem

Batch `pixel_to_world` and `world_to_pixel` spend ~80% of their runtime in the
per-column projection + celestial-rotation loop.  Each column calls
`intermediate_to_native` (TAN: `hypot`, `atan`) and `native_to_celestial`
(three `sincos`, one `asin`, one `atan`).  These trig functions are compiled
to scalar x86-64 libm calls (`fsin`, `fcos`, `fpatan`, etc.) and LLVM cannot
auto-vectorise them because they cross `ccall` boundaries.

The linear part (`CD * pixels` via BLAS) takes <3% of total runtime for N=2
celestial WCS and is not a bottleneck.

Hoisting `wcs.projection` to a typed local, `@simd`, and `@inbounds` were
applied in a prior iteration.  They are harmless housekeeping but provide no
measurable speedup — the libm trig calls dominate.

## Why SLEEFPirates.jl alone does not work

[SLEEFPirates.jl](https://github.com/JuliaSIMD/SLEEFPirates.jl) provides
SLEEF-based SIMD `sin`, `cos`, `atan`, `asin`, `hypot`, and `sincos`.  It
extends `Base` trig for `AbstractSIMD` types (e.g. `Vec{4,Float64}`) from
VectorizationBase, but it does **not** override `Base.sin(::Float64)` or
`Base.sin(::Float32)`.  A `@simd` loop over scalar `Float64` values continues
to emit libm calls through `ccall`.  LLVM cannot auto-vectorise across
`ccall` boundaries, so the loop stays scalar regardless of annotations.

To get SIMD trig, the loop body must explicitly operate on `Vec{N,Float64}`
elements.  SLEEFPirates provides `SLEEFPirates.sin(::Vec{N,Float64})` which
compiles to SLEEF vector intrinsics.  But manually packing columns into
`Vec` elements and writing SIMD-aware loop logic is verbose and brittle —
that is exactly the problem LoopVectorization.jl solves.

## Approach: LoopVectorization.jl `@turbo`

[LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl)
provides the `@turbo` macro, which:

1. **Auto-packs scalars into SIMD vectors** — reads consecutive array elements
   into `Vec{N,T}` registers.
2. **Replaces trig calls with SLEEF** — detects `sin`, `cos`, `atan`, `asin`,
   `hypot`, `sincos` inside the loop body and lowers them to SLEEF vector
   intrinsics (via its dependency on VectorizationBase + SLEEFPirates).
3. **Handles predication** — branches inside the loop body (e.g. the TAN
   `iszero(r)` check, the world_to_pixel fiducial-point early return) are
   lowered to masked SIMD operations, keeping the vector pipeline full.
4. **Unrolls and schedules** — optimises loop unrolling and instruction
   scheduling for the target microarchitecture (AVX2, AVX-512, etc.).

The transformation is:

```
# Before:
for k in 1:N
    phi, theta = intermediate_to_native(proj, im[li,k], im[la,k])
    alpha, delta = native_to_celestial(phi, theta, ap, dp, pp)
    world[li,k] = mod(alpha * r2d, 360)
    world[la,k] = delta * r2d
end

# After:
@turbo for k in 1:N
    phi, theta = intermediate_to_native(proj, im[li,k], im[la,k])
    alpha, delta = native_to_celestial(phi, theta, ap, dp, pp)
    world[li,k] = mod(alpha * r2d, 360)
    world[la,k] = delta * r2d
end
```

`@turbo` inlines the projection and celestial functions into the loop body,
replaces all trig calls with SLEEF vector intrinsics, and lowers the result
to packed AVX2/AVX-512 instructions.  This works for **any projection** whose
forward/inverse is closed-form trig — no projection-specific inlining needed.

| Vectorisable | Not vectorisable (iterative inverse) |
|---|---|
| TAN, SIN, STG, ARC, ZEA | ZPN |
| CAR, CEA, CYP, MER, SFL | AIR |
| PAR, MOL, AIT, PCO | COP, COD, COE, COO |
| BON, TSC, CSC, QSC, HPX, XPH | TPV (when inverse is solved) |

For iterative projections, `@turbo` cannot vectorise the inner convergence
loop and falls back to scalar execution silently.

### Dependencies to add

LoopVectorization.jl pulls in several packages.  All are already present in
the Manifest from the SLEEFPirates experiment (they were its indirect
dependencies):

| Package | Role |
|---------|------|
| `LoopVectorization` | `@turbo` macro, loop restructuring |
| `VectorizationBase` | `Vec` SIMD types, low-level intrinsics |
| `SLEEFPirates` | SLEEF vector trig (sin, cos, atan, etc.) |
| `CPUSummary` | CPU feature detection (AVX2/AVX-512) |
| `LayoutPointers` | Strided array access patterns |
| `StaticArrayInterface` | Static array support |
| `ArrayInterface` | Array trait queries |
| `ManualMemory` | Manual memory management helpers |
| `SIMDTypes` | SIMD type definitions |
| `BitTwiddlingConvenienceFunctions` | Bit utilities |
| `HostCPUFeatures` | Runtime CPU feature query |
| `IfElse` | Predicated select operations |
| `Static` | Compile-time static values |

Only `LoopVectorization` itself needs to be added to `[deps]`; the rest are
pulled in automatically.  Compat bound: `LoopVectorization = "0.12"` (supports
Julia ≥ 1.10).

## Changes to `src/transforms.jl`

### 1. Add `using LoopVectorization` at module level

In `src/FITSWCS.jl`:
```julia
using LoopVectorization
```

### 2. Replace `@simd` with `@turbo` in batch projection loops

In `pixel_to_world(wcs, pixels::AbstractMatrix)`:

```julia
proj = wcs.projection
@turbo for k in 1:N
    phi, theta = intermediate_to_native(proj, intermediate[lon_idx, k], intermediate[lat_idx, k])
    alpha, delta = native_to_celestial(phi, theta, alpha_p, delta_p, phi_p)
    world[lon_idx, k] = mod(alpha * r2d, T(360))
    world[lat_idx, k] = delta * r2d
end
```

Same replacement in `world_to_pixel(wcs, worlds::AbstractMatrix)`:

```julia
proj = wcs.projection
@turbo for k in 1:N
    lon_delta = mod(T(worlds[lon_idx, k]) - T(wcs.crval[lon_idx]) + T(180), T(360)) - T(180)
    if abs(lon_delta) <= T(1e-10) && abs(T(worlds[lat_idx, k]) - T(wcs.crval[lat_idx])) <= T(1e-10)
        intermediate[lon_idx, k] = zero(T)
        intermediate[lat_idx, k] = zero(T)
    else
        alpha = T(worlds[lon_idx, k]) * d2r
        delta = T(worlds[lat_idx, k]) * d2r
        phi, theta = celestial_to_native(alpha, delta, alpha_p, delta_p, phi_p)
        x_lon, x_lat = native_to_intermediate(proj, phi, theta)
        intermediate[lon_idx, k] = T(x_lon)
        intermediate[lat_idx, k] = T(x_lat)
    end
end
```

The `@turbo` macro handles the `if`/`else` branch via masked predication.
The `mod` operation is also lowered to SLEEF vector intrinsics where
available.

### 3. `@inbounds` remains, `@simd` is replaced

`@turbo` implies `@inbounds` (it skips bounds checks automatically) and
replaces `@simd` (it is a superset).  The `@inbounds` annotations on other
batch helper loops (`_add_crval_rows!`, `_world_offsets_batch`, etc.) remain
— those loops do only arithmetic and already benefit from generic LLVM
auto-vectorisation.

## `@turbo` requirements and constraints

`@turbo` requires:
- Loop bounds known at entry (they are: `1:N` with `N` from `size(pixels, 2)`)
- Array arguments accessed with constant stride (they are: columns of a dense
  `Matrix{Float64}` with stride `N`)
- No function calls that cannot be inlined — `intermediate_to_native` and
  `native_to_celestial` dispatch on concrete projection types (thanks to the
  `proj` hoisting) and are inlinable for closed-form projections

Caveats:
- `@turbo` does not support `mod(x, y)` where `y` is non-integer in some
  versions.  `mod(alpha * r2d, T(360))` with integer `360` is fine.
- The `clamp` call inside `native_to_celestial` may not lower to SIMD in
  all cases; if it proves problematic, replace with
  `max(-one(T), min(one(T), x))`.
- First call to a `@turbo`-annotated function compiles target-specific code;
  warm-up is needed in benchmarks.

## Performance expectations

| Scenario | Current | With @turbo | Speedup |
|----------|---------|-------------|---------|
| TAN batch 100k coords | 7.8 ms | ~2–3 ms (est.) | 2.5–4× |
| SIN batch 100k coords | ~8 ms | ~2–3 ms (est.) | 2.5–4× |
| TAN scalar SVector | 35 ns | unchanged | 1× |
| ZPN batch (iterative) | ~100 ms | unchanged | 1× |

The speedup ceiling is ~4× on AVX2 (4-wide `Float64` SIMD) minus loop
overhead.  This brings batch throughput from ~78 ns/coord to ~20–30 ns/coord,
making the batch path faster per-element than the scalar `SVector` path
(~35 ns/coord).  AVX-512 (`Float64 × 8`) would give a further 1.5–2×.

## Tests

- All existing regression tests must pass — SLEEF guarantees <1 ULP accuracy.
- The existing "Batch transforms agree with scalar" tests (runtests.jl lines
  899 and 1340) already verify batch-vs-scalar identity at tight tolerances.
- No new consistency tests needed unless `@turbo` introduces numerical
  differences beyond floating-point noise.

## Verification

### 1. `@code_native` inspection

After implementation, inspect the generated code for the TAN batch loop:
```julia
@code_native debuginfo=:none batch_tan_loop!(world, intermediate, wcs, N, T)
```
Output should contain packed AVX instructions (`vaddpd`, `vmulpd`,
`vfnmadd213pd`) and SLEEF vector calls (e.g. `Sleef_atan2d4_u35avx2` for
4-wide `Float64` atan2 on AVX2).

### 2. Benchmark against scalar

```julia
K = 100_000
pix_mat = rand(2, K)
@btime pixel_to_world($wcs, $pix_mat)        # batch with @turbo
@btime pixel_to_world.(Ref($wcs), eachcol($pix_mat))  # broadcast
```

The batch should be 2.5–4× faster than broadcast.

### 3. Benchmark against pre-@turbo baseline

Run `benchmark/benchmarks.jl` before and after the change and compare
"TAN/batch-1M" and "TAN/batch-100" timings.

## Implementation order

1. Add `LoopVectorization` to `Project.toml` deps and compat.
2. Add `using LoopVectorization` to `src/FITSWCS.jl`.
3. Replace `@simd` with `@turbo` in both batch projection loops.
4. Run full test suite — all existing tests must pass.
5. Verify `@code_native` contains packed instructions and SLEEF vector calls.
6. Benchmark and compare against pre-`@turbo` baseline.
