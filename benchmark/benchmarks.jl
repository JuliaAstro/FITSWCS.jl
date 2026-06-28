"""
Benchmarks for FITSWCS.jl.

Organized by operation type. Each suite group covers one conceptual category
of operations. Run with AirSpeedVelocity or BenchmarkTools directly.

Usage (manual):
    julia --project=benchmark benchmark/benchmarks.jl

Usage (AirSpeedVelocity, from repo root):
    avs run

Groups:
    "pixel_to_world" – forward pixel → world transforms
    "world_to_pixel" – inverse world → pixel transforms
    "parsing"        – from_header construction
"""

using BenchmarkTools
using FITSWCS
using StaticArrays
using PrettyTables

function show_benchmarks(results)
    # Collect results — results may be a flat Dict or a nested BenchmarkGroup;
    # flatten first so that every value is a Trial.
    flat = flatten_results(results)
    # Collect results
    sorted = sort(collect(flat), by=first)
    names = [k for (k,_) in sorted]
    trials = [v for (_,v) in sorted]

    # Pack into matrix
    data = hcat(
        names,
        [BenchmarkTools.prettytime(median(t).time) for t in trials],
        [BenchmarkTools.prettymemory(median(t).memory) for t in trials],
        [median(t).allocs for t in trials]
    )

    # Make pretty table
    pretty_table(data;
        column_labels = ["Benchmark", "Median Time", "Memory", "Allocs"],
        alignment     = [:l, :r, :r, :r]
    )
end

function flatten_results(group)
    # Recursively flatten a (possibly nested) BenchmarkGroup of Trial results
    # into a flat Dict{String, Trial}.  Dict inputs are returned as-is so that
    # the function is idempotent.
    flat = Dict{String, Any}()
    _flatten_results!(flat, group, "")
    return flat
end

function _flatten_results!(flat, group, prefix)
    for (k, v) in group
        fullname = isempty(prefix) ? string(k) : "$prefix/$k"
        if v isa BenchmarkGroup
            _flatten_results!(flat, v, fullname)
        else
            flat[fullname] = v
        end
    end
end

#  ────────────────────────────────────────────────────────────────────────────

const SUITE = BenchmarkGroup()

# ── Shared WCS objects ───────────────────────────────────────────────────────

const _hdr_tan = Dict(
    "NAXIS"  => 2,
    "CTYPE1" => "RA---TAN",  "CTYPE2" => "DEC--TAN",
    "CRPIX1" => 512.0,       "CRPIX2" => 512.0,
    "CRVAL1" => 83.8221,     "CRVAL2" => -5.3911,
    "CDELT1" => -2.7778e-4,  "CDELT2" =>  2.7778e-4,
)
const WCS_TAN = from_header(_hdr_tan)

const _hdr_ait = Dict(
    "NAXIS"  => 2,
    "CTYPE1" => "GLON-AIT",  "CTYPE2" => "GLAT-AIT",
    "CRPIX1" => 360.5,       "CRPIX2" => 180.5,
    "CRVAL1" => 0.0,         "CRVAL2" => 0.0,
    "CDELT1" => -0.5,        "CDELT2" => 0.5,
)
const WCS_AIT = from_header(_hdr_ait)

const _hdr_sip = Dict(
    "NAXIS"   => 2,
    "CTYPE1"  => "RA---TAN-SIP",  "CTYPE2"  => "DEC--TAN-SIP",
    "CRPIX1"  => 512.0,            "CRPIX2"  => 512.0,
    "CRVAL1"  => 150.0,            "CRVAL2"  => 2.5,
    "CDELT1"  => -2.7778e-4,       "CDELT2"  =>  2.7778e-4,
    "A_ORDER" => 2,
    "A_2_0"   => 5.0e-6,           "A_0_2"   => 2.0e-6,  "A_1_1" => 0.0,
    "B_ORDER" => 2,
    "B_2_0"   => 1.0e-6,           "B_0_2"   => 0.0,     "B_1_1" => 3.0e-6,
)
const WCS_SIP = from_header(_hdr_sip)

const _hdr_cube = Dict(
    "NAXIS"  => 3,
    "CTYPE1" => "RA---TAN",  "CTYPE2" => "DEC--TAN",  "CTYPE3" => "FREQ",
    "CRPIX1" => 50.0,        "CRPIX2" => 50.0,         "CRPIX3" => 1.0,
    "CRVAL1" => 10.0,        "CRVAL2" => 25.0,          "CRVAL3" => 1.42e9,
    "CDELT1" => -0.01,       "CDELT2" =>  0.01,         "CDELT3" => 1.0e6,
)
const WCS_CUBE = from_header(_hdr_cube)

# Sample pixels
const _pix_tan   = [400.0, 300.0]
const _pix_ait   = [300.0, 150.0]
const _pix_sip   = [400.0, 300.0]
const _pix_cube  = [40.0, 60.0, 5.0]

const _world_tan  = pixel_to_world(WCS_TAN, _pix_tan)
const _world_ait  = pixel_to_world(WCS_AIT, _pix_ait)
const _world_sip  = pixel_to_world(WCS_SIP, _pix_sip)
const _world_cube = pixel_to_world(WCS_CUBE, _pix_cube)

# Batch of 100 pixels for TAN
const _batch_pix = [p .+ [i*0.5, i*0.3] for i in 1:100 for p in [_pix_tan]] |>
                   (x -> reduce(hcat, x))  # 2×100 matrix
const _batch_world = [pixel_to_world(WCS_TAN, _batch_pix[:,i]) for i in 1:100]

# Batch of 1M pixels for TAN
const _batch_pix_1M = [p .+ [i*0.5, i*0.3] for i in 1:1_000_000 for p in [_pix_tan]] |>
                      (x -> reduce(hcat, x))  # 2×1_000_000 matrix
const _batch_world_1M = pixel_to_world(WCS_TAN, _batch_pix_1M)

# ── pixel_to_world ───────────────────────────────────────────────────────────

SUITE["pixel_to_world"] = BenchmarkGroup()
let g = SUITE["pixel_to_world"]
    g["TAN/scalar"]     = @benchmarkable pixel_to_world($WCS_TAN, $_pix_tan)
    g["TAN/scalar/SVector Float64"] = @benchmarkable pixel_to_world($WCS_TAN, $(SVector{2,Float64}(_pix_tan)))
    g["TAN/scalar/SVector Float32"] = @benchmarkable pixel_to_world($WCS_TAN, $(SVector{2,Float32}(Float32.(_pix_tan))))
    g["TAN/scalar/Tuple"] = @benchmarkable pixel_to_world($WCS_TAN, $(Tuple(_pix_tan)))
    g["AIT/scalar"]     = @benchmarkable pixel_to_world($WCS_AIT, $_pix_ait)
    g["TAN-SIP/scalar"] = @benchmarkable pixel_to_world($WCS_SIP, $_pix_sip)
    g["3D-cube/scalar"] = @benchmarkable pixel_to_world($WCS_CUBE, $_pix_cube)
    g["TAN/batch-100"]  = @benchmarkable pixel_to_world($WCS_TAN, $_batch_pix)
    g["TAN/batch-1M"]  = @benchmarkable pixel_to_world($WCS_TAN, $_batch_pix_1M)
end

# ── world_to_pixel ───────────────────────────────────────────────────────────

SUITE["world_to_pixel"] = BenchmarkGroup()
let g = SUITE["world_to_pixel"]
    g["TAN/scalar"]     = @benchmarkable world_to_pixel($WCS_TAN, $_world_tan)
    g["AIT/scalar"]     = @benchmarkable world_to_pixel($WCS_AIT, $_world_ait)
    g["TAN-SIP/scalar"] = @benchmarkable world_to_pixel($WCS_SIP, $_world_sip)
    g["3D-cube/scalar"] = @benchmarkable world_to_pixel($WCS_CUBE, $_world_cube)
    g["TAN/batch-100"]  = @benchmarkable world_to_pixel($WCS_TAN, $_batch_world)
    g["TAN/batch-1M"]  = @benchmarkable world_to_pixel($WCS_TAN, $_batch_world_1M)
end

# ── parsing ──────────────────────────────────────────────────────────────────

SUITE["parsing"] = BenchmarkGroup()
let g = SUITE["parsing"]
    g["from_header/TAN"]     = @benchmarkable from_header($_hdr_tan)
    g["from_header/AIT"]     = @benchmarkable from_header($_hdr_ait)
    g["from_header/TAN-SIP"] = @benchmarkable from_header($_hdr_sip)
    g["from_header/3D-cube"] = @benchmarkable from_header($_hdr_cube)
end


# ── If not on CI, show a nice table ──────────────────────────────────────────
if get(ENV, "CI", "false") == "false"
    # Run the requested benchmarks and print a table for each suite.
    # run_selected_suites(ARGS)
    # Run the benchmarks
    for name in keys(SUITE)
        results = run(SUITE[name], verbose=true)
        println("\nBenchmark suite: $name")
        show_benchmarks(results)
    end
    # results = run(SUITE, verbose=true)
    # show_benchmarks(results)
end
