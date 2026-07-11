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
    "parsing"        – WCS construction
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
const WCS_TAN = WCS(_hdr_tan)

const _hdr_tan_preserved = Dict(
    "NAXIS"  => 2,
    "CTYPE1" => "RA---TAN",  "CTYPE2" => "DEC--TAN",
    "CRPIX1" => 512.0,       "CRPIX2" => 512.0,
    "CRVAL1" => 83.8221 * 3600,     "CRVAL2" => -5.3911 * 3600,
    "CDELT1" => -2.7778e-4 * 3600,  "CDELT2" =>  2.7778e-4 * 3600,
    "CUNIT1" => "arcsec",  "CUNIT2" => "arcsec",
)
const WCS_TAN_PRESERVED = WCS(_hdr_tan_preserved; preserve_units = true)

const WCS_TAN_SLICED = slice_wcs(WCS_TAN, 1:512, 1:512)

const _hdr_ait = Dict(
    "NAXIS"  => 2,
    "CTYPE1" => "GLON-AIT",  "CTYPE2" => "GLAT-AIT",
    "CRPIX1" => 360.5,       "CRPIX2" => 180.5,
    "CRVAL1" => 0.0,         "CRVAL2" => 0.0,
    "CDELT1" => -0.5,        "CDELT2" => 0.5,
)
const WCS_AIT = WCS(_hdr_ait)

const _hdr_sip = Dict(
    "NAXIS"   => 2,
    "CTYPE1"  => "RA---TAN-SIP",  "CTYPE2"  => "DEC--TAN-SIP",
    "CRPIX1"  => 512.0,           "CRPIX2"  => 512.0,
    "CRVAL1"  => 150.0,           "CRVAL2"  => 2.5,
    "CDELT1"  => -2.7778e-4,      "CDELT2"  =>  2.7778e-4,
    "A_ORDER" => 2,
    "A_2_0"   => 5.0e-6,          "A_0_2"   => 2.0e-6,  "A_1_1" => 0.0,
    "B_ORDER" => 2,
    "B_2_0"   => 1.0e-6,          "B_0_2"   => 0.0,     "B_1_1" => 3.0e-6,
)
const WCS_SIP = WCS(_hdr_sip)

const _hdr_sip_paperiv = merge(
    copy(_hdr_sip),
    Dict(
        "D2IMDIS1" => "LOOKUP", "D2IM1.AXIS.1" => 1,
        "D2IMDIS2" => "LOOKUP", "D2IM2.AXIS.2" => 2,
        "CPDIS1" => "LOOKUP", "DP1.AXIS.1" => 1,
        "CPDIS2" => "LOOKUP", "DP2.AXIS.2" => 2,
    ),
)

"""Synthetic in-memory auxiliary-data source for TAN-SIP + Paper IV benchmarks."""
struct BenchmarkPaperIVFobj end

const _d2im_x = FITSWCS.LookupTable2D(
    [1.0e-3 * (i - 1) + 2.0e-4 * (j - 1) for i in 1:16, j in 1:16];
    crpix = (1.0, 1.0),
    crval = (1.0, 1.0),
    cdelt = (64.0, 64.0),
)
const _d2im_y = FITSWCS.LookupTable2D(
    [-3.0e-4 * (i - 1) + 8.0e-4 * (j - 1) for i in 1:16, j in 1:16];
    crpix = (1.0, 1.0),
    crval = (1.0, 1.0),
    cdelt = (64.0, 64.0),
)
const _cpdis_x = FITSWCS.LookupTable2D(
    [7.0e-4 * sin(0.2 * i) + 2.0e-4 * cos(0.15 * j) for i in 1:16, j in 1:16];
    crpix = (1.0, 1.0),
    crval = (1.0, 1.0),
    cdelt = (64.0, 64.0),
)
const _cpdis_y = FITSWCS.LookupTable2D(
    [5.0e-4 * cos(0.18 * i) - 3.0e-4 * sin(0.12 * j) for i in 1:16, j in 1:16];
    crpix = (1.0, 1.0),
    crval = (1.0, 1.0),
    cdelt = (64.0, 64.0),
)

function FITSWCS._auxiliary_wcs_data(header::AbstractDict, ::BenchmarkPaperIVFobj; alt::Char = ' ', minerr::Real = 0.0)
    # Return prebuilt backend-neutral lookup tables so benchmarks avoid FITS I/O.
    return FITSWCS.AuxiliaryWCSData(
        det2im = (_d2im_x, _d2im_y),
        cpdis = (_cpdis_x, _cpdis_y),
    )
end

const WCS_SIP_PAPERIV = WCS(_hdr_sip_paperiv; fobj = BenchmarkPaperIVFobj())

const _hdr_cube = Dict(
    "NAXIS"  => 3,
    "CTYPE1" => "RA---TAN",  "CTYPE2" => "DEC--TAN",  "CTYPE3" => "FREQ",
    "CRPIX1" => 50.0,        "CRPIX2" => 50.0,         "CRPIX3" => 1.0,
    "CRVAL1" => 10.0,        "CRVAL2" => 25.0,          "CRVAL3" => 1.42e9,
    "CDELT1" => -0.01,       "CDELT2" =>  0.01,         "CDELT3" => 1.0e6,
)
const WCS_CUBE = WCS(_hdr_cube)

# ── Spectral cube with -TAB ──────────────────────────────────────────────────

"""Synthetic in-memory auxiliary-data source for 3D TAB spectral-cube benchmarks."""
struct BenchmarkTabularFobj end

function FITSWCS._auxiliary_wcs_data(header::AbstractDict, ::BenchmarkTabularFobj; alt::Char = ' ', minerr::Real = 0.0)
    # Return a prebuilt backend-neutral TAB payload so benchmarks avoid FITS I/O.
    return FITSWCS.AuxiliaryWCSData(
        tabular = FITSWCS._tabular_auxiliary_data(
            header,
            (extname, extver, extlev, column) -> begin
                column == "FREQS" && return Float64[10.0, 20.0, 40.0]
                throw(KeyError(column))
            end;
            alt = alt,
        ),
    )
end

const _hdr_cube_tab = Dict(
    "NAXIS"  => 3,
    "CTYPE1" => "RA---TAN",  "CTYPE2" => "DEC--TAN",  "CTYPE3" => "FREQ-TAB",
    "CRPIX1" => 512.0,       "CRPIX2" => 512.0,        "CRPIX3" => 1.0,
    "CRVAL1" => 83.8221,     "CRVAL2" => -5.3911,      "CRVAL3" => 1.0,
    "CD1_1"  => -2.7778e-4,  "CD1_2" => 5.5556e-5,    "CD1_3" => 0.0,
    "CD2_1"  => 5.5556e-5,   "CD2_2" => 2.7778e-4,    "CD2_3" => 0.0,
    "CD3_1"  => 0.0,         "CD3_2" => 0.0,           "CD3_3" => 1.0,
    "PS3_0"  => "WCS-TABLE",
    "PS3_1"  => "FREQS",
    "LONPOLE" => 180.0,
)
const WCS_CUBE_TAB = WCS(_hdr_cube_tab; fobj = BenchmarkTabularFobj())

# ── Spectral cube (FREQ-LOG) ─────────────────────────────────────────────────

const _hdr_cube_spec = Dict(
    "NAXIS"  => 3,
    "CTYPE1" => "RA---TAN",  "CTYPE2" => "DEC--TAN",  "CTYPE3" => "FREQ-LOG",
    "CRPIX1" => 512.0,       "CRPIX2" => 512.0,        "CRPIX3" => 1.0,
    "CRVAL1" => 83.8221,     "CRVAL2" => -5.3911,      "CRVAL3" => 1.42e9,
    "CD1_1"  => -2.7778e-4,  "CD1_2" => 5.5556e-5,    "CD1_3" => 0.0,
    "CD2_1"  => 5.5556e-5,   "CD2_2" => 2.7778e-4,    "CD2_3" => 0.0,
    "CD3_1"  => 0.0,         "CD3_2" => 0.0,           "CD3_3" => 1.0,
    "CUNIT3" => "Hz",
    "LONPOLE" => 180.0,
)
const WCS_CUBE_SPEC = WCS(_hdr_cube_spec)

# ── Grism (AWAV-GRA, KPNO Coude Feed) ───────────────────────────────────────

const _hdr_grism = Dict(
    "NAXIS"  => 1,
    "CTYPE1" => "AWAV-GRA",
    "CRPIX1" => 1801.7,
    "CRVAL1" => 5225.2,
    "CDELT1" => -0.4334,
    "CUNIT1" => "Angstrom",
    "PV1_0"  => 3.16e5,
    "PV1_1"  => 1.0,
    "PV1_2"  => 13.9,
)
const WCS_GRISM = WCS(_hdr_grism)

# ── 2D coupled TAB ───────────────────────────────────────────────────────────

"""Synthetic in-memory auxiliary-data source for 2D coupled-TAB benchmarks."""
struct BenchmarkTabular2DFobj end

function FITSWCS._auxiliary_wcs_data(header::AbstractDict, ::BenchmarkTabular2DFobj; alt::Char = ' ', minerr::Real = 0.0)
    coords = Array{Float64}(undef, 2, 2, 2)
    for k1 in 1:2, k2 in 1:2
        coords[1, k1, k2] = 100.0 + 10.0 * k1 + k2
        coords[2, k1, k2] = 200.0 + k1 + 10.0 * k2
    end
    return FITSWCS.AuxiliaryWCSData(
        tabular = FITSWCS._tabular_auxiliary_data(
            header,
            (extname, extver, extlev, column) -> begin
                column == "COORDS" && return coords
                column == "XINDEX" && return Float64[1.0, 2.0]
                column == "YINDEX" && return Float64[1.0, 2.0]
                throw(KeyError(column))
            end;
            alt = alt,
        ),
    )
end

const _hdr_coupled_tab = Dict(
    "NAXIS"  => 2,
    "CTYPE1" => "RA---TAB",  "CTYPE2" => "DEC--TAB",
    "CRPIX1" => 1.0,         "CRPIX2" => 1.0,
    "CRVAL1" => 1.0,         "CRVAL2" => 1.0,
    "CD1_1"  => 0.5,         "CD1_2" => 0.1,
    "CD2_1"  => -0.1,        "CD2_2" => 0.5,
    "PS1_0"  => "WCS-TABLE", "PS2_0" => "WCS-TABLE",
    "PS1_1"  => "COORDS",    "PS2_1" => "COORDS",
    "PS1_2"  => "XINDEX",    "PS2_2" => "YINDEX",
    "PV1_3"  => 1,           "PV2_3" => 2,
    "LONPOLE" => 180.0,
)
const WCS_COUPLED_TAB = WCS(_hdr_coupled_tab; fobj = BenchmarkTabular2DFobj())

# Sample pixels
const _pix_tan   = [400.0, 300.0]
const _pix_ait   = [300.0, 150.0]
const _pix_sip   = [400.0, 300.0]
const _pix_sip_paperiv = [400.0, 300.0]
const _pix_cube  = [40.0, 60.0, 5.0]
const _pix_cube_tab = [612.0, 412.0, 1.5]
const _pix_coupled_tab = [1.5, 1.5]
const _pix_cube_spec = [612.0, 412.0, 3.0]
const _pix_grism = [1900.0]

const _world_tan  = pixel_to_world(WCS_TAN, _pix_tan)
const _world_ait  = pixel_to_world(WCS_AIT, _pix_ait)
const _world_sip  = pixel_to_world(WCS_SIP, _pix_sip)
const _world_sip_paperiv = pixel_to_world(WCS_SIP_PAPERIV, _pix_sip_paperiv)
const _world_cube = pixel_to_world(WCS_CUBE, _pix_cube)
const _world_cube_tab = pixel_to_world(WCS_CUBE_TAB, _pix_cube_tab)
const _world_coupled_tab = pixel_to_world(WCS_COUPLED_TAB, _pix_coupled_tab)
const _world_cube_spec = pixel_to_world(WCS_CUBE_SPEC, _pix_cube_spec)
const _world_grism = pixel_to_world(WCS_GRISM, _pix_grism)

# Batch of 100 pixels for TAN
const _batch_pix = [p .+ [i*0.5, i*0.3] for i in 1:100 for p in [_pix_tan]] |>
                   (x -> reduce(hcat, x))  # 2×100 matrix
const _batch_world = pixel_to_world(WCS_TAN, _batch_pix)

# Batch of 1M pixels for TAN
const _batch_pix_1M = [p .+ [i*0.5, i*0.3] for i in 1:1_000_000 for p in [_pix_tan]] |>
                      (x -> reduce(hcat, x))  # 2×1_000_000 matrix
const _batch_world_1M = pixel_to_world(WCS_TAN, _batch_pix_1M)

# Batch of 100 pixels for 3D-TAB cube
const _batch_pix_cube_tab = reduce(hcat, [[612.0 + i*0.5, 412.0 + i*0.3, 1.0 + i*0.02] for i in 1:100])
const _batch_world_cube_tab = pixel_to_world(WCS_CUBE_TAB, _batch_pix_cube_tab)

# ── pixel_to_world ───────────────────────────────────────────────────────────

SUITE["pixel_to_world"] = BenchmarkGroup()
let g = SUITE["pixel_to_world"]
    g["TAN/scalar"] = @benchmarkable pixel_to_world($WCS_TAN, $_pix_tan) evals=100
    g["TAN/scalar/preserve_units"] = @benchmarkable pixel_to_world($WCS_TAN_PRESERVED, $_pix_tan) evals=100
    g["TAN/scalar/SVector Float64"] = @benchmarkable pixel_to_world($WCS_TAN, $(SVector{2,Float64}(_pix_tan))) evals=100
    g["TAN/scalar/SVector Float32"] = @benchmarkable pixel_to_world($WCS_TAN, $(SVector{2,Float32}(Float32.(_pix_tan)))) evals=100
    g["TAN/scalar/Tuple"] = @benchmarkable pixel_to_world($WCS_TAN, $(Tuple(_pix_tan))) evals=100
    g["TAN/scalar/sliced"] = @benchmarkable pixel_to_world($WCS_TAN_SLICED, $_pix_tan) evals=100
    g["AIT/scalar"] = @benchmarkable pixel_to_world($WCS_AIT, $_pix_ait) evals=100
    g["TAN-SIP/scalar"] = @benchmarkable pixel_to_world($WCS_SIP, $_pix_sip) evals=100
    g["TAN-SIP-PaperIV/scalar"] = @benchmarkable pixel_to_world($WCS_SIP_PAPERIV, $_pix_sip_paperiv) evals=100
    g["3D-cube/scalar"] = @benchmarkable pixel_to_world($WCS_CUBE, $_pix_cube) evals=100
    g["3D-cube-TAB/scalar"] = @benchmarkable pixel_to_world($WCS_CUBE_TAB, $_pix_cube_tab) evals=100
    g["3D-cube-spec/scalar"] = @benchmarkable pixel_to_world($WCS_CUBE_SPEC, $_pix_cube_spec) evals=100
    g["2D-coupled-TAB/scalar"] = @benchmarkable pixel_to_world($WCS_COUPLED_TAB, $_pix_coupled_tab) evals=100
    g["TAN/batch-100/Float64"] = @benchmarkable pixel_to_world($WCS_TAN, $_batch_pix) evals=1
    g["TAN/batch-100/Float32"] = @benchmarkable pixel_to_world($WCS_TAN, $(Float32.(_batch_pix))) evals=1
    g["TAN/batch-1M/Float64"] = @benchmarkable pixel_to_world($WCS_TAN, $_batch_pix_1M) evals=2 samples=3
    g["TAN/batch-1M/Float32"] = @benchmarkable pixel_to_world($WCS_TAN, $(Float32.(_batch_pix_1M))) evals=2 samples=3
    g["3D-cube-TAB/batch-100"] = @benchmarkable pixel_to_world($WCS_CUBE_TAB, $_batch_pix_cube_tab) evals=1
    g["grism/AWAV-GRA/scalar"] = @benchmarkable pixel_to_world($WCS_GRISM, $_pix_grism) evals=100
end

# ── world_to_pixel ───────────────────────────────────────────────────────────

SUITE["world_to_pixel"] = BenchmarkGroup()
let g = SUITE["world_to_pixel"]
    g["TAN/scalar"] = @benchmarkable world_to_pixel($WCS_TAN, $_world_tan) evals=100
    g["TAN/scalar/sliced"] = @benchmarkable world_to_pixel($WCS_TAN_SLICED, $_world_tan) evals=100
    g["TAN/scalar/preserve_units"] = @benchmarkable world_to_pixel($WCS_TAN_PRESERVED, $(_world_tan * 3600)) evals=100
    g["AIT/scalar"] = @benchmarkable world_to_pixel($WCS_AIT, $_world_ait) evals=100
    g["TAN-SIP/scalar"] = @benchmarkable world_to_pixel($WCS_SIP, $_world_sip) evals=100
    g["TAN-SIP-PaperIV/scalar"] = @benchmarkable world_to_pixel($WCS_SIP_PAPERIV, $_world_sip_paperiv) evals=100
    g["3D-cube/scalar"] = @benchmarkable world_to_pixel($WCS_CUBE, $_world_cube) evals=100
    g["3D-cube-TAB/scalar"] = @benchmarkable world_to_pixel($WCS_CUBE_TAB, $_world_cube_tab) evals=100
    g["3D-cube-spec/scalar"] = @benchmarkable world_to_pixel($WCS_CUBE_SPEC, $_world_cube_spec) evals=100
    g["2D-coupled-TAB/scalar"] = @benchmarkable world_to_pixel($WCS_COUPLED_TAB, $_world_coupled_tab) evals=100
    g["TAN/batch-100/Float64"] = @benchmarkable world_to_pixel($WCS_TAN, $_batch_world) evals=1
    g["TAN/batch-100/Float32"] = @benchmarkable world_to_pixel($WCS_TAN, $(Float32.(_batch_world))) evals=1
    g["TAN/batch-1M/Float64"] = @benchmarkable world_to_pixel($WCS_TAN, $_batch_world_1M) evals=2 samples=3
    g["TAN/batch-1M/Float32"] = @benchmarkable world_to_pixel($WCS_TAN, $(Float32.(_batch_world_1M))) evals=2 samples=3
    g["3D-cube-TAB/batch-100"] = @benchmarkable world_to_pixel($WCS_CUBE_TAB, $_batch_world_cube_tab) evals=1
    g["grism/AWAV-GRA/scalar"] = @benchmarkable world_to_pixel($WCS_GRISM, $_world_grism) evals=100
end

# ── parsing ──────────────────────────────────────────────────────────────────

SUITE["parsing"] = BenchmarkGroup()
let g = SUITE["parsing"]
    g["WCS/TAN"] = @benchmarkable WCS($_hdr_tan) evals=1 samples=10
    g["WCS/AIT"] = @benchmarkable WCS($_hdr_ait) evals=1 samples=10
    g["WCS/TAN-SIP"] = @benchmarkable WCS($_hdr_sip) evals=1 samples=10
    g["WCS/3D-cube"] = @benchmarkable WCS($_hdr_cube) evals=1 samples=10
    g["WCS/3D-cube-TAB"] = @benchmarkable WCS($_hdr_cube_tab; fobj = $(BenchmarkTabularFobj())) evals=5 samples=3
    g["WCS/3D-cube-spec"] = @benchmarkable WCS($_hdr_cube_spec) evals=5 samples=3
    g["WCS/grism/AWAV-GRA"] = @benchmarkable WCS($_hdr_grism) evals=5 samples=3
end

# ── slicing ──────────────────────────────────────────────────────────────────

SUITE["slicing"] = BenchmarkGroup()
let g = SUITE["slicing"]
    g["slice/2-D spatial"] = @benchmarkable slice_wcs($WCS_TAN, $(10:512), $(15:512)) evals=100
    g["slice/3-D drop spectral"] = @benchmarkable slice_wcs($WCS_CUBE, $(10:512), $(15:512), $5) evals=100
    g["slice/3-D drop spatial"] = @benchmarkable slice_wcs($WCS_CUBE, $10, $15) evals=100
    g["slice/2-D coupled-TAB"] = @benchmarkable slice_wcs($WCS_COUPLED_TAB, $(10:20), $(15:20)) evals=100
    g["slice/2-D spatial recursive"] = @benchmarkable slice_wcs($WCS_TAN_SLICED, $(10:20), $(15:30)) evals=100
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
