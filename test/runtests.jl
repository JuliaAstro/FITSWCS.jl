using Test
using LinearAlgebra
using Random
using FITSWCS
using FITSFiles
using FITSIO
using StaticArrays
using FITSWCS: pixel_to_intermediate, intermediate_to_pixel,
               parse_ctype, projection_from_code,
               native_to_celestial, celestial_to_native,
               compute_native_pole,
               unit_to_deg, build_cd_matrix,
               evaluate_sip_polynomial, sip_pixel_to_focal, sip_focal_to_pixel,
               NoDistortionPipeline, DistortionPipeline,
               pixel_to_focal, focal_to_pixel, has_distortion,
               LookupTable2D, interpolate_lookup_table,
               TabularAxisSpec, TabularWCSData, NoTabularWCSData,
               TabularWCSTable,
               NoAuxiliaryWCSData, AuxiliaryWCSData,
               has_auxiliary, _auxiliary_wcs_data,
               _tabular_forward

# Convenience shorthand
const D2R = π / 180.0
const R2D = 180.0 / π

"""Test that two angles (in radians) are equal modulo 2π."""
function angle_approx(a::Real, b::Real; atol=1e-12)
    d = mod(a - b + π, 2π) - π   # wrap difference to (-π, π]
    return abs(d) <= atol
end

"""Fake FITS container used to exercise auxiliary-data resolver dispatch."""
struct FakeAuxiliaryFobj end

function FITSWCS._auxiliary_wcs_data(header::AbstractDict, ::FakeAuxiliaryFobj; alt::Char=' ', minerr::Real=0.0)
    # Return a recognizable payload without requiring a real FITS backend.
    return AuxiliaryWCSData(det2im=:fake_det2im, cpdis=:fake_cpdis, tabular=:fake_tabular)
end

"""Fake FITS container that supplies concrete Paper IV lookup tables."""
struct FakeLookupFobj end

function FITSWCS._auxiliary_wcs_data(header::AbstractDict, ::FakeLookupFobj; alt::Char=' ', minerr::Real=0.0)
    # Make CPDIS1 depend on the D2IM-corrected y coordinate to test stage ordering.
    d2im_y = LookupTable2D(fill(0.5, 2, 2); crpix=(1, 1), crval=(1, 1), cdelt=(1, 1))
    cpdis_x = LookupTable2D([0.0 2.0; 0.0 2.0]; crpix=(1, 1), crval=(1, 1), cdelt=(1, 1))
    cpdis_y = LookupTable2D(fill(-0.25, 2, 2); crpix=(1, 1), crval=(1, 1), cdelt=(1, 1))

    # Return the same backend-independent payload that real FITS extensions produce.
    return AuxiliaryWCSData(det2im=(nothing, d2im_y), cpdis=(cpdis_x, cpdis_y))
end

"""Fake FITS container that supplies sparse CPDIS-only interpolation tables."""
struct FakeCPDISInterpolationFobj end

function FITSWCS._auxiliary_wcs_data(header::AbstractDict, ::FakeCPDISInterpolationFobj; alt::Char=' ', minerr::Real=0.0)
    x_data = zeros(Float32, 25, 25)
    y_data = zeros(Float32, 25, 25)

    # Place isolated offsets at known lookup-table coordinates.
    x_data[21, 11] = 0.5
    y_data[6, 11] = 0.7

    cpdis_x = LookupTable2D(x_data; crpix=(5, 10), crval=(10, 20), cdelt=(2, 2))
    cpdis_y = LookupTable2D(y_data; crpix=(5, 10), crval=(10, 20), cdelt=(3, 3))
    return AuxiliaryWCSData(cpdis=(cpdis_x, cpdis_y))
end

"""Fake FITS container that supplies a simple one-dimensional TAB table."""
struct FakeTabularFobj end

function FITSWCS._auxiliary_wcs_data(header::AbstractDict, ::FakeTabularFobj; alt::Char=' ', minerr::Real=0.0)
    # Load a backend-neutral TAB payload without requiring a FITS table reader.
    coords = Array{Float64}(undef, 2, 2, 2)
    for k1 in 1:2, k2 in 1:2
        coords[1, k1, k2] = 100.0 + 10.0 * k1 + k2
        coords[2, k1, k2] = 200.0 + k1 + 10.0 * k2
    end
    tabular = FITSWCS._tabular_auxiliary_data(
        header,
        (extname, extver, extlev, column) -> begin
            column == "FREQS" && return [10.0, 20.0, 40.0]
            column == "PIXELS" && return [1.0, 2.0, 3.0]
            column == "COORDS" && return coords
            column == "XINDEX" && return [1.0, 2.0]
            column == "YINDEX" && return [1.0, 2.0]
            throw(KeyError(column))
        end;
        alt = alt,
    )
    return AuxiliaryWCSData(tabular = tabular)
end

@testset "FITSWCS" begin

# ──────────────────────────────────────────────────────────────────────────────
@testset "CTYPE parsing" begin
    @testset "Standard celestial CTYPEs" begin
        sys, ct, pc = parse_ctype("RA---TAN")
        @test sys == "RA"
        @test ct  == :lon
        @test pc  == "TAN"

        sys, ct, pc = parse_ctype("DEC--TAN")
        @test sys == "DEC"
        @test ct  == :lat
        @test pc  == "TAN"

        sys, ct, pc = parse_ctype("GLON-AIT")
        @test sys == "GLON"
        @test ct  == :lon
        @test pc  == "AIT"

        sys, ct, pc = parse_ctype("GLAT-AIT")
        @test sys == "GLAT"
        @test ct  == :lat
        @test pc  == "AIT"
    end

    @testset "Linear / spectral CTYPEs" begin
        sys, ct, pc = parse_ctype("WAVE")
        @test ct == :linear
        @test pc == ""

        sys, ct, pc = parse_ctype("")
        @test ct == :linear
        @test pc == ""

        sys, ct, pc = parse_ctype("FREQ")
        @test ct == :linear
        @test pc == ""
    end

    @testset "Strict 4-3 CTYPE parsing" begin
        # Paper I treats non-4-3 CTYPE values as linear and suffixes as conventions.
        sys, ct, pc = parse_ctype("RA---TAN-SIP")
        @test sys == "RA"
        @test ct == :lon
        @test pc == "TAN"

        sys, ct, pc = parse_ctype("RA-TAN")
        @test sys == "RA-TAN"
        @test ct == :linear
        @test pc == ""

        sys, ct, pc = parse_ctype("FREQ-TAB")
        @test sys == "FREQ"
        @test ct == :linear
        @test pc == "TAB"

        sys, ct, pc = parse_ctype("FREQ-LOG")
        @test sys == "FREQ"
        @test ct == :linear
        @test pc == "LOG"
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Unit conversion" begin
    @test unit_to_deg("deg")    == 1.0
    @test unit_to_deg("")       == 1.0
    @test unit_to_deg("arcmin") ≈ 1.0/60.0
    @test unit_to_deg("arcsec") ≈ 1.0/3600.0
    @test unit_to_deg("rad")    ≈ R2D
    @test isnan(unit_to_deg("pixel"))
    # Case insensitivity
    @test unit_to_deg("DEG")    == 1.0
    @test unit_to_deg("ARCSEC") ≈ 1.0/3600.0
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Header parsing" begin

    @testset "Minimal 2D header with NAXIS" begin
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---TAN",
            "CTYPE2" => "DEC--TAN",
            "CRPIX1" => 100.0, "CRPIX2" => 200.0,
            "CRVAL1" => 10.0,  "CRVAL2" =>  20.0,
            "CDELT1" => -1e-4, "CDELT2" =>  1e-4,
        )
        wcs = WCS(hdr)
        @test wcs.naxis     == 2
        @test wcs.crpix     == [100.0, 200.0]
        @test wcs.crval     == [10.0, 20.0]
        @test wcs.lon_axis  == 1
        @test wcs.lat_axis  == 2
        @test wcs.projection isa TAN
    end

    @testset "CD matrix takes precedence over PC/CDELT" begin
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---TAN",
            "CTYPE2" => "DEC--TAN",
            "CRPIX1" => 1.0, "CRPIX2" => 1.0,
            "CRVAL1" => 0.0, "CRVAL2" => 0.0,
            "CDELT1" => 1.0, "CDELT2" => 1.0,   # would give identity
            "CD1_1"  => 2e-4, "CD1_2" => 0.0,
            "CD2_1"  => 0.0,  "CD2_2" => 2e-4,
        )
        wcs = WCS(hdr)
        @test wcs.cd[1, 1] ≈ 2e-4
        @test wcs.cd[2, 2] ≈ 2e-4
    end

    @testset "PC matrix + CDELT" begin
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---TAN",
            "CTYPE2" => "DEC--TAN",
            "CRPIX1" => 1.0, "CRPIX2" => 1.0,
            "CRVAL1" => 0.0, "CRVAL2" => 0.0,
            "CDELT1" => -2e-4, "CDELT2" => 2e-4,
            "PC1_1"  => 1.0, "PC1_2" => 0.0,
            "PC2_1"  => 0.0, "PC2_2" => 1.0,
        )
        wcs = WCS(hdr)
        @test wcs.cd[1, 1] ≈ -2e-4
        @test wcs.cd[2, 2] ≈  2e-4
        @test wcs.cd[1, 2] == 0.0
        @test wcs.cd[2, 1] == 0.0
    end

    @testset "SIN projection parameters parse from latitude-axis PV keywords" begin
        # Slant orthographic SIN stores xi/eta as PVi_1/PVi_2 on latitude axis i.
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---SIN",
            "CTYPE2" => "DEC--SIN",
            "PV2_1"  => 0.1,
            "PV2_2"  => -0.2,
        )
        wcs = WCS(hdr)
        @test wcs.projection == SIN(0.1, -0.2)
    end

    @testset "CEA projection parameter parses from latitude-axis PV keyword" begin
        # Cylindrical equal-area CEA stores lambda as PVi_1 on the latitude axis.
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---CEA",
            "CTYPE2" => "DEC--CEA",
            "PV2_1"  => 0.75,
        )
        wcs = WCS(hdr)
        @test wcs.projection == CEA(0.75)

        bad = copy(hdr)
        bad["PV2_1"] = 0.0
        @test_throws ArgumentError WCS(bad)
    end

    @testset "Error: PC and CD matrix forms cannot be mixed" begin
        # Paper I prohibits mixing PCi_ja and CDi_ja in one WCS description.
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "X",
            "CTYPE2" => "Y",
            "PC1_1"  => 1.0,
            "CD1_1"  => 1.0,
        )
        @test_throws ArgumentError WCS(hdr)
    end

    @testset "CROTA2 legacy rotation" begin
        # A 45° rotation with CDELT = 1 arcsec
        cdelt = 1.0 / 3600.0     # 1 arcsec in degrees
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---TAN",
            "CTYPE2" => "DEC--TAN",
            "CRPIX1" => 1.0, "CRPIX2" => 1.0,
            "CRVAL1" => 0.0, "CRVAL2" => 0.0,
            "CDELT1" => -cdelt, "CDELT2" => cdelt,
            "CROTA2" => 45.0,
        )
        wcs = WCS(hdr)
        rho = 45.0 * D2R
        @test wcs.cd[1, 1] ≈ -cdelt * cos(rho)   atol=1e-15
        @test wcs.cd[1, 2] ≈ -cdelt * sin(rho)   atol=1e-15
        @test wcs.cd[2, 1] ≈ -cdelt * sin(rho)   atol=1e-15
        @test wcs.cd[2, 2] ≈  cdelt * cos(rho)   atol=1e-15
    end

    @testset "Purely linear WCS (no celestial CTYPE)" begin
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "X",
            "CTYPE2" => "Y",
            "CRPIX1" => 1.0, "CRPIX2" => 1.0,
            "CRVAL1" => 0.0, "CRVAL2" => 0.0,
            "CDELT1" => 1.0, "CDELT2" => 1.0,
        )
        wcs = WCS(hdr)
        @test wcs.projection === nothing
        @test wcs.lon_axis   == 0
        @test wcs.lat_axis   == 0
    end

    @testset "1D WCS" begin
        hdr = Dict(
            "NAXIS"  => 1,
            "CTYPE1" => "FREQ",
            "CRPIX1" => 1.0,
            "CRVAL1" => 1.4e9,
            "CDELT1" => 1e6,
        )
        wcs = WCS(hdr)
        @test wcs.naxis == 1
        @test wcs.projection === nothing
    end

    @testset "WCSAXES takes precedence over NAXIS" begin
        hdr = Dict(
            "NAXIS"   => 3,
            "WCSAXES" => 2,
            "CTYPE1" => "RA---TAN",
            "CTYPE2" => "DEC--TAN",
            "CRPIX1" => 1.0, "CRPIX2" => 1.0,
            "CRVAL1" => 0.0, "CRVAL2" => 0.0,
            "CDELT1" => -1e-4, "CDELT2" => 1e-4,
        )
        wcs = WCS(hdr)
        @test wcs.naxis == 2
    end

    @testset "WCSAXES defaults to the largest WCS keyword index" begin
        # Paper I allows WCS dimensionality to exceed the image dimensionality.
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "X",
            "CTYPE2" => "Y",
            "CTYPE3" => "FREQ",
            "CRVAL3" => 1.4e9,
            "CDELT3" => 1e6,
        )
        wcs = WCS(hdr)
        @test wcs.naxis == 3
        @test wcs.crpix == [0.0, 0.0, 0.0]
    end

    @testset "Error: explicit WCSAXES bounds indexed keywords" begin
        # An explicit WCSAXES value makes higher-index WCS keywords malformed.
        hdr = Dict(
            "NAXIS"   => 2,
            "WCSAXES" => 2,
            "CTYPE1"  => "X",
            "CTYPE2"  => "Y",
            "CTYPE3"  => "FREQ",
        )
        @test_throws ArgumentError WCS(hdr)
    end

    @testset "Error: explicit WCSAXES bounds PV keywords" begin
        # Projection parameters are indexed WCS keywords and must fit WCSAXES.
        hdr = Dict(
            "NAXIS"   => 2,
            "WCSAXES" => 2,
            "CTYPE1"  => "RA---SIN",
            "CTYPE2"  => "DEC--SIN",
            "PV3_1"   => 0.1,
        )
        @test_throws ArgumentError WCS(hdr)
    end

    @testset "LONPOLE defaults" begin
        # zenithal (TAN): crval[lat] = 45 < theta0=90 → lonpole = 180
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---TAN",
            "CTYPE2" => "DEC--TAN",
            "CRPIX1" => 1.0, "CRPIX2" => 1.0,
            "CRVAL1" => 0.0, "CRVAL2" => 45.0,
            "CDELT1" => -1e-4, "CDELT2" => 1e-4,
        )
        wcs = WCS(hdr)
        @test wcs.lonpole == 180.0

        # CAR/AIT: crval[lat] = 0 = theta0=0 → lonpole = 180 (standard orientation)
        hdr_car = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---CAR",
            "CTYPE2" => "DEC--CAR",
            "CRPIX1" => 1.0, "CRPIX2" => 1.0,
            "CRVAL1" => 0.0, "CRVAL2" => 0.0,
            "CDELT1" => -1e-4, "CDELT2" => 1e-4,
        )
        wcs_car = WCS(hdr_car)
        @test wcs_car.lonpole == 180.0

        # CAR: crval[lat] = 30 > theta0=0 → lonpole = 0
        hdr_car2 = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---CAR",
            "CTYPE2" => "DEC--CAR",
            "CRPIX1" => 1.0, "CRPIX2" => 1.0,
            "CRVAL1" => 0.0, "CRVAL2" => 30.0,
            "CDELT1" => -1e-4, "CDELT2" => 1e-4,
        )
        wcs_car2 = WCS(hdr_car2)
        @test wcs_car2.lonpole == 0.0
    end

    @testset "Error: missing NAXIS/WCSAXES" begin
        @test_throws ArgumentError WCS(Dict("CRPIX1" => 1.0))
    end

    @testset "Alternate WCS keyword suffix selects independent solution" begin
        # Alternate descriptions should use suffixed keywords without changing primary WCS behavior.
        hdr = Dict(
            "NAXIS"   => 1,
            "CTYPE1"  => "FREQ",
            "CRPIX1"  => 1.0,
            "CRVAL1"  => 1.0e9,
            "CDELT1"  => 1.0e6,
            "CTYPE1A" => "WAVE",
            "CRPIX1A" => 5.0,
            "CRVAL1A" => 500.0,
            "CDELT1A" => 0.5,
        )

        primary = WCS(hdr)
        alternate = WCS(hdr; alt='A')

        @test pixel_to_world(primary, [2.0]) ≈ [1.001e9]
        @test pixel_to_world(alternate, [7.0]) ≈ [501.0]
        @test_throws ArgumentError WCS(hdr; alt='a')
    end

    @testset "Error: mismatched projection codes" begin
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---TAN",
            "CTYPE2" => "DEC--SIN",
            "CRPIX1" => 1.0, "CRPIX2" => 1.0,
            "CRVAL1" => 0.0, "CRVAL2" => 0.0,
            "CDELT1" => 1e-4, "CDELT2" => 1e-4,
        )
        @test_throws ArgumentError WCS(hdr)
    end

    @testset "TAB lookup axes require external data" begin
        # Paper III -TAB axes require table arrays, so header-only parsing must reject them.
        hdr = Dict(
            "NAXIS"  => 1,
            "CTYPE1" => "FREQ-TAB",
            "CRPIX1" => 1.0,
            "CRVAL1" => 1.0,
            "CDELT1" => 1.0,
            "PS1_0"  => "WCS-TABLE",
            "PS1_1"  => "FREQS",
        )
        @test_throws ArgumentError WCS(hdr)
        @test_throws ArgumentError WCS(hdr; fobj=:unsupported)

        alt_hdr = Dict(
            "NAXIS"   => 1,
            "CTYPE1"  => "FREQ",
            "CTYPE1A" => "WAVE-TAB",
            "PS1_0A"  => "WCS-TABLE",
            "PS1_1A"  => "FREQS",
        )
        @test pixel_to_world(WCS(alt_hdr), [2.0]) ≈ [2.0]
        @test_throws ArgumentError WCS(alt_hdr; alt='A')

        # TAB axis numbers must not exceed NAXIS.
        bad_axis_hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---TAN",  "CTYPE2" => "DEC--TAN",
            "CTYPE3" => "FREQ-TAB",
            "CRPIX1" => 1.0, "CRPIX2" => 1.0, "CRPIX3" => 1.0,
            "CRVAL1" => 0.0, "CRVAL2" => 0.0, "CRVAL3" => 1.0,
            "CDELT1" => 1.0, "CDELT2" => 1.0, "CDELT3" => 1.0,
            "PS3_0"  => "WCS-TABLE",
            "PS3_1"  => "FREQS",
        )
        @test_throws ArgumentError WCS(bad_axis_hdr)
    end

    @testset "Error: non-linear spectral algorithms are explicitly unsupported" begin
        # Paper III algorithms such as LOG and F2W are not equivalent to linear axes.
        hdr = Dict(
            "NAXIS"  => 1,
            "CTYPE1" => "FREQ-LOG",
            "CRPIX1" => 1.0,
            "CRVAL1" => 1.0e9,
            "CDELT1" => 1.0e6,
        )
        @test_throws ArgumentError WCS(hdr)

        alt_hdr = Dict(
            "NAXIS"   => 1,
            "CTYPE1"  => "FREQ",
            "CTYPE1A" => "WAVE-F2W",
            "CRPIX1A" => 1.0,
            "CRVAL1A" => 500.0,
            "CDELT1A" => 1.0,
        )
        @test pixel_to_world(WCS(alt_hdr), [2.0]) ≈ [2.0]
        @test_throws ArgumentError WCS(alt_hdr; alt='A')
    end

    @testset "Celestial units are converted without changing linear axes" begin
        # The public celestial API uses degrees while unrelated axes keep header units.
        hdr = Dict(
            "NAXIS"  => 3,
            "CTYPE1" => "RA---TAN",
            "CTYPE2" => "DEC--TAN",
            "CTYPE3" => "FREQ",
            "CRPIX1" => 0.0, "CRPIX2" => 0.0, "CRPIX3" => 0.0,
            "CRVAL1" => 0.0, "CRVAL2" => 0.0, "CRVAL3" => 1.0,
            "CDELT1" => 3600.0, "CDELT2" => 3600.0, "CDELT3" => 2.0,
            "CUNIT1" => "arcsec", "CUNIT2" => "arcsec", "CUNIT3" => "pixel",
        )
        wcs = WCS(hdr)
        @test wcs.cd[1, 1] ≈ 1.0
        @test wcs.cd[2, 2] ≈ 1.0
        @test wcs.cd[3, 3] ≈ 2.0
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "SIP distortion" begin

    @testset "Header parsing builds coefficient matrices" begin
        # SIP coefficient keywords should become triangular polynomial matrices.
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---TAN-SIP",
            "CTYPE2" => "DEC--TAN-SIP",
            "CRPIX1" => 10.0, "CRPIX2" => 20.0,
            "A_ORDER" => 2, "B_ORDER" => 2,
            "A_2_0" => 1e-3, "B_0_2" => -2e-3,
        )
        wcs = WCS(hdr)
        @test wcs.pipeline isa DistortionPipeline
        @test wcs.pipeline.sip isa SIPDistortion
        @test wcs.pipeline.sip.a[3, 1] == 1e-3
        @test wcs.pipeline.sip.b[1, 3] == -2e-3
        @test wcs.projection isa TAN
    end

    @testset "Distortion pipeline identity and SIP dispatch" begin
        # NoDistortionPipeline should be a no-op, while SIP pipelines correct axes 1–2.
        hdr_plain = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "X",
            "CTYPE2" => "Y",
            "CRPIX1" => 1.0, "CRPIX2" => 1.0,
            "CDELT1" => 1.0, "CDELT2" => 1.0,
        )
        plain = WCS(hdr_plain)
        @test plain.pipeline isa NoDistortionPipeline
        @test !has_distortion(plain.pipeline)
        @test pixel_to_focal(plain.pipeline, [2.0, 3.0], Val(2)) ≈ [2.0, 3.0]
        @test focal_to_pixel(plain.pipeline, [2.0, 3.0], Val(2)) ≈ [2.0, 3.0]
        static_pixel = SVector(2.0, 3.0)
        @test pixel_to_focal(plain.pipeline, static_pixel, Val(2)) === static_pixel
        @test focal_to_pixel(plain.pipeline, static_pixel, Val(2)) === static_pixel
        @test pixel_to_focal(plain.pipeline, [2.0, 3.0], Val(2)) isa SVector{2,Float64}

        hdr_sip = copy(hdr_plain)
        hdr_sip["A_ORDER"] = 2
        hdr_sip["B_ORDER"] = 2
        hdr_sip["A_2_0"] = 0.1
        sip_wcs = WCS(hdr_sip)
        @test sip_wcs.pipeline isa DistortionPipeline
        @test has_distortion(sip_wcs.pipeline)
        @test pixel_to_focal(sip_wcs.pipeline, [3.0, 3.0], Val(2)) ≈ [3.4, 3.0]
    end

    @testset "Polynomial evaluation uses CRPIX-relative powers" begin
        # The polynomial evaluator should sum terms up to total SIP order.
        coeff = zeros(3, 3)
        coeff[3, 1] = 0.5
        coeff[2, 2] = -0.25
        @test evaluate_sip_polynomial(coeff, 2.0, 3.0) ≈ 0.5 * 2.0^2 - 0.25 * 2.0 * 3.0
    end

    @testset "Forward correction feeds the linear transform" begin
        # A forward SIP offset should change intermediate coordinates before CD.
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "X",
            "CTYPE2" => "Y",
            "CRPIX1" => 10.0, "CRPIX2" => 20.0,
            "CRVAL1" => 0.0, "CRVAL2" => 0.0,
            "CDELT1" => 1.0, "CDELT2" => 1.0,
            "A_ORDER" => 2, "B_ORDER" => 2,
            "A_2_0" => 0.1, "B_0_2" => -0.2,
        )
        wcs = WCS(hdr)
        pix = [12.0, 23.0]
        focal = sip_pixel_to_focal(wcs.pipeline.sip, pix)
        @test focal ≈ [12.4, 21.2]
        @test pixel_to_intermediate(wcs, pix) ≈ [2.4, 1.2]
    end

    @testset "Inverse coefficients are used when present" begin
        # AP/BP coefficients provide the direct focal-to-pixel approximation.
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "X",
            "CTYPE2" => "Y",
            "CRPIX1" => 0.0, "CRPIX2" => 0.0,
            "A_ORDER" => 2, "B_ORDER" => 2,
            "A_2_0" => 0.0, "B_0_2" => 0.0,
            "AP_ORDER" => 2, "BP_ORDER" => 2,
            "AP_1_0" => -0.1, "BP_0_1" => 0.2,
        )
        wcs = WCS(hdr)
        @test sip_focal_to_pixel(wcs.pipeline.sip, [10.0, 5.0]) ≈ [9.0, 6.0]
    end

    @testset "Iterative inverse round-trips when AP/BP are absent" begin
        # Without inverse coefficients, the fixed-point solver should recover pixels.
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "X",
            "CTYPE2" => "Y",
            "CRPIX1" => 100.0, "CRPIX2" => 100.0,
            "A_ORDER" => 2, "B_ORDER" => 2,
            "A_2_0" => 1e-4, "B_0_2" => -2e-4,
        )
        wcs = WCS(hdr)
        pix = [103.0, 98.0]
        focal = sip_pixel_to_focal(wcs.pipeline.sip, pix)
        @test sip_focal_to_pixel(wcs.pipeline.sip, focal) ≈ pix atol=1e-8
    end

    @testset "SIP public transforms round-trip through linear WCS" begin
        # The high-level pixel/world API should include SIP in both directions.
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "X",
            "CTYPE2" => "Y",
            "CRPIX1" => 10.0, "CRPIX2" => 20.0,
            "CRVAL1" => 100.0, "CRVAL2" => 200.0,
            "CDELT1" => 2.0, "CDELT2" => 3.0,
            "A_ORDER" => 2, "B_ORDER" => 2,
            "A_2_0" => 1e-3, "B_0_2" => -1e-3,
        )
        wcs = WCS(hdr)
        pix = [12.0, 23.0]
        world = pixel_to_world(wcs, pix)
        @test world ≈ [104.008, 208.973]
        @test world_to_pixel(wcs, world) ≈ pix atol=1e-8
    end

    @testset "Malformed SIP headers throw clear errors" begin
        # SIP requires explicit CRPIX and matched forward/inverse order pairs.
        @test_throws ArgumentError WCS(Dict(
            "NAXIS" => 2,
            "A_ORDER" => 2,
            "B_ORDER" => 2,
        ))
        @test_throws ArgumentError WCS(Dict(
            "NAXIS" => 2,
            "CRPIX1" => 0.0, "CRPIX2" => 0.0,
            "A_ORDER" => 2,
        ))
        @test_throws ArgumentError WCS(Dict(
            "NAXIS" => 2,
            "CRPIX1" => 0.0, "CRPIX2" => 0.0,
            "A_ORDER" => 2, "B_ORDER" => 2,
            "AP_ORDER" => 2,
        ))
    end

    @testset "Unsupported lookup distortion keywords throw clear errors" begin
        # Paper IV lookup-table distortions must not be ignored as plain WCS.
        base = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---TAN",
            "CTYPE2" => "DEC--TAN",
            "CRPIX1" => 10.0, "CRPIX2" => 10.0,
            "CRVAL1" => 0.0,  "CRVAL2" => 0.0,
            "CDELT1" => -1e-3, "CDELT2" => 1e-3,
        )

        for (key, value) in [
            ("CPDIS1", "LOOKUP"),
            ("D2IMDIS1", "LOOKUP"),
            ("D2IMERR1", 0.01),
            ("AXISCORR", "OMIT"),
            ("DP1.NAXES", 2),
            ("DQ2.AXIS.1", 1),
        ]
            hdr = copy(base)
            hdr[key] = value
            @test_throws ArgumentError WCS(hdr)
        end

        # Alternate WCS lookup metadata should affect only the selected alternate.
        alt_base = copy(base)
        alt_base["CTYPE1A"] = "RA---TAN"
        alt_base["CTYPE2A"] = "DEC--TAN"
        for (key, value) in [
            ("CPDIS1A", "LOOKUP"),
            ("D2IMDIS1A", "LOOKUP"),
            ("D2IMERR1A", 0.01),
            ("AXISCORRA", "OMIT"),
            ("DP1A.NAXES", 2),
            ("DQ2A.AXIS.1", 1),
        ]
            hdr = copy(alt_base)
            hdr[key] = value
            @test WCS(hdr) isa WCSTransform
            @test_throws ArgumentError WCS(hdr; alt='A')
        end

        # SIP is implemented and should still parse through the distortion path.
        sip_hdr = copy(base)
        sip_hdr["CTYPE1"] = "RA---TAN-SIP"
        sip_hdr["CTYPE2"] = "DEC--TAN-SIP"
        sip_hdr["A_ORDER"] = 2
        sip_hdr["B_ORDER"] = 2
        @test WCS(sip_hdr).pipeline.sip isa SIPDistortion
    end

    @testset "SIP with 3D cube (RA/DEC/FREQ)" begin
        # SIP only corrects pixel axes 1–2; the spectral axis (3) must pass
        # through unchanged in both forward and inverse directions.
        hdr = Dict{String,Any}(
            "NAXIS"  => 3,
            "CTYPE1" => "RA---TAN-SIP", "CTYPE2" => "DEC--TAN-SIP", "CTYPE3" => "FREQ",
            "CRPIX1" => 256.0, "CRPIX2" => 256.0, "CRPIX3" => 1.0,
            "CRVAL1" => 83.8221, "CRVAL2" => -5.3911, "CRVAL3" => 1.42e9,
            "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4, "CDELT3" => 1.0e6,
            "A_ORDER" => 2, "B_ORDER" => 2,
            "A_0_1" => 1e-5, "B_1_0" => -2e-5,
        )
        wcs = WCS(hdr)
        @test wcs.naxis == 3
        @test wcs.pipeline.sip isa SIPDistortion

        # Forward: spectral axis passes through unchanged.
        pix = [300.0, 280.0, 5.0]  # off-reference in all axes
        world = pixel_to_world(wcs, pix)
        @test world[3] ≈ 1.42e9 + 1.0e6 * (5.0 - 1.0)  # linear spectral

        # Round-trip: spectral axis must be preserved.
        pix2 = world_to_pixel(wcs, world)
        @test pix2 ≈ pix  atol=1e-5
        @test pix2[3] ≈ pix[3]  # spectral axis unchanged by SIP
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Higher-dimensional WCS" begin
    @testset "Split celestial and spectral cube axes" begin
        # Mixed cubes can place a linear spectral axis between celestial longitude and latitude.
        hdr = Dict(
            "NAXIS"  => 3,
            "CTYPE1" => "RA---TAN", "CTYPE2" => "FREQ", "CTYPE3" => "DEC--TAN",
            "CRPIX1" => 30.0,       "CRPIX2" => 40.0,   "CRPIX3" => 45.0,
            "CRVAL1" => 10.0,       "CRVAL2" => 1.42e9, "CRVAL3" => 25.0,
            "CDELT1" => -0.01,      "CDELT2" => 1.0e6,  "CDELT3" => 0.01,
        )
        wcs = WCS(hdr)

        @test wcs.lon_axis == 1
        @test wcs.lat_axis == 3
        @test pixel_to_world(wcs, [30.0, 40.0, 45.0]) ≈ [10.0, 1.42e9, 25.0] atol=1e-12
        @test pixel_to_world(wcs, [30.0, 43.0, 45.0]) ≈ [10.0, 1.423e9, 25.0] atol=1e-12

        for pix in ([29.0, 43.0, 44.0], [31.5, 37.0, 47.0])
            world = pixel_to_world(wcs, pix)
            @test world[2] ≈ 1.42e9 + 1.0e6 * (pix[2] - 40.0)
            @test world_to_pixel(wcs, world) ≈ pix atol=1e-8
        end

        # Same split axis ordering with SIP distortion on pixel axes 1–2.
        # SIP always corrects the first two pixel axes; the CD matrix then
        # maps those corrections to world axes 1 (lon) and 3 (lat).
        hdr_sip = Dict{String,Any}(
            "NAXIS"  => 3,
            "CTYPE1" => "RA---TAN-SIP", "CTYPE2" => "FREQ", "CTYPE3" => "DEC--TAN-SIP",
            "CRPIX1" => 30.0,           "CRPIX2" => 40.0,    "CRPIX3" => 45.0,
            "CRVAL1" => 10.0,           "CRVAL2" => 1.42e9,  "CRVAL3" => 25.0,
            "CDELT1" => -0.01,          "CDELT2" => 1.0e6,   "CDELT3" => 0.01,
            "A_ORDER" => 2, "B_ORDER" => 2,
            "A_0_1" => 1e-5, "B_1_0" => -2e-5,
        )
        wcs_sip = WCS(hdr_sip)

        @test wcs_sip.lon_axis == 1
        @test wcs_sip.lat_axis == 3
        @test wcs_sip.pipeline.sip isa SIPDistortion

        # Spectral axis must pass through SIP unchanged (SIP only touches axes 1–2).
        pix = [30.0, 43.0, 45.0]
        world = pixel_to_world(wcs_sip, pix)
        @test world[2] ≈ 1.42e9 + 1.0e6 * (pix[2] - 40.0)

        # Full round-trip: split axes + SIP must recover the original pixel.
        # Tolerance is dominated by the SIP iterative inverse (1e-10 convergence
        # in focal-plane space, scaled by CD⁻¹ → ~1e-8 pixels; 1e-4 is safe).
        pix2 = world_to_pixel(wcs_sip, world)
        @test pix2 ≈ pix  atol=1e-4
    end

    @testset "Batch mixed-axis transforms agree with scalar transforms" begin
        # Batch transforms must preserve axis order for mixed celestial/non-celestial WCS.
        hdr = Dict(
            "NAXIS"  => 3,
            "CTYPE1" => "RA---TAN", "CTYPE2" => "FREQ", "CTYPE3" => "DEC--TAN",
            "CRPIX1" => 30.0,       "CRPIX2" => 40.0,   "CRPIX3" => 45.0,
            "CRVAL1" => 10.0,       "CRVAL2" => 1.42e9, "CRVAL3" => 25.0,
            "CDELT1" => -0.01,      "CDELT2" => 1.0e6,  "CDELT3" => 0.01,
        )
        wcs = WCS(hdr)
        pixels = [30.0 29.0 31.5;
                  40.0 43.0 37.0;
                  45.0 44.0 47.0]

        worlds = pixel_to_world(wcs, pixels)
        for k in axes(pixels, 2)
            @test worlds[:, k] ≈ pixel_to_world(wcs, pixels[:, k])
        end
        @test world_to_pixel(wcs, worlds) ≈ pixels atol=1e-8

        pixel_vectors = [pixels[:, k] for k in axes(pixels, 2)]
        world_vectors = pixel_to_world(wcs, pixel_vectors)
        @test world_vectors == [worlds[:, k] for k in axes(worlds, 2)]
        @test world_to_pixel(wcs, world_vectors) ≈ pixel_vectors atol=1e-8
    end

    @testset "Time and Stokes axes remain linear" begin
        # Paper I linear axes should work before adding physical time/Stokes interpretation.
        hdr = Dict(
            "NAXIS"  => 4,
            "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN",
            "CTYPE3" => "TIME",     "CTYPE4" => "STOKES",
            "CRPIX1" => 10.0,       "CRPIX2" => 20.0,
            "CRPIX3" => 1.0,        "CRPIX4" => 1.0,
            "CRVAL1" => 30.0,       "CRVAL2" => -5.0,
            "CRVAL3" => 59000.0,    "CRVAL4" => 1.0,
            "CDELT1" => -0.001,     "CDELT2" => 0.001,
            "CDELT3" => 0.5,        "CDELT4" => 1.0,
        )
        wcs = WCS(hdr)
        pix = [10.0, 20.0, 3.0, 4.0]
        world = pixel_to_world(wcs, pix)

        @test world[1:2] ≈ [30.0, -5.0] atol=1e-12
        @test world[3] ≈ 59001.0
        @test world[4] ≈ 4.0
        @test world_to_pixel(wcs, world) ≈ pix atol=1e-8
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "API compatibility helpers" begin
    # WCS.jl-style names should delegate to the canonical FITS 1-based API.
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "X", "CTYPE2" => "Y",
        "CRPIX1" => 1.0, "CRPIX2" => 1.0,
        "CRVAL1" => 0.0, "CRVAL2" => 0.0,
        "CDELT1" => 2.0, "CDELT2" => 3.0,
    )

    wcs = WCS(hdr)
    pixels = [2.0, 3.0]
    pixel_batch = [1.0 2.0;
                   1.0 3.0]

    @test wcs isa WCSTransform
    @test pix_to_world(wcs, pixels) ≈ pixel_to_world(wcs, pixels)
    @test pix_to_world(wcs, 2.0, 3.0) ≈ [2.0, 6.0]
    @test pix_to_world(wcs, 2.0, 3.0) isa SVector{2,Float64}
    @test pix_to_world(wcs, (2.0, 3.0)) isa SVector{2,Float64}
    @test world_to_pix(wcs, [2.0, 6.0]) ≈ [2.0, 3.0]
    @test world_to_pix(wcs, 2.0, 6.0) ≈ [2.0, 3.0]
    @test world_to_pix(wcs, 2.0, 6.0) isa SVector{2,Float64}
    @test world_to_pix(wcs, (2.0, 6.0)) isa SVector{2,Float64}
    @test pix_to_world(wcs, pixel_batch) ≈ pixel_to_world(wcs, pixel_batch)

    # Constructor compatibility should match the same transform built from a header.
    constructed = WCSTransform(2;
        ctype=["X", "Y"],
        crpix=[1.0, 1.0],
        crval=[0.0, 0.0],
        cdelt=[2.0, 3.0],
    )
    @test constructed isa WCSTransform
    @test pixel_to_world(constructed, pixels) ≈ [2.0, 6.0]

    rotated = WCSTransform(2;
        ctype=["X", "Y"],
        crpix=[1.0, 1.0],
        crval=[0.0, 0.0],
        pc=[0.0 -1.0; 1.0 0.0],
        cdelt=[2.0, 3.0],
    )
    @test pixel_to_world(rotated, [2.0, 1.0]) ≈ [0.0, 3.0]

    # Mutating WCS.jl-style aliases should fill caller-provided output arrays.
    world_out = similar(pixels)
    pixel_out = similar(pixels)
    batch_out = similar(pixel_batch)
    @test pix_to_world!(wcs, pixels, world_out) === world_out
    @test world_out ≈ [2.0, 6.0]
    @test world_to_pix!(wcs, world_out, pixel_out) === pixel_out
    @test pixel_out ≈ pixels
    @test pix_to_world!(wcs, pixel_batch, batch_out) === batch_out
    @test batch_out ≈ pixel_to_world(wcs, pixel_batch)

    @test_throws DimensionMismatch pix_to_world(wcs, 1.0)
    @test_throws DimensionMismatch world_to_pix(wcs, [1.0])
    @test_throws DimensionMismatch pix_to_world!(wcs, pixels, zeros(3))
    @test_throws DimensionMismatch world_to_pix!(wcs, world_out, zeros(3))
    @test_throws DimensionMismatch WCSTransform(2; crpix=[1.0])
    @test_throws DimensionMismatch WCSTransform(2; pc=ones(3, 3))
    @test_throws ArgumentError WCSTransform(2; restfrq=1.42e9)
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Auxiliary WCS data plumbing" begin
    # Auxiliary data should be resolved and stored without affecting ordinary headers.
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "X", "CTYPE2" => "Y",
        "CRPIX1" => 1.0, "CRPIX2" => 1.0,
        "CDELT1" => 1.0, "CDELT2" => 1.0,
    )

    plain = WCS(hdr)
    @test plain.aux isa NoAuxiliaryWCSData
    @test !has_auxiliary(plain.aux)

    generic_fobj = WCS(hdr; fobj=:unused)
    @test generic_fobj.aux isa NoAuxiliaryWCSData

    fake = WCS(hdr; fobj=FakeAuxiliaryFobj())
    @test fake.aux isa AuxiliaryWCSData
    @test has_auxiliary(fake.aux)
    @test fake.aux.det2im == :fake_det2im
    @test fake.aux.cpdis == :fake_cpdis
    @test fake.aux.tabular == :fake_tabular

    lookup_hdr = copy(hdr)
    lookup_hdr["CPDIS1"] = "LOOKUP"
    @test_throws ArgumentError WCS(lookup_hdr)
    @test_throws ArgumentError WCS(lookup_hdr; fobj=:unsupported)
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "TAB auxiliary data and transforms" begin
    # TAB axes should parse concrete specs and round-trip through interpolation.
    hdr = Dict(
        "NAXIS"  => 1,
        "CTYPE1" => "FREQ-TAB",
        "CRPIX1" => 1.0,
        "CRVAL1" => 1.0,
        "CDELT1" => 1.0,
        "PS1_0"  => "WCS-TABLE",
        "PS1_1"  => "FREQS",
    )
    indexed_hdr = copy(hdr)
    indexed_hdr["PS1_2"] = "PIXELS"

    specs = FITSWCS._tabular_axis_specs(hdr)
    indexed_specs = FITSWCS._tabular_axis_specs(indexed_hdr)
    @test only(specs) isa TabularAxisSpec{Nothing}
    @test only(indexed_specs) isa TabularAxisSpec{String}

    wcs = WCS(hdr; fobj=FakeTabularFobj())
    @test wcs.aux.tabular isa TabularWCSData
    @test pixel_to_world(wcs, [1.0]) ≈ [10.0]
    @test pixel_to_world(wcs, [1.5]) ≈ [15.0]
    @test pixel_to_world(wcs, [2.5]) ≈ [30.0]
    @test world_to_pixel(wcs, [10.0]) ≈ [1.0]
    @test world_to_pixel(wcs, [15.0]) ≈ [1.5]
    @test world_to_pixel(wcs, [30.0]) ≈ [2.5]

    pixels = reshape([1.0, 1.5, 2.5, 3.0], 1, 4)
    worlds = reshape([10.0, 15.0, 30.0, 40.0], 1, 4)
    @test pixel_to_world(wcs, pixels) ≈ worlds
    @test world_to_pixel(wcs, worlds) ≈ pixels

    indexed_wcs = WCS(indexed_hdr; fobj=FakeTabularFobj())
    @test pixel_to_world(indexed_wcs, [2.5]) ≈ [30.0]
    @test world_to_pixel(indexed_wcs, [30.0]) ≈ [2.5]

    # 2D coupled TAB with a non-trivial CD matrix that includes scaling
    # and rotation; the CD transform is applied before TAB lookup on both axes.
    coupled_hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---TAB", "CTYPE2" => "DEC--TAB",
        "CRPIX1" => 1.0, "CRPIX2" => 1.0,
        "CRVAL1" => 1.0, "CRVAL2" => 1.0,
        "CD1_1"  => 0.5, "CD1_2" => 0.1,
        "CD2_1"  => -0.1, "CD2_2" => 0.5,
        "PS1_0"  => "WCS-TABLE", "PS2_0" => "WCS-TABLE",
        "PS1_1"  => "COORDS", "PS2_1" => "COORDS",
        "PS1_2"  => "XINDEX", "PS2_2" => "YINDEX",
        "PV1_3"  => 1, "PV2_3" => 2,
        "LONPOLE" => 180.0,
    )
    coupled_wcs = WCS(coupled_hdr; fobj=FakeTabularFobj())
    c_fw = pixel_to_world(coupled_wcs, [1.5, 1.5])
    @test c_fw ≈ [114.2, 213.3]  atol=1e-10
    c_bw = world_to_pixel(coupled_wcs, SVector{2,Float64}(c_fw[1], c_fw[2]))
    @test c_bw ≈ [1.5, 1.5]  atol=1e-10

    coupled_pixels = [1.0 1.5 2.0; 1.0 1.5 2.0]
    coupled_worlds = [111.0 114.2 117.4; 211.0 213.3 215.6]
    @test pixel_to_world(coupled_wcs, coupled_pixels) ≈ coupled_worlds  atol=1e-10
    @test world_to_pixel(coupled_wcs, coupled_worlds) ≈ coupled_pixels  atol=1e-10

    # 2D coupled TAB with default 1-based indexing (no PS?_2 keywords).
    default_idx_hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---TAB", "CTYPE2" => "DEC--TAB",
        "CRPIX1" => 1.0, "CRPIX2" => 1.0,
        "CRVAL1" => 1.0, "CRVAL2" => 1.0,
        "CDELT1" => 1.0, "CDELT2" => 1.0,
        "PS1_0"  => "WCS-TABLE", "PS2_0" => "WCS-TABLE",
        "PS1_1"  => "COORDS", "PS2_1" => "COORDS",
        # No PS1_2 / PS2_2 — default 1-based indexing.
        "PV1_3"  => 1, "PV2_3" => 2,
        "LONPOLE" => 180.0,
    )
    default_wcs = WCS(default_idx_hdr; fobj=FakeTabularFobj())
    dfw = pixel_to_world(default_wcs, [1.5, 1.5])
    @test dfw ≈ [116.5, 216.5]  atol=1e-10
    @test world_to_pixel(default_wcs, SVector{2,Float64}(dfw[1], dfw[2])) ≈ [1.5, 1.5]  atol=1e-10

    # Mixed TAB + non-TAB: one ordinary linear axis and one TAB axis sharing a WCS.
    mixed_hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "LINEAR",
        "CTYPE2" => "FREQ-TAB",
        "CRPIX1" => 1.0, "CRPIX2" => 1.0,
        "CRVAL1" => 0.0, "CRVAL2" => 1.0,
        "CDELT1" => 1.0, "CDELT2" => 1.0,
        "PS2_0"  => "WCS-TABLE",
        "PS2_1"  => "FREQS",
    )
    mixed_wcs = WCS(mixed_hdr; fobj=FakeTabularFobj())
    @test pixel_to_world(mixed_wcs, [1.5, 1.0]) ≈ [0.5, 10.0]
    @test world_to_pixel(mixed_wcs, [0.5, 10.0]) ≈ [1.5, 1.0]

    # Batch path for mixed TAB + non-TAB.
    mixed_pixels = [1.0 1.5; 1.0 2.5]
    mixed_worlds = [0.0 0.5; 10.0 30.0]
    @test pixel_to_world(mixed_wcs, mixed_pixels) ≈ mixed_worlds
    @test world_to_pixel(mixed_wcs, mixed_worlds) ≈ mixed_pixels

    # Spectral cube: RA---TAN + DEC--TAN celestial projection with a FREQ-TAB
    # spectral axis.  This is the dominant real-world use case for -TAB.
    cube_hdr = Dict(
        "NAXIS"  => 3,
        "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CTYPE3" => "FREQ-TAB",
        "CRPIX1" => 512.0, "CRPIX2" => 512.0, "CRPIX3" => 1.0,
        "CRVAL1" => 83.8221, "CRVAL2" => -5.3911, "CRVAL3" => 1.0,
        "CD1_1"  => -2.7778e-4, "CD1_2" => 5.5556e-5, "CD1_3" => 0.0,
        "CD2_1"  => 5.5556e-5, "CD2_2" => 2.7778e-4, "CD2_3" => 0.0,
        "CD3_1"  => 0.0, "CD3_2" => 0.0, "CD3_3" => 1.0,
        "PS3_0"  => "WCS-TABLE",
        "PS3_1"  => "FREQS",
        "LONPOLE" => 180.0,
    )
    cube_wcs = WCS(cube_hdr; fobj=FakeTabularFobj())
    # At the reference pixel: spatial coords return CRVAL, spectral returns
    # the first coordinate-array entry.  Matches astropy all_pix2world exactly.
    @test pixel_to_world(cube_wcs, [512.0, 512.0, 1.0]) ≈ [83.8221, -5.3911, 10.0]
    # Spectral interpolation between coordinate-array entries.
    @test pixel_to_world(cube_wcs, [512.0, 512.0, 1.5]) ≈ [83.8221, -5.3911, 15.0]
    # Spatial offset with rotation in the CD matrix; the spectral axis is unchanged.
    fw_spatial = pixel_to_world(cube_wcs, [612.0, 512.0, 1.0])
    @test fw_spatial[1] ≈ 83.79419883756202  atol=1e-10
    @test fw_spatial[2] ≈ -5.385543765216007  atol=1e-10
    @test fw_spatial[3] ≈ 10.0
    # Round-trip at the reference point.
    cube_bw = world_to_pixel(cube_wcs, SVector{3,Float64}(83.8221, -5.3911, 10.0))
    @test cube_bw ≈ [512.0, 512.0, 1.0]
    # Round-trip with spectral interpolation.
    cube_bw2 = world_to_pixel(cube_wcs, SVector{3,Float64}(83.8221, -5.3911, 15.0))
    @test cube_bw2 ≈ [512.0, 512.0, 1.5]
    # Batch round-trip.
    batch_pix = [512.0 512.0 612.0; 512.0 512.0 512.0; 1.0 1.5 1.0]
    batch_world = pixel_to_world(cube_wcs, batch_pix)
    batch_round = world_to_pixel(cube_wcs, batch_world)
    @test batch_round ≈ batch_pix  atol=1e-10

    # Float32 input must be preserved through the TAB pipeline.
    @test eltype(pixel_to_world(wcs, Float32[1.5])) == Float32
    @test eltype(pixel_to_world(coupled_wcs, Float32[1.5, 1.5])) == Float32
    @test eltype(world_to_pixel(wcs, SVector{1,Float32}(15.0))) == Float32
    @test eltype(world_to_pixel(coupled_wcs, SVector{2,Float32}(116.5, 216.5))) == Float32
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "TAB high-dimensional fallback (M=5)" begin
    # The runtime-M construction path must produce a correctly typed table that
    # _tabular_forward can interpolate over.
    data = Array{Float64}(undef, 5, 2, 2, 2, 2, 2)
    for i1 in 1:2, i2 in 1:2, i3 in 1:2, i4 in 1:2, i5 in 1:2
        data[1, i1, i2, i3, i4, i5] = Float64(i1 + 2*i2 + 3*i3 + 4*i4 + 5*i5)
        for c in 2:5
            data[c, i1, i2, i3, i4, i5] = 0.0
        end
    end
    indices = [Float64[1.0, 2.0] for _ in 1:5]
    M = 5
    axes = SVector{M, Int}(i for i in 1:5)
    table_axes = SVector{M, Int}(i for i in 1:5)
    table = TabularWCSTable{M, typeof(data), typeof(indices)}(axes, table_axes, data, indices)
    @test typeof(table).parameters[1] == 5

    # Evaluate at the centre of the coordinate array.
    result = _tabular_forward(table, Float64[0.0, 0.0, 0.0, 0.0, 0.0],
                              Float64[1.0, 1.0, 1.0, 1.0, 1.0])
    @test result isa NTuple{5, Float64}
    @test result[1] ≈ 15.0
    # Components 2-5 are all zero in the test data.
    @test all(iszero, result[2:5])
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "LookupTable2D interpolation" begin
    # Paper IV lookup tables need predictable bilinear interpolation and clamping.
    data = Float64[0.0 10.0; 20.0 30.0]
    table = LookupTable2D(data; crpix=(1, 1), crval=(1, 1), cdelt=(1, 1))

    @test interpolate_lookup_table(table, 1.0, 1.0) == 0.0
    @test interpolate_lookup_table(table, 2.0, 1.0) == 20.0
    @test interpolate_lookup_table(table, SVector(1.5, 1.5)) == 15.0
    @test interpolate_lookup_table(table, 1.25, 1.75) ≈ 12.5
    @test interpolate_lookup_table(table, -100.0, 1.0) == 0.0
    @test interpolate_lookup_table(table, 100.0, 100.0) == 30.0
    @test table(1.25, 1.75) == interpolate_lookup_table(table, 1.25, 1.75)
    @test table(SVector(1.5, 1.5)) == interpolate_lookup_table(table, SVector(1.5, 1.5))

    metadata_table = LookupTable2D(
        reshape(Float64.(1:12), 3, 4);
        crpix=(2, 3),
        crval=(100, 200),
        cdelt=(10, -20),
    )
    @test interpolate_lookup_table(metadata_table, 100.0, 200.0) == metadata_table.data[2, 3]
    @test metadata_table(100.0, 200.0) == interpolate_lookup_table(metadata_table, 100.0, 200.0)

    typed_table = LookupTable2D(Float32[1 3; 5 7]; crpix=(1.0, 1.0), crval=(0.0, 0.0), cdelt=(2.0, 2.0))
    @test typed_table.crpix isa SVector{2,Float32}
    @test typed_table.crval isa SVector{2,Float32}
    @test typed_table.cdelt isa SVector{2,Float32}

    singleton_table = LookupTable2D(reshape(Float64[2.0, 6.0], 1, 2))
    @test interpolate_lookup_table(singleton_table, 100.0, 0.5) == 4.0

    @test_throws ArgumentError LookupTable2D(zeros(Float64, 0, 2))
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Linear transforms" begin

    @testset "Identity (CRPIX=1, CDELT=1, no rotation)" begin
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "X",
            "CTYPE2" => "Y",
            "CRPIX1" => 1.0, "CRPIX2" => 1.0,
            "CRVAL1" => 0.0, "CRVAL2" => 0.0,
            "CDELT1" => 1.0, "CDELT2" => 1.0,
        )
        wcs = WCS(hdr)
        # pixel [1,1] → world [0,0]
        @test pixel_to_world(wcs, [1.0, 1.0]) ≈ [0.0, 0.0]
        # world [0,0] → pixel [1,1]
        @test world_to_pixel(wcs, [0.0, 0.0]) ≈ [1.0, 1.0]
        # pixel [5,3] → world [4,2]
        @test pixel_to_world(wcs, [5.0, 3.0]) ≈ [4.0, 2.0]
        @test world_to_pixel(wcs, [4.0, 2.0]) ≈ [5.0, 3.0]
    end

    @testset "Paper I defaults make pixel values equal world values" begin
        # With default CRPIX=0, CRVAL=0, CDELT=1, coordinates follow pixel values.
        hdr = Dict("NAXIS" => 2)
        wcs = WCS(hdr)
        @test pixel_to_world(wcs, [1.0, 1.0]) ≈ [1.0, 1.0]
        @test world_to_pixel(wcs, [3.0, 4.0]) ≈ [3.0, 4.0]
    end

    @testset "Pure offset (CRPIX != 1)" begin
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "X",
            "CTYPE2" => "Y",
            "CRPIX1" => 10.0, "CRPIX2" => 20.0,
            "CRVAL1" => 100.0, "CRVAL2" => 200.0,
            "CDELT1" => 1.0, "CDELT2" => 1.0,
        )
        wcs = WCS(hdr)
        @test pixel_to_world(wcs, [10.0, 20.0]) ≈ [100.0, 200.0]
        @test world_to_pixel(wcs, [100.0, 200.0]) ≈ [10.0, 20.0]
        @test pixel_to_world(wcs, [11.0, 20.0]) ≈ [101.0, 200.0]
    end

    @testset "Pure scale (CDELT)" begin
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "X",
            "CTYPE2" => "Y",
            "CRPIX1" => 1.0, "CRPIX2" => 1.0,
            "CRVAL1" => 0.0, "CRVAL2" => 0.0,
            "CDELT1" => 0.5, "CDELT2" => 2.0,
        )
        wcs = WCS(hdr)
        @test pixel_to_world(wcs, [3.0, 2.0]) ≈ [1.0, 2.0]
        @test world_to_pixel(wcs, [1.0, 2.0]) ≈ [3.0, 2.0]
    end

    @testset "Rotation via PC matrix" begin
        # 90° rotation: x → -y, y → x
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "X",
            "CTYPE2" => "Y",
            "CRPIX1" => 0.0, "CRPIX2" => 0.0,
            "CRVAL1" => 0.0, "CRVAL2" => 0.0,
            "CDELT1" => 1.0, "CDELT2" => 1.0,
            "PC1_1"  => 0.0, "PC1_2" => -1.0,
            "PC2_1"  => 1.0, "PC2_2" =>  0.0,
        )
        wcs = WCS(hdr)
        # pixel [1,0] → world [0, 1]   (90° CCW rotation)
        w = pixel_to_world(wcs, [1.0, 0.0])
        @test w ≈ [0.0, 1.0]  atol=1e-14
        # round-trip
        @test world_to_pixel(wcs, w) ≈ [1.0, 0.0]  atol=1e-14
    end

    @testset "Round-trip random linear WCS" begin
        rng = (v -> v) ∘ identity   # deterministic via fixed data
        # Fixed random-ish WCS
        hdr = Dict(
            "NAXIS"  => 3,
            "CTYPE1" => "X", "CTYPE2" => "Y", "CTYPE3" => "Z",
            "CRPIX1" => 13.0, "CRPIX2" => 7.5, "CRPIX3" => 22.1,
            "CRVAL1" => 5.0,  "CRVAL2" => -3.0, "CRVAL3" => 0.1,
            "CDELT1" => 2.0,  "CDELT2" => -0.5, "CDELT3" => 0.01,
        )
        wcs = WCS(hdr)
        pixels = [1.0, 10.0, 50.0]
        world  = pixel_to_world(wcs, pixels)
        @test world_to_pixel(wcs, world) ≈ pixels  atol=1e-10
    end

    @testset "Batch transforms agree with scalar" begin
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "X", "CTYPE2" => "Y",
            "CRPIX1" => 5.0, "CRPIX2" => 5.0,
            "CRVAL1" => 0.0, "CRVAL2" => 0.0,
            "CDELT1" => 1.0, "CDELT2" => 1.0,
        )
        wcs = WCS(hdr)
        pix_mat = [1.0 5.0 10.0;
                   1.0 5.0 10.0]
        w_batch = pixel_to_world(wcs, pix_mat)
        for k in 1:3
            @test w_batch[:, k] ≈ pixel_to_world(wcs, pix_mat[:, k])
        end
        # Batch round-trip
        p_back = world_to_pixel(wcs, w_batch)
        @test p_back ≈ pix_mat  atol=1e-12
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Spherical rotation" begin
    # Identity rotation: alpha_p = 0, delta_p = 90°, phi_p = 0
    # → (phi, theta) should equal (alpha, delta)
    @testset "Identity rotation (pole at equatorial N pole)" begin
        alpha_p = 0.0 * D2R
        delta_p = 90.0 * D2R
        phi_p   = 0.0 * D2R
        alpha_in = 45.0 * D2R
        delta_in = 30.0 * D2R
        phi, theta = celestial_to_native(alpha_in, delta_in, alpha_p, delta_p, phi_p)
        alpha_out, delta_out = native_to_celestial(phi, theta, alpha_p, delta_p, phi_p)
        @test alpha_out ≈ alpha_in  atol=1e-12
        @test delta_out ≈ delta_in  atol=1e-12
    end

    @testset "Round-trip arbitrary rotation" begin
        alpha_p = 83.8221 * D2R
        delta_p = -5.3911 * D2R
        phi_p   = 180.0   * D2R
        for (a, d) in [(10.0, -5.0), (83.8, -5.4), (0.0, 0.0), (180.0, 45.0)]
            ai = a * D2R
            di = d * D2R
            phi, theta = celestial_to_native(ai, di, alpha_p, delta_p, phi_p)
            ao, do_ = native_to_celestial(phi, theta, alpha_p, delta_p, phi_p)
            @test ao ≈ ai  atol=1e-12
            @test do_ ≈ di  atol=1e-12
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "TAN projection" begin
    proj = TAN()

    @testset "Reference point (x=y=0)" begin
        # At (0,0) intermediate coords, native coords should be (0, 90°)
        phi, theta = intermediate_to_native(proj, 0.0, 0.0)
        @test isapprox(theta, π/2, atol=1e-12)
        # phi is undefined at the pole but atan(0,0) gives 0
    end

    @testset "Forward/inverse round-trip" begin
        # phi in radians, theta in radians, well away from singularity
        for (phi_d, theta_d) in [(0.0, 60.0), (45.0, 70.0), (-90.0, 80.0), (180.0, 45.0)]
            phi_r   = phi_d   * D2R
            theta_r = theta_d * D2R
            x, y    = native_to_intermediate(proj, phi_r, theta_r)
            phi2, theta2 = intermediate_to_native(proj, x, y)
            # phi is only defined mod 2π
            @test angle_approx(phi2, phi_r)
            @test theta2 ≈ theta_r  atol=1e-12
        end
    end

    @testset "TAN: singularity at theta <= 0 throws" begin
        @test_throws ErrorException native_to_intermediate(proj, 0.0, 0.0)
        @test_throws ErrorException native_to_intermediate(proj, 0.0, -0.1)
    end

    @testset "Hand-computed: x=0, y=-1 degree" begin
        # R_θ = 1, phi = 0, theta = atan(180/π, 1) ≈ atan(57.2958, 1)
        x = 0.0; y = -1.0
        phi, theta = intermediate_to_native(proj, x, y)
        @test phi ≈ 0.0  atol=1e-12
        expected_theta = atan(180.0 / π, 1.0)  # atan(R2D, R_theta_deg)
        @test theta ≈ expected_theta  atol=1e-12
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "SIN projection" begin
    proj = SIN()

    @testset "Forward/inverse round-trip" begin
        for (phi_d, theta_d) in [(0.0, 90.0), (45.0, 70.0), (-90.0, 45.0), (180.0, 60.0)]
            phi_r   = phi_d   * D2R
            theta_r = theta_d * D2R
            x, y    = native_to_intermediate(proj, phi_r, theta_r)
            phi2, theta2 = intermediate_to_native(proj, x, y)
            @test angle_approx(phi2, phi_r)
            @test theta2 ≈ theta_r  atol=1e-12
        end
    end

    @testset "SIN: singularity outside unit circle" begin
        # x=y=1 → R_θ² = 2 > 1
        @test_throws ErrorException intermediate_to_native(proj, 1.0/D2R, 1.0/D2R)
    end

    @testset "Slant parameters round-trip through inverse formula" begin
        # Nonzero xi/eta exercise the quadratic inverse and slant offset terms.
        slant = SIN(0.1, -0.05)
        phi_r = 20.0 * D2R
        theta_r = 70.0 * D2R
        x, y = native_to_intermediate(slant, phi_r, theta_r)
        phi2, theta2 = intermediate_to_native(slant, x, y)
        @test angle_approx(phi2, phi_r)
        @test theta2 ≈ theta_r atol=1e-12
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Full slant SIN celestial WCS" begin
    # Slant SIN PV parameters should participate in both transform directions.
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---SIN",
        "CTYPE2" => "DEC--SIN",
        "CRPIX1" => 64.0, "CRPIX2" => 64.0,
        "CRVAL1" => 10.0, "CRVAL2" => 20.0,
        "CDELT1" => -0.01, "CDELT2" => 0.01,
        "PV2_1"  => 0.1, "PV2_2" => -0.05,
    )
    wcs = WCS(hdr)

    @testset "Round-trip off reference pixel" begin
        # The slant terms make this distinct from the standard orthographic SIN path.
        pix_in = [70.0, 58.0]
        world = pixel_to_world(wcs, pix_in)
        pix_out = world_to_pixel(wcs, world)
        @test pix_out ≈ pix_in atol=1e-8
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "STG projection" begin
    proj = STG()
    @testset "Forward/inverse round-trip" begin
        for (phi_d, theta_d) in [(0.0, 45.0), (60.0, 30.0), (-90.0, 70.0)]
            phi_r   = phi_d   * D2R
            theta_r = theta_d * D2R
            x, y    = native_to_intermediate(proj, phi_r, theta_r)
            phi2, theta2 = intermediate_to_native(proj, x, y)
            @test angle_approx(phi2, phi_r)
            @test theta2 ≈ theta_r  atol=1e-12
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "ARC projection" begin
    proj = ARC()
    @testset "Forward/inverse round-trip" begin
        for (phi_d, theta_d) in [(0.0, 90.0), (30.0, 60.0), (-45.0, 80.0), (180.0, 30.0)]
            phi_r   = phi_d   * D2R
            theta_r = theta_d * D2R
            x, y    = native_to_intermediate(proj, phi_r, theta_r)
            phi2, theta2 = intermediate_to_native(proj, x, y)
            @test angle_approx(phi2, phi_r)
            @test theta2 ≈ theta_r  atol=1e-12
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "ZEA projection" begin
    proj = ZEA()
    @testset "Forward/inverse round-trip" begin
        for (phi_d, theta_d) in [(0.0, 90.0), (30.0, 60.0), (-45.0, 45.0), (180.0, 20.0)]
            phi_r   = phi_d   * D2R
            theta_r = theta_d * D2R
            x, y    = native_to_intermediate(proj, phi_r, theta_r)
            phi2, theta2 = intermediate_to_native(proj, x, y)
            @test angle_approx(phi2, phi_r)
            @test theta2 ≈ theta_r  atol=1e-12
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "CAR projection" begin
    proj = CAR()
    @testset "Forward/inverse round-trip" begin
        for (phi_d, theta_d) in [(0.0, 0.0), (45.0, 30.0), (-90.0, -45.0)]
            phi_r   = phi_d   * D2R
            theta_r = theta_d * D2R
            x, y    = native_to_intermediate(proj, phi_r, theta_r)
            phi2, theta2 = intermediate_to_native(proj, x, y)
            @test phi2   ≈ phi_r    atol=1e-12
            @test theta2 ≈ theta_r  atol=1e-12
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "CEA projection" begin
    proj = CEA(0.75)

    @testset "Forward equations follow equal-area formula" begin
        # CEA uses x=phi and y=sin(theta)/lambda in angular-degree units.
        x, y = native_to_intermediate(proj, 30.0 * D2R, 45.0 * D2R)
        @test x ≈ 30.0
        @test y ≈ R2D * sin(45.0 * D2R) / 0.75
    end

    @testset "Forward/inverse round-trip" begin
        # Round trips should recover native coordinates throughout the valid domain.
        for (phi_d, theta_d) in [(0.0, 0.0), (45.0, 30.0), (-90.0, -45.0), (120.0, 60.0)]
            phi_r = phi_d * D2R
            theta_r = theta_d * D2R
            x, y = native_to_intermediate(proj, phi_r, theta_r)
            phi2, theta2 = intermediate_to_native(proj, x, y)
            @test phi2 ≈ phi_r atol=1e-12
            @test theta2 ≈ theta_r atol=1e-12
        end
    end

    @testset "Inverse rejects outside latitude domain" begin
        # The inverse sine argument must stay within [-1, 1].
        @test_throws ErrorException intermediate_to_native(proj, 0.0, R2D / 0.75 + 1.0)
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Expanded cylindrical and pseudocylindrical projections" begin
    # These projections fill in a WCSLIB-compatible slice with formulas checked
    # by Astropy reference fixtures in regression_astropy_values.jl.
    for (proj, theta_range) in [
        (CYP(), (-70.0, 70.0)),
        (MER(), (-70.0, 70.0)),
        (SFL(), (-70.0, 70.0)),
        (PAR(), (-70.0, 70.0)),
        (MOL(), (-70.0, 70.0)),
    ]
        for (phi_d, theta_d) in [(0.0, 0.0), (30.0, 20.0), (-45.0, -25.0), (120.0, 50.0)]
            theta_range[1] <= theta_d <= theta_range[2] || continue
            phi_r = phi_d * D2R
            theta_r = theta_d * D2R
            x, y = native_to_intermediate(proj, phi_r, theta_r)
            phi2, theta2 = intermediate_to_native(proj, x, y)
            @test angle_approx(phi2, phi_r; atol=1e-11)
            @test theta2 ≈ theta_r atol=1e-11
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Additional WCSLIB projection defaults" begin
    # AZP and SZP default to central perspective; PCO exercises the WCSLIB inverse.
    for (proj, samples) in [
        (AZP(), [(0.0, 80.0), (30.0, 60.0), (-45.0, 70.0)]),
        (SZP(), [(0.0, 80.0), (30.0, 60.0), (-45.0, 70.0)]),
        (PCO(), [(0.0, 80.0), (30.0, 60.0), (-45.0, -20.0), (120.0, 30.0),
                 (155.95791089690434, 39.71570316255433)]),
    ]
        for (phi_d, theta_d) in samples
            phi_r = phi_d * D2R
            theta_r = theta_d * D2R
            x, y = native_to_intermediate(proj, phi_r, theta_r)
            phi2, theta2 = intermediate_to_native(proj, x, y)
            @test angle_approx(phi2, phi_r; atol=1e-11)
            @test theta2 ≈ theta_r atol=1e-11
        end
    end

    # Non-default AZP/SZP parameters parse and roundtrip correctly.
    azp_wcs = WCS(Dict(
        "NAXIS" => 2, "CTYPE1" => "RA---AZP", "CTYPE2" => "DEC--AZP",
        "CRPIX1" => 128.0, "CRPIX2" => 96.0, "CRVAL1" => 120.0, "CRVAL2" => 35.0,
        "CDELT1" => -0.05, "CDELT2" => 0.05, "PV2_1" => 1.0,
    ))
    @test azp_wcs.projection == AZP(1.0, 0.0)
    szp_wcs = WCS(Dict(
        "NAXIS" => 2, "CTYPE1" => "RA---SZP", "CTYPE2" => "DEC--SZP",
        "CRPIX1" => 128.0, "CRPIX2" => 96.0, "CRVAL1" => 120.0, "CRVAL2" => 35.0,
        "CDELT1" => -0.05, "CDELT2" => 0.05, "PV2_3" => 45.0,
    ))
    @test szp_wcs.projection == SZP(0.0, 0.0, 45.0)
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "WCSLIB projection edge semantics" begin
    # WCSLIB defines deterministic native coordinates at otherwise ambiguous
    # projection centers and boundaries; these cases guard those conventions.
    for proj in (AZP(), SZP(), TAN(), STG(), ARC(), ZEA())
        phi, theta = intermediate_to_native(proj, 0.0, 0.0)
        @test phi == 0.0
        @test theta ≈ π/2 atol=1e-12
    end

    @testset "SIN center longitude" begin
        # Standard and slant SIN both use phi=0 at the projection center.
        for proj in (SIN(), SIN(0.1, -0.05))
            phi, theta = intermediate_to_native(proj, 0.0, 0.0)
            @test phi == 0.0
            @test theta ≈ π/2 atol=1e-12
        end

        phi32, theta32 = intermediate_to_native(SIN(), Float32(0), Float32(0))
        @test phi32 isa Float32
        @test theta32 isa Float32
    end

    @testset "PAR projected pole" begin
        # PAR accepts the projected pole only when x is effectively zero.
        phi, theta = intermediate_to_native(PAR(), 0.0, 90.0)
        @test phi == 0.0
        @test theta ≈ π/2 atol=1e-12
        @test_throws ErrorException intermediate_to_native(PAR(), 1e-8, 90.0)
    end

    @testset "MOL projected pole" begin
        # MOL also treats nonzero x at the projected pole as outside the valid map.
        ypole = sqrt(2.0) * R2D
        phi, theta = intermediate_to_native(MOL(), 0.0, ypole)
        @test phi == 0.0
        @test theta ≈ π/2 atol=1e-12
        @test_throws ErrorException intermediate_to_native(MOL(), 1e-8, ypole)
    end

    @testset "Boundary tolerance" begin
        # WCSLIB allows tiny floating-point overshoot at closed projection boundaries.
        _, zea_theta = intermediate_to_native(ZEA(), 2.0 * R2D * (1.0 + 1e-13), 0.0)
        @test zea_theta ≈ -π/2 atol=1e-12

        _, cea_theta = intermediate_to_native(CEA(0.75), 0.0, (R2D / 0.75) * (1.0 + 5e-14))
        @test cea_theta ≈ π/2 atol=1e-12

        ait_phi, ait_theta = intermediate_to_native(AIT(), 2.0 * sqrt(2.0) * R2D * (1.0 + 1e-14), 0.0)
        @test angle_approx(ait_phi, π; atol=1e-10)
        @test ait_theta ≈ 0.0 atol=1e-12
    end

    @testset "Non-default CYP parameters" begin
        # CYP PV parameters should parse from the latitude axis and round-trip.
        hdr = Dict(
            "NAXIS" => 2, "CTYPE1" => "RA---CYP", "CTYPE2" => "DEC--CYP",
            "PV2_1" => 0.75, "PV2_2" => 1.5,
        )
        @test WCS(hdr).projection == CYP(0.75, 1.5)

        proj = CYP(0.75, 1.5)
        phi_r = 30.0 * D2R
        theta_r = 20.0 * D2R
        x, y = native_to_intermediate(proj, phi_r, theta_r)
        phi2, theta2 = intermediate_to_native(proj, x, y)
        @test angle_approx(phi2, phi_r; atol=1e-12)
        @test theta2 ≈ theta_r atol=1e-12
    end

    @testset "Float32 edge precision" begin
        # Projection edge paths should not silently promote Float32 inputs.
        for proj in (ZEA(), CEA(0.75), PAR(), MOL(), PCO(), AIT())
            phi, theta = intermediate_to_native(proj, Float32(0), Float32(0))
            @test phi isa Float32
            @test theta isa Float32
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "AIT projection" begin
    proj = AIT()
    @testset "Forward/inverse round-trip" begin
        for (phi_d, theta_d) in [(0.0, 0.0), (90.0, 30.0), (-90.0, -30.0), (45.0, 45.0)]
            phi_r   = phi_d   * D2R
            theta_r = theta_d * D2R
            x, y    = native_to_intermediate(proj, phi_r, theta_r)
            phi2, theta2 = intermediate_to_native(proj, x, y)
            @test phi2   ≈ phi_r    atol=1e-12
            @test theta2 ≈ theta_r  atol=1e-12
        end
    end

    @testset "Wrapped longitude keeps x-coordinate branch" begin
        # AIT inverse calls can receive longitudes outside [-180, 180] after rotation.
        x1, y1 = native_to_intermediate(proj, -20.0 * D2R, 5.0 * D2R)
        x2, y2 = native_to_intermediate(proj, 340.0 * D2R, 5.0 * D2R)
        @test x2 ≈ x1 atol=1e-12
        @test y2 ≈ y1 atol=1e-12
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Full TAN celestial WCS" begin
    # Reference: Orion Nebula region, typical HST-like header
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---TAN",
        "CTYPE2" => "DEC--TAN",
        "CRPIX1" => 512.0,   "CRPIX2" => 512.0,
        "CRVAL1" => 83.8221, "CRVAL2" => -5.3911,
        "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4,
    )
    wcs = WCS(hdr)

    @testset "Reference pixel → CRVAL" begin
        world = pixel_to_world(wcs, [512.0, 512.0])
        @test world[1] ≈ wcs.crval[1]  atol=1e-6
        @test world[2] ≈ wcs.crval[2]  atol=1e-6
    end

    @testset "Round-trip center" begin
        pix_in = [512.0, 512.0]
        world  = pixel_to_world(wcs, pix_in)
        pix_out = world_to_pixel(wcs, world)
        @test pix_out ≈ pix_in  atol=1e-8
    end

    @testset "Round-trip at corners" begin
        corners = [[1.0, 1.0], [1024.0, 1.0], [1.0, 1024.0], [1024.0, 1024.0]]
        for pix_in in corners
            world   = pixel_to_world(wcs, pix_in)
            pix_out = world_to_pixel(wcs, world)
            @test pix_out ≈ pix_in  atol=1e-6
        end
    end

    @testset "Small offset pixel" begin
        # Move 1 pixel in RA (negative CDELT1 → RA increases as pixel decreases)
        pix_center = [512.0, 512.0]
        pix_offset = [511.0, 512.0]
        w_center = pixel_to_world(wcs, pix_center)
        w_offset = pixel_to_world(wcs, pix_offset)
        # RA should have increased by approximately |CDELT1| degrees
        dra = w_offset[1] - w_center[1]
        @test dra ≈ -wcs.cd[1,1]  atol=1e-5  # cd[1,1] = CDELT1 ≈ -2.7778e-4
    end

    @testset "Batch matches scalar" begin
        pix_mat = hcat([1.0, 1.0], [512.0, 512.0], [1024.0, 1024.0])
        w_batch = pixel_to_world(wcs, pix_mat)
        for k in 1:3
            @test w_batch[:, k] ≈ pixel_to_world(wcs, pix_mat[:, k])  atol=1e-14
        end
        p_batch = world_to_pixel(wcs, w_batch)
        for k in 1:3
            @test p_batch[:, k] ≈ world_to_pixel(wcs, w_batch[:, k])  atol=1e-10
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "AIT full celestial WCS" begin
    # All-sky AIT with GALACTIC coordinates.
    # For standard all-sky maps, CRVAL gives the celestial coordinates at
    # native phi=180°, theta=0° (the anti-meridian of the native frame).
    # With LONPOLE=180° (default for delta0=0), the native phi=180°, theta=0°
    # maps to (alpha_p, 0°) = (CRVAL_lon, 0°).
    # A typical full-sky map uses CRVAL1=0 (or 180) with LONPOLE=180.
    # Here we test the round-trip: the projection center is at native (0, 0),
    # which maps to some celestial coords we can verify by round-trip.
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "GLON-AIT",
        "CTYPE2" => "GLAT-AIT",
        "CRPIX1" => 360.5, "CRPIX2" => 180.5,
        "CRVAL1" => 0.0,   "CRVAL2" => 0.0,
        "CDELT1" => -0.5,  "CDELT2" => 0.5,
    )
    wcs = WCS(hdr)
    @test wcs.projection isa AIT

    @testset "Round-trip center" begin
        pix_in  = [360.5, 180.5]
        world   = pixel_to_world(wcs, pix_in)
        pix_out = world_to_pixel(wcs, world)
        @test pix_out ≈ pix_in  atol=1e-8
    end

    @testset "Round-trip off-center" begin
        # Use pixels that map to valid AIT domain.
        # AIT domain: (x_deg/162)² + (y_deg/81)² < 1 (approximately)
        # where x = CDELT1*(p1-CRPIX1), y = CDELT2*(p2-CRPIX2).
        # |x| < 162° → |p1-360.5| < 324.  |y| < 81° → |p2-180.5| < 162.
        for pix_in in ([300.0, 150.0], [400.0, 190.0], [360.5, 120.0])
            world   = pixel_to_world(wcs, pix_in)
            pix_out = world_to_pixel(wcs, world)
            @test pix_out ≈ pix_in  atol=1e-6
        end
    end

    @testset "Latitude at equator pixel row" begin
        # Any pixel in the equatorial row (CRPIX2) should have latitude = 0°
        world = pixel_to_world(wcs, [100.0, 180.5])
        @test world[2] ≈ 0.0  atol=1e-6
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "CAR full celestial WCS" begin
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---CAR",
        "CTYPE2" => "DEC--CAR",
        "CRPIX1" => 181.0, "CRPIX2" => 91.0,
        "CRVAL1" => 0.0,   "CRVAL2" => 0.0,
        "CDELT1" => -1.0,  "CDELT2" => 1.0,
    )
    wcs = WCS(hdr)
    @test wcs.projection isa CAR

    @testset "Equatorial pixel row has latitude = 0" begin
        # Any pixel on the central row (y=0 intermediate) → DEC = 0°
        world = pixel_to_world(wcs, [100.0, 91.0])
        @test world[2] ≈ 0.0  atol=1e-10
    end

    @testset "Round-trip" begin
        pix_in  = [100.0, 50.0]
        world   = pixel_to_world(wcs, pix_in)
        pix_out = world_to_pixel(wcs, world)
        @test pix_out ≈ pix_in  atol=1e-8
    end

    @testset "Round-trip at multiple points" begin
        # Avoid poles (Dec = ±90°) where longitude is undefined and
        # the round-trip is degenerate.
        for pix_in in ([50.0, 20.0], [100.0, 50.0], [181.0, 91.0], [50.0, 120.0])
            world   = pixel_to_world(wcs, pix_in)
            pix_out = world_to_pixel(wcs, world)
            @test pix_out ≈ pix_in  atol=1e-8
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "CEA full celestial WCS" begin
    # CEA headers from real all-sky products commonly include PV latitude parameter 1.
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---CEA",
        "CTYPE2" => "DEC--CEA",
        "CRPIX1" => 181.0, "CRPIX2" => 91.0,
        "CRVAL1" => 0.0,   "CRVAL2" => 0.0,
        "CDELT1" => -1.0,  "CDELT2" => 1.0,
        "PV2_1"  => 1.0,
    )
    wcs = WCS(hdr)
    @test wcs.projection == CEA(1.0)

    @testset "Equatorial pixel row has latitude = 0" begin
        # The CEA equator lies on the row whose intermediate latitude ordinate is zero.
        world = pixel_to_world(wcs, [100.0, 91.0])
        @test world[2] ≈ 0.0 atol=1e-10
    end

    @testset "Round-trip away from poles" begin
        # Avoid poles where longitude is undefined.
        for pix_in in ([50.0, 50.0], [100.0, 70.0], [181.0, 91.0], [50.0, 120.0])
            world = pixel_to_world(wcs, pix_in)
            pix_out = world_to_pixel(wcs, world)
            @test pix_out ≈ pix_in atol=1e-8
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Unknown projection raises error on transform" begin
    proj = UnknownProjection("XXX")
    @test_throws ErrorException intermediate_to_native(proj, 0.0, 0.0)
    @test_throws ErrorException native_to_intermediate(proj, 0.0, 1.0)
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "FITSIO extension" begin
    # FITSIO headers and HDUs should parse without making FITSIO a hard dependency.
    header = FITSIO.FITSHeader(
        ["NAXIS", "CTYPE1", "CRPIX1", "CRVAL1", "CDELT1"],
        [1, "FREQ", 1.0, 100.0, 2.0],
        fill("", 5),
    )

    wcs = WCS(header)
    @test pixel_to_world(wcs, [3.0]) ≈ [104.0]

    mktempdir() do dir
        path = joinpath(dir, "linear.fits")
        FITSIO.FITS(path, "w") do file
            FITSIO.write(file, zeros(Float32, 4); header=header)
            hdu_wcs = WCS(file[1])
            @test pixel_to_world(hdu_wcs, [3.0]) ≈ [104.0]

            fobj_wcs = WCS(header; fobj=file, minerr=0.1)
            @test fobj_wcs.aux isa NoAuxiliaryWCSData

            hdu_fobj_wcs = WCS(file[1]; fobj=file, minerr=0.1)
            @test hdu_fobj_wcs.aux isa NoAuxiliaryWCSData

            lookup_header = FITSIO.FITSHeader(
                [
                    "NAXIS", "CTYPE1", "CTYPE2",
                    "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2",
                    "CPDIS1", "DP1.EXTVER", "DP1.AXIS.1",
                    "D2IMDIS2", "D2IM2.EXTVER", "D2IM2.AXIS.2",
                ],
                [
                    2, "X", "Y",
                    1.0, 1.0, 0.0, 0.0, 1.0, 1.0,
                    "LOOKUP", 1, 1,
                    "LOOKUP", 2, 2,
                ],
                fill("", 15),
            )
            table_header = FITSIO.FITSHeader(
                ["CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2"],
                [0.0, 0.0, 0.0, 0.0, 1.0, 1.0],
                fill("", 6),
            )
            FITSIO.write(file, fill(Float32(0.5), 2, 2); header=table_header, name="WCSDVARR", ver=1)
            FITSIO.write(file, fill(Float32(0.25), 2, 2); header=table_header, name="D2IMARR", ver=2)

            lookup_wcs = WCS(lookup_header; fobj=file)
            @test lookup_wcs.aux isa AuxiliaryWCSData
            @test lookup_wcs.aux.cpdis[1] isa LookupTable2D
            @test lookup_wcs.aux.det2im[2] isa LookupTable2D
            @test pixel_to_world(lookup_wcs, [1.0, 1.0]) ≈ [0.5, 0.25]
            @test world_to_pixel(lookup_wcs, [0.5, 0.25]) ≈ [1.0, 1.0]
            @test world_to_pixel(lookup_wcs, reshape([0.5, 0.25], 2, 1)) ≈ reshape([1.0, 1.0], 2, 1)

            tab_header = FITSIO.FITSHeader(
                [
                    "NAXIS", "CTYPE1", "CRPIX1", "CRVAL1", "CDELT1",
                    "PS1_0", "PS1_1", "PS1_2",
                ],
                [
                    1, "FREQ-TAB", 1.0, 1.0, 1.0,
                    "WCS-TABLE", "FREQS", "PIXELS",
                ],
                fill("", 8),
            )
            FITSIO.write(
                file,
                ["FREQS", "PIXELS"],
                Any[
                    reshape([10.0, 20.0, 40.0], 3, 1),
                    reshape([1.0, 2.0, 3.0], 3, 1),
                ];
                name = "WCS-TABLE",
            )

            tab_wcs = WCS(tab_header; fobj=file)
            @test tab_wcs.aux.tabular isa TabularWCSData
            @test pixel_to_world(tab_wcs, [2.5]) ≈ [30.0]
            @test world_to_pixel(tab_wcs, [30.0]) ≈ [2.5]
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "FITSFiles extension" begin
    # FITSFiles card vectors and HDUs should parse through the same WCS path.
    cards = FITSFiles.Card[
        FITSFiles.Card("SIMPLE", true),
        FITSFiles.Card("BITPIX", -32),
        FITSFiles.Card("NAXIS", 1),
        FITSFiles.Card("NAXIS1", 4),
        FITSFiles.Card("CTYPE1", "FREQ"),
        FITSFiles.Card("CRPIX1", 1.0),
        FITSFiles.Card("CRVAL1", 100.0),
        FITSFiles.Card("CDELT1", 2.0),
    ]

    wcs = WCS(cards)
    @test pixel_to_world(wcs, [3.0]) ≈ [104.0]

    hdu = FITSFiles.HDU(cards)
    hdu_wcs = WCS(hdu)
    @test pixel_to_world(hdu_wcs, [3.0]) ≈ [104.0]

    hdus = [hdu]

    cards_fobj_wcs = WCS(cards; fobj=hdus, minerr=0.1)
    @test cards_fobj_wcs.aux isa NoAuxiliaryWCSData

    hdu_fobj_wcs = WCS(hdu; fobj=hdus, minerr=0.1)
    @test hdu_fobj_wcs.aux isa NoAuxiliaryWCSData

    lookup_cards = FITSFiles.Card[
        FITSFiles.Card("SIMPLE", true),
        FITSFiles.Card("BITPIX", -32),
        FITSFiles.Card("NAXIS", 2),
        FITSFiles.Card("NAXIS1", 4),
        FITSFiles.Card("NAXIS2", 4),
        FITSFiles.Card("CTYPE1", "X"),
        FITSFiles.Card("CTYPE2", "Y"),
        FITSFiles.Card("CRPIX1", 1.0),
        FITSFiles.Card("CRPIX2", 1.0),
        FITSFiles.Card("CRVAL1", 0.0),
        FITSFiles.Card("CRVAL2", 0.0),
        FITSFiles.Card("CDELT1", 1.0),
        FITSFiles.Card("CDELT2", 1.0),
        FITSFiles.Card("CPDIS1", "LOOKUP"),
        FITSFiles.Card("HIERARCH", "DP1.EXTVER", 1, ""),
        FITSFiles.Card("HIERARCH", "DP1.AXIS.1", 1, ""),
        FITSFiles.Card("D2IMDIS2", "LOOKUP"),
        FITSFiles.Card("HIERARCH", "D2IM2.EXTVER", 2, ""),
        FITSFiles.Card("HIERARCH", "D2IM2.AXIS.2", 2, ""),
    ]
    lookup_hdu = FITSFiles.HDU(lookup_cards)
    table_cards_1 = FITSFiles.Card[
        FITSFiles.Card("EXTNAME", "WCSDVARR"),
        FITSFiles.Card("EXTVER", 1),
        FITSFiles.Card("CRPIX1", 0.0),
        FITSFiles.Card("CRPIX2", 0.0),
        FITSFiles.Card("CRVAL1", 0.0),
        FITSFiles.Card("CRVAL2", 0.0),
        FITSFiles.Card("CDELT1", 1.0),
        FITSFiles.Card("CDELT2", 1.0),
    ]
    table_cards_2 = FITSFiles.Card[
        FITSFiles.Card("EXTNAME", "D2IMARR"),
        FITSFiles.Card("EXTVER", 2),
        FITSFiles.Card("CRPIX1", 0.0),
        FITSFiles.Card("CRPIX2", 0.0),
        FITSFiles.Card("CRVAL1", 0.0),
        FITSFiles.Card("CRVAL2", 0.0),
        FITSFiles.Card("CDELT1", 1.0),
        FITSFiles.Card("CDELT2", 1.0),
    ]
    lookup_hdus = [
        lookup_hdu,
        FITSFiles.HDU(FITSFiles.Image, fill(Float32(0.5), 2, 2), table_cards_1),
        FITSFiles.HDU(FITSFiles.Image, fill(Float32(0.25), 2, 2), table_cards_2),
    ]

    lookup_wcs = WCS(lookup_hdu; fobj=lookup_hdus)
    @test lookup_wcs.aux isa AuxiliaryWCSData
    @test lookup_wcs.aux.cpdis[1] isa LookupTable2D
    @test lookup_wcs.aux.det2im[2] isa LookupTable2D
    @test pixel_to_world(lookup_wcs, [1.0, 1.0]) ≈ [0.5, 0.25]
    @test world_to_pixel(lookup_wcs, [0.5, 0.25]) ≈ [1.0, 1.0]
    @test world_to_pixel(lookup_wcs, reshape([0.5, 0.25], 2, 1)) ≈ reshape([1.0, 1.0], 2, 1)

    tab_cards = FITSFiles.Card[
        FITSFiles.Card("SIMPLE", true),
        FITSFiles.Card("BITPIX", -32),
        FITSFiles.Card("NAXIS", 1),
        FITSFiles.Card("NAXIS1", 3),
        FITSFiles.Card("CTYPE1", "FREQ-TAB"),
        FITSFiles.Card("CRPIX1", 1.0),
        FITSFiles.Card("CRVAL1", 1.0),
        FITSFiles.Card("CDELT1", 1.0),
        FITSFiles.Card("PS1_0", "WCS-TABLE"),
        FITSFiles.Card("PS1_1", "FREQS"),
        FITSFiles.Card("PS1_2", "PIXELS"),
    ]
    tab_hdu = FITSFiles.HDU(tab_cards)
    table_hdu = FITSFiles.HDU(
        (FREQS = [10.0, 20.0, 40.0], PIXELS = [1.0, 2.0, 3.0]),
        [
            FITSFiles.Card("EXTNAME", "WCS-TABLE"),
            FITSFiles.Card("EXTVER", 1),
        ],
    )
    tab_wcs = WCS(tab_hdu; fobj=[tab_hdu, table_hdu])
    @test tab_wcs.aux.tabular isa TabularWCSData
    @test pixel_to_world(tab_wcs, [2.5]) ≈ [30.0]
    @test world_to_pixel(tab_wcs, [30.0]) ≈ [2.5]
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Paper IV forward distortion pipeline" begin
    # Loaded D2IM and CPDIS tables should affect scalar and batched pixel-to-world transforms.
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "X", "CTYPE2" => "Y",
        "CRPIX1" => 1.0, "CRPIX2" => 1.0,
        "CDELT1" => 1.0, "CDELT2" => 1.0,
        "CPDIS1" => "LOOKUP", "DP1.AXIS.1" => 1,
        "D2IMDIS2" => "LOOKUP", "D2IM2.AXIS.2" => 2,
    )

    wcs = WCS(hdr; fobj=FakeLookupFobj())

    @test wcs.pipeline isa DistortionPipeline
    @test pixel_to_focal(wcs.pipeline, [1.0, 1.0], Val(2)) ≈ [2.0, 1.25]
    @test focal_to_pixel(wcs.pipeline, [2.0, 1.25], Val(2)) ≈ [1.0, 1.0]
    @test pixel_to_world(wcs, [1.0, 1.0]) ≈ [1.0, 0.25]
    @test world_to_pixel(wcs, [1.0, 0.25]) ≈ [1.0, 1.0]

    pixels = [1.0 2.0; 1.0 1.0]
    worlds = pixel_to_world(wcs, pixels)
    @test worlds[:, 1] ≈ pixel_to_world(wcs, pixels[:, 1])
    @test worlds[:, 2] ≈ pixel_to_world(wcs, pixels[:, 2])
    @test world_to_pixel(wcs, worlds) ≈ pixels
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Paper IV and SIP composition" begin
    # Nonzero SIP and Paper IV lookup stages should compose in detector, SIP, CPDIS order.
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "X", "CTYPE2" => "Y",
        "CRPIX1" => 1.0, "CRPIX2" => 1.0,
        "CDELT1" => 1.0, "CDELT2" => 1.0,
        "A_ORDER" => 2, "B_ORDER" => 2,
        "A_2_0" => 0.1, "B_0_2" => 0.2,
        "CPDIS1" => "LOOKUP", "DP1.AXIS.1" => 1,
        "D2IMDIS2" => "LOOKUP", "D2IM2.AXIS.2" => 2,
    )

    wcs = WCS(hdr; fobj=FakeLookupFobj())
    pixel = [1.0, 1.0]
    focal_ref = [2.0, 1.3]
    world_ref = [1.0, 0.3]

    @test wcs.pipeline.sip isa SIPDistortion
    @test pixel_to_focal(wcs.pipeline, pixel, Val(2)) ≈ focal_ref
    @test focal_to_pixel(wcs.pipeline, focal_ref, Val(2)) ≈ pixel atol=1e-10
    @test pixel_to_world(wcs, pixel) ≈ world_ref
    @test world_to_pixel(wcs, world_ref) ≈ pixel atol=1e-10

    pixels = [1.0 1.5; 1.0 1.0]
    worlds = pixel_to_world(wcs, pixels)
    @test world_to_pixel(wcs, worlds) ≈ pixels atol=1e-10
    @test pixel_to_world(wcs, pixel) ≈ world_ref atol=1e-6
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Paper IV CPDIS-only interpolation" begin
    # Sparse CPDIS tables should apply exact and bilinear offsets using table metadata.
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "X", "CTYPE2" => "Y",
        "CRPIX1" => 1.0, "CRPIX2" => 1.0,
        "CRVAL1" => 0.0, "CRVAL2" => 0.0,
        "CDELT1" => 1.0, "CDELT2" => 1.0,
        "CPDIS1" => "LOOKUP", "DP1.AXIS.1" => 1,
        "CPDIS2" => "LOOKUP", "DP2.AXIS.2" => 2,
    )

    wcs = WCS(hdr; fobj=FakeCPDISInterpolationFobj())
    @test wcs.aux.cpdis[1] isa LookupTable2D
    @test wcs.aux.cpdis[2] isa LookupTable2D

    @test pixel_to_world(wcs, [42.0, 22.0]) ≈ [41.5, 21.0]
    @test pixel_to_world(wcs, [13.0, 23.0]) ≈ [12.0, 22.7]
    @test focal_to_pixel(wcs.pipeline, [42.5, 22.0], Val(2)) ≈ [42.0, 22.0]
    @test world_to_pixel(wcs, [41.5, 21.0]) ≈ [42.0, 22.0]
    @test world_to_pixel(wcs, [12.0, 22.7]) ≈ [13.0, 23.0]

    for (pixel, world_ref) in [
        ([43.0, 22.0], [42.25, 21.0]),
        ([41.0, 22.0], [40.25, 21.0]),
        ([42.0, 23.0], [41.25, 22.0]),
        ([42.0, 21.0], [41.25, 20.0]),
    ]
        @test pixel_to_world(wcs, pixel) ≈ world_ref
        @test world_to_pixel(wcs, world_ref) ≈ pixel
    end

    for (pixel, world_ref) in [
        ([14.5, 23.0], [13.5, 22.35]),
        ([11.5, 23.0], [10.5, 22.35]),
        ([13.0, 24.5], [12.0, 23.85]),
        ([13.0, 21.5], [12.0, 20.85]),
    ]
        @test pixel_to_world(wcs, pixel) ≈ world_ref
        @test world_to_pixel(wcs, world_ref) ≈ pixel
    end

    for (pixel, world_ref) in [
        ([46.0, 22.0], [45.0, 21.0]),
        ([38.0, 22.0], [37.0, 21.0]),
        ([42.0, 26.0], [41.0, 25.0]),
        ([42.0, 18.0], [41.0, 17.0]),
        ([19.0, 23.0], [18.0, 22.0]),
        ([7.0, 23.0], [6.0, 22.0]),
        ([13.0, 29.0], [12.0, 28.0]),
        ([13.0, 17.0], [12.0, 16.0]),
    ]
        @test pixel_to_world(wcs, pixel) ≈ world_ref
        @test world_to_pixel(wcs, world_ref) ≈ pixel
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Paper IV forward distortion Astropy regression" begin
    # Astropy's dist_lookup fixture exercises real SIP, D2IM, and CPDIS tables.
    refs = [
        ([1.0, 1.0], [5.5264578963294078, -72.051718954260323]),
        ([24.0, 24.0], [5.5278318183240547, -72.051898631854385]),
        ([512.0, 512.0], [5.5571048564879018, -72.055692801825842]),
        ([1024.0, 1024.0], [5.5880253054822395, -72.059638092645642]),
        ([100.0, 200.0], [5.5343152229426584, -72.053752245702228]),
        ([700.25, 300.75], [5.5607016650989554, -72.05211216803022]),
        ([2048.0, 1024.0], [5.6305686380276061, -72.054571792078278]),
    ]
    focal_refs = [
        ([1.0, 1.0], [34.071093418170157, 0.62674039300433648]),
        ([24.0, 24.0], [56.301206127480675, 23.608765313128124]),
        ([512.0, 512.0], [530.78963861021305, 511.09249280708974]),
        ([1024.0, 1024.0], [1033.4930976678936, 1022.0527435313898]),
        ([100.0, 200.0], [130.29194796774041, 199.32780993366652]),
        ([700.25, 300.75], [713.51391137340056, 300.39852642546708]),
        ([2048.0, 1024.0], [2048.0122948031303, 1024.0008731349378]),
    ]

    FITSIO.FITS(joinpath(@__DIR__, "dist_lookup.fits.gz")) do file
        wcs = WCS(file[2]; fobj=file)

        @test wcs.aux isa AuxiliaryWCSData
        @test wcs.aux.det2im[1] isa LookupTable2D
        @test wcs.aux.cpdis[1] isa LookupTable2D
        @test wcs.aux.cpdis[2] isa LookupTable2D

        for (pixel, focal_ref) in focal_refs
            @test pixel_to_focal(wcs.pipeline, pixel, Val(2)) ≈ focal_ref atol=1e-10
            @test focal_to_pixel(wcs.pipeline, focal_ref, Val(2)) ≈ pixel atol=1e-10
        end

        for (pixel, world_ref) in refs
            @test pixel_to_world(wcs, pixel) ≈ world_ref atol=1e-10
        end

        for (pixel, world_ref) in refs
            @test world_to_pixel(wcs, world_ref) ≈ pixel atol=1e-8
        end

        pixels = hcat((ref[1] for ref in refs)...)
        worlds = pixel_to_world(wcs, pixels)
        for (k, (_, world_ref)) in pairs(refs)
            @test worlds[:, k] ≈ world_ref atol=1e-10
        end

        # Test batch version
        inverse_pixels = hcat((ref[1] for ref in refs)...)
        inverse_worlds = hcat((ref[2] for ref in refs)...)
        @test maximum(abs.(world_to_pixel(wcs, inverse_worlds) .- inverse_pixels)) <= 1e-8
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Type stability" begin
    # Plan Milestone 8: type stability where practical.
    # pixel_to_world and world_to_pixel must return Vector{Float64} regardless of
    # which code path is taken (linear, celestial, SIP).

    hdr_lin = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "X", "CTYPE2" => "Y",
        "CRPIX1" => 1.0, "CRPIX2" => 1.0,
        "CRVAL1" => 0.0, "CRVAL2" => 0.0,
        "CDELT1" => 1.0, "CDELT2" => 1.0,
    )
    hdr_tan = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN",
        "CRPIX1" => 512.0,      "CRPIX2" => 512.0,
        "CRVAL1" => 83.8221,    "CRVAL2" => -5.3911,
        "CDELT1" => -2.7778e-4, "CDELT2" =>  2.7778e-4,
    )
    hdr_sip = Dict(
        "NAXIS"   => 2,
        "CTYPE1"  => "RA---TAN-SIP", "CTYPE2"  => "DEC--TAN-SIP",
        "CRPIX1"  => 512.0,           "CRPIX2"  => 512.0,
        "CRVAL1"  => 150.0,           "CRVAL2"  => 2.5,
        "CDELT1"  => -2.7778e-4,      "CDELT2"  =>  2.7778e-4,
        "A_ORDER" => 2,  "A_2_0" => 5.0e-6,  "A_0_2" => 2.0e-6,
        "B_ORDER" => 2,  "B_2_0" => 1.0e-6,  "B_1_1" => 3.0e-6,
    )
    hdr_3d = Dict(
        "NAXIS"  => 3,
        "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CTYPE3" => "FREQ",
        "CRPIX1" => 50.0,       "CRPIX2" => 50.0,        "CRPIX3" => 1.0,
        "CRVAL1" => 10.0,       "CRVAL2" => 25.0,         "CRVAL3" => 1.42e9,
        "CDELT1" => -0.01,      "CDELT2" =>  0.01,        "CDELT3" => 1.0e6,
    )

    wcs_lin = WCS(hdr_lin)
    wcs_tan = WCS(hdr_tan)
    wcs_sip = WCS(hdr_sip)
    wcs_3d  = WCS(hdr_3d)

    # Parsed transforms should store small numeric WCS state in static arrays.
    @test WCS(hdr_lin) isa WCSTransform
    @test WCS(hdr_tan) isa WCSTransform
    @test WCS(hdr_sip) isa WCSTransform
    @test wcs_tan.crpix isa SVector{2,Float64}
    @test wcs_tan.crval isa SVector{2,Float64}
    @test wcs_tan.cd isa SMatrix{2,2,Float64,4}

    # pixel_to_world always returns SVector regardless of input container type.
    @test @inferred(pixel_to_world(wcs_lin, [2.0, 3.0]))  isa SVector{2,Float64}
    @test @inferred(pixel_to_world(wcs_tan, [400.0, 300.0])) isa SVector{2,Float64}
    @test @inferred(pixel_to_world(wcs_sip, [600.0, 500.0])) isa SVector{2,Float64}
    @test @inferred(pixel_to_world(wcs_3d,  [40.0, 60.0, 5.0])) isa SVector{3,Float64}
    @test @inferred(pixel_to_world(wcs_lin, SVector(2.0, 3.0))) isa SVector{2,Float64}
    @test @inferred(pixel_to_world(wcs_lin, SVector(Float32(2), Float32(3)))) isa SVector{2,Float32}

    # world_to_pixel always returns SVector regardless of input container type.
    w_lin = pixel_to_world(wcs_lin, [2.0, 3.0])
    w_tan = pixel_to_world(wcs_tan, [400.0, 300.0])
    w_sip = pixel_to_world(wcs_sip, [600.0, 500.0])
    w_3d  = pixel_to_world(wcs_3d,  [40.0, 60.0, 5.0])
    @test @inferred(world_to_pixel(wcs_lin, w_lin)) isa SVector{2,Float64}
    @test @inferred(world_to_pixel(wcs_tan, w_tan)) isa SVector{2,Float64}
    @test @inferred(world_to_pixel(wcs_sip, w_sip)) isa SVector{2,Float64}
    @test @inferred(world_to_pixel(wcs_3d,  w_3d))  isa SVector{3,Float64}

    # Vector input also returns SVector (input type does not dictate output type).
    @test @inferred(world_to_pixel(wcs_lin, collect(w_lin))) isa SVector{2,Float64}
    @test @inferred(world_to_pixel(wcs_lin, SVector(2.0, 3.0))) isa SVector{2,Float64}

    # Batch (matrix) forms return Matrix{Float64}.
    pix_mat = [1.0 512.0; 1.0 512.0]
    @test @inferred(pixel_to_world(wcs_lin, pix_mat)) isa Matrix{Float64}
    @test @inferred(pixel_to_world(wcs_tan, pix_mat)) isa Matrix{Float64}
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Randomized property tests" begin
    # Plan Testing Strategy section 5: use well-conditioned randomized WCS objects
    # with a fixed RNG seed so results are deterministic.
    rng = MersenneTwister(0x5f3759df)

    @testset "Linear WCS round-trip (random)" begin
        # Generate well-conditioned random linear WCS objects and verify that
        # pixel → world → pixel recovers the original pixel within tolerance.
        for _ in 1:20
            naxis = rand(rng, 2:4)
            crpix = 50.0 .+ 200.0 .* rand(rng, naxis)
            crval = -180.0 .+ 360.0 .* rand(rng, naxis)
            cdelt = 1e-3 .+ 1e-1 .* rand(rng, naxis)

            hdr = Dict{String,Any}("NAXIS" => naxis)
            for i in 1:naxis
                hdr["CTYPE$i"] = string('A' - 1 + i)
                hdr["CRPIX$i"] = crpix[i]
                hdr["CRVAL$i"] = crval[i]
                hdr["CDELT$i"] = cdelt[i]
            end
            wcs = WCS(hdr)

            # Random pixel point well away from any edge effects.
            pix = crpix .+ 10.0 .* (rand(rng, naxis) .- 0.5)
            world = pixel_to_world(wcs, pix)
            pix2 = world_to_pixel(wcs, world)
            @test pix2 ≈ pix  atol=1e-8
        end
    end

    @testset "PC matrix inverse consistency (random)" begin
        # Random well-conditioned rotation matrices via QR decomposition of a random matrix.
        # Verify CD = CDELT * PC elementwise, and that the transform round-trips.
        for _ in 1:10
            Q, _ = qr(randn(rng, 2, 2))
            pc = Matrix(Q)   # orthogonal rotation matrix
            cdelt = [-1e-4, 1e-4] .* (1.0 .+ 0.5 .* rand(rng, 2))

            hdr = Dict{String,Any}(
                "NAXIS"  => 2,
                "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN",
                "CRPIX1" => 256.0,      "CRPIX2" => 256.0,
                "CRVAL1" => 45.0,       "CRVAL2" => 20.0,
                "CDELT1" => cdelt[1],   "CDELT2" => cdelt[2],
                "PC1_1"  => pc[1,1],    "PC1_2"  => pc[1,2],
                "PC2_1"  => pc[2,1],    "PC2_2"  => pc[2,2],
            )
            wcs = WCS(hdr)

            # CD matrix should equal CDELT * PC elementwise (up to float round-off).
            for i in 1:2, j in 1:2
                @test wcs.cd[i,j] ≈ cdelt[i] * pc[i,j]  atol=1e-12
            end

            # Round-trip near reference pixel.
            pix = [256.0 + randn(rng), 256.0 + randn(rng)]
            @test world_to_pixel(wcs, pixel_to_world(wcs, pix)) ≈ pix  atol=1e-4
        end
    end

    @testset "Projection inverse consistency (random, away from singularities)" begin
        # Each projection should invert exactly for points well within its domain.
        # Singularity-sensitive projections use restricted ranges.
        for (proj, phi_range, theta_range) in [
            (AZP(),    (-π, π), (0.1, π/2)),   # AZP default: central perspective
            (SZP(),    (-π, π), (0.1, π/2)),   # SZP default: central perspective
            (TAN(),    (-π, π), (0.1, π/2)),   # TAN: theta > 0
            (SIN(),    (-π, π), (0.1, π/2)),   # SIN: R_theta <= 1
            (STG(),    (-π, π), (0.1, π/2)),   # STG: theta > -90
            (ARC(),    (-π, π), (0.0, π/2)),   # ARC: all theta valid
            (ZEA(),    (-π, π), (0.0, π/2)),   # ZEA: all theta valid
            (CAR(),    (-π, π), (-π/2, π/2)),  # CAR: full sky
            (CEA(0.8), (-π, π), (-π/3, π/3)),  # CEA: restricted latitude
            (CYP(),    (-π, π), (-π/3, π/3)),  # CYP: avoid perspective singularities
            (MER(),    (-π, π), (-π/3, π/3)),  # MER: avoid polar singularities
            (SFL(),    (-π, π), (-π/3, π/3)),  # SFL: avoid undefined polar longitude
            (PAR(),    (-π, π), (-π/3, π/3)),  # PAR: avoid domain edges
            (MOL(),    (-π, π), (-π/3, π/3)),  # MOL: avoid polar auxiliary degeneracy
            (PCO(),    (-π, π), (-π/3, π/3)),      # PCO: avoid polar branch degeneracy
            (AIT(),    (-π, π), (-π/2, π/2)),  # AIT: full sky
        ]
            for _ in 1:10
                # Sample a point in the valid domain.
                phi_lo, phi_hi = phi_range
                th_lo, th_hi   = theta_range
                phi   = phi_lo + (phi_hi - phi_lo) * rand(rng)
                theta = th_lo  + (th_hi  - th_lo)  * rand(rng)

                # SIN projection: ensure the point is inside the unit circle.
                if proj isa SIN
                    while true
                        Rth = abs(cos(theta))   # = |cos(theta)|
                        Rth < 0.99 && break
                        theta = th_lo + (th_hi - th_lo) * rand(rng)
                    end
                end

                x, y = FITSWCS.native_to_intermediate(proj, phi, theta)
                phi2, theta2 = FITSWCS.intermediate_to_native(proj, x, y)

                # Angles are equal modulo 2π.
                @test angle_approx(phi2, phi;   atol=1e-10)
                @test theta2 ≈ theta  atol=1e-10
            end
        end
    end

    @testset "TAN celestial WCS round-trip (random)" begin
        # Random TAN celestial WCS with well-conditioned CDELT and near-pole-free
        # CRVAL should round-trip to within floating-point projection tolerance.
        for _ in 1:15
            crval1 = -170.0 + 340.0 * rand(rng)      # RA away from 0/360 boundary
            crval2 = -60.0  + 120.0 * rand(rng)      # Dec away from ±90
            cdelt  = 1e-4   + 9e-4  * rand(rng)      # 0.1–1 arcsec/pixel range

            hdr = Dict{String,Any}(
                "NAXIS"  => 2,
                "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN",
                "CRPIX1" => 256.0,      "CRPIX2" => 256.0,
                "CRVAL1" => crval1,     "CRVAL2" => crval2,
                "CDELT1" => -cdelt,     "CDELT2" => cdelt,
            )
            wcs = WCS(hdr)

            # Pixel near the reference point.
            pix = [256.0 + 10.0*randn(rng), 256.0 + 10.0*randn(rng)]
            world = pixel_to_world(wcs, pix)
            @test world_to_pixel(wcs, world) ≈ pix  atol=1e-5
        end
    end
end

    @testset "compute_native_pole type stability" begin
        halfpi64 = π / 2
        halfpi32 = Float32(π / 2)

        # ── Float64 ──────────────────────────────────────────
        # Branch 1: zenithal (theta0 ≈ 90°)
        @test compute_native_pole(0.0, 0.0, 0.0, halfpi64, 0.0, halfpi64) isa Tuple{Float64,Float64}
        @inferred compute_native_pole(0.0, 0.0, 0.0, halfpi64, 0.0)

        # Branch 2: degenerate R (sin(theta0) ≈ 0, cos(dphi) ≈ 0)
        @test compute_native_pole(0.0, 0.0, halfpi64, 0.0, 0.0, halfpi64) isa Tuple{Float64,Float64}
        @inferred compute_native_pole(0.0, 0.0, halfpi64, 0.0, 0.0)

        # Branch 4: normal path
        @test compute_native_pole(0.0, deg2rad(45.0), 0.0, deg2rad(45.0), 0.0, halfpi64) isa Tuple{Float64,Float64}
        @inferred compute_native_pole(0.0, deg2rad(45.0), 0.0, deg2rad(45.0), 0.0)

        # ── Float32 ──────────────────────────────────────────
        # Branch 1: zenithal
        @test compute_native_pole(0f0, 0f0, 0f0, halfpi32, 0f0, halfpi32) isa Tuple{Float32,Float32}
        @inferred compute_native_pole(0f0, 0f0, 0f0, halfpi32, 0f0)

        # Branch 2: degenerate R
        @test compute_native_pole(0f0, 0f0, halfpi32, 0f0, 0f0, halfpi32) isa Tuple{Float32,Float32}
        @inferred compute_native_pole(0f0, 0f0, halfpi32, 0f0, 0f0)

        # Branch 4: normal path
        @test compute_native_pole(0f0, 1f0, 0f0, 1f0, 0f0, halfpi32) isa Tuple{Float32,Float32}
        @inferred compute_native_pole(0f0, 1f0, 0f0, 1f0, 0f0)

        # ── Default latpole with Float32 inputs ──────────────
        # Float64 default arg must not infect the return type.
        ap, dp = compute_native_pole(0f0, 1f0, 0f0, 1f0, 0f0)
        @test ap isa Float32
        @test dp isa Float32

        # ── Int inputs → Float64 promotion ───────────────────
        # Standard Julia promotion: float(Int) = Float64.
        @test compute_native_pole(0, 0, 0, 90, 0, 90) isa Tuple{Float64,Float64}
        @inferred compute_native_pole(0, 0, 0, 90, 0)

        # ── _reduce_lat ──────────────────────────────────────
        @inferred FITSWCS._reduce_lat(1.0)
        @inferred FITSWCS._reduce_lat(1.0f0)
        @test FITSWCS._reduce_lat(1.0) isa Float64
        @test FITSWCS._reduce_lat(1.0f0) isa Float32
        @test FITSWCS._reduce_lat(1) isa Float64
    end

    # ──────────────────────────────────────────────────────────────────────────
    # TPV/TPD distortion
    # ──────────────────────────────────────────────────────────────────────────

    @testset "TPD term table" begin
        terms = FITSWCS._TPD_TERMS
        @test length(terms) == 60
        # Spot-check the table entries (m is 0-based, Julia index is m+1).
        @test terms[1]  == (:mono, 0, 0)   # m=0: constant
        @test terms[2]  == (:mono, 1, 0)   # m=1: x
        @test terms[3]  == (:mono, 0, 1)   # m=2: y
        @test terms[4]  == (:radial, 1)    # m=3: r
        @test terms[5]  == (:mono, 2, 0)   # m=4: x²
        @test terms[6]  == (:mono, 1, 1)   # m=5: xy
        @test terms[7]  == (:mono, 0, 2)   # m=6: y²
        @test terms[8]  == (:mono, 3, 0)   # m=7: x³
        @test terms[9]  == (:mono, 2, 1)   # m=8: x²y
        @test terms[10] == (:mono, 1, 2)   # m=9: xy²
        @test terms[11] == (:mono, 0, 3)   # m=10: y³
        @test terms[12] == (:radial, 3)    # m=11: r³
    end

    @testset "TPV polynomial evaluation" begin
        # Identity coefficients: should return input unchanged.
        xcoeff = Float64[0.0, 1.0]          # x' = x
        ycoeff = Float64[0.0, 0.0, 1.0]     # y' = y
        @test FITSWCS._evaluate_tpv_polynomial(xcoeff, 1.5, 2.5) ≈ 1.5
        @test FITSWCS._evaluate_tpv_polynomial(ycoeff, 1.5, 2.5) ≈ 2.5

        # Quadratic terms: xcoeff[m=4] = x² coefficient.
        xc2 = Float64[0.0, 1.0, 0.0, 0.0, 1e-6]  # m=0..4: x' = x + 1e-6·x²
        yc2 = Float64[0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1e-6]  # m=0..6: y' = y + 1e-6·y²
        @test FITSWCS._evaluate_tpv_polynomial(xc2, 2.0, 0.0) ≈ 2.0 + 4e-6
        @test FITSWCS._evaluate_tpv_polynomial(yc2, 0.0, 3.0) ≈ 3.0 + 9e-6

        # Gapped coefficients: skip m indices with zero fill.
        xc_gap = zeros(Float64, 20)
        xc_gap[2] = 1.0           # m=1: x
        xc_gap[12] = 1e-6         # m=11: r³ = (x²+y²)^(3/2)
        val = FITSWCS._evaluate_tpv_polynomial(xc_gap, 3.0, 4.0)
        @test val ≈ 3.0 + 1e-6 * 125.0  # r³ = 5³ = 125
    end

    @testset "TPV inverse" begin
        xcoeff = Float64[0.0, 1.0]
        ycoeff = Float64[0.0, 0.0, 1.0]
        # Identity: inverse recovers target.
        u, v = FITSWCS._tpv_inverse(xcoeff, ycoeff, 5.0, -3.0)
        @test u ≈ 5.0
        @test v ≈ -3.0

        # Cubic distortion on x: x' = x + 1e-5·x³
        xc = zeros(Float64, 8)
        xc[2] = 1.0      # m=1: x
        xc[8] = 1e-5     # m=7: x³
        yc = Float64[0.0, 0.0, 1.0]
        target = 2.0 + 1e-5 * 8.0  # 2 + 1e-5*8
        u, v = FITSWCS._tpv_inverse(xc, yc, target, 0.0)
        @test u ≈ 2.0  atol=1e-10
        @test v ≈ 0.0  atol=1e-10

        # Divergent polynomial (negative identity: x' = -x).
        # Fixed-point iteration u_{k+1} = u_k - (-u_k - target) = 2*u_k + target
        # diverges monotonically.  Should warn and return the best estimate.
        xc_div = Float64[0.0, -1.0]  # x' = -x
        @test_warn "diverging" FITSWCS._tpv_inverse(xc_div, yc, 10.0, 0.0)
    end

    @testset "TPV projection round-trip" begin
        # Identity TPV matches TAN.
        tpv0 = TPV()
        @test FITSWCS.native_theta0(tpv0) == 90.0
        @test FITSWCS.native_phi0(tpv0) == 0.0

        # Round-trip: intermediate → native → intermediate.
        x, y = native_to_intermediate(tpv0, 0.3, deg2rad(80.0))
        phi, theta = intermediate_to_native(tpv0, x, y)
        @test phi ≈ 0.3 atol=1e-10
        @test theta ≈ deg2rad(80.0) atol=1e-10

        # Non-trivial TPV: x² term on x, y² term on y.
        xc = Float64[0.0, 1.0, 0.0, 0.0, 5e-6]   # x' = x + 5e-6·x²
        yc = Float64[0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -3e-6]  # y' = y - 3e-6·y²
        tpv1 = TPV(xc, yc)
        for _ in 1:5
            phi = (rand() - 0.5) * 0.1
            theta = deg2rad(75.0 + 10.0 * rand())
            x, y = native_to_intermediate(tpv1, phi, theta)
            phi2, theta2 = intermediate_to_native(tpv1, x, y)
            @test phi2 ≈ phi atol=1e-10
            @test theta2 ≈ theta atol=1e-10
        end
    end

    @testset "TPV header parsing" begin
        # Minimal TPV header (no PV keywords → identity).
        hdr = Dict{String,Any}(
            "NAXIS" => 2,
            "CTYPE1" => "RA---TPV", "CTYPE2" => "DEC--TPV",
            "CRPIX1" => 256.0, "CRPIX2" => 256.0,
            "CRVAL1" => 83.8221, "CRVAL2" => -5.3911,
            "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4,
        )
        wcs = WCS(hdr)
        @test wcs.projection isa TPV
        # Should behave like TAN with no extra distortion.
        pix = [256.0, 256.0]
        world = pixel_to_world(wcs, pix)
        @test world[1] ≈ 83.8221  atol=1e-4
        @test world[2] ≈ -5.3911  atol=1e-4

        # TPV with PV coefficients.
        hdr2 = Dict{String,Any}(
            "NAXIS" => 2,
            "CTYPE1" => "RA---TPV", "CTYPE2" => "DEC--TPV",
            "CRPIX1" => 128.0, "CRPIX2" => 96.0,
            "CRVAL1" => 120.0, "CRVAL2" => 35.0,
            "CDELT1" => -0.05, "CDELT2" => 0.05,
            "PV1_1" => 1.0,     # identity x term
            "PV1_4" => 1e-6,    # x² term on axis 1
            "PV2_2" => 1.0,     # identity y term
            "PV2_6" => -2e-6,   # y² term on axis 2
        )
        wcs2 = WCS(hdr2)
        t = wcs2.projection
        @test t isa TPV
        @test t.xcoeff[1] == 0.0     # m=0 (constant), zero-filled gap
        @test t.xcoeff[2] == 1.0     # m=1 (x)
        @test t.xcoeff[5] == 1e-6    # m=4 (x²)
        @test t.ycoeff[3] == 1.0     # m=2 (y)
        @test t.ycoeff[7] == -2e-6   # m=6 (y²)

        # Round-trip at reference pixel.
        pix2 = [128.0, 96.0]
        world2 = pixel_to_world(wcs2, pix2)
        pix2b = world_to_pixel(wcs2, world2)
        @test pix2b ≈ pix2  atol=1e-10

        # TPD CTYPE also works.
        hdr3 = Dict{String,Any}(
            "NAXIS" => 2,
            "CTYPE1" => "RA---TPD", "CTYPE2" => "DEC--TPD",
            "CRPIX1" => 256.0, "CRPIX2" => 256.0,
            "CRVAL1" => 45.0, "CRVAL2" => 30.0,
            "CDELT1" => -0.01, "CDELT2" => 0.01,
        )
        wcs3 = WCS(hdr3)
        @test wcs3.projection isa TPV
    end

    @testset "SCAMP compatibility" begin
        # Pre-2012 SCAMP: CTYPE=-TAN with PVi_j (j≥5) → -TPV.
        hdr_pre = Dict{String,Any}(
            "NAXIS" => 2,
            "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN",
            "CRPIX1" => 256.0, "CRPIX2" => 256.0,
            "CRVAL1" => 83.8221, "CRVAL2" => -5.3911,
            "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4,
            "PV1_1" => 1.0,     # TPV coefficients hiding as TAN
            "PV1_5" => 1e-6,    # j=5 signals SCAMP
            "PV2_2" => 1.0,
        )
        wcs_pre = WCS(hdr_pre)
        @test wcs_pre.projection isa TPV

        # TAN with only low-index PV keywords (j < 5) is NOT misdetected.
        hdr_tan_low = Dict{String,Any}(
            "NAXIS" => 2,
            "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN",
            "CRPIX1" => 256.0, "CRPIX2" => 256.0,
            "CRVAL1" => 83.8221, "CRVAL2" => -5.3911,
            "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4,
            "PV1_1" => 0.0,    # j=1 — low, should not trigger SCAMP
            "PV1_2" => 0.0,    # j=2
            "PV2_3" => 0.0,    # j=3
            "PV2_4" => 0.0,    # j=4 — highest sub-threshold index
        )
        wcs_tan_low = WCS(hdr_tan_low)
        @test wcs_tan_low.projection isa TAN

        # TPV + SIP: SIP keywords stripped (SCAMP rule).
        hdr_sip = Dict{String,Any}(
            "NAXIS" => 2,
            "CTYPE1" => "RA---TPV", "CTYPE2" => "DEC--TPV",
            "CRPIX1" => 256.0, "CRPIX2" => 256.0,
            "CRVAL1" => 83.8221, "CRVAL2" => -5.3911,
            "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4,
            "A_ORDER" => 2, "B_ORDER" => 2,
            "A_0_1" => 1e-5, "B_0_1" => -2e-5,
        )
        wcs_sip = WCS(hdr_sip)
        @test wcs_sip.projection isa TPV
        @test wcs_sip.pipeline isa NoDistortionPipeline   # SIP was stripped
    end

    @testset "TPV Float32 type stability" begin
        # Polynomial eval in Float32.
        xc32 = Float32[0.0, 1.0]
        yc32 = Float32[0.0, 0.0, 1.0]
        @test @inferred(FITSWCS._evaluate_tpv_polynomial(xc32, 1.0f0, 2.0f0)) isa Float32

        # Inverse in Float32.
        u, v = @inferred FITSWCS._tpv_inverse(xc32, yc32, 3.0f0, 4.0f0)
        @test u isa Float32
        @test v isa Float32

        # Projection functions in Float32.
        tpv32 = TPV(xc32, yc32)
        phi, theta = @inferred intermediate_to_native(tpv32, 0.0f0, 0.0f0)
        @test phi isa Float32
        @test theta isa Float32
        x32, y32 = @inferred native_to_intermediate(tpv32, 0.5f0, Float32(π/3))
        @test x32 isa Float32
        @test y32 isa Float32

        # Full-pipeline type preservation with Float32 pixel input.
        # Forward pass should preserve Float32 through the pipeline.
        hdr32 = Dict{String,Any}(
            "NAXIS" => 2,
            "CTYPE1" => "RA---TPV", "CTYPE2" => "DEC--TPV",
            "CRPIX1" => 256.0, "CRPIX2" => 256.0,
            "CRVAL1" => 83.8221, "CRVAL2" => -5.3911,
            "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4,
        )
        wcs32 = WCS(hdr32)
        pix32 = Float32[256.0, 256.0]
        world32 = pixel_to_world(wcs32, pix32)
        @test eltype(world32) == Float32
        # Float32 forward → Float64 world → Float64 pixel is consistent.
        world_copy = Float64[world32[1], world32[2]]
        pix_from_copy = world_to_pixel(wcs32, world_copy)
        world2 = pixel_to_world(wcs32, pix_from_copy)
        @test world2 ≈ world_copy  atol=1e-7
    end

    @testset "TPV full pipeline round-trip" begin
        # Build a TPV WCS with a moderate cubic distortion and verify that
        # pixel→world→pixel round-trips converge at multiple positions.
        hdr = Dict{String,Any}(
            "NAXIS" => 2,
            "CTYPE1" => "RA---TPV", "CTYPE2" => "DEC--TPV",
            "CRPIX1" => 256.0, "CRPIX2" => 256.0,
            "CRVAL1" => 83.8221, "CRVAL2" => -5.3911,
            "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4,
            "PV1_1" => 1.0,     # identity x
            "PV1_4" => 1e-7,    # small x² correction
            "PV2_2" => 1.0,     # identity y
            "PV2_6" => -1e-7,   # small y² correction
        )
        wcs = WCS(hdr)
        @test wcs.projection isa TPV

        test_pixels = [
            [256.0, 256.0],   # reference pixel
            [300.0, 256.0],   # off-reference in x
            [256.0, 300.0],   # off-reference in y
            [200.0, 200.0],   # lower-left quadrant
            [350.0, 350.0],   # upper-right quadrant
            [150.0, 100.0],   # corner
        ]
        for pix0 in test_pixels
            world = pixel_to_world(wcs, pix0)
            pix1 = world_to_pixel(wcs, world)
            @test pix1 ≈ pix0  atol=1e-7
            # World round-trip error is bounded by the _tpv_inverse convergence
            # tolerance (1e-10 degrees).
            pix2 = world_to_pixel(wcs, world)
            world2 = pixel_to_world(wcs, pix2)
            @test world2 ≈ world  atol=1e-10
        end
    end
end  # @testset "FITSWCS"

# Regression tests against wcslib reference values (Milestone 5).
include("regression_wcslib.jl")

# Regression tests against stored Astropy values for projections without wcslib fixtures.
include("regression_astropy_values.jl")
