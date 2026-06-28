using Test
using LinearAlgebra
using Random
using FITSWCS
using FITSFiles
using FITSIO
using FITSWCS: pixel_to_intermediate, intermediate_to_pixel,
               parse_ctype, projection_from_code,
               native_to_celestial, celestial_to_native,
               unit_to_deg, build_cd_matrix,
               evaluate_sip_polynomial, sip_pixel_to_focal, sip_focal_to_pixel

# Convenience shorthand
const D2R = π / 180.0
const R2D = 180.0 / π

"""Test that two angles (in radians) are equal modulo 2π."""
function angle_approx(a::Real, b::Real; atol=1e-12)
    d = mod(a - b + π, 2π) - π   # wrap difference to (-π, π]
    return abs(d) <= atol
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
        wcs = from_header(hdr)
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
        wcs = from_header(hdr)
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
        wcs = from_header(hdr)
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
        wcs = from_header(hdr)
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
        wcs = from_header(hdr)
        @test wcs.projection == CEA(0.75)

        bad = copy(hdr)
        bad["PV2_1"] = 0.0
        @test_throws ArgumentError from_header(bad)
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
        @test_throws ArgumentError from_header(hdr)
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
        wcs = from_header(hdr)
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
        wcs = from_header(hdr)
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
        wcs = from_header(hdr)
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
        wcs = from_header(hdr)
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
        wcs = from_header(hdr)
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
        @test_throws ArgumentError from_header(hdr)
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
        @test_throws ArgumentError from_header(hdr)
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
        wcs = from_header(hdr)
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
        wcs_car = from_header(hdr_car)
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
        wcs_car2 = from_header(hdr_car2)
        @test wcs_car2.lonpole == 0.0
    end

    @testset "Error: missing NAXIS/WCSAXES" begin
        @test_throws ArgumentError from_header(Dict("CRPIX1" => 1.0))
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

        primary = from_header(hdr)
        alternate = from_header(hdr; alt='A')

        @test pixel_to_world(primary, [2.0]) ≈ [1.001e9]
        @test pixel_to_world(alternate, [7.0]) ≈ [501.0]
        @test_throws ArgumentError from_header(hdr; alt='a')
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
        @test_throws ArgumentError from_header(hdr)
    end

    @testset "Error: TAB lookup axes are explicitly unsupported" begin
        # Paper III -TAB axes require table arrays, so they must not be treated as linear axes.
        hdr = Dict(
            "NAXIS"  => 1,
            "CTYPE1" => "FREQ-TAB",
            "CRPIX1" => 1.0,
            "CRVAL1" => 1.0,
            "CDELT1" => 1.0,
        )
        @test_throws ArgumentError from_header(hdr)

        alt_hdr = Dict(
            "NAXIS"   => 1,
            "CTYPE1"  => "FREQ",
            "CTYPE1A" => "WAVE-TAB",
        )
        @test pixel_to_world(from_header(alt_hdr), [2.0]) ≈ [2.0]
        @test_throws ArgumentError from_header(alt_hdr; alt='A')
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
        @test_throws ArgumentError from_header(hdr)

        alt_hdr = Dict(
            "NAXIS"   => 1,
            "CTYPE1"  => "FREQ",
            "CTYPE1A" => "WAVE-F2W",
            "CRPIX1A" => 1.0,
            "CRVAL1A" => 500.0,
            "CDELT1A" => 1.0,
        )
        @test pixel_to_world(from_header(alt_hdr), [2.0]) ≈ [2.0]
        @test_throws ArgumentError from_header(alt_hdr; alt='A')
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
        wcs = from_header(hdr)
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
        wcs = from_header(hdr)
        @test wcs.sip isa SIPDistortion
        @test wcs.sip.a[3, 1] == 1e-3
        @test wcs.sip.b[1, 3] == -2e-3
        @test wcs.projection isa TAN
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
        wcs = from_header(hdr)
        pix = [12.0, 23.0]
        focal = sip_pixel_to_focal(wcs.sip, pix)
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
        wcs = from_header(hdr)
        @test sip_focal_to_pixel(wcs.sip, [10.0, 5.0]) ≈ [9.0, 6.0]
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
        wcs = from_header(hdr)
        pix = [103.0, 98.0]
        focal = sip_pixel_to_focal(wcs.sip, pix)
        @test sip_focal_to_pixel(wcs.sip, focal) ≈ pix atol=1e-8
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
        wcs = from_header(hdr)
        pix = [12.0, 23.0]
        world = pixel_to_world(wcs, pix)
        @test world ≈ [104.008, 208.973]
        @test world_to_pixel(wcs, world) ≈ pix atol=1e-8
    end

    @testset "Malformed SIP headers throw clear errors" begin
        # SIP requires explicit CRPIX and matched forward/inverse order pairs.
        @test_throws ArgumentError from_header(Dict(
            "NAXIS" => 2,
            "A_ORDER" => 2,
            "B_ORDER" => 2,
        ))
        @test_throws ArgumentError from_header(Dict(
            "NAXIS" => 2,
            "CRPIX1" => 0.0, "CRPIX2" => 0.0,
            "A_ORDER" => 2,
        ))
        @test_throws ArgumentError from_header(Dict(
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
            @test_throws ArgumentError from_header(hdr)
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
            @test from_header(hdr) isa WCSTransform
            @test_throws ArgumentError from_header(hdr; alt='A')
        end

        # SIP is implemented and should still parse through the distortion path.
        sip_hdr = copy(base)
        sip_hdr["CTYPE1"] = "RA---TAN-SIP"
        sip_hdr["CTYPE2"] = "DEC--TAN-SIP"
        sip_hdr["A_ORDER"] = 2
        sip_hdr["B_ORDER"] = 2
        @test from_header(sip_hdr).sip isa SIPDistortion
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
        wcs = from_header(hdr)

        @test wcs.lon_axis == 1
        @test wcs.lat_axis == 3
        @test pixel_to_world(wcs, [30.0, 40.0, 45.0]) ≈ [10.0, 1.42e9, 25.0] atol=1e-12
        @test pixel_to_world(wcs, [30.0, 43.0, 45.0]) ≈ [10.0, 1.423e9, 25.0] atol=1e-12

        for pix in ([29.0, 43.0, 44.0], [31.5, 37.0, 47.0])
            world = pixel_to_world(wcs, pix)
            @test world[2] ≈ 1.42e9 + 1.0e6 * (pix[2] - 40.0)
            @test world_to_pixel(wcs, world) ≈ pix atol=1e-8
        end
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
        wcs = from_header(hdr)
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
        wcs = from_header(hdr)
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
    @test world_to_pix(wcs, [2.0, 6.0]) ≈ [2.0, 3.0]
    @test world_to_pix(wcs, 2.0, 6.0) ≈ [2.0, 3.0]
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
        wcs = from_header(hdr)
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
        wcs = from_header(hdr)
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
        wcs = from_header(hdr)
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
        wcs = from_header(hdr)
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
        wcs = from_header(hdr)
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
        wcs = from_header(hdr)
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
        wcs = from_header(hdr)
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
    wcs = from_header(hdr)

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
    wcs = from_header(hdr)

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
    wcs = from_header(hdr)
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
    wcs = from_header(hdr)
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
    wcs = from_header(hdr)
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
    proj = UnknownProjection("TPV")
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

    wcs = from_header(header)
    @test pixel_to_world(wcs, [3.0]) ≈ [104.0]

    mktempdir() do dir
        path = joinpath(dir, "linear.fits")
        FITSIO.FITS(path, "w") do file
            FITSIO.write(file, zeros(Float32, 4); header=header)
            hdu_wcs = from_header(file[1])
            @test pixel_to_world(hdu_wcs, [3.0]) ≈ [104.0]
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

    wcs = from_header(cards)
    @test pixel_to_world(wcs, [3.0]) ≈ [104.0]

    hdu = FITSFiles.HDU(cards)
    hdu_wcs = from_header(hdu)
    @test pixel_to_world(hdu_wcs, [3.0]) ≈ [104.0]
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

    wcs_lin = from_header(hdr_lin)
    wcs_tan = from_header(hdr_tan)
    wcs_sip = from_header(hdr_sip)
    wcs_3d  = from_header(hdr_3d)

    # from_header return type is stable.
    @test @inferred(from_header(hdr_lin)) isa WCSTransform
    @test @inferred(from_header(hdr_tan)) isa WCSTransform
    @test @inferred(from_header(hdr_sip)) isa WCSTransform

    # pixel_to_world return type is stable across linear, celestial, SIP, and 3D paths.
    @test @inferred(pixel_to_world(wcs_lin, [2.0, 3.0]))  isa Vector{Float64}
    @test @inferred(pixel_to_world(wcs_tan, [400.0, 300.0])) isa Vector{Float64}
    @test @inferred(pixel_to_world(wcs_sip, [600.0, 500.0])) isa Vector{Float64}
    @test @inferred(pixel_to_world(wcs_3d,  [40.0, 60.0, 5.0])) isa Vector{Float64}

    # world_to_pixel return type is stable across the same paths.
    w_lin = pixel_to_world(wcs_lin, [2.0, 3.0])
    w_tan = pixel_to_world(wcs_tan, [400.0, 300.0])
    w_sip = pixel_to_world(wcs_sip, [600.0, 500.0])
    w_3d  = pixel_to_world(wcs_3d,  [40.0, 60.0, 5.0])
    @test @inferred(world_to_pixel(wcs_lin, w_lin)) isa Vector{Float64}
    @test @inferred(world_to_pixel(wcs_tan, w_tan)) isa Vector{Float64}
    @test @inferred(world_to_pixel(wcs_sip, w_sip)) isa Vector{Float64}
    @test @inferred(world_to_pixel(wcs_3d,  w_3d))  isa Vector{Float64}

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
            wcs = from_header(hdr)

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
            wcs = from_header(hdr)

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
            wcs = from_header(hdr)

            # Pixel near the reference point.
            pix = [256.0 + 10.0*randn(rng), 256.0 + 10.0*randn(rng)]
            world = pixel_to_world(wcs, pix)
            @test world_to_pixel(wcs, world) ≈ pix  atol=1e-5
        end
    end
end

end  # @testset "FITSWCS"

# Regression tests against wcslib reference values (Milestone 5).
include("regression_wcslib.jl")

# Regression tests against stored Astropy values for projections without wcslib fixtures.
include("regression_astropy_values.jl")
