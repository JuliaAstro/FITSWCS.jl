"""
Reference-comparison regression tests against WCS.jl/wcslib.

Expected values were generated using WCS.jl v0.6.3 (wrapping wcslib), which is
the authoritative reference implementation for FITS WCS.

Each test asserts that FITSWCS.jl pixel_to_world and world_to_pixel produce
results consistent with wcslib to within floating-point tolerances.

Tolerances are tight (1e-8 degrees ≈ 36 µas) unless the geometry near a
singularity requires otherwise.
"""

using Test
using FITSWCS

# ---------------------------------------------------------------------------
# TAN projection – CDELT form (Orion Nebula-style)
# ctype = ["RA---TAN", "DEC--TAN"]
# crpix = [512.0, 512.0]
# crval = [83.8221, -5.3911]
# cdelt = [-2.7778e-4, 2.7778e-4]
#
# Reference values generated with WCS.jl 0.6.3 / wcslib.
# ---------------------------------------------------------------------------

@testset "TAN CDELT-form (wcslib comparison)" begin
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---TAN",  "CTYPE2" => "DEC--TAN",
        "CRPIX1" => 512.0,       "CRPIX2" => 512.0,
        "CRVAL1" => 83.8221,     "CRVAL2" => -5.3911,
        "CDELT1" => -2.7778e-4,  "CDELT2" =>  2.7778e-4,
    )
    wcs = WCS(hdr)

    # Reference pixel maps to CRVAL
    w = pixel_to_world(wcs, [512.0, 512.0])
    @test w[1] ≈ 83.8221         atol=1e-8
    @test w[2] ≈ -5.3911         atol=1e-8

    # Corners
    for (px, py, ra_ref, dec_ref) in [
        (  1.0,    1.0, 83.96470930312101, -5.533028256996351),
        (1024.0,   1.0, 83.67921161916392, -5.533028190267628),
        (  1.0, 1024.0, 83.96464257062250, -5.248860779320227),
        (1024.0,1024.0, 83.67927848225278, -5.248860716038349),
    ]
        w = pixel_to_world(wcs, [px, py])
        @test w[1] ≈ ra_ref  atol=1e-7
        @test w[2] ≈ dec_ref atol=1e-7
        # Round-trip
        p = world_to_pixel(wcs, w)
        @test p[1] ≈ px  atol=1e-6
        @test p[2] ≈ py  atol=1e-6
    end

    # Interior point
    w = pixel_to_world(wcs, [100.0, 200.0])
    @test w[1] ≈ 83.93707010775421  atol=1e-8
    @test w[2] ≈ -5.477756332944454 atol=1e-8
    p = world_to_pixel(wcs, w)
    @test p[1] ≈ 100.0  atol=1e-6
    @test p[2] ≈ 200.0  atol=1e-6
end

# ---------------------------------------------------------------------------
# TAN projection – PC matrix form (45° rotation, 1 arcsec/pixel)
# crpix = [256.0, 256.0]
# crval = [45.0, 20.0]
# cdelt = [-1/3600, 1/3600]
# pc    = [cos(45°), -sin(45°); sin(45°), cos(45°)]
# ---------------------------------------------------------------------------

@testset "TAN PC-matrix 45° rotation (wcslib comparison)" begin
    cdelt_v = 1.0 / 3600.0
    rho     = 45.0 * π / 180.0
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---TAN",   "CTYPE2" => "DEC--TAN",
        "CRPIX1" => 256.0,        "CRPIX2" => 256.0,
        "CRVAL1" => 45.0,         "CRVAL2" => 20.0,
        "CDELT1" => -cdelt_v,     "CDELT2" =>  cdelt_v,
        "PC1_1"  => cos(rho),     "PC1_2"  => -sin(rho),
        "PC2_1"  => sin(rho),     "PC2_2"  =>  cos(rho),
    )
    wcs = WCS(hdr)

    # Reference pixel
    w = pixel_to_world(wcs, [256.0, 256.0])
    @test w[1] ≈ 45.0  atol=1e-10
    @test w[2] ≈ 20.0  atol=1e-10

    # Off-reference pixels (wcslib reference values)
    w2 = pixel_to_world(wcs, [300.0, 256.0])
    @test w2[1] ≈ 44.99080242788927   atol=1e-8
    @test w2[2] ≈ 20.008642178799985  atol=1e-8

    w3 = pixel_to_world(wcs, [256.0, 300.0])
    @test w3[1] ≈ 45.00919757211073   atol=1e-8
    @test w3[2] ≈ 20.008642178799985  atol=1e-8

    # Round-trip
    for pix in ([256.0, 256.0], [300.0, 256.0], [256.0, 300.0], [200.0, 310.0])
        @test world_to_pixel(wcs, pixel_to_world(wcs, pix)) ≈ pix  atol=1e-7
    end
end

# ---------------------------------------------------------------------------
# TAN projection – CD matrix form (HST ACS/WFC-like)
# cd11 = -scale*cos(θ)  cd12 = -scale*sin(θ)
# cd21 = -scale*sin(θ)  cd22 =  scale*cos(θ)
# where scale = 0.05/3600 (0.05 arcsec/pixel), θ = 0.1°
# crpix = [2048.0, 1024.0], crval = [150.0, 2.5]
# ---------------------------------------------------------------------------

@testset "TAN CD-matrix HST ACS/WFC-like (wcslib comparison)" begin
    theta = 0.1 * π / 180.0
    scale = 0.05 / 3600.0
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---TAN",   "CTYPE2" => "DEC--TAN",
        "CRPIX1" => 2048.0,       "CRPIX2" => 1024.0,
        "CRVAL1" => 150.0,        "CRVAL2" => 2.5,
        "CD1_1"  => -scale * cos(theta),  "CD1_2" => -scale * sin(theta),
        "CD2_1"  => -scale * sin(theta),  "CD2_2" =>  scale * cos(theta),
    )
    wcs = WCS(hdr)

    # Reference pixel
    w = pixel_to_world(wcs, [2048.0, 1024.0])
    @test w[1] ≈ 150.0  atol=1e-10
    @test w[2] ≈   2.5  atol=1e-10

    # Corner comparison against wcslib
    w1 = pixel_to_world(wcs, [1.0, 1.0])
    @test w1[1] ≈ 150.02848210976424   atol=1e-7
    @test w1[2] ≈   2.485841002491482  atol=1e-7

    w2 = pixel_to_world(wcs, [4096.0, 2048.0])
    @test w2[1] ≈ 149.9715033488136   atol=1e-7
    @test w2[2] ≈   2.514172244812679 atol=1e-7

    # Round-trip at interior points
    for pix in ([2048.0, 1024.0], [1000.0, 500.0], [3000.0, 1500.0])
        @test world_to_pixel(wcs, pixel_to_world(wcs, pix)) ≈ pix  atol=1e-5
    end
end

# ---------------------------------------------------------------------------
# AIT projection – galactic all-sky
# ctype = ["GLON-AIT", "GLAT-AIT"]
# crpix = [360.5, 180.5], crval = [0.0, 0.0], cdelt = [-0.5, 0.5]
# ---------------------------------------------------------------------------

@testset "AIT galactic all-sky (wcslib comparison)" begin
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "GLON-AIT",  "CTYPE2" => "GLAT-AIT",
        "CRPIX1" => 360.5,       "CRPIX2" => 180.5,
        "CRVAL1" => 0.0,         "CRVAL2" => 0.0,
        "CDELT1" => -0.5,        "CDELT2" => 0.5,
    )
    wcs = WCS(hdr)

    # Center
    w = pixel_to_world(wcs, [360.5, 180.5])
    @test w[1] ≈ 0.0  atol=1e-8
    @test w[2] ≈ 0.0  atol=1e-8

    # Off-center against wcslib
    w2 = pixel_to_world(wcs, [300.0, 150.0])
    @test w2[1] ≈ 31.167480197179597   atol=1e-7
    @test w2[2] ≈ -15.155848828445626  atol=1e-7

    w3 = pixel_to_world(wcs, [400.0, 190.0])
    @test w3[1] ≈ 340.1744667748081   atol=1e-7
    @test w3[2] ≈   4.733615120590103 atol=1e-7

    # Equatorial row: latitude = 0°
    w4 = pixel_to_world(wcs, [100.0, 180.5])
    @test w4[1] ≈ 138.5334189628714  atol=1e-7
    @test w4[2] ≈ 0.0                atol=1e-7

    # Round-trip
    for pix in ([360.5, 180.5], [300.0, 150.0], [400.0, 190.0], [200.0, 160.0])
        @test world_to_pixel(wcs, pixel_to_world(wcs, pix)) ≈ pix  atol=1e-6
    end
end

# ---------------------------------------------------------------------------
# SIN projection
# ctype = ["RA---SIN", "DEC--SIN"]
# crpix = [100.0, 100.0], crval = [180.0, 30.0], cdelt = [-0.01, 0.01]
# ---------------------------------------------------------------------------

@testset "SIN projection (wcslib comparison)" begin
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---SIN",  "CTYPE2" => "DEC--SIN",
        "CRPIX1" => 100.0,       "CRPIX2" => 100.0,
        "CRVAL1" => 180.0,       "CRVAL2" => 30.0,
        "CDELT1" => -0.01,       "CDELT2" => 0.01,
    )
    wcs = WCS(hdr)

    # Reference pixel
    w = pixel_to_world(wcs, [100.0, 100.0])
    @test w[1] ≈ 180.0  atol=1e-10
    @test w[2] ≈  30.0  atol=1e-10

    # Off-reference against wcslib
    w2 = pixel_to_world(wcs, [150.0, 80.0])
    @test w2[1] ≈ 179.42380496785952  atol=1e-7
    @test w2[2] ≈  29.79874251298449  atol=1e-7

    w3 = pixel_to_world(wcs, [60.0, 120.0])
    @test w3[1] ≈ 180.46281699820287  atol=1e-7
    @test w3[2] ≈  30.199192628819894 atol=1e-7

    # Round-trip
    for pix in ([100.0, 100.0], [150.0, 80.0], [60.0, 120.0])
        @test world_to_pixel(wcs, pixel_to_world(wcs, pix)) ≈ pix  atol=1e-7
    end
end

# ---------------------------------------------------------------------------
# CAR projection
# ctype = ["RA---CAR", "DEC--CAR"]
# crpix = [181.0, 91.0], crval = [0.0, 0.0], cdelt = [-1.0, 1.0]
# ---------------------------------------------------------------------------

@testset "CAR projection (wcslib comparison)" begin
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---CAR",  "CTYPE2" => "DEC--CAR",
        "CRPIX1" => 181.0,       "CRPIX2" => 91.0,
        "CRVAL1" => 0.0,         "CRVAL2" => 0.0,
        "CDELT1" => -1.0,        "CDELT2" => 1.0,
    )
    wcs = WCS(hdr)

    # Reference pixel
    w = pixel_to_world(wcs, [181.0, 91.0])
    @test w[1] ≈ 0.0  atol=1e-10
    @test w[2] ≈ 0.0  atol=1e-10

    # Off-reference against wcslib
    w2 = pixel_to_world(wcs, [100.0, 50.0])
    @test w2[1] ≈ 81.0   atol=1e-10
    @test w2[2] ≈ -41.0  atol=1e-10

    w3 = pixel_to_world(wcs, [50.0, 20.0])
    @test w3[1] ≈ 131.0  atol=1e-10
    @test w3[2] ≈ -71.0  atol=1e-10

    # Round-trip uses the local longitude branch, matching Astropy/wcslib.
    for pix in ([181.0, 91.0], [100.0, 50.0], [50.0, 20.0], [160.0, 100.0], [200.0, 100.0])
        @test world_to_pixel(wcs, pixel_to_world(wcs, pix)) ≈ pix  atol=1e-8
    end
end

# ---------------------------------------------------------------------------
# STG projection
# ctype = ["RA---STG", "DEC--STG"]
# crpix = [64.0, 64.0], crval = [270.0, -45.0], cdelt = [-0.1, 0.1]
# ---------------------------------------------------------------------------

@testset "STG projection (wcslib comparison)" begin
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---STG",  "CTYPE2" => "DEC--STG",
        "CRPIX1" => 64.0,        "CRPIX2" => 64.0,
        "CRVAL1" => 270.0,       "CRVAL2" => -45.0,
        "CDELT1" => -0.1,        "CDELT2" => 0.1,
    )
    wcs = WCS(hdr)

    # Reference pixel
    w = pixel_to_world(wcs, [64.0, 64.0])
    @test w[1] ≈ 270.0  atol=1e-10
    @test w[2] ≈ -45.0  atol=1e-10

    # Off-reference against wcslib
    w2 = pixel_to_world(wcs, [70.0, 64.0])
    @test w2[1] ≈ 269.1515106303453    atol=1e-7
    @test w2[2] ≈ -44.996858579590125  atol=1e-7

    w3 = pixel_to_world(wcs, [64.0, 70.0])
    @test w3[1] ≈ 270.0                atol=1e-8
    @test w3[2] ≈ -44.400005483023364  atol=1e-7

    # Round-trip
    for pix in ([64.0, 64.0], [70.0, 64.0], [64.0, 70.0], [58.0, 58.0])
        @test world_to_pixel(wcs, pixel_to_world(wcs, pix)) ≈ pix  atol=1e-7
    end
end

# ---------------------------------------------------------------------------
# 3D cube – celestial + spectral
# ctype = ["RA---TAN", "DEC--TAN", "FREQ"]
# crpix = [50.0, 50.0, 1.0], crval = [10.0, 25.0, 1.42e9]
# cdelt = [-0.01, 0.01, 1.0e6]
# ---------------------------------------------------------------------------

@testset "3D cube RA+DEC+FREQ (wcslib comparison)" begin
    hdr = Dict(
        "NAXIS"  => 3,
        "CTYPE1" => "RA---TAN",  "CTYPE2" => "DEC--TAN",  "CTYPE3" => "FREQ",
        "CRPIX1" => 50.0,        "CRPIX2" => 50.0,         "CRPIX3" => 1.0,
        "CRVAL1" => 10.0,        "CRVAL2" => 25.0,          "CRVAL3" => 1.42e9,
        "CDELT1" => -0.01,       "CDELT2" =>  0.01,         "CDELT3" => 1.0e6,
    )
    wcs = WCS(hdr)

    # Reference pixel
    w = pixel_to_world(wcs, [50.0, 50.0, 1.0])
    @test w[1] ≈ 10.0    atol=1e-10
    @test w[2] ≈ 25.0    atol=1e-10
    @test w[3] ≈ 1.42e9  atol=1e-3

    # Off-reference against wcslib
    w2 = pixel_to_world(wcs, [60.0, 40.0, 5.0])
    @test w2[1] ≈ 9.8897520707026      atol=1e-7
    @test w2[2] ≈ 24.89995959401809    atol=1e-7
    @test w2[3] ≈ 1.424e9              atol=1e3

    # Round-trip
    for pix in ([50.0, 50.0, 1.0], [60.0, 40.0, 5.0], [45.0, 55.0, 10.0])
        @test world_to_pixel(wcs, pixel_to_world(wcs, pix)) ≈ pix  atol=1e-7
    end
end

# ---------------------------------------------------------------------------
# SIP distortion – TAN-SIP with 2nd-order A and B polynomials
# ctype = ["RA---TAN-SIP", "DEC--TAN-SIP"]
# crpix = [512.0, 512.0], crval = [150.0, 2.5], cdelt = [-2.7778e-4, 2.7778e-4]
# A_ORDER = 2, A_2_0 = 5e-6, A_0_2 = 2e-6
# B_ORDER = 2, B_2_0 = 1e-6, B_1_1 = 3e-6
#
# Reference values generated with WCS.jl 0.6.3 / wcslib.
# ---------------------------------------------------------------------------

@testset "TAN-SIP distortion (wcslib comparison)" begin
    hdr = Dict(
        "NAXIS"   => 2,
        "CTYPE1"  => "RA---TAN-SIP",  "CTYPE2"  => "DEC--TAN-SIP",
        "CRPIX1"  => 512.0,            "CRPIX2"  => 512.0,
        "CRVAL1"  => 150.0,            "CRVAL2"  => 2.5,
        "CDELT1"  => -2.7778e-4,       "CDELT2"  =>  2.7778e-4,
        "A_ORDER" => 2,
        "A_2_0"   => 5.0e-6,           "A_0_2"   => 2.0e-6,   "A_1_1" => 0.0,
        "B_ORDER" => 2,
        "B_2_0"   => 1.0e-6,           "B_0_2"   => 0.0,      "B_1_1" => 3.0e-6,
    )
    wcs = WCS(hdr)

    # Reference pixel maps to CRVAL
    w = pixel_to_world(wcs, [512.0, 512.0])
    @test w[1] ≈ 150.0  atol=1e-9
    @test w[2] ≈   2.5  atol=1e-9

    # Off-reference against wcslib
    w2 = pixel_to_world(wcs, [600.0, 500.0])
    @test w2[1] ≈ 149.97552128963324  atol=1e-7
    @test w2[2] ≈   2.4966676835562676 atol=1e-7

    w3 = pixel_to_world(wcs, [400.0, 600.0])
    @test w3[1] ≈ 150.03111983052793  atol=1e-7
    @test w3[2] ≈   2.524439537711708  atol=1e-7

    w4 = pixel_to_world(wcs, [300.0, 300.0])
    @test w4[1] ≈ 150.05885532834      atol=1e-7
    @test w4[2] ≈   2.4411593124884754 atol=1e-7

    # Round-trip (uses iterative inverse or AP/BP if present)
    for pix in ([512.0, 512.0], [600.0, 500.0], [400.0, 600.0], [300.0, 300.0])
        @test world_to_pixel(wcs, pixel_to_world(wcs, pix)) ≈ pix  atol=1e-5
    end

    # Verify SIP path differs from no-SIP path for the same base transform
    hdr_nosip = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---TAN",  "CTYPE2" => "DEC--TAN",
        "CRPIX1" => 512.0,       "CRPIX2" => 512.0,
        "CRVAL1" => 150.0,       "CRVAL2" => 2.5,
        "CDELT1" => -2.7778e-4,  "CDELT2" =>  2.7778e-4,
    )
    wcs_nosip = WCS(hdr_nosip)
    # At the reference pixel SIP adds zero distortion, so results match
    @test pixel_to_world(wcs, [512.0, 512.0]) ≈ pixel_to_world(wcs_nosip, [512.0, 512.0])  atol=1e-10
    # Off-center the SIP result must differ from the non-SIP result
    w_sip   = pixel_to_world(wcs,       [600.0, 500.0])
    w_nosip = pixel_to_world(wcs_nosip, [600.0, 500.0])
    @test !isapprox(w_sip[1], w_nosip[1]; atol=1e-8)
end
