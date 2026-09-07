@testset "TPV sequent distortion ordering" begin
    # A nontrivial PC/CDELT split must match WCSLIB's PC -> TPV -> CDELT pipeline.
    header = Dict{String,Any}(
        "NAXIS" => 2,
        "CTYPE1" => "RA---TPV", "CTYPE2" => "DEC--TPV",
        "CRPIX1" => 100.0, "CRPIX2" => 200.0,
        "CRVAL1" => 123.4, "CRVAL2" => -22.5,
        "CDELT1" => -0.02, "CDELT2" => 0.035,
        "PC1_1" => 0.8, "PC1_2" => 0.3,
        "PC2_1" => -0.2, "PC2_2" => 1.1,
        "PV1_0" => 0.01, "PV1_1" => 1.02, "PV1_2" => -0.03,
        "PV1_4" => 0.0007, "PV1_5" => -0.0002, "PV1_6" => 0.0005,
        "PV2_0" => -0.02, "PV2_1" => 0.5, "PV2_2" => 0.5,
        "PV2_4" => 0.0004, "PV2_5" => -0.0003, "PV2_6" => 0.0004,
    )
    wcs = WCS(header)
    pixel = [112.3, 176.4]
    expected_intermediate = [-0.08205336719999999, -0.43751199240000005]
    expected_world = [123.31090434225322, -22.937478626425648]
    @test pixel_to_intermediate(wcs, pixel) ≈ expected_intermediate atol=2e-15
    @test pixel_to_world(wcs, pixel) ≈ expected_world atol=2e-14
    @test world_to_pixel(wcs, expected_world) ≈ pixel atol=1e-11

    # Explicit CD is interpreted as PC with unit CDELT, as WCSLIB requires for TPV.
    cd_header = copy(header)
    for key in ("CDELT1", "CDELT2", "PC1_1", "PC1_2", "PC2_1", "PC2_2")
        delete!(cd_header, key)
    end
    cd_header["CD1_1"] = -0.016
    cd_header["CD1_2"] = -0.006
    cd_header["CD2_1"] = -0.007
    cd_header["CD2_2"] = 0.0385
    cd_wcs = WCS(cd_header)
    @test pixel_to_intermediate(cd_wcs, pixel) ≈ [-0.01597713451499998, -0.5445694821799998] atol=2e-15
    @test pixel_to_world(cd_wcs, pixel) ≈ [123.38263812267824, -23.0445521374377] atol=2e-14
end

@testset "TPV axis-local coefficients" begin
    # Standard identity coefficients use term 1 on both attached axes.
    header = Dict{String,Any}(
        "NAXIS" => 2,
        "CTYPE1" => "RA---TPV", "CTYPE2" => "DEC--TPV",
        "CRPIX1" => 0.0, "CRPIX2" => 0.0,
        "CRVAL1" => 0.0, "CRVAL2" => 0.0,
        "CDELT1" => 1.0, "CDELT2" => 1.0,
        "PV1_1" => 1.0, "PV2_1" => 1.0,
    )
    wcs = WCS(header)
    @test pixel_to_intermediate(wcs, [2.0, 3.0]) == [2.0, 3.0]
    @test pixel_to_world(wcs, [2.0, 3.0]) ≈ [1.9991882802315373, 2.9954418975929276] atol=2e-14
    @test world_to_pixel(wcs, pixel_to_world(wcs, [2.0, 3.0])) ≈ [2.0, 3.0] atol=2e-14

    # TPV has only the first 40 TPD terms.
    header["PV1_40"] = 1.0
    @test_throws "parameter numbers 0 through 39" WCS(header)
end

@testset "Paper IV prior TPD" begin
    # Prior additive TPD must run before the linear WCS matrix and invert cleanly.
    header = Dict{String,Any}(
        "NAXIS" => 2, "CTYPE1" => "LINEAR", "CTYPE2" => "LINEAR",
        "CRPIX1" => 0.0, "CRPIX2" => 0.0,
        "CRVAL1" => 0.0, "CRVAL2" => 0.0,
        "CDELT1" => 1.0, "CDELT2" => 1.0,
        "CPDIS1" => "TPD", "CPDIS2" => "TPD",
        "DP1.NAXES" => 2, "DP1.AXIS.1" => 1, "DP1.AXIS.2" => 2,
        "DP1.TPD.FWD.4" => 0.01,
        "DP2.NAXES" => 2, "DP2.AXIS.1" => 2, "DP2.AXIS.2" => 1,
        "DP2.TPD.FWD.4" => -0.02,
    )
    wcs = WCS(header)
    @test pixel_to_world(wcs, [2.0, 3.0]) ≈ [2.04, 2.82] atol=2e-15
    @test world_to_pixel(wcs, [2.04, 2.82]) ≈ [2.0, 3.0] atol=1e-13

    # This nonzero-CRPIX case verifies WCSLIB's prior TPD -> CRPIX -> PC -> CDELT order.
    rich = Dict{String,Any}(
        "NAXIS" => 2, "CTYPE1" => "LINEAR", "CTYPE2" => "LINEAR",
        "CRPIX1" => 10.0, "CRPIX2" => 20.0,
        "CRVAL1" => 5.0, "CRVAL2" => -7.0,
        "CDELT1" => -0.02, "CDELT2" => 0.035,
        "PC1_1" => 0.8, "PC1_2" => 0.3,
        "PC2_1" => -0.2, "PC2_2" => 1.1,
        "CPDIS1" => "TPD", "CPDIS2" => "TPD",
        "DP1.NAXES" => 2, "DP1.AXIS.1" => 1, "DP1.AXIS.2" => 2,
        "DP1.OFFSET.1" => 1.0, "DP1.OFFSET.2" => -2.0,
        "DP1.SCALE.1" => 0.5, "DP1.SCALE.2" => 2.0,
        "DP1.TPD.FWD.0" => 0.01, "DP1.TPD.FWD.1" => 1.02, "DP1.TPD.FWD.2" => -0.03,
        "DP1.TPD.FWD.4" => 0.0007, "DP1.TPD.FWD.5" => -0.0002, "DP1.TPD.FWD.6" => 0.0005,
        "DP2.NAXES" => 2, "DP2.AXIS.1" => 2, "DP2.AXIS.2" => 1,
        "DP2.TPD.FWD.4" => -0.002, "DP2.TPD.FWD.5" => 0.0003,
    )
    rich_wcs = WCS(rich)
    pixel = [12.3, -3.4]
    expected_world = [5.00963112, -7.95957814925]
    @test pixel_to_intermediate(rich_wcs, pixel) ≈ [0.009631119999999993, -0.9595781492500003] atol=2e-15
    @test pixel_to_world(rich_wcs, pixel) ≈ expected_world atol=2e-14
    @test world_to_pixel(rich_wcs, expected_world) ≈ pixel atol=2e-12
end

@testset "Paper IV sequent TPD" begin
    # Sequent TPD exercises axis maps, normalization, direct output, auxiliary variables, and reverse seeds.
    header = Dict{String,Any}(
        "NAXIS" => 2, "CTYPE1" => "LINEAR", "CTYPE2" => "LINEAR",
        "CRPIX1" => 10.0, "CRPIX2" => 20.0,
        "CRVAL1" => 5.0, "CRVAL2" => -7.0,
        "CDELT1" => -0.02, "CDELT2" => 0.035,
        "PC1_1" => 0.8, "PC1_2" => 0.3,
        "PC2_1" => -0.2, "PC2_2" => 1.1,
        "CQDIS1" => "TPD", "CQDIS2" => "TPD",
        "DQ1.NAXES" => 2, "DQ1.AXIS.1" => 1, "DQ1.AXIS.2" => 2,
        "DQ1.OFFSET.1" => 1.0, "DQ1.OFFSET.2" => -2.0,
        "DQ1.SCALE.1" => 0.5, "DQ1.SCALE.2" => 2.0, "DQ1.DOCORR" => 0,
        "DQ1.AUX.1.COEFF.0" => 0.1, "DQ1.AUX.1.COEFF.1" => 1.1, "DQ1.AUX.1.COEFF.2" => -0.2,
        "DQ1.AUX.2.COEFF.0" => -0.3, "DQ1.AUX.2.COEFF.1" => 0.4, "DQ1.AUX.2.COEFF.2" => 0.9,
        "DQ1.TPD.FWD.0" => 0.01, "DQ1.TPD.FWD.1" => 1.02, "DQ1.TPD.FWD.2" => -0.03,
        "DQ1.TPD.FWD.4" => 0.0007, "DQ1.TPD.FWD.5" => -0.0002, "DQ1.TPD.FWD.6" => 0.0005,
        "DQ1.TPD.REV.1" => 1.0,
        "DQ2.NAXES" => 2, "DQ2.AXIS.1" => 2, "DQ2.AXIS.2" => 1,
        "DQ2.TPD.FWD.4" => -0.002, "DQ2.TPD.FWD.5" => 0.0003,
    )
    wcs = WCS(header)
    pixel = [12.3, -3.4]
    expected_world = [4.820512437282, -7.963625782]
    @test pixel_to_intermediate(wcs, pixel) ≈ [-0.17948756271800004, -0.9636257820000002] atol=2e-15
    @test pixel_to_world(wcs, pixel) ≈ expected_world atol=2e-14
    @test world_to_pixel(wcs, expected_world) ≈ pixel atol=2e-12

    # minerr suppresses each complete TPD axis before constructing the pipeline.
    header["CQERR1"] = 0.01
    header["CQERR2"] = 0.01
    skipped = WCS(header; minerr = 0.02)
    @test pixel_to_world(skipped, pixel) ≈ [5.1036, -7.917] atol=2e-15
end

@testset "One-dimensional ninth-degree TPD" begin
    # One-dimensional TPD keeps pure attached-axis powers and ignores radial terms.
    header = Dict{String,Any}(
        "NAXIS" => 2, "CTYPE1" => "LINEAR", "CTYPE2" => "LINEAR",
        "CRPIX1" => 0.0, "CRPIX2" => 0.0,
        "CRVAL1" => 0.0, "CRVAL2" => 0.0,
        "CDELT1" => 1.0, "CDELT2" => 1.0,
        "CPDIS1" => "TPD", "CPDIS2" => "TPD",
        "DP1.NAXES" => 1, "DP1.AXIS.1" => 1,
        "DP1.TPD.FWD.49" => 1e-9, "DP1.TPD.FWD.59" => 10.0,
        "DP2.NAXES" => 1, "DP2.AXIS.1" => 2,
    )
    @test pixel_to_world(WCS(header), [2.0, 3.0]) ≈ [2.000000512, 3.0] atol=2e-15

    # Two-dimensional radial terms use both mapped coordinates.
    radial = copy(header)
    radial["DP1.NAXES"] = 2
    radial["DP1.AXIS.2"] = 2
    delete!(radial, "DP1.TPD.FWD.49")
    radial["DP1.TPD.FWD.59"] = 1e-9
    @test pixel_to_world(WCS(radial), [3.0, 4.0]) ≈ [3.001953125, 4.0] atol=2e-15
end

@testset "FITSFiles repeated TPD parameter cards" begin
    # Repeated DP cards must survive backend conversion as distinct flattened parameters.
    cards = FITSFiles.Card[
        FITSFiles.Card("NAXIS", 2),
        FITSFiles.Card("CTYPE1", "LINEAR"), FITSFiles.Card("CTYPE2", "LINEAR"),
        FITSFiles.Card("CRPIX1", 0.0), FITSFiles.Card("CRPIX2", 0.0),
        FITSFiles.Card("CRVAL1", 0.0), FITSFiles.Card("CRVAL2", 0.0),
        FITSFiles.Card("CDELT1", 1.0), FITSFiles.Card("CDELT2", 1.0),
        FITSFiles.Card("CPDIS1", "TPD"), FITSFiles.Card("CPDIS2", "TPD"),
        FITSFiles.Card("DP1", "NAXES: 2"), FITSFiles.Card("DP1", "AXIS.1: 1"),
        FITSFiles.Card("DP1", "AXIS.2: 2"), FITSFiles.Card("DP1", "TPD.FWD.4: 0.01"),
        FITSFiles.Card("DP2", "NAXES: 2"), FITSFiles.Card("DP2", "AXIS.1: 2"),
        FITSFiles.Card("DP2", "AXIS.2: 1"), FITSFiles.Card("DP2", "TPD.FWD.4: -0.02"),
    ]
    @test pixel_to_world(WCS(cards), [2.0, 3.0]) ≈ [2.04, 2.82] atol=2e-15
end

@testset "FITSIO repeated TPD parameter cards" begin
    # FITSIO's repeated string-card representation must retain every TPD field.
    extension = Base.get_extension(FITSWCS, :FITSWCSFITSIOExt)
    header = Dict{String,Any}()
    for value in ("NAXES: 2", "AXIS.1: 1", "AXIS.2: 2", "DOCORR: 0",
                  "OFFSET.1: 3", "SCALE.2: 4", "AUX.1.COEFF.0: 5",
                  "TPD.FWD.4: 0.01", "TPD.REV.4: -0.01")
        extension._set_fitsio_header_card!(header, "DP1", value)
    end
    @test header["DP1.NAXES"] == 2.0
    @test header["DP1.DOCORR"] == 0.0
    @test header["DP1.OFFSET.1"] == 3.0
    @test header["DP1.SCALE.2"] == 4.0
    @test header["DP1.AUX.1.COEFF.0"] == 5.0
    @test header["DP1.TPD.FWD.4"] == 0.01
    @test header["DP1.TPD.REV.4"] == -0.01
end
