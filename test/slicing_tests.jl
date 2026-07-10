@testset "WCS Slicing" begin

# ──────────────────────────────────────────────────────────────────────────────
@testset "axis_correlation_matrix" begin
    @testset "Diagonal linear WCS" begin
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "X", "CTYPE2" => "Y",
            "CDELT1" => 1.0, "CDELT2" => 1.0,
        )
        wcs = WCS(hdr)
        corr = axis_correlation_matrix(wcs)
        @test corr == Bool[1 0; 0 1]
    end

    @testset "Rotated linear WCS" begin
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "X", "CTYPE2" => "Y",
            "CD1_1"  => 1.0, "CD1_2" => 0.5,
            "CD2_1"  => 0.5, "CD2_2" => 1.0,
        )
        wcs = WCS(hdr)
        corr = axis_correlation_matrix(wcs)
        @test corr == Bool[1 1; 1 1]
    end

    @testset "Celestial pair shares dependencies" begin
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN",
            "CRPIX1" => 1.0, "CRPIX2" => 1.0,
            "CRVAL1" => 0.0, "CRVAL2" => 0.0,
            "CDELT1" => -1e-4, "CDELT2" => 1e-4,
        )
        wcs = WCS(hdr)
        corr = axis_correlation_matrix(wcs)
        # Celestial axes share dependencies: both depend on both pixel axes.
        @test corr == Bool[1 1; 1 1]
    end

    @testset "3D cube with independent spectral axis" begin
        hdr = Dict(
            "NAXIS"  => 3,
            "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CTYPE3" => "FREQ",
            "CRPIX1" => 512.0, "CRPIX2" => 512.0, "CRPIX3" => 1.0,
            "CRVAL1" => 83.8221, "CRVAL2" => -5.3911, "CRVAL3" => 1.42e9,
            "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4, "CDELT3" => 1.0e6,
        )
        wcs = WCS(hdr)
        corr = axis_correlation_matrix(wcs)
        # Celestial pair (axes 1-2) coupled; spectral (axis 3) independent.
        @test corr == Bool[1 1 0; 1 1 0; 0 0 1]
    end

    @testset "SIP distortion forces all-true" begin
        hdr = Dict{String,Any}(
            "NAXIS"  => 2,
            "CTYPE1" => "X", "CTYPE2" => "Y",
            "CRPIX1" => 10.0, "CRPIX2" => 20.0,
            "A_ORDER" => 2, "B_ORDER" => 2,
            "A_2_0" => 1e-3,
        )
        wcs = WCS(hdr)
        corr = axis_correlation_matrix(wcs)
        @test corr == Bool[1 1; 1 1]
    end

    @testset "SlicedWCSTransform delegates to parent submatrix" begin
        hdr = Dict(
            "NAXIS"  => 3,
            "CTYPE1" => "X", "CTYPE2" => "Y", "CTYPE3" => "Z",
            "CDELT1" => 1.0, "CDELT2" => 1.0, "CDELT3" => 1.0,
        )
        wcs = WCS(hdr)
        swcs = slice_wcs(wcs, 1:10, 1:10, 5)  # drop axis 3
        corr = axis_correlation_matrix(swcs)
        @test corr == Bool[1 0; 0 1]
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "pixel_n_dim / world_n_dim" begin
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "X", "CTYPE2" => "Y",
        "CDELT1" => 1.0, "CDELT2" => 1.0,
    )
    wcs = WCS(hdr)
    @test pixel_n_dim(wcs) == 2
    @test world_n_dim(wcs) == 2

    swcs = slice_wcs(wcs, 1:5, 1:10)
    @test pixel_n_dim(swcs) == 2
    @test world_n_dim(swcs) == 2

    # Coupled case: drop one pixel axis, both world axes survive.
    hdr_rot = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "X", "CTYPE2" => "Y",
        "CD1_1"  => 1.0, "CD1_2" => 0.5,
        "CD2_1"  => 0.5, "CD2_2" => 1.0,
    )
    wcs_rot = WCS(hdr_rot)
    swcs_rot = slice_wcs(wcs_rot, 10, 1:20)
    @test pixel_n_dim(swcs_rot) == 1
    @test world_n_dim(swcs_rot) == 2
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Linear 2D WCS slicing" begin
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "X", "CTYPE2" => "Y",
        "CRPIX1" => 10.0, "CRPIX2" => 20.0,
        "CRVAL1" => 100.0, "CRVAL2" => 200.0,
        "CDELT1" => 2.0, "CDELT2" => 3.0,
    )
    wcs = WCS(hdr)

    @testset "UnitRange slice → round-trip" begin
        swcs = slice_wcs(wcs, 3:7, 5:10)
        pix = [1.0, 1.0]  # first pixel of sliced image
        world = pixel_to_world(swcs, pix)
        @test world ≈ [100.0 + 2*(3 - 10), 200.0 + 3*(5 - 20)]
        @test world_to_pixel(swcs, world) ≈ pix
    end

    @testset "StepRange slice → round-trip" begin
        swcs = slice_wcs(wcs, 1:2:10, 3:3:15)
        pix = [2.0, 2.0]
        world = pixel_to_world(swcs, pix)
        # pixel 2 → original 1 + 2*(2-1) = 3; original 3 + 3*(2-1) = 6
        @test world ≈ [100.0 + 2*(3 - 10), 200.0 + 3*(6 - 20)]
        @test world_to_pixel(swcs, world) ≈ pix
    end

    @testset "Drop one axis via Integer" begin
        swcs = slice_wcs(wcs, 5, 1:10)
        @test pixel_n_dim(swcs) == 1
        @test world_n_dim(swcs) == 1
        # pixel axis 1 dropped at 5, pixel axis 2 = 1 (kept)
        world = pixel_to_world(swcs, [1.0])
        @test world ≈ [200.0 + 3*(1 - 20)]  # only Y world axis kept
        @test world_to_pixel(swcs, world) ≈ [1.0]
    end

    @testset "Multiple points round-trip" begin
        swcs = slice_wcs(wcs, 3:7, 5:10)
        for pix_sub in ([1.0, 1.0], [3.0, 4.0], [5.0, 6.0])
            world = pixel_to_world(swcs, pix_sub)
            @test world_to_pixel(swcs, world) ≈ pix_sub
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Linear 3D WCS slicing" begin
    hdr = Dict(
        "NAXIS"  => 3,
        "CTYPE1" => "X", "CTYPE2" => "Y", "CTYPE3" => "Z",
        "CRPIX1" => 1.0, "CRPIX2" => 1.0, "CRPIX3" => 1.0,
        "CRVAL1" => 10.0, "CRVAL2" => 20.0, "CRVAL3" => 30.0,
        "CDELT1" => 1.0, "CDELT2" => 2.0, "CDELT3" => 3.0,
    )
    wcs = WCS(hdr)

    @testset "Drop middle axis" begin
        swcs = slice_wcs(wcs, 1:10, 5, 1:10)
        @test pixel_n_dim(swcs) == 2
        # Axis 2 dropped; axes 1 and 3 kept.
        world = pixel_to_world(swcs, [1.0, 1.0])
        @test world ≈ [10.0, 30.0]  # at reference pixels
        @test world_to_pixel(swcs, world) ≈ [1.0, 1.0]
    end

    @testset "Step on all axes" begin
        swcs = slice_wcs(wcs, 1:2:10, 1:3:15, 1:1:20)
        pix = [1.0, 1.0, 1.0]
        world = pixel_to_world(swcs, pix)
        # all pixels = 1 → original 1
        @test world ≈ [10.0, 20.0, 30.0]
        @test world_to_pixel(swcs, world) ≈ pix
    end

    @testset "Drop first and step on remaining" begin
        swcs = slice_wcs(wcs, 3, 2:4:20, 1:2:15)
        @test pixel_n_dim(swcs) == 2
        @test world_n_dim(swcs) == 2
        # Axis 1 dropped at 3, kept axes 2 and 3.
        world = pixel_to_world(swcs, [1.0, 1.0])
        # axis 2: pixel 1 → original 2; axis 3: pixel 1 → original 1
        @test world ≈ [20.0 + 2*(2-1), 30.0 + 3*(1-1)]
        @test world_to_pixel(swcs, world) ≈ [1.0, 1.0]
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Celestial WCS slicing" begin
    @testset "3D cube: slice spectral axis" begin
        # Most common real-world use case.
        hdr = Dict(
            "NAXIS"  => 3,
            "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN", "CTYPE3" => "FREQ",
            "CRPIX1" => 512.0, "CRPIX2" => 512.0, "CRPIX3" => 1.0,
            "CRVAL1" => 83.8221, "CRVAL2" => -5.3911, "CRVAL3" => 1.42e9,
            "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4, "CDELT3" => 1.0e6,
        )
        wcs = WCS(hdr)

        # Slice the spectral axis at pixel 1, keep all spatial.
        swcs = slice_wcs(wcs, 1:1024, 1:1024, 1)
        @test pixel_n_dim(swcs) == 2
        @test world_n_dim(swcs) == 2

        # At reference pixel: should get CRVAL lon/lat.
        world = pixel_to_world(swcs, [512.0, 512.0])
        @test world ≈ [83.8221, -5.3911]
        @test world_to_pixel(swcs, world) ≈ [512.0, 512.0] atol=1e-10
    end

    @testset "2D celestial: spatial cutout" begin
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN",
            "CRPIX1" => 512.0, "CRPIX2" => 512.0,
            "CRVAL1" => 83.8221, "CRVAL2" => -5.3911,
            "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4,
        )
        wcs = WCS(hdr)
        swcs = slice_wcs(wcs, 400:600, 400:600)

        # At reference pixel in sliced image.
        pix_ref = [113.0, 113.0]  # original (512, 512) = (400+113-1, 400+113-1)
        world = pixel_to_world(swcs, pix_ref)
        @test world ≈ [83.8221, -5.3911] atol=1e-10
        @test world_to_pixel(swcs, world) ≈ pix_ref atol=1e-10
    end

    @testset "Coupled rotated 2D: drop one axis → both world axes survive" begin
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "X", "CTYPE2" => "Y",
            "CRPIX1" => 10.0, "CRPIX2" => 10.0,
            "CRVAL1" => 100.0, "CRVAL2" => 200.0,
            "CD1_1"  => 1.0, "CD1_2" => 0.5,
            "CD2_1"  => 0.5, "CD2_2" => 1.0,
        )
        wcs = WCS(hdr)
        swcs = slice_wcs(wcs, 10, 1:20)
        @test pixel_n_dim(swcs) == 1
        @test world_n_dim(swcs) == 2

        # At reference: pixel axis 2 = 10, pixel axis 1 = 10 (fixed).
        world = pixel_to_world(swcs, [10.0])
        @test world ≈ [100.0, 200.0]
        @test world_to_pixel(swcs, world) ≈ [10.0]
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "SIP distortion slicing" begin
    hdr = Dict{String,Any}(
        "NAXIS"  => 2,
        "CTYPE1" => "X", "CTYPE2" => "Y",
        "CRPIX1" => 10.0, "CRPIX2" => 20.0,
        "CRVAL1" => 100.0, "CRVAL2" => 200.0,
        "CDELT1" => 2.0, "CDELT2" => 3.0,
        "A_ORDER" => 2, "B_ORDER" => 2,
        "A_2_0" => 1e-3, "B_0_2" => -1e-3,
    )
    wcs = WCS(hdr)

    @testset "Range slice → round-trip" begin
        swcs = slice_wcs(wcs, 8:15, 18:25)
        # Pixel (3, 3) in sliced = original (10, 20) = SIP reference.
        pix = [3.0, 3.0]
        world = pixel_to_world(swcs, pix)
        @test world ≈ [100.0, 200.0] atol=1e-10
        @test world_to_pixel(swcs, world) ≈ pix atol=1e-8
    end

    @testset "Step slice → round-trip" begin
        swcs = slice_wcs(wcs, 8:2:14, 18:2:24)
        pix = [2.0, 2.0]  # original (10, 20)
        world = pixel_to_world(swcs, pix)
        @test world ≈ [100.0, 200.0] atol=1e-10
        @test world_to_pixel(swcs, world) ≈ pix atol=1e-8
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "TPV distortion slicing" begin
    hdr = Dict{String,Any}(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---TPV", "CTYPE2" => "DEC--TPV",
        "CRPIX1" => 512.0, "CRPIX2" => 512.0,
        "CRVAL1" => 83.8221, "CRVAL2" => -5.3911,
        "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4,
        "PV1_0" => 0.0, "PV1_1" => 1.0,
        "PV2_0" => 0.0, "PV2_1" => 0.0, "PV2_2" => 1.0,
    )
    wcs = WCS(hdr)

    @testset "Range slice → round-trip" begin
        swcs = slice_wcs(wcs, 400:600, 400:600)
        # Pixel (113, 113) in sliced = original (512, 512) = reference.
        pix = [113.0, 113.0]
        world = pixel_to_world(swcs, pix)
        @test world ≈ [83.8221, -5.3911] atol=1e-10
        @test world_to_pixel(swcs, world) ≈ pix atol=1e-8
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Float32 preservation" begin
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "X", "CTYPE2" => "Y",
        "CRPIX1" => 1.0, "CRPIX2" => 1.0,
        "CDELT1" => 2.0, "CDELT2" => 3.0,
    )
    wcs = WCS(hdr)
    swcs = slice_wcs(wcs, 3:7, 5:10)

    world32 = pixel_to_world(swcs, Float32[1.0, 1.0])
    @test eltype(world32) == Float32
    @test world32 ≈ Float32[4.0, 12.0]

    pix32 = world_to_pixel(swcs, Float32[4.0, 12.0])
    @test eltype(pix32) == Float32
    @test pix32 ≈ Float32[1.0, 1.0]
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Batch transforms" begin
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "X", "CTYPE2" => "Y",
        "CRPIX1" => 1.0, "CRPIX2" => 1.0,
        "CDELT1" => 1.0, "CDELT2" => 1.0,
    )
    wcs = WCS(hdr)
    swcs = slice_wcs(wcs, 3:7, 5:9)

    @testset "Matrix batch → round-trip" begin
        pix_batch = [1.0 3.0 5.0; 1.0 2.0 4.0]
        world_batch = pixel_to_world(swcs, pix_batch)
        @test size(world_batch) == (2, 3)
        @test world_to_pixel(swcs, world_batch) ≈ pix_batch
    end

    @testset "Vector-of-vectors batch → round-trip" begin
        pix_vec = [[1.0, 1.0], [3.0, 2.0], [5.0, 4.0]]
        world_vec = pixel_to_world(swcs, pix_vec)
        @test world_to_pixel(swcs, world_vec) ≈ hcat(pix_vec...)
    end

    @testset "Tuple and varargs convenience" begin
        @test pixel_to_world(swcs, (1.0, 1.0)) ≈ [2.0, 4.0]
        @test pixel_to_world(swcs, 1.0, 1.0) ≈ [2.0, 4.0]
        @test world_to_pixel(swcs, (2.0, 4.0)) ≈ [1.0, 1.0]
        @test world_to_pixel(swcs, 2.0, 4.0) ≈ [1.0, 1.0]
    end

    @testset "Batch with dropped axis" begin
        hdr3 = Dict(
            "NAXIS"  => 3,
            "CTYPE1" => "X", "CTYPE2" => "Y", "CTYPE3" => "Z",
            "CDELT1" => 1.0, "CDELT2" => 1.0, "CDELT3" => 1.0,
        )
        wcs3 = WCS(hdr3)
        swcs3 = slice_wcs(wcs3, 1:5, 3, 1:5)

        pix_batch = [1.0 3.0; 1.0 4.0]
        world_batch = pixel_to_world(swcs3, pix_batch)
        @test size(world_batch) == (2, 2)
        @test world_to_pixel(swcs3, world_batch) ≈ pix_batch
    end
end

# ──────────────────────────────────────────────────────────────────────────────
@testset "Edge cases and errors" begin
    hdr = Dict(
        "NAXIS"  => 2,
        "CRPIX1" => 0.0, "CRPIX2" => 0.0,
        "CRVAL1" => 0.0, "CRVAL2" => 0.0,
        "CTYPE1" => "X", "CTYPE2" => "Y",
        "CDELT1" => 1.0, "CDELT2" => 1.0,
    )
    wcs = WCS(hdr)

    @testset "All axes dropped → error" begin
        @test_throws ArgumentError slice_wcs(wcs, 1, 2)
    end

    @testset "Too many slice arguments → error" begin
        @test_throws ArgumentError slice_wcs(wcs, 1, 2, 3)
    end

    @testset "Missing trailing arguments default to Colon" begin
        # slice_wcs(wcs, 1:10) with a 2D WCS implicitly keeps axis 2.
        swcs = slice_wcs(wcs, 1:10)
        @test pixel_n_dim(swcs) == 2
        @test world_n_dim(swcs) == 2
        world = pixel_to_world(swcs, [1.0, 1.0])
        @test world ≈ [1.0, 1.0]
    end

    @testset "Explicit Colon keeps axis unchanged" begin
        swcs = slice_wcs(wcs, 1:5, :)
        @test pixel_n_dim(swcs) == 2
        # axis 1: pixel 1 → original 1, axis 2: pixel 3 → original 3 (identity).
        world = pixel_to_world(swcs, [1.0, 3.0])
        @test world ≈ [1.0, 3.0]
        @test world_to_pixel(swcs, world) ≈ [1.0, 3.0]
    end

    @testset "Invalid slice type → error" begin
        @test_throws ArgumentError slice_wcs(wcs, 1:10, "invalid")
    end

    @testset "Dimension mismatch for SlicedWCSTransform" begin
        swcs = slice_wcs(wcs, 3:7, 5:10)
        @test_throws DimensionMismatch pixel_to_world(swcs, [1.0])
        @test_throws DimensionMismatch pixel_to_world(swcs, [1.0, 2.0, 3.0])
        @test_throws DimensionMismatch world_to_pixel(swcs, [1.0])
    end

    @testset "preserve_units round-trip through sliced WCS" begin
        hdr_arcsec = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---TAN", "CTYPE2" => "DEC--TAN",
            "CRPIX1" => 1.0, "CRPIX2" => 1.0,
            "CRVAL1" => 7200.0, "CRVAL2" => 3600.0,
            "CDELT1" => 1.0, "CDELT2" => 1.0,
            "CUNIT1" => "arcsec", "CUNIT2" => "arcsec",
        )
        wcs_pu = WCS(hdr_arcsec; preserve_units = true)
        swcs = slice_wcs(wcs_pu, 1:10, 1:10)
        world = pixel_to_world(swcs, [1.0, 1.0])
        @test world ≈ [7200.0, 3600.0]
        @test world_to_pixel(swcs, world) ≈ [1.0, 1.0]
    end
end

end # @testset "WCS Slicing"
