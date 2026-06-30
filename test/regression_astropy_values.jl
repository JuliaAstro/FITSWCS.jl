"""
Reference-comparison regression tests against Astropy.

Expected values were generated with Astropy 6.1.7 using
`WCS(header, relax=true).all_pix2world(pixels, 1)`, where `origin=1` matches
the FITS 1-based pixel convention used by FITSWCS.jl.
"""

using Test
using FITSWCS

@testset "ARC projection (Astropy comparison)" begin
    # Validate ARC against stored Astropy values because round-trip tests alone
    # can miss a symmetric projection-scale error.
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---ARC",  "CTYPE2" => "DEC--ARC",
        "CRPIX1" => 128.0,       "CRPIX2" => 96.0,
        "CRVAL1" => 120.0,       "CRVAL2" => 35.0,
        "CDELT1" => -0.05,       "CDELT2" => 0.05,
    )
    wcs = from_header(hdr)

    for (pix, world_ref) in [
        ([128.0,  96.0], [120.0, 35.0]),
        ([140.0,  96.0], [119.267548373357940, 34.997800282811028]),
        ([128.0, 110.0], [120.0, 35.699999999999989]),
        ([100.0,  80.0], [121.692487124705124, 34.188218938603562]),
        ([160.0, 120.0], [118.017663701139099, 36.183966405908677]),
    ]
        world = pixel_to_world(wcs, pix)
        @test world ≈ world_ref atol=1e-10
        @test world_to_pixel(wcs, world_ref) ≈ pix atol=1e-8
    end
end

@testset "ZEA projection (Astropy comparison)" begin
    # Validate ZEA against stored Astropy values to cover absolute projection
    # coordinates, not just inverse consistency.
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---ZEA",  "CTYPE2" => "DEC--ZEA",
        "CRPIX1" => 128.0,       "CRPIX2" => 96.0,
        "CRVAL1" => 120.0,       "CRVAL2" => 35.0,
        "CDELT1" => -0.05,       "CDELT2" => 0.05,
    )
    wcs = from_header(hdr)

    for (pix, world_ref) in [
        ([128.0,  96.0], [120.0, 35.0]),
        ([140.0,  96.0], [119.267545026673716, 34.997800262709013]),
        ([128.0, 110.0], [120.0, 35.700004353563713]),
        ([100.0,  80.0], [121.692542438183011, 34.188191764846081]),
        ([160.0, 120.0], [118.017561544720678, 36.184025689320833]),
    ]
        world = pixel_to_world(wcs, pix)
        @test world ≈ world_ref atol=1e-10
        @test world_to_pixel(wcs, world_ref) ≈ pix atol=1e-8
    end
end

@testset "CEA projection (Astropy comparison)" begin
    # Validate CEA with a non-default lambda parameter so PV parsing and
    # cylindrical equal-area scaling are both externally checked.
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---CEA",  "CTYPE2" => "DEC--CEA",
        "CRPIX1" => 181.0,       "CRPIX2" => 91.0,
        "CRVAL1" => 40.0,        "CRVAL2" => 0.0,
        "CDELT1" => -1.0,        "CDELT2" => 1.0,
        "PV2_1"  => 0.75,
    )
    wcs = from_header(hdr)

    for (pix, world_ref) in [
        ([181.0,  91.0], [40.0, 0.0]),
        ([170.0, 100.0], [51.0, 6.765712354999467]),
        ([200.0,  80.0], [21.0, -8.278777211051736]),
        ([150.0,  60.0], [71.0, -23.940582595236908]),
        ([210.0, 120.0], [11.0, 22.309472280415399]),
    ]
        world = pixel_to_world(wcs, pix)
        @test world ≈ world_ref atol=1e-10
        @test world_to_pixel(wcs, world_ref) ≈ pix atol=1e-8
    end
end

@testset "Additional WCSLIB projection codes (Astropy comparison)" begin
    # Validate the first projection-expansion slice against stored Astropy
    # values so the new formulas are checked outside round-trip symmetry.
    refs = Dict(
        "AZP" => [
            ([101.0, 81.0], [30.0, -10.0]),
            ([110.0, 90.0], [25.50225799973873, -5.492356509625904]),
            ([80.0, 60.0], [40.88513357349926, -20.047444732794045]),
            ([120.0, 70.0], [20.28164106558495, -15.271451202182682]),
        ],
        "SZP" => [
            ([101.0, 81.0], [30.0, -10.0]),
            ([110.0, 90.0], [25.502257999738685, -5.492356509625858]),
            ([80.0, 60.0], [40.88513357349926, -20.047444732794045]),
            ([120.0, 70.0], [20.28164106558495, -15.271451202182682]),
        ],
        "CYP" => [
            ([101.0, 81.0], [30.0, -10.0]),
            ([110.0, 90.0], [25.493309718808764, -5.471594360589324]),
            ([80.0, 60.0], [41.01485018705852, -20.29598959498479]),
            ([120.0, 70.0], [20.19064257183051, -15.35488544402523]),
        ],
        "MER" => [
            ([101.0, 81.0], [30.0, -10.0]),
            ([110.0, 90.0], [25.493277982512954, -5.473902349204117]),
            ([80.0, 60.0], [41.01380556079776, -20.266961601401988]),
            ([120.0, 70.0], [20.190772321298176, -15.350677488836727]),
        ],
        "SFL" => [
            ([101.0, 81.0], [30.0, -10.0]),
            ([110.0, 90.0], [25.479407365355677, -5.469093373275083]),
            ([80.0, 60.0], [41.203287790375555, -20.319228679179997]),
            ([120.0, 70.0], [20.14517683081226, -15.357798258796976]),
        ],
        "PAR" => [
            ([101.0, 81.0], [30.0, -10.0]),
            ([110.0, 90.0], [25.479257388932428, -5.671479180738607]),
            ([80.0, 60.0], [41.150731126139675, -19.853100979351957]),
            ([120.0, 70.0], [20.1614031505245, -15.111086198565213]),
        ],
        "MOL" => [
            ([101.0, 81.0], [30.0, -10.0]),
            ([110.0, 90.0], [24.979836656715705, -5.9092224081239495]),
            ([80.0, 60.0], [42.29667953176584, -19.25145201212617]),
            ([120.0, 70.0], [19.098842743047673, -14.779881183544061]),
        ],
        "PCO" => [
            ([101.0, 81.0], [30.0, -10.0]),
            ([110.0, 90.0], [25.47927349034995, -5.482957694859953]),
            ([80.0, 60.0], [41.19267068916983, -20.144144812286417]),
            ([120.0, 70.0], [20.148314355587686, -15.283052955995876]),
        ],
    )

    for code in ("AZP", "SZP", "CYP", "MER", "SFL", "PAR", "MOL", "PCO")
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---$code", "CTYPE2" => "DEC--$code",
            "CRPIX1" => 101.0,        "CRPIX2" => 81.0,
            "CRVAL1" => 30.0,         "CRVAL2" => -10.0,
            "CDELT1" => -0.5,         "CDELT2" => 0.5,
        )
        wcs = from_header(hdr)

        for (pix, world_ref) in refs[code]
            world = pixel_to_world(wcs, pix)
            @test world ≈ world_ref atol=1e-10
            @test world_to_pixel(wcs, world_ref) ≈ pix atol=1e-7
        end
    end
end

@testset "Celestial CUNIT arcsec (Astropy comparison)" begin
    # Astropy normalizes celestial arcsecond units to degree-valued world
    # coordinates; FITSWCS should expose the same public convention.
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---TAN",  "CTYPE2" => "DEC--TAN",
        "CRPIX1" => 100.0,       "CRPIX2" => 100.0,
        "CRVAL1" => 36000.0,     "CRVAL2" => 72000.0,
        "CDELT1" => -1.0,        "CDELT2" => 1.0,
        "CUNIT1" => "arcsec",    "CUNIT2" => "arcsec",
    )
    wcs = from_header(hdr)

    for (pix, world_ref) in [
        ([100.0, 100.0], [10.0, 20.0]),
        ([110.0,  95.0], [9.997043976715300, 19.998611086605219]),
        ([ 80.0, 120.0], [10.005912307369485, 20.005555457476497]),
    ]
        world = pixel_to_world(wcs, pix)
        @test world ≈ world_ref atol=1e-10
        @test world_to_pixel(wcs, world_ref) ≈ pix atol=1e-6
    end
end

@testset "Celestial CUNIT rad (Astropy comparison)" begin
    # Astropy normalizes celestial radian units to degree-valued world
    # coordinates; this checks both CRVAL and CDELT conversion.
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---TAN",  "CTYPE2" => "DEC--TAN",
        "CRPIX1" => 100.0,       "CRPIX2" => 100.0,
        "CRVAL1" => 0.5,         "CRVAL2" => 0.25,
        "CDELT1" => -1.0e-4,     "CDELT2" => 1.0e-4,
        "CUNIT1" => "rad",       "CUNIT2" => "rad",
    )
    wcs = from_header(hdr)

    for (pix, world_ref) in [
        ([100.0, 100.0], [28.647889756541161, 14.323944878270581]),
        ([110.0,  95.0], [28.588763210302314, 14.295289691157484]),
        ([ 80.0, 120.0], [28.766218248792566, 14.438506780437534]),
    ]
        world = pixel_to_world(wcs, pix)
        @test world ≈ world_ref atol=1e-10
        @test world_to_pixel(wcs, world_ref) ≈ pix atol=1e-6
    end
end

@testset "Split celestial and spectral axes (Astropy comparison)" begin
    # Validate mixed-axis ordering when a linear spectral axis separates the
    # celestial longitude and latitude axes.
    hdr = Dict(
        "NAXIS"  => 3,
        "CTYPE1" => "RA---TAN",  "CTYPE2" => "FREQ",  "CTYPE3" => "DEC--TAN",
        "CRPIX1" => 30.0,        "CRPIX2" => 40.0,    "CRPIX3" => 45.0,
        "CRVAL1" => 10.0,        "CRVAL2" => 1.42e9,  "CRVAL3" => 25.0,
        "CDELT1" => -0.01,       "CDELT2" => 1.0e6,   "CDELT3" => 0.01,
    )
    wcs = from_header(hdr)

    for (pix, world_ref) in [
        ([30.0, 40.0, 45.0], [10.0, 1.42e9, 25.0]),
        ([29.0, 43.0, 44.0], [10.011032881130753, 1.423e9, 24.989999593356771]),
        ([31.5, 37.0, 47.0], [9.983446637250788, 1.417e9, 25.019999082760339]),
        ([20.0, 50.0, 35.0], [10.110247929297399, 1.43e9, 24.899959594018089]),
    ]
        world = pixel_to_world(wcs, pix)
        @test world ≈ world_ref atol=1e-7
        @test world_to_pixel(wcs, world_ref) ≈ pix atol=1e-6
    end
end

@testset "Mixed celestial, time, and Stokes axes (Astropy comparison)" begin
    # Validate the basic linear TIME/STOKES subset without claiming physical
    # interpretation beyond Paper I linear axes.
    hdr = Dict(
        "NAXIS"  => 4,
        "CTYPE1" => "RA---TAN",  "CTYPE2" => "DEC--TAN",
        "CTYPE3" => "TIME",      "CTYPE4" => "STOKES",
        "CRPIX1" => 10.0,        "CRPIX2" => 20.0,
        "CRPIX3" => 1.0,         "CRPIX4" => 1.0,
        "CRVAL1" => 30.0,        "CRVAL2" => -5.0,
        "CRVAL3" => 59000.0,     "CRVAL4" => 1.0,
        "CDELT1" => -0.001,      "CDELT2" => 0.001,
        "CDELT3" => 0.5,         "CDELT4" => 1.0,
    )
    wcs = from_header(hdr)

    for (pix, world_ref) in [
        ([10.0, 20.0,  1.0, 1.0], [30.0, -5.0, 59000.0, 1.0]),
        ([10.0, 20.0,  3.0, 4.0], [30.0, -5.0, 59001.0, 4.0]),
        ([12.0, 18.0,  5.0, 2.0], [29.997992354194512, -5.001999996944045, 59002.0, 2.0]),
        ([ 7.5, 22.5, -1.0, 3.0], [30.002509540012301, -4.997499995232223, 58999.0, 3.0]),
    ]
        world = pixel_to_world(wcs, pix)
        @test world ≈ world_ref atol=1e-7
        @test world_to_pixel(wcs, world_ref) ≈ pix atol=1e-6
    end
end

@testset "ZPN projection (Astropy comparison)" begin
    hdr = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---ZPN",  "CTYPE2" => "DEC--ZPN",
        "CRPIX1" => 128.0,       "CRPIX2" => 96.0,
        "CRVAL1" => 120.0,       "CRVAL2" => 35.0,
        "CDELT1" => -0.05,       "CDELT2" => 0.05,
        "PV2_0"  => 0.0,         "PV2_1" => 1.0,
    )
    wcs = from_header(hdr)
    for (pix, world_ref) in [
        ([128.0,  96.0], [120.0, 35.0]),
        ([140.0,  96.0], [119.26754837335794, 34.99780028281103]),
        ([128.0, 110.0], [120.0, 35.69999999999999]),
        ([100.0,  80.0], [121.69248712470512, 34.18821893860356]),
        ([160.0, 120.0], [118.0176637011391, 36.18396640590868]),
    ]
        world = pixel_to_world(wcs, pix)
        @test world ≈ world_ref atol=1e-10
        @test world_to_pixel(wcs, world_ref) ≈ pix atol=1e-8
    end

    # Quadratic ZPN
    hdr2 = Dict(
        "NAXIS"  => 2,
        "CTYPE1" => "RA---ZPN",  "CTYPE2" => "DEC--ZPN",
        "CRPIX1" => 128.0,       "CRPIX2" => 96.0,
        "CRVAL1" => 120.0,       "CRVAL2" => 35.0,
        "CDELT1" => -0.05,       "CDELT2" => 0.05,
        "PV2_0"  => 0.0,         "PV2_1" => 1.0,  "PV2_2" => 0.1,
    )
    wcs2 = from_header(hdr2)
    for (pix, world_ref) in [
        ([128.0,  96.0], [120.0, 35.0]),
        ([140.0,  96.0], [119.26831376526451, 34.99780487775835]),
        ([128.0, 110.0], [120.0, 35.69914687198124]),
        ([100.0,  80.0], [121.68779664466189, 34.19052306871118]),
        ([160.0, 120.0], [118.02463752599317, 36.17991893853489]),
    ]
        world = pixel_to_world(wcs2, pix)
        @test world ≈ world_ref atol=1e-10
        @test world_to_pixel(wcs2, world_ref) ≈ pix atol=1e-8
    end
end

@testset "AIR projection (Astropy comparison)" begin
    for (theta_b, refs) in [
        (90.0, [
            ([128.0,  96.0], [120.0, 35.0]),
            ([140.0,  96.0], [119.26755004668192, 34.99780029286188]),
            ([128.0, 110.0], [120.0, 35.699997823255245]),
            ([100.0,  80.0], [121.69245946850039, 34.18823252520399]),
            ([160.0, 120.0], [118.01771477773288, 36.183936765073874]),
        ]),
        (45.0, [
            ([128.0,  96.0], [120.0, 35.0]),
            ([140.0,  96.0], [119.23818894479591, 34.99762039997232]),
            ([128.0, 110.0], [120.0, 35.72805914003979]),
            ([100.0,  80.0], [121.75962048333146, 34.15520728913019]),
            ([160.0, 120.0], [117.93701469091603, 36.23071372079401]),
        ]),
        (30.0, [
            ([128.0,  96.0], [120.0, 35.0]),
            ([140.0,  96.0], [119.21370912218947, 34.99746500968105]),
            ([128.0, 110.0], [120.0, 35.75145537251257]),
            ([100.0,  80.0], [121.81557645996834, 34.1276442398319]),
            ([160.0, 120.0], [117.86965611191242, 36.269673163902375]),
        ]),
    ]
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---AIR",  "CTYPE2" => "DEC--AIR",
            "CRPIX1" => 128.0,       "CRPIX2" => 96.0,
            "CRVAL1" => 120.0,       "CRVAL2" => 35.0,
            "CDELT1" => -0.05,       "CDELT2" => 0.05,
        )
        if theta_b != 90.0
            hdr["PV2_1"] = theta_b
        end
        wcs = from_header(hdr)
        for (pix, world_ref) in refs
            world = pixel_to_world(wcs, pix)
            @test world ≈ world_ref atol=1e-10
            @test world_to_pixel(wcs, world_ref) ≈ pix atol=1e-8
        end
    end
end

@testset "Conic projections (Astropy comparison)" begin
    for (code, sigma, delta, refs) in [
        ("COP", 50.0, 20.0, [
            ([128.0,  96.0], [120.0, 35.0]),
            ([140.0,  96.0], [119.22054311078567, 34.997508840329246]),
            ([128.0, 110.0], [120.0, 35.744882471976084]),
            ([100.0,  80.0], [121.79984541747162, 34.13530271405697]),
            ([160.0, 120.0], [117.88876360949762, 36.25882207068821]),
        ]),
        ("COD", 50.0, 20.0, [
            ([128.0,  96.0], [120.0, 35.0]),
            ([140.0,  96.0], [119.22054352690223, 34.99784494509101]),
            ([128.0, 110.0], [120.0, 35.7]),
            ([100.0,  80.0], [121.79970244406476, 34.1883699953639]),
            ([160.0, 120.0], [117.88848534551155, 36.184469862072056]),
        ]),
        ("COE", 50.0, 20.0, [
            ([128.0,  96.0], [120.0, 34.999999999999986]),
            ([140.0,  96.0], [119.22281223848552, 34.99822040224724]),
            ([128.0, 110.0], [120.0, 35.65945458076792]),
            ([100.0,  80.0], [121.79512970531837, 34.23604832834036]),
            ([160.0, 120.0], [117.89575988235443, 36.117457204567046]),
        ]),
        ("COO", 50.0, 20.0, [
            ([128.0,  96.0], [120.0, 35.0]),
            ([140.0,  96.0], [119.21970166898389, 34.99740790590618]),
            ([128.0, 110.0], [120.0, 35.745839642927315]),
            ([100.0,  80.0], [121.80104633350241, 34.13398818748572]),
            ([160.0, 120.0], [117.88518463787781, 36.25992655890607]),
        ]),
    ]
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---$code",  "CTYPE2" => "DEC--$code",
            "CRPIX1" => 128.0,         "CRPIX2" => 96.0,
            "CRVAL1" => 120.0,         "CRVAL2" => 35.0,
            "CDELT1" => -0.05,         "CDELT2" => 0.05,
            "PV2_1"  => sigma,         "PV2_2"  => delta,
        )
        wcs = from_header(hdr)
        for (pix, world_ref) in refs
            world = pixel_to_world(wcs, pix)
            @test world ≈ world_ref atol=1e-10
            @test world_to_pixel(wcs, world_ref) ≈ pix atol=1e-8
        end
    end
end

@testset "BON projection (Astropy comparison)" begin
    for (theta1, refs) in [
        (45.0, [
            ([128.0,  96.0], [120.0, 35.0]),
            ([140.0,  96.0], [119.26755992296864, 34.99604071663145]),
            ([128.0, 110.0], [120.0, 35.7]),
            ([100.0,  80.0], [121.69240176434272, 34.17863386316033]),
            ([160.0, 120.0], [118.01775259718073, 36.17146088488652]),
        ]),
        (30.0, [
            ([128.0,  96.0], [120.0, 35.0]),
            ([140.0,  96.0], [119.26755820821381, 34.996407545680896]),
            ([128.0, 110.0], [120.0, 35.69999999999999]),
            ([100.0,  80.0], [121.69242195840391, 34.1806031776717]),
            ([160.0, 120.0], [118.01771634524357, 36.174124748366154]),
        ]),
        (0.0, [
            ([128.0,  96.0], [120.0, 35.0]),
            ([140.0,  96.0], [119.26754837335795, 34.99780028281102]),
            ([128.0, 110.0], [120.0, 35.7]),
            ([100.0,  80.0], [121.6925404912451, 34.18813856862629]),
            ([160.0, 120.0], [118.01751492655592, 36.18412003847978]),
        ]),
    ]
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---BON",  "CTYPE2" => "DEC--BON",
            "CRPIX1" => 128.0,       "CRPIX2" => 96.0,
            "CRVAL1" => 120.0,       "CRVAL2" => 35.0,
            "CDELT1" => -0.05,       "CDELT2" => 0.05,
            "PV2_1"  => theta1,
        )
        wcs = from_header(hdr)
        for (pix, world_ref) in refs
            world = pixel_to_world(wcs, pix)
            @test world ≈ world_ref atol=1e-10
            @test world_to_pixel(wcs, world_ref) ≈ pix atol=1e-8
        end
    end
end

@testset "Quadcube projections (Astropy comparison)" begin
    for (code, refs) in [
        ("TSC", [
            ([181.0,  91.0], [0.0, 0.0]),
            ([200.0,  91.0], [337.1094483437517, 0.0]),
            ([150.0,  91.0], [34.562524648881826, 0.0]),
            ([181.0, 110.0], [0.0, 22.890551656248327]),
            ([181.0,  70.0], [0.0, -25.01689347810002]),
            ([181.0, 140.0], [0.0, 47.66300076606715]),
            ([181.0,  40.0], [0.0, -49.085616779974885]),
        ]),
        ("CSC", [
            ([181.0,  91.0], [0.0, 0.0]),
            ([200.0,  91.0], [342.1944963868453, 0.0]),
            ([150.0,  91.0], [29.70172368252648, 0.0]),
            ([181.0, 110.0], [0.0, 17.805503640889842]),
            ([181.0,  70.0], [0.0, -19.729070839390136]),
            ([181.0, 140.0], [0.0, 49.49898534337126]),
            ([181.0,  40.0], [0.0, -51.72912697934343]),
        ]),
        ("QSC", [
            ([181.0,  91.0], [0.0, 0.0]),
            ([200.0,  91.0], [341.4030622394675, 0.0]),
            ([150.0,  91.0], [30.57069424698958, 0.0]),
            ([181.0, 110.0], [0.0, 18.59693776053253]),
            ([181.0,  70.0], [0.0, -20.57477188521203]),
            ([181.0, 140.0], [0.0, 49.18837536120323]),
            ([181.0,  40.0], [0.0, -51.261001060330926]),
        ]),
    ]
        hdr = Dict(
            "NAXIS"  => 2,
            "CTYPE1" => "RA---$code",  "CTYPE2" => "DEC--$code",
            "CRPIX1" => 181.0,         "CRPIX2" => 91.0,
            "CRVAL1" => 0.0,           "CRVAL2" => 0.0,
            "CDELT1" => -1.0,          "CDELT2" => 1.0,
        )
        wcs = from_header(hdr)
        # CSC uses float32 coefficients in WCSLIB; allow 2e-6 tolerance.
        tol = code == "CSC" ? 2e-6 : 1e-10
        for (pix, world_ref) in refs
            world = pixel_to_world(wcs, pix)
            @test world ≈ world_ref atol=tol
        end
    end
end
