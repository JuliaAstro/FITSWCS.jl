"""
    FITSWCS

Pure-Julia implementation of the FITS World Coordinate System standard.

## References

- Paper I:  Greisen & Calabretta (2002), A&A, 395, 1061.
            Linear WCS, core keywords.
- Paper II: Calabretta & Greisen (2002), A&A, 395, 1077.
            Celestial projections and spherical rotation.

## Quick start

```julia
using FITSWCS

# Build a WCS from a header dictionary (keys are FITS keyword names).
hdr = Dict(
    "NAXIS"  => 2,
    "CTYPE1" => "RA---TAN",
    "CTYPE2" => "DEC--TAN",
    "CRPIX1" => 512.0, "CRPIX2" => 512.0,
    "CRVAL1" => 83.8221, "CRVAL2" => -5.3911,
    "CDELT1" => -2.7778e-4, "CDELT2" => 2.7778e-4,
)
wcs = WCS(hdr)

# Convert pixel -> world (RA/Dec in degrees)
world = pixel_to_world(wcs, [512.0, 512.0])   # -> [83.8221, -5.3911]

# Convert world -> pixel
pix   = world_to_pixel(wcs, world)             # -> [512.0, 512.0]

# Return world coordinates in original header units (e.g., arcsec, Angstrom)
wcs_preserved = WCS(hdr; preserve_units = true)
world_as = pixel_to_world(wcs_preserved, [512.0, 512.0])

# Batch transform: each column of pixels is one coordinate.
pix_batch = [1.0 512.0 1024.0;
             1.0 512.0 1024.0]
world_batch = pixel_to_world(wcs, pix_batch)
```
"""
module FITSWCS

using LinearAlgebra: I, \
using StaticArrays: SMatrix, SVector, StaticVector, MVector, MMatrix

include("utilities.jl")
include("projections.jl")
include("lookup_tables.jl")
include("tabular.jl")
include("spectral.jl")
include("auxiliary_data.jl")
include("distortion.jl")
include("celestial.jl")
include("parsing.jl")
include("linear.jl")
include("transforms.jl")
include("api.jl")

export
    # Types
    WCSTransform,
    SIPDistortion,
    AbstractProjection,
    AZP, SZP, TAN, TPV, SIN, STG, ARC, ZEA, CAR, CEA, CYP, MER, SFL, PAR, MOL, PCO, AIT, UnknownProjection,

    # Parsing
    WCS,

    # Transforms
    pixel_to_world,
    world_to_pixel,
    pix_to_world,
    world_to_pix,
    pix_to_world!,
    world_to_pix!,

    # Projection primitives (useful for extension and testing)
    intermediate_to_native,
    native_to_intermediate,

    # Spherical rotation primitives
    native_to_celestial,
    celestial_to_native

end # module FITSWCS
