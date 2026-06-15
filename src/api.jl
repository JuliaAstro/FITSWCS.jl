"""
Small public API compatibility helpers.

These names mirror common WCS.jl workflows, but keep FITSWCS.jl's FITS
1-based pixel coordinate convention.
"""

"""
    WCS(header; alt=' ') -> WCSTransform

Construct a `WCSTransform` from a FITS-like header object.
"""
function WCS(header; alt::Char=' ')
    # Keep construction behavior centralized in the header parser.
    return from_header(header; alt=alt)
end

"""
    pix_to_world(wcs, pixels)
    pix_to_world(wcs, x, y, ...)

Compatibility alias for `pixel_to_world`.
"""
function pix_to_world(wcs::WCSTransform, pixels::Union{AbstractVector, AbstractMatrix, Tuple})
    # Delegate to the canonical transform implementation.
    return pixel_to_world(wcs, pixels)
end

function pix_to_world(wcs::WCSTransform, coords::Real...)
    # Collect scalar coordinate arguments into the vector form expected internally.
    return pixel_to_world(wcs, collect(Float64, coords))
end

"""
    world_to_pix(wcs, worlds)
    world_to_pix(wcs, x, y, ...)

Compatibility alias for `world_to_pixel`.
"""
function world_to_pix(wcs::WCSTransform, worlds::Union{AbstractVector, AbstractMatrix, Tuple})
    # Delegate to the canonical transform implementation.
    return world_to_pixel(wcs, worlds)
end

function world_to_pix(wcs::WCSTransform, coords::Real...)
    # Collect scalar coordinate arguments into the vector form expected internally.
    return world_to_pixel(wcs, collect(Float64, coords))
end
