"""
Linear pixel-to-intermediate and intermediate-to-pixel transforms.

Paper I (Greisen & Calabretta 2002), Section 2.

The intermediate world coordinate xᵢ is related to pixel coordinate pⱼ by

    xᵢ = Σⱼ CDᵢⱼ (pⱼ − CRPIXⱼ)

where `CDᵢⱼ = CDELTᵢ × PCᵢⱼ`.  All intermediate coordinates share the
units of their corresponding `CRVALᵢ` (typically degrees for celestial axes).

The inverse is

    pⱼ = CRPIXⱼ + Σᵢ (CD⁻¹)ⱼᵢ xᵢ
"""

"""
    pixel_to_intermediate(wcs, pixel) -> intermediate

Apply the linear pixel-to-intermediate transform.

`pixel` is a length-`naxis` vector of 1-based FITS pixel coordinates.
Returns a length-`naxis` vector of intermediate world coordinates in degrees
(for celestial axes) or in whatever units are implied by the CD matrix.
"""
function pixel_to_intermediate(wcs::WCSTransform, pixel::AbstractVector)
    # Apply pre-linear SIP distortion when present.
    focal = wcs.sip === nothing ? pixel : sip_pixel_to_focal(wcs.sip, pixel)

    # Convert focal/image-plane pixel coordinates to intermediate coordinates.
    return wcs.cd * (focal .- wcs.crpix)
end

"""
    intermediate_to_pixel(wcs, intermediate) -> pixel

Inverse linear transform: intermediate world coordinates → pixel coordinates.

Requires the CD matrix to be invertible.  Throws a `LinearAlgebra.SingularException`
if the matrix is singular.
"""
function intermediate_to_pixel(wcs::WCSTransform, intermediate::AbstractVector)
    # Undo the linear matrix to recover focal/image-plane pixel coordinates.
    focal = wcs.crpix .+ (wcs.cd \ collect(intermediate))

    # Convert focal/image-plane coordinates back to detector pixels if needed.
    return wcs.sip === nothing ? focal : sip_focal_to_pixel(wcs.sip, focal)
end
