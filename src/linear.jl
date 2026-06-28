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
function pixel_to_intermediate(wcs::WCSTransform{N}, pixel::StaticVector{N}) where {N}
    # Apply pre-linear SIP distortion when present.
    focal = wcs.sip === nothing ? pixel : sip_pixel_to_focal(wcs.sip, pixel)

    # Build a static offset vector before applying the static CD matrix.
    T = _coordinate_float_type(pixel)
    delta = SVector{N,T}(ntuple(i -> T(focal[i]) - T(wcs.crpix[i]), N))
    return wcs.cd * delta
end

function pixel_to_intermediate(wcs::WCSTransform{N}, pixel::AbstractVector) where {N}
    length(pixel) == N ||
        throw(DimensionMismatch("pixel has length $(length(pixel)), expected $N"))

    # Apply pre-linear SIP distortion when present.
    focal = wcs.sip === nothing ? pixel : sip_pixel_to_focal(wcs.sip, pixel)

    # Use a static temporary for the matrix product, then return ordinary storage.
    T = _coordinate_float_type(pixel)
    delta = SVector{N,T}(ntuple(i -> T(focal[i]) - T(wcs.crpix[i]), N))
    return Vector{T}(wcs.cd * delta)
end

"""
    intermediate_to_pixel(wcs, intermediate) -> pixel

Inverse linear transform: intermediate world coordinates → pixel coordinates.

Requires the CD matrix to be invertible.  Throws a `LinearAlgebra.SingularException`
if the matrix is singular.
"""
function intermediate_to_pixel(wcs::WCSTransform{N}, intermediate::StaticVector{N}) where {N}
    # Undo the linear matrix to recover focal/image-plane pixel coordinates.
    focal = wcs.crpix .+ (wcs.cd \ intermediate)

    # Convert focal/image-plane coordinates back to detector pixels if needed.
    pixel = wcs.sip === nothing ? focal : sip_focal_to_pixel(wcs.sip, focal)
    return SVector{N,_coordinate_float_type(intermediate)}(pixel)
end

function intermediate_to_pixel(wcs::WCSTransform{N}, intermediate::AbstractVector) where {N}
    length(intermediate) == N ||
        throw(DimensionMismatch("intermediate has length $(length(intermediate)), expected $N"))

    # Use a static temporary for the matrix solve, then return ordinary storage.
    T = _coordinate_float_type(intermediate)
    x = SVector{N,T}(ntuple(i -> T(intermediate[i]), N))
    focal = wcs.crpix .+ (wcs.cd \ x)

    # Convert focal/image-plane coordinates back to detector pixels if needed.
    pixel = wcs.sip === nothing ? focal : sip_focal_to_pixel(wcs.sip, focal)
    return Vector{T}(pixel)
end
