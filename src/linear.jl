"""
Linear pixel-to-intermediate and intermediate-to-pixel transforms.

Paper I (Greisen & Calabretta 2002), Section 2.

Without a sequent distortion, the intermediate world coordinate xᵢ is related
to pixel coordinate pⱼ by

    xᵢ = Σⱼ CDᵢⱼ (pⱼ − CRPIXⱼ)

where `CDᵢⱼ = CDELTᵢ × PCᵢⱼ`.  All intermediate coordinates share the
units of their corresponding `CRVALᵢ` (typically degrees for celestial axes).

The inverse is

    pⱼ = CRPIXⱼ + Σᵢ (CD⁻¹)ⱼᵢ xᵢ

With TPV or a Paper IV sequent distortion D, WCSLIB instead evaluates

    q = PC (p − CRPIX),    x = CDELT ⊙ D(q)

so the PC/CDELT decomposition must be preserved.
"""

"""
    pixel_to_intermediate(wcs, pixel) -> intermediate

Apply the linear pixel-to-intermediate transform.

`pixel` is a length-`naxis` vector of 1-based FITS pixel coordinates.
Returns a length-`naxis` vector of intermediate world coordinates in degrees
(for celestial axes) or in whatever units are implied by the CD matrix.
"""
function pixel_to_intermediate(wcs::WCSTransform{N}, pixel::AbstractVector) where {N}
    length(pixel) == N ||
        throw(DimensionMismatch("pixel has length $(length(pixel)), expected $N"))

    T = _coordinate_float_type(pixel)
    focal = pixel_to_focal(wcs.pipeline, pixel, Val(N))
    delta = SVector{N,T}(ntuple(i -> T(focal[i]) - T(wcs.crpix[i]), N))

    # The ordinary affine path can use the combined CD matrix directly.
    if !(wcs.projection isa TPV) && !has_sequent_distortion(wcs.pipeline)
        cd_T = SMatrix{N, N, T}(wcs.cd)
        return cd_T * delta
    end

    # Sequent distortions operate after PC and before the per-axis CDELT scale.
    pc_T = SMatrix{N, N, T}(wcs.pc)
    coord = pc_T * delta
    if wcs.projection isa TPV
        lon = wcs.lon_axis
        lat = wcs.lat_axis
        x = _evaluate_tpv_polynomial(wcs.projection.xcoeff, coord[lon], coord[lat])
        y = _evaluate_tpv_polynomial(wcs.projection.ycoeff, coord[lat], coord[lon])
        coord = SVector{N, T}(ntuple(i -> i == lon ? T(x) : i == lat ? T(y) : coord[i], N))
    end
    coord = apply_sequent_distortion(wcs.pipeline, coord)

    # Apply the final scale in the input coordinate's floating-point type.
    return coord .* SVector{N, T}(wcs.cdelt)
end


"""
    intermediate_to_pixel(wcs, intermediate) -> pixel

Inverse linear transform: intermediate world coordinates → pixel coordinates.

Requires the CD matrix to be invertible.  Throws a `LinearAlgebra.SingularException`
if the matrix is singular.
"""
function intermediate_to_pixel(wcs::WCSTransform{N}, intermediate::StaticVector{N}) where {N}
    T = _coordinate_float_type(intermediate)
    cpix_T = SVector{N, T}(wcs.crpix)

    # The ordinary affine path can invert the combined CD matrix directly.
    if !(wcs.projection isa TPV) && !has_sequent_distortion(wcs.pipeline)
        cd_T = SMatrix{N, N, T}(wcs.cd)
        focal = cpix_T .+ (cd_T \ intermediate)
        return focal_to_pixel(wcs.pipeline, focal, Val(N))
    end

    # Undo CDELT, sequent TPD, and TPV before solving the PC matrix.
    coord = intermediate ./ SVector{N, T}(wcs.cdelt)
    coord = invert_sequent_distortion(wcs.pipeline, coord)
    if wcs.projection isa TPV
        lon = wcs.lon_axis
        lat = wcs.lat_axis
        x, y = _tpv_inverse(wcs.projection.xcoeff, wcs.projection.ycoeff, coord[lon], coord[lat])
        coord = SVector{N, T}(ntuple(i -> i == lon ? T(x) : i == lat ? T(y) : coord[i], N))
    end
    pc_T = SMatrix{N, N, T}(wcs.pc)
    focal = cpix_T .+ (pc_T \ coord)

    return focal_to_pixel(wcs.pipeline, focal, Val(N))
end

function intermediate_to_pixel(wcs::WCSTransform{N}, intermediate::AbstractVector) where {N}
    length(intermediate) == N ||
        throw(DimensionMismatch("intermediate has length $(length(intermediate)), expected $N"))

    # Use a static temporary for the matrix solve, then return ordinary storage.
    T = _coordinate_float_type(intermediate)
    x = SVector{N,T}(ntuple(i -> T(intermediate[i]), N))
    return intermediate_to_pixel(wcs, x)
end
