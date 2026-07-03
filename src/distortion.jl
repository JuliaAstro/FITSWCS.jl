"""
SIP distortion parsing and evaluation.

The Simple Imaging Polynomial convention represents image-plane distortion as
polynomial offsets in pixel coordinates relative to `CRPIX1` and `CRPIX2`.
"""

"""
    SIPDistortion

Simple Imaging Polynomial distortion model.

SIP applies only to the first two pixel axes.  The forward coefficients `a`
and `b` map detector pixel coordinates to focal/image-plane pixel coordinates.
The optional inverse coefficients `ap` and `bp` map focal/image-plane
coordinates back to detector pixel coordinates.
"""
struct SIPDistortion
    crpix::SVector{2, Float64}
    a::Matrix{Float64}
    b::Matrix{Float64}
    ap::Union{Nothing, Matrix{Float64}}
    bp::Union{Nothing, Matrix{Float64}}
end

# ──────────────────────────────────────────────────────────────────────────────

"""Abstract supertype for pre-linear pixel/focal-plane distortion pipelines."""
abstract type AbstractDistortionPipeline end

"""Identity distortion pipeline for WCS transforms with no pre-linear distortion."""
struct NoDistortionPipeline <: AbstractDistortionPipeline end

"""
    DistortionPipeline

Pre-linear distortion pipeline.  Detector-to-image lookup tables are applied
first; SIP and CPDIS offsets are then evaluated at the detector-corrected
coordinate before the linear WCS matrix.
"""
struct DistortionPipeline{S <: Union{Nothing, SIPDistortion}, D <: Tuple, C <: Tuple} <: AbstractDistortionPipeline
    det2im::D
    sip::S
    cpdis::C
end

DistortionPipeline(sip::SIPDistortion) =
    DistortionPipeline{typeof(sip), Tuple{Nothing, Nothing}, Tuple{Nothing, Nothing}}((nothing, nothing), sip, (nothing, nothing))

function distortion_pipeline(sip::Union{Nothing, SIPDistortion}, ::NoAuxiliaryWCSData)
    # Header-only WCS uses the existing no-op or SIP-only pipeline.
    return distortion_pipeline(sip)
end

function distortion_pipeline(sip::Union{Nothing, SIPDistortion}, aux::AuxiliaryWCSData)
    det2im = aux.det2im isa Tuple && length(aux.det2im) == 2 ? aux.det2im : (nothing, nothing)
    cpdis = aux.cpdis isa Tuple && length(aux.cpdis) == 2 ? aux.cpdis : (nothing, nothing)

    # Avoid allocating a pipeline for auxiliary payloads that contain only TAB data.
    if isnothing(sip) && all(isnothing, det2im) && all(isnothing, cpdis)
        return NoDistortionPipeline()
    end

    # Preserve exact tuple types so no-lookup and partial-lookup stages stay concrete.
    return DistortionPipeline{typeof(sip), typeof(det2im), typeof(cpdis)}(det2im, sip, cpdis)
end

function has_sip_keywords(header::AbstractDict, alt_str::AbstractString)
    # Detect any SIP order keyword for this WCS version.
    for prefix in ("A", "B", "AP", "BP")
        haskey(header, "$(prefix)_ORDER$(alt_str)") && return true
    end
    return false
end

function read_sip_matrix(header::AbstractDict, prefix::AbstractString, order::Int, alt_str::AbstractString)
    coeff = zeros(Float64, order + 1, order + 1)

    # Fill only the triangular coefficient region defined by total order.
    for i in 0:order, j in 0:(order - i)
        key = "$(prefix)_$(i)_$(j)$(alt_str)"
        if haskey(header, key)
            coeff[i + 1, j + 1] = Float64(header[key])
        end
    end

    return coeff
end

function parse_sip_distortion(header::AbstractDict, crpix::Vector{Float64}, naxis::Int, alt::Char)
    alt_str = alt == ' ' ? "" : string(alt)
    has_sip_keywords(header, alt_str) || return nothing

    # SIP is defined for two image axes and requires explicit reference pixels.
    naxis >= 2 || throw(ArgumentError("SIP distortion requires at least two WCS axes"))
    for key in ("CRPIX1$(alt_str)", "CRPIX2$(alt_str)")
        haskey(header, key) || throw(ArgumentError("SIP distortion requires explicit $key"))
    end

    # Forward coefficients must be present as a matched A/B pair.
    has_a = haskey(header, "A_ORDER$(alt_str)")
    has_b = haskey(header, "B_ORDER$(alt_str)")
    has_a == has_b || throw(ArgumentError("SIP A_ORDER and B_ORDER must be provided together"))
    has_a || throw(ArgumentError("SIP distortion requires A_ORDER and B_ORDER"))

    a_order = Int(header["A_ORDER$(alt_str)"])
    b_order = Int(header["B_ORDER$(alt_str)"])
    a_order >= 0 || throw(ArgumentError("A_ORDER must be non-negative, got $a_order"))
    b_order >= 0 || throw(ArgumentError("B_ORDER must be non-negative, got $b_order"))
    a = read_sip_matrix(header, "A", a_order, alt_str)
    b = read_sip_matrix(header, "B", b_order, alt_str)

    # Inverse coefficients are optional but must be provided as an AP/BP pair.
    has_ap = haskey(header, "AP_ORDER$(alt_str)")
    has_bp = haskey(header, "BP_ORDER$(alt_str)")
    has_ap == has_bp || throw(ArgumentError("SIP AP_ORDER and BP_ORDER must be provided together"))
    ap = nothing
    bp = nothing
    if has_ap
        ap_order = Int(header["AP_ORDER$(alt_str)"])
        bp_order = Int(header["BP_ORDER$(alt_str)"])
        ap_order >= 0 || throw(ArgumentError("AP_ORDER must be non-negative, got $ap_order"))
        bp_order >= 0 || throw(ArgumentError("BP_ORDER must be non-negative, got $bp_order"))
        ap = read_sip_matrix(header, "AP", ap_order, alt_str)
        bp = read_sip_matrix(header, "BP", bp_order, alt_str)
    end

    return SIPDistortion(SVector{2, Float64}(crpix[1:2]), a, b, ap, bp)
end

distortion_pipeline(::Nothing) = NoDistortionPipeline()
distortion_pipeline(sip::SIPDistortion) = DistortionPipeline(sip)

has_distortion(::NoDistortionPipeline) = false
has_distortion(::DistortionPipeline) = true

function _lookup_stage_offset(tables::Tuple, coord::StaticVector{N, T}) where {N, T}
    if length(tables) != 2 || N < 2
        throw(ArgumentError("lookup stage requires two tables and at least two pixel axes, got $(length(tables)) tables and $N axes"))
    end
    x = coord[1]
    y = coord[2]

    # Paper IV image arrays store additive offsets for each corrected axis.
    dx = isnothing(tables[1]) ? zero(T) : T(tables[1](x, y))
    dy = isnothing(tables[2]) ? zero(T) : T(tables[2](x, y))
    return SVector{N, T}(ntuple(i -> i == 1 ? dx :
                                     i == 2 ? dy :
                                     zero(T), N))
end

function evaluate_sip_polynomial(coeff::AbstractMatrix, u::Real, v::Real)
    T = _promote_float_type(u, v)
    order = size(coeff, 1) - 1
    value = zero(T)

    # Sum terms whose total degree is within the SIP polynomial order.
    for i in 0:order, j in 0:(order - i)
        c = coeff[i + 1, j + 1]
        iszero(c) && continue
        value += T(c) * _smallpow(u, i) * _smallpow(v, j)
    end

    return value
end

sip_pixel_to_focal(::Nothing, pixel::AbstractVector) = SVector{2, _coordinate_float_type(pixel)}(pixel[1], pixel[2])
sip_pixel_to_focal(::Nothing, pixel::StaticVector) = pixel
function sip_pixel_to_focal(sip::SIPDistortion, pixel::AbstractVector)
    length(pixel) >= 2 || throw(DimensionMismatch("SIP distortion requires at least two pixel axes"))

    # Evaluate forward offsets relative to the SIP reference pixel.
    # FITS SIP convention: f_i = p_i + Σ A_ij (p_1−CRPIX1)^i (p_2−CRPIX2)^j
    u = pixel[1] - sip.crpix[1]
    v = pixel[2] - sip.crpix[2]
    fx = pixel[1] + evaluate_sip_polynomial(sip.a, u, v)
    fy = pixel[2] + evaluate_sip_polynomial(sip.b, u, v)

    return SVector{2, _coordinate_float_type(pixel)}(fx, fy)
end

function sip_focal_to_pixel(sip::SIPDistortion, focal::AbstractVector)
    length(focal) >= 2 || throw(DimensionMismatch("SIP distortion requires at least two pixel axes"))
    T = _coordinate_float_type(focal)

    # Prefer explicit inverse SIP coefficients when the header provides them.
    # FITS SIP convention: p_i = f_i + Σ AP_ij (f_1−CRPIX1)^i (f_2−CRPIX2)^j
    if sip.ap !== nothing && sip.bp !== nothing
        u = focal[1] - sip.crpix[1]
        v = focal[2] - sip.crpix[2]
        px = focal[1] + evaluate_sip_polynomial(sip.ap, u, v)
        py = focal[2] + evaluate_sip_polynomial(sip.bp, u, v)
        return SVector{2, T}(px, py)
    end

    # Otherwise solve forward(pixel) = focal with a fixed-point correction.
    # TODO: add a keyword argument (e.g. `error::Bool = false`) that raises
    # a `NoConvergence`-style exception carrying the best solution.
    pixel = SVector{2, T}(T(focal[1]), T(focal[2]))  # initial guess
    target = SVector{2, T}(T(focal[1]), T(focal[2]))
    max_iter = 64
    tol = _convergence_tol(T)
    prev_r = T(Inf)
    div_count = 0

    for k in 1:max_iter
        corrected = sip_pixel_to_focal(sip, pixel)
        dx = corrected[1] - target[1]
        dy = corrected[2] - target[2]
        pixel = SVector{2, T}(pixel[1] - dx, pixel[2] - dy)
        r = sum(abs2, (dx, dy))
        r <= tol^2 && return pixel

        if r >= prev_r
            div_count += 1
            if div_count >= 3
                @warn "SIP inverse is diverging at iteration $k " *
                    "(residual $(sqrt(prev_r)) → $(sqrt(r)) > tolerance $tol); " *
                    "returning best estimate so far"
                return pixel
            end
        else
            div_count = 0
        end
        prev_r = r
    end

    @warn "SIP inverse failed to converge after $max_iter iterations " *
        "(final residual $sqrt(prev_r) > tolerance $tol); " *
        "returning best estimate"
    return SVector{2, T}(pixel)
end

function pixel_to_focal(::NoDistortionPipeline, pixel::AbstractVector, ::Val{N}) where {N}
    length(pixel) == N ||
        throw(DimensionMismatch("pixel has length $(length(pixel)), expected $N"))

    # Materialize the identity focal coordinate in stable static storage.
    T = _coordinate_float_type(pixel)
    return SVector{N, T}(ntuple(i -> T(pixel[i]), N))
end

function pixel_to_focal(::NoDistortionPipeline, pixel::StaticVector{N}, ::Val{N}) where {N}
    # Static coordinates are already fixed-size, so identity distortion can return them directly.
    return pixel
end

function pixel_to_focal(pipeline::DistortionPipeline, pixel::AbstractVector, v::Val{N}) where {N}
    length(pixel) == N ||
        throw(DimensionMismatch("pixel has length $(length(pixel)), expected $N"))

    T = _coordinate_float_type(pixel)
    coord = SVector{N, T}(ntuple(i -> T(pixel[i]), N)) # Forward to function below
    return pixel_to_focal(pipeline, coord, v)
end

function pixel_to_focal(pipeline::DistortionPipeline, pixel::StaticVector{N}, ::Val{N}) where {N}
    T = _coordinate_float_type(pixel)

    # Evaluate all prior distortion offsets at the detector-corrected coordinate.
    detector = pixel + _lookup_stage_offset(pipeline.det2im, pixel)
    coord = detector
    if !isnothing(pipeline.sip)
        fx, fy = sip_pixel_to_focal(pipeline.sip, detector)
        coord = SVector{N, T}(ntuple(i -> i == 1 ? T(fx) :
                                          i == 2 ? T(fy) :
                                          detector[i], N))
    end
    return coord + _lookup_stage_offset(pipeline.cpdis, detector)
end

function focal_to_pixel(::NoDistortionPipeline, focal::AbstractVector, ::Val{N}) where {N}
    length(focal) == N ||
        throw(DimensionMismatch("focal coordinate has length $(length(focal)), expected $N"))

    # Materialize the identity pixel coordinate in stable static storage.
    T = _coordinate_float_type(focal)
    return SVector{N, T}(ntuple(i -> T(focal[i]), N))
end

function focal_to_pixel(::NoDistortionPipeline, focal::StaticVector{N}, ::Val{N}) where {N}
    # Static coordinates are already fixed-size, so identity inversion can return them directly.
    return focal
end

function focal_to_pixel(pipeline::DistortionPipeline, focal::AbstractVector, v::Val{N}) where {N}
    length(focal) == N ||
        throw(DimensionMismatch("focal coordinate has length $(length(focal)), expected $N"))

    T = _coordinate_float_type(focal)
    coord = SVector{N, T}(ntuple(i -> T(focal[i]), N)) # Forward to function below
    return focal_to_pixel(pipeline, coord, v)
end

# No Paper IV lookup stage is present, so we invert only SIP stage.
# This dispatch is needed to avoid allocations in the common case of SIP-only distortion.
function focal_to_pixel(pipeline::DistortionPipeline{S, Tuple{Nothing, Nothing}, Tuple{Nothing, Nothing}}, focal::StaticVector{N}, ::Val{N}) where {S, N}
    T = _coordinate_float_type(focal)

    # Preserve identity behavior for SIP-free pipeline variants.
    if isnothing(pipeline.sip)
        return SVector{N, T}(ntuple(i -> T(focal[i]), N))
    end

    # Invert SIP-only pipelines through the existing SIP inverse path.
    px, py = sip_focal_to_pixel(pipeline.sip, focal)
    return SVector{N,T}(ntuple(i ->
        i == 1 ? T(px) :
        i == 2 ? T(py) :
        T(focal[i]), N))
end

# Has Paper IV lookup stage, so we must iterate to invert the full pipeline.
function focal_to_pixel(pipeline::DistortionPipeline, focal::StaticVector{N}, ::Val{N}) where {N}
    T = _coordinate_float_type(focal)

    target = SVector{N, T}(ntuple(i -> T(focal[i]), N))
    pixel = target

    # Use the SIP inverse as a better starting point when it is available with fast inverse coefficients.
    if !isnothing(pipeline.sip) && !isnothing(pipeline.sip.ap) && !isnothing(pipeline.sip.bp)
        px, py = sip_focal_to_pixel(pipeline.sip, target)
        pixel = SVector{N, T}(ntuple(i -> i == 1 ? T(px) :
                                          i == 2 ? T(py) :
                                          target[i], N))
    end

    max_iter = 64
    # TODO: Consider reducing tolerance, iterations here are expensive
    tol = _convergence_tol(T)
    prev_r = T(Inf)
    div_count = 0

    for k in 1:max_iter
        # Correct the current estimate using the full forward distortion model.
        residual = pixel_to_focal(pipeline, pixel, Val(N)) - target
        r = sum(abs2, residual)
        r <= tol^2 && return pixel

        pixel = pixel - residual

        # Match SIP inverse behavior: warn and return a best effort if the solve diverges.
        if r > prev_r
            div_count += 1
            if div_count >= 3
                @warn "Paper IV lookup inverse is diverging at iteration $k " *
                    "(residual $(sqrt(prev_r)) → $(sqrt(r)) > tolerance $tol); " *
                    "returning best estimate so far"
                return pixel
            end
        else
            div_count = 0
        end
        prev_r = r
    end

    @warn "Paper IV lookup inverse failed to converge after $max_iter iterations " *
        "(final residual $(sqrt(prev_r)) > tolerance $tol); " *
        "returning best estimate"
    return pixel
end
