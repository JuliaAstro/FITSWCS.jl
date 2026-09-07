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

"""
    TPDDistortion

One axis of a Paper IV Template Polynomial Distortion.  `axmap`, `offset`, and
`scale` define the independent variables.  Forward and optional reverse
coefficients use the 60-term TPD basis.  By default the polynomial is an
additive correction; `docorr=false` makes it return the corrected coordinate
directly.
"""
struct TPDDistortion{M, A <: Union{Nothing, SMatrix{2, 3, Float64, 6}}}
    axmap::SVector{M, Int}
    offset::SVector{M, Float64}
    scale::SVector{M, Float64}
    aux::A
    forward::Vector{Float64}
    reverse::Union{Nothing, Vector{Float64}}
    docorr::Bool
end

# ──────────────────────────────────────────────────────────────────────────────

"""Abstract supertype for distortion pipelines surrounding the linear transform."""
abstract type AbstractDistortionPipeline end

"""Identity distortion pipeline for WCS transforms with no distortion stages."""
struct NoDistortionPipeline <: AbstractDistortionPipeline end

"""
    DistortionPipeline

Distortion pipeline surrounding the linear WCS matrix.  Detector-to-image
lookup tables, SIP, CPDIS lookup offsets, and prior TPD are applied before PC.
Sequent TPD is applied after PC and before CDELT.
"""
struct DistortionPipeline{S <: Union{Nothing, SIPDistortion}, D <: Tuple, C <: Tuple, P <: Tuple, Q <: Tuple} <: AbstractDistortionPipeline
    det2im::D
    sip::S
    cpdis::C
    tpd_pre::P
    tpd_seq::Q
end

DistortionPipeline(sip::SIPDistortion) =
    DistortionPipeline((nothing, nothing), sip, (nothing, nothing), (), ())

function distortion_pipeline(sip::Union{Nothing, SIPDistortion}, aux::AbstractAuxiliaryWCSData, tpd_pre::Tuple, tpd_seq::Tuple)
    # External lookup tables are absent from header-only auxiliary payloads.
    det2im = aux isa AuxiliaryWCSData && aux.det2im isa Tuple && length(aux.det2im) == 2 ? aux.det2im : (nothing, nothing)
    cpdis = aux isa AuxiliaryWCSData && aux.cpdis isa Tuple && length(aux.cpdis) == 2 ? aux.cpdis : (nothing, nothing)

    # Avoid allocating a pipeline when every distortion stage is empty.
    if isnothing(sip) && all(isnothing, det2im) && all(isnothing, cpdis) &&
            all(isnothing, tpd_pre) && all(isnothing, tpd_seq)
        return NoDistortionPipeline()
    end

    # Preserve empty TPD stages explicitly for cheap common-case dispatch.
    prior = all(isnothing, tpd_pre) ? () : tpd_pre
    sequent = all(isnothing, tpd_seq) ? () : tpd_seq
    return DistortionPipeline(det2im, sip, cpdis, prior, sequent)
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

function _collect_tpd_coefficients(header::AbstractDict, prefix::AbstractString, direction::AbstractString)
    coeff = Float64[]

    # Preserve coefficient numbering while zero-filling omitted terms.
    for m in 0:59
        key = "$(prefix).TPD.$(direction).$(m)"
        haskey(header, key) || continue
        while length(coeff) < m
            push!(coeff, 0.0)
        end
        push!(coeff, Float64(header[key]))
    end
    return coeff
end

function _tpd_auxiliary_matrix(header::AbstractDict, prefix::AbstractString)
    has_aux = any(haskey(header, "$(prefix).AUX.$(k).COEFF.$(m)") for k in 1:2, m in 0:2)
    has_aux || return nothing

    # WCSLIB defaults the auxiliary transform to the two-dimensional identity.
    aux = zeros(Float64, 2, 3)
    aux[1, 2] = 1.0
    aux[2, 3] = 1.0
    for k in 1:2, m in 0:2
        key = "$(prefix).AUX.$(k).COEFF.$(m)"
        haskey(header, key) && (aux[k, m + 1] = Float64(header[key]))
    end
    return SMatrix{2, 3, Float64, 6}(aux)
end

function _valid_tpd_parameter_field(field::AbstractString)
    field in ("DOCORR", "NAXES") && return true
    occursin(r"^(AXIS|OFFSET|SCALE)\.[1-9][0-9]*$", field) && return true
    occursin(r"^TPD\.(FWD|REV)\.([0-9]|[1-5][0-9])$", field) && return true
    return occursin(r"^AUX\.[12]\.COEFF\.[0-2]$", field)
end

function _build_tpd_distortion(header::AbstractDict, prefix::AbstractString, naxis::Int, ::Val{M}) where {M}
    # Resolve the independent-axis map and its normalization parameters.
    axmap = SVector{M, Int}(ntuple(k -> Int(get(header, "$(prefix).AXIS.$(k)", k)), M))
    all(i -> 1 <= i <= naxis, axmap) || throw(ArgumentError("$prefix axis map contains an axis outside 1:$naxis"))
    length(unique(axmap)) == M || throw(ArgumentError("$prefix axis map contains duplicate axes"))
    offset = SVector{M, Float64}(ntuple(k -> Float64(get(header, "$(prefix).OFFSET.$(k)", 0.0)), M))
    scale = SVector{M, Float64}(ntuple(k -> Float64(get(header, "$(prefix).SCALE.$(k)", 1.0)), M))

    # TPD forward coefficients are required structurally but may all be zero.
    forward = _collect_tpd_coefficients(header, prefix, "FWD")
    reverse_raw = _collect_tpd_coefficients(header, prefix, "REV")
    reverse = isempty(reverse_raw) ? nothing : reverse_raw
    aux = _tpd_auxiliary_matrix(header, prefix)
    docorr = Int(get(header, "$(prefix).DOCORR", 1)) != 0
    return TPDDistortion(axmap, offset, scale, aux, forward, reverse, docorr)
end

function _parse_tpd_axis(header::AbstractDict, axis::Int, naxis::Int, alt_str::AbstractString,
                         dist_prefix::AbstractString, param_prefix::AbstractString,
                         error_prefix::AbstractString, minerr::Real)
    dist_key = "$(dist_prefix)$(axis)$(alt_str)"
    haskey(header, dist_key) || return nothing
    dtype = uppercase(strip(String(header[dist_key])))
    if dtype == "LOOKUP"
        dist_prefix == "CPDIS" && return nothing
        throw(ArgumentError("sequent LOOKUP distortion $dist_key is not supported"))
    end
    dtype == "TPD" || throw(ArgumentError("unsupported Paper IV distortion type $(header[dist_key]) in $dist_key"))

    # Error thresholds suppress the complete distortion on this axis.
    error_key = "$(error_prefix)$(axis)$(alt_str)"
    Float64(get(header, error_key, 0.0)) < Float64(minerr) && return nothing

    prefix = "$(param_prefix)$(axis)$(alt_str)"
    naxes_key = "$(prefix).NAXES"
    haskey(header, naxes_key) || throw(ArgumentError("TPD distortion $dist_key requires $naxes_key"))
    nhat = Int(header[naxes_key])
    nhat in (1, 2) || throw(ArgumentError("$naxes_key must be 1 or 2, got $nhat"))

    # Reject misspelled or unsupported TPD parameter records.
    dotted_prefix = "$(prefix)."
    for key in keys(header)
        key isa AbstractString || continue
        startswith(key, dotted_prefix) || continue
        field = key[length(dotted_prefix) + 1:end]
        _valid_tpd_parameter_field(field) || throw(ArgumentError("unrecognized TPD parameter $key"))
        mapped = match(r"^(AXIS|OFFSET|SCALE)\.([1-9][0-9]*)$", field)
        mapped !== nothing && parse(Int, mapped.captures[2]) > nhat &&
            throw(ArgumentError("$key exceeds the $nhat independent axes declared by $naxes_key"))
    end

    return nhat == 1 ? _build_tpd_distortion(header, prefix, naxis, Val(1)) :
                       _build_tpd_distortion(header, prefix, naxis, Val(2))
end

function parse_tpd_distortions(header::AbstractDict, naxis::Int, alt::Char, minerr::Real)
    alt_str = alt == ' ' ? "" : string(alt)

    # CPDIS/DP is prior to PC; CQDIS/DQ is sequent to PC and prior to CDELT.
    prior = ntuple(axis -> _parse_tpd_axis(header, axis, naxis, alt_str, "CPDIS", "DP", "CPERR", minerr), naxis)
    sequent = ntuple(axis -> _parse_tpd_axis(header, axis, naxis, alt_str, "CQDIS", "DQ", "CQERR", minerr), naxis)
    return prior, sequent
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

function _evaluate_tpd_coefficients(coeff::AbstractVector, u::T, v::T, ::Val{M}) where {T, M}
    result = zero(T)

    # One-dimensional TPD uses only pure powers of its attached coordinate.
    @inbounds for i in eachindex(coeff)
        c = coeff[i]
        iszero(c) && continue
        term = _TPD_TERMS[i]
        if M == 1
            term[1] === :mono && term[3] == 0 || continue
        end
        result += T(c) * _tpd_term_value(i - 1, u, v)
    end
    return result
end

function _evaluate_tpd(model::TPDDistortion{M}, raw::StaticVector, coeff::AbstractVector) where {M}
    T = _coordinate_float_type(raw)

    # Normalize and reorder the independent coordinates through the axis map.
    vars = SVector{M, T}(ntuple(k -> (T(raw[model.axmap[k]]) - T(model.offset[k])) * T(model.scale[k]), M))
    u = vars[1]
    v = M == 1 ? zero(T) : vars[2]

    # Optional auxiliary variables are a two-dimensional affine transform.
    if !isnothing(model.aux)
        aux = SMatrix{2, 3, T, 6}(model.aux)
        transformed = aux * SVector{3, T}(one(T), u, v)
        u, v = transformed
    end
    return _evaluate_tpd_coefficients(coeff, u, v, Val(M))
end

apply_tpd_stage(::Tuple{}, raw::StaticVector) = raw

function apply_tpd_stage(models::Tuple, raw::StaticVector{N, T}) where {N, T}
    length(models) == N || throw(DimensionMismatch("TPD stage has $(length(models)) axes, expected $N"))

    # Every output is evaluated from the same undistorted input coordinate.
    return SVector{N, T}(ntuple(j -> begin
        model = models[j]
        if isnothing(model)
            raw[j]
        else
            value = _evaluate_tpd(model, raw, model.forward)
            model.docorr ? raw[j] + value : value
        end
    end, N))
end

function _tpd_reverse_guess(models::Tuple, target::StaticVector{N, T}) where {N, T}
    # Reverse polynomials provide only the starting point for forward iteration.
    return SVector{N, T}(ntuple(j -> begin
        model = models[j]
        if isnothing(model) || isnothing(model.reverse)
            target[j]
        else
            value = _evaluate_tpd(model, target, model.reverse)
            model.docorr ? target[j] + value : value
        end
    end, N))
end

invert_tpd_stage(::Tuple{}, target::StaticVector) = target

function invert_tpd_stage(models::Tuple, target::StaticVector{N, T}) where {N, T}
    raw = _tpd_reverse_guess(models, target)
    tol = T(_convergence_tol(T))

    # Refine the reverse-polynomial estimate using the forward model and a numerical Jacobian.
    for _ in 1:30
        distorted = apply_tpd_stage(models, raw)
        residual = distorted - target
        all(j -> abs(residual[j]) <= tol * max(one(T), abs(target[j])), 1:N) && return raw

        jacobian = MMatrix{N, N, T}(undef)
        for column in 1:N
            step = clamp(abs(residual[column]) / 2, T(1e-6), one(T))
            trial = Base.setindex(raw, raw[column] + step, column)
            shifted = apply_tpd_stage(models, trial)
            for row in 1:N
                jacobian[row, column] = (shifted[row] - distorted[row]) / step
            end
        end
        raw -= SMatrix{N, N, T}(jacobian) \ residual
    end

    residual = apply_tpd_stage(models, raw) - target
    @warn "TPD inverse failed to converge after 30 iterations (residual $(sqrt(sum(abs2, residual))) > tolerance $tol); returning best estimate"
    return raw
end

has_sequent_distortion(::NoDistortionPipeline) = false
has_sequent_distortion(pipeline::DistortionPipeline) = any(!isnothing, pipeline.tpd_seq)

apply_sequent_distortion(::NoDistortionPipeline, coord::StaticVector) = coord
apply_sequent_distortion(pipeline::DistortionPipeline, coord::StaticVector) = apply_tpd_stage(pipeline.tpd_seq, coord)

invert_sequent_distortion(::NoDistortionPipeline, coord::StaticVector) = coord
invert_sequent_distortion(pipeline::DistortionPipeline, coord::StaticVector) = invert_tpd_stage(pipeline.tpd_seq, coord)

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
    lookup_coord = coord + _lookup_stage_offset(pipeline.cpdis, detector)
    isempty(pipeline.tpd_pre) && return lookup_coord

    # CPDIS functions on different output axes are evaluated from the same input.
    tpd_coord = apply_tpd_stage(pipeline.tpd_pre, coord)
    return SVector{N, T}(ntuple(i -> isnothing(pipeline.tpd_pre[i]) ? lookup_coord[i] : tpd_coord[i], N))
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
function focal_to_pixel(pipeline::DistortionPipeline{S, Tuple{Nothing, Nothing}, Tuple{Nothing, Nothing}, Tuple{}, Q}, focal::StaticVector{N}, ::Val{N}) where {S, Q, N}
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

# A TPD-only prior stage can use its reverse polynomial and Newton refinement directly.
function focal_to_pixel(pipeline::DistortionPipeline{Nothing, Tuple{Nothing, Nothing}, Tuple{Nothing, Nothing}, P, Q}, focal::StaticVector{N}, ::Val{N}) where {P <: Tuple{Any, Vararg{Any}}, Q, N}
    T = _coordinate_float_type(focal)
    target = SVector{N, T}(focal)
    return invert_tpd_stage(pipeline.tpd_pre, target)
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
