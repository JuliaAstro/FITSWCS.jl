"""
SIP distortion parsing and evaluation.

The Simple Imaging Polynomial convention represents image-plane distortion as
polynomial offsets in pixel coordinates relative to `CRPIX1` and `CRPIX2`.
"""


function has_sip_keywords(header::AbstractDict, alt_str::AbstractString)
    # Detect any SIP order keyword for this WCS version.
    for prefix in ("A", "B", "AP", "BP")
        haskey(header, "$(prefix)_ORDER$(alt_str)") && return true
    end
    return false
end

function is_lookup_distortion_keyword(key::AbstractString, alt_str::AbstractString="")
    ukey = uppercase(String(key))
    suffix = uppercase(String(alt_str))

    # Astropy exposes Paper IV lookup tables as CPDIS and detector-to-image D2IM.
    occursin(Regex("^CPDIS[1-9][0-9]*$(suffix)\$"), ukey) && return true
    occursin(Regex("^D2IMDIS[1-9][0-9]*$(suffix)\$"), ukey) && return true
    occursin(Regex("^D2IMERR[1-9][0-9]*$(suffix)\$"), ukey) && return true
    ukey == "AXISCORR$(suffix)" && return true

    # wcslib Paper IV distortion parameters use DPja/DQia keyword families.
    occursin(Regex("^D[QP][1-9][0-9]*$(suffix)\\."), ukey) && return true

    return false
end

function reject_lookup_distortion_keywords(header::AbstractDict, alt_str::AbstractString="")
    # Lookup-table distortions affect only the selected WCS version's pixel pipeline.
    for key in keys(header)
        key isa AbstractString || continue
        if is_lookup_distortion_keyword(key, alt_str)
            throw(ArgumentError(
                "distortion lookup keyword $key is not implemented yet"
            ))
        end
    end

    return nothing
end

function read_sip_matrix(header::AbstractDict, prefix::AbstractString,
                          order::Int, alt_str::AbstractString)
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

function parse_sip_distortion(header::AbstractDict, crpix::Vector{Float64},
                               naxis::Int, alt::Char)
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

    return SIPDistortion(SVector{2,Float64}(crpix[1:2]), a, b, ap, bp)
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

function sip_pixel_to_focal(sip::SIPDistortion, pixel::AbstractVector)
    length(pixel) >= 2 || throw(DimensionMismatch("SIP distortion requires at least two pixel axes"))

    # Evaluate forward offsets relative to the SIP reference pixel.
    # FITS SIP convention: f_i = p_i + Σ A_ij (p_1−CRPIX1)^i (p_2−CRPIX2)^j
    u = pixel[1] - sip.crpix[1]
    v = pixel[2] - sip.crpix[2]
    fx = pixel[1] + evaluate_sip_polynomial(sip.a, u, v)
    fy = pixel[2] + evaluate_sip_polynomial(sip.b, u, v)

    return SVector{2,_coordinate_float_type(pixel)}(fx, fy)
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
        return SVector{2,T}(px, py)
    end

    # Otherwise solve forward(pixel) = focal with a fixed-point correction.
    # TODO: add a keyword argument (e.g. `error::Bool = false`) that raises
    # a `NoConvergence`-style exception carrying the best solution.
    pixel = SVector{2,T}(T(focal[1]), T(focal[2]))  # initial guess
    target = SVector{2,T}(T(focal[1]), T(focal[2]))
    max_iter = 64
    tol = 1e-10
    prev_r = Inf
    div_count = 0

    for k in 1:max_iter
        corrected = sip_pixel_to_focal(sip, pixel)
        dx = corrected[1] - target[1]
        dy = corrected[2] - target[2]
        pixel = SVector{2,T}(pixel[1] - dx, pixel[2] - dy)
        r = hypot(dx, dy)
        r <= tol && return pixel

        if r >= prev_r
            div_count += 1
            if div_count >= 3
                @warn "SIP inverse is diverging at iteration $k " *
                      "(residual $prev_r → $r > tolerance $tol); " *
                      "returning best estimate so far"
                return pixel
            end
        else
            div_count = 0
        end
        prev_r = r
    end

    @warn "SIP inverse failed to converge after $max_iter iterations " *
          "(final residual $prev_r > tolerance $tol); " *
          "returning best estimate"
    return SVector{2,T}(pixel)
end
