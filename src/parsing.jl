"""
FITS WCS header parsing.

Converts a dictionary-like FITS header into a `WCSTransform`.

Supported keyword forms (Paper I, Section 2):
- `CTYPEia`, `CRPIXia`, `CRVALia`, `CDELTia`, `CUNITia`
- `PCi_ja` (linear transform matrix; multiplied by CDELT to form CD)
- `CDi_ja` (explicit CD matrix; takes precedence over PC/CDELT)
- `CROTAi` (legacy rotation keyword; converted to PC form)
- `WCSAXES`, `NAXIS`
- `LONPOLE`, `LATPOLE`
- `EQUINOX`, `RADESYS`, `MJD-OBS`, `DATE-OBS`

The suffix `a` denotes an alternate WCS description character ('A'–'Z').
The primary description uses a blank suffix.
"""

"""
    WCSTransform

Parsed, validated FITS World Coordinate System transform.

## Coordinate conventions

Pixel coordinates follow the **FITS 1-based** convention: pixel 1 is the
centre of the first array element.  This matches the FITS standard and the
values stored in `CRPIX` header keywords.  When working with Julia 1-based
array indices the values are numerically identical; no offset is required.

## Fields

- `naxis`      – number of WCS axes.
- `crpix`      – reference pixel position (FITS 1-based), length `naxis`.
- `crval`      – world coordinate value at the reference pixel, length `naxis`.
- `cd`         – combined CD matrix (naxis × naxis), where
                 `cd[i,j] = CDELT_i * PC_i_j` (or the explicit `CD_i_j` value).
                 Units match `cunit`.
- `ctype`      – FITS `CTYPEi` strings, length `naxis`.
- `cunit`      – FITS `CUNITi` strings (empty string means degrees for
                 celestial axes), length `naxis`.
- `lonpole`    – native longitude of the celestial pole (degrees), φₚ in Paper II.
- `latpole`    – native latitude of the celestial pole hint (degrees), θₚ in
                 Paper II; used to resolve the `delta_p` ambiguity.
- `alpha_p`    – celestial longitude of the native north pole (degrees), αₚ in
                 Paper II.  Precomputed during construction; not a direct FITS
                 keyword.
- `delta_p`    – celestial latitude of the native north pole (degrees), δₚ in
                 Paper II.  Precomputed during construction.
- `projection` – spherical projection for the celestial axes, or `nothing` for
                 purely linear WCS.
- `pipeline`   – pre-linear pixel/focal-plane distortion pipeline.
- `aux`        – resolved auxiliary WCS data, or `NoAuxiliaryWCSData`.
- `lon_axis`   – 1-based index of the longitude axis; 0 if no celestial axes.
- `lat_axis`   – 1-based index of the latitude axis; 0 if no celestial axes.
"""
struct WCSTransform{N, L, P <: Union{Nothing, AbstractProjection}, D <: AbstractDistortionPipeline, A <: AbstractAuxiliaryWCSData}
    naxis::Int
    crpix::SVector{N, Float64}
    crval::SVector{N, Float64}
    cd::SMatrix{N, N, Float64, L}       # naxis × naxis
    ctype::Vector{String}
    cunit::Vector{String}
    lonpole::Float64          # degrees; φₚ (Paper II)
    latpole::Float64          # degrees; used during construction
    alpha_p::Float64          # degrees; celestial lon of native N pole
    delta_p::Float64          # degrees; celestial lat of native N pole
    projection::P
    pipeline::D
    aux::A
    lon_axis::Int             # 1-based index; 0 = no celestial lon axis
    lat_axis::Int             # 1-based index; 0 = no celestial lat axis
end

# ──────────────────────────────────────────────────────────────────────────────
# CTYPE parsing helpers
# ──────────────────────────────────────────────────────────────────────────────

"""
Known longitude axis coordinate system prefixes.
"""
const _LON_SYSTEMS = Set(["RA", "GLON", "ELON", "HLON", "SLON", "HPLN"])

"""
Known latitude axis coordinate system prefixes.
"""
const _LAT_SYSTEMS = Set(["DEC", "GLAT", "ELAT", "HLAT", "SLAT", "HPLT"])

"""
    parse_ctype(ctype_str) -> (system, coord_type, proj_code)

Parse a FITS `CTYPEi` string into its components.

The FITS WCS convention (Paper I, Section 2.4) uses the format:
`'AAAA-BBB'` where `AAAA` is the coordinate type (padded to 4 characters
with hyphens) and `BBB` is the three-character projection code.

Returns:
- `system`     – the coordinate system string (e.g. `"RA"`, `"DEC"`, `"WAVE"`)
- `coord_type` – `:lon`, `:lat`, `:linear`, or `:unknown`
- `proj_code`  – the projection code string, or `""` for non-celestial axes
"""
function parse_ctype(ctype_str::AbstractString)
    s = uppercase(strip(ctype_str))
    if isempty(s)
        return "", :linear, ""
    end

    # Recognize the strict FITS 4-3 form, preserving the base code for suffixes.
    if ncodeunits(s) >= 8 && s[5] == '-'
        coord_part = rstrip(s[1:4], '-')
        proj_code = s[6:8]
    else
        coord_part = s
        proj_code = ""
    end

    coord_type = if coord_part in _LON_SYSTEMS
        :lon
    elseif coord_part in _LAT_SYSTEMS
        :lat
    else
        :linear
    end

    return coord_part, coord_type, proj_code
end

"""
    projection_from_code(code::AbstractString) -> AbstractProjection

Create a projection object from a three-character FITS projection code.
Returns `UnknownProjection(code)` for unrecognised codes rather than throwing.
"""
function projection_from_code(code::AbstractString)
    c = uppercase(strip(code))
    c == "AZP" && return AZP()
    c == "SZP" && return SZP()
    c == "TAN" && return TAN()
    c == "SIN" && return SIN()
    c == "STG" && return STG()
    c == "ARC" && return ARC()
    c == "ZEA" && return ZEA()
    c == "CYP" && return CYP()
    c == "MER" && return MER()
    c == "SFL" && return SFL()
    c == "PAR" && return PAR()
    c == "MOL" && return MOL()
    c == "PCO" && return PCO()
    c == "CAR" && return CAR()
    c == "CEA" && return CEA()
    c == "AIT" && return AIT()
    c == "ZPN" && return ZPN()
    c == "AIR" && return AIR()
    c == "BON" && return BON(45.0)   # will be overridden by projection_from_header
    c == "COP" && return COP(45.0, 0.0)
    c == "COD" && return COD(45.0, 0.0)
    c == "COE" && return COE(45.0, 0.0)
    c == "COO" && return COO(45.0, 0.0)
    c == "TSC" && return TSC()
    c == "CSC" && return CSC()
    c == "QSC" && return QSC()
    c == "HPX" && return HPX()
    c == "XPH" && return XPH()
    c == "TPV" && return TPV()
    c == "TPD" && return TPV()
    return UnknownProjection(c)
end

function projection_from_header(code::AbstractString, header::AbstractDict,
                                 lon_axis::Int, lat_axis::Int, alt::Char)
    c = uppercase(strip(code))

    if c == "AZP"
        alt_str = alt == ' ' ? "" : string(alt)
        mu = Float64(get(header, "PV$(lat_axis)_1$(alt_str)", 0.0))
        gamma = Float64(get(header, "PV$(lat_axis)_2$(alt_str)", 0.0))
        return AZP(mu, gamma)
    end

    if c == "SZP"
        alt_str = alt == ' ' ? "" : string(alt)
        mu = Float64(get(header, "PV$(lat_axis)_1$(alt_str)", 0.0))
        phi_c = Float64(get(header, "PV$(lat_axis)_2$(alt_str)", 0.0))
        theta_c = Float64(get(header, "PV$(lat_axis)_3$(alt_str)", 90.0))
        return SZP(mu, phi_c, theta_c)
    end

    # SIN uses slant orthographic parameters on the latitude-like axis.
    if c == "SIN"
        alt_str = alt == ' ' ? "" : string(alt)
        xi = Float64(get(header, "PV$(lat_axis)_1$(alt_str)", 0.0))
        eta = Float64(get(header, "PV$(lat_axis)_2$(alt_str)", 0.0))
        return SIN(xi, eta)
    end

    # CEA uses PV latitude-axis parameter 1 as the cylindrical equal-area lambda.
    if c == "CEA"
        alt_str = alt == ' ' ? "" : string(alt)
        lambda = Float64(get(header, "PV$(lat_axis)_1$(alt_str)", 1.0))
        (0.0 < lambda <= 1.0) ||
            throw(ArgumentError("CEA PV$(lat_axis)_1 must be in (0, 1], got $lambda"))
        return CEA(lambda)
    end

    # CYP uses latitude-axis parameters 1 and 2 for lambda and mu.
    if c == "CYP"
        alt_str = alt == ' ' ? "" : string(alt)
        lambda = Float64(get(header, "PV$(lat_axis)_1$(alt_str)", 1.0))
        mu = Float64(get(header, "PV$(lat_axis)_2$(alt_str)", 1.0))
        lambda != 0.0 ||
            throw(ArgumentError("CYP PV$(lat_axis)_1 must be non-zero"))
        mu + lambda != 0.0 ||
            throw(ArgumentError("CYP PV$(lat_axis)_1 + PV$(lat_axis)_2 must be non-zero"))
        return CYP(lambda, mu)
    end

    # ZPN uses latitude-axis parameters 0..N as polynomial coefficients.
    if c == "ZPN"
        alt_str = alt == ' ' ? "" : string(alt)
        # Collect coefficients up to m=30 (WCSLIB limit).
        pvs = Float64[]
        for m in 0:30
            key = "PV$(lat_axis)_$(m)$(alt_str)"
            if haskey(header, key)
                # Extend vector if needed, fill gaps with 0.
                while length(pvs) < m
                    push!(pvs, 0.0)
                end
                push!(pvs, Float64(header[key]))
            end
        end
        isempty(pvs) && (pvs = [0.0, 1.0])  # default: r = zd (same as ARC)
        return ZPN(pvs)
    end

    # AIR uses latitude-axis parameter 1 as theta_b (break latitude, degrees).
    if c == "AIR"
        alt_str = alt == ' ' ? "" : string(alt)
        theta_b = Float64(get(header, "PV$(lat_axis)_1$(alt_str)", 90.0))
        return AIR(theta_b)
    end

    # BON uses latitude-axis parameter 1 as theta1 (standard parallel, degrees).
    if c == "BON"
        alt_str = alt == ' ' ? "" : string(alt)
        theta1 = Float64(get(header, "PV$(lat_axis)_1$(alt_str)", 45.0))
        return BON(theta1)
    end

    # Conic projections: PVlat_1 = sigma (deg), PVlat_2 = delta (deg, default 0).
    if c in ("COP", "COD", "COE", "COO")
        alt_str = alt == ' ' ? "" : string(alt)
        sigma = Float64(get(header, "PV$(lat_axis)_1$(alt_str)", 45.0))
        delta = Float64(get(header, "PV$(lat_axis)_2$(alt_str)", 0.0))
        c == "COP" && return COP(sigma, delta)
        c == "COD" && return COD(sigma, delta)
        c == "COE" && return COE(sigma, delta)
        c == "COO" && return COO(sigma, delta)
    end

    # HPX uses PVlat_1 = H (default 4) and PVlat_2 = K (default 3).
    if c == "HPX"
        alt_str = alt == ' ' ? "" : string(alt)
        H = Int(get(header, "PV$(lat_axis)_1$(alt_str)", 4))
        K = Int(get(header, "PV$(lat_axis)_2$(alt_str)", 3))
        return HPX(H, K)
    end

    # TPV/TPD: TAN + sequent polynomial distortion encoded in PVi_m keywords.
    # Coefficients are collected from PV keywords on both celestial axes.
    if c == "TPV" || c == "TPD"
        alt_str = alt == ' ' ? "" : string(alt)
        xcoeff_raw = _collect_tpv_coeffs(header, lon_axis, alt_str)
        ycoeff_raw = _collect_tpv_coeffs(header, lat_axis, alt_str)
        # If both empty, return identity TPV (= plain TAN).
        xcoeff = isempty(xcoeff_raw) ? Float64[0.0, 1.0] : xcoeff_raw
        ycoeff = isempty(ycoeff_raw) ? Float64[0.0, 0.0, 1.0] : ycoeff_raw
        return TPV(xcoeff, ycoeff)
    end

    # Other supported projections currently have no parsed PV parameters.
    return projection_from_code(c)
end

"""
    _collect_tpv_coeffs(header, axis, alt_str) -> Vector{Float64}

Collect TPV/TPD polynomial coefficients from `PV{axis}_m` keywords for
``m = 0..59`` (the TPD coefficient range).  Gaps are zero-filled so the
returned vector is indexed directly by ``m`` (i.e. `result[m+1] = PVm`).
Returns an empty vector if no PV keywords are found on this axis.
"""
function _collect_tpv_coeffs(header::AbstractDict, axis::Int, alt_str::AbstractString)
    coeff = Float64[]
    for m in 0:59
        key = "PV$(axis)_$(m)$(alt_str)"
        if haskey(header, key)
            while length(coeff) < m
                push!(coeff, 0.0)   # fill gap with zero
            end
            push!(coeff, Float64(header[key]))
        end
    end
    return coeff
end

# ── SCAMP TPV compatibility ────────────────────────────────────────────────────

"""
    _remove_sip_keywords!(header, alt_str)

Delete all SIP-distortion keywords (`A_ORDER`, `B_ORDER`, `A_i_j`, `B_i_j`,
`AP_*`, `BP_*`) for the given alternate WCS suffix from `header`.
"""
function _remove_sip_keywords!(header::AbstractDict, alt_str::AbstractString)
    sip_prefixes = ("A_", "B_", "AP_", "BP_")
    to_remove = String[]
    for (key, _) in header
        key isa AbstractString || continue
        u = uppercase(key)
        for pfx in sip_prefixes
            startswith(u, pfx) || continue
            (occursin(r"^(A|B|AP|BP)_(ORDER|[0-9]+_[0-9]+)$" * alt_str, u) ||
             occursin(r"^(A|B|AP|BP)_ORDER$" * alt_str, u)) &&
                push!(to_remove, key)
        end
    end
    for k in to_remove
        delete!(header, k)
    end
end

"""
    _detect_pre2012_scamp_tpv(header, alt_str) -> Bool

Detect pre-2012 SCAMP headers: CTYPE ends in `-TAN` on celestial axes AND
PVi_j keywords with j ≥ 5 are present on those axes.
"""
function _detect_pre2012_scamp_tpv(header::AbstractDict, alt_str::AbstractString)
    tan_axes = Int[]
    pv_axes_high = Int[]

    for (key, value) in header
        key isa AbstractString || continue
        ukey = uppercase(key)

        # Detect CTYPE with -TAN suffix
        m = match(Regex("^CTYPE([1-9][0-9]*)$(alt_str)\$"), ukey)
        if m !== nothing
            sys, ct, pc = parse_ctype(String(value))
            if pc == "TAN"
                push!(tan_axes, parse(Int, m.captures[1]))
            end
            continue
        end

        # Detect PVi_j with j ≥ 5
        m = match(Regex("^PV([1-9][0-9]*)_([0-9]+)$(alt_str)\$"), ukey)
        if m !== nothing
            param = parse(Int, m.captures[2])
            if param >= 5
                push!(pv_axes_high, parse(Int, m.captures[1]))
            end
        end
    end

    !isempty(tan_axes) && !isempty(pv_axes_high) && return true
    return false
end

"""
    _fix_scamp_compatibility!(header, alt_str)

Apply SCAMP WCS compatibility fixes to `header` before parsing.

Pre-2012 SCAMP wrote TPV coefficients in `PVi_m` keywords but used
`CTYPE=-TAN`.  When detected, the CTYPE is mutated from `-TAN` to
`-TPV` and any SIP keywords are stripped so that TPV is parsed as
the active distortion model.
"""
function _fix_scamp_compatibility!(header::AbstractDict, alt_str::AbstractString)
    if !_detect_pre2012_scamp_tpv(header, alt_str)
        return
    end

    # Mutate CTYPE from -TAN to -TPV on celestial axes.
    for (key, value) in header
        key isa AbstractString || continue
        ukey = uppercase(key)
        m = match(Regex("^CTYPE([1-9][0-9]*)$(alt_str)\$"), ukey)
        m === nothing && continue
        sys, ct, pc = parse_ctype(String(value))
        if pc == "TAN"
            # Preserve the original 4-char prefix including hyphens (e.g. "RA--").
            prefix = String(value)[1:4]
            header[key] = "$(prefix)-TPV"
        end
    end

    # Remove SIP keywords (TPV takes precedence).
    _remove_sip_keywords!(header, alt_str)
end

function reject_unsupported_nonlinear_axes(ctype::Vector{String})
    # Non-celestial algorithm-coded axes need Paper III or domain-specific logic.
    for (i, ctype_i) in pairs(ctype)
        _, coord_type, pc = parse_ctype(ctype_i)
        if coord_type == :linear && !isempty(pc) && pc != "TAB"
            throw(ArgumentError(
                "CTYPE$(i)=$(repr(ctype_i)) uses unsupported algorithm code " *
                "$(repr(pc)); only plain linear non-celestial axes are implemented"
            ))
        end
    end

    return nothing
end

function _wcs_axis_indices(key::AbstractString, alt_str::AbstractString)
    if !endswith(key, alt_str)
        return Int[]
    end
    if !isempty(alt_str) && ncodeunits(key) <= ncodeunits(alt_str)
        return Int[]
    end

    # Strip the alternate suffix before matching indexed WCS keyword forms.
    core = isempty(alt_str) ? key : key[1:(lastindex(key) - ncodeunits(alt_str))]

    # Single-index keywords carry one WCS axis number.
    m = match(r"^(?:CRPIX|CRVAL|CDELT|CTYPE|CUNIT)([1-9][0-9]*)$", core)
    if m !== nothing
        return [parse(Int, m.captures[1])]
    end

    # Matrix keywords carry a world-axis row and pixel-axis column.
    m = match(r"^(?:PC|CD)([1-9][0-9]*)_([1-9][0-9]*)$", core)
    if m !== nothing
        return [parse(Int, m.captures[1]), parse(Int, m.captures[2])]
    end

    # Numeric projection parameters are attached to one WCS axis.
    m = match(r"^PV([1-9][0-9]*)_([0-9]|[1-9][0-9])$", core)
    if m !== nothing
        return [parse(Int, m.captures[1])]
    end

    # Paper III tabular / spectral parameters also carry one WCS axis.
    m = match(r"^PS([1-9][0-9]*)_([0-9]+)$", core)
    if m !== nothing
        return [parse(Int, m.captures[1])]
    end

    return Int[]
end

function inferred_wcs_axes(header::AbstractDict, alt::Char)
    alt_str = alt == ' ' ? "" : string(alt)
    max_index = 0

    # Use indexed WCS keywords when WCSAXES is omitted, as specified by Paper I.
    for key in keys(header)
        key isa AbstractString || continue
        for idx in _wcs_axis_indices(String(key), alt_str)
            max_index = max(max_index, idx)
        end
    end

    return max_index
end

function validate_wcs_keyword_dimensions(header::AbstractDict, naxis::Int, alt::Char)
    alt_str = alt == ' ' ? "" : string(alt)

    # Explicit dimensionality must bound every indexed keyword for this version.
    for key in keys(header)
        key isa AbstractString || continue
        for idx in _wcs_axis_indices(String(key), alt_str)
            if idx > naxis
                throw(ArgumentError(
                    "keyword $key uses axis $idx but WCSAXES/NAXIS is $naxis"
                ))
            end
        end
    end
end

function has_indexed_wcs_keyword(header::AbstractDict, prefix::AbstractString,
                                  naxis::Int, alt_str::AbstractString)
    # Detect whether a matrix family is present before reading its entries.
    for i in 1:naxis, j in 1:naxis
        haskey(header, "$(prefix)$(i)_$(j)$(alt_str)") && return true
    end
    return false
end

# ──────────────────────────────────────────────────────────────────────────────
# CD-matrix construction
# ──────────────────────────────────────────────────────────────────────────────

"""
    build_cd_matrix(header, naxis, cdelt, alt) -> Matrix{Float64}

Construct the naxis × naxis CD matrix from the header keywords.

Priority (Paper I, Appendix A):
1. `CDi_ja` – explicit CD matrix; CDELT/PC are ignored.
2. `PCi_ja` + `CDELTia` – linear transform matrix.
3. `CROTAi` – legacy rotation keyword; converted to equivalent CD.
4. Diagonal `CDELTia` – no rotation.

The returned matrix has entry `[i, j]` = `CD_ij` in the standard sense
(`x_i = Σ_j CD_ij (p_j − CRPIX_j)`).
"""
function build_cd_matrix(header::AbstractDict, naxis::Int,
                          cdelt::Vector{Float64}, alt::Char)
    alt_str = alt == ' ' ? "" : string(alt)
    has_cd_keywords = has_indexed_wcs_keyword(header, "CD", naxis, alt_str)
    has_pc_keywords = has_indexed_wcs_keyword(header, "PC", naxis, alt_str)

    if has_cd_keywords && has_pc_keywords
        throw(ArgumentError("CDi_ja and PCi_ja keywords must not be mixed"))
    end

    # ── Try CDi_j keywords ────────────────────────────────────────────────────
    cd = zeros(Float64, naxis, naxis)
    has_cd = false
    for i in 1:naxis, j in 1:naxis
        key = "CD$(i)_$(j)$(alt_str)"
        if haskey(header, key)
            cd[i, j] = Float64(header[key])
            has_cd = true
        end
    end
    has_cd && return cd

    # ── Try PCi_j keywords (multiply by CDELTi) ───────────────────────────────
    pc = Matrix{Float64}(I, naxis, naxis)  # default: identity
    has_pc = false
    for i in 1:naxis, j in 1:naxis
        key = "PC$(i)_$(j)$(alt_str)"
        if haskey(header, key)
            pc[i, j] = Float64(header[key])
            has_pc = true
        end
    end
    if has_pc
        # CD_ij = CDELT_i * PC_ij
        for i in 1:naxis, j in 1:naxis
            cd[i, j] = cdelt[i] * pc[i, j]
        end
        return cd
    end

    # ── Try CROTA2 legacy keyword ─────────────────────────────────────────────
    crota_key = "CROTA2$(alt_str)"
    if naxis >= 2 && haskey(header, crota_key)
        crota = Float64(header[crota_key])
        # Paper I, Appendix A, Eq. A2:
        # CD11 = CDELT1*cos(rho)   CD12 = -CDELT2*sin(rho)
        # CD21 = CDELT1*sin(rho)   CD22 =  CDELT2*cos(rho)
        rho = crota * (π / 180.0)
        cr, sr = cos(rho), sin(rho)
        # Fill in only 2D part for CROTA; off-diagonal higher axes stay 0.
        for k in 1:naxis
            cd[k, k] = cdelt[k]
        end
        cd[1, 1] = cdelt[1] * cr
        cd[1, 2] = -cdelt[2] * sr
        cd[2, 1] = cdelt[1] * sr
        cd[2, 2] = cdelt[2] * cr
        return cd
    end

    # ── Fall back: diagonal CDELT ─────────────────────────────────────────────
    for k in 1:naxis
        cd[k, k] = cdelt[k]
    end
    return cd
end

# ──────────────────────────────────────────────────────────────────────────────
# Main parsing entry point
# ──────────────────────────────────────────────────────────────────────────────

"""
    WCS(header; fobj=nothing, alt=' ', minerr=0.0) -> WCSTransform

Parse a FITS WCS header into a `WCSTransform`.

`header` is any `AbstractDict{<:AbstractString, Any}` containing FITS keyword
values.  FITSIO.jl `FITSHeader` objects can be converted to `Dict` with
`Dict(header)`, or use the FITSIO extension if available.

`fobj` optionally supplies the owning FITS container for external WCS arrays.
`alt` selects the alternate WCS description character (`' '` for the primary,
`'A'`–`'Z'` for alternates).
`minerr` is reserved for Paper IV distortion support: when auxiliary lookup
tables are implemented, distortion components whose declared error estimate is
below `minerr` may be skipped.  It is currently passed through to auxiliary-data
resolver methods but is not used by the core header-only parser.

Throws `ArgumentError` for headers that are clearly malformed (e.g., axis
count mismatch in a CD or PC matrix keyword).

## Example

```julia
hdr = Dict(
    "NAXIS"  => 2,
    "CTYPE1" => "RA---TAN",
    "CTYPE2" => "DEC--TAN",
    "CRPIX1" => 512.0,
    "CRPIX2" => 512.0,
    "CRVAL1" => 83.8221,
    "CRVAL2" => -5.3911,
    "CDELT1" => -2.7778e-4,
    "CDELT2" =>  2.7778e-4,
)
wcs = WCS(hdr)
```
"""
function WCS(header::AbstractDict; fobj = nothing, alt::Char = ' ', minerr::Real = 0.0)
    if alt != ' ' && !(('A' <= alt <= 'Z'))
        throw(ArgumentError("alt must be ' ' or 'A'–'Z', got $(repr(alt))"))
    end
    alt_str = alt == ' ' ? "" : string(alt)
    aux = _auxiliary_wcs_data(header, fobj; alt = alt, minerr = minerr)

    # ── Number of WCS axes ────────────────────────────────────────────────────
    # WCSAXES takes precedence over NAXIS (Paper I, Sec. 2.1).
    naxis_key = "WCSAXES$(alt_str)"
    naxis_key2 = "WCSAXES"
    naxis = if haskey(header, naxis_key)
        Int(header[naxis_key])
    elseif alt == ' ' && haskey(header, naxis_key2)
        Int(header[naxis_key2])
    elseif haskey(header, "NAXIS")
        max(Int(header["NAXIS"]), inferred_wcs_axes(header, alt))
    else
        throw(ArgumentError("header contains neither WCSAXES nor NAXIS"))
    end

    naxis >= 1 || throw(ArgumentError("NAXIS/WCSAXES must be ≥ 1, got $naxis"))
    validate_wcs_keyword_dimensions(header, naxis, alt)

    # ── SCAMP compatibility: pre-2012 -TAN + PV≥5 → -TPV ─────────────────────
    # Must run before CTYPE parsing so the mutated header values are read.
    _fix_scamp_compatibility!(header, alt_str)

    # ── Core per-axis keywords ────────────────────────────────────────────────
    crpix = Vector{Float64}(undef, naxis)
    crval = Vector{Float64}(undef, naxis)
    cdelt = Vector{Float64}(undef, naxis)
    ctype = Vector{String}(undef, naxis)
    cunit = Vector{String}(undef, naxis)

    for i in 1:naxis
        crpix[i] = Float64(get(header, "CRPIX$(i)$(alt_str)", 0.0))
        crval[i] = Float64(get(header, "CRVAL$(i)$(alt_str)", 0.0))
        cdelt[i] = Float64(get(header, "CDELT$(i)$(alt_str)", 1.0))
        ctype[i] = String(get(header, "CTYPE$(i)$(alt_str)", ""))
        cunit[i] = String(get(header, "CUNIT$(i)$(alt_str)", ""))
    end
    reject_unsupported_nonlinear_axes(ctype)

    # ── Parse CTYPE to determine projection and axis roles ────────────────────
    lon_axis = 0
    lat_axis = 0
    proj_code = ""

    for i in 1:naxis
        _, coord_type, pc = parse_ctype(ctype[i])
        pc == "TAB" && continue
        if coord_type == :lon
            if lon_axis != 0
                throw(ArgumentError("header has multiple longitude axes (CTYPE$(lon_axis) and CTYPE$(i))"))
            end
            lon_axis = i
            proj_code = pc
        elseif coord_type == :lat
            if lat_axis != 0
                throw(ArgumentError("header has multiple latitude axes (CTYPE$(lat_axis) and CTYPE$(i))"))
            end
            lat_axis = i
            # Consistency: both axes should have the same projection code.
            lat_pc = parse_ctype(ctype[i])[3]
            if !isempty(proj_code) && !isempty(lat_pc) && proj_code != lat_pc
                throw(ArgumentError(
                    "longitude and latitude CTYPEs have different projection codes: " *
                    "\"$proj_code\" vs \"$lat_pc\""
                ))
            end
            if isempty(proj_code)
                proj_code = lat_pc
            end
        end
    end

    # Validate paired celestial axes.
    if (lon_axis == 0) != (lat_axis == 0)
        throw(ArgumentError(
            "header has a celestial longitude axis without a matching latitude axis, " *
            "or vice versa"
        ))
    end

    projection = isempty(proj_code) ? nothing :
        projection_from_header(proj_code, header, lon_axis, lat_axis, alt)

    # ── Parse optional SIP distortion before building transform output ────────
    # TPV/TPD takes precedence over SIP (SCAMP convention).
    if proj_code in ("TPV", "TPD")
        _remove_sip_keywords!(header, alt_str)
    end
    sip = parse_sip_distortion(header, crpix, naxis, alt)
    pipeline = distortion_pipeline(sip, aux)

    # ── Build the CD matrix ───────────────────────────────────────────────────
    cd = build_cd_matrix(header, naxis, cdelt, alt)

    # ── Normalize celestial angular units to degrees at the public API boundary
    # Linear and spectral axes remain in the units encoded by their headers.
    for i in (lon_axis, lat_axis)
        i == 0 && continue
        f = unit_to_deg(cunit[i])
        isnan(f) && throw(ArgumentError("unsupported celestial CUNIT$(i): $(repr(cunit[i]))"))
        crval[i] *= f
        cd[i, :] .*= f
    end

    # ── Celestial pole parameters ─────────────────────────────────────────────
    lonpole_raw = get(header, "LONPOLE$(alt_str)", get(header, "LONPOLE", nothing))
    latpole_raw = get(header, "LATPOLE$(alt_str)", get(header, "LATPOLE", nothing))

    lonpole::Float64 = if !isnothing(lonpole_raw)
        Float64(lonpole_raw)
    elseif projection !== nothing
        theta0 = native_theta0(projection)
        delta0 = (lat_axis > 0) ? crval[lat_axis] : 0.0
        default_lonpole(delta0, theta0)
    else
        0.0
    end

    latpole::Float64 = !isnothing(latpole_raw) ? Float64(latpole_raw) : 90.0

    # ── Compute native pole position in celestial coords (alpha_p, delta_p) ──
    # These are the celestial coordinates of the native north pole (theta=90°).
    # For zenithal projections (theta0=90°), alpha_p=CRVAL_lon, delta_p=CRVAL_lat.
    # For other projections, compute from the constraint that the fiducial native
    # point (phi0, theta0) maps to celestial (alpha0, delta0) = CRVAL.
    # Paper II, Section 2.4, Eq. 11.
    alpha_p::Float64 = 0.0
    delta_p::Float64 = 90.0

    if projection !== nothing && lon_axis > 0 && lat_axis > 0
        alpha0 = crval[lon_axis] * (π / 180.0)
        delta0 = crval[lat_axis] * (π / 180.0)
        phi0 = native_phi0(projection) * (π / 180.0)
        theta0 = native_theta0(projection) * (π / 180.0)
        phi_p = lonpole * (π / 180.0)

        ap, dp = compute_native_pole(
            alpha0, delta0, phi0, theta0, phi_p,
            latpole * (π / 180.0)
        )
        alpha_p = ap * (180.0 / π)
        delta_p = dp * (180.0 / π)
    end

    return WCSTransform(
        naxis,
        SVector{naxis, Float64}(crpix),
        SVector{naxis, Float64}(crval),
        SMatrix{naxis, naxis, Float64, naxis * naxis}(cd),
        ctype, cunit,
        lonpole, latpole, alpha_p, delta_p,
        projection, pipeline, aux, lon_axis, lat_axis,
    )
end
