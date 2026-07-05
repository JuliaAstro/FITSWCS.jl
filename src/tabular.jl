"""
Tabular WCS data structures and header parsing.
"""

"""Abstract supertype for resolved FITS `-TAB` coordinate payloads."""
abstract type AbstractTabularWCSData end

"""No-op tabular payload for WCS transforms with no `-TAB` axes."""
struct NoTabularWCSData <: AbstractTabularWCSData end

"""
    TabularAxisSpec

Parse-time reference from one WCS `-TAB` axis to an external binary-table
coordinate array and optional index-vector column.
"""
struct TabularAxisSpec{A}
    axis::Int
    table_axis::Int
    extname::String
    extver::Int
    extlev::Int
    coord_column::String
    index_column::A
end

"""
    TabularWCSTable

Backend-independent resolved `-TAB` table used by coordinate transforms.

The type parameter `M` encodes the number of table axes, enabling the compiler to
stack-allocate intermediate buffers and unroll interpolation loops.
"""
struct TabularWCSTable{M, C, I}
    axes::SVector{M, Int}
    table_axes::SVector{M, Int}
    coord::C
    indices::I
end

"""Resolved collection of all `-TAB` tables associated with one WCS."""
struct TabularWCSData{T} <: AbstractTabularWCSData
    tables::T
end

has_tabular(::NoTabularWCSData) = false
has_tabular(::TabularWCSData) = true

function _tabular_axis_numbers(header::AbstractDict, alt_str::AbstractString)
    naxis = Int(get(header, "NAXIS", 0))
    axes = Int[]

    # Identify selected CTYPE axes whose algorithm code is TAB.
    suffix = uppercase(String(alt_str))
    for (key, value) in header
        key isa AbstractString || continue
        m = match(Regex("^CTYPE([1-9][0-9]*)$(suffix)\$"), uppercase(String(key)))
        m === nothing && continue
        ctype = uppercase(strip(String(value)))
        ncodeunits(ctype) >= 8 && ctype[5] == '-' && ctype[6:8] == "TAB" ||
            continue
        axis = parse(Int, m.captures[1])
        axis > naxis &&
            throw(ArgumentError("CTYPE$(axis) references axis $(axis), " *
                                "but NAXIS is $(naxis)"))
        push!(axes, axis)
    end

    return sort!(axes)
end

function _tabular_axis_spec(header::AbstractDict, axis::Int, alt_str::AbstractString)
    ext_key = "PS$(axis)_0$(alt_str)"
    coord_key = "PS$(axis)_1$(alt_str)"
    index_key = "PS$(axis)_2$(alt_str)"

    # WCSLIB requires the extension name and coordinate-array column.
    haskey(header, ext_key) ||
        throw(ArgumentError("TAB axis CTYPE$(axis) requires $ext_key"))
    haskey(header, coord_key) ||
        throw(ArgumentError("TAB axis CTYPE$(axis) requires $coord_key"))

    table_axis = Int(get(header, "PV$(axis)_3$(alt_str)", 1))
    table_axis >= 1 ||
        throw(ArgumentError("PV$(axis)_3$(alt_str) must be positive, got $table_axis"))
    extver = Int(get(header, "PV$(axis)_1$(alt_str)", 1))
    extlev = Int(get(header, "PV$(axis)_2$(alt_str)", 1))
    extname = String(header[ext_key])
    coord_column = String(header[coord_key])

    # Keep the optional index column concrete by constructing distinct spec types.
    if haskey(header, index_key)
        return TabularAxisSpec(axis, table_axis, extname, extver, extlev, coord_column, String(header[index_key]))
    end
    return TabularAxisSpec(axis, table_axis, extname, extver, extlev, coord_column, nothing)
end

function _tabular_axis_specs(header::AbstractDict; alt::Char = ' ')
    alt_str = alt == ' ' ? "" : string(alt)
    axes = _tabular_axis_numbers(header, alt_str)

    # Build parse-time specs only for axes selected by the active alternate WCS.
    return [_tabular_axis_spec(header, axis, alt_str) for axis in axes]
end

function _tabular_groups(specs::AbstractVector{<:TabularAxisSpec})
    groups = Dict{Tuple{String, Int, Int, String}, Vector{TabularAxisSpec}}()

    # Group axes that share one coordinate-array column.
    for spec in specs
        key = (spec.extname, spec.extver, spec.extlev, spec.coord_column)
        push!(get!(groups, key, TabularAxisSpec[]), spec)
    end

    return collect(values(groups))
end

function _validate_tabular_group(specs::AbstractVector{<:TabularAxisSpec})
    isempty(specs) && throw(ArgumentError("TAB table group must not be empty"))
    seen = Set{Int}()

    # Each WCS axis in a group must map to a distinct table dimension.
    for spec in specs
        spec.table_axis in seen &&
            throw(ArgumentError("duplicate TAB table axis $(spec.table_axis) for coordinate column $(spec.coord_column)"))
        push!(seen, spec.table_axis)
    end

    max_axis = maximum(seen)
    all(in(seen), 1:max_axis) ||
        throw(ArgumentError("TAB table axes must define every dimension from 1 through $max_axis"))
    return specs
end

function _tabular_coordinate_array(coord, M::Int)
    M >= 1 || throw(ArgumentError("TAB coordinate dimensionality must be positive"))

    # One-dimensional TAB accepts any vector-like coordinate payload.
    if M == 1
        data = Float64.(vec(coord))
        isempty(data) && throw(ArgumentError("TAB coordinate array must not be empty"))
        return reshape(data, 1, length(data))
    end

    nd = ndims(coord)
    nd >= 2 || throw(ArgumentError("TAB coordinate array for $M axes must have at least two dimensions"))

    # Normalize to coord[component, table_axis_1, table_axis_2, ...].
    if size(coord, 1) == M
        return Array{Float64}(coord)
    elseif size(coord, nd) == M
        order = (nd, (1:(nd - 1))...)
        return Array{Float64}(permutedims(coord, order))
    end

    throw(ArgumentError("TAB coordinate array must have a coordinate component dimension of length $M"))
end

function _tabular_index_vector(index, K::Int)
    K >= 1 || throw(ArgumentError("TAB coordinate axes must not be empty"))

    # Missing index vectors use FITS default indexing in the 1-relative TAB coordinate.
    if isnothing(index)
        return collect(range(1.0, step = 1.0, length = K))
    end

    values = Float64.(vec(index))
    length(values) == K ||
        throw(ArgumentError("TAB index vector has length $(length(values)), expected $K"))
    return values
end

function _tabular_table_from_group(specs::AbstractVector{<:TabularAxisSpec}, coord, raw_indices)
    ordered = sort(collect(_validate_tabular_group(specs)); by = spec -> spec.table_axis)
    M = length(ordered)
    data = _tabular_coordinate_array(coord, M)
    ndims(data) == M + 1 ||
        throw(ArgumentError("TAB coordinate array has $(ndims(data) - 1) table axes, expected $M"))

    # Convert optional index arrays into concrete vectors aligned with table axes.
    indices = Vector{Vector{Float64}}(undef, M)
    for (m, spec) in pairs(ordered)
        K = size(data, m + 1)
        raw = get(raw_indices, spec.axis, nothing)
        indices[m] = _tabular_index_vector(raw, K)
    end

    # Dispatch on M to construct the appropriately-typed TabularWCSTable{M}.
    # This is a runtime-to-compile-time bridge: each branch returns a concrete
    # TabularWCSTable{M} so that the hot path (_tabular_forward) can dispatch
    # on M and stack-allocate intermediate buffers.
    if M == 1
        return TabularWCSTable{1, typeof(data), typeof(indices)}(
            SVector{1, Int}(ordered[1].axis),
            SVector{1, Int}(ordered[1].table_axis),
            data, indices)
    elseif M == 2
        return TabularWCSTable{2, typeof(data), typeof(indices)}(
            SVector{2, Int}(ordered[1].axis, ordered[2].axis),
            SVector{2, Int}(ordered[1].table_axis, ordered[2].table_axis),
            data, indices)
    elseif M == 3
        return TabularWCSTable{3, typeof(data), typeof(indices)}(
            SVector{3, Int}(ordered[1].axis, ordered[2].axis, ordered[3].axis),
            SVector{3, Int}(ordered[1].table_axis, ordered[2].table_axis, ordered[3].table_axis),
            data, indices)
    elseif M == 4
        return TabularWCSTable{4, typeof(data), typeof(indices)}(
            SVector{4, Int}(ordered[1].axis, ordered[2].axis, ordered[3].axis, ordered[4].axis),
            SVector{4, Int}(ordered[1].table_axis, ordered[2].table_axis, ordered[3].table_axis, ordered[4].table_axis),
            data, indices)
    else
        # Fallback for high-dimensional TAB tables: construct with runtime M.
        # SVector{length(v)}(v) works with a runtime length, producing a
        # properly-typed SVector.  This path is type-unstable but only runs
        # during WCS construction, not in the hot transform path.
        return TabularWCSTable{M, typeof(data), typeof(indices)}(
            SVector{M, Int}(spec.axis for spec in ordered),
            SVector{M, Int}(spec.table_axis for spec in ordered),
            data, indices)
    end
end

function _tabular_auxiliary_data(header::AbstractDict, loader; alt::Char = ' ')
    specs = _tabular_axis_specs(header; alt = alt)
    isempty(specs) && return NoTabularWCSData()

    # Load each shared coordinate-array group through the backend callback.
    tables = map(_tabular_groups(specs)) do group
        ordered = sort(collect(group); by = spec -> spec.table_axis)
        first_spec = first(ordered)
        coord = loader(first_spec.extname, first_spec.extver, first_spec.extlev, first_spec.coord_column)
        raw_indices = Dict{Int, Any}()
        for spec in ordered
            isnothing(spec.index_column) && continue
            raw_indices[spec.axis] = loader(spec.extname, spec.extver, spec.extlev, spec.index_column)
        end
        _tabular_table_from_group(ordered, coord, raw_indices)
    end

    return TabularWCSData(tuple(tables...))
end

function _tabular_fractional_index(index::AbstractVector{<:Real}, value::Real)
    T = promote_type(_promote_float_type(value), _float_type(eltype(index)))
    K = length(index)
    K >= 1 || throw(ArgumentError("TAB index vector must not be empty"))
    K == 1 && return one(T)

    first_value = T(index[1])
    last_value = T(index[end])
    increasing = last_value >= first_value

    # Locate the bracketing index-vector cell.
    for k in 1:(K - 1)
        a = T(index[k])
        b = T(index[k + 1])
        lo = min(a, b)
        hi = max(a, b)
        if lo <= T(value) <= hi
            iszero(b - a) && return T(k)
            return T(k) + (T(value) - a) / (b - a)
        end
    end

    # Match WCSLIB's half-cell tolerance at the edges.
    if increasing
        before = first_value - (T(index[2]) - first_value) / 2
        after = last_value + (last_value - T(index[end - 1])) / 2
        before <= T(value) < first_value && return one(T)
        last_value < T(value) <= after && return T(K)
    else
        before = first_value + (first_value - T(index[2])) / 2
        after = last_value - (T(index[end - 1]) - last_value) / 2
        first_value < T(value) <= before && return one(T)
        after <= T(value) < last_value && return T(K)
    end

    throw(ArgumentError("TAB coordinate $value is outside index-vector bounds"))
end

function _tabular_index_value(index::AbstractVector{<:Real}, upsilon::Real)
    T = promote_type(_promote_float_type(upsilon), _float_type(eltype(index)))
    K = length(index)
    K >= 1 || throw(ArgumentError("TAB index vector must not be empty"))
    K == 1 && return T(index[1])

    clamped = clamp(T(upsilon), one(T), T(K))
    lo = min(floor(Int, clamped), K - 1)
    hi = lo + 1
    weight = clamped - T(lo)
    return (one(T) - weight) * T(index[lo]) + weight * T(index[hi])
end

function _tabular_cell(upsilon::Real, K::Int)
    T = _promote_float_type(upsilon)
    K >= 1 || throw(ArgumentError("TAB coordinate axes must not be empty"))
    K == 1 && return 1, 1, zero(T)

    clamped = clamp(T(upsilon), one(T), T(K))
    lo = min(floor(Int, clamped), K - 1)
    hi = lo + 1
    weight = clamped - T(lo)
    return lo, hi, weight
end

"""
    _tabular_forward(table::TabularWCSTable, intermediate::AbstractVector, crval::AbstractVector)

Evaluate one `-TAB` coordinate array at the given intermediate world coordinates.

For each axis controlled by `table`, the coordinate lookup key is
``psi_m = intermediate[axis_m] + crval[axis_m]``, i.e. the intermediate coordinate
with the WCS reference value folded back in.  `psi_m` is converted to a fractional
index position via `table.indices[m]` (or 1-based default indexing when no index
vector is supplied), then multilinear interpolation is performed over the enclosing
``2^M``-vertex cell of `table.coord`.

Returns a tuple of `M` world-coordinate values, one per table axis, in the same
order as `table.axes`.
"""
function _tabular_forward(table::TabularWCSTable{M}, intermediate::AbstractVector, crval::AbstractVector) where {M}
    cells = ntuple(m -> begin
        axis = table.axes[m]
        psi = intermediate[axis] + crval[axis]
        upsilon = _tabular_fractional_index(table.indices[m], psi)
        _tabular_cell(upsilon, size(table.coord, m + 1))
    end, Val(M))

    # Multilinearly interpolate each coordinate component over the enclosing cell.
    return ntuple(component -> begin
        total = 0.0
        for corner in 0:(2^M - 1)
            weight = 1.0
            idx = MVector{M, Int}(undef)
            for m in 1:M
                lo, hi, t = cells[m]
                if iszero((corner >> (m - 1)) & 1)
                    idx[m] = lo
                    weight *= 1 - t
                else
                    idx[m] = hi
                    weight *= t
                end
            end
            total += weight * table.coord[component, Tuple(idx)...]
        end
        total
    end, Val(M))
end

function _monotonic_sense(values::AbstractVector{<:Real})
    length(values) <= 1 && return 1
    all(values[i] <= values[i + 1] for i in 1:(length(values) - 1)) && return 1
    all(values[i] >= values[i + 1] for i in 1:(length(values) - 1)) && return -1
    return 0
end

function _tabular_inverse_1d(table::TabularWCSTable{M}, world_value::Real, crval::Real) where {M}
    M == 1 ||
        throw(ArgumentError("coupled multidimensional TAB inverse is not implemented yet"))
    coord_values = vec(view(table.coord, 1, :))
    _monotonic_sense(coord_values) != 0 ||
        throw(ArgumentError("1D TAB inverse requires monotonic coordinate values"))

    upsilon = _tabular_fractional_index(coord_values, world_value)
    psi = _tabular_index_value(table.indices[1], upsilon)
    return psi - crval
end

# ── M=2 bilinear cell-scan ──────────────────────────────────────────────────

"""
    _tabular_inverse_bilinear(table::TabularWCSTable{2}, world, crval)

Invert a two-dimensional `-TAB` lookup by scanning coordinate-array cells
for the one whose bilinear map contains the target world coordinate, then
solving the bilinear system in closed form.

The bilinear map within cell ``(k_1, k_2)`` is
``f(u,v) = a + b·u + c·v + d·u·v`` where ``u,v ∈ [0,1]``.
Substituting ``u = (r₁ - c₁·v) / (b₁ + d₁·v)`` into the second component
yields a quadratic in ``v``, solved by the quadratic formula.  Two degenerate
cases are handled explicitly: purely linear cells (``d₁ = d₂ = 0``) and
cells where one cross-term vanishes.
"""
function _tabular_inverse_bilinear(table::TabularWCSTable{2}, world::AbstractVector, crval::AbstractVector)
    T = _coordinate_float_type(world)
    axis1 = table.axes[1]
    axis2 = table.axes[2]
    t1 = T(world[axis1])
    t2 = T(world[axis2])
    K1 = size(table.coord, 2)
    K2 = size(table.coord, 3)
    tol = _convergence_tol(T)
    idx1 = table.indices[1]
    idx2 = table.indices[2]

    for k1 in 1:K1-1, k2 in 1:K2-1
        a1 = table.coord[1, k1, k2]
        a2 = table.coord[2, k1, k2]
        b1 = table.coord[1, k1+1, k2] - a1
        b2 = table.coord[2, k1+1, k2] - a2
        c1 = table.coord[1, k1, k2+1] - a1
        c2 = table.coord[2, k1, k2+1] - a2
        d1 = table.coord[1, k1+1, k2+1] - a1 - b1 - c1
        d2 = table.coord[2, k1+1, k2+1] - a2 - b2 - c2
        r1 = t1 - a1
        r2 = t2 - a2

        uv = _solve_bilinear_cell(b1, b2, c1, c2, d1, d2, r1, r2, tol)
        uv === nothing && continue
        u, v = uv

        psi1 = _tabular_index_value(idx1, k1 + u)
        psi2 = _tabular_index_value(idx2, k2 + v)
        return SVector{2, T}(
            psi1 - T(crval[axis1]),
            psi2 - T(crval[axis2]))
    end

    throw(ArgumentError("TAB inverse: target not found in any cell"))
end

"""
    _solve_bilinear_cell(b1, b2, c1, c2, d1, d2, r1, r2, tol)

Solve ``r₁ = b₁·u + c₁·v + d₁·u·v``, ``r₂ = b₂·u + c₂·v + d₂·u·v``
for ``(u,v)`` in closed form.

`tol` is the threshold below which cross-term coefficients ``d₁, d₂`` are
treated as zero, routing to the purely-linear 2×2 solve.  It also guards
against near-singular denominators (``b₁ + d₁·v ≈ 0``).
We recommend using `_convergence_tol` to set this.

Returns `(u, v)` if a solution exists within the cell (with minor
extrapolation tolerance ±0.5), or `nothing` otherwise.
"""
function _solve_bilinear_cell(b1, b2, c1, c2, d1, d2, r1, r2, tol)
    # Purely linear cell: no cross-terms.
    if abs(d1) < tol && abs(d2) < tol
        det = b1*c2 - b2*c1
        abs(det) < tol && return nothing
        u = (c2*r1 - c1*r2) / det
        v = (b1*r2 - b2*r1) / det
        return (-0.5 <= u <= 1.5 && -0.5 <= v <= 1.5) ? (u, v) : nothing
    end

    # d₁ ≠ 0: solve  u = (r₁ − c₁·v) / (b₁ + d₁·v)  → quadratic in v.
    #   A·v² + B·v + C = 0
    #   A = c₂·d₁ − d₂·c₁
    #   B = c₂·b₁ + d₂·r₁ − r₂·d₁ − b₂·c₁
    #   C = b₂·r₁ − r₂·b₁
    if abs(d1) >= tol
        A = c2*d1 - d2*c1
        B = c2*b1 + d2*r1 - r2*d1 - b2*c1
        C = b2*r1 - r2*b1

        if abs(A) < tol
            # Degenerate: quadratic → linear.
            abs(B) < tol && return nothing
            v = -C / B
            denom = b1 + d1*v
            abs(denom) < tol && return nothing
            u = (r1 - c1*v) / denom
            return (-0.5 <= u <= 1.5 && -0.5 <= v <= 1.5) ? (u, v) : nothing
        end

        disc = B*B - 4*A*C
        disc < 0 && return nothing
        sqrt_disc = sqrt(disc)

        for v in ((-B + sqrt_disc) / (2*A), (-B - sqrt_disc) / (2*A))
            denom = b1 + d1*v
            if abs(denom) > tol
                u = (r1 - c1*v) / denom
                if -0.5 <= u <= 1.5 && -0.5 <= v <= 1.5
                    return (u, v)
                end
            end
        end
        return nothing
    end

    # d₁ = 0, d₂ ≠ 0: solve  v = (r₂ − b₂·u) / (c₂ + d₂·u)  → quadratic in u.
    # Symmetric to above with (b ↔ c, 1 ↔ 2):
    A = b1*d2
    B = b1*c2 - r1*d2 - c1*b2
    C = c1*r2 - r1*c2

    if abs(A) < tol
        abs(B) < tol && return nothing
        u = -C / B
        denom = c2 + d2*u
        abs(denom) < tol && return nothing
        v = (r2 - b2*u) / denom
        return (-0.5 <= u <= 1.5 && -0.5 <= v <= 1.5) ? (u, v) : nothing
    end

    disc = B*B - 4*A*C
    disc < 0 && return nothing
    sqrt_disc = sqrt(disc)

    for u in ((-B + sqrt_disc) / (2*A), (-B - sqrt_disc) / (2*A))
        denom = c2 + d2*u
        if abs(denom) > tol
            v = (r2 - b2*u) / denom
            if -0.5 <= u <= 1.5 && -0.5 <= v <= 1.5
                return (u, v)
            end
        end
    end
    return nothing
end

# ── General-M Newton solver ─────────────────────────────────────────────────

"""
    _tabular_inverse_newton(table::TabularWCSTable{M}, world, crval)

Invert a multi-dimensional (M ≥ 3) `-TAB` lookup via Newton-Raphson iteration
with a finite-difference Jacobian.

A starting point is obtained by scanning the coordinate-array vertices for the
one nearest to the target world coordinate in Euclidean distance.  The Newton
iteration then refines this estimate using the piecewise-multilinear forward
map.  All work arrays are stack-allocated via `MVector` / `MMatrix`.
"""
function _tabular_inverse_newton(table::TabularWCSTable{M}, world::AbstractVector, crval::AbstractVector) where {M}
    T = _coordinate_float_type(world)
    target = SVector{M, T}(T(world[axis]) for axis in table.axes)

    # Nearest-vertex initial guess.
    vars = MVector{M, T}(undef)
    best_r = typemax(T)
    for index in CartesianIndices(size(table.coord)[2:end])
        r = zero(T)
        for component in 1:M
            residual = T(table.coord[component, Tuple(index)...] - target[component])
            r += residual^2
        end
        if r < best_r
            best_r = r
            for m in 1:M
                vars[m] = T(table.indices[m][index[m]] - crval[table.axes[m]])
            end
        end
    end

    max_axis = maximum(table.axes)
    trial = zeros(T, max_axis)
    tol = _convergence_tol(T)
    h_base = sqrt(eps(T))
    residual = MVector{M, T}(undef)
    step = MVector{M, T}(undef)

    for _ in 1:32
        # Evaluate forward at current estimate.
        @inbounds for m in 1:M
            trial[table.axes[m]] = vars[m]
        end
        values = _tabular_forward(table, trial, crval)
        @inbounds for m in 1:M
            residual[m] = T(values[m]) - target[m]
        end
        sum(abs2, residual) <= tol^2 && return vars

        # Finite-difference Jacobian.
        jac = MMatrix{M, M, T}(undef)
        @inbounds for m in 1:M
            h = h_base * max(abs(vars[m]), one(T))
            psi = vars[m] + T(crval[table.axes[m]])
            if psi >= maximum(table.indices[m])
                h = -h
            elseif psi <= minimum(table.indices[m])
                h = abs(h)
            end
            trial[table.axes[m]] = vars[m] + h
            shifted = _tabular_forward(table, trial, crval)
            for row in 1:M
                jac[row, m] = (T(shifted[row]) - T(values[row])) / h
            end
            trial[table.axes[m]] = vars[m]
        end

        step_vec = try
            jac \ residual
        catch
            break
        end
        @inbounds for m in 1:M
            step[m] = step_vec[m]
        end
        vars .-= step
        sum(abs2, step) <= tol^2 && return vars
    end

    @warn "coupled TAB inverse failed to converge; returning best estimate"
    return SVector{M, T}(vars)
end

"""
    _tabular_inverse(table::TabularWCSTable, world::AbstractVector, crval::AbstractVector)

Invert the `-TAB` coordinate lookup to recover intermediate world coordinates.

For a one-dimensional table this searches the coordinate array for the bracket
containing each world value, then inverts through the index vector to recover
``psi_m``; the returned value is ``psi_m - crval[axis_m]``, i.e. the intermediate
coordinate without the reference offset.

For two-dimensional tables a bilinear cell scan solves the inverse in closed
form via the quadratic formula.

For higher-dimensional tables a Newton-Raphson iteration with a finite-difference
Jacobian is used, initialised from the nearest coordinate-array vertex.

Returns a vector of intermediate coordinate values, one per table axis, in the
same order as `table.axes`.
"""
function _tabular_inverse(table::TabularWCSTable{M}, world::AbstractVector, crval::AbstractVector) where {M}
    T = _coordinate_float_type(world)
    if M == 1
        axis = only(table.axes)
        return SVector{1, T}(_tabular_inverse_1d(table, world[axis], crval[axis]))
    elseif M == 2
        return _tabular_inverse_bilinear(table, world, crval)
    end
    return _tabular_inverse_newton(table, world, crval)
end
