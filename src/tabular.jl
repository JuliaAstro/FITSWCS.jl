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

function _tabular_initial_guess(table::TabularWCSTable{M}, target::AbstractVector, crval::AbstractVector) where {M}
    table_shape = size(table.coord)[2:end]
    best = fill(0.0, M)
    best_r = Inf

    # Use the nearest coordinate-array vertex as a robust starting point.
    for index in CartesianIndices(table_shape)
        r = 0.0
        for component in 1:M
            residual = table.coord[component, Tuple(index)...] - target[component]
            r += residual^2
        end
        if r < best_r
            best_r = r
            for m in 1:M
                best[m] = table.indices[m][index[m]] - crval[table.axes[m]]
            end
        end
    end

    return best
end

function _tabular_inverse_coupled(table::TabularWCSTable{M}, world::AbstractVector, crval::AbstractVector) where {M}
    target = [Float64(world[axis]) for axis in table.axes]
    vars = _tabular_initial_guess(table, target, crval)
    trial = zeros(Float64, maximum(table.axes))
    max_iter = 32
    tol = _convergence_tol(Float64)

    for _ in 1:max_iter
        # Evaluate the residual at the current intermediate-coordinate estimate.
        for (m, axis) in pairs(table.axes)
            trial[axis] = vars[m]
        end
        values = collect(_tabular_forward(table, trial, crval))
        residual = values .- target
        sqrt(sum(abs2, residual)) <= tol && return vars

        # Build a finite-difference Jacobian for the coupled TAB forward map.
        jac = MMatrix{M, M, Float64}(undef)
        for m in 1:M
            h = sqrt(eps(Float64)) * max(abs(vars[m]), 1.0)
            psi = vars[m] + crval[table.axes[m]]
            if psi >= maximum(table.indices[m])
                h = -h
            elseif psi <= minimum(table.indices[m])
                h = abs(h)
            end
            trial[table.axes[m]] = vars[m] + h
            shifted = collect(_tabular_forward(table, trial, crval))
            for row in 1:M
                jac[row, m] = (shifted[row] - values[row]) / h
            end
            trial[table.axes[m]] = vars[m]
        end

        step = try
            jac \ residual
        catch
            break
        end
        vars .-= step
        sqrt(sum(abs2, step)) <= tol && return vars
    end

    @warn "coupled TAB inverse failed to converge; returning best estimate"
    return vars
end

"""
    _tabular_inverse(table::TabularWCSTable, world::AbstractVector, crval::AbstractVector)

Invert the `-TAB` coordinate lookup to recover intermediate world coordinates.

For a one-dimensional table this searches the coordinate array for the bracket
containing each world value, then inverts through the index vector to recover
``psi_m``; the returned value is ``psi_m - crval[axis_m]``, i.e. the intermediate
coordinate without the reference offset.

For coupled multi-dimensional tables a Newton-Raphson iteration with a
finite-difference Jacobian is used, initialised from the nearest coordinate-array
vertex.  Convergence is tested against `_convergence_tol(Float64)` with a maximum
of 32 iterations; if the iteration stalls a warning is emitted and the best
estimate seen so far is returned.

Returns a vector of intermediate coordinate values, one per table axis, in the
same order as `table.axes`.
"""
function _tabular_inverse(table::TabularWCSTable{M}, world::AbstractVector, crval::AbstractVector) where {M}
    if M == 1
        axis = only(table.axes)
        return [_tabular_inverse_1d(table, world[axis], crval[axis])]
    end
    return _tabular_inverse_coupled(table, world, crval)
end
