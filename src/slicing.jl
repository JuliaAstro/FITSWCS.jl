"""
WCS slicing wrapper.

Provides `slice_wcs` for producing sub-image WCS transforms and the
`SlicedWCSTransform` type that implements the sliced coordinate mapping.
"""

# ── Slice descriptor types ──────────────────────────────────────────────────────

"""
    KeepAll

Sentinel indicating a pixel axis is kept unchanged (no offset, no stride).
Equivalent to `:` in array-slicing notation.
"""
struct KeepAll end
_is_kept(::KeepAll) = true

"""
    KeepRange{R}(range::AbstractRange)

A kept pixel axis trimmed to `range`.  `R` is the concrete range type
(`UnitRange`, `StepRange`, etc.).  The slice maps pixel 1 in the sliced
image to `first(range)` in the original image, with `step(range)` stride.
"""
struct KeepRange{R <: AbstractRange}
    range::R
end
_is_kept(::KeepRange) = true

"""
    DropAxis{Float64}(pixel::Float64)

A dropped pixel axis fixed at `pixel` in the original image. `pixel` may be fractional.
Unlike array slicing, freezing a WCS axis at a fractional pixel position
is well defined (the transform is continuous in pixel space).
"""
struct DropAxis
    pixel::Float64
end
_is_kept(::DropAxis) = false

"""Normalize a slice argument (Colon, Real, or AbstractRange) into `KeepAll`, `KeepRange`, or `DropAxis`."""
_normalize_slice(::Colon) = KeepAll()
_normalize_slice(s::Real) = DropAxis(Float64(s))
_normalize_slice(s::AbstractRange) = KeepRange(s)
_normalize_slice(s) = throw(ArgumentError(
    "slice must be a Real, AbstractRange, or Colon, got $(typeof(s))"
))

# ── Slice composition (for recursive slicing of SlicedWCSTransform) ─────────────

"""
    _compose_slices(old, new) -> slice_spec

Combine an existing slice spec `old` (KeepAll, KeepRange, or DropAxis) with a
new raw slice argument `new` (Colon, Integer, or AbstractRange) on the
sub-image coordinate space.  Returns a new slice spec on the parent coordinate
space.

Used by `slice_wcs(::SlicedWCSTransform, ...)` to compose nested slices
without wrapper nesting.
"""
_compose_slices(::KeepAll, ::Colon) = KeepAll()
_compose_slices(::KeepAll, k::Real) = DropAxis(Float64(k))
_compose_slices(::KeepAll, r::AbstractRange) = KeepRange(r)

# old=KeepRange: sub-coordinate p_sub maps to parent p = a1 + s1·(p_sub - 1).
_compose_slices(old::KeepRange, ::Colon) = old
function _compose_slices(old::KeepRange, k::Real)
    a1 = first(old.range)
    s1 = step(old.range)
    return DropAxis(Float64(a1 + s1 * (k - 1)))
end
function _compose_slices(old::KeepRange, new::AbstractRange)
    a1 = first(old.range); s1 = step(old.range)
    a2 = first(new); s2 = step(new)
    a3 = a1 + s1 * (a2 - 1)
    s3 = s1 * s2
    return KeepRange(a3:s3:(a3 + s3 * (length(new) - 1)))
end

# old=DropAxis: unreachable — dropped axes are not in pixel_keep.

# ── Sliced WCS transform ────────────────────────────────────────────────────────

"""
    SlicedWCSTransform{N, Np, Nw, S, T} <: AbstractWCSTransform

WCS transform representing a sliced view of a parent `WCSTransform`.
Constructed with `slice_wcs`.

# Type parameters

- `N`  – original number of pixel/world axes in the parent.
- `Np` – number of kept pixel axes.
- `Nw` – number of kept world axes (may differ from `Np` when axes are coupled).
- `S`  – tuple type of the N slice descriptors (each `KeepRange` or `DropAxis`).
- `T`  – concrete type of the parent `WCSTransform`.

# Fields

- `parent` – the original full `WCSTransform`.
- `slices` – length-N tuple of `KeepRange` or `DropAxis`, one per parent axis.
- `pixel_keep` – indices of kept pixel axes (1-based, in parent axis order).
- `world_keep` – indices of kept world axes (1-based, in parent axis order).
- `dropped_world_values` – precomputed world-coordinate values for dropped axes.
  These are exact because a world axis is only dropped if it does not depend on
  any kept pixel axis (see `axis_correlation_matrix`).
"""
struct SlicedWCSTransform{N, Np, Nw, S <: Tuple, T <: WCSTransform} <: AbstractWCSTransform
    parent::T
    slices::S
    pixel_keep::SVector{Np, Int}
    world_keep::SVector{Nw, Int}
    dropped_world_values::SVector{N, Float64}
end

# ── Constructor ─────────────────────────────────────────────────────────────────

"""
    _slice_normalized(wcs::WCSTransform, normalized)
Core construction logic shared by `slice_wcs` on both concrete types."""
function _slice_normalized(wcs::WCSTransform{N}, normalized) where {N}
    # Determine which pixel axes are kept.
    pixel_keep = Int[i for (i, s) in enumerate(normalized) if _is_kept(s)]
    Np = length(pixel_keep)
    Np > 0 || throw(ArgumentError("All pixel axes were dropped; a WCS must have at least one pixel axis"))

    # Determine which world axes to keep using the correlation matrix.
    corr = axis_correlation_matrix(wcs)
    world_keep = Int[w for w in 1:N if any(corr[w, pixel_keep])]
    Nw = length(world_keep)
    Nw > 0 || throw(ArgumentError("No world axes depend on the kept pixel axes"))

    # Precompute world-coordinate values for dropped world axes.
    # We evaluate the forward transform at a nominal pixel where kept axes are
    # at the reference pixel (crpix) and dropped axes are at their fixed values.
    # Because a world axis is only dropped if it does not depend on any kept
    # pixel axis, these values are exact for all kept-pixel positions.
    nominal_pixel = MVector{N, Float64}(undef)
    @inbounds for i in 1:N
        s = normalized[i]
        if s isa DropAxis
            nominal_pixel[i] = s.pixel
        else
            nominal_pixel[i] = wcs.crpix[i]
        end
    end
    dropped_world_values = pixel_to_world(wcs, SVector{N, Float64}(nominal_pixel))

    # Function barrier: Val(Np) and Val(Nw) let the compiler specialize
    # SVector construction for the exact sizes, avoiding type instability.
    return _construct_sliced(wcs, normalized, pixel_keep, world_keep,
                             dropped_world_values, Val(Np), Val(Nw))
end

function _construct_sliced(wcs::WCSTransform{N}, normalized::S,
                           pixel_keep::Vector{Int},
                           world_keep::Vector{Int},
                           dropped_world_values::AbstractVector,
                           ::Val{Np}, ::Val{Nw}) where {N, Np, Nw, S}
    pk = SVector{Np, Int}(pixel_keep)
    wk = SVector{Nw, Int}(world_keep)
    dwv = SVector{N, Float64}(ntuple(i -> dropped_world_values[i], N))
    return SlicedWCSTransform{N, Np, Nw, S, typeof(wcs)}(
        wcs, normalized, pk, wk, dwv
    )
end

# ── Public constructors ──────────────────────────────────────────────────────────

"""
    slice_wcs(wcs::WCSTransform, slices...) -> SlicedWCSTransform

Slice a `WCSTransform` along its pixel axes.

Each positional argument corresponds to one pixel axis in FITS/WCSTransform
axis order (axis 1, axis 2, ..., axis N) and is one of:

- `a:b` (`AbstractUnitRange`): keep the axis, trimmed to `[a, b]`.
  New CRPIX = CRPIX - (a - 1).
- `a:s:b` (`StepRange`): keep the axis with stride `s`.
  CRPIX and the CD matrix column rescale by `s`.
- `k` (`Real`): drop the axis, fixing it at pixel `k`. Fractional positions
  are allowed.

Returns a `SlicedWCSTransform` with dimensionality equal to the number of
range arguments.

Pixel ``1`` in the sliced image corresponds to the first element of each
range argument.  For example, ``pixel_to_world(swcs, [1, 1])`` returns
the same world coordinates as ``pixel_to_world(wcs, [a, b])`` when
slicing with ``a:b`` ranges.

# Examples

```jldoctest
julia> hdr = Dict("NAXIS" => 2, "CTYPE1" => "X", "CTYPE2" => "Y",
                  "CRPIX1" => 512.0, "CRPIX2" => 512.0,
                  "CRVAL1" => 0.0, "CRVAL2" => 0.0,
                  "CDELT1" => 1.0, "CDELT2" => 1.0);

julia> wcs = WCS(hdr);

julia> swcs = slice_wcs(wcs, 400:600, 500:600);

julia> pixel_to_world(swcs, [1, 1]) == pixel_to_world(wcs, [400, 500])
true
```
"""
function slice_wcs(wcs::WCSTransform{N}, slices::Vararg{Any, M}) where {N, M}
    # Pad missing trailing axes with Colon (keep all).
    if M < N
        padded = ntuple(i -> i <= M ? slices[i] : Colon(), N)
    else
        M == N || throw(ArgumentError(
            "too many slice arguments: got $M, expected at most $N"))
        padded = slices
    end
    # Normalize each slice to KeepAll, KeepRange, or DropAxis.
    normalized = _normalize_slice.(padded)
    return _slice_normalized(wcs, normalized)
end

"""
    slice_wcs(swcs::SlicedWCSTransform, slices...) -> SlicedWCSTransform

Slice an already-sliced WCS.  Maps the new slice arguments (applied to the
sub-image pixel axes) back to the original parent axes, composes them with the
existing slice specs, and constructs a new `SlicedWCSTransform` that wraps the
original parent — never nests wrappers.
"""
function slice_wcs(swcs::SlicedWCSTransform{N, Np}, slices::Vararg{Any, M}) where {N, Np, M}
    # Pad/validate: M new slices for Np kept pixel axes.
    if M < Np
        new_slices = ntuple(i -> i <= M ? slices[i] : Colon(), Np)
    elseif M == Np
        new_slices = slices
    else
        throw(ArgumentError(
            "too many slice arguments: got $M, expected at most $Np"))
    end

    # Map new slices to parent axes and compose with existing specs.
    old_slices = swcs.slices
    pk = swcs.pixel_keep
    composed = ntuple(N) do i
        k = something(findfirst(==(i), pk), 0)
        if k == 0
            old_slices[i]  # already dropped, pass through unchanged
        else
            _compose_slices(old_slices[i], new_slices[k])
        end
    end

    # Construct directly from the composed (already normalized) specs.
    return _slice_normalized(swcs.parent, composed)
end

"""
    slice_wcs(dict::Dict{Char, <:AbstractWCSTransform}, slices...) -> Dict{Char, SlicedWCSTransform}

Slice every WCS in a dictionary of the type returned by `WCS_all`.
"""
function slice_wcs(dict::Dict{Char, T}, slices...) where {T <: AbstractWCSTransform}
    return Dict(c => slice_wcs(wcs, slices...) for (c, wcs) in dict)
end

# ── Disambiguation: batch vector-of-vectors vs single-coordinate AbstractVector ─

# A Vector{<:AbstractVector} matches both the batch method on
# AbstractWCSTransform and the single-coordinate method on SlicedWCSTransform.
# Resolve by forwarding to the generic batch implementation.
function pixel_to_world(swcs::SlicedWCSTransform,
                        pixels::AbstractVector{<:AbstractVector})
    return invoke(pixel_to_world,
                  Tuple{AbstractWCSTransform, typeof(pixels)}, swcs, pixels)
end
function world_to_pixel(swcs::SlicedWCSTransform,
                        worlds::AbstractVector{<:AbstractVector})
    return invoke(world_to_pixel,
                  Tuple{AbstractWCSTransform, typeof(worlds)}, swcs, worlds)
end

# ── pixel_to_world ──────────────────────────────────────────────────────────────

function pixel_to_world(swcs::SlicedWCSTransform, pixel_sub::AbstractVector)
    length(pixel_sub) == pixel_n_dim(swcs) ||
        throw(DimensionMismatch("pixel_sub has length $(length(pixel_sub)), expected $(pixel_n_dim(swcs))"))

    # Expand the sliced pixel coordinate to full N-dim pixel space.
    full_pixel = _expand_pixel(swcs, pixel_sub)

    # Transform using the parent WCS (handles all distortions automatically).
    full_world = pixel_to_world(swcs.parent, full_pixel)

    # Extract world coordinates for kept world axes only.
    return _extract_kept(full_world, swcs.world_keep)
end

"""
Expand a sliced pixel coordinate to the full N-dimensional pixel space of the
parent WCS.

Uses `pixel_keep` to determine the mapping: for each parent axis `i`,
`findfirst(==(i), pixel_keep)` gives the index into `pixel_sub` (or
`nothing` if the axis was dropped).
"""
function _expand_pixel(swcs::SlicedWCSTransform{N}, pixel_sub::AbstractVector) where {N}
    T = _coordinate_float_type(pixel_sub)
    slices = swcs.slices
    pk = swcs.pixel_keep
    return SVector{N, T}(ntuple(Val(N)) do i
        s = slices[i]
        k = something(findfirst(==(i), pk), 0)
        if k == 0  # dropped axis
            T(s.pixel)
        elseif s isa KeepAll
            T(pixel_sub[k])
        else  # KeepRange
            a = T(first(s.range))
            stp = T(step(s.range))
            T(a) + T(stp) * (T(pixel_sub[k]) - one(T))
        end
    end)
end

# ── world_to_pixel ──────────────────────────────────────────────────────────────

function world_to_pixel(swcs::SlicedWCSTransform, world_sub::AbstractVector)
    length(world_sub) == world_n_dim(swcs) ||
        throw(DimensionMismatch("world_sub has length $(length(world_sub)), expected $(world_n_dim(swcs))"))

    # Build full N-dim world vector: user values for kept axes, precomputed
    # values for dropped axes.
    full_world = _build_full_world(swcs, world_sub)

    # Invert through parent WCS (single shot, no iteration).
    full_pixel = world_to_pixel(swcs.parent, full_world)

    # Extract kept pixel axes and undo the slice offset.
    pixel_sub = _extract_kept_pixel(swcs, full_pixel)

    return pixel_sub
end

"""
Build a full N-dimensional world vector from a sliced world vector.

Kept world axes get the user-supplied values; dropped world axes use the
precomputed `dropped_world_values` (exact because they don't depend on
kept pixel axes).
"""
function _build_full_world(swcs::SlicedWCSTransform{N}, world_sub::AbstractVector) where {N}
    T = _coordinate_float_type(world_sub)
    wk = swcs.world_keep
    dwv = swcs.dropped_world_values
    return SVector{N, T}(ntuple(Val(N)) do i
        k = something(findfirst(==(i), wk), 0)
        k == 0 ? T(dwv[i]) : T(world_sub[k])
    end)
end

"""
Extract kept pixel coordinates from a full pixel vector and undo slice offsets.

For each kept pixel axis: pixel_sub = (pixel_old - a) / s + 1
where a and s come from the slice range.
"""
function _extract_kept_pixel(swcs::SlicedWCSTransform{N, Np}, full_pixel::AbstractVector) where {N, Np}
    T = _coordinate_float_type(full_pixel)
    slices = swcs.slices

    # Build the converted N-element vector in one pass, then index with
    # pixel_keep (integer indexing is ~0 ns vs ~27 ns for boolean mask).
    all_pix = SVector{N, T}(ntuple(Val(N)) do i
        s = slices[i]
        if s isa KeepAll
            T(full_pixel[i])
        elseif s isa KeepRange
            a = T(first(s.range))
            stp = T(step(s.range))
            (T(full_pixel[i]) - a) / stp + one(T)
        else  # DropAxis — pass through, filtered out by pixel_keep indexing
            T(full_pixel[i])
        end
    end)

    return _extract_kept(all_pix, swcs.pixel_keep)
end

# ── Generic extraction helper ───────────────────────────────────────────────────

"""
    _extract_kept(vec, keep) -> SVector

Return a static vector containing elements of `vec` at indices `keep`.
"""
@inline function _extract_kept(vec::AbstractVector, keep::SVector{K, Int}) where {K}
    T = _coordinate_float_type(vec)
    return SVector{K, T}(ntuple(j -> T(vec[keep[j]]), K))
end
