"""
Ancillary API functions for WCS transforms, including dimensionality accessors and axis correlation matrix computation.
"""

# ── Dimensionality accessors ────────────────────────────────────────────────────

"""
    pixel_n_dim(wcs) -> Int

Number of pixel axes.  A compile-time constant (extracted from the type parameter).
"""
pixel_n_dim(::WCSTransform{N}) where {N} = N
pixel_n_dim(::SlicedWCSTransform{N, Np}) where {N, Np} = Np

"""
    world_n_dim(wcs) -> Int

Number of world axes.  A compile-time constant (extracted from the type parameter).
"""
world_n_dim(::WCSTransform{N}) where {N} = N
world_n_dim(::SlicedWCSTransform{N, Np, Nw}) where {N, Np, Nw} = Nw

# ── Axis correlation matrix ─────────────────────────────────────────────────────

"""
    axis_correlation_matrix(wcs) -> BitMatrix

Return a `(world_n_dim, pixel_n_dim)` boolean matrix where entry `[w, p]` is
`true` if world axis `w` depends on pixel axis `p`.

For `WCSTransform`: derived from the CD matrix and projection structure,
matching the APE-14 convention.  For `SlicedWCSTransform`: the parent's
matrix indexed by the kept axes.
"""
function axis_correlation_matrix(wcs::WCSTransform{N}) where {N}
    # If any pre-linear distortion is present (SIP, D2IM, CPDIS), we must be
    # conservative: distortions can introduce arbitrary cross-axis coupling
    # that isn't captured by the CD matrix alone.
    if has_distortion(wcs.pipeline)
        return fill(true, SMatrix{N, N, Bool})
    end

    li = wcs.lon_axis
    la = wcs.lat_axis
    has_celestial = li > 0 && la > 0

    # Build the correlation SMatrix in one pass.  For celestial axes the
    # longitude and latitude share pixel dependencies (spherical rotation
    # couples them), so those rows take the union of the CD-based booleans.
    return SMatrix{N, N, Bool}(ntuple(k -> begin
        i = (k - 1) ÷ N + 1
        j = (k - 1) % N + 1
        if has_celestial && (i == li || i == la)
            wcs.cd[li, j] != 0 || wcs.cd[la, j] != 0
        else
            wcs.cd[i, j] != 0
        end
    end, Val(N * N)))
end

function axis_correlation_matrix(swcs::SlicedWCSTransform)
    full = axis_correlation_matrix(swcs.parent)
    return full[swcs.world_keep, swcs.pixel_keep]
end
