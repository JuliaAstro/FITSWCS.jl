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

    # Base coupling from the CD matrix: world axis i depends on pixel axis j
    # if CD[i,j] is non-zero.  Convert to mutable BitMatrix so we can mutate
    # rows for celestial sharing below.
    matrix = MMatrix{N, N, Bool}(ntuple(k -> begin
        i = (k - 1) ÷ N + 1
        j = (k - 1) % N + 1
        wcs.cd[i, j] != 0
    end, N * N))

    # Celestial longitude and latitude share pixel dependencies because
    # spherical rotation couples them.  Union their dependencies.
    if wcs.lon_axis > 0 && wcs.lat_axis > 0
        li = wcs.lon_axis
        la = wcs.lat_axis
        coupled = matrix[li, :] .| matrix[la, :]
        matrix[li, :] .= coupled
        matrix[la, :] .= coupled
    end

    return SMatrix{N, N, Bool}(matrix)
end

function axis_correlation_matrix(swcs::SlicedWCSTransform)
    full = axis_correlation_matrix(swcs.parent)
    return full[swcs.world_keep, swcs.pixel_keep]
end
