"""
Cells `(nrow, ncol)` of the upper triangle are 2-D integer points. Single-linkage
agglomerative clustering with a distance cut-off is exactly the connected
components of the graph that joins every pair of points closer than `distance`
(sklearn merges only when the linkage distance is *below* the threshold). Since
the points sit on an integer grid, each point only has to look at the few grid
cells within `distance` of it, found through a hash of occupied cells and merged
with a union-find, so the cost is O(m) for `m` outliers rather than the O(m²) of
a pairwise distance matrix.
"""
@inline function findroot!(parent::Vector{Int}, i::Int)
    @inbounds while parent[i] != i
        parent[i] = parent[parent[i]]   # path halving
        i = parent[i]
    end
    return i
end

"""
    outlierClusters(nrow, ncol; distance=3) -> Vector{Int}

Cluster id (`1…k`, numbered by first appearance) of every outlier, grouping cells
whose Euclidean distance in `(ncol, nrow)` space is `< distance`, transitively
(single linkage). Matches `AgglomerativeClustering(n_clusters=None,
distance_threshold=distance, linkage="single")` from the original Python.
"""
function outlierClusters(nrow::AbstractVector{<:Integer}, ncol::AbstractVector{<:Integer}; distance::Real=3)
    m = length(nrow)
    length(ncol) == m || throw(DimensionMismatch("nrow and ncol must have the same length"))
    m == 0 && return Int[]
    R = ceil(Int, distance) - 1        # largest integer offset that can still be < distance
    d² = distance^2
    # half of the neighbourhood suffices, because every pair is seen from one side
    offsets = [(dx, dy) for dx in 0:R for dy in -R:R if (dx > 0 || dy > 0) && dx^2 + dy^2 < d²]

    minr, minc = minimum(nrow), minimum(ncol)
    stride = Int(maximum(nrow)) - Int(minr) + 2R + 1
    key(r, c) = (Int(c) - Int(minc) + R) * stride + (Int(r) - Int(minr) + R)

    cells = Dict{Int,Int}()
    sizehint!(cells, m)
    parent = collect(1:m)
    @inbounds for p in 1:m
        k = key(nrow[p], ncol[p])
        if haskey(cells, k)             # repeated cell: same location, same cluster
            parent[findroot!(parent, p)] = findroot!(parent, cells[k])
            continue
        end
        cells[k] = p
    end
    @inbounds for p in 1:m
        r, c = Int(nrow[p]), Int(ncol[p])
        for (dx, dy) in offsets
            q = get(cells, key(r + dy, c + dx), 0)
            q == 0 && continue
            rp, rq = findroot!(parent, p), findroot!(parent, q)
            rp == rq || (parent[max(rp, rq)] = min(rp, rq))
        end
    end

    ids = zeros(Int, m)
    labels = Vector{Int}(undef, m)
    k = 0
    @inbounds for p in 1:m
        root = findroot!(parent, p)
        ids[root] == 0 && (ids[root] = (k += 1))
        labels[p] = ids[root]
    end
    return labels
end

"""
    detectSVs(outliers; winsize=1, distance=3) -> NamedTuple

Group the outliers returned by [`findOutliers`](@ref) (or `detectOutliers(...).outliers`;
only its `nrow` and `ncol` are used) into putative structural variants, a port of
WRATH's `sv_detection.py`. Outlier cells closer than `distance` windows are
clustered (single linkage), and each cluster gives one SV spanning
`minrow … maxcol`, so `length = maxcol - minrow` windows. `winsize` converts
window indices to base pairs for `start`, `stop` and `length`.

Returns equal-length vectors, sorted by `length` (longest first):
`(id, start, stop, length, minrow, maxrow, mincol, maxcol)`, where the last four
are the cluster bounds in window units. No outliers gives empty vectors.
"""
function detectSVs(outliers; winsize::Integer=1, distance::Real=3)
    (; nrow, ncol) = outliers
    labels = outlierClusters(nrow, ncol; distance)
    k = maximum(labels; init=0)
    minrow = fill(typemax(Int), k); maxrow = fill(typemin(Int), k)
    mincol = fill(typemax(Int), k); maxcol = fill(typemin(Int), k)
    @inbounds for p in eachindex(labels)
        g, r, c = labels[p], Int(nrow[p]), Int(ncol[p])
        minrow[g] = min(minrow[g], r); maxrow[g] = max(maxrow[g], r)
        mincol[g] = min(mincol[g], c); maxcol[g] = max(maxcol[g], c)
    end
    len = maxcol .- minrow
    o = sortperm(len; rev=true, alg=MergeSort)    # stable, so ties keep first-seen order
    return (id=o, start=minrow[o] .* winsize, stop=maxcol[o] .* winsize, length=len[o] .* winsize,
        minrow=minrow[o], maxrow=maxrow[o], mincol=mincol[o], maxcol=maxcol[o])
end

"""
    writeSVs(path, svs) -> path

Write the result of [`detectSVs`](@ref) as CSV with the columns
`SV_id,start,end,length`, like the original script.
"""
function writeSVs(path::AbstractString, svs)
    open(path, "w") do io
        println(io, "SV_id,start,end,length")
        for i in eachindex(svs.id)
            println(io, svs.id[i], ',', svs.start[i], ',', svs.stop[i], ',', svs.length[i])
        end
    end
    return path
end
