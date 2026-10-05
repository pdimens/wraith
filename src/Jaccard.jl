using Statistics, LsqFit, Distributions, LinearAlgebra


"""
Return length(intersection(x,y)), length(union(x,y)) for Sets `x` and `y`
without creating (allocating) the union and intersection sets.
"""
function jaccardIdent(x::Set{String}, y::Set{String})::Float64
    _intersect = count(i -> i ∈ y, x)
    _union = length(x) + length(y) - _intersect
    return _intersect / _union
end

"""
Return a generator traversing an upper triangle, omitting the diagonal.
Returns (i,j) indices, intended to be used as:
```
for (i,j) in uppertriangle(mat)
        ...
    end
```
"""
uppertriangle(A) = ((i, j) for j in axes(A, 2) for i in 1:j-1)

"""
Given a column index `col` to start, creates a generator
that traverses the diagonal of Matrix `A`
"""
diagonal(A, col::Int) = (A[i, col+i] for i in 1:size(A, 1)-1)

# ---------------------------------------------------------------------------
# Everything downstream of the Jaccard matrix depends only on the distance `d`
# from the diagonal, never on the individual cell. A matrix with `n` windows has
# n(n-1)/2 cells but only n-1 distinct distances, so all fitting and band
# calculations below work on O(n) per-diagonal summaries instead of O(n²) cells.
# ---------------------------------------------------------------------------

"""
Per-diagonal summary of the upper triangle of an `n × n` matrix. Every vector
has length `n-1` and is indexed by the distance `d = j - i` from the diagonal.

- `μ[d]`: mean of diagonal `d`
- `ss[d]`: within-diagonal sum of squares, `Σ(y - μ[d])²`
- `k[d]`: number of cells on diagonal `d` (`n - d`)
"""
struct DiagonalStats
    n::Int
    μ::Vector{Float64}
    ss::Vector{Float64}
    k::Vector{Int}
end

"""
    diagonalStats(mat) -> DiagonalStats

Compute the mean and sum of squares of every diagonal of the upper triangle in a
single column-major (cache-friendly) pass using O(n) extra memory. Sums are
accumulated relative to the first element of each diagonal to avoid
catastrophic cancellation.
"""
function diagonalStats(mat::AbstractMatrix{<:Real})::DiagonalStats
    n = size(mat, 1)
    size(mat, 2) == n || throw(DimensionMismatch("matrix must be square"))
    nd = n - 1
    shift = Vector{Float64}(undef, nd)
    s1 = zeros(nd)
    s2 = zeros(nd)
    @inbounds for d in 1:nd
        shift[d] = mat[1, 1+d]
    end
    @inbounds for j in 2:n
        for i in 1:j-1
            d = j - i
            δ = mat[i, j] - shift[d]
            s1[d] += δ
            s2[d] += δ * δ
        end
    end
    k = [n - d for d in 1:nd]
    μ = shift .+ s1 ./ k
    ss = max.(s2 .- s1 .^ 2 ./ k, 0.0)
    return DiagonalStats(n, μ, ss, k)
end

"""
Per-diagonal sample standard deviation (`n-1` denominator, like R's `scale()`).
`NaN` for the single-cell diagonal.
"""
diagonalSD(s::DiagonalStats) = sqrt.(s.ss ./ (s.k .- 1))

@inline zscore(y, μ, σ) = σ == 0 ? 0.0 : (y - μ) / σ

"""
Using the Jaccard identity matrix (intersect/union) as input, calculate the
Z-scores of each upper-triangle cell relative to its diagonal. Returns a
`Float32` matrix (lower triangle and diagonal are zero). If you only need to
flag outliers, prefer [`detectOutliers`](@ref), which never stores this matrix.
"""
function jaccardScores(mat::AbstractMatrix{<:Real}, stats::DiagonalStats=diagonalStats(mat))::Matrix{Float32}
    n = stats.n
    σ = diagonalSD(stats)
    z_scores = zeros(Float32, n, n)
    @inbounds for j in 2:n, i in 1:j-1
        d = j - i
        z_scores[i, j] = zscore(mat[i, j], stats.μ[d], σ[d])
    end
    return z_scores
end

# Much less code, about the same speed, significantly more allocations
function jaccardScores2(mat::Matrix{Float32})::Matrix{Float64}
    z_scores = similar(mat, axes(mat))
    @inbounds for col in 1:(size(mat, 1)-1)
        idx = diagind(mat, col)
        diag = @view mat[idx]
        μ = mean(diag)
        σ = stdm(diag, μ)
        for i in idx
            z_scores[i] = σ == 0 ? 0.0 : (mat[i] - μ) / σ
        end
    end
    return z_scores
end

"""
    gompertz(x, a, b, c)

The model `exp(a + b·exp(-c·x))`.
"""
@inline gompertz(x, a, b, c) = exp(a + b * exp(-c * x))

"""
    gompertzJacobianRow(x, a, b, c) -> (f, ∂a, ∂b, ∂c)

Model value and analytic gradient at `x`:
`∂a = f`, `∂b = f·e`, `∂c = -f·b·x·e`, where `e = exp(-c·x)`.
"""
@inline function gompertzJacobianRow(x, a, b, c)
    e = exp(-c * x)
    f = exp(a + b * e)
    return (f, f, f * e, -f * b * x * e)
end

"""
    fitGompertz(stats; p0=[1.0, 1.0, 1.0]) -> Vector{Float64}

Fit `y = exp(a + b·exp(-x·c))` to every cell of the matrix, using only its
[`DiagonalStats`](@ref).

For a model that depends on `x` alone,
`Σ(y - f(x))² = Σ(y - μ[d])² + k[d]·(μ[d] - f(d))²` per diagonal, and the first
term does not depend on the parameters. The least-squares fit over all
n(n-1)/2 cells is therefore *exactly* a fit of the n-1 diagonal means weighted
by the diagonal lengths `k`, so the parameters match a fit on every cell while
costing O(n) time and memory. The analytic Jacobian is supplied to avoid
finite differences.
"""
function fitGompertz(stats::DiagonalStats; p0::AbstractVector{Float64}=[1.0, 1.0, 1.0])::Vector{Float64}
    nd = length(stats.μ)
    xs = collect(1.0:nd)
    sw = sqrt.(Float64.(stats.k))   # weights enter the cost as sw²
    ys = sw .* stats.μ
    model(x, p) = [sw[i] * gompertz(x[i], p[1], p[2], p[3]) for i in eachindex(x)]
    function jacobian(x, p)
        J = Matrix{Float64}(undef, length(x), 3)
        @inbounds for i in eachindex(x)
            _, ja, jb, jc = gompertzJacobianRow(x[i], p[1], p[2], p[3])
            J[i, 1] = sw[i] * ja
            J[i, 2] = sw[i] * jb
            J[i, 3] = sw[i] * jc
        end
        return J
    end
    fit = curve_fit(model, jacobian, xs, ys, collect(p0))
    fit.converged || @warn "Gompertz fit did not converge; bands may be unreliable"
    return fit.param
end

"""
    predictionBands(stats, params, level=0.95)

Pointwise prediction bands for the Gompertz model with parameters
`params = [a, b, c]`, returned per diagonal (vectors of length `n-1` indexed by
distance `d`) as a NamedTuple `(fitted, lower, upper, est_error, rss, mse, df, tcrit)`.

Each band is `ŷ(d) ± t · sqrt(MSE · (1 + h(d)))` where
- `MSE = RSS / (m - 3)` with `m = n(n-1)/2` cells,
- `h(d) = J(d)ᵀ (JᵀJ)⁻¹ J(d)` is the leverage, with
  `JᵀJ = Σ_d k[d] · J(d) J(d)ᵀ` (every cell on a diagonal shares the same Jacobian row),
- `RSS = Σ ss[d] + Σ k[d]·(μ[d] - ŷ(d))²`.

`h` is computed through a hand-rolled 3×3 Cholesky factorisation (`JᵀJ = LLᵀ`,
`h = ‖L⁻¹J‖²`), which is better conditioned than forming `(JᵀJ)⁻¹` since the
`b` and `c` columns of J are nearly collinear. No heap allocation besides the
output vectors.
"""
function predictionBands(stats::DiagonalStats, params::AbstractVector{Float64}, level::Float64=0.95)
    nd = length(stats.μ)
    nd >= 3 || throw(ArgumentError("need at least 4 windows (3 distinct diagonals) to fit 3 parameters"))
    a, b, c = params
    m = sum(stats.k)
    df = m - 3

    # Pass 1: JᵀJ (upper triangle, weighted by diagonal length), fitted values, RSS
    s11 = s12 = s13 = s22 = s23 = s33 = 0.0
    rss = sum(stats.ss)
    fitted = Vector{Float64}(undef, nd)
    @inbounds for d in 1:nd
        f, j1, j2, j3 = gompertzJacobianRow(d, a, b, c)
        fitted[d] = f
        w = stats.k[d]
        s11 += w * j1 * j1
        s12 += w * j1 * j2
        s13 += w * j1 * j3
        s22 += w * j2 * j2
        s23 += w * j2 * j3
        s33 += w * j3 * j3
        rss += w * (stats.μ[d] - f)^2
    end
    mse = rss / df

    # Cholesky JᵀJ = L Lᵀ
    l11 = sqrt(s11)
    l21 = s12 / l11
    l31 = s13 / l11
    p22 = s22 - l21 * l21
    p22 > 0 || throw(ArgumentError("JᵀJ is not positive definite; check the fitted parameters"))
    l22 = sqrt(p22)
    l32 = (s23 - l31 * l21) / l22
    p33 = s33 - l31 * l31 - l32 * l32
    p33 > 0 || throw(ArgumentError("JᵀJ is not positive definite; check the fitted parameters"))
    l33 = sqrt(p33)

    tcrit = quantile(TDist(df), 1.0 - (1.0 - level) / 2.0)

    # Pass 2: leverage and bands, one value per diagonal
    lower = Vector{Float64}(undef, nd)
    upper = Vector{Float64}(undef, nd)
    est_error = Vector{Float64}(undef, nd)
    @inbounds for d in 1:nd
        f, j1, j2, j3 = gompertzJacobianRow(d, a, b, c)
        w1 = j1 / l11
        w2 = (j2 - l21 * w1) / l22
        w3 = (j3 - l31 * w1 - l32 * w2) / l33
        h = w1 * w1 + w2 * w2 + w3 * w3
        se = sqrt(mse * (1.0 + h))
        est_error[d] = se
        lower[d] = f - tcrit * se
        upper[d] = f + tcrit * se
    end
    return (; fitted, lower, upper, est_error, rss, mse, df, tcrit)
end

"""
    findOutliers(mat, stats, bands; zthreshold=2.0)

Stream over the upper triangle and keep only cells that are outside the
prediction bands of their diagonal **and** have `|z| > zthreshold` (the WRATH
criterion). Nothing proportional to n² is allocated beyond the (small) result.

Returns a NamedTuple of equal-length vectors:
`(nrow, ncol, x, value, z_score, fitted, lower, upper, est_error, is_upper, is_lower)`.
`NaN` cells are never flagged.
"""
function findOutliers(mat::AbstractMatrix{<:Real}, stats::DiagonalStats, bands; zthreshold::Float64=2.0)
    n = stats.n
    σ = diagonalSD(stats)
    out = (nrow=Int[], ncol=Int[], x=Int[], value=Float64[], z_score=Float64[],
        fitted=Float64[], lower=Float64[], upper=Float64[], est_error=Float64[],
        is_upper=Bool[], is_lower=Bool[])
    @inbounds for j in 2:n, i in 1:j-1
        d = j - i
        y = Float64(mat[i, j])
        up = y > bands.upper[d]
        lo = y < bands.lower[d]
        (up || lo) || continue
        z = zscore(y, stats.μ[d], σ[d])
        abs(z) > zthreshold || continue
        push!(out.nrow, i); push!(out.ncol, j); push!(out.x, d)
        push!(out.value, y); push!(out.z_score, z)
        push!(out.fitted, bands.fitted[d]); push!(out.lower, bands.lower[d])
        push!(out.upper, bands.upper[d]); push!(out.est_error, bands.est_error[d])
        push!(out.is_upper, up); push!(out.is_lower, lo)
    end
    return out
end

"""
    detectOutliers(mat; level=0.95, zthreshold=2.0) -> (outliers, bands, params)

End-to-end WRATH outlier detection on a Jaccard matrix: per-diagonal statistics,
Gompertz fit, prediction bands, then a streaming scan for outliers. Peak extra
memory is O(n) plus the outliers themselves, instead of O(n²).
"""
function detectOutliers(mat::AbstractMatrix{<:Real}; level::Float64=0.95, zthreshold::Float64=2.0)
    stats = diagonalStats(mat)
    params = fitGompertz(stats)
    bands = predictionBands(stats, params, level)
    return (outliers=findOutliers(mat, stats, bands; zthreshold), bands=bands, params=params)
end
