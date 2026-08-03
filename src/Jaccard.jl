using Statistics, LsqFit, Distributions, LinearAlgebra, StaticArrays


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

"""
Using the Jaccard identity matrix (intersect/union) as input,
calculate the Z-scores. Returns a Z score matrix.
"""
function jaccardScores(mat::Matrix{Float32})::Matrix{Float64}
    n = size(mat, 1)
    z_scores = similar(mat, axes(mat))
    @inbounds for diag in 1:(n-1)
        diag_start = diag + 1
        diag_len = n - diag
        μ = 0.0
        for j in diag_start:n
            μ += mat[j-diag, j]
        end
        μ /= diag_len

        σ² = 0.0
        for j in diag_start:n
            σ² += (mat[j-diag, j] - μ)^2
        end
        σ = sqrt(σ² / (diag_len - 1))

        for j in diag_start:n
            z_scores[j-diag, j] = σ == 0 ? 0.0 : (mat[j-diag, j] - μ) / σ
        end
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
Fit the model y = exp(a + b*exp(-x*c)) to (xs, ys) using nonlinear least squares.
Returns the fitted parameter vector [a, b, c].
"""
function fitGompertz(xs::Vector{Float64}, ys::Vector{Float64})::Vector{Float64}
    model(x, p) = exp.(p[1] .+ p[2] .* exp.(-x .* p[3]))
    p0 = [1.0, 1.0, 1.0]
    fit = curve_fit(model, xs, ys, p0)
    return fit.param
end


"""
    predictionBands(xs, ys, params, level=0.95) -> (ŷ, lower, upper)

Compute pointwise prediction bands for the Gompertz model
`y = exp(a + b·exp(-x·c))` fitted to observations `(xs, ys)`.

Each band is centred on the fitted value `ŷᵢ` with half-width
`t_{α/2, n-3} · sqrt(MSE · (1 + hᵢ))`, where:
- `MSE = RSS / (n - 3)` is the residual mean squared error
- `hᵢ = Jᵢᵀ (JᵀJ)⁻¹ Jᵢ` is the leverage of point `i`
- `Jᵢ` is the analytic Jacobian row `[∂f/∂a, ∂f/∂b, ∂f/∂c]` evaluated at `xᵢ`

# Arguments
- `xs::Vector{Float64}`: predictor values (diagonal distances from the matrix diagonal)
- `ys::Vector{Float64}`: response values (Jaccard similarities)
- `params::Vector{Float64}`: fitted parameters `[a, b, c]` from [`fitGompertz`](@ref)
- `level::Float64`: confidence level for the prediction bands (default: `0.95`)

# Returns
A `Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}` of equal-length vectors:
- `ŷ`: fitted values
- `lower`: lower prediction bound
- `upper`: upper prediction bound

# Performance notes
All intermediate values (`Jᵢ`, `JᵀJ`, `(JᵀJ)⁻¹`) are stack-allocated via
`StaticArrays`. The only heap allocations are the three returned output vectors.
"""
function predictionBands(
    xs::Vector{Float64},
    ys::Vector{Float64},
    params::Vector{Float64},
    level::Float64=0.95
)::Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}

    n = length(ys)
    a, b, c = params[1], params[2], params[3]
    α = 1.0 - level

    # First pass: accumulate JᵀJ, RSS, and ŷ together
    JtJ = @MMatrix zeros(3, 3)   # heap-allocated, but fast(er) mutable 3×3
    rss = 0.0
    ŷ = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        x = xs[i]
        e = exp(-x * c)
        fi = exp(a + b * e)         # ŷᵢ
        ŷ[i] = fi

        # Jacobian row
        Ji = SVector{3,Float64}(fi, fi * e, fi * (-b * x * e))

        # Accumulate JᵀJ in-place: rank-1 update
        for r in 1:3, s in r:3
            v = Ji[r] * Ji[s]
            JtJ[r, s] += v
            if r != s
                JtJ[s, r] += v
            end
        end

        r = ys[i] - fi
        rss += r * r
    end

    mse = rss / (n - 3)
    JtJ_inv = inv(SMatrix{3,3}(JtJ))   # static 3×3 inversion, no heap alloc
    t_crit = quantile(TDist(n - 3), 1.0 - α / 2.0)

    # Second pass: compute per-point leverage and bands
    lower = Vector{Float64}(undef, n)
    upper = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        x = xs[i]
        e = exp(-x * c)
        fi = ŷ[i]
        Ji = SVector{3,Float64}(fi, fi * e, fi * (-b * x * e))

        # hᵢ = Jᵢᵀ (JᵀJ)⁻¹ Jᵢ
        v = JtJ_inv * Ji
        hi = dot(Ji, v)
        half_w = t_crit * sqrt(mse * (1.0 + hi))
        lower[i] = fi - half_w
        upper[i] = fi + half_w
    end

    return (ŷ, lower, upper)
end

# this one allocates a lot, perhaps don't use it, see function above
function predictionBands(xs::Vector{Float64}, ys::Vector{Float64}, params::Vector{Float64}, level::Float64=0.95)::Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}
    model(x, p) = exp.(p[1] .+ p[2] .* exp.(-x .* p[3]))

    n = length(ys)
    p = length(params)
    ŷ = model(xs, params)

    # Residual mean squared error
    residuals = ys .- ŷ
    mse = dot(residuals, residuals) / (n - p)

    # Jacobian of the model at each point: shape (n × p)
    # ∂/∂a = exp(a + b*exp(-xc))          = ŷ
    # ∂/∂b = exp(a + b*exp(-xc))*exp(-xc) = ŷ .* exp.(-xs.*params[3])
    # ∂/∂c = exp(a + b*exp(-xc))*(-b*x*exp(-xc)) = ŷ .* (-params[2] .* xs .* exp.(-xs.*params[3]))
    e = exp.(-xs .* params[3])
    J = hcat(ŷ, ŷ .* e, ŷ .* (-params[2] .* xs .* e))   # n × 3

    # (JᵀJ)⁻¹  — small 3×3 matrix, direct inversion is fine
    JtJ_inv = inv(J' * J)

    # Leverage scores: hᵢ = Jᵢᵀ (JᵀJ)⁻¹ Jᵢ  (one scalar per point)
    leverage = [dot(J[i, :], JtJ_inv * J[i, :]) for i in 1:n]

    # Prediction std: sqrt(MSE * (1 + hᵢ))
    pred_std = sqrt.(mse .* (1.0 .+ leverage))

    # Two-sided t critical value
    α = 1.0 - level
    t_crit = quantile(TDist(n - p), 1.0 - α / 2.0)

    @inbounds begin
        lower = ŷ .- t_crit .* pred_std
        upper = ŷ .+ t_crit .* pred_std
    end
    return (ŷ, lower, upper)
end

"""
Given the Jaccard matrix, return a NamedTuple of flat vectors ready for
outlier detection, mirroring the R `dataset` dataframe.

Fields: nrow, ncol, diag_dist, value, z_score, abs_z, fitted, lower, upper
"""
function buildOutlierDataset(
    mat::Matrix{Float32},
    z_scores::Matrix{Float64},
    prediction_level::Float64=0.95
)
    n = size(mat, 1)
    m = n * (n - 1) ÷ 2

    # pre-allocate output vectors because we know what the size will be based on the input matrix
    rows = Vector{Int}(undef, m)
    cols = Vector{Int}(undef, m)
    dists = Vector{Int}(undef, m)
    values = Vector{Float64}(undef, m)
    zs = Vector{Float64}(undef, m)

    k = 0
    @inbounds for j in 1:n, i in 1:(j-1)
        k += 1
        rows[k] = i
        cols[k] = j
        dists[k] = j - i
        values[k] = mat[i, j]
        zs[k] = z_scores[i, j]
    end

    xs = Float64.(dists)
    ys = values

    params = fitGompertz(xs, ys)
    ŷ, lower, upper = predictionBands(xs, ys, params, prediction_level)

    return (
        nrow=rows,
        ncol=cols,
        x=dists,
        y=values,
        z_score=zs,
        abs_z=abs.(zs),
        fitted=ŷ,
        lower=lower,
        upper=upper,
    )
end