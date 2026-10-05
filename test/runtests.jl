using WraithSV
using Test
using LinearAlgebra, LsqFit, Distributions

# deterministic pseudo-noise in [-0.5, 0.5) so tests need no extra dependencies
noise(i, j) = mod(sin(12.9898i + 78.233j) * 43758.5453, 1.0) - 0.5

function synthetic(n; a=-1.0, b=2.0, c=0.08, σ=0.01)
    mat = zeros(Float32, n, n)
    for j in 2:n, i in 1:j-1
        mat[i, j] = gompertz(j - i, a, b, c) + σ * noise(i, j)
    end
    return mat
end

# Straightforward O(n²) reference: fit and bands over every individual cell
function bruteforce(mat, level)
    n = size(mat, 1)
    xs = Float64[]; ys = Float64[]
    for j in 2:n, i in 1:j-1
        push!(xs, j - i); push!(ys, mat[i, j])
    end
    model(x, p) = exp.(p[1] .+ p[2] .* exp.(-x .* p[3]))
    fit = curve_fit(model, xs, ys, [1.0, 1.0, 1.0])
    p = fit.param
    ŷ = model(xs, p)
    m = length(ys)
    mse = sum(abs2, ys .- ŷ) / (m - 3)
    e = exp.(-xs .* p[3])
    J = hcat(ŷ, ŷ .* e, ŷ .* (-p[2] .* xs .* e))
    JtJinv = inv(J' * J)
    h = [dot(J[i, :], JtJinv * J[i, :]) for i in 1:m]
    t = quantile(TDist(m - 3), 1 - (1 - level) / 2)
    hw = t .* sqrt.(mse .* (1 .+ h))
    return (; xs, p, ŷ, lower=ŷ .- hw, upper=ŷ .+ hw, mse)
end

@testset "WraithSV.jl" begin
    @testset "analytic Jacobian matches finite differences" begin
        a, b, c = -1.0, 2.0, 0.08
        for x in (1.0, 7.0, 40.0)
            _, ja, jb, jc = gompertzJacobianRow(x, a, b, c)
            h = 1e-6
            fd(f) = (f(h) - f(-h)) / 2h
            @test ja ≈ fd(δ -> gompertz(x, a + δ, b, c)) rtol = 1e-5
            @test jb ≈ fd(δ -> gompertz(x, a, b + δ, c)) rtol = 1e-5
            @test jc ≈ fd(δ -> gompertz(x, a, b, c + δ)) rtol = 1e-5
        end
    end

    @testset "diagonalStats" begin
        mat = synthetic(30)
        s = diagonalStats(mat)
        for d in (1, 5, 29)
            vals = [mat[i, i+d] for i in 1:30-d]
            @test s.k[d] == 30 - d
            @test s.μ[d] ≈ sum(vals) / length(vals)
            @test s.ss[d] ≈ sum(abs2, vals .- s.μ[d]) atol = 1e-9
        end
        z = jaccardScores(mat, s)
        @test z isa Matrix{Float32}
        d = 4
        vals = [mat[i, i+d] for i in 1:30-d]
        zref = (vals .- sum(vals) / length(vals)) ./ sqrt(sum(abs2, vals .- sum(vals) / length(vals)) / (length(vals) - 1))
        @test [z[i, i+d] for i in 1:30-d] ≈ zref rtol = 1e-4
    end

    @testset "reduced fit and bands equal the all-cells fit" begin
        mat = synthetic(60)
        level = 0.95
        ref = bruteforce(mat, level)
        stats = diagonalStats(mat)
        p = fitGompertz(stats)
        @test p ≈ ref.p rtol = 1e-4
        bands = predictionBands(stats, p, level)
        @test bands.mse ≈ ref.mse rtol = 1e-4
        # every cell on a diagonal shares the same fitted value and bands
        for (k, x) in enumerate(ref.xs)
            d = Int(x)
            @test bands.fitted[d] ≈ ref.ŷ[k] rtol = 1e-4
            @test bands.lower[d] ≈ ref.lower[k] rtol = 1e-3 atol = 1e-6
            @test bands.upper[d] ≈ ref.upper[k] rtol = 1e-3 atol = 1e-6
        end
    end

    @testset "detectOutliers flags an injected block" begin
        mat = synthetic(60)
        for j in 40:42, i in 10:12
            mat[i, j] += 0.3f0
        end
        res = detectOutliers(mat)
        @test !isempty(res.outliers.nrow)
        flagged = Set(zip(res.outliers.nrow, res.outliers.ncol))
        @test all((i, j) in flagged for i in 10:12, j in 40:42)
        @test all(res.outliers.is_upper .| res.outliers.is_lower)
        @test all(abs.(res.outliers.z_score) .> 2.0)
    end
end
