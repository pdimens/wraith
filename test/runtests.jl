using WraithSV
using Test
using LinearAlgebra, LsqFit, Distributions
import Mmap

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
    @testset "SV detection" begin
        # clustering equals brute-force single linkage (distance < 3), incl. diagonal-only neighbours
        pts = unique([(rand(1:60), rand(1:60)) for _ in 1:300])
        nr, nc = first.(pts), last.(pts)
        lab = WraithSV.outlierClusters(nr, nc)
        m = length(pts)
        comp = collect(1:m)
        for i in 1:m, j in i+1:m
            (nr[i] - nr[j])^2 + (nc[i] - nc[j])^2 < 9 || continue
            a, b = comp[i], comp[j]
            a == b || (comp .= ifelse.(comp .== b, a, comp))
        end
        @test all((lab[i] == lab[j]) == (comp[i] == comp[j]) for i in 1:m, j in 1:m)
        @test sort(unique(lab)) == 1:maximum(lab)

        # hand-checked: (1,10)-(2,12) are √5 apart (joined); (20,30)-(20,33) are exactly 3 apart (not joined)
        o = (nrow=[1, 2, 20, 20], ncol=[10, 12, 30, 33])
        sv = detectSVs(o; winsize=100)
        @test sv.length == [1300, 1100, 1000]       # sorted, longest first
        @test sv.start == [2000, 100, 2000] && sv.stop == [3300, 1200, 3000]
        @test sv.minrow == [20, 1, 20] && sv.maxcol == [33, 12, 30]

        # nothing in, nothing out
        e = detectSVs((nrow=Int[], ncol=Int[]))
        @test isempty(e.id) && isempty(e.length)

        # end to end: a block of excess barcode sharing is reported as one SV
        mat = synthetic(60)
        for j in 40:48, i in 10:14
            mat[i, j] += 0.3f0
        end
        svs = detectSVs(detectOutliers(mat).outliers)
        @test length(svs.id) >= 1
        @test svs.minrow[1] <= 10 + 2 && svs.maxcol[1] >= 48 - 2

        mktempdir() do dir
            f = writeSVs(joinpath(dir, "sv.csv"), detectSVs(o; winsize=100))
            lines = readlines(f)
            @test lines[1] == "SV_id,start,end,length" && length(lines) == 4
        end
    end

    @testset "non-finite cells are skipped" begin
        mat = synthetic(30)
        clean = diagonalStats(mat)
        mat[3, 7] = NaN
        mat[1, 2] = Inf
        mat[1, 30] = NaN            # the whole single-cell diagonal
        st = diagonalStats(mat)
        @test st.k[4] == clean.k[4] - 1
        @test st.k[1] == clean.k[1] - 1
        @test st.k[29] == 0 && isnan(st.μ[29])
        @test all(isfinite, st.μ[1:28]) && all(isfinite, st.ss)
        z = jaccardScores(mat)
        @test count(!isfinite, z) == 2   # (3,7) NaN and (1,2) Inf; (1,30) sits on an empty diagonal so scores 0
        r = detectOutliers(mat)
        @test all(isfinite, r.params)
        @test r.params ≈ detectOutliers(synthetic(30)).params rtol=0.05
        @test_throws ArgumentError detectOutliers(fill(NaN, 5, 5))
    end

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

    @testset "PackedTriangle" begin
        n = 37
        mat = synthetic(n)
        p = packTriangle(mat)
        @test p isa PackedTriangle
        @test length(p.data) == packedlength(n) == n * (n - 1) ÷ 2
        @test size(p) == (n, n)

        # indexing: column-major packing, zeros on/below the diagonal
        k = 0
        for j in 2:n, i in 1:j-1
            k += 1
            @test p.data[k] == mat[i, j]
            @test p[i, j] == mat[i, j]
        end
        @test all(p[i, i] == 0f0 for i in 1:n)
        @test all(p[j, i] == 0f0 for i in 1:n for j in i+1:n)
        @test_throws ArgumentError (p[5, 3] = 1f0)
        @test_throws BoundsError p[0, 1]
        @test unpack(p) == mat

        # downstream functions give the same answer on packed and dense input
        @test diagonalStats(p).μ ≈ diagonalStats(mat).μ
        zp = jaccardScores(p)
        @test zp isa PackedTriangle
        @test unpack(zp) ≈ jaccardScores(mat)
        rd = detectOutliers(mat)
        rp = detectOutliers(p)
        @test rp.params ≈ rd.params
        @test rp.outliers.nrow == rd.outliers.nrow && rp.outliers.ncol == rd.outliers.ncol

        mktempdir() do dir
            # compressed round trip, in memory and memory-mapped
            f = joinpath(dir, "tri.wpt")
            @test writePacked(f, p) == f
            @test filesize(f) < 4length(p.data)   # smooth synthetic data compresses
            @test unpack(readPacked(f)) == mat
            q = readPacked(f; mmap=true, mmap_path=joinpath(dir, "q.bin"))
            @test q.mapped
            @test unpack(q) == mat

            # construct directly into a memory-mapped file, flush, re-open it
            raw = joinpath(dir, "raw.bin")
            m = packTriangle(mat; mmap=true, path=raw)
            @test m.mapped
            @test unpack(m) == mat
            Mmap.sync!(m)
            @test filesize(raw) == 4length(p.data)
            @test unpack(openMapped(raw, n; readonly=true)) == mat
            @test_throws ArgumentError openMapped(raw, n + 1)

            # temp-file mmap (no path)
            t = packTriangle(mat; mmap=true)
            @test t.mapped && unpack(t) == mat

            # malformed files
            bad = joinpath(dir, "bad.wpt")
            write(bad, "not a packed file")
            @test_throws ArgumentError readPacked(bad)
            trunc = joinpath(dir, "trunc.wpt")
            write(trunc, read(f, filesize(f) - 20))
            @test_throws Exception readPacked(trunc)
        end

        # degenerate sizes
        @test length(PackedTriangle(0).data) == 0
        @test length(PackedTriangle(1; mmap=true).data) == 0
    end
end
