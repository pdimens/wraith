module WraithSV

include("Jaccard.jl")
include("PackedTriangle.jl")
include("DetectSV.jl")
include("BAM.jl")
include("Windows.jl")

export jaccardIdent, jaccardScores, DiagonalStats, diagonalStats, diagonalSD,
    gompertz, gompertzJacobianRow, fitGompertz, predictionBands, findOutliers,
    detectOutliers, PackedTriangle, packedlength, packTriangle, unpack,
    writePacked, readPacked, openMapped, detectSVs, writeSVs, wrath

"""
    wraithsv(bams, genome, chromosome; winsize=50_000, start=nothing, stop=nothing,
          outdir="wrath_out", threads=1, detect=false, plot=true, step=nothing)

Entry point for the whole WRATH workflow, the Julia counterpart of the `wrath`
bash script: build a barcode-sharing (Jaccard) matrix between genomic windows of
one chromosome, optionally flag outliers and call SVs from it.

- `bams`: BAM files (with `.bai` indexes) of the population/phenotype of interest
- `genome`: reference FASTA (its `.fai` supplies the chromosome lengths)
- `chromosome`: contig to analyse; `start`/`stop` optionally subset it
- `winsize`: window size in bp
- `detect`: also run outlier detection and SV calling (`-l` in the script)
- `plot`: draw the heatmap(s) (disabled with `-p` in the script)
- `step`: resume from `:windows`, `:barcodes`, `:matrix`, `:outliers` or `:plot`,
  reusing the files an earlier run left in `outdir`

This is a scaffold: stages marked `TODO` are not implemented yet.
"""
function wraithsv(bams::AbstractVector{<:AbstractString}, genome::AbstractString, chromosome::AbstractString;
    winsize::Integer=50_000, start::Union{Nothing,Integer}=nothing, stop::Union{Nothing,Integer}=nothing,
    outdir::AbstractString="wrath_out", threads::Integer=1, detect::Bool=false, plot::Bool=true,
    step::Union{Nothing,Symbol}=nothing)

    step in (nothing, :windows, :barcodes, :matrix, :outliers, :plot) ||
        throw(ArgumentError("unknown step $(repr(step)); use :windows, :barcodes, :matrix, :outliers or :plot"))
    # a stage runs when no `step` was given, or when it is the requested one or comes after it
    stages = (:windows, :barcodes, :matrix, :outliers, :plot)
    runs(s) = step === nothing || findfirst(==(s), stages) >= findfirst(==(step), stages)

    mkpath.(joinpath.(outdir, ("beds", "matrices", "outliers", "SVs", "plots")))
    tag = "$(winsize)_$(chromosome)_$(something(start, 1))_$(something(stop, "end"))"

    # 1. Windows -------------------------------------------------------------
    # Read the chromosome length from `genome * ".fai"`, then tile
    # `start:winsize:stop` (default: the whole contig) into windows. Windows.jl
    # has `GenomicWindows` for this. Keep the windows in memory and also write
    # `outdir/beds/windows_<tag>.bed` so `step=:matrix` can resume from it.
    runs(:windows) && error("TODO: build windows")

    # 2. Barcodes ------------------------------------------------------------
    # For every BAM, stream the reads on `chromosome:start-stop` (BAM.jl:
    # `BamReader`, `fetch`, `getValidBX`), keep MAPQ >= 20 and a valid BX tag,
    # and collect `(position, barcode)` pairs. Combine all samples into one
    # position-sorted table per window; if the barcodes are not kept in memory,
    # write `outdir/beds/barcodes_<tag>.bed.gz` (+ index) as the script does.
    runs(:barcodes) && error("TODO: collect barcodes per window")

    # 3. Jaccard matrix ------------------------------------------------------
    # For every pair of windows (i < j) compute |A ∩ B| / |A ∪ B| of their barcode
    # sets with `jaccardIdent`, parallelised over `threads` (rows are independent).
    # Windows with no barcodes must give 0, not NaN. Store the strict upper
    # triangle in a `PackedTriangle` (`mmap=true` for large chromosomes) and save
    # it with `writePacked` to `outdir/matrices/jaccard_matrix_<tag>.wpt`.
    runs(:matrix) && error("TODO: build the Jaccard matrix")
    # mat = PackedTriangle(length(windows); mmap = ...) ; fill it ; writePacked(...)
    # on resume: mat = readPacked(joinpath(outdir, "matrices", "jaccard_matrix_$tag.wpt"))

    # 4. Outliers and SVs ----------------------------------------------------
    # Per-diagonal z-scores, the Gompertz fit and its prediction bands, then keep
    # the cells outside the bands with |z| > 2, and cluster them into SVs.
    # res = detectOutliers(mat; level = 0.95, zthreshold = 2.0)
    # svs = detectSVs(res.outliers; winsize)
    # writeSVs(joinpath(outdir, "SVs", "sv_$tag.csv"), svs)
    # Also write `res.outliers` to `outdir/outliers/outliers_<tag>.csv`.
    detect && runs(:outliers) && error("TODO: outliers and SV calling (detectOutliers, detectSVs, writeSVs)")

    # 5. Plots ---------------------------------------------------------------
    # Heatmap of log(jaccard + 1e-4) * 100 with the windows on the axes; when SVs
    # were called, mark them on top. Written to `outdir/plots/heatmap_<tag>.png`.
    plot && runs(:plot) && error("TODO: heatmap")

    return nothing
end

end
