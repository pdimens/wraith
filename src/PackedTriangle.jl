import Mmap
using CodecZlib

"""
    PackedTriangle <: AbstractMatrix{Float32}

The strict upper triangle of an `n × n` matrix (diagonal excluded) stored as a
flat `Vector{Float32}` of length `n(n-1)/2`, about half the memory of a dense
matrix. Cell `(i, j)` with `i < j` lives at index `(j-1)(j-2)÷2 + i`, i.e. in
column-major order, so walking `for j in 2:n, i in 1:j-1` reads the data
sequentially. Reads of the diagonal or lower triangle return `0f0`; writing
there throws.

The backing vector can be an ordinary in-memory `Vector` or a memory-mapped
file (see [`PackedTriangle(n; mmap, path)`](@ref)). Because it is an
`AbstractMatrix`, it works with [`diagonalStats`](@ref), [`findOutliers`](@ref),
[`detectOutliers`](@ref), etc.
"""
struct PackedTriangle <: AbstractMatrix{Float32}
    n::Int
    data::Vector{Float32}
    mapped::Bool
end

"Length of the packed vector for an `n × n` matrix: `n(n-1)/2`."
packedlength(n::Integer) = n * (n - 1) ÷ 2

"Index into the packed vector of cell `(i, j)`, requires `1 ≤ i < j`."
@inline trilindex(i::Integer, j::Integer) = ((j - 1) * (j - 2)) ÷ 2 + i

Base.size(p::PackedTriangle) = (p.n, p.n)

@inline function Base.getindex(p::PackedTriangle, i::Int, j::Int)
    @boundscheck checkbounds(p, i, j)
    return i < j ? (@inbounds p.data[trilindex(i, j)]) : 0f0
end

@inline function Base.setindex!(p::PackedTriangle, v, i::Int, j::Int)
    @boundscheck checkbounds(p, i, j)
    i < j || throw(ArgumentError("only the strict upper triangle (i < j) is stored; cannot set ($i, $j)"))
    @inbounds p.data[trilindex(i, j)] = v
    return v
end

# map a file as a Vector{Float32} of length `len`; `create` truncates/creates it
function mapfile(path::AbstractString, len::Int; create::Bool, readonly::Bool=false)
    io = open(path, create ? "w+" : (readonly ? "r" : "r+"))
    try
        return Mmap.mmap(io, Vector{Float32}, len)
    finally
        close(io)   # the mapping stays valid after the file handle is closed
    end
end

"""
    PackedTriangle(n; mmap=false, path=nothing)

Allocate a zero-filled packed triangle for an `n × n` matrix.

- `mmap=false`: backed by an ordinary in-memory vector (`path` is ignored).
- `mmap=true`: backed by a memory-mapped file at `path`, so the OS pages it in
  and out and it can exceed available RAM. An existing file at `path` is
  **overwritten**. If `path` is omitted a temporary file is used and, on Unix,
  unlinked immediately (the mapping stays valid and the disk space is released
  when the object is garbage collected). Use [`openMapped`](@ref) to re-open a
  file you kept.
"""
function PackedTriangle(n::Integer; mmap::Bool=false, path::Union{Nothing,AbstractString}=nothing)
    n >= 0 || throw(ArgumentError("n must be non-negative"))
    len = packedlength(n)
    (mmap && len > 0) || return PackedTriangle(n, zeros(Float32, len), false)
    tmp = path === nothing
    file = tmp ? tempname() : String(path)
    data = mapfile(file, len; create=true)
    tmp && Sys.isunix() && rm(file; force=true)
    return PackedTriangle(n, data, true)
end

"""
    openMapped(path, n; readonly=false) -> PackedTriangle

Memory-map an existing raw packed file (e.g. one created with
`PackedTriangle(n; mmap=true, path=path)`). With `readonly=true` writes throw.
"""
function openMapped(path::AbstractString, n::Integer; readonly::Bool=false)
    len = packedlength(n)
    filesize(path) == 4len || throw(ArgumentError("$path is $(filesize(path)) bytes; expected $(4len) for n = $n"))
    len == 0 && return PackedTriangle(n, Float32[], false)
    return PackedTriangle(n, mapfile(path, len; create=false, readonly), true)
end

"""
    Mmap.sync!(p::PackedTriangle)

Flush a memory-mapped packed triangle to disk (no-op for in-memory ones).
"""
function Mmap.sync!(p::PackedTriangle)
    p.mapped && Mmap.sync!(p.data)
    return p
end

"""
    packTriangle(mat; mmap=false, path=nothing) -> PackedTriangle

Pack the strict upper triangle of a square matrix. `mmap` and `path` are as in
[`PackedTriangle(n; mmap, path)`](@ref).
"""
function packTriangle(mat::AbstractMatrix{<:Real}; mmap::Bool=false, path::Union{Nothing,AbstractString}=nothing)
    n = size(mat, 1)
    size(mat, 2) == n || throw(DimensionMismatch("matrix must be square"))
    p = PackedTriangle(n; mmap, path)
    k = 0
    @inbounds for j in 2:n, i in 1:j-1
        k += 1
        p.data[k] = mat[i, j]
    end
    return p
end

"""
    unpack(p::PackedTriangle) -> Matrix{Float32}

Expand to a dense `n × n` matrix with zeros on and below the diagonal. This
allocates the full n² matrix, so it is meant for tests and sanity checks on
small inputs.
"""
function unpack(p::PackedTriangle)::Matrix{Float32}
    n = p.n
    mat = zeros(Float32, n, n)
    k = 0
    @inbounds for j in 2:n, i in 1:j-1
        k += 1
        mat[i, j] = p.data[k]
    end
    return mat
end

const PACKED_MAGIC = b"WRAITHPT"
const PACKED_VERSION = UInt8(1)
const PACKED_CHUNK = 1 << 20   # Float32 elements per (de)compression chunk (4 MB)

"""
    writePacked(path, p::PackedTriangle; level=6) -> path

Write `p` to `path` gzip-compressed. File layout: the 8-byte magic `WRAITHPT`,
a `UInt8` format version, the `Int64` matrix size `n`, then a gzip stream
holding the packed `Float32` values in the machine's (little-endian)
byte order. Data is compressed in 4 MB chunks, so memory use stays flat even
for memory-mapped triangles. `level` is the gzip level (1 fastest, 9 smallest).
"""
function writePacked(path::AbstractString, p::PackedTriangle; level::Integer=6)
    open(path, "w") do io
        write(io, PACKED_MAGIC)
        write(io, PACKED_VERSION)
        write(io, Int64(p.n))
        gz = GzipCompressorStream(io; level=level)
        try
            chunk = Vector{Float32}(undef, min(PACKED_CHUNK, length(p.data)))
            for start in 1:PACKED_CHUNK:length(p.data)
                len = min(PACKED_CHUNK, length(p.data) - start + 1)
                buf = len == length(chunk) ? chunk : Vector{Float32}(undef, len)
                copyto!(buf, 1, p.data, start, len)
                write(gz, buf)
            end
        finally
            close(gz)   # flushes the gzip trailer
        end
    end
    return path
end

"""
    readPacked(path; mmap=false, mmap_path=nothing) -> PackedTriangle

Read a file written by [`writePacked`](@ref). With `mmap=true` the decompressed
values are streamed straight into a memory-mapped file at `mmap_path` (see
[`PackedTriangle(n; mmap, path)`](@ref)), so the whole triangle never has to fit
in RAM.
"""
function readPacked(path::AbstractString; mmap::Bool=false, mmap_path::Union{Nothing,AbstractString}=nothing)
    open(path, "r") do io
        read(io, length(PACKED_MAGIC)) == PACKED_MAGIC || throw(ArgumentError("$path is not a WRAITHPT packed file"))
        version = read(io, UInt8)
        version == PACKED_VERSION || throw(ArgumentError("unsupported packed file version $version"))
        n = Int(read(io, Int64))
        p = PackedTriangle(n; mmap, path=mmap_path)
        gz = GzipDecompressorStream(io)
        chunk = Vector{Float32}(undef, min(PACKED_CHUNK, length(p.data)))
        for start in 1:PACKED_CHUNK:length(p.data)
            len = min(PACKED_CHUNK, length(p.data) - start + 1)
            buf = len == length(chunk) ? chunk : Vector{Float32}(undef, len)
            read!(gz, buf)   # throws EOFError on a truncated file
            copyto!(p.data, start, buf, 1, len)
        end
        eof(gz) || throw(ArgumentError("$path contains more data than n = $n requires"))
        return p
    end
end

"""
Z-scores for a [`PackedTriangle`](@ref), returned packed as well (`mmap` and
`path` as in [`PackedTriangle(n; mmap, path)`](@ref)).
"""
function jaccardScores(mat::PackedTriangle, stats::DiagonalStats=diagonalStats(mat);
    mmap::Bool=false, path::Union{Nothing,AbstractString}=nothing)::PackedTriangle
    σ = diagonalSD(stats)
    z = PackedTriangle(stats.n; mmap, path)
    k = 0
    @inbounds for j in 2:stats.n, i in 1:j-1
        k += 1
        d = j - i
        z.data[k] = zscore(mat.data[k], stats.μ[d], σ[d])
    end
    return z
end
