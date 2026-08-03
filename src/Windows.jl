
"""
Struct that stores a `StepRange`. This is for an implementation of `iterate()` and
`next!()` such that we can use this in `for` loops, enumeration, zip, etc. the
struct accepts a start value for the range, a stop value, and a step, just like 
`range(start, stop, step=...)` (a.k.a) `start:step:stop`. The additional things here
make sure that a UnitRange is returned for use in `eachoverlap` and since StepRanges must have
equal-length elements, that the final genomic window, which is usually shorter than the others,
is included as well. The range in `GenomicWindows` is stateful, meaning that it gets consumed
upon use.
"""
mutable struct GenomicWindows
    range::Base.Iterators.Stateful{StepRange{Int64,Int64},Union{Nothing,Tuple{Int64,Int64}}}
    current::Int
    stop::Int
end

"""
Instantiate a `GenomicWindows` object from a `start`, `stop`, and `step`
like you would a UnitRange. This is intended to be used with a known chromosome
length or series of positions. e.g.

    contig_len = 432234426
    gws = GenomicWindows(1, contig_len, 50_000)
"""
function GenomicWindows(start::Int, stop::Int, step::Int)::GenomicWindows
    stateful = Iterators.Stateful(start:step:stop)
    effective_start = first(iterate(stateful)) - 1
    GenomicWindows(stateful, effective_start, stop)
end

"""
This is an interface for manual iteration, mutating the state
of the underlying `GenomicWindows.range` and updating `GenomicWindows.current`.
It returns a `UnitRange` of x:y where `x` is the value of the range
when this function was invoked and `y` is the value of the range after iterating
the range once (y-1 becoming the new `GenomicWindows.current` value). 
"""
function next!(gws::GenomicWindows)::UnitRange{Int64}
    if gws.current == gws.stop
        error("Nothing left to iterate")
    end
    start = gws.current + 1
    state = iterate(gws.range)
    if isnothing(state)
        gws.current = gws.stop
        return start:gws.stop
    end
    gws.current = first(state) - 1
    return start:gws.current
end

function Base.iterate(gws::GenomicWindows, _=nothing)
    if gws.current == gws.stop
        return nothing
    end
    return (next!(gws), nothing)
end
