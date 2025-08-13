

export Segment



"""
    Segment(start::Int, stop::Int)

A genomic segment defined by its start and stop positions.
The start and end positions are inclusive.
"""
struct Segment
    first::Int
    last::Int
end

Base.first(s::Segment) = s.first
Base.last(s::Segment) = s.last
Base.length(s::Segment) = s.last - s.first + 1

Base.show(io::IO, s::Segment) = print(io, "Segment(first=$(s.first), last=$(s.last))")

