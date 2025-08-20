

export ARGsegment

mutable struct ARGsegment{T}
    segment::Segment
    tree::T
end

DataStructures.@delegate ARGsegment.segment (Base.first, Base.last, Base.length)
DataStructures.@delegate ARGsegment.tree (timespan, iscoalescent)

ARGsegment(segment::Segment, tree::T) where T = ARGsegment{T}(segment, tree)


