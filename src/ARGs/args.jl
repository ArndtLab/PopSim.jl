

export ARG

mutable struct ARG{T}
    segment::Segment
    tree::T
end

DataStructures.@delegate ARG.segment (Base.first, Base.last, Base.length)
DataStructures.@delegate ARG.tree (timespan, iscoalescent)

ARG(segment::Segment, tree::T) where T = ARG{T}(segment, tree)


