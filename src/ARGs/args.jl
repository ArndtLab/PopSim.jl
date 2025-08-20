

export ARGsegment

mutable struct ARGsegment{T, Tr}
    segment::Segment{T}
    tree::Tr
end

DataStructures.@delegate ARGsegment.segment (Base.first, Base.last, Base.length)
DataStructures.@delegate ARGsegment.tree (timespan, iscoalescent)

# ARGsegment(segment::Segment{T}, tree::Tr) where {T, Tr} = ARGsegment{T, Tr}(segment, tree)


