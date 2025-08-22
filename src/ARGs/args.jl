

export ARGsegment

mutable struct ARGsegment{T, D}
    segment::Segment{T}
    data::D
end

data(as::ARGsegment) = as.data

DataStructures.@delegate ARGsegment.segment (Base.first, Base.last, Base.length)


timespan(as::ARGsegment{T,D}) where {T,D <: AbstractCoalescentTree} = timespan(as.data)
iscoalescent(as::ARGsegment{T,D}) where {T,D <: AbstractCoalescentTree} = iscoalescent(as.data)


