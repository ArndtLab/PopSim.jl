
export AbstractCoalescentTree, AbstractMutatedCoalescentTree,
    CoalescentTreeTwoLineages, CoalescentTree,
    timespan, iscoalescent



abstract type AbstractCoalescentTree end
abstract type AbstractMutatedCoalescentTree end


struct CoalescentTreeTwoLineages <: AbstractCoalescentTree
    root_id::Int64
    timespan::Float64
end


timespan(tree::CoalescentTreeTwoLineages) = tree.timespan
iscoalescent(tree::CoalescentTreeTwoLineages) = tree.timespan >= 0.0



struct CoalescentTree{T} <: AbstractCoalescentTree
    ids::Vector{Int64}
    root_id::Int64
    start_time::Float64
    end_time::Float64
    tree::T
end

CoalescentTree(ids::Vector{Int64}, root_id::Int64, start_time::Float64, end_time::Float64) = CoalescentTree{Nothing}(ids, root_id, start_time, end_time, nothing)

start_time(ct::CoalescentTree{T}) where {T} = ct.start_time
end_time(ct::CoalescentTree{T}) where {T} = ct.end_time
timespan(ct::CoalescentTree{T}) where {T} = end_time(ct) - start_time(ct)
root_id(ct::CoalescentTree{T}) where {T}   = ct.root_id
iscoalescent(ct::CoalescentTree{T}) where {T} = root_id(ct) > 0

Base.show(io::IO, ct::CoalescentTree{T}) where {T} = print(io, "CoalescentTree starting at $(ct.start_time) in $(ct.root_id) for $(length(ct.ids)) individuals")

