
export AbstractCoalescentTree, AbstractMutatedCoalescentTree,
    CoalescentTreeTwoLineages, CoalescentTree, Branch, MutatedBranch,
    timespan, iscoalescent



abstract type AbstractCoalescentTree end
abstract type AbstractMutatedCoalescentTree end


struct CoalescentTreeTwoLineages <: AbstractCoalescentTree
    root_id::Int64
    timespan::Float64
end


timespan(tree::CoalescentTreeTwoLineages) = tree.timespan
iscoalescent(tree::CoalescentTreeTwoLineages) = tree.timespan >= 0.0



mutable struct Branch
    id::Int64
    time::Float64
    ancestor_k::Int64
end

mutable struct MutatedBranch
    id::Int64
    time::Float64
    ancestor_k::Int64
    mutations::Vector{Int64}
end


struct CoalescentTree{T} <: AbstractCoalescentTree
    ids::Vector{Int64}
    root_id::Int64
    start_time::Float64
    end_time::Float64
    branches::T
end

CoalescentTree(ids::Vector{Int64}, root_id::Int64, start_time::Float64, end_time::Float64) = CoalescentTree{Nothing}(ids, root_id, start_time, end_time, nothing)

start_time(ct::CoalescentTree{T}) where {T} = ct.start_time
end_time(ct::CoalescentTree{T}) where {T} = ct.end_time
timespan(ct::CoalescentTree{T}) where {T} = end_time(ct) - start_time(ct)
root_id(ct::CoalescentTree{T}) where {T}   = ct.root_id
iscoalescent(ct::CoalescentTree{T}) where {T} = root_id(ct) > 0

Base.show(io::IO, ct::CoalescentTree{T}) where {T} = print(io, "CoalescentTree starting at $(ct.start_time) in $(ct.root_id) for $(length(ct.ids)) individuals")

