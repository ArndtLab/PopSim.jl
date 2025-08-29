
export AbstractCoalescentTree,
    CoalescentTreeTwoLineages, CoalescentTree, Branch, MutatedBranch,
    timespan, iscoalescent



abstract type AbstractCoalescentTree end


struct CoalescentTreeTwoLineages <: AbstractCoalescentTree
    root_id::Int64
    timespan::Float64
end


timespan(tree::CoalescentTreeTwoLineages) = tree.timespan
iscoalescent(tree::CoalescentTreeTwoLineages) = tree.timespan >= 0.0

Base.show(io::IO, ct::CoalescentTreeTwoLineages) = print(io, "CoalescentTreeTwoLineages starting in $(ct.root_id) with timespan $(ct.timespan)")


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


struct CoalescentTree{T, F} <: AbstractCoalescentTree
    ids::Vector{Int64}
    root_id::Int64
    start_time::F
    end_time::F
    branches::T
end


start_time(ct::CoalescentTree{T, F}) where {T, F} = ct.start_time
end_time(ct::CoalescentTree{T, F}) where {T, F} = ct.end_time
timespan(ct::CoalescentTree{T, F}) where {T, F} = end_time(ct) - start_time(ct)
root_id(ct::CoalescentTree{T, F}) where {T, F}   = ct.root_id
iscoalescent(ct::CoalescentTree{T, F}) where {T, F} = root_id(ct) > 0

Base.show(io::IO, ct::CoalescentTree{T, F}) where {T, F} = print(io, "CoalescentTree starting at $(ct.start_time) in $(ct.root_id) for $(length(ct.ids)) individuals")

