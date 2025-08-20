
export AbstractCoalescentTree, AbstractMutatedCoalescentTree,
    CoalescentTreeTwoLineages,
    timespan, iscoalescent



abstract type AbstractCoalescentTree end
abstract type AbstractMutatedCoalescentTree end


struct CoalescentTreeTwoLineages <: AbstractCoalescentTree
    root_id::Int64
    timespan::Float64
end


timespan(tree::CoalescentTreeTwoLineages) = tree.timespan
iscoalescent(tree::CoalescentTreeTwoLineages) = tree.timespan >= 0.0


