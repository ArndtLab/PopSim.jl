module SMCapprox

using Distributions
using Random
using ..APop

export SMC

struct SMC <: AbstractEvolutionaryModel end



mutable struct IBDIterator
    pop::StationaryPopulation
    tau_recombination::Int
end

IBDIterator(pop::StationaryPopulation) = IBDIterator(pop, 1)


Base.IteratorSize(::Type{IBDIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{IBDIterator}) = Base.HasEltype()
Base.eltype(::Type{IBDIterator}) = ARGsegment{Int64, CoalescentTreeTwoLineages}


function Base.iterate(ti::IBDIterator, pos = 1)
    if pos > genome_length(ti.pop)
        return nothing
    end

    tau = rand(Geometric(1 / (2 * size(ti.pop)))) + ti.tau_recombination
    ti.tau_recombination = rand(DiscreteUniform(1, tau))
    recombination_rate_bp = 2 * recombination_rate(ti.pop) * tau
    if recombination_rate_bp < 1
        len = rand(Geometric(recombination_rate_bp))
    else
        len = 1
    end

    tree = CoalescentTreeTwoLineages(0, tau) # for two individuals no actual tree is returned

    stop = min(pos + len - 1, genome_length(ti.pop))
    return ARGsegment(Segment(pos, stop), tree), stop + 1
end



end # module SMCapprox