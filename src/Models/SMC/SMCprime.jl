module SMCprimeapprox

using ..PopSim
using Distributions
using Random


import ..PopSim: sim_ancestry

export SMCprime

struct SMCprime <: AbstractEvolutionaryModel end

function sim_ancestry(
    model::SMCprime,
    demography::Demography,
    genome::Genome{R,M},
    pop_sample = nothing;
) where {R<:AbstractRateDistribution,M<:AbstractRateDistribution}

    @assert genome.recombination isa PopSim.UniformRate "SMCprime only implemented for uniform recombination rate"
    @assert genome.mutation isa PopSim.UniformRate "SMCprime only implemented for uniform mutation rate"

    tnv = TNvector(demography, length(genome))
    recombination_rate=PopSim.rate(recombination(genome))
    mutation_rate=PopSim.rate(mutation(genome))
    if length(tnv) == 2
        pop = PopSim.StationaryPopulation(;size = tnv[2], genome_length = tnv[1],
            mutation_rate, recombination_rate)
    else
        pop = PopSim.VaryingPopulation(; TNvector=tnv, mutation_rate, recombination_rate)
    end

    PopSim.SimulatedAncestry(model, demography, genome, pop_sample, pop)
end



mutable struct IBDIterator{T}
    pop::T
    tau_recombination::Int
    tau_previous::Int
    max_length::Int
end

IBDIterator(pop::StationaryPopulation, max_length = genome_length(pop)) = IBDIterator(pop, 1, 1, max_length)
IBDIterator(pop::VaryingPopulation{T}, max_length = genome_length(pop)) where {T} = IBDIterator(pop, 1, 1, max_length)


Base.IteratorSize(::Type{IBDIterator{T}}) where {T} = Base.SizeUnknown()
Base.IteratorEltype(::Type{IBDIterator{T}}) where {T} = Base.HasEltype()
Base.eltype(::Type{IBDIterator{T}}) where {T} = ARGsegment{Int64, CoalescentTreeTwoLineages}


function Base.iterate(ti::IBDIterator{StationaryPopulation}, pos = 1)
    if pos > ti.max_length
        return nothing
    end

    tau_back = rand(Geometric(1 / size(ti.pop)))
    if (ti.tau_recombination + tau_back) > ti.tau_previous
        tau = rand(Geometric(1 / (2 * size(ti.pop)))) + ti.tau_previous
    else
        equal = rand(Bool)
        if equal
            tau = ti.tau_previous
        else
            tau = ti.tau_recombination + tau_back
        end
    end
    ti.tau_previous = tau
    ti.tau_recombination = rand(DiscreteUniform(1, tau))
    recombination_rate_bp = 2 * recombination_rate(ti.pop) * tau
    if recombination_rate_bp < 1
        len = rand(Geometric(recombination_rate_bp))
    else
        len = 1
    end

    tree = CoalescentTreeTwoLineages(0, tau) # for two individuals no actual tree is returned

    stop = min(pos + len - 1, ti.max_length)
    return ARGsegment(Segment(pos, stop), tree), stop + 1
end


function find_last_epoch(t::T, times::Vector{S}) where {T, S}
    for i in length(times):-1:1
        if t >= times[i]
            return i
        end
    end
    return 1
end


function Base.iterate(ti::IBDIterator{VaryingPopulation{T}}, pos = 1) where {T}
    if pos > ti.max_length
        return nothing
    end

    (; tau_recombination, tau_previous) = ti
    (; population_sizes, times) = ti.pop

    epoch = find_last_epoch(tau_recombination, times)
    tau_back = tau_recombination + rand(Geometric(1 / population_sizes[epoch]))
    while (epoch < length(times)) && (tau_back >= times[epoch + 1]) && (tau_previous > times[epoch + 1])
        epoch += 1
        tau_back = times[epoch] + rand(Geometric(1 / population_sizes[epoch]))
    end
    if tau_back > tau_previous
        epoch_pr = find_last_epoch(tau_previous, times)
        tau = tau_previous + rand(Geometric(1 / (2 * population_sizes[epoch_pr])))
        while (epoch_pr < length(times)) && (tau >= times[epoch_pr + 1])
            epoch_pr += 1
            tau = times[epoch_pr] + rand(Geometric(1 / (2 * population_sizes[epoch_pr])))
        end
    else
        equal = rand(Bool)
        if equal
            tau = tau_previous
        else
            tau = tau_back
        end
    end
    ti.tau_previous = ceil(Int, tau)
    ti.tau_recombination = rand(DiscreteUniform(1, ti.tau_previous))
    recombination_rate_bp = 2 * recombination_rate(ti.pop) * tau
    if recombination_rate_bp < 1
        len = rand(Geometric(recombination_rate_bp))
    else
        len = 1
    end

    tree = CoalescentTreeTwoLineages(0, ti.tau_previous) # for two individuals no actual tree is returned

    stop = min(pos + len - 1, ti.max_length)
    return ARGsegment(Segment(pos, stop), tree), stop + 1
end


end # module SMCprimeapprox
