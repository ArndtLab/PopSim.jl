
module WrightFisherForwardModel

export WrightFisher

using ..PopSim
import ..PopSim: sim_ancestry

include("CrossOverStores.jl")
using .CrossoverStores









struct WrightFisher <: AbstractEvolutionaryModel
    allow_selfing::Bool
end

WrightFisher(;allow_selfing::Bool = true) = WrightFisher(allow_selfing)

randomswap(a, b) = rand(Bernoulli()) ? (a,b) : (b,a)





# ============================================================
# Individual
# ============================================================

using StaticArrays
using DataStructures
using Distributions


export Individual, randallele, ploidy

struct Individual{P} 
    alleles::SVector{P, Int64}
end

Individual(alleles...) = Individual(SVector{length(alleles)}(alleles...))
Individual(alleles::Vector) = Individual(SVector{length(alleles)}(alleles))


Base.show(io::IO, ind::Individual{P}) where P = print(io, "Individual (ploidy: $P, $(ind.alleles))")

DataStructures.@delegate Individual.alleles (Base.getindex, Base.setindex!, Base.length)

ploidy(::Individual{P}) where P = P

randallele(ind::Individual) = rand(ind.alleles)







function sim_ancestry(model::WrightFisher, demography::Demography, genome::Genome, popsample = nothing)
    # prechecks
    @assert demography.ploidy == 2 "Not implemented"

    # println(PopSim.summary(demography))


    # setup
    PopSim.fix_population_sizes!(demography)


    cos = CrossoverStores.MemoryCrossoverStore(length(genome))
    t = demography.start_time
    
    # initialize individuals
    alives = map(enumerate(demography.populations)) do (i, pop)
        N = demography.population_sizes[i, t]
        alive = map(1:N) do i
            ids = map(1:demography.ploidy) do j
                CrossoverStores.newid!(cos, float(t))
            end    
            Individual(ids...)
        end
        alive
    end    
    
    # run loop 
    nextevent = 1
    for t in demography.start_time + 1 : demography.end_time

        # handle events
        while nextevent <= length(demography.events) && demography.events[nextevent].time <= t
            e = demography.events[nextevent]
            nextevent += 1
            if e isa PopulationSizeEvent
                continue  # already taken care of in fix_population_sizes!
                throw(ArgumentError("Not implemented event type: $(typeof(e))"))
            elseif e isa PopulationSplitEvent
                si = get_population_index_by_id(demography, e.source_population_id)
                for target_population_id in e.target_population_ids
                    ti = get_population_index_by_id(demography, target_population_id)
                    alives[ti] = copy(alives[si])
                end
                alives[si] = Vector{Individual}()
            elseif e isa PopulationMergeEvent
                ti = get_population_index_by_id(demography, e.target_population_id)
                alives[ti] = Vector{Individual}()
                for source_population_id in e.source_population_ids
                    si = get_population_index_by_id(demography, source_population_id)
                    append!(alives[ti], alives[si])
                    alives[si] = Vector{Individual}()
                end
            else
                throw(ArgumentError("Unknown event type: $(typeof(e))"))
            end
        end

        # reproduction
        nalives = map(enumerate(alives)) do (i, alive)
            targetN = demography.population_sizes[i, t]
            alive1 = map(1:targetN) do k
                parentpool = alives[PopSim.get_rand_parentpool(demography, i)]
                a1, a2 = sample(parentpool, 2, replace=model.allow_selfing)
                
                a11, a12 = randomswap(a1[1], a1[2])
                breaks1 = sample(genome.recombination, 1.0, 1, length(genome))

                a21, a22 = randomswap(a2[1], a2[2])
                breaks2 = sample(genome.recombination, 1.0, 1, length(genome))

                Individual(CrossoverStores.newid!(cos, a11, a12, float(t), breaks1), CrossoverStores.newid!(cos, a21, a22, float(t), breaks2))
            end
            alive1
        end
        
        alives .= nalives

    end

    # sample individuals
    if isnothing(popsample)
        cos_ids = mapreduce(vcat, alives) do alive
            mapreduce(vcat, alive, init = Vector{Int}()) do indv
                indv.alleles
            end
        end
    end
    if popsample isa Int
        popsample = Sample(demography, popsample)
    end
    if popsample isa Sample
        cos_ids = Int[]
        for id in popsample.ids
            pop_id = get_population_index_by_id(demography, id)
            if isnothing(pop_id) || length(alives[pop_id]) == 0
                throw(ArgumentError("Population name '$id' not found in demography populations or population has no individuals at the end of the simulation."))
            end
            n = 10000
            while n > 0
                indv = rand(alives[pop_id])
                cos_id = randallele(indv)
                if !in(cos_id, cos_ids)
                    push!(cos_ids, cos_id)
                    break
                end
                n -= 1
            end
            n == 0 && throw(ArgumentError("Could not find a sample for "))
        end
    end

    return PopSim.SimulatedAncestry(model, demography, genome, popsample, (;cos, alives, cos_ids))
end


include("arg.jl")

end