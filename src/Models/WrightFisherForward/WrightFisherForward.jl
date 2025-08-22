

abstract type AbstractEvolutionaryModel end


export WrightFisher, AbstractEvolutionaryModel, sim_ancestry, get_ARGsegments




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

function get_migration_parent_pop_sampler(demography::Demography)
    P = length(demography.populations)
    migration_parent_pop_sampler = map(1:P) do to_id
        if demography.migration[to_id, to_id] == 1.0
            () -> to_id
        else
            dist = Categorical([demography.migration[to_id, from_id] for from_id in 1:P])
            () -> rand(dist)
        end
    end
    return migration_parent_pop_sampler
end



using .CrossoverStores



struct SimulatedAncestry{M <: WrightFisher}
    model::M
    demography::Demography
    genome::Genome
    cos::CrossoverStores.AbstractCrossoverStore
    alives::Vector{Vector{Individual}}
end






function sim_ancestry(model::WrightFisher, demography::Demography, genome::Genome)
    # prechecks
    @assert demography.ploidy == 2 "Not implemented"

    # println(APop.summary(demography))


    # setup
    fix_population_sizes!(demography)
    P = length(demography.populations)

    migration_parent_pop_sampler = get_migration_parent_pop_sampler(demography)

    cos = CrossoverStores.MemoryCrossoverStore(length(genome))
    t = demography.start_time
    
    # initialize individuals
    alives = map(demography.population_sizes) do sizes
        N = sizes[t]
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
                t1i = get_population_index_by_id(demography, e.target_population_id1)
                t2i = get_population_index_by_id(demography, e.target_population_id2)
                alives[t1i] = alives[si]
                alives[t2i] = copy(alives[si])
                alives[si] = Vector{Individual}()
            elseif e isa PopulationMergeEvent
                si1 = get_population_index_by_id(demography, e.source_population_id1)
                si2 = get_population_index_by_id(demography, e.source_population_id2)
                ti = get_population_index_by_id(demography, e.target_population_id)
                alives[ti] = vcat(alives[si1], alives[si2])
                alives[si1] = Vector{Individual}()
                alives[si2] = Vector{Individual}()
            else
                throw(ArgumentError("Unknown event type: $(typeof(e))"))
            end
        end

        # reproduction
        nalives = map(enumerate(alives)) do (i, alive)
            targetN = demography.population_sizes[i][t]
            alive1 = map(1:targetN) do k
                parentpool = alives[migration_parent_pop_sampler[i]()]
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

    return SimulatedAncestry{typeof(model)}(model, demography, genome, cos, alives)
end

