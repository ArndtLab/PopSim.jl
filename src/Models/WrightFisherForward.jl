

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

    # setup
    fix_population_sizes!(demography)
    
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
    for t in demography.start_time : demography.end_time
        # reproduction
        nalives = map(enumerate(alives)) do (i, alive)
            targetN = demography.population_sizes[i][t]
            parentpool = alive
            
            alive1 = map(1:targetN) do i
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

        # handle events
        while nextevent <= length(demography.events) && demography.events[nextevent].time == t
            e = demography.events[nextevent]
            nextevent += 1
            if e isa ParameterChangeEvent
                throw(ArgumentError("Not implemented event type: $(typeof(e))"))
            elseif e isa PopulationSplitEvent
                si = get_population_index_by_id(demography, e.source_population_id)
                t1i = get_population_index_by_id(demography, e.target_population_id1)
                t2i = get_population_index_by_id(demography, e.target_population_id2)
                throw(ArgumentError("Not implemented event type: $(typeof(e))"))

            elseif e isa PopulationMergeEvent
                si1 = get_population_index_by_id(demography, e.source_population_id1)
                si2 = get_population_index_by_id(demography, e.source_population_id2)
                ti = get_population_index_by_id(demography, e.target_population_id)
                throw(ArgumentError("Not implemented event type: $(typeof(e))"))

            else
                throw(ArgumentError("Unknown event type: $(typeof(e))"))
            end
        end
    end

    return SimulatedAncestry{typeof(model)}(model, demography, genome, cos, alives)
end


function get_ARGsegments(sa::SimulatedAncestry, ids::Vector{Int64}) 
    @assert length(ids) == 2 "Only two lineages are supported for now"
    IBDIterator(sa.cos, ids, float(sa.demography.end_time))
end


struct IBDIterator{C}
    cos::C
    ids::Vector{Int64}
    end_time::Float64

    function IBDIterator{C}(cos, ids::Vector{Int64}, end_time::Float64) where{C <: AbstractCrossoverStore} 
        length(ids) == 2 || throw(ArgumentError("The ids must be a vector of length 2"))
        new(cos, ids, end_time)
    end
end

IBDIterator(cos::AbstractCrossoverStore, ids::Vector{Int64}, end_time::Float64 = NaN) = 
    IBDIterator{typeof(cos)}(cos, ids, end_time)




Base.IteratorSize(::Type{IBDIterator{T}}) where {T} = Base.SizeUnknown()
Base.IteratorEltype(::Type{IBDIterator{T}}) where {T} = Base.HasEltype()
Base.eltype(::Type{IBDIterator{T}}) where {T} = ARGsegment{Int64, CoalescentTreeTwoLineages}




function Base.iterate(ti::IBDIterator)
    pos = 1
    iterate(ti, pos)
end


function Base.iterate(ti::IBDIterator, pos)
    if pos > ti.cos.genome_length
        return nothing
    end
    pid1, nb1 = ti.ids[1], ti.cos.genome_length
    pid2, nb2 = ti.ids[2], ti.cos.genome_length

    while ((pid1 > 0) || (pid2 > 0)) && (pid1 != pid2)
        if pid1 < pid2
            pid2, nb = CrossoverStores.getparentat(ti.cos, pid2, pos)
            nb2 = min(nb2, nb)
        else
            pid1, nb = CrossoverStores.getparentat(ti.cos, pid1, pos)
            nb1 = min(nb1, nb)
        end
    end

    pid = pid1 == pid2 ? pid1 : 0    # pid == 0 signals a non-coalescent event
    time = pid > 0 ? gettime(ti.cos, pid) : -Inf


    tree = CoalescentTreeTwoLineages(pid, ti.end_time - time)

    nb = min(nb1, nb2)
    return ARGsegment(Segment(pos, nb), tree), nb + 1
end

