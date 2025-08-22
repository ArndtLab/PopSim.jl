
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

