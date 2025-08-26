
function get_ARGsegments(sa::SimulatedAncestry, ids::Vector{Int64}) 
    if length(ids) == 2
        IBDIterator(sa.cos, ids, float(sa.demography.end_time))
    else
        IBDIteratorMulti(sa.cos, ids, float(sa.demography.end_time))
    end
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




function Base.iterate(ti::IBDIterator, pos::Int = 1)
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



# ---------------------------------



struct IBDIteratorMulti{C<: AbstractCrossoverStore}
    cos::C
    ids::Vector{Int64}
    end_time::Float64

    function IBDIteratorMulti{C}(cos, ids::Vector{Int64}, end_time::Float64) where{C <: AbstractCrossoverStore} 
        length(ids) >= 2 || throw(ArgumentError("The ids must be a vector of length >= 2"))
        new(cos, ids, end_time)
    end
end

IBDIteratorMulti(cos::C, ids::Vector{Int64}, end_time::Float64 = NaN) where {C <: AbstractCrossoverStore} = 
    IBDIteratorMulti{C}(cos, ids, end_time)




Base.IteratorSize(::Type{IBDIteratorMulti{T}}) where {T} = Base.SizeUnknown()
Base.IteratorEltype(::Type{IBDIteratorMulti{T}}) where {T} = Base.HasEltype()
Base.eltype(::Type{IBDIteratorMulti{T}}) where {T} = ARGsegment{Int64, CoalescentTree{Vector{Branch}}}




mutable struct BranchToMerge
    child_id::Int
    child_k::Int
    ancestor_id::Int
    nb::Int
end



function Base.iterate(ti::IBDIteratorMulti, pos::Int = 1)
    if pos > ti.cos.genome_length
        return nothing
    end

    branches = map(ti.ids) do id
        Branch(id, gettime(ti.cos, id), -1)
    end


    tmp = map(enumerate(ti.ids)) do (i, p)
        BranchToMerge(p, i, p, ti.cos.genome_length)
    end

    while true
        sort!(tmp, by = x -> x.ancestor_id, rev = true)

        tmp[1].ancestor_id == 0 && break # no more coalescent events

        if tmp[1].ancestor_id != tmp[2].ancestor_id 
            # move first backwards in time
            npid, nnb = getparentat(ti.cos, tmp[1].ancestor_id, pos)
            nnb = min(tmp[1].nb, nnb)
            tmp[1].ancestor_id = npid
            tmp[1].nb = nnb
        else 
            # coalesce at least two lineages
            i = 1
            while i+1 <= length(tmp)
                if tmp[i+1].ancestor_id == tmp[1].ancestor_id
                    i += 1
                else
                    break
                end
            end
            nb = minimum(t -> t.nb, tmp[1:i])

            for k in 1:i
                branches[tmp[k].child_k].ancestor_k = length(branches) + 1
            end
            push!(branches, Branch(tmp[1].ancestor_id, gettime(ti.cos, tmp[1].ancestor_id), -1))

            tmp[1].child_id = tmp[1].ancestor_id
            tmp[1].child_k = length(branches)
            tmp[1].nb = nb
            deleteat!(tmp, 2:i)
        end

        length(tmp) == 1 && break # no more coalescent events
    end

    nb = minimum(t -> t.nb, tmp)


    root_id = length(tmp) == 1 ? tmp[1].child_id : 0    # root_id == 0 signals a non-coalescent event
    root_time = root_id > 0 ? gettime(ti.cos, root_id) : -Inf
    data = CoalescentTree(ti.ids, root_id, root_time, ti.end_time, branches)

    return ARGsegment(Segment(pos, nb), data), nb + 1
end

