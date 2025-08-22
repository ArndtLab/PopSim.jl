
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
# Base.IteratorEltype(::Type{IBDIteratorMulti{T}}) where {T} = Base.HasEltype()
# Base.eltype(::Type{IBDIteratorMulti{T}}) where {T} = 
#     SegItem{Int64, CoalescentTree{Vector{@NamedTuple{pid::Int64, time::Float64, vid_anc::Int64}}}}

function Base.iterate(ti::IBDIteratorMulti)
    pos = 1
    iterate(ti, pos)
end


function Base.iterate(ti::IBDIteratorMulti, pos)
    if pos > ti.cos.genome_length
        return nothing
    end

    branches = map(1: 2 * length(ti.ids)-1) do i
        (pid = i <= length(ti.ids) ? ti.ids[i] : -1, 
            time = 0.0,
            vid_anc = -1)
    end
    nbranch = length(ti.ids) + 1
    tmp = map(enumerate(ti.ids)) do (i, p)
        (vid_child = i, pid_child = p, pid_anc = p, nb = ti.cos.genome_length)
    end

    while true
        sort!(tmp, by = x -> x.pid_anc, rev = true)

        tmp[1].pid_anc == 0 && break # no more coalescent events

        if tmp[1].pid_anc != tmp[2].pid_anc # move
            npid, nnb = getparentat(ti.cos, tmp[1].pid_anc, pos)
            nnb = min(tmp[1].nb, nnb)
            tmp[1] = (; tmp[1]..., pid_anc = npid, nb = nnb)
        else # merge
            i = 1
            while i+1 <= length(tmp)
                if tmp[i+1].pid_anc == tmp[1].pid_anc
                    i += 1
                else
                    break
                end
            end
            nb = minimum(getindex.(tmp[1:i], :nb))
            for j in 1:i
                branches[tmp[j].vid_child] = (pid = tmp[j].pid_child, time = gettime(ti.cos, tmp[j].pid_child), vid_anc = nbranch)
            end
            branches[nbranch] = (pid = tmp[1].pid_anc, time = gettime(ti.cos, tmp[1].pid_anc), vid_anc = -1)
            tmp[1] = ( vid_child = nbranch, pid_child = tmp[1].pid_anc, pid_anc = tmp[1].pid_anc, nb = nb)
            nbranch += 1
            deleteat!(tmp, 2:i)
        end

        length(tmp) == 1 && break # no more coalescent events
    end
    for i in 1:length(ti.ids)
        branches[i] = (; branches[i] ..., time = ti.end_time)
    end

    nb = minimum(getindex.(tmp, :nb))
    filter!(x -> x.pid >= 0, branches)

    pid = length(tmp) == 1 ? tmp[1].pid_child : 0    # pid == 0 signals a non-coalescent event
    time = pid > 0 ? gettime(ti.cos, pid) : -Inf

    data = CoalescentTree(ti.ids, pid, time, ti.end_time, branches)

    return ARGsegment(Segment(pos, nb), data), nb + 1
end

