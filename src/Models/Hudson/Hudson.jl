module HudsonModel

export Hudson, sim_ancestry

using ..APop

struct Hudson end


mutable struct StatefulWithDefaultIterator{T, D}
    v::Vector{T}
    k::Int
    default::D
end

StatefulWithDefaultIterator(v::Vector{T}, default::D) where {T, D} = 
    StatefulWithDefaultIterator(v, 1, default)

function StatefulWithDefaultIterator(v::Vector{Segment{T}}, default::D) where {T, D} 
    vt = map(v -> (first(v), last(v)), v)
    StatefulWithDefaultIterator(vt, 1, (default, default))
end

function StatefulWithDefaultIterator(v::Vector{ARGsegment{T, D}}, default::S) where {T, D, S} 
    vt = map(v -> (first(v), last(v), v.data), v)
    StatefulWithDefaultIterator(vt, 1, (default, default, nothing))
end


function nextitem(si::StatefulWithDefaultIterator)
    if si.k > length(si.v)
        return si.default
    else
        si.k += 1
        return si.v[si.k-1]
    end
end




struct HudsonARG{T}
    id::Int
    time::T
    nleaves::Int
    child1::HudsonARG{T}
    child2::HudsonARG{T}

    HudsonARG{T}(id::Int, time::T) where {T} = new(id, time, 1)

    function HudsonARG{T}(id::Int, time::T, child1::HudsonARG{T}, child2::HudsonARG{T}) where {T}
        new(id, time, child1.nleaves + child2.nleaves, child1, child2)
    end
end

Base.show(io::IO, h::HudsonARG) = print(io, "HudsonARG(id=$(h.id), time=$(h.time), nleaves=$(h.nleaves))")




function distribute(vi::Vector{Segment{T}}, bps::Vector{T}) where {T<:Integer}
    v1 = Vector{Segment{T}}()
    v2 = Vector{Segment{T}}()
    length(vi) == 0 && return (v1, v2)
    length(bps) == 0 && return (vi, v2)

    posmax = last(vi[end]) + 1
    vis = StatefulWithDefaultIterator(vi, posmax)
    nextpushv1 = true

    bpsi = StatefulWithDefaultIterator(bps, posmax)
    nextbppos = nextitem(bpsi)

    nextintervalstart, nextintervalstop = nextitem(vis)

    pos = min(nextintervalstop, nextbppos)
    k = 1
    while pos < posmax
        # @show (k, pos) (nextintervalstart, nextintervalstop, nextbppos)

        if (pos == nextintervalstop || pos == nextbppos) && nextintervalstart <= pos
            if nextpushv1
                push!(v1, Segment(nextintervalstart, pos))
            else
                push!(v2, Segment(nextintervalstart, pos))
            end
            nextintervalstart = pos + 1
        end
        if pos == nextintervalstop
            nextintervalstart, nextintervalstop = nextitem(vis)
        end
        if pos == nextbppos
            nextbppos = nextitem(bpsi)
            nextpushv1 = !nextpushv1
        end
        # @show (nextintervalstart, nextintervalstop, nextbppos)

        # @show v1 v2

        pos = min(nextintervalstop, nextbppos)
        k += 1
        # k > 20 && break
    end

    v1, v2
end



function distribute(vi::Vector{ARGsegment{T, D}}, bps::Vector{T}) where {D, T<:Integer}
    v1 = similar(vi, 0)
    v2 = similar(vi, 0)
    length(vi) == 0 && return (v1, v2)
    length(bps) == 0 && return (vi, v2)

    posmax = last(vi[end]) + 1
    vis = StatefulWithDefaultIterator(vi, posmax)
    nextpushv1 = true

    bpsi = StatefulWithDefaultIterator(bps, posmax)
    nextbppos = nextitem(bpsi)

    nextintervalstart, nextintervalstop, nextintervaldata = nextitem(vis)

    pos = min(nextintervalstop, nextbppos)
    k = 1
    while pos < posmax

        if (pos == nextintervalstop || pos == nextbppos) && nextintervalstart <= pos
            if nextpushv1
                push!(v1, ARGsegment(Segment(nextintervalstart, pos), nextintervaldata))
            else
                push!(v2, ARGsegment(Segment(nextintervalstart, pos), nextintervaldata))
            end
            nextintervalstart = pos + 1
        end
        if pos == nextintervalstop
            nextintervalstart, nextintervalstop, nextintervaldata = nextitem(vis)
        end
        if pos == nextbppos
            nextbppos = nextitem(bpsi)
            nextpushv1 = !nextpushv1
        end

        pos = min(nextintervalstop, nextbppos)
        k += 1
    end

    v1, v2
end



function distribute(vi::Vector, mutation::AbstractRateDistribution)

    length(vi) == 0 && return (similar(vi, 0), similar(vi, 0))

    bps = APop.sample(mutation, 1.0, first(vi[1]), last(vi[end]))
    length(bps) == 0 && return (vi, similar(vi, 0))
   
    distribute(vi, bps)
end




function coalesce(
    v1::Vector{Segment{T}},
    v2::Vector{Segment{T}},
    vc::Vector{ARGsegment{T, CoalescentTreeTwoLineages}},
    t::F, tmax::F,
    n::Int
) where {T<:Integer,F<:Real}
    
    @assert n == 2
    delta_t = tmax - t
    length(v1) == 0 && return v2
    length(v2) == 0 && return v1

    v = Vector{Segment{T}}()

    posmax = max(last(v1[end]), last(v2[end])) + 1

    v1s = StatefulWithDefaultIterator(v1, posmax)
    v2s = StatefulWithDefaultIterator(v2, posmax)
    next1start, next1stop = nextitem(v1s)
    next2start, next2stop = nextitem(v2s)

    pos = min(next1start, next2start)

    while pos < posmax

        if pos == next1start && next2start == pos
            nextpos = min(next1stop, next2stop)
        elseif pos == next1start && next2start > pos
            nextpos = min(next1stop, next2start)
        elseif pos == next2start && next1start > pos
            nextpos = min(next2stop, next1start)
        elseif pos == next1start && next2start < pos
            push!(v, Segment(next2start, pos - 1))
            next2start = pos
            nextpos = min(next1stop, next2stop)
        elseif pos == next2start && next1start < pos
            push!(v, Segment(next1start, pos - 1))
            next1start = pos
            nextpos = min(next1stop, next2stop)
        end

        if pos == next1stop && pos == next2stop
            @assert next1start == next2start
            tree = CoalescentTreeTwoLineages(0, delta_t)
            push!(vc, ARGsegment(Segment(next1start, pos), tree))
            next1start, next1stop = nextitem(v1s)
            next2start, next2stop = nextitem(v2s)
            nextpos = min(next1start, next2start)
        elseif pos == next1stop && next2start > pos
            push!(v, Segment(next1start, pos))
            next1start, next1stop = nextitem(v1s)
            nextpos = min(next1start, next2start)
        elseif pos == next1stop && next2stop > pos
            @assert next1start == next2start
            tree = CoalescentTreeTwoLineages(0, delta_t)
            push!(vc, ARGsegment(Segment(next1start, pos), tree))
            next1start, next1stop = nextitem(v1s)
            next2start = pos + 1
            nextpos = min(next1start, next2stop)
        elseif pos == next2stop && next1start > pos
            push!(v, Segment(next2start, pos))
            next2start, next2stop = nextitem(v2s)
            nextpos = min(next1start, next2start)
        elseif pos == next2stop && next1stop > pos
            @assert next1start == next2start
            tree = CoalescentTreeTwoLineages(0, delta_t)
            push!(vc, ARGsegment(Segment(next1start, pos), tree))
            next2start, next2stop = nextitem(v2s)
            next1start = pos + 1
            nextpos = min(next1stop, next2start)
        end

        pos = nextpos
    end

    return v
end


function coalesce(
    v1::Vector{ARGsegment{T, HudsonARG{F}}} ,
    v2::Vector{ARGsegment{T, HudsonARG{F}}} ,
    vc::Vector{ARGsegment{T, HudsonARG{F}}},
    t::F, tmax::F,
    n::Int
) where {T<:Integer ,F<:Real}


    length(v1) == 0 && return v2
    length(v2) == 0 && return v1

    v = similar(v1, 0)

    posmax = max(last(v1[end]), last(v2[end])) + 1

    v1s = StatefulWithDefaultIterator(v1, posmax)
    v2s = StatefulWithDefaultIterator(v2, posmax)
    next1start, next1stop, next1data = nextitem(v1s)
    next2start, next2stop, next2data = nextitem(v2s)

    pos = min(next1start, next2start)

    while pos < posmax

        if pos == next1start && next2start == pos
            nextpos = min(next1stop, next2stop)
        elseif pos == next1start && next2start > pos
            nextpos = min(next1stop, next2start)
        elseif pos == next2start && next1start > pos
            nextpos = min(next2stop, next1start)
        elseif pos == next1start && next2start < pos
            push!(v, ARGsegment(Segment(next2start, pos - 1), next2data))
            next2start = pos
            nextpos = min(next1stop, next2stop)
        elseif pos == next2start && next1start < pos
            push!(v, ARGsegment(Segment(next1start, pos - 1), next1data))
            next1start = pos
            nextpos = min(next1stop, next2stop)
        end

        if pos == next1stop && pos == next2stop
            @assert next1start == next2start
            tree = HudsonARG{T}(0, t, next1data, next2data)
            if tree.nleaves == n
                push!(vc, ARGsegment(Segment(next1start, pos), tree))
            else
                push!(v, ARGsegment(Segment(next1start, pos), tree))
            end
            next1start, next1stop, next1data = nextitem(v1s)
            next2start, next2stop, next2data = nextitem(v2s)
            nextpos = min(next1start, next2start)
        elseif pos == next1stop && next2start > pos
            push!(v, ARGsegment(Segment(next1start, pos), next1data))
            next1start, next1stop, next1data = nextitem(v1s)
            nextpos = min(next1start, next2start)
        elseif pos == next1stop && next2stop > pos
            @assert next1start == next2start
            tree = HudsonARG{T}(0, t, next1data, next2data)
            if tree.nleaves == n
                push!(vc, ARGsegment(Segment(next1start, pos), tree))
            else
                push!(v, ARGsegment(Segment(next1start, pos), tree))
            end
            next1start, next1stop, next1data = nextitem(v1s)
            next2start = pos + 1
            nextpos = min(next1start, next2stop)
        elseif pos == next2stop && next1start > pos
            push!(v, ARGsegment(Segment(next2start, pos), next2data))
            next2start, next2stop, next2data = nextitem(v2s)
            nextpos = min(next1start, next2start)
        elseif pos == next2stop && next1stop > pos
            @assert next1start == next2start
            tree = HudsonARG{T}(0, t, next1data, next2data)
            if tree.nleaves == n
                push!(vc, ARGsegment(Segment(next1start, pos), tree))
            else
                push!(v, ARGsegment(Segment(next1start, pos), tree))
            end
            next2start, next2stop, next2data = nextitem(v2s)
            next1start = pos + 1
            nextpos = min(next1stop, next2start)
        end
        
        pos = nextpos
    end
    
    return v
end



struct SimulatedAncestry{M <: Hudson, T}
    model::M
    demography::Demography
    genome::Genome
    sample::Dict{String, Int}
    treedata::T
end

function APop.sim_ancestry(model::Hudson, demography::Demography, genome::Genome,
    sample::Int)
    
    if length(demography.populations) != 1
        throw(ArgumentError("You need to specify sample as a Dict{String, Int}"))
    end
    if sample < 2
        throw(ArgumentError("Sample size must be at least 2"))
    end
    sample = Dict(demography.populations[1].id => sample)
    APop.sim_ancestry(model, demography, genome, sample)
end

function APop.sim_ancestry(model::Hudson, demography::Demography, genome::Genome,
    sample::Dict{String, Int})
    
    
    @assert length(demography.populations) == 1 "not implemented"

    tmax = demography.end_time
    tmin = demography.start_time
    L = length(genome)

    n = sum(values(sample))

    # initial segments
    if n == 2
        v1 = [[Segment{Int}(1, L)], [Segment{Int}(1, L)]]
        v2 = similar(v1, 0)
        vc = Vector{ARGsegment{Int, CoalescentTreeTwoLineages}}()
    else
        v1 = map(1:n) do i
            [ARGsegment(Segment{Int}(1, L), HudsonARG{Int}(i, tmax))]
        end
        v2 = similar(v1, 0)
        vc = similar(v1[1], 0)
    end



    t = tmax - 1
    while true
        empty!(v2)

        N = demography.population_sizes[1][t < tmin ? tmin : t]
        for vi in v1
            vi1, vi2 = distribute(vi, recombination(genome))

            if !isempty(vi1)
                k = rand(1:2*N)
                if k > length(v2)
                    push!(v2, vi1)
                else
                    v2[k] = coalesce(v2[k], vi1, vc, t, tmax, n)
                end
            end

            if !isempty(vi2)
                k = rand(1:2*N)
                if k > length(v2)
                    push!(v2, vi2)
                else
                    v2[k] = coalesce(v2[k], vi2, vc, t, tmax, n)
                end
            end
        end

        v1 = filter(x -> length(x) > 0, v2)
        if isempty(v1)
            break
        end

        t -= 1

    end
    sort!(vc, by = first)
    @assert sum(length, vc) == L

    SimulatedAncestry(model, demography, genome, sample, vc)
end




function createbranches!(arg, branches, nextinternal, idc, parent_idc, last_time)
    branches[idc] = Branch(idc, arg.time, parent_idc)
    if arg.nleaves == 1
        return
    else 
        if arg.child1.nleaves == 1
            branches[arg.child1.id] = Branch(arg.child1.id, last_time, idc)
        else
            nextinternal -= 1
            createbranches!(arg.child1, branches, nextinternal, nextinternal + 1, idc, last_time)
        end
        if arg.child2.nleaves == 1
            branches[arg.child2.id] = Branch(arg.child2.id, last_time, idc)
        else
            nextinternal -= 1
            createbranches!(arg.child2, branches, nextinternal, nextinternal + 1, idc, last_time)
        end
        return
    end
end

function get_ARGsegments(sa::SimulatedAncestry{M, Vector{ARGsegment{Int, CoalescentTreeTwoLineages}}}) where {M<:Hudson}
    return sa.treedata
end

function get_ARGsegments(sa::SimulatedAncestry{M, Vector{ARGsegment{Int, HudsonARG{F}}}}) where {M<:Hudson, F}
    n = sum(values(sa.sample))

    map(sa.treedata) do vci
        ids = collect(1:n)
        first_id = 2 * n - 1 
        idc = first_id
        first_time = vci.data.time
        last_time = sa.demography.end_time
        branches = map(i -> Branch(-1, 0.0, -1), 1:idc)

        nextinternal = idc - 1
        parent_idc = -1
        createbranches!(vci.data, branches, nextinternal, idc, parent_idc, last_time)

        tree = CoalescentTree(ids, first_id, first_time, last_time, branches)
        ARGsegment(vci.segment, tree)
    end
end

end # module Hudson