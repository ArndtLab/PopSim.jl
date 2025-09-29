
export IBSIterator, IBMIterator, IBD2Iterator



mutable struct IBSIteratorTwoLineages{T,RD <: AbstractRateDistribution} 
    ibxs::T
    mutation::RD
    multiple_hits::Symbol
    breaks::Vector{Int64}
    nbi::Int64
    lastibxstop::Int64
end

IBSIteratorTwoLineages(ibds, mutation::AbstractRateDistribution; multiple_hits::Symbol = :ignore) = 
    IBSIteratorTwoLineages(Iterators.Stateful(ibds), mutation, multiple_hits, Int[], 0, 0)

IBSIteratorTwoLineages(ibds, mutation_rate::Number; mut_kwargs...) = 
    IBSIteratorTwoLineages(ibds, UniformRate(mutation_rate); mut_kwargs...)

Base.IteratorSize(::Type{IBSIteratorTwoLineages{T, RD}}) where {T, RD} = Base.SizeUnknown()
Base.IteratorEltype(::Type{IBSIteratorTwoLineages{T}}) where {T} = Base.HasEltype()
Base.eltype(::Type{IBSIteratorTwoLineages{T, RD}}) where {T, RD} = ARGsegment{Int64, Int64}




function Base.iterate(si::IBSIteratorTwoLineages, pos = 1)

    isnothing(pos) && return nothing # no position to iterate

    if si.nbi > 0 && si.nbi <= length(si.breaks)
        mystart = pos
        mystop = si.breaks[si.nbi]
        si.nbi += 1
        return ARGsegment(Segment(mystart, mystop), 0), mystop + 1
    else
        r = 0
        while true
            if Iterators.isdone(si.ibxs) # last ibd reached, emit last interval
                return  ARGsegment(Segment(pos, si.lastibxstop), r), nothing
            end
            
            ibx = popfirst!(si.ibxs)
            r += 1
            si.lastibxstop = last(ibx)

            dt = timespan(ibx)
            si.breaks = PopSim.sample(si.mutation, 2 * dt, first(ibx), last(ibx); multiple_hits = si.multiple_hits)

            if !isempty(si.breaks)
                si.nbi = 1
                mystart = pos
                mystop = si.breaks[si.nbi]
                si.nbi += 1
                return ARGsegment(Segment(mystart, mystop), r), mystop + 1
            end
        end
    end
end








# -----------------------------------------------------------------------------


function sprinckle_mutations(s::ARGsegment{Int64, CoalescentTree{Vector{Branch},F}}, mut::AbstractRateDistribution; kwargs...)::ARGsegment{Int64, CoalescentTree{Vector{MutatedBranch},F}} where {F}

    tree = data(s)
    branches = tree.branches
    mutbranches = map(branches) do b
        if b.ancestor_k > 0
            dt = b.time > -Inf ? b.time - branches[b.ancestor_k].time : 0.0
            @assert dt >= 0.0 "dt=$dt < 0.0 for branch $(b.id) with ancestor $(b.ancestor_k): time1: $(b.time) time2: $(branches[b.ancestor_k].time)"
            MutatedBranch(b.id, b.time, b.ancestor_k,
                PopSim.sample(mut, dt, first(s), last(s); kwargs...)
            )
        else
            MutatedBranch(b.id, b.time, b.ancestor_k, Int64[])
        end
    end
    muttree = CoalescentTree(tree.ids, tree.root_id, tree.start_time, tree.end_time, mutbranches)

    ARGsegment(s.segment, muttree)::ARGsegment{Int64, CoalescentTree{Vector{MutatedBranch},F}}
end


function IBMIterator(iter, mutation; kwargs...)
    Iterators.map(s->sprinckle_mutations(s, mutation; kwargs...),iter)
end



# -----------------------------------------------------------------------------


function collect_sprinckled_mutations(ct::CoalescentTree{Vector{MutatedBranch},F}, id_ks::Vector{Int64}) where {F}
    breaks = Int64[]
    ks = sort(id_ks)
    while true
        length(ks) == 1 && break

        if (length(ks) >= 2) && (ks[1] == ks[2])
            ks = ks[2:end]
            continue 
        end
        
        # move along branch to ancestor
        ancestor_k = ct.branches[ks[1]].ancestor_k

        if ancestor_k > 0
            append!(breaks, ct.branches[ks[1]].mutations)
            ks[1] = ancestor_k
        else
            ks = ks[2:end]
        end
        sort!(ks)
    end
    sort!(breaks)
end







mutable struct IBSIteratorMutated{T, V}
    ibxs::Iterators.Stateful{T}
    vids::V
    breaks::Vector{Int64}
    nbi::Int64
    lastibxstop::Int64
end


IBSIteratorMutated(ibx, vids::Vector{Int64}) = 
    IBSIteratorMutated(Iterators.Stateful(ibx), vids, Int[], 0, 0)

IBSIteratorMutated(ibx, vid1::Int64, vid2::Int64) = 
    IBSIteratorMutated(ibx, [vid1, vid2])


Base.IteratorSize(::Type{IBSIteratorMutated{T,V}}) where {T,V} = Base.SizeUnknown()
Base.IteratorEltype(::Type{IBSIteratorMutated{T,V}}) where {T,V} = Base.HasEltype()
Base.eltype(::Type{IBSIteratorMutated{T,V}}) where {T,V} = ARGsegment{Int64, Int64}



function Base.iterate(si::IBSIteratorMutated, pos = 1)

    isnothing(pos) && return nothing # no position to iterate

    if si.nbi > 0 && si.nbi <= length(si.breaks)
        mystart = pos
        mystop = si.breaks[si.nbi]
        si.nbi += 1
        return ARGsegment(Segment(mystart, mystop), 0), mystop + 1
    else
        r = 0
        while true
            if Iterators.isdone(si.ibxs) # last ibd reached, emit last interval
                return  ARGsegment(Segment(pos, si.lastibxstop), r), nothing
            end
            
            ibx, _ = iterate(si.ibxs)
            r += 1
            si.lastibxstop = last(ibx)
            si.breaks = collect_sprinckled_mutations(data(ibx), si.vids)

            if !isempty(si.breaks)
                si.nbi = 1
                mystart = pos
                mystop = si.breaks[si.nbi]
                si.nbi += 1
                return ARGsegment(Segment(mystart, mystop), r), mystop + 1
            end
        end
    end
end




IBSIterator(collection, args...; kwargs...) = _IBSIterator(collection, Base.IteratorEltype(collection), args...; kwargs...)
_IBSIterator(collection, ::Base.EltypeUnknown, args...; kwargs...) = throw(ArgumentError("The eltype of the collection is unknown"))
_IBSIterator(collection, ::Base.HasEltype, args...; kwargs...) =  _IBSIterator(collection, eltype(collection), args...; kwargs...)


_IBSIterator(collection, ::Type{ARGsegment{Int64, CoalescentTreeTwoLineages}}, args...; kwargs...) = IBSIteratorTwoLineages(collection, args...; kwargs...)
_IBSIterator(collection, ::Type{ARGsegment{Int64, CoalescentTree{Vector{MutatedBranch}, F}}}, args...; kwargs...) where {F} = IBSIteratorMutated(collection, args...; kwargs...)





function find_last_common_ancestor(ct::CoalescentTree{Vector{Branch},F}, id_ks::Vector{Int64}) where {F}
    ks = sort(id_ks)
    while true
        length(ks) == 1 && break

        if (length(ks) >= 2) && (ks[1] == ks[2])
            ks = ks[2:end]
            continue 
        end
        
        # move along branch to ancestor
        ancestor_k = ct.branches[ks[1]].ancestor_k

        if ancestor_k > 0
            ks[1] = ancestor_k
        else
            ks = ks[2:end]
        end
        sort!(ks)
    end
    ks[1]
end




function IBD2Iterator(iter, ids::Vector{Int64})
    Iterators.map(iter) do a
        s = data(a)
        lca = find_last_common_ancestor(s, ids)
        ARGsegment{Int64, CoalescentTreeTwoLineages}(a.segment, CoalescentTreeTwoLineages(s.root_id, end_time(s) - s.branches[lca].time))
    end
end



