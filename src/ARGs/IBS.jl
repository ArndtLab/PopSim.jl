
export IBSIterator, IBMIterator



mutable struct IBSIterator{T,RD <: AbstractRateDistribution} 
    ibxs::Iterators.Stateful{T}
    mutation::RD
    mut_kwargs
    breaks::Vector{Int64}
    nbi::Int64
    lastibxstop::Int64
end

IBSIterator(ibds, mutation; mut_kwargs...) = 
    IBSIterator(Iterators.Stateful(ibds), mutation, mut_kwargs, Int[], 0, 0)


Base.IteratorSize(::Type{IBSIterator{T, RD}}) where {T, RD} = Base.SizeUnknown()
Base.IteratorEltype(::Type{IBSIterator{T}}) where {T} = Base.HasEltype()
Base.eltype(::Type{IBSIterator{T, RD}}) where {T, RD} = ARGsegment{Int64, Int64}




function Base.iterate(si::IBSIterator, pos = 1)

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

            dt = timespan(ibx)
            si.breaks = APop.sample(si.mutation, 2 * dt, first(ibx), last(ibx); si.mut_kwargs...)

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


function sprinckle_mutations(s::ARGsegment{Int64, CoalescentTree{Vector{Branch}}}, mut::AbstractRateDistribution; kwargs...) 

    tree = data(s)
    branches = tree.branches
    mutbranches = map(branches) do b
        if b.ancestor_k > 0
            dt = b.time - branches[b.ancestor_k].time
            @assert dt >= 0.0 "dt=$dt < 0.0 for branch $(b.id) with ancestor $(b.ancestor_k)"
            MutatedBranch(b.id, b.time, b.ancestor_k,
                APop.sample(mut, dt, first(s), last(s); kwargs...)
            )
        else
            MutatedBranch(b.id, b.time, b.ancestor_k, Int64[])
        end
    end
    muttree = CoalescentTree(tree.ids, tree.root_id, tree.start_time, tree.end_time, mutbranches)

    ARGsegment(s.segment, muttree)
end


function IBMIterator(iter, mutation; kwargs...)
    Iterators.map(s->sprinckle_mutations(s, mutation; kwargs...),iter)
end


# -----------------------------------------------------------------------------


function collect_sprinckled_mutations(ct::CoalescentTree{Vector{MutatedBranch}}, id_ks::Vector{Int64})
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
    ibxs::T
    vids::V
    breaks::Vector{Int64}
    lastibxstop::Int64
end


IBSIteratorMutated(ibx, vid1::Int64, vid2::Int64) = 
    IBSIteratorMutated(ibx, [vid1, vid2],Int[], 0)


Base.IteratorSize(::Type{IBSIteratorMutated{T,V}}) where {T,V} = Base.SizeUnknown()
Base.IteratorEltype(::Type{IBSIteratorMutated{T,V}}) where {T,V} = Base.HasEltype()
Base.eltype(::Type{IBSIteratorMutated{T,V}}) where {T,V} = ARGsegment{Int64, Int64}


mutable struct IBSIteratorMutatedState{T}
    ibxstate::Union{Nothing,T}
    laststop::Int64
    b::Int64
end



function Base.iterate(si::IBSIteratorMutated)
    ibx = iterate(si.ibxs)
    isnothing(ibx) && return nothing # no ibx to iterate

    seg = ibx[1]
    @assert all(1 .<= si.vids .<= length(data(seg).ids))
    si.breaks = collect_sprinckled_mutations(data(seg), si.vids)
    si.lastibxstop = last(seg)
    state = IBSIteratorMutatedState(ibx[2], 0, 1)
    iterate(si, state)
end


function Base.iterate(si::IBSIteratorMutated, state)
    isnothing(state) && return nothing

    mystart = state.laststop + 1
    # @show mystart, state.b, si.breaks

    if state.b <= length(si.breaks)
        mystop = si.breaks[state.b]
        state.b += 1
        state.laststop = mystop
        return ARGsegment(Segment(mystart, mystop), 0), state
    else
        r = 0
        while true
            ibx = iterate(si.ibxs, state.ibxstate)
            if isnothing(ibx) # last ibd reached, emit last interval
                return ARGsegment(Segment(mystart, si.lastibxstop), r), nothing
            end

            state.ibxstate = ibx[2]
            seg = ibx[1]
            r += 1
            # println("seg: ", seg)

            si.breaks = collect_sprinckled_mutations(data(seg), si.vids)
            si.lastibxstop = last(seg)

            if isempty(si.breaks) # no breaks in ibx
                continue
            end

            # @show si.breaks
            mystop = si.breaks[1]
            state.b = 2
            state.laststop = mystop

            return ARGsegment(Segment(mystart, mystop), r), state
        end
    end
end




# IBSIterator(collection, args...; kwargs...) = _IBSIterator(collection, Base.IteratorEltype(collection), args...; kwargs...)
# _IBSIterator(collection, ::Base.EltypeUnknown, args...; kwargs...) = throw(ArgumentError("The eltype of the collection is unknown"))
# _IBSIterator(collection, ::Base.HasEltype, args...; kwargs...) =  _IBSIterator(collection, eltype(collection), args...; kwargs...)


# _IBSIterator(collection, ::Type{SegItem{Int64, PopSimBase.CoalescentTrees.SimpleCoalescentTree}}, args...; kwargs...)  = IBSIteratorNonMutated(collection, args...; kwargs...)
# _IBSIterator(collection, ::Type{SegItem{Int64, T}}, args...; kwargs...) where {T<:PopSimBase.CoalescentTrees.MutatedCoalescentTree} = IBSIteratorMutated(collection, args...; kwargs...)
# _IBSIterator(collection, ::Type{SegItem{Int64, T}}, args...; kwargs...) where {T<:PopSimBase.CoalescentTrees.MutatedSimpleCoalescentTree} = IBSIteratorMutated(collection, args...; kwargs...)

