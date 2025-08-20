



mutable struct IBSIteratorNonMutated{T,RD <: AbstractRateDistribution} 
    ibxs::T
    breaks::Vector{Int64}
    lastibxstop::Int64
    mutation::RD
    kwargs
end

IBSIteratorNonMutated(ibx, mutation; kwargs...) = 
    IBSIteratorNonMutated(ibx, Int[], 0, mutation, kwargs)


Base.IteratorSize(::Type{IBSIteratorNonMutated{T, RD}}) where {T, RD} = Base.SizeUnknown()
# Base.IteratorEltype(::Type{IBSIteratorNonMutated{T}}) where {T} = Base.HasEltype()
# Base.eltype(::Type{IBSIteratorNonMutated{T}}) where {T} = Segments.SegItem{Int64, Int64}


mutable struct IBSIteratorNonMutatedState{T, S}
    ibxstate::Union{Nothing,T}
    tree::S
    laststop::Int64
    b::Int64
end



function Base.iterate(si::IBSIteratorNonMutated)
    ibx = iterate(si.ibxs)
    isnothing(ibx) && return nothing # no ibx to iterate

    seg = ibx[1]

    dt = timespan(seg)
    si.breaks = APop.sample(si.mutation, 2 * dt, first(seg), last(seg); si.kwargs...)

    si.lastibxstop = last(seg)
    state = IBSIteratorNonMutatedState(ibx[2], seg.tree, 0, 1)
    iterate(si, state)
end


function Base.iterate(si::IBSIteratorNonMutated, state)
    isnothing(state) && return nothing

    mystart = state.laststop + 1
    # @show mystart, state.b, si.breaks

    if state.b <= length(si.breaks)
        mystop = si.breaks[state.b]
        state.b += 1
        state.laststop = mystop
        return ARGsegment(Segment(mystart, mystop), state.tree), state
    else
        r = 0
        while true
            ibx = iterate(si.ibxs, state.ibxstate)
            if isnothing(ibx) # last ibd reached, emit last interval
                return ARGsegment(Segment(mystart, si.lastibxstop), state.tree), nothing
            end

            state.ibxstate = ibx[2]
            seg = ibx[1]
            r += 1
            # println("seg: ", seg)

            dt = timespan(seg)
            si.breaks = APop.sample(si.mutation, 2 * dt, first(seg), last(seg); si.kwargs...)

            si.lastibxstop = last(seg)

            if isempty(si.breaks) # no breaks in ibx
                continue
            end

            # @show si.breaks
            mystop = si.breaks[1]
            state.b = 2
            state.tree = seg.tree
            state.laststop = mystop

            return ARGsegment(Segment(mystart, mystop), state.tree), state
        end
    end
end

