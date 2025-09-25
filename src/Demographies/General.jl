
# ============================================================
# Population
# ============================================================

export Population

const TimeType = Int64

struct Population
    id::String
    description::String
    size::Int
    growth_rate::Float64
    time_offset::TimeType
end



function Population(;id::String = "", 
        description::String = "", 
        size::Int = 0, 
        growth_rate::Float64 = 0.0, 
        time_offset::TimeType = 0
        )
    Population(id, description, size, growth_rate, time_offset)
end




Base.size(p::Population) = p.size

growth_rate(p::Population) = p.growth_rate

time_offset(p::Population) = p.time_offset



# ============================================================
# Event
# ============================================================

export AbstractAncestryEvent, PopulationSizeEvent, PopulationSplitEvent, PopulationMergeEvent

abstract type AbstractAncestryEvent end

struct PopulationSizeEvent <: AbstractAncestryEvent
    time::TimeType
    population_id::AbstractString
    value
end

struct PopulationSplitEvent <: AbstractAncestryEvent
    time::TimeType
    source_population_id::AbstractString
    target_population_ids::Vector{AbstractString}
end

PopulationSplitEvent(time::TimeType, source_population_id::AbstractString, 
    target_population_id1::AbstractString, target_population_id2::AbstractString) =
    PopulationSplitEvent(time, source_population_id, [target_population_id1, target_population_id2])

struct PopulationMergeEvent <: AbstractAncestryEvent
    time::TimeType
    source_population_ids::Vector{AbstractString}
    target_population_id::AbstractString
end

PopulationMergeEvent(time::TimeType, source_population_id1::AbstractString, 
    source_population_id2::AbstractString, target_population_id::AbstractString) =
    PopulationMergeEvent(time, [source_population_id1, source_population_id2], target_population_id)

export Demography, add_population!, add_event!, set_start_time!, set_end_time!, set_migration!,
    get_population_index_by_id, fix_population_sizes!, summary, TNvector, get_migration_parent_pop_sampler, test_population_sizes

using OffsetArrays

mutable struct Demography
    populations::Vector{Population}
    events::Vector{AbstractAncestryEvent}
    
    start_time::TimeType
    end_time::TimeType
    ploidy::Int64
    
    population_sizes::Vector{OffsetVector{Int64}}
    migration::Array{Float64, 2}
end



Demography() = Demography(Vector{Population}(), Vector{AbstractAncestryEvent}(), 0.0, 0.0, 2, Vector{OffsetVector{Int64}}(), zeros(Float64, 0, 0))


function add_population!(d::Demography, population::Population)
    if any(p -> p.id == population.id, d.populations)
        throw(ArgumentError("Population with id $(population.id) already exists"))
    end
    push!(d.populations, population)

    # Update migration matrix
    np = length(d.populations)
    m = zeros(Float64, np, np)
    for i in 1:np-1, j in 1:np-1
        m[i, j] = d.migration[i, j]  # copy existing migration rates
    end
    m[np, np] = 1.0  # self-migration rate
    d.migration = m 

    fix_population_sizes!(d)
end

function add_event!(d::Demography, event::AbstractAncestryEvent)
    push!(d.events, event)
    fix_population_sizes!(d)
end

function set_start_time!(d::Demography, time::TimeType)
    d.start_time = time
    fix_population_sizes!(d)
end

function set_end_time!(d::Demography, time::TimeType)
    d.end_time = time
    fix_population_sizes!(d)
end

function set_via_TNvector!(d::Demography, tnv::Vector)
    add_population!(d, Population(id="pop", size=round(Int, tnv[2])))
    t = 0
    set_end_time!(d, t)
    for i in length(tnv):-2:3
        t -= round(Int, tnv[i-1])
        add_event!(d, PopulationSizeEvent(t, "pop", round(Int, tnv[i])))
        set_start_time!(d, t-1)
    end
    fix_population_sizes!(d)
    return d
end


function set_migration!(d::Demography, fromid::AbstractString, toid::AbstractString, rate::Float64)
    @assert rate >= 0.0 "Migration rate must be non-negative"
    @assert fromid != toid "Cannot set migration rate from a population to itself"

    fromindex = get_population_index_by_id(d, fromid)
    toindex = get_population_index_by_id(d, toid)
    d.migration[toindex, fromindex] = rate
    s = sum(d.migration[toindex, 1:toindex-1]) + sum(d.migration[toindex, toindex+1:end])
    d.migration[toindex, toindex] = 1.0 - s
    @assert d.migration[toindex, toindex] >= 0.0 "Self-migration rate must be non-negative"
    d
end

function get_population_index_by_id(d::Demography, id::AbstractString)
    for (i, p) in enumerate(d.populations)
        if p.id == id
            return i
        end
    end
    throw(ArgumentError("Population with id $id not found"))
    return nothing
end


function fix_population_sizes!(d::Demography)
    sort!(d.events, by = e -> e.time)
    nextevent = 1


    T = d.end_time - d.start_time + 1
    T <= 0 && return d

    d.population_sizes = map(p -> OffsetVector(fill(p.size, T), d.start_time : d.end_time), d.populations)
    @assert all(p -> p.growth_rate == 0.0, d.populations)

    for t in d.start_time+1 : d.end_time
        for i in 1:length(d.populations)
            d.population_sizes[i][t] = d.population_sizes[i][t-1]
        end

        while nextevent <= length(d.events) && d.events[nextevent].time <= t
            e = d.events[nextevent]
            nextevent += 1
            if e isa PopulationSizeEvent
                i = get_population_index_by_id(d, e.population_id)
                d.population_sizes[i][t] = e.value
            elseif e isa PopulationSplitEvent
                si = get_population_index_by_id(d, e.source_population_id)
                d.population_sizes[si][t : end] .= 0
                for target_population_id in e.target_population_ids
                    ti = get_population_index_by_id(d, target_population_id)
                    d.population_sizes[ti][begin : t-1] .= 0
                    @assert d.population_sizes[ti][t] > 0
                end
            elseif e isa PopulationMergeEvent
                for source_population_id in e.source_population_ids
                    si = get_population_index_by_id(d, source_population_id)
                    d.population_sizes[si][t : end] .= 0
                end
                ti = get_population_index_by_id(d, e.target_population_id)
                d.population_sizes[ti][begin : t-1] .= 0
                @assert d.population_sizes[ti][t] > 0
            else
                throw(ArgumentError("Unknown event type: $(typeof(e))"))
            end
        end
    end
    return d
end

function summary(d::Demography)
    io = IOBuffer()
    write(io, "Demography with $(length(d.populations)) populations, " *
           "$(length(d.events)) events, start time $(d.start_time), end time $(d.end_time)")

    ts = map(enumerate(d.start_time : d.end_time)) do (i, t)
        p = map(p -> lpad(string(d.population_sizes[p][t]), 7), 1:length(d.populations))
        (i, t, join(p, " "))
    end
    ks = map(t->t[1], unique(i -> i[3], ts))
    push!(ks, length(ts))

    lastprintedk = 1
    for k in ks
        if k > 1
            (ip, tp, sp) = ts[k-1]
            if ip > lastprintedk + 1
                write(io, "\n ...")
            end
            write(io, "\n", lpad(string(tp), 7), ": $sp")
        end
        write(io, "\n", lpad(string(ts[k][2]), 7), ": ", ts[k][3])
        lastprintedk = k
    end
        
    return String(take!(io))
end



function TNvector(d::Demography, sequence_length::Int) 
    if length(d.populations) != 1
        throw(ArgumentError("TNvector is only defined for demographies with a single population"))
    end
    sizes = d.population_sizes[1]
    t = tlast = d.start_time
    Nlast = sizes[t]
    TN = [sequence_length, Nlast]

    while t <= d.end_time
        if sizes[t] != Nlast || t == d.end_time
            if tlast != d.start_time
                push!(TN, t - tlast)
                push!(TN, sizes[tlast])
            end
            tlast = t
            Nlast = sizes[t]
        end
        t += 1
    end

    return TN
end



function get_rand_parentpool(demography::Demography, pop_id::Int64)
    if demography.migration[pop_id, pop_id] == 1.0
        return pop_id
    else
        dist = Categorical([demography.migration[pop_id, from_id] for from_id in 1:length(demography.populations)])
        return rand(dist)
    end
end

function test_population_sizes(d::Demography)
    for ns in d.population_sizes
        length(ns) == 0 && continue
        i1 = findfirst(!=(0), ns)
        i2 = findlast(!=(0), ns)
        if isnothing(i1)
            @warn "Population sizes are all zero."
        end
        if !all(ns[i1:i2] .> 0)
            @warn "Population sizes reach zero."
        end
    end
end

