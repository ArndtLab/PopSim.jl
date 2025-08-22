
# ============================================================
# Population
# ============================================================

export Population, setparam!

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

time0(p::Population) = p.time0




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

    target_population_id1::AbstractString
    target_population_id2::AbstractString
end

struct PopulationMergeEvent <: AbstractAncestryEvent
    time::TimeType
    source_population_id1::AbstractString
    source_population_id2::AbstractString
    target_population_id::AbstractString
end



export Demography, add_population!, add_event!, set_start_time!, set_end_time!, set_migration!

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
                t1i = get_population_index_by_id(d, e.target_population_id1)
                t2i = get_population_index_by_id(d, e.target_population_id2)
                d.population_sizes[si][t : end] .= 0
                d.population_sizes[t1i][begin : t-1] .= 0
                d.population_sizes[t2i][begin : t-1] .= 0
                @assert d.population_sizes[t1i][t] > 0
                @assert d.population_sizes[t2i][t] > 0
            elseif e isa PopulationMergeEvent
                si1 = get_population_index_by_id(d, e.source_population_id1)
                si2 = get_population_index_by_id(d, e.source_population_id2)
                ti = get_population_index_by_id(d, e.target_population_id)
                d.population_sizes[ti][begin : t-1] .= 0
                d.population_sizes[si1][t : end] .= 0
                d.population_sizes[si2][t : end] .= 0
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
    
    ts = map(d.start_time : d.end_time) do t
        p = map(p -> lpad(string(d.population_sizes[p][t]), 7), 1:length(d.populations))
        (t, join(p, " "))
    end

    for (t, s) in ts
        write(io, "\n", lpad(string(t), 7), ": $s")
    end

    
    return String(take!(io))
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

