
# ============================================================
# Population
# ============================================================

export Population, setparam!

const TimeType = Int64

mutable struct Population
    id::String
    description::String
    size::Int
    growth_rate::Float64
    time0::TimeType
end



function Population(;id::String = "", 
        description::String = "", 
        size::Int = 0, 
        growth_rate::Float64 = 0.0, 
        time0::TimeType = 0
        )
    Population(id, description, size, growth_rate, time0)
end


function setparam!(p::Population, param::Symbol, value)
    if param == :size
        p.size = value
    elseif param == :growth_rate
        p.growth_rate = value
    elseif param == :time0
        p.time0 = value
    else
        throw(ArgumentError("Unknown or immutable parameter $param"))
    end
end

Base.size(p::Population) = p.size

growth_rate(p::Population) = p.growth_rate

time0(p::Population) = p.time0




# ============================================================
# Event
# ============================================================

export AbstractAncestryEvent, ParameterChangeEvent, PopulationSplitEvent, PopulationMergeEvent

abstract type AbstractAncestryEvent end

struct ParameterChangeEvent <: AbstractAncestryEvent
    time::TimeType
    population_id::AbstractString
    parameter::Symbol
    value
end

struct PopulationSplitEvent <: AbstractAncestryEvent
    time::TimeType
    source_population_id::AbstractString

    target_population_id1::AbstractString
    target_size1::Int64

    target_population_id2::AbstractString
    target_size2::Int64
end

struct PopulationMergeEvent <: AbstractAncestryEvent
    time::TimeType
    source_population_id1::AbstractString
    source_population_id2::AbstractString
    target_population_id::AbstractString
    target_size::Int64
end



export Demography, add_population!, add_event!, set_start_time!, set_end_time!, set_migration!

using OffsetArrays

mutable struct Demography
    populations::Vector{Population}
    events::Vector{AbstractAncestryEvent}
    
    start_time::TimeType
    end_time::TimeType
    
    population_sizes::Vector{OffsetVector{Int64}}
    migration::Array{Float64, 2}
end



Demography() = Demography(Vector{Population}(), Vector{AbstractAncestryEvent}(), 0.0, 0.0, Vector{OffsetVector{Int64}}(), zeros(Float64, 0, 0))


function add_population!(d::Demography, population::Population)
    if any(p -> p.id == population.id, d.populations)
        throw(ArgumentError("Population with id $(population.id) already exists"))
    end
    push!(d.populations, population)
    if time0(population) < d.start_time
        d.start_time = time0(population)
    end
    m = zeros(Float64, length(d.populations), length(d.populations))
    for i in 1:length(d.populations)-1, j in 1:length(d.populations)-1
        m[i, j] = d.migration[i, j]  # copy existing migration rates
    end
    m[length(d.populations), length(d.populations)] = 1.0  # self-migration rate
    d.migration = m 
    fix_population_sizes!(d)
end

function add_event!(d::Demography, event::AbstractAncestryEvent)
    push!(d.events, event)
    if event.time < d.start_time
        d.start_time = event.time
    end
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

    d.population_sizes = map(p -> OffsetVector(zeros(Int64, T), d.start_time:d.end_time), d.populations)
    for t in d.start_time:d.end_time
        for (i, p) in enumerate(d.populations)
            @assert p.growth_rate == 0.0
            d.population_sizes[i][t] = p.size
        end
        while nextevent <= length(d.events) && d.events[nextevent].time == t
            e = d.events[nextevent]
            nextevent += 1
            if e isa ParameterChangeEvent
                i = get_population_index_by_id(d, e.population_id)
                setparam!(d.populations[i], e.parameter, e.value)
            elseif e isa PopulationSplitEvent
                si = get_population_index_by_id(d, e.source_population_id)
                t1i = get_population_index_by_id(d, e.target_population_id1)
                t2i = get_population_index_by_id(d, e.target_population_id2)

                setparam!(d.populations[si], :size, 0)
                setparam!(d.populations[t1i], :size, e.target_size1)
                setparam!(d.populations[t2i], :size, e.target_size2)
            elseif e isa PopulationMergeEvent
                si1 = get_population_index_by_id(d, e.source_population_id1)
                si2 = get_population_index_by_id(d, e.source_population_id2)
                ti = get_population_index_by_id(d, e.target_population_id)
                setparam!(d.populations[si1], :size, 0)
                setparam!(d.populations[si2], :size, 0)
                setparam!(d.populations[ti], :size, e.target_size)
            else
                throw(ArgumentError("Unknown event type: $(typeof(e))"))
            end
        end
    end
    return d
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

