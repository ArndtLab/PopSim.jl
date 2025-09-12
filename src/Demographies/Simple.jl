


export AbstractDemography, StationaryPopulation, VaryingPopulation,
    genome_length, 
    ploidy, recombination_rate, mutation_rate,
    TNvector


abstract type AbstractDemography end

"""
    StationaryPopulation(ploidy::Int, size::Int)

    A simple demography for population modeling.
    A stationary population with a fixed ploidy and size.
"""
struct StationaryPopulation <: AbstractDemography
    ploidy::Int
    population_size::Int
    genome_length::Int
    recombination_rate::Float64
    mutation_rate::Float64


    function StationaryPopulation(ploidy::Int, 
        size::Int, 
        genome_length::Int, 
        recombination_rate::Float64, 
        mutation_rate::Float64 )
        if ploidy < 1 || ploidy > 2
            throw(ArgumentError("Ploidy must be 1 or 2."))
        end
        if size < 1
            throw(ArgumentError("Size must be a positive integer."))
        end
        return new(ploidy, size, genome_length, recombination_rate, mutation_rate)
    end
end




StationaryPopulation(;ploidy::Int = 2, 
    size::Int = 1000, 
    genome_length::Int = 1000000, 
    recombination_rate::Float64 = 1e-8, 
    mutation_rate::Float64 = 1e-8) = StationaryPopulation(ploidy, size, genome_length, recombination_rate, mutation_rate)


ploidy(d::StationaryPopulation) = d.ploidy
Base.size(d::StationaryPopulation) = d.population_size
genome_length(d::StationaryPopulation) = d.genome_length
recombination_rate(d::StationaryPopulation) = d.recombination_rate
mutation_rate(d::StationaryPopulation) = d.mutation_rate

Base.show(io::IO, d::StationaryPopulation) = print(io, "StationaryPopulation(ploidy=$(d.ploidy), population_size=$(d.population_size))")


TNvector(d::StationaryPopulation) = [genome_length(d), size(d)]





struct VaryingPopulation
    ploidy::Int
    genome_length::Int
    recombination_rate::Float64
    mutation_rate::Float64

    population_sizes::Vector{Int}
    times::Vector{Float64}
end


"""
    VaryingPopulation(; ploidy = 2, population_sizes = [1_000], times = [0.0],
    genome_length = 1_000_000, recombination_rate = 1.0e-8, mutation_rate = 1.0e-8)
Create a varying population with the given parameters.
# Arguments
- `ploidy::Int`: The ploidy of the population. Default is 2.
- `genome_length::Int`: The length of the genome. Default is 1_000_000.
- `recombination_rate::Float64`: The recombination rate. Default is 1.0e-8.
- `mutation_rate::Float64`: The mutation rate. Default is 1.0e-8.
- `population_sizes::Vector{Int}`: The sizes of the population at different times. Default is [1_000].
- `times::Vector{Float64}`: The times at which the population sizes are specified. Default is [0.0].

# Notes
- The population sizes and times must have the same length.
- The times must be sorted in ascending order.
- The first time must be 0.0.

# Example
To model a population which starts at N=1000, goes throw a bottleneck (N=200) 
from 150 genertions ago to 120 generations ago and then recovers to N=2000 individuals, 
you can use the following:

```julia
times = [0.0, 120.0, 150.0]
population_sizes = [2000, 200, 1000]
```
"""
function VaryingPopulation(;
        ploidy = 2,
        population_sizes = nothing,
        times = nothing,
        TNvector = nothing,
        genome_length = nothing,
        recombination_rate = 1.0e-8,
        mutation_rate = 1.0e-8
    )
    if !isnothing(TNvector)
        @assert isnothing(population_sizes) "Cannot provide both TNvector and population_sizes"
        @assert isnothing(times) "Cannot provide both TNvector and times"
        @assert isnothing(genome_length) "Cannot provide both TNvector and genome_length"

        t = 0
        times = [0.0]
        population_sizes = Int[]
        i = length(TNvector)
        while i > 2
            push!(population_sizes, TNvector[i])
            deltat = TNvector[i-1]
            t += deltat
            push!(times, t)
            i -= 2
        end
        @assert i == 2 "TNvector must have even length"
        push!(population_sizes, TNvector[2])
        genome_length = TNvector[1]
    end
    if isnothing(TNvector)
        @assert !isnothing(population_sizes) "Must provide either TNvector or population_sizes"
        @assert !isnothing(times) "Must provide either TNvector or times"
        @assert !isnothing(genome_length) "Must provide either TNvector or genome_length"
    end


    if length(population_sizes) != length(times)
        throw(ArgumentError("population_sizes and times must have the same length"))
    end
    if !issorted(times)
        throw(ArgumentError("times must be sorted"))
    end
    if times[1] != 0.0
        throw(ArgumentError("times must start at 0.0"))
    end
    VaryingPopulation(ploidy, genome_length, recombination_rate, mutation_rate, population_sizes, times)
end





genome_length(pop::VaryingPopulation) = pop.genome_length
recombination_rate(pop::VaryingPopulation) = pop.recombination_rate
mutation_rate(pop::VaryingPopulation) = pop.mutation_rate
ploidy(pop::VaryingPopulation) = pop.ploidy

function TNvector(pop::VaryingPopulation)
    N = length(pop.population_sizes)
    tnv = Vector{Float64}(undef, 2*N)
    
    tnv[1] = pop.genome_length
    for i in 1:N
        tnv[2*(N-i+1)] = pop.population_sizes[i]
    end
    for i in 2:N
        tnv[2*(N - i + 1)+1] = pop.times[i] - pop.times[i-1]
    end
    return tnv
end

Base.show(io::IO, d::VaryingPopulation) = print(io, "VaryingPopulation(ploidy=$(d.ploidy), population_sizes=$(d.population_sizes))")
