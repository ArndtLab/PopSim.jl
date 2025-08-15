
export AbstractRateDistribution, UniformRate, NonUniformRate,
    average_rate

using StatsBase, Distributions, Random

abstract type AbstractRateDistribution end


"""
    UniformRate(rate::Float64)

    A uniform rate distribution for genetic processes.
"""
struct UniformRate <: AbstractRateDistribution
    rate::Float64
end

rate(r::UniformRate) = r.rate

Base.show(io::IO, r::UniformRate) = print(io, "Uniform(rate=$(r.rate))")

average_rate(r::UniformRate) = r.rate


function StatsBase.sample(u::UniformRate, dt::Float64, n1::Int, n2::Int)
    n = n2 - n1 + 1
    k = rand(Poisson(n * u.rate * dt))
    sort!(rand(n1:n2, k))
end

StatsBase.sample(u::UniformRate, dt::Float64, n::Int) = sample(u, dt, 1, n)



struct NonUniformRate <: AbstractRateDistribution
    cprob::Vector{Float64}
    
    function NonUniformRate(ratesperposition::Vector{Float64})
        cprob = cumsum(ratesperposition)
        new(cprob)
    end
end
    
Base.show(io::IO, nu::NonUniformRate) = print(io, "NonUniformRateMap of length $(length(nu.cprob))")

average_rate(nu::NonUniformRate) = nu.cprob[end] / length(nu.cprob)

Base.length(nu::NonUniformRate) = length(nu.cprob)



function StatsBase.sample(nu::NonUniformRate, dt::Float64, n1::Int, n2::Int)
    L = n2 - n1 + 1
    R2 = nu.cprob[n2] 
    R1 = (n1 > 1 ? nu.cprob[n1 - 1] : 0.0)

    k = rand(Poisson((R2 - R1) * dt))


    xs = sort!((R2 - R1) .* rand(k) .+ R1)
    return [searchsortedlast(nu.cprob, x)+1 for x in xs]
end


