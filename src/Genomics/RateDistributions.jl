
export AbstractRateDistribution, UniformRate, NonUniformRate,
    average_rate, sample

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


function StatsBase.sample(u::UniformRate, dt::Float64, n1::Int, n2::Int;
    multiple_hits=:ignore)::Vector{Int}
    n = n2 - n1 + 1
    n == 0 && return Int64[]
    if dt == Inf && multiple_hits != :JCcorrect
        return collect(n1:n2)
    end
    if multiple_hits == :ignore
        k = rand(Poisson(n * u.rate * dt))
        k == 0 && return Int64[]
        return sort!(rand(n1:n2, k))
    elseif multiple_hits == :as_one
        k = rand(Poisson(n * u.rate * dt))
        k == 0 && return Int64[]
        return unique!(sort!(rand(n1:n2, k)))
    elseif multiple_hits == :JCcorrect
        prob_corrected = (3 / 4) * (1 - exp(-4 * u.rate * dt / 3))
        k = rand(Poisson(n * prob_corrected))
        k == 0 && return Int64[]
        k >= n && return collect(n1:n2)
        return sort!(sample(n1:n2, k; replace=false))
    end
end





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


