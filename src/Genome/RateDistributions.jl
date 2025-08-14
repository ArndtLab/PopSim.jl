
export AbstractRateDistribution, UniformRate

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



function StatsBase.sample(u::UniformRate, n0::Int, n1::Int)
    n = n1 - n0 + 1
    k = rand(Poisson(n * u.rate))
    sort!(rand(n0:n1, k))
end

StatsBase.sample(u::UniformRate, n::Int) = sample(u, 1, n)
