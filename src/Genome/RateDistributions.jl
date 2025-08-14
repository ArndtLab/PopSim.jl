
export AbstractRateDistribution, UniformRate


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




