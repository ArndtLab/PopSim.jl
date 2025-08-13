


export AbstractDemography, StationaryPopulation


abstract type AbstractDemography end

"""
    StationaryPopulation(ploidy::Int, size::Int)

    A simple demography for population modeling.
    A stationary population with a fixed ploidy and size.
"""
struct StationaryPopulation <: AbstractDemography
    ploidy::Int
    size::Int

    function StationaryPopulation(ploidy::Int, size::Int)
        if ploidy < 1 || ploidy > 2
            throw(ArgumentError("Ploidy must be 1 or 2."))
        end
        if size < 1
            throw(ArgumentError("Size must be a positive integer."))
        end
        return new(ploidy, size)
    end
end



StationaryPopulation(size::Int) = StationaryPopulation(2, size)

ploidy(d::StationaryPopulation) = d.ploidy
Base.size(d::StationaryPopulation) = d.size

Base.show(io::IO, d::StationaryPopulation) = print(io, "StationaryPopulation(ploidy=$(d.ploidy), size=$(d.size))")

