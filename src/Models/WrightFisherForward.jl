

abstract type AbstractEvolutionaryModel end


export WrightFisher, AbstractEvolutionaryModel


struct WrightFisher <: AbstractEvolutionaryModel
    allow_selfing::Bool
end

WrightFisher(;allow_selfing::Bool = true) = WrightFisher(allow_selfing)






# ============================================================
# Individual
# ============================================================

using StaticArrays
using DataStructures


export Individual, randallele, ploidy

struct Individual{P} 
    alleles::SVector{P, Int64}
end

Individual(alleles...) = Individual(SVector{length(alleles)}(alleles...))
Individual(alleles::Vector) = Individual(SVector{length(alleles)}(alleles))


Base.show(io::IO, ind::Individual{P}) where P = print(io, "Individual (ploidy: $P, $(ind.alleles))")

DataStructures.@delegate Individual.alleles (Base.getindex, Base.setindex!, Base.length)

ploidy(::Individual{P}) where P = P

randallele(ind::Individual) = rand(ind.alleles)
