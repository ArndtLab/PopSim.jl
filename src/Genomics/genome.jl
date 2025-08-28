
export Genome, recombination, mutation




struct Genome{R <: AbstractRateDistribution, M <: AbstractRateDistribution}
    recombination::R
    mutation::M
    length::Int
end

function Genome(; recombination_rate::Float64 = 0.0, mutation_rate::Float64 = 0.0, length::Int = 0)
    recombination = UniformRate(recombination_rate)
    mutation = UniformRate(mutation_rate)
    return Genome(recombination, mutation, length)
end



Base.length(g::Genome) = g.length
recombination(g::Genome) = g.recombination
mutation(g::Genome) = g.mutation   

Base.show(io::IO, g::Genome) = print(io, "Genome(recombination=$(g.recombination), mutation=$(g.mutation), length=$(g.length))")




