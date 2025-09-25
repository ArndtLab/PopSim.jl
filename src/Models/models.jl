abstract type AbstractEvolutionaryModel end

export AbstractEvolutionaryModel, sim_ancestry

function sim_ancestry end


include("SMC/SMC.jl")
include("SMC/SMCprime.jl")


include("WrightFisherForward/WrightFisherForward.jl")

include("Hudson/Hudson.jl")
