abstract type AbstractEvolutionaryModel end

export AbstractEvolutionaryModel, sim_ancestry, get_ARGsegments, get_IBDsegments

function sim_ancestry end
# function get_ARGsegments end
# function get_IBDsegments end


include("SMC/SMC.jl")
include("SMC/SMCprime.jl")


include("WrightFisherForward/WrightFisherForward.jl")

include("Hudson/Hudson.jl")
