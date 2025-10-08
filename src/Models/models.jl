abstract type AbstractEvolutionaryModel end

export AbstractEvolutionaryModel, sim_ancestry,
    # SimulatedAncestry,  
    Hudson, SMC, SMCprime, WrightFisher

function sim_ancestry end


struct SimulatedAncestry{M <: AbstractEvolutionaryModel,S,D}
    model::M
    demography::Demography
    genome::Genome
    sample::S
    data::D
end




include("SMC/SMC.jl")
include("SMC/SMCprime.jl")
include("WrightFisherForward/WrightFisherForward.jl")
include("Hudson/Hudson.jl")

import .HudsonModel: Hudson  # needs to be imported for export to work
import .SMCapprox: SMC  # needs to be imported for export to work
import .SMCprimeapprox: SMCprime  # needs to be imported for export to work
import .WrightFisherForwardModel: WrightFisher





