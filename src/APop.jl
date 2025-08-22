module APop

using DataStructures
using Random
using Distributions

include("Genome/RateDistributions.jl")
include("Genome/genome.jl")



include("ARGs/segments.jl")
include("ARGs/trees.jl")
include("ARGs/args.jl")
include("ARGs/IBS.jl")

include("Demographies/Simple.jl")
include("Demographies/General.jl")



include("Models/WrightFisherForward/CrossOverStores.jl")
include("Models/WrightFisherForward/WrightFisherForward.jl")
include("Models/WrightFisherForward/arg.jl")


end
