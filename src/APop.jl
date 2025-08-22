module APop

using DataStructures
using Random
using Distributions

include("Genome/RateDistributions.jl")
include("Genome/genome.jl")

include("Models/CrossOverStores.jl")


include("ARGs/segments.jl")
include("ARGs/trees.jl")
include("ARGs/args.jl")
include("ARGs/IBS.jl")

include("Demographies/Simple.jl")
include("Demographies/General.jl")



include("Models/WrightFisherForward.jl")


end
