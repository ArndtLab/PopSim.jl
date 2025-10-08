module PopSim

using DataStructures
using Random
using Distributions

include("Genomics/RateDistributions.jl")
include("Genomics/genome.jl")

include("ARGs/segments.jl")
include("ARGs/trees.jl")
include("ARGs/args.jl")

include("Demographies/Simple.jl")
include("Demographies/General.jl")
include("Demographies/Samples.jl")

include("Models/models.jl")


include("ARGs/IBS.jl")


end
