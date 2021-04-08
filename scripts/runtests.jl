using DrWatson
@quickactivate "CARM"

using Random, Distributions
include("../src/CARM.jl")
using Random

include("ExactProjectionQuadratic.jl")
include("InexactProjectionQuadratic.jl")
include("InexactProjectionEllipsis.jl")