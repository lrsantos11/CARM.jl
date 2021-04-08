using DrWatson
@quickactivate "CARM"

using Random, Distributions
include(srcdir("CARM.jl"))

include("ExactProjectionQuadratic.jl")
include("InexactProjectionQuadratic.jl")
include("InexactProjectionEllipsis.jl")