using DrWatson
@quickactivate "CARM"

include("../src/CARM.jl")
using Random

include("ExactProjectionQuadratic.jl")
include("InexactProjectionQuadratic.jl")
include("InexactProjectionEllipsis.jl")