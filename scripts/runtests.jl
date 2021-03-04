using DrWatson
@quickactivate "CARM"

using Random, Distributions
include("../src/CARM.jl")
#Random.seed!(1234)

#n = 1
#A = [zeros(n);1]
#b =  0.
#Affine = IndAffine(A,b)
#ProjectU(x) = ProjectIndicator(Affine,x)
#ReflectU(x) = Reflection(x,ProjectU)
#u₀ = StartingPoint(n)
#s₀ = StartingPoint(1)
#x₀ = ProjectU([u₀;s₀])

#include("Exact_Epi_Quadratic.jl")
#include("Approx_Epi_Quadratic.jl")
include("makeplots.jl")