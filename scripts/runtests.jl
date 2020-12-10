using DrWatson
@quickactivate "CARM"

include("../src/CARM.jl")
using Random
Random.seed!(1234)

n = 1
A = [zeros(n);1]
b =  0.
Affine = IndAffine(A,b)
ProjectU(x) = ProjectIndicator(Affine,x)
ReflectU(x) = Reflection(x,ProjectU)
u₀ = StartingPoint(n)
s₀ = StartingPoint(1)
x₀ = ProjectU([u₀;s₀])

# Exact Projection onto convex $K = \\operatorname(epi)(\alpha x^Tx)$
α = 1.00
ProjectK(x) = ProjectEpiQuadratic(x,α=α)
ReflectK(x) = Reflection(x,ProjectK)


xCRM, tol, k = CRM(x₀,ReflectK,ReflectU)
@show tol, k
xMAP, tol, k = MAP(x₀,ProjectK,ProjectU,itmax=2000)
@show tol, k

# Outer Projection onto convex $K = \\operatorname(epi)(\alpha x^Tx)$

g(x) = α*dot(x[1:end-1],x[1:end-1]) - x[end]
∂g(x) = [2*α*x[1:end-1]; -1]
ApProjectK(x) = ApproxProject(x,g,∂g)
ApReflectK(x) = Reflection(x,ApProjectK)

filedirCARM = joinpath(datadir(),"sims","CARM_Epi.dat")
xCARM, tol, k = CRM(x₀,ApReflectK,ReflectU,filedir=filedirCARM)
@show tol, k
filedirAMAP = joinpath(datadir(),"sims","AMAP_Epi.dat")
xAMAP, tol, k = MAP(x₀,ApProjectK,ProjectU,itmax=2000,filedir=filedirAMAP)
@show tol, k;
