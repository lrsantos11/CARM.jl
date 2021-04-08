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

# Outer Projection onto convex $K = \\operatorname(epi)(\alpha x^Tx)$

g(x) = α*dot(x[1:end-1],x[1:end-1]) - x[end]
∂g(x) = [2*α*x[1:end-1]; -1]
λ=1.75
ApProjectK(x) = RelaxApproxProject(x,g,∂g,λ=λ)
ApReflectK(x) = Reflection(x,ApProjectK)

filedirCARM = joinpath(datadir(),"sims","CARM_Epi.dat")
xCARM, tol, k = CRM(x₀,ApReflectK,ReflectU,filedir=filedirCARM)
@show "CARM", tol, k, λ
filedirAMAP = joinpath(datadir(),"sims","AMAP_Epi.dat")
xAMAP, tol, k = MAP(x₀,ApProjectK,ProjectU,itmax=2000,filedir=filedirAMAP)
@show "AMAP", tol, k, λ;
