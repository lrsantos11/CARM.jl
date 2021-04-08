Random.seed!(1234)

n = 200
A = [zeros(n);1]
b =  0.
Affine = IndAffine(A,b)
ProjectU(x) = ProjectIndicator(Affine,x)
ReflectU(x) = Reflection(x,ProjectU)
u₀ = StartingPoint(n)
s₀ = StartingPoint(1)
x₀ = ProjectU([u₀;s₀])

# Outer Projection onto convex $K = \\operatorname(epi)(\alpha x^Tx)$
α = 1.00
v = rand(n)
D = diagm(v)
g(x) = α*dot(x[1:end-1],D*x[1:end-1]) - x[end]
∂g(x) = [2*α*D*x[1:end-1]; -1]
λ=1.
ApProjectK(x) = ApproxProject(x,g,∂g,λ=λ)
ApReflectK(x) = Reflection(x,ApProjectK)

filedirCARM = datadir("sims","CARM_Epi.dat")
xCARM, tol, k = CRM(x₀,ApReflectK,ReflectU,filedir=filedirCARM,itmax=2000)
@show "CARM", tol, k, λ
filedirAMAP = datadir("sims","AMAP_Epi.dat")
xAMAP, tol, k = MAP(x₀,ApProjectK,ProjectU,itmax=7000,filedir=filedirAMAP)
@show "AMAP", tol, k, λ;
