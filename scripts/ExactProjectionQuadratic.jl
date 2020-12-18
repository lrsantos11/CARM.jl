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
@show "CRM", tol, k
xMAP, tol, k = MAP(x₀,ProjectK,ProjectU,itmax=2000)
@show "MAP", tol, k