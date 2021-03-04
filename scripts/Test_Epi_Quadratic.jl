"""
This script builds the results and plots presented in Section 4.1 of [^Behling2020]

[^Behling2020] R. Behling, J.-Y. Bello-Cruz, e L.-R. Santos, “On the Circumcentered-Reflection Method 
for the Convex Feasibility Problem”, Numer. Algorithms, jul. 2020, doi: [10.1007/s11075-020-00941-6](https://doi.org/10.1007/s11075-020-00941-6). 
"""

include("../scripts/plots_util.jl")

function TestEpiQuadratic(;n::Int64 = 200,samples::Int64 = 10, restarts::Int64 = 1,
        EPSVAL::Float64 = 1e-6,itmax::Int64 = 5_000,λ::Float64=1.)
    # Fix Random
    Random.seed!(1234)
    # Defines DataFrame for Results
    dfResults= DataFrame(Problem=String[],CARMit=Int[],AMAPit=Int[],CRMit=Int[],MAPit=Int[])

    for j in 1:samples
        
        α = rand(Uniform(0,10))
        A = [zeros(n);1]
        b =  abs(rand(Normal(0,10)))
        
        Affine = IndAffine(A,b)
        ProjectU(x) = ProjectIndicator(Affine,x)
        ReflectU(x) = Reflection(x,ProjectU)
        
        ProjectK(x) = ProjectEpiQuadratic(x,α=α)
        ReflectK(x) = Reflection(x,ProjectK)
        
        ApProjectK(x) = RelaxApproxProject(x,g,∂g,λ=λ)
        ApReflectK(x) = Reflection(x,ApProjectK)        
        
        g(x) = α*dot(x[1:end-1],x[1:end-1]) - x[end]
        ∂g(x) = [2*α*x[1:end-1]; -1]
        # Restarts
        for i = 1:restarts
            u₀ = StartingPoint(n)
            s₀ = StartingPoint(1)
            x₀ = ProjectU([u₀;s₀])
            # @show SOC(xzero)
            while g(x₀) < EPSVAL
                u₀ = StartingPoint(n)
                s₀ = StartingPoint(1)
                x₀ = ProjectU([u₀;s₀])
            end
            prob_name  = String("Problem$j"*"Restart$i")
            # println(prob_name)
            resultCARM  = CRM(x₀,ApReflectK,ReflectU,itmax=itmax,EPSVAL=EPSVAL)
            resultAMAP  = MAP(x₀,ApProjectK,ProjectU,itmax=itmax,EPSVAL=EPSVAL)
            resultCRM  = CRM(x₀,ReflectK,ReflectU,itmax=itmax,EPSVAL=EPSVAL)
            resultMAP  = MAP(x₀,ProjectK,ProjectU,itmax=itmax,EPSVAL=EPSVAL)
            push!(dfResults,(prob_name, resultCARM[3], resultAMAP[3], resultCRM[3], resultMAP[3]))
        end
    end
    return dfResults
end

n = 200
samples = 100
restarts = 10
ε = 1e-5
itmax = 5000
dfResults = TestEpiQuadratic(n = n, samples = samples,itmax=itmax, EPSVAL = ε, restarts = restarts, λ=1.)
dfResults25 = TestEpiQuadratic(n = n, samples = samples,itmax=itmax, EPSVAL = ε, restarts = restarts,λ =1.25)
dfResults50 = TestEpiQuadratic(n = n, samples = samples,itmax=itmax, EPSVAL = ε, restarts = restarts,λ =1.50)
dfResults75 = TestEpiQuadratic(n = n, samples = samples,itmax=itmax, EPSVAL = ε, restarts = restarts,λ =1.75)

perprof = performance_profile(hcat(dfResults.CARMit, dfResults.AMAPit), ["CARM","AMAP"],
    title=L"Performance Profile -- Gap error -- $\varepsilon = 10^{-6}$",
    legend = :bottomright, framestyle = :box, linestyles=[:solid, :dash],
    logscale=true)
ylabel!("Percentage of problems solved")
savefig(perprof,("./plots/fig.pdf"))
@show describe(dfResults)
@show describe(dfResults25)
@show describe(dfResults50)
@show describe(dfResults75)
perprof


