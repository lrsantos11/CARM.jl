"""
This script builds the results and plots presented in Section 4.1 of [^Behling2020]

[^Behling2020] R. Behling, J.-Y. Bello-Cruz, e L.-R. Santos, “On the Circumcentered-Reflection Method 
for the Convex Feasibility Problem”, Numer. Algorithms, jul. 2020, doi: [10.1007/s11075-020-00941-6](https://doi.org/10.1007/s11075-020-00941-6). 
"""

include("../scripts/plots_util.jl")

function TestEpiQuadratic(;n::Int64 = 200,samples::Int64 = 10, restarts::Int64 = 1,
        EPSVAL::Float64 = 1e-6,itmax::Int64 = 5_000,λ::Float64=1.,σ::Float64=0.)
    # Fix Random
    Random.seed!(310595)
    # Defines DataFrame for Results
    dfResults= DataFrame(Problem=String[],CARMit=Int[],AMAPit=Int[],CRMit=Int[],MAPit=Int[])

    for j in 1:samples
        
        α = rand(Uniform(0,10))
        A = [zeros(n);1]
        b =  abs(rand(Normal(0,σ)))
    
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
σ=5.

dfResults = TestEpiQuadratic(n = n, samples = samples,itmax=itmax, EPSVAL = ε, restarts = restarts, λ=1.,σ=0.)
dfResults25 = TestEpiQuadratic(n = n, samples = samples,itmax=itmax, EPSVAL = ε, restarts = restarts,λ =1.25,σ=0.)
dfResults50 = TestEpiQuadratic(n = n, samples = samples,itmax=itmax, EPSVAL = ε, restarts = restarts,λ =1.50,σ=0.)
dfResults75 = TestEpiQuadratic(n = n, samples = samples,itmax=itmax, EPSVAL = ε, restarts = restarts,λ =1.75,σ=0.)

dfResultsEB = TestEpiQuadratic(n = n, samples = samples,itmax=itmax, EPSVAL = ε, restarts = restarts, λ=1.,σ=σ)
dfResults25EB = TestEpiQuadratic(n = n, samples = samples,itmax=itmax, EPSVAL = ε, restarts = restarts,λ =1.25,σ=σ)
dfResults50EB = TestEpiQuadratic(n = n, samples = samples,itmax=itmax, EPSVAL = ε, restarts = restarts,λ =1.50,σ=σ)
dfResults75EB = TestEpiQuadratic(n = n, samples = samples,itmax=itmax, EPSVAL = ε, restarts = restarts,λ =1.75,σ=σ)

@show describe(dfResults)
@show describe(dfResults25)
@show describe(dfResults50)
@show describe(dfResults75)

@show describe(dfResultsEB)
@show describe(dfResults25EB)
@show describe(dfResults50EB)
@show describe(dfResults75EB)

perprof1 = performance_profile(hcat(dfResults.CARMit, dfResults.AMAPit, dfResults.CRMit, dfResults.MAPit), ["CARM","AMAP","CRM","MAP"],
    title=L"Performance Profile (no EB) -- Gap error -- $\varepsilon = 10^{-6}$", 
#    linestyles=[:solid, :dash, :dot, :dashdot],
    legend = :bottomright, framestyle = :box,
    logscale=true)
ylabel!("Percentage of problems solved")
savefig(perprof1,("./plots/fig1.pdf"))
perprof1

perprof1EB = performance_profile(hcat(dfResultsEB.CARMit, dfResultsEB.AMAPit, dfResultsEB.CRMit, dfResultsEB.MAPit), ["CARM","AMAP","CRM","MAP"],
    title=L"Performance Profile (with EB) -- Gap error -- $\varepsilon = 10^{-6}$",
 #   linestyles=[:solid, :dash, :dot, :dashdot],
    legend = :bottomright, framestyle = :box,
    logscale=true)
ylabel!("Percentage of problems solved")
savefig(perprof1EB,("./plots/fig1EB.pdf"))
perprof1EB

perprof2 = performance_profile(hcat(dfResults.CARMit, dfResults.AMAPit, dfResults.MAPit), ["CARM","AMAP","MAP"],
    title=L"Performance Profile (no EB) -- Gap error -- $\varepsilon = 10^{-6}$",
  #  linestyles=[:solid, :dash, :dashdot],
    legend = :bottomright, framestyle = :box,
    logscale=true)
ylabel!("Percentage of problems solved")
savefig(perprof2,("./plots/fig2.pdf"))
perprof2

perprof2EB = performance_profile(hcat(dfResultsEB.CARMit, dfResultsEB.AMAPit, dfResultsEB.MAPit), ["CARM","AMAP","MAP"],
    title=L"Performance Profile (with EB) -- Gap error -- $\varepsilon = 10^{-6}$",
   # linestyles=[:solid, :dash, :dashdot],
    legend = :bottomright, framestyle = :box,
    logscale=true)
ylabel!("Percentage of problems solved")
savefig(perprof2EB,("./plots/fig2EB.pdf"))
perprof2EB

perprof3 = performance_profile(hcat(dfResults.CARMit, dfResults.AMAPit), ["CARM","AMAP"],
    title=L"Performance Profile (no EB) -- Gap error -- $\varepsilon = 10^{-6}$",
    #linestyles=[:solid, :dash],
    legend = :bottomright, framestyle = :box,
    logscale=true)
ylabel!("Percentage of problems solved")
savefig(perprof3,("./plots/fig3.pdf"))
perprof3

perprof3EB = performance_profile(hcat(dfResultsEB.CARMit, dfResultsEB.AMAPit), ["CARM","AMAP"],
    title=L"Performance Profile (with EB) -- Gap error -- $\varepsilon = 10^{-6}$",
    #linestyles=[:solid, :dash],
    legend = :bottomright, framestyle = :box,
    logscale=true)
ylabel!("Percentage of problems solved")
savefig(perprof3EB,("./plots/fig3EB.pdf"))
perprof3EB

perprof4 = performance_profile(hcat(dfResults.CARMit, dfResults25.CARMit, dfResults50.CARMit, dfResults.CRMit), ["CARM","CARM (λ=1.25)","CARM (λ=1.5)","CRM"],
    title=L"Performance Profile (no EB) -- Gap error -- $\varepsilon = 10^{-6}$",
    #linestyles=[:solid, :solid, :solid, :dot],
    legend = :bottomright, framestyle = :box,
    logscale=true)
ylabel!("Percentage of problems solved")
savefig(perprof4,("./plots/fig4.pdf"))
perprof4

perprof4EB = performance_profile(hcat(dfResultsEB.CARMit, dfResults25EB.CARMit, dfResults50EB.CARMit, dfResultsEB.CRMit), ["CARM","CARM (λ=1.25)","CARM (λ=1.5)","CRM"],
    title=L"Performance Profile (with EB) -- Gap error -- $\varepsilon = 10^{-6}$",
    #linestyles=[:solid, :solid, :solid, :dot],
    legend = :bottomright, framestyle = :box,
    logscale=true)
ylabel!("Percentage of problems solved")
savefig(perprof4EB,("./plots/fig4EB.pdf"))
perprof4EB

perprof5 = performance_profile(hcat(dfResults50.CARMit, dfResults75.CARMit, dfResults.CRMit), ["CARM (λ=1.5)","CARM (λ=1.75)","CRM"],
    title=L"Performance Profile (no EB) -- Gap error -- $\varepsilon = 10^{-6}$",
    #linestyles=[:solid, :solid, :dot],
    legend = :bottomright, framestyle = :box,
    logscale=true)
ylabel!("Percentage of problems solved")
savefig(perprof5,("./plots/fig5.pdf"))
perprof5

perprof5EB = performance_profile(hcat(dfResults50EB.CARMit, dfResults75EB.CARMit, dfResultsEB.CRMit), ["CARM (λ=1.5)","CARM (λ=1.75)","CRM"],
    title=L"Performance Profile (with EB) -- Gap error -- $\varepsilon = 10^{-6}$",
    #linestyles=[:solid, :solid, :dot],
    legend = :bottomright, framestyle = :box,
    logscale=true)
ylabel!("Percentage of problems solved")
savefig(perprof5EB,("./plots/fig5EB.pdf"))
perprof5EB



