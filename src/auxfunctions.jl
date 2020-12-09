__precompile__()
using LinearAlgebra
using Printf
using Random
using ProximalOperators
using PolynomialRoots

####################################

"""jldoctest
FindCircumcentermSet(X)

Finds the Circumcenter of linearly independent vectors ``x_0,x_1,…,x_m``, columns of matrix ``X``,
as described in [^Behling2018a] and [^Behling2018b].

[^Behling2018a]: Behling, R., Bello Cruz, J.Y., Santos, L.-R.:
Circumcentering the Douglas–Rachford method. Numer. Algorithms. 78(3), 759–776 (2018).
[doi:10.1007/s11075-017-0399-5](https://doi.org/10.1007/s11075-017-0399-5)
[^Behling2018b]: Behling, R., Bello Cruz, J.Y., Santos, L.-R.:
On the linear convergence of the circumcentered-reflection method. Oper. Res. Lett. 46(2), 159-162 (2018).
[doi:10.1016/j.orl.2017.11.018](https://doi.org/10.1016/j.orl.2017.11.018)
"""
function FindCircumcentermSet(X)
    # Finds the Circumcenter of  linearly independent points  X = [X1, X2, X3, ... Xn]
        # println(typeof(X))
        lengthX = length(X)
        if lengthX  == 1
            return X[1]
        elseif lengthX == 2
            return .5*(X[1] + X[2])
        end
        V = []
        b = Float64[]
        # Forms V = [X[2] - X[1] ... X[n]-X[1]]
        # and b = [dot(V[1],V[1]) ... dot(V[n-1],V[n-1])]
        for ind in 2:lengthX
            difXnX1 = X[ind]-X[1]
            push!(V,difXnX1)
            push!(b,dot(difXnX1,difXnX1))
        end

       # Forms Gram Matrix
        dimG = lengthX-1
        G = diagm(b)

        for irow in 1:(dimG-1)
            for icol in  (irow+1):dimG
                G[irow,icol] = dot(V[irow],V[icol])
                G[icol,irow] = G[irow,icol]
            end
        end
        # println(rank(G))
        y = G\b
        # if isposdef(G)
        #     L = cholesky(G)
        #     y = L\b
        # else
        #     @warn "Gram matrix is not SPD"
        #     L = qr(G)
        #     y=L\b
        # end
        CC = X[1]
        for ind in 1:dimG
            CC += .5*y[ind]*V[ind]
        end
        return CC
    end
####################################
"""
    proj = ProjectIndicator(indicator,x)
    Projection using Indicator Function from `ProximalOperators.jl`
    """
    function ProjectIndicator(indicator,x)
        proj, fproj = prox(indicator,x)
        return proj
    end

####################################
    """
    reflec = ReflectIndicator(indicator,x)
    Reflection using Indicator Function from `ProximalOperators.jl`
    """
    function ReflectIndicator(indicator,x)
        proj = ProjectIndicator(indicator,x)
        reflec = 2*proj - x
    end
####################################
function StartingPoint(n::Int64)
    ## Creates a random point in R^n
    x=zeros(n);
    while norm(x)<2
        x = randn(n);
    end
    # norm between 5 and 15
    foonorm = (15-5)*rand() + 5
    return foonorm*x/norm(x);
end

####################################
function StartingPoint(m::Int64,n::Int64)
    ## Creates a random point in R^{m\times n}
    X=zeros(m,n)
    for j in 1:n
        X[:,j] = StartingPoint(m)
    end
    return X
end

"""
        ProjectEpigraphofQuadratic(v,s)
Project ``(v,t) ∈ R^{n+1}`` onto the epigraph of ``f(x) = αx^Tx``, such that ``f(v) '≤ t``.
"""
function ProjectEpigraphofQuadratic(v::Union{AbstractArray,Number},t::Number, α::Float64=1.0)
        #PolynomialCoefficient
        if α*dot(v,v) <= t
                return v, t
        end
        a0 = t-α*dot(v,v)
        a1 = 4*α*t + 1.
        a2 = 4*α^2*t + 4*α
        a3 = 4*α^2
        r = (roots([a0,a1,a2,a3]))
        indexreal =  findall(x->abs.(x)<1e-12,imag.(r))
        λ = (maximum(real.(r[indexreal])))
        x = 1/(1+2*α*λ) * v
        t += λ
        return x, t
end
