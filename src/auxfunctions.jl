__precompile__()
using LinearAlgebra
using Printf
using Random
using ProximalOperators

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
    y = G\b
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
