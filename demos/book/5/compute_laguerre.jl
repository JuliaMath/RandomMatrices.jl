## This function computes the Generalized Laguerre polynomials using the Tridiagonal matrix
include("orthopoly_evaluate.jl")
function compute_laguerre(n, a, x, all=true)
    # Evaluate Orthonormalized Laguerre polynomials at x
    #  Laguerre:
    #     w(x)=x^a exp(-x) with integral gamma(a+1)  (a>-1)
    #     Equivalent Mathematica:LaguerreL[n,a,x]/Sqrt[Gamma[n+a+1]/Gamma[n+1]] 
    ## Laguerre Bidiagonal matrix
    B = diagm(sqrt((a + 1 +(0:n)))) # diagonal
    B = B - diagm(sqrt((1:n)), 1)   # superdiagonal
    T = B'*B
    c = gamma(a+1) #Integral of weight function
    ## Default to 0 through n computation
    if all
        phi = orthopoly_evaluate_all(T, x)
    else
        phi = orthopoly_evaluate_n(T, x)
    end
    phi/sqrt(c)
end
