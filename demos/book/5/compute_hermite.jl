include("orthopoly_evaluate.jl")
function compute_hermite(n, x, all=true);
    # Evaluate Orthonormalized Hermite polynomials at x
    # Physicists Hermite:
    #       w(x)=exp(-x^2) with integral sqrt(pi)
    #       Equivalent Mathematica: HermiteH[n,x]/Sqrt[Sqrt[Pi]*n!*2^n] 
    # Probabilists Hermite:
    #       w(x)=exp(-x^2/2) with integral sqrt(2*pi)
    #       Equivalent Mathematica: HermiteH[n,x/Sqrt[2]]/Sqrt[Sqrt[Pi]*n!*2^(1/2+n)]
    
    ## Physicists Hermite Tridiagonal matrix (w=exp(-x^2))
    T = diagm(sqrt((1:n)/2),1) # n+1 by n+1 matrix 
    T = T+T'
    c = sqrt(pi) #Integral of weight function
    ## Probabilists Hermite Tridiagonal matrix (w=exp(-x^2/2))
    # T = diagm(sqrt((1:n)),1) # n+1 by n+1 matrix
    # T = T+T'
    # c = sqrt(2*pi)
    ## Default to 0 through n computation
    if all
        phi = orthopoly_evaluate_all(T, x)
    else
        phi = orthopoly_evaluate_n(T, x)
    end
    phi/sqrt(c)
end
