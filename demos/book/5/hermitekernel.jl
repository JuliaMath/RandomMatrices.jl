include("compute_hermite.jl")
function hermitekernel(n,x)
    #Compute the Hermite kernel matrix at points x
    #including degrees from 0 to n inclusive
    
    # # Method 1: Direct Computation
       phi=compute_hermite(n,x);
       K=phi'*phi;
       
    # # Method 2: Christoffel-Darboux
    #   phinp1 = compute_hermite(n+1,x,false)
    #   phin = compute_hermite(n,x,false)
    #   K = phinp1'*phin*sqrt((n+1)/2)
    #   x = x[:]*ones(1,length(x))
    #   K = (K-K')./(x-x')

    K
end
