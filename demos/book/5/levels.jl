#levels.jl
#Algorithm 5.7 of Random Eigenvalues by Alan Edelman

function levels(n)
    # Plot exact semicircle formula for GUE
    # This is given in formula (5.2.16) in Mehta as the sum of phi_j^2
    # but Christoffel-Darboux seems smarter
    ## Parameter
    # n: Dimension of matrix for which to determine exact semicircle
    x = (-1:0.001:1) * sqrt(2n) * 1.3
    pold = zeros(size(x))   # -1st Hermite Polynomial
    p = ones(size(x))       #  0th Hermite Polynomial
    k = p
    for j = 1:n # Compute the three-term recurrence
        pnew = (sqrt(2) * x .* p - sqrt(j - 1) * pold) / sqrt(j)
        pold = p
        p = pnew
    end
    pnew = (sqrt(2) * x .* p - sqrt(n) * pold)/sqrt(n + 1)

    k = n * p .^ 2 - sqrt(n * (n+1)) * pnew .* pold # Use p.420 of Mehta
    k = k .* exp(-x .^ 2) / sqrt(pi) # Normalize
    
    # Rescale so that "semicircle" is on [-1, 1] and area is pi/2
    (x / sqrt(2n), k * pi / sqrt(2n))
end

