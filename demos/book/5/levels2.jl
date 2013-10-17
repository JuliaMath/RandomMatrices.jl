#levels2.jl
#Algorithm 5.8 of Random Eigenvalues by Alan Edelman

function levels2(n)
    # Plot exact semicircle formula for GUE
    xfull = (-1:0.001:1) * sqrt(2n) * 1.3

    # Form the T matrix
    T = diagm(sqrt(1:n-1),1)
    T = T + T'
    
    # Do the eigen-decomposition of T, T = UVU'
    V, U = eig(T)

    # Precompute U'*e_n
    # tmp_en = U' * ((0:n-1) == n-1)'
    tmp_en = U[end, :]'

    y = Array(Float64, length(xfull))
    for i = 1:length(xfull)
        x = xfull[i]
        # generate the v vector as in (2.5)
        v = U * (tmp_en ./ (V - sqrt(2) * x))
        # multiply the normalization term
        y[i] = norm((sqrt(pi)) ^ (-1/2) * exp(-x^2/2) * v/v[1])^2
    end

    # Rescale so that "semicircle" is on [-1, 1] and area is pi/2
    (xfull/sqrt(2n), y * pi / sqrt(2n))
end
