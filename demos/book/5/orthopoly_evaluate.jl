function orthopoly_evaluate_all(T, x)
    ##  Evaluate 0th through nth orthogonal polynomial at x
    #   Parameters
    #     T: Tridiagonal matrix encoding the three term recurrence
    #        size is (n+1)x(n+1)
    #     x: evaluation points
    #
    #  Returns
    #    phi: Polynomials 0 through n at x
    #         Method:  Spectrally solve (T-XI)u = constant * (last column of I)
    #         WARNING: potential for dividing by 0
    #         1st row of phi is 1 -- for w not a probability measure
    #                divide by 1/sqrt(c), where c is the total weight
    #                of the real line
    n = size(T,1)
    Lambda, H = eig(T)
    Hn = H[end, :]'
    phi = Array(Float64, length(Lambda), length(x))
    for m = 1:length(x)
        xi = x[m]
        # generate the u vector as in (2.5)
        u = H * (Hn./(Lambda - xi))
        phi[:,m] = u/u[1]
    end
    phi
end

function orthopoly_evaluate_n(T, x)
    ##  Evaluate nth orthogonal polynomial at x
    #   Parameters
    #     T: Tridiagonal matrix encoding the three term recurrence
    #        size is (n+1)x(n+1)
    #     x: evaluation points
    #
    #  Returns
    #    phi: nth at x
    #         normalized so that 0'th is 1 (w has unit weight)
    #         Method: characteristic polynomial of n x n T
    #                 with leading coefficient taken from the 
    #                 product of the superdiagonal of the
    #                 (n+1)x(n+1) T
    
    n = size(T,1)-1
    ## compute the characteristic polynomial of T[1:n,1:n]
    Lambda = eigvals(T[1:n, 1:n])
    phi = ones(size(x))
    for i = 1:n
        phi = phi.*(x - Lambda[i])
    end
    ## normalization
    if n > 0
        phi = phi / prod(diag(T,1))
    end
    phi
end
