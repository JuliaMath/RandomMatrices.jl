# Generates samples of dense random matrices that are distributed according
# to the classical Gaussian random matrix ensembles
#
# The notation follows closely that in:
#
# Ioana Dumitriu and Alan Edelman, "Matrix Models for Beta Ensembles"
# Journal of Mathematical Physics, vol. 43 no. 11 (2002), pp. 5830--5547
# doi: 10.1063/1.1507823
# arXiv: math-ph/0206043

#Generates a NxN symmetric Wigner matrix
#Hermite ensemble
function GaussianHermiteMatrix(n :: Integer, beta :: Integer)
    if beta == 1 #real
        A = randn(n, n)
        normalization = sqrt(2*n)
    elseif beta == 2 #complex
        A = randn(n, n) + im*randn(n, n)
        normalization = sqrt(4*n)
    elseif beta == 4 #quaternion
        #Employs 2x2 matrix representation of quaternions
        X = randn(n, n) + im*randn(n, n)
        Y = randn(n, n) + im*randn(n, n)
        A = [X Y; -conj(Y) conj(X)]
        normalization = sqrt(8*n) #TODO check normalization
    else
        error(@sprintf("beta = %d is not implemented", beta))
    end
    return (A + A') / normalization
end

#Generates a NxN Hermitian Wishart matrix
#Laguerre ensemble
#TODO Check - the eigenvalue distribution looks funky
function GaussianLaguerreMatrix(m :: Integer, n :: Integer, beta :: Integer)
    if beta == 1 #real
        A = randn(m, n)
    elseif beta == 2 #complex
        A = randn(m, n) + im*randn(m, n)
    elseif beta == 4 #quaternion
        #Employs 2x2 matrix representation of quaternions
        X = randn(m, n) + im*randn(m, n)
        Y = randn(m, n) + im*randn(m, n)
        A = [X Y; -conj(Y) conj(X)]
        error(@sprintf("beta = %d is not implemented", beta))
    end
    return (A * A') / m
end

#Generates a NxN self-dual MANOVA Matrix
#Jacobi ensemble
function GaussianJacobiMatrix(m :: Integer, n1 :: Integer, n2 :: Integer, beta :: Integer)
    w1 = Wishart(m, n1, beta)
    w2 = Wishart(m, n2, beta)
    return (w1 + w2) \ w1
end


#A convenience function
using Distributions
chi(df) = sqrt(rand(Chisq(df)))

#Generates a NxN tridiagonal Wigner matrix
#Hermite ensemble
function GaussianHermiteTridiagonalMatrix(n :: Integer, beta :: FloatingPoint)
    @assert beta > 0
    Hdiag = randn(n)/sqrt(n)
    Hsup = [chi(beta*i)/sqrt(2*n) for i=n-1:-1:1]
    return SymTridiagonal(Hdiag, Hsup)
end


#TODO Check normalization - I guessed this
function GaussianLaguerreTridiagonalMatrix(m :: Integer, a :: FloatingPoint, beta :: FloatingPoint)
    if a <= beta*(m-1)/2.0
        error(@sprintf("Given your choice of m and beta, a must be at least %f (You said a = %f)", beta*(m-1)/2.0, a))
    end
    Hdiag = [chi(2*a-i*beta) for i=0:m-1]
    Hsub = [chi(beta*i) for i=m-1:-1:1]
    #Julia has no bidiagonal type... yet
    B = Tridiagonal(Hsub, Hdiag, zeros(m-1))
    L = B * B' #Currently returns a dense matrix
    return SymTridiagonal(diag(L)/sqrt(m), diag(L,1)/sqrt(m))
end

