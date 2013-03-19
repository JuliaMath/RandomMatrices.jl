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


#Generates a NxN tridiagonal Wishart matrix
#Laguerre ensemble
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


# A helper function for Jacobi samples
function SampleCSValues(n :: Integer, a :: FloatingPoint, b :: FloatingPoint, beta :: FloatingPoint)
    if beta == Inf
        c=sqrt((a+[1:n])./(a+b+2*[1:n]))
        s=sqrt((b+[1:n])./(a+b+2*[1:n]))
        cp=sqrt([1:n-1]./(a+b+1+2*[1:n-1]))
        sp=sqrt((a+b+1+[1:n-1])./(a+b+1+2*[1:n-1]))
    else
        #Generate cosine-squared values
        csq = [rand(Beta(beta*(a+i)/2,beta*(b+i)/2)) for i=1:n]
        cpsq = [rand(Beta(beta*i/2,beta*(a+b+1+i)/2)) for i=1:n]
        #Cosine-sine pairs
        c, s = sqrt(csq), sqrt(1-csq)
        cp, sp = sqrt(cpsq), sqrt(1-cpsq)
    end
    return c, s, cp, sp
end
#Generates a 2Mx2M sparse MANOVA matrix
#Jacobi ensemble
#
# Reference:
#     A. Edelman and B. D. Sutton, "The beta-Jacobi matrix model, the CS decomposition,
#     and generalized singular value problems", Foundations of Computational Mathematics,
#     vol. 8 iss. 2 (2008), pp 259-285.
#TODO check normalization
function GaussianJacobiSparseMatrix(n :: Integer, a :: FloatingPoint, b :: FloatingPoint, beta :: FloatingPoint)
    CoordI = zeros(8n-4)
    CoordJ = zeros(8n-4)
    Values = zeros(8n-4)

    c, s, cp, sp = SampleCSValues(n, a, b, beta)

    #Diagonals of each block
    for i=1:n
        CoordI[i], CoordJ[i] = i, i
        Values[i] = i==1 ? c[n] : c[n+1-i] * sp[n+1-i]
    end    
    for i=1:n
        CoordI[n+i], CoordJ[n+i] = i, n+i
        Values[n+i] = i==n ? s[1] : s[n+1-i] * sp[n-i]
    end
    for i=1:n
        CoordI[2n+i], CoordJ[2n+i] = n+i, i
        Values[2n+i] = i==1 ? -s[n] : -s[n+1-i] * sp[n+1-i]
    end
    for i=1:n
        CoordI[3n+i], CoordJ[3n+i] = n+i, n+i
        Values[3n+i] = i==n ? c[1] : c[n+1-i] * sp[n-i]
    end
    #Off-diagonals of each block
    for i=1:n+1
        CoordI[4n+i], CoordJ[4n+i] = i,i+1
        Values[4n+i] = -s[n+1-i]*cp[n-i]
    end
    for i=1:n+1
        CoordI[5n-1+i], CoordJ[5n-1+i] = i+1,n+i
        Values[5n-1+i] = c[n-i]*cp[n-i]
    end
    for i=1:n+1
        CoordI[6n-2+i], CoordJ[6n-2+i] = n+i,i+1
        Values[6n-2+i] = -c[n+1-i]*cp[n-i]
    end
    for i=1:n+1
        CoordI[7n-3+i], CoordJ[7n-3+i] = n+i,i+1
        Values[7n-3+i] = -s[n-i]*cp[n-i]
    end
    
    return sparse(CoordI, CoordJ, Values)
end

#Return n eigenvalues distributed according to the Hermite ensemble
function GaussianHermiteSamples(n :: Integer, beta :: FloatingPoint)
    eigvals(GaussianHermiteTridiagonalMatrix(n, beta))
end

#Return n eigenvalues distributed according to the Laguerre ensemble
function GaussianLaguerreSamples(m :: Integer, a :: FloatingPoint, beta :: FloatingPoint)
    eigvals(GaussianLaguerreTridiagonalMatrix(m, a, beta))
end

#Return n eigenvalues distributed according to the Jacobi ensemble
function GaussianJacobiSamples(n :: Integer, a :: FloatingPoint, b :: FloatingPoint, beta :: FloatingPoint)
    #Generate just the upper left quadrant of the matrix
    c, s, cp, sp = SampleCSValues(n, a, b, beta)
    dv = [i==1 ? c[n] : c[n+1-i] * sp[n+1-i] for i=1:n]
    ev = [-s[n+1-i]*cp[n-i] for i=1:n-1]
    M = Tridiagonal(zeros(n-1), dv, ev)
    return svdvals(full(M)) #No SVD for tridiagonal matrices... yet
end
    
