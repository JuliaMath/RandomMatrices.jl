# Generates samples of dense random matrices that are distributed according
# to the classical Gaussian random matrix ensembles
#
# The notation follows closely that in:
#
# Ioana Dumitriu and Alan Edelman, "Matrix Models for Beta Ensembles"
# Journal of Mathematical Physics, vol. 43 no. 11 (2002), pp. 5830--5547
# doi: 10.1063/1.1507823
# arXiv: math-ph/0206043
#
# References
#
# Alan Edelman, Per-Olof Persson and Brian D Sutton, "The fourfold way"
# http://www-math.mit.edu/~edelman/homepage/papers/ffw.pdf

using Distributions
export GaussianHermiteMatrix, GaussianLaguerreMatrix, GaussianJacobiMatrix,
       GaussianHermiteTridiagonalMatrix, GaussianLaguerreTridiagonalMatrix,
       GaussianJacobiSparseMatrix,
       GaussianHermiteSamples, GaussianLaguerreSamples, GaussianJacobiSamples


#########################
# Convenience functions #
#########################

#Produces a random variate of the chi distribution
chi(df) = df==0? 0.0 : sqrt(rand(Chisq(df)))

#############################
# Gaussian Wigner  ensemble #
# Gaussian Hermite ensemble #
#############################

#Generates a NxN symmetric Wigner matrix
function GaussianHermiteMatrix(n::Integer, beta::Integer)
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

#Generates a NxN tridiagonal Wigner matrix
#The beta=infinity case is defined in Edelman, Persson and Sutton, 2012
function GaussianHermiteTridiagonalMatrix(n::Integer, beta::Real)
    if beta<=0 error("beta must be positive") end
    if beta==Inf return SymTridiagonal(zeros(n), [sqrt(x/2) for x=n-1:-1:1]) end
    Hdiag = randn(n)/sqrt(n)
    Hsup = [chi(beta*i)/sqrt(2*n) for i=n-1:-1:1]
    return SymTridiagonal(Hdiag, Hsup)
end

##############################
# Gaussian Wishart  ensemble #
# Gaussian Laguerre ensemble #
##############################

#Generates a NxN Hermitian Wishart matrix
#TODO Check - the eigenvalue distribution looks funky
function GaussianLaguerreMatrix(m::Integer, n::Integer, beta::Integer)
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


#Generates a NxN bidiagonal Wishart matrix
#Laguerre ensemble
function GaussianLaguerreBidiagonalMatrix(m::Integer, a::Real, beta::Real)
    min_a = beta*(m-1)/2
    a<min_a ? error(@sprintf("Given your choice of m and beta, a must be at least %f (You said a = %f)", min_a, a)) : nothing
    Hdiag = [chi(2*a-i*beta) for i=0:m-1]
    Hsub = [chi(beta*i) for i=m-1:-1:1]
    Bidiagonal(Hsub, Hdiag, false)/(m^0.25)
end


#Generates a NxN tridiagonal Wishart matrix
#Laguerre ensemble
function GaussianLaguerreTridiagonalMatrix(m::Integer, a::Real, beta::Real)
    B = GaussianLaguerreBidiagonalMatrix(m, a, beta)
    B * B'
end


#Return n eigenvalues distributed according to the Hermite ensemble
function GaussianHermiteSamples(n::Integer, beta::Real)
    eigvals(GaussianHermiteTridiagonalMatrix(n, beta))
end


############################
# Gaussian MANOVA ensemble #
# Jacobi ensemble          #
############################

#Generates a NxN self-dual MANOVA Matrix
function GaussianJacobiMatrix(m::Integer, n1::Integer, n2::Integer, beta::Integer)
    w1 = Wishart(m, n1, beta)
    w2 = Wishart(m, n2, beta)
    return (w1 + w2) \ w1
end


# A helper function for Jacobi samples
function SampleCSValues(n::Integer, a::Real, b::Real, beta::Real)
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


#Return n eigenvalues distributed according to the Laguerre ensemble
#Compute the singular values of the bidiagonal matrix
function GaussianLaguerreSamples(m::Integer, a::Real, beta::Real)
    svdvals(GaussianLaguerreBidiagonalMatrix(m, a, beta))
end


#Generates a 2Mx2M sparse MANOVA matrix
#Jacobi ensemble
#
# Reference:
#     A. Edelman and B. D. Sutton, "The beta-Jacobi matrix model, the CS decomposition,
#     and generalized singular value problems", Foundations of Computational Mathematics,
#     vol. 8 iss. 2 (2008), pp 259-285.
#TODO check normalization
function GaussianJacobiSparseMatrix(n::Integer, a::Real, b::Real, beta::Real)
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


#Return n-dimensional bidiagonal matrix representation of the Jacobi ensemble
function GaussianJacobiBidiagonalMatrix(n::Integer, a::Real, b::Real, beta::Real)
    #Generate just the upper left quadrant of the matrix
    c, s, cp, sp = SampleCSValues(n, a, b, beta)
    dv = [i==1 ? c[n] : c[n+1-i] * sp[n+1-i] for i=1:n]
    ev = [-s[n+1-i]*cp[n-i] for i=1:n-1]
    Bidiagonal(dv, ev, true)
end
    

#Return n eigenvalues distributed according to the Jacobi ensemble
function GaussianJacobiSamples(n::Integer, a::Real, b::Real, beta::Real)
    svdvals(GaussianJacobiBidiagonalMatrix(n, a, b))
end
    

