export NeedPiecewiseCorrection


# Computes samplex of real or complex Haar matrices of size nxn
#
# For beta=1,2,4, generates random orthogonal, unitary and symplectic matrices
# of uniform Haar measure. 
# These matrices are distributed with uniform Haar measure over the
# classical orthogonal, unitary and symplectic groups O(N), U(N) and
# Sp(N)~USp(2N) respectively.
#
# The last parameter specifies whether or not the piecewise correction
# is applied to ensure that it truly of Haar measure
# This addresses an inconsistency in the Householder reflections as
# implemented in most versions of LAPACK
# Method 0: No correction
# Method 1: Multiply rows by uniform random phases
# Method 2: Multiply rows by phases of diag(R)
# References:
#    Edelman and Rao, 2005
#    Mezzadri, 2006, math-ph/0609050
#TODO implement O(n^2) method
function rand(W::Haar, n::Integer, doCorrection::Integer)
    beta = W.beta
    M=rand(Ginibre(beta,n))
    q,r=qr(M)
    if doCorrection==0
        q
    elseif doCorrection==1 
        if beta==1
            L = sign(rand(n).-0.5)
        elseif beta==2
            L = exp(im*rand(n)*2pi)
        elseif beta==4
            L = exp(im*rand(2n)*2pi)
        else
            error(string("beta = ",beta, " not implemented."))
        end
        q*diagm(L)
    elseif doCorrection==2
        if beta==1
            L=sign(diag(r))
        elseif (beta==2 || beta==4)
            L=diag(r)
            L=L./abs(L)
        else
            error(string("beta = ",beta, " not implemented."))
        end
        q*diagm(L)
    end
end

#By default, always do piecewise correction
#For most applications where you use the HaarMatrix as a similarity transform
#it doesn't matter, but better safe than sorry... let the user choose else
rand(W::Haar,n::Integer) = rand(W,n, 1)

#A utility method to check if the piecewise correction is needed
#This checks the R part of the QR factorization; if correctly done,
#the diagonals are all chi variables so are non-negative
function NeedPiecewiseCorrection()
    n=20
    R=qr(randn(Ginibre(2,n)))[2]
    return any([x<0 for x in diag(R)])
end


#TODO maybe, someday
#Haar measure on U(N) in terms of local coordinates
#Zyczkowski and Kus, Random unitary matrices, J. Phys. A: Math. Gen. 27,
#4235–4245 (1994).

