# Computes samples of real or complex Haar matrices of size nxn
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
#By default, always do piecewise correction
#For most applications where you use the HaarMatrix as a similarity transform
#it doesn't matter, but better safe than sorry... let the user choose else
function rand(W::Haar, n::Int, doCorrection::Int=1)
    beta = W.beta
    M=rand(Ginibre(beta,n))
    q,r=qr(M)
    if doCorrection==0
        q
    elseif doCorrection==1
        if beta==1
            L = sign.(rand(n).-0.5)
        elseif beta==2
            L = exp.(im*rand(n)*2pi)
        elseif beta==4
            L = exp.(im*rand(2n)*2pi)
        else
            error(string("beta = ",beta, " not implemented."))
        end
        q*diagm(0 => L)
    elseif doCorrection==2
        if beta==1
            L=sign.(diag(r))
        elseif (beta==2 || beta==4)
            L=diag(r)
            L=L./abs.(L)
        else
            error(string("beta = ",beta, " not implemented."))
        end
        q*diagm(0 => L)
    end
end

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


import LinearAlgebra.BLAS: @blasfunc


using LinearAlgebra: BlasInt
for (s, elty) in (("dlarfg_", Float64),
                  ("zlarfg_", ComplexF64))
    @eval begin
        function larfg!(n::Int, α::Ptr{$elty}, x::Ptr{$elty}, incx::Int, τ::Ptr{$elty})
	    ccall((@blasfunc($s), LAPACK.liblapack), Nothing,
		  (Ptr{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty}),
                  Ref(n), α, x, Ref(incx), τ)
        end
    end
end

"""
Stewarts algorithm for n^2 orthogonal random matrix
"""
function Stewart(::Type{Float64}, n)
    τ = Array{Float64}(undef, n)
    H = randn(n, n)

    pτ = pointer(τ)
    pβ = pointer(H)
    pH = pointer(H, 2)

    for i = 0:n-2
        larfg!(n - i, pβ + (n + 1)*8i, pH + (n + 1)*8i, 1, pτ + 8i)
    end
    LinearAlgebra.QRPackedQ(H,τ)
end

function Stewart(::Type{ComplexF64}, n)
    τ = Array{ComplexF64}(undef, n)
    H = complex.(randn(n, n), randn(n,n))

    pτ = pointer(τ)
    pβ = pointer(H)
    pH = pointer(H, 2)

    for i = 0:n-2
        larfg!(n - i, pβ + (n + 1)*16i, pH + (n + 1)*16i, 1, pτ + 16i)
    end
    LinearAlgebra.QRPackedQ(H,τ)
end

export randfast
function randfast(W::Haar, n::Int)
    if W.beta==1
        Stewart(Float64, n)
    elseif W.beta==2
        Stewart(ComplexF64, n)
    else
        error("beta = $beta not implemented")
    end
end
