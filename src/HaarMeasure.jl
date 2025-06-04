#TODO implement O(n^2) method
"""
    rand(W::Haar, n::Int)

Computes samples of real or complex Haar matrices of size `n`×`n`.

For `beta = 1,2,4`, generates random orthogonal, unitary and symplectic matrices
of uniform Haar measure.
These matrices are distributed with uniform Haar measure over the
classical orthogonal, unitary and symplectic groups `O(n)`, `U(n)` and
`Sp(n)~USp(2n)` respectively.

    rand(W::Haar, n::Int, doCorrection::Int = 1)

The additional argument `doCorrection` specifies whether or not the piecewise correction
is applied to ensure that it is truly of Haar measure.
This addresses an inconsistency in the Householder reflections as
implemented in most versions of LAPACK.
- Method 0: No correction
- Method 1: Multiply rows by uniform random phases
- Method 2: Multiply rows by phases of diag(R)

## References:
- Edelman and Rao, 2005
- Mezzadri, 2006, math-ph/0609050
"""
function rand(W::Haar, n::Int, doCorrection::Int=1)
    beta = W.beta
    M=rand(Ginibre(beta), n)
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
        q*Diagonal(L)
    elseif doCorrection==2
        if beta==1
            L=sign.(diag(r))
        elseif (beta==2 || beta==4)
            L=diag(r)
            L=L./abs.(L)
        else
            error(string("beta = ",beta, " not implemented."))
        end
        q*Diagonal(L)
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


import Base.size
import LinearAlgebra: lmul!, rmul!, QRPackedQ, Diagonal, Adjoint



struct StewartQ{T,S<:AbstractMatrix{T},C<:AbstractVector{T},D<:AbstractVector{T}} <: LinearAlgebra.AbstractQ{T}
    factors::S
    τ::C
    signs::D
end
size(Q::StewartQ, dim::Integer) = size(Q.factors, dim == 2 ? 1 : dim)
size(Q::StewartQ) = (n = size(Q.factors, 1); (n, n))

lmul!(A::StewartQ, B::AbstractVecOrMat) = lmul!(QRPackedQ(A.factors,A.τ), lmul!(Diagonal(A.signs), B))
rmul!(A::AbstractVecOrMat, B::StewartQ) = rmul!(rmul!(A, QRPackedQ(B.factors,B.τ)), Diagonal(B.signs))

@static if VERSION ≥ v"1.10"
    import LinearAlgebra.AdjointQ
    lmul!(adjA::AdjointQ{<:Any,<:StewartQ}, B::AbstractVecOrMat) =
        lmul!(Diagonal(adjA.Q.signs), lmul!(QRPackedQ(adjA.Q.factors, adjA.Q.τ)', B))
    rmul!(A::AbstractVecOrMat, adjB::AdjointQ{<:Any,<:StewartQ}) =
        rmul!(rmul!(A, Diagonal(adjB.Q.signs)), QRPackedQ(adjB.Q.factors, adjB.Q.τ)')
    else
    lmul!(adjA::Adjoint{<:Any,<:StewartQ}, B::AbstractVecOrMat) =
        lmul!(Diagonal(adjA.parent.signs), lmul!(QRPackedQ(adjA.parent.factors, adjA.parent.τ)', B))
    rmul!(A::AbstractVecOrMat, adjB::Adjoint{<:Any,<:StewartQ}) =
        rmul!(rmul!(A, Diagonal(adjB.parent.signs)), QRPackedQ(adjB.parent.factors,adjB.parent.τ)')
end

StewartQ{T}(Q::StewartQ) where {T} = StewartQ(convert(AbstractMatrix{T}, Q.factors), convert(Vector{T}, Q.τ), convert(Vector{T}, Q.signs))
AbstractMatrix{T}(Q::StewartQ{T}) where {T} = Q
AbstractMatrix{T}(Q::StewartQ) where {T} = StewartQ{T}(Q)

"""
Stewart's algorithm for sampling orthogonal/unitary random matrices in time O(n^2)
"""
function Stewart(::Type{T}, n) where {T<:Union{Float64,ComplexF64}}
    τ = Array{T}(undef, n)
    signs = Vector{T}(undef, n)
    H = randn(T, n, n)

    pτ = pointer(τ)
    pβ = pointer(H)
    pH = pointer(H, 2)

    incr = (T <: Real ? 1 : 2)
    for i = 0:n-1
        larfg!(n - i, pβ + (n + 1)*8*incr*i, pH + (n + 1)*8*incr*i, 1, pτ + 8*incr*i)
        signs[i+1] = sign(real(H[i+1, i+1]))
    end
    return StewartQ(H, τ, signs)
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
