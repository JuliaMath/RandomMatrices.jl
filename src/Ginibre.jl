export rand, Ginibre
import Base.rand

#Samples a matrix from the Ginibre ensemble
#This ensemble lives in GL(N, F), the set of all invertible N x N matrices
#over the field F
#For beta=1,2,4, F=R, C, H respectively
struct Ginibre <: ContinuousMatrixDistribution
   beta::Float64
   N::Integer
end

function rand(W::Ginibre)
    beta, n = W.beta, W.N
    if beta==1
        randn(n,n)
    elseif beta==2
        randn(n,n)+im*randn(n,n)
    elseif beta==4
        Q0=randn(n,n)
        Q1=randn(n,n)
        Q2=randn(n,n)
        Q3=randn(n,n)
        [Q0+im*Q1 Q2+im*Q3;-Q2+im*Q3 Q0-im*Q1]
    else 
        error(string("beta = ", beta, " not implemented"))
    end
end

function jpdf(Z::AbstractMatrix{z}) where {z<:Complex}
    pi^(size(Z,1)^2)*exp(-trace(Z'*Z))
end

#Dyson Circular orthogonal, unitary and symplectic ensembles
struct CircularOrthogonal
    N :: Int64
end

function rand(M::CircularOrthogonal)
    U = rand(Ginibre(2, M.N))
    U * U'
end

struct CircularUnitary
    N :: Int64
end

rand(M::CircularUnitary) = rand(Ginibre(2, M.N))

struct CircularSymplectic
    N :: Int64
end

rand(M::CircularSymplectic) = error("Not implemented")

