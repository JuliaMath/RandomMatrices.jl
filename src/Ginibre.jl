export rand, Ginibre
import Base.rand

"""
    Ginibre(β::Int, N::Int) <: ContinuousMatrixDistribution

Represents a Ginibre ensemble with Dyson index `β` living in `GL(N, F)`, the set
of all invertible `N × N` matrices over the field `F`. 

## Fields
- `beta`: Dyson index
- `N`: Matrix dimension over the field `F`.

## Examples

```jldoctest
julia> Random.seed!(1234);

julia> rand(Ginibre(2, 3))
3×3 Matrix{ComplexF64}:
  0.970656-0.763689im  -0.0328031-0.0998909im    2.70742+0.942733im
 -0.979218-0.534709im   -0.600792-0.726142im     1.52445-0.00991272im
  0.901861-0.837116im    -1.44518-0.00420647im  -0.20563-0.66748im
```

## References:
- Edelman and Rao, 2005
"""
struct Ginibre <: ContinuousMatrixDistribution
   beta::Float64
   N::Integer
end

"""
    rand(W::Ginibre)

Samples a matrix from the Ginibre ensemble.

For `β = 1,2,4`, generates matrices randomly sampled from the real, complex, and quaternion
Ginibre ensemble, respectively.
"""
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

