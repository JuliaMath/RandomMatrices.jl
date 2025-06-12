export rand, Ginibre
import Base.rand

"""
    Ginibre{β} <: ContinuousMatrixDistribution
    Ginibre(β::Real) -> Ginibre{β}()

Represents a Ginibre ensemble with Dyson index `β` living in `GL(N, F)`, the set
of all invertible `N × N` matrices over the field `F`. 

## Examples

```@example
julia> rand(Ginibre(2), 3)
3×3 Matrix{ComplexF64}:
 0.781329+2.00346im   0.0595122+0.488652im  -0.323494-0.35966im
  1.11089+0.935174im  -0.384457+1.71419im    0.114358-0.360676im
  1.54119+0.362003im  -0.693623-2.50141im    -1.42383-1.06341im
```

## References:
- Edelman and Rao, 2005
"""
struct Ginibre{β} <: ContinuousMatrixDistribution end
Ginibre(β::B) where {B} = Ginibre{β}()

"""
    rand(W::Ginibre{β}, n::Int)

Samples a matrix from the Ginibre ensemble.

For `β = 1,2,4`, generates matrices randomly sampled from the real, complex, and quaternion
Ginibre ensemble, respectively.
"""
rand(W::Ginibre{1}, n::Int) = randn(n, n)
rand(W::Ginibre{2}, n::Int) = randn(ComplexF64, n, n)
function rand(W::Ginibre{4}, n::Int)
    Q0=randn(n,n)
    Q1=randn(n,n)
    Q2=randn(n,n)
    Q3=randn(n,n)
    return [Q0+im*Q1 Q2+im*Q3;-Q2+im*Q3 Q0-im*Q1]
end
rand(W::Ginibre{β}, n::Int) where {β} = throw(ArgumentError("Cannot sample random matrix of size $n x $n for β=$β"))


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

