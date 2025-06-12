# Generates samples of dense random matrices that are distributed according
# to the classical Gaussian random matrix ensembles
#
# The notation follows closely that in:
#
# Ioana Dumitriu and Alan Edelman, " Models for Beta s"
# Journal of Mathematical Physics, vol. 43 no. 11 (2002), pp. 5830--5547
# doi: 10.1063/1.1507823
# arXiv: math-ph/0206043
#
# References
#
# Alan Edelman, Per-Olof Persson and Brian D Sutton, "The fourfold way"
# http://www-math.mit.edu/~edelman/homepage/papers/ffw.pdf

export GaussianHermite, GaussianLaguerre, GaussianJacobi,
       Wigner, Wishart, MANOVA #Synonyms for Hermite, Laguerre and Jacobi

#A convenience function to define a chi scalar random variable
χ(df::Real) = rand(Chi(df))

#####################
# Hermite ensemble #
#####################

"""
    GaussianHermite{β} <: ContinuousMatrixDistribution
    GaussianHermite(β::Real) -> GaussianHermite{β}()

Represents a Gaussian Hermite ensemble with Dyson index `β`.

`Wigner{β}` is a synonym.

## Examples

```@example
julia> rand(Wigner(2), 3)
3×3 LinearAlgebra.Hermitian{ComplexF64, Matrix{ComplexF64}}:
   0.383322+0.0im       -0.0452508+0.491032im   -0.313208-0.330435im
 -0.0452508-0.491032im   -0.264521+0.0im        -0.131337+0.0904235im
  -0.313208+0.330435im   -0.131337-0.0904235im  -0.481758+0.0im
```
"""
struct GaussianHermite{β} <: ContinuousMatrixDistribution end
GaussianHermite(β::Real) = GaussianHermite{β}()

"""
Synonym for GaussianHermite{β}
"""
const Wigner{β} = GaussianHermite{β}

"""
    rand(d::Wigner{β}, n::Int)

Generates an `n × n` matrix randomly sampled from the Gaussian-Hermite ensemble (also known as the Wigner ensemble).

The Dyson index `β` is restricted to `β = 1,2` or `4`, for real, complex, and quaternionic fields, respectively.
"""
rand(d::Type{Wigner{β}}, dims...) where {β} = rand(d(), dims...)

function rand(d::Wigner{1}, n::Int)
    A = randn(n, n)
    normalization = 1 / √(2n)
    return Symmetric((A + A') / 2 * normalization)
end

function rand(d::Wigner{2}, n::Int)
    A = randn(ComplexF64, n, n)
    normalization = √(4*n)
    return Hermitian((A + A') / normalization)
end

function rand(d::Wigner{4}, n::Int)
    #Employs 2x2 matrix representation of quaternions
    X = randn(ComplexF64, n, n)
    Y = randn(ComplexF64, n, n)
    A = Matrix{ComplexF64}(undef, 2n, 2n)
    @inbounds for j in 1:n, i in 1:n
        x = X[i, j]
        y = Y[i, j]
        A[i, j] = x
        A[i+n, j] = -conj(y)
        A[i, j+n] = y
        A[i+n, j+n] = conj(x)
    end
    normalization = √(8*n)
    return Hermitian((A + A') / normalization)
end

rand(d::Wigner{β}, n::Int) where {β} =
    throw(ArgumentError("Cannot sample random matrix of size $n x $n for β=$β"))

function rand(d::Wigner{β}, dims::Int...) where {β}
    if length(dims)==2 && dims[1] == dims[2]
	return rand(d, dims[1])
    else
        throw(ArgumentError("Cannot sample random matrix of size $dims for β=$β"))
    end
end

"""
    tridand(d::Wigner{β}, n::Int)

Generates an `n × n` symmetric tridiagonal matrix from the Gaussian-Hermite ensemble (also known as the Wigner ensemble).

Unlike for `rand(Wigner(β), n)`, which is restricted to `β = 1,2` or `4`,
the call `trirand(Wigner(β), n)` will generate a tridiagonal Wigner matrix for any positive
value of `β`, including infinity.

The `β == ∞` case is defined in Edelman, Persson and Sutton, 2012.
"""
function tridrand(d::Wigner{β}, n::Int) where {β}
    χ(df::Real) = rand(Distributions.Chi(df))
    if β≤0
        throw(ArgumentError("β = $β cannot be nonpositive"))
    elseif isinf(β)
        return tridrand(Wigner{Inf}, n)
    else
        normalization = 1 / √(2n)
        Hd = rand(Distributions.Normal(0,2), n)./√2
        He = [χ(β*i)/√2 for i=n-1:-1:1]
        return normalization * SymTridiagonal(Hd, He)
    end
end

function tridrand(d::Wigner{β}, dims...) where {β}
    if length(dims)==2 && dims[1] == dims[2]
	return rand(d, dims[1])
    else
        throw(ArgumentError("Cannot sample random matrix of size $dims for β=$β"))
    end
end

"""
Return n eigenvalues distributed according to the Hermite ensemble
"""
eigvalrand(d::Wigner, n::Int) = eigvals(tridrand(d, n))

"""
Calculate Vandermonde determinant term
"""
function VandermondeDeterminant(λ::AbstractVector{Eigenvalue}, β::Real) where {Eigenvalue<:Number}
    n = length(λ)
    Vandermonde = one(Eigenvalue)^β
    for j=1:n, i=1:j-1
        Vandermonde *= abs(λ[i] - λ[j])^β
    end
    Vandermonde
end

"""
Calculate joint eigenvalue density
"""
function eigvaljpdf(d::Wigner{β}, λ::AbstractVector{Eigenvalue}) where {β,Eigenvalue<:Number}
    n = length(λ)
    #Calculate normalization constant
    c = (2π)^(-n/2)
    for j=1:n
        c *= gamma(1 + β/2)/gamma(1 + β*j/2)
    end
    Energy = sum(λ.^2/2) #Calculate argument of exponential
    VandermondeDeterminant(λ, β) * exp(-Energy)
end

#####################
# Laguerre ensemble #
#####################

"""
    GaussianLaguerre{β}(a::Real)` <: ContinuousMatrixDistribution
    GaussianLaguerre(β::Real, a::Real) -> GaussianLaguerre{β}(a)

Represents a Gaussian-Laguerre ensemble with Dyson index `β` and `a` parameter
used to control the density of eigenvalues near `λ = 0`.

`Wishart{β}(a)` is a synonym.

## Fields
- `a`: Parameter used for weighting the joint probability density function of the ensemble

## Examples
```@example
julia> rand(GaussianLaguerre(1, 2), 2)
2×2 Matrix{Float64}:
  5.08544   -0.156801
 -0.156801   3.17245

julia> rand(GaussianLaguerre(4, 8), 2)
4×4 Matrix{ComplexF64}:
      2.66965+0.0im            0.616044-0.336025im     -8.48419e-18-2.17815e-17im    0.661922-0.190965im
     0.616044+0.336025im        3.27444+0.0im             -0.661922+0.190965im     -1.239e-17+8.84953e-17im
 -8.48419e-18+2.17815e-17im   -0.661922-0.190965im          2.66965+0.0im            0.616044+0.336025im
     0.661922+0.190965im     -1.239e-17-8.84953e-17im      0.616044-0.336025im        3.27444+0.0im
```
## References:
- Edelman and Rao, 2005
"""
struct GaussianLaguerre{β} <: ContinuousMatrixDistribution
  a::Real
end
GaussianLaguerre(β::Real, a::Real) = GaussianLaguerre{β}(a::Real)
const Wishart{β} = GaussianLaguerre{β}

#TODO Check - the eigenvalue distribution looks funky
#TODO The appropriate matrix size should be calculated from a and one matrix dimension
"""
    rand(d::GaussianLaguerre{β}, n::Int)

Generate a random matrix sampled from the Gaussian Laguerre ensemble (also known as the Wishart ensemble)
with parameters defined in `d`.

The Dyson index `β` is restricted to `β = 1,2` (`n × n` matrix) or `4` (`2n × 2n` block matrix representation),
for real, complex, and quaternionic fields, respectively.
"""
function rand(d::GaussianLaguerre{1}, n::Int)
    a = d.a
    a >=  n / 2 || throw(ArgumentError("the minimum value of `a` must be `βn/2`."))
    m = Int(2a)
    A = randn(m, n)
    return (A' * A) / n
end
function rand(d::GaussianLaguerre{2}, n::Int)
    a = d.a
    a >= n || throw(ArgumentError("the minimum value of `a` must be `βn/2`."))
    m = Int(2a)
    A = randn(ComplexF64, m, n)
    return (A' * A) / n
end
function rand(d::GaussianLaguerre{4}, n::Int)
    a = d.a
    a >= 2n || throw(ArgumentError("the minimum value of `a` must be `βn/2`."))
    m = Int(2a)
    # employs 2x2 matrix representation of quaternions
    X = randn(ComplexF64, m, n)
    Y = randn(ComplexF64, m, n)
    A = Matrix{ComplexF64}(undef, 2m, 2n)
    @inbounds for j in 1:n, i in 1:m
        x = X[i, j]
        y = Y[i, j]
        A[i, j] = x
        A[i+m, j] = -conj(y)
        A[i, j+n] = y
        A[i+m, j+n]   = conj(x)
    end
    return (A' * A) / n
end
rand(d::GaussianLaguerre{β}, n::Int) where {β} = throw(ArgumentError("beta = $(β) is not implemented"))

"""
    bidrand(d::GaussianLaguerre{β}, n::Int)

Generate an `n × n` bidiagonal matrix sampled from the Gaussian Laguerre ensemble (also known as the Wishart ensemble).
"""
function bidrand(d::GaussianLaguerre{β}, m::Integer) where {β}
  if d.a <= β*(m-1)/2.0
      error("Given your choice of m and beta, a must be at least $(β*(m-1)/2.0) (You said a = $(d.a))")
  end
  Bidiagonal([chi(2*d.a-i*β) for i=0:m-1], [chi(β*i) for i=m-1:-1:1], true)
end

"""
    tridrand(d::GaussianLaguerre{β}, n::Int)

Generate an `n × n` tridiagonal matrix sampled from the Gaussian Laguerre ensemble (also known as the Wishart ensemble).
"""
function tridrand(d::GaussianLaguerre, m::Integer)
  B = bidrand(d, m)
  L = B * B'
  SymTridiagonal(diag(L)/sqrt(m), diag(L,1)/sqrt(m))
end

#Return n eigenvalues distributed according to the Laguerre ensemble
eigvalrand(d::GaussianLaguerre, m::Integer) = eigvals(tridrand(d, m))

#TODO Check m and ns
function eigvaljpdf(d::GaussianLaguerre{β}, lambda::Vector{Eigenvalue}) where {β,Eigenvalue<:Number}
  m = length(lambda)
  #Laguerre parameters
  p = 0.5*β*(m-1) + 1.0
  #Calculate normalization constant
  c = 2.0^-(m*d.a)
  z = (d.a - β*(m)*0.5)
  for j=1:m
     z += 0.5*β
     if z < 0 && (int(z) - z) < eps()
         #Pole of gamma function, there is no density here no matter what
         return 0.0
     end
     c *= gamma(1 + β/2)/(gamma(1 + β*j/2)*gamma(z))
  end

  Prod = prod(lambda.^(a-p)) #Calculate Laguerre product term
  Energy = sum(lambda)/2 #Calculate argument of exponential
  return c * VandermondeDeterminant(lambda, β) * Prod * exp(-Energy)
end



###################
# Jacobi ensemble #
###################

"""
    GaussianJacobi{β}(a::Real, b::Real) <: ContinuousMatrixDistribution
    GaussianJacobi(β::Real, a::Real, b::Real) -> GaussianJacobi{β}(a, b)

Represents a Gaussian-Jacobi ensemble with Dyson index `β`, while
`a`and `b` are parameters used to weight the joint probability density function of the ensemble.

`MANOVA{β}(a, b)` is a synonym.

## Fields
- `a`: Parameter used for shaping the joint probability density function near `λ = 0`
- `b`: Parameter used for shaping the joint probability density function near `λ = 1`

## References:
- Edelman and Rao, 2005
"""
struct GaussianJacobi{β} <: ContinuousMatrixDistribution
  a::Real
  b::Real
end
const MANOVA{β} = GaussianJacobi{β}

"""
    rand(d::GaussianJacobi{β}, n::Int)

Generate an `n × n` random matrix sampled from the Gaussian-Jacobi ensemble (also known as the MANOVA ensemble)
with parameters defined in `d`.
"""
function rand(d::GaussianJacobi{β}, m::Integer) where {β}
  w1 = Wishart(m, int(2.0*d.a/β), β)
  w2 = Wishart(m, int(2.0*d.b/β), β)
  return (w1 + w2) \ w1
end


# A helper function for Jacobi samples
function SampleCSValues(n::Integer, a::Real, b::Real, beta::Real)
  if beta == Inf
    c =sqrt((a+[1:n])./(a+b+2*[1:n]))
    s =sqrt((b+[1:n])./(a+b+2*[1:n]))
    cp=sqrt([1:n-1]./(a+b+1+2*[1:n-1]))
    sp=sqrt((a+b+1+[1:n-1])./(a+b+1+2*[1:n-1]))
  else
    #Generate cosine-squared values
    csq  = [rand(Beta(beta*(a+i)/2,beta*(b+i)/2)) for i=1:n]
    cpsq = [rand(Beta(beta*i/2,beta*(a+b+1+i)/2)) for i=1:n]
    #Cosine-sine pairs
    c , s  = sqrt(csq) , sqrt(1-csq)
    cp, sp = sqrt(cpsq), sqrt(1-cpsq)
  end
  return c, s, cp, sp
end

#Generates a 2Mx2M sparse MANOVA matrix
#Jacobi ensemble
#
# Reference:
#   A. Edelman and B. D. Sutton, "The beta-Jacobi matrix model, the CS decomposition,
#   and generalized singular value problems", Foundations of Computational Mathematics,
#   vol. 8 iss. 2 (2008), pp 259-285.
#TODO check normalization
function sprand(d::GaussianJacobi{β}, n::Integer, a::Real, b::Real) where {β}
  CoordI = zeros(8n-4)
  CoordJ = zeros(8n-4)
  Values = zeros(8n-4)

  c, s, cp, sp = SampleCSValues(n, a, b, β)

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

#Return n eigenvalues distributed according to the Jacobi ensemble
function eigvalrand(d::GaussianJacobi{β}, n::Integer) where {β}
  #Generate just the upper left quadrant of the matrix
  c, s, cp, sp = SampleCSValues(n, d.a, d.b, β)
  dv = [i==1 ? c[n] : c[n+1-i] * sp[n+1-i] for i=1:n]
  ev = [-s[n+1-i]*cp[n-i] for i=1:n-1]

  ##TODO: understand why dv and ev are returned as Array{Any,1}
  M = Bidiagonal(convert(Array{Float64,1},dv), convert(Array{Float64,1},ev), false)
  return svdvals(M)
end

#TODO Check m and ns
function eigvaljpdf(d::GaussianJacobi{β}, lambda::Vector{Eigenvalue}) where {β,Eigenvalue<:Number}
  m = length(lambda)
  #Jacobi parameters
  a1, a2 = d.a, d.b
  p = 1.0 + β*(m-1)/2.0
  #Calculate normalization constant
  c = 1.0
  for j=1:m
    z1 = (a1 - β*(m-j)/2.0)
    if z1 < 0 && (int(z1) - z1) < eps()
      return 0.0 #Pole of gamma function, there is no density here no matter what
    end
    z2 = (a2 - β*(m-j)/2.0)
    if z2 < 0 && (int(z2) - z2) < eps()
      return 0.0 #Pole of gamma function, there is no density here no matter what
    end
    c *= gamma(1 + β/2)*gamma(a1+a2-β*(m-j)/2)
    c /= gamma(1 + β*j/2)*gamma(z1)*gamma(z2)
  end

  Prod = prod(lambda.^(a1-p))*prod((1-lambda).^(a2-p)) #Calculate Laguerre product term
  Energy = sum(lambda/2) #Calculate argument of exponential

  return c * VandermondeDeterminant(lambda, β) * Prod * exp(-Energy)
end
