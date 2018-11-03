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
GaussianHermite{β} represents a Gaussian Hermite ensemble with parameter β.

Wigner{β} is a synonym.

Example of usage:

    β = 2 #β = 1, 2, 4 generates real, complex and quaternionic matrices respectively.
    d = Wigner{β} #same as GaussianHermite{β}

    n = 20 #Generate square matrices of this size

    S = rand(d, n)       #Generate a 20x20 symmetric Wigner matrix
    T = tridrand(d, n)   #Generate the symmetric tridiagonal form
    v = eigvalrand(d, n) #Generate a sample of eigenvalues
"""
struct GaussianHermite{β} <: ContinuousMatrixDistribution end
GaussianHermite(β) = GaussianHermite{β}()

"""
Synonym for GaussianHermite{β}
"""
const Wigner{β} = GaussianHermite{β}

rand(d::Type{Wigner{β}}, dims...) where {β} = rand(d(), dims...)

function rand(d::Wigner{1}, n::Int)
    A = randn(n, n)
    normalization = 1 / √(2n)
    return Symmetric((A + A') / 2 * normalization)
end

function rand(d::Wigner{2}, n::Int)
    A = randn(n, n) + im*randn(n, n)
    normalization = √(4*n)
    return Hermitian((A + A') / normalization)
end

function rand(d::Wigner{4}, n::Int)
    #Employs 2x2 matrix representation of quaternions
    X = randn(n, n) + im*randn(n, n)
    Y = randn(n, n) + im*randn(n, n)
    A = [X Y; -conj(Y) conj(X)]
    normalization = √(8*n)
    return Hermitian((A + A') / normalization)
end

rand(d::Wigner{β}, n::Int) where {β} =
    throw(ValueError("Cannot sample random matrix of size $n x $n for β=$β"))

function rand(d::Wigner{β}, dims...) where {β}
    if length(dims)==2 && dims[1] == dims[2]
	return rand(d, dims[1])
    else
        throw(ValueError("Cannot sample random matrix of size $dims for β=$β"))
    end
end

"""
Generates a nxn symmetric tridiagonal Wigner matrix

Unlike for `rand(Wigner{β}, n)`, which is restricted to β=1,2 or 4,
`trirand(Wigner{β}, n)` will generate a

The β=∞ case is defined in Edelman, Persson and Sutton, 2012
"""
function tridrand(d::Wigner{β}, n::Int) where {β}
    χ(df::Real) = rand(Distributions.Chi(df))
    if β≤0
        throw(ValueError("β = $β cannot be nonpositive"))
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
        throw(ValueError("Cannot sample random matrix of size $dims for β=$β"))
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

mutable struct GaussianLaguerre <: ContinuousMatrixDistribution
  beta::Real
  a::Real
end
const Wishart = GaussianLaguerre


#Generates a NxN Hermitian Wishart matrix
#TODO Check - the eigenvalue distribution looks funky
#TODO The appropriate matrix size should be calculated from a and one matrix dimension
function rand(d::GaussianLaguerre, dims::Dim2)
  n = 2.0*a/d.beta
  if d.beta == 1 #real
    A = randn(dims)
  elseif d.beta == 2 #complex
    A = randn(dims) + im*randn(dims)
  elseif d.beta == 4 #quaternion
    #Employs 2x2 matrix representation of quaternions
    X = randn(dims) + im*randn(dims)
    Y = randn(dims) + im*randn(dims)
    A = [X Y; -conj(Y) conj(X)]
    error("beta = $(d.beta) is not implemented")
  end
  return (A * A') / dims[1]
end

#Generates a NxN bidiagonal Wishart matrix
function bidrand(d::GaussianLaguerre, m::Integer)
  if d.a <= d.beta*(m-1)/2.0
      error("Given your choice of m and beta, a must be at least $(d.beta*(m-1)/2.0) (You said a = $(d.a))")
  end
  Bidiagonal([chi(2*d.a-i*d.beta) for i=0:m-1], [chi(d.beta*i) for i=m-1:-1:1], true)
end

#Generates a NxN tridiagonal Wishart matrix
function tridrand(d::GaussianLaguerre, m::Integer)
  B = bidrand(d, m)
  L = B * B'
  SymTridiagonal(diag(L)/sqrt(m), diag(L,1)/sqrt(m))
end

#Return n eigenvalues distributed according to the Laguerre ensemble
eigvalrand(d::GaussianLaguerre, m::Integer) = eigvals(tridrand(d, m))

#TODO Check m and ns
function eigvaljpdf(d::GaussianLaguerre, lambda::Vector{Eigenvalue}) where {Eigenvalue<:Number}
  m = length(lambda)
  #Laguerre parameters
  p = 0.5*d.beta*(m-1) + 1.0
  #Calculate normalization constant
  c = 2.0^-(m*d.a)
  z = (d.a - d.beta*(m)*0.5)
  for j=1:m
     z += 0.5*d.beta
     if z < 0 && (int(z) - z) < eps()
         #Pole of gamma function, there is no density here no matter what
         return 0.0
     end
     c *= gamma(1 + beta/2)/(gamma(1 + beta*j/2)*gamma(z))
  end

  Prod = prod(lambda.^(a-p)) #Calculate Laguerre product term
  Energy = sum(lambda)/2 #Calculate argument of exponential
  return c * VandermondeDeterminant(lambda, beta) * Prod * exp(-Energy)
end



###################
# Jacobi ensemble #
###################

#Generates a NxN self-dual MANOVA matrix
mutable struct GaussianJacobi <: ContinuousMatrixDistribution
  beta::Real
  a::Real
  b::Real
end
const MANOVA = GaussianJacobi

function rand(d::GaussianJacobi, m::Integer)
  w1 = Wishart(m, int(2.0*d.a/d.beta), d.beta)
  w2 = Wishart(m, int(2.0*d.b/d.beta), d.beta)
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
function sprand(d::GaussianJacobi, n::Integer, a::Real, b::Real)
  CoordI = zeros(8n-4)
  CoordJ = zeros(8n-4)
  Values = zeros(8n-4)

  c, s, cp, sp = SampleCSValues(n, a, b, d.beta)

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
function eigvalrand(d::GaussianJacobi, n::Integer)
  #Generate just the upper left quadrant of the matrix
  c, s, cp, sp = SampleCSValues(n, d.a, d.b, d.beta)
  dv = [i==1 ? c[n] : c[n+1-i] * sp[n+1-i] for i=1:n]
  ev = [-s[n+1-i]*cp[n-i] for i=1:n-1]

  ##TODO: understand why dv and ev are returned as Array{Any,1}
  M = Bidiagonal(convert(Array{Float64,1},dv), convert(Array{Float64,1},ev), false)
  return svdvals(M)
end

#TODO Check m and ns
function eigvaljpdf(d::GaussianJacobi, lambda::Vector{Eigenvalue}) where {Eigenvalue<:Number}
  m = length(lambda)
  #Jacobi parameters
  a1, a2 = d.a, d.b
  p = 1.0 + d.beta*(m-1)/2.0
  #Calculate normalization constant
  c = 1.0
  for j=1:m
    z1 = (a1 - beta*(m-j)/2.0)
    if z1 < 0 && (int(z1) - z1) < eps()
      return 0.0 #Pole of gamma function, there is no density here no matter what
    end
    z2 = (a2 - beta*(m-j)/2.0)
    if z2 < 0 && (int(z2) - z2) < eps()
      return 0.0 #Pole of gamma function, there is no density here no matter what
    end
    c *= gamma(1 + beta/2)*gamma(a1+a2-beta*(m-j)/2)
    c /= gamma(1 + beta*j/2)*gamma(z1)*gamma(z2)
  end

  Prod = prod(lambda.^(a1-p))*prod((1-lambda).^(a2-p)) #Calculate Laguerre product term
  Energy = sum(lambda/2) #Calculate argument of exponential

  return c * VandermondeDeterminant(lambda, beta) * Prod * exp(-Energy)
end
