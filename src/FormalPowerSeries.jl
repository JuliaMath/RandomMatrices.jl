# Excursions in formal power series (fps)
#
# Jiahao Chen <jiahao@mit.edu> 2013
#
# The primary reference is
# [H] Peter Henrici, "Applied and Computational Complex Analysis", Volume I:
#     Power Series---Integration---Conformal Mapping---Location of Zeros,
#     Wiley-Interscience: New York, 1974

import Base: zero, one, inv, length, ==, +, -, *, ^
import Base.Broadcast: broadcasted
export FormalPowerSeries, fps, tovector, trim, isunit, isnonunit,
       MatrixForm, reciprocal, derivative, isconstant, compose,
       isalmostunit, FormalLaurentSeries

#############################
# Formal power series (fps) #
#############################

# Definition [H, p.9] - we limit this to finitely many coefficients
mutable struct FormalPowerSeries{Coefficients}
  c :: Dict{BigInt, Coefficients}
  FormalPowerSeries{Coefficients}(c::Dict{BigInt, Coefficients}) where Coefficients = new(c)
  function FormalPowerSeries{Coefficients}(v::Vector{Coefficients}) where Coefficients
    c=Dict{BigInt, Coefficients}()
    for i in eachindex(v)
      if v[i] != 0
        c[i-1] = v[i] #Note this off by 1 - allows constant term c[0] to be set
      end
    end
    FormalPowerSeries{Coefficients}(c)
  end
end

FormalPowerSeries(v::Vector{T}) where T = FormalPowerSeries{T}(v)

#Convenient abbreviation for floats
fps = FormalPowerSeries{Float64}

zero(P::FormalPowerSeries{T}) where {T} = FormalPowerSeries(T[])
one(P::FormalPowerSeries{T}) where {T} = FormalPowerSeries(T[1])

#Return truncated vector with c[i] = P[n[i]]
function tovector(P::FormalPowerSeries{T}, n::AbstractVector{Index}) where {T,Index<:Integer}
  nn = length(n)
  c = zeros(nn)
  for (k,v) in P.c, i in eachindex(n)
      n[i]==k && (c[i]=v)
  end
  c
end

tovector(P::FormalPowerSeries, n::Integer)=tovector(P,1:n)


# Basic housekeeping and properties

# Remove extraneous zeros
function trim(P::FormalPowerSeries{T}) where {T}
  for (k,v) in P.c
    if v==0
      delete!(P.c, k)
    end
  end
  return P
end

length(P::FormalPowerSeries{T}) where {T} = maximum([k for (k,v) in P.c])

function ==(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T}) where {T}
  for (k,v) in P.c
    if v==0 #ignore explicit zeros
      continue
    elseif !haskey(Q.c, k)
      return false
    elseif Q.c[k] != v
      return false
    end
  end
  for (k,v) in Q.c
    if v==0 #ignore explicit zeros
      continue
    elseif !haskey(P.c, k)
      return false
    elseif P.c[k] != v
      return false
    end
  end
  return true
end

# Basic arithmetic [H, p.10]
function +(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T}) where {T}
  c = Dict{BigInt, T}()
  for (k,v) in P.c
    haskey(c,k) ? (c[k]+=v) : (c[k]=v)
  end
  for (k,v) in Q.c
    haskey(c,k) ? (c[k]+=v) : (c[k]=v)
  end
  FormalPowerSeries{T}(c)
end

function -(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T}) where {T}
  c = Dict{BigInt, T}()
  for (k,v) in P.c
    c[k] = get(c,k,zero(T)) + v
  end
  for (k,v) in Q.c
    c[k] = get(c,k,zero(T)) - v
  end
  FormalPowerSeries{T}(c)
end

#negation
function -(P::FormalPowerSeries{T}) where {T}
  c = Dict{BigInt, T}()
  for (k,v) in P.c
    c[k] = -v
  end
  FormalPowerSeries{T}(c)
end

#multiplication by scalar
function *(P::FormalPowerSeries{T}, n::Number) where {T}
  c = Dict{BigInt, T}()
  for (k,v) in P.c
    c[k] = n*v
  end
  FormalPowerSeries{T}(c)
end
*(n::Number, P::FormalPowerSeries{T}) where {T} = *(P, n)

#Cauchy product [H, p.10]
*(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T}) where {T} = CauchyProduct(P, Q)
function CauchyProduct(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T}) where {T}
  c = Dict{BigInt, T}()
  for (k1, v1) in P.c, (k2, v2) in Q.c
      c[k1+k2]=get(c, k1+k2, zero(T))+v1*v2
  end
  FormalPowerSeries{T}(c)
end

#Hadamard product [H, p.10] - the elementwise product
broadcasted(::typeof(*), P::FormalPowerSeries{T}, Q::FormalPowerSeries{T}) where T =
  HadamardProduct(P, Q)
function HadamardProduct(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T}) where {T}
  c = Dict{BigInt, T}()
  for (k,v) in P.c
    if v!=0 && haskey(Q.c,k) && Q.c[k] !=0
      c[k] = v * Q.c[k]
    end
  end
  FormalPowerSeries{T}(c)
end

#The identity element over the complex numbers
#replacement for previous function eye(P::FormalPowerSeries{T})
function FormalPowerSeries{T}(s::UniformScaling) where {T}
  c = Dict{BigInt, T}()
  c[0] = 1
  return FormalPowerSeries{T}(c)
end

isunit(P::FormalPowerSeries{T}) where {T <: Number} = P==FormalPowerSeries{Float64}(I)

# [H, p.12]
isnonunit(P::FormalPowerSeries{T}) where {T} = (!haskey(P.c, 0) || P.c[0]==0) && !isunit(P)

#Constructs the top left m x m block of the (infinite) semicirculant matrix
#associated with the fps [H, Sec.1.3, p.14]
#[H] calls it the semicirculant, but in contemporary nomenclature this is an
#upper triangular Toeplitz matrix
#This constructs the dense matrix - Toeplitz matrices don't exist in Julia yet
function MatrixForm(P::FormalPowerSeries{T}, m :: Integer) where {T}
  m<0 && error("Invalid matrix dimension $m requested")
  M=zeros(T, m, m)
  for (k,v) in P.c
    if k < m
      for i=1:m-k
        M[i, i+k]=v
      end
    end
  end
  M
end

#Reciprocal of power series [H, Sec.1.3, p.15]
#
#This algorithm is NOT that in [H], since it uses determinantal inversion
#formulae for the inverse. It is much quicker to find the inverse of the
#associated infinite triangular Toeplitz matrix!
#
#Every fps is isomorphic to such a triangular Toeplitz matrix [H. Sec.1.3, p.15]
#
#As most reciprocals of finite series are actually infinite,
#allow the inverse to have a finite truncation
#
#TODO implement the FFT-based algorithm of one of the following
#  doi:10.1109/TAC.1984.1103499
#  https://cs.uwaterloo.ca/~glabahn/Papers/tamir-george-1.pdf
function reciprocal(P::FormalPowerSeries, n::Integer)
  n<0 && error(sprintf("Invalid inverse truncation length %d requested", n))

  a = tovector(P, 0:n-1) #Extract original coefficients in vector
  a[1]==0 && error("Inverse does not exist")
  inv_a1 = inv(a[1])

  b = zeros(n) #Inverse
  b[1] = inv_a1
  for k=2:n
    b[k] = -inv_a1 * sum([b[k+1-i] * a[i] for i=2:k])
  end
  FormalPowerSeries(b)
end

#Derivative of fps [H. Sec.1.4, p.18]
function derivative(P::FormalPowerSeries{T}) where T
    c = Dict{BigInt, T}()
    for (k,v) in P.c
        if k != 0 && v != 0
            c[k-1] = k*v
        end
    end
    FormalPowerSeries{T}(c)
end

#[H, Sec.1.4, p.18]
function isconstant(P::FormalPowerSeries)
  for (k,v) in P.c
    if k!=0 && v!=0
      return false
    end
  end
  return true
end



# Power of fps [H, Sec.1.5, p.35]
function ^(P::FormalPowerSeries{T}, n::Integer) where T
    if n == 0
        return FormalPowerSeries{T}([one(T)])
    elseif n > 1
        #The straightforward recursive way
        return P^(n-1) * P
        #TODO implement the non-recursive formula
    elseif n==1
        return P
    else
        error("I don't know what to do for power n = $n")
    end
end



# Composition of fps [H, Sec.1.5, p.35]
# This is the quick and dirty version
function compose(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T}) where {T}
  sum([v * Q^k for (k, v) in P.c])
end


#[H, p.45]
function isalmostunit(P::FormalPowerSeries{T}) where {T}
  (haskey(P.c, 0) && P.c[0]!=0) ? (return false) :
  (haskey(P.c, 1) && P.c[1]!=0) ? (return true) : (return false)
end


# Reversion of a series (inverse with respect to composition)
# P^[-1]
# [H. Sec.1.7, p.47, but more succinctly stated on p.55]
# Constructs the upper left nxn subblock of the matrix representation
# and inverts it
function reversion(P::FormalPowerSeries{T}, n :: Integer) where {T}
  n>0 ? error("Need non-negative dimension") : nothing

  Q = P
  A = zeros(n, n)
  #Construct the matrix representation (1.9-1), p.55
  for i=1:n
      Q = P
      ai = tovector(Q, n) #Extract coefficients P[1]...P[n]
      A[i,:] = ai
      i<n && (Q *= P)
  end

  #TODO I just need the first row of the inverse
  B = inv(A)
  FormalPowerSeries{T}(B[1, :]'[:,1])
end

inv(P::FormalPowerSeries, n::Integer) = reversion(P, n)


###############################
# Formal Laurent series (fLs) #
# [H., Sec. 1.8, p. 51]       #
###############################

mutable struct FormalLaurentSeries{Coefficients}
  c :: Dict{BigInt, Coefficients}
  FormalLaurentSeries{Coefficients}(c::Dict{BigInt, Coefficients}) where Coefficients = new(c)
  function FormalLaurentSeries{Coefficients}(v::Vector{Coefficients}) where Coefficients
    c = new(c)
    for i=1:length(v)
      c[i] = v[i]
    end
    c
  end
end
