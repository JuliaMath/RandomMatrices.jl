#!/usr/bin/env julia
#
# Excursions in formal power series (fps)
#
# Jiahao Chen <jiahao@mit.edu> 2013
#
# The primary reference is
# [H] Peter Henrici, "Applied and Computational Complex Analysis", Volume I:
#     Power Series---Integration---Conformal Mapping---Location of Zeros,
#     Wiley-Interscience: New York, 1974

import Base.inv, Base.length

#############################
# Formal power series (fps) #
#############################

# Definition [H, p.9] - we limit this to finitely many coefficients
type FormalPowerSeries{Coefficients}
    c :: Dict{BigInt, Coefficients}
    FormalPowerSeries{Coefficients}(c::Dict{BigInt, Coefficients}) = new(c)
    function FormalPowerSeries{Coefficients}(v::Vector{Coefficients})
        c=Dict{BigInt, Coefficients}()
        for i=1:length(v)
            if v[i] != 0
                c[i-1] = v[i] #Note this off by 1 - allows constant term c[0] to be set
            end     
        end
        FormalPowerSeries{Coefficients}(c)
    end
end



#Return truncated vector with c[i+1] = P[i]
function tovector{T}(P::FormalPowerSeries{T}, n :: Integer)
    c = zeros(n)
    for (k,v) in P.c
        if 1<=k+1<=n
            c[k+1]=v
        end
    end
    c
end


function =={T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T})
    for (k,v) in P.c
        if v==0 #ignore explicit zeros
            continue
        elseif !has(Q.c, k)
            return false
        elseif Q.c[k] != v
            return false
        end
    end
    for (k,v) in Q.c
        if v==0 #ignore explicit zeros
            continue
        elseif !has(P.c, k)
            return false
        elseif P.c[k] != v
            return false
        end
    end
    return true
end



# Basic housekeeping and properties

# Remove extraneous zeros
function trim{T}(P::FormalPowerSeries{T})
    for (k,v) in P.c
        if v==0
            delete!(P.c, k)
        end
    end
    return P
end

length{T}(P::FormalPowerSeries{T})=max([k for (k,v) in P.c])

# Basic arithmetic [H, p.10]
function +{T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T})
    c = Dict{BigInt, T}()
    for (k,v) in P.c
        has(c,k) ? (c[k]+=v) : (c[k]=v) 
    end
    for (k,v) in Q.c
        has(c,k) ? (c[k]+=v) : (c[k]=v) 
    end
    FormalPowerSeries{T}(c)
end

function -{T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T})
    c = Dict{BigInt, T}()
    for (k,v) in P.c
        has(c,k) ? (c[k]+=v) : (c[k]=v) 
    end
    for (k,v) in Q.c
        has(c,k) ? (c[k]-=v) : (c[k]=-v) 
    end
    FormalPowerSeries{T}(c)
end

#negation
function -{T}(P::FormalPowerSeries{T})
    c = Dict{BigInt, T}()
    for (k,v) in P.c
        c[k] = -v
    end
    FormalPowerSeries{T}(c)
end

#multiplication by scalar
function *{T}(P::FormalPowerSeries{T}, n::Number)
    c = Dict{BigInt, T}()
    for (k,v) in P.c
        c[k] = n*v
    end
    FormalPowerSeries{T}(c)
end
*{T}(n::Number, P::FormalPowerSeries{T}) = *(P, n)

#Cauchy product [H, p.10]
*{T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T}) = CauchyProduct(P, Q)
function CauchyProduct{T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T})
    c = Dict{BigInt, T}()
    for (k1, v1) in P.c
        for (k2, v2) in Q.c
            has(c, k1+k2) ? (c[k1+k2]+=v1*v2) : (c[k1+k2]=v1*v2)
        end
    end
    FormalPowerSeries{T}(c)
end

#Hadamard product [H, p.10] - the elementwise product
.*{T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T}) = HadamardProduct(P, Q)
function HadamardProduct{T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T})
    c = Dict{BigInt, T}()
    for (k,v) in P.c
        if v!=0 && has(Q.c,k) && Q.c[k]==0
            c[k] = v * Q.c[k]
        end
    end
    FormalPowerSeries{T}(c)
end

#The identity element over the complex numbers
function eye{T <: Number}(P::FormalPowerSeries{T})
    c = Dict{BigInt, T}()
    c[0] = 1
    FormalPowerSeries{T}(c)
end

isunit{T <: Number}(P::FormalPowerSeries{T}) = P==eye(P)

# [H, p.12]
isnonunit{T}(P::FormalPowerSeries{T}) = (!has(P.c, 0) || P.c[0]==0) && !isunit(P)

#Constructs the top left m x m block of the (infinite) semicirculant matrix
#associated with the fps [H, Sec.1.3, p.14]
#[H] calls it the semicirculant, but in contemporary nomenclature this is an
#upper triangular Toeplitz matrix
#This constructs the dense matrix - Toeplitz matrices don't exist in Julia yet
function MatrixForm{T}(P::FormalPowerSeries{T}, m :: Integer)
    m<0 ? error(sprintf("Invalid matrix dimension %d requested", m)) : true
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
function reciprocal{T}(P::FormalPowerSeries{T}, n :: Integer)
    n<0 ? error(sprintf("Invalid inverse truncation length %d requested", n)) : true

    a = zeros(n) #Extract original coefficients in vector
    for (k,v) in P.c
        if 1<=k+1<=n
            a[k+1] = v
        end
    end
    a[1]==0 ? (error("Inverse does not exist")) : true
    inv_a1 = T<:Number ? 1.0/a[1] : inv(a[1])
    
    b = zeros(n) #Inverse
    b[1] = inv_a1
    for k=2:n
        b[k] = -inv_a1 * sum([b[k+1-i] * a[i] for i=2:k])
    end
    FormalPowerSeries{T}(b)
end


#Derivative of fps [H. Sec.1.4, p.18] 
function derivative{T}(P::FormalPowerSeries{T})
    c = Dict{BigInt, T}()
    for (k,v) in P.c
        if k != 0 && v != 0 
            c[k-1] = k*v
        end
    end
    FormalPowerSeries{T}(c)
end

#[H, Sec.1.4, p.18]
function isconstant{T}(P::FormalPowerSeries{T})
    for (k,v) in P.c
        if k!=0 && v!=0
            return false
        end
    end
    return true
end



# Power of fps [H, Sec.1.5, p.35]
function ^{T}(P::FormalPowerSeries{T}, n :: Integer)
    if n == 0
        return FormalPowerSeries{T}([convert(T, 1)])
    elseif n > 1
        #The straightforward recursive way
        return P^(n-1) * P
        #TODO implement the non-recursive formula 
    elseif n==1
        return P
    else
        error(@sprintf("I don't know what to do for power n = %d", n))
    end
end



# Composition of fps [H, Sec.1.5, p.35]
# This is the quick and dirty version
function compose{T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T})
    sum([v * Q^k for (k, v) in P.c]) 
end


#[H, p.45]
function isalmostunit{T}(P::FormalPowerSeries{T})
    (has(P.c, 0) && P.c[0]!=0) ? (return false) : true
    (has(P.c, 1) && P.c[1]!=0) ? (return true) : (return false)
end

###############################
# Formal Laurent series (fLs) #
###############################

type FormalLaurentSeries{Coefficients}
    c :: Dict{BigInt, Coefficients}
    FormalLaurentSeries{Coefficients}(c::Dict{BigInt, Coefficients}) = new(c)
    function FormalLaurentSeries{Coefficients}(v::Vector{Coefficients})
       c = new(c)
       for i=1:length(v)
           c[i] = v[i]
       end
       c
    end
end
    

###########
# Sandbox #
###########

MaxSeriesSize=10
MaxRange = 50
MatrixSize=15
TT=Int64

(n1, n2, n3) = int(rand(3)*MaxSeriesSize)

X = FormalPowerSeries{TT}(int(rand(n1)*MaxRange))
Y = FormalPowerSeries{TT}(int(rand(n2)*MaxRange))
Z = FormalPowerSeries{TT}(int(rand(n3)*MaxRange))

c = int(rand()*MatrixSize) #Size of matrix representation to generate

nzeros = int(rand()*MaxSeriesSize)
@assert X == trim(X)
XX = deepcopy(X)
for i=1:nzeros  
    idx = int(rand()*MaxRange)
    if !has(XX.c, idx)
        XX.c[idx] = convert(TT, 0)
    end
end
@assert trim(XX) == X

#Test addition, p.15, (1.3-4)
@assert X+X == 2X
@assert X+Y == Y+X
@assert MatrixForm(X+Y,c) == MatrixForm(X,c)+MatrixForm(Y,c)

#Test subtraction, p.15, (1.3-4)
@assert X-X == 0X
@assert X-Y == -(Y-X)
@assert MatrixForm(X-Y,c) == MatrixForm(X,c)-MatrixForm(Y,c)

#Test multiplication, p.15, (1.3-5)
@assert X*Y == Y*X
@assert MatrixForm(X*X,c) == MatrixForm(X,c)*MatrixForm(X,c)
@assert MatrixForm(X*Y,c) == MatrixForm(X,c)*MatrixForm(Y,c)
@assert MatrixForm(Y*X,c) == MatrixForm(Y,c)*MatrixForm(X,c)
@assert MatrixForm(Y*Y,c) == MatrixForm(Y,c)*MatrixForm(Y,c)

@assert X.*Y == Y.*X

#The reciprocal series has associated matrix that is the matrix inverse
#of the original series
#Force reciprocal to exist
X.c[0] = 1
discrepancy = (norm(inv(float(MatrixForm(X,c)))[1, :]'[:, 1] - tovector(reciprocal(X, c),c)))
if discrepancy > c*sqrt(eps())
    error(sprintf("Error %f exceeds tolerance %f", discrepancy, c*sqrt(eps())))
end
#@assert norm(inv(float(MatrixForm(X,c)))[1, :]'[:, 1] - tovector(reciprocal(X, c),c)) < c*sqrt(eps())

#Test differentiation
XX = derivative(X)
for (k, v) in XX.c
    k==0 ? continue : true
    @assert X.c[k+1] == v/(k+1)
end

#Test product rule [H, Sec.1.4, p.19]
@assert derivative(X*Y) == derivative(X)*Y + X*derivative(Y)

#Test right distributive law of composition [H, Sec.1.6, p.38]
@assert compose(X,Z)*compose(Y,Z) == compose(X*Y, Z)

#Test chain rule [H, Sec.1.6, p.40]
@assert derivative(compose(X,Y)) == compose(derivative(X),Y)*derivative(Y)
