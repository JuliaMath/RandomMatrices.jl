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

# Definition [H, p.9] - we limit this to finitely many coefficients
type FormalPowerSeries{Coefficients}
    #TODO this really should be a sparse vector
    c :: Vector{Coefficients}
    FormalPowerSeries{Coefficients <: Number}(c::Vector{Coefficients}) = new(c)
end

#Comparisons allowed between fps and Vector
function =={T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T})
    a, b = P.c, Q.c
    n1, n2 = length(a), length(b)
    #Pad shorter coefficient list with zeros
    (n1 > n2) ? (b = [b; zeros(T, n1 - n2)]) : (a = [a; zeros(T, n2 - n1)])
    return all([a[i] == b[i] for i=1:max(n1, n2)])
end

function =={T}(P::Vector{T}, Q::FormalPowerSeries{T})
    a, b = P, Q.c
    n1, n2 = length(a), length(b)
    #Pad shorter coefficient list with zeros
    (n1 > n2) ? (b = [b; zeros(T, n1 - n2)]) : (a = [a; zeros(T, n2 - n1)])
    return all([a[i] == b[i] for i=1:max(n1, n2)])
end

function =={T}(P::FormalPowerSeries{T}, Q::Vector{T})
    a, b = P.c, Q
    n1, n2 = length(a), length(b)
    #Pad shorter coefficient list with zeros
    (n1 > n2) ? (b = [b; zeros(T, n1 - n2)]) : (a = [a; zeros(T, n2 - n1)])
    return all([a[i] == b[i] for i=1:max(n1, n2)])
end

# Basic housekeeping and properties

# Remove extraneous zeros
function trim{T}(P::FormalPowerSeries{T})
    NumExtraZeros::Integer = 0
    for k=length(P.c):-1:2
        P.c[k]==convert(T, 0) ? (NumExtraZeros+=1) : break
    end
    return FormalPowerSeries{T}(copy(P.c[1:end-NumExtraZeros]))
end

length{T}(P::FormalPowerSeries{T})=length(P.c)

# Basic arithmetic [H, p.10]
function +{T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T})
    a, b = P.c, Q.c
    n1, n2 = length(a), length(b)
    #Pad shorter coefficient list with zeros
    (n1 > n2) ? (b = [b; zeros(T, n1 - n2)]) : (a = [a; zeros(T, n2 - n1)])
    FormalPowerSeries{T}(a + b)
end

function -{T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T})
    a, b = P.c, Q.c
    n1, n2 = length(a), length(b)
    #Pad shorter coefficient list with zeros
    (n1 > n2) ? (b = [b; zeros(T, n1 - n2)]) : (a = [a; zeros(T, n2 - n1)])
    FormalPowerSeries{T}(a - b)
end

#negation
-{T}(P::FormalPowerSeries{T}) = FormalPowerSeries{T}(-P.c)

#multiplication by scalar
*{T}(P::FormalPowerSeries{T}, k::Number) = FormalPowerSeries{T}(P.c * k)
*{T}(k::Number, P::FormalPowerSeries{T}) = FormalPowerSeries{T}(P.c * k)

#Cauchy product [H, p.10]
#
#TODO This is actually a discrete convolution so it should be amenable
#      to an FFT algorithm
*{T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T}) = CauchyProduct(P, Q)
function CauchyProduct{T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T})
    a, b = P.c, Q.c
    n1, n2 = length(a), length(b)
    c = zeros(T, n1+n2-1)
    for i1=1:n1
        for i2=1:n2
            c[i1+i2-1] += a[i1]*b[i2]
        end
    end
    FormalPowerSeries{T}(c)
end

#Hadamard product [H, p.10] - the elementwise product
.*{T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T}) = HadamardProduct(P, Q)
function HadamardProduct{T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T})
    a, b = P.c, Q.c
    n1, n2 = length(a), length(b)
    n1>n2 ? (a=a[1:n2]) : (b=b[1:n1])
    FormalPowerSeries{T}(a.*b)
end

#The identity element over the complex numbers
function eye{T <: Number}(P::FormalPowerSeries{T})
    FormalPowerSeries{T}([convert(T, 0), convert(T, 1)])
end

isunit{T <: Number}(P::FormalPowerSeries{T}) = P==eye(P)

# [H, p.12]
isnonunit{T}(P::FormalPowerSeries{T}) = P.c[1]==0 && !isunit(P)

#Constructs the top left m x m block of the (infinite) semicirculant matrix
#associated with the fps [H, Sec.1.3, p.14]
#[H] calls it the semicirculant, but in contemporary nomenclature this is an
#upper triangular Toeplitz matrix
#This constructs the dense matrix - Toeplitz matrices don't exist in Julia yet
function MatrixForm{T}(P::FormalPowerSeries{T}, m :: Integer)
    m<0 ? error(sprintf("Invalid matrix dimension %d requested", m)) : true
    
    a = P.c
    n = length(a)
    M = sum([diagm(repmat(ones(T,1,1)*a[i], 1, m+1-i), i-1) for i=1:min(m, n)])
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
function inv{T}(P::FormalPowerSeries{T}, Size :: Integer)
    Size<0 ? error(sprintf("Invalid inverse truncation length %d requested", Size)) : true
    
    a = Size>=length(P) ? [P.c; zeros(T, Size-length(P))] : P.c[1:Size]
    n = length(a)
    a[1] == 0 ? (error("Inverse does not exist")) : true
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
    a = P.c
    n = length(a)
    FormalPowerSeries{T}(n==1 ? [convert(T, 0)]: [(k-1)*a[k] for k=2:n])
end

function isconstant{T <: Number}(P::FormalPowerSeries{T})
    length(P.c)==1 || all([c==0 for c in P.c[2:end]])
end

#[H, Sec.1.4, p.18]
isconst{T}(P::FormalPowerSeries{T}) = length(P)==1 || all([x==0 for x in P[2:end]])


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
#TODO this doesn't work in infix notation yet
#(⋅){T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T}) = compose(P, Q)
function compose{T}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T})
    a, b = P.c, Q.c
    n, m = length(a), length(b)

    l = (m-1)*(n-1)+1
    bb = zeros(T, n, l)
    QQ = FormalPowerSeries{T}([convert(T, 1)])
    for i=1:n
        bb[i, 1:length(QQ.c)] = QQ.c
        QQ *= Q
    end
    for i=1:n
        bb[i, :] *= a[i]
    end
    
    return FormalPowerSeries{T}(sum(bb, 1)'[:,1])
end


#[H, p.45]
isalmostunit{T}(P::FormalPowerSeries{T}) = length(P) > 1 && P.c[1]==0 && P.c[2] != 0
###########
# Sandbox #
###########

MaxSeriesSize=50
MaxRange = 1000
MatrixSize=150

(n1, n2, n3) = int(rand(3)*MaxSeriesSize)

X = FormalPowerSeries{Int64}(int(rand(n1)*MaxRange))
Y = FormalPowerSeries{Int64}(int(rand(n2)*MaxRange))
Z = FormalPowerSeries{Int64}(int(rand(n3)*MaxRange))

c = int(rand()*MatrixSize) #Size of matrix representation to generate


#Test very basic things
@assert length(X) == length(X.c)
@assert length(Y) == length(Y.c)

TT = typeof(X.c[1])
nzeros = int(rand()*MaxSeriesSize)
@assert trim(FormalPowerSeries{TT}([X.c; zeros(TT, nzeros)])) == X
@assert trim(FormalPowerSeries{TT}([Y.c; zeros(TT, nzeros)])) == Y

#padded vectors for later testing
Xv = [X.c; zeros(TT, max(0, length(Y.c)-length(X.c)))]
Yv = [Y.c; zeros(TT, max(0, length(X.c)-length(Y.c)))]

#Test addition, p.15, (1.3-4)
@assert (X+X).c == 2X.c
@assert (X+Y).c == Xv + Yv
@assert (Y+X).c == Yv + Xv
@assert (Y+Y).c == 2Y.c
@assert MatrixForm(X+Y,c) == MatrixForm(X,c)+MatrixForm(Y,c)

#Test subtraction, p.15, (1.3-4)
@assert (X-X).c == 0X.c
@assert (X-Y).c == Xv - Yv
@assert (Y-X).c == Yv - Xv
@assert (Y-Y).c == 0Y.c
@assert MatrixForm(X-Y,c) == MatrixForm(X,c)-MatrixForm(Y,c)

#Test multiplication, p.15, (1.3-5)
@assert MatrixForm(X*X,c) == MatrixForm(X,c)*MatrixForm(X,c)
@assert MatrixForm(X*Y,c) == MatrixForm(X,c)*MatrixForm(Y,c)
@assert MatrixForm(Y*X,c) == MatrixForm(Y,c)*MatrixForm(X,c)
@assert MatrixForm(Y*Y,c) == MatrixForm(Y,c)*MatrixForm(Y,c)

@assert X.*X == Xv.*Xv
@assert X.*Y == Xv.*Yv
@assert Y.*X == Yv.*Xv
@assert Y.*Y == Yv.*Yv

#Test reversion of series
#The reciprocal series has associated matrix that is the matrix inverse
#of the original series
#@assert inv(float(MatrixForm(X,c)))[1, :]'[:, 1] == inv(X, c)

#Test differentiation
XX = X
for i =1:length(X)-1
    XX = derivative(XX)
end
@assert XX == [X.c[end]*factorial(length(X)-1)]

#Test product rule [H, Sec.1.4, p.19]
@assert derivative(X*Y) == derivative(X)*Y + X*derivative(Y)

#Test composition against the one-liner
compose_quickanddirty{T <: Number}(P::FormalPowerSeries{T}, Q::FormalPowerSeries{T}) = sum([P.c[i] * Q^(i-1) for i=1:length(P.c)]) 

#@assert ⋅(X,Y) == compose(X, Y)
@assert compose(X, X) == compose_quickanddirty(X, X)
@assert compose(X, Y) == compose_quickanddirty(X, Y)
@assert compose(Y, X) == compose_quickanddirty(Y, X)
@assert compose(Y, Y) == compose_quickanddirty(Y, Y)

#Test right distributive law of composition [H, Sec.1.6, p.38]
@assert compose(X,Z)*compose(Y,Z) == compose(X*Y, Z)

#Test chain rule [H, Sec.1.6, p.40]
@assert derivative(compose(X,Y)) == compose(derivative(X),Y)*derivative(Y)
