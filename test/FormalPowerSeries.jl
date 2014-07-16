#!/usr/bin/env julia
#
# Excursions in formal power series (fps)
#
# Jiahao Chen <jiahao@mit.edu> 2013

using Base.Test
using RandomMatrices

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
discrepancy = (norm(inv(float(MatrixForm(X,c)))[1, :]'[:, 1] - tovector(reciprocal(X, c),[0:c-1])))
if discrepancy > c*sqrt(eps())
    error(@sprintf("Error %f exceeds tolerance %f", discrepancy, c*sqrt(eps())))
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
