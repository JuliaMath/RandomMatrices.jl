using RandomMatrices
using LinearAlgebra: norm
using Test

@testset "FormalPowerSeries" begin
seed!(4)

# Excursions in formal power series (fps)
MaxSeriesSize=8
MaxRange = 20
MatrixSize=10
TT=Int64

(n1, n2, n3) = rand(1:MaxSeriesSize, 3)

X = FormalPowerSeries{TT}(rand(1:MaxRange, n1))
Y = FormalPowerSeries{TT}(rand(1:MaxRange, n2))
Z = FormalPowerSeries{TT}(rand(1:MaxRange, n3))

c = rand(1:MatrixSize) #Size of matrix representation to generate

nzeros = rand(1:MaxSeriesSize)
@test X == trim(X)
XX = deepcopy(X)
for i=1:nzeros
    idx = rand(1:MaxRange)
    if !haskey(XX.c, idx)
        XX.c[idx] = convert(TT, 0)
    end
end
@test trim(XX) == X

#Test addition, p.15, (1.3-4)
@test X+X == 2X
@test X+Y == Y+X
@test MatrixForm(X+Y,c) == MatrixForm(X,c)+MatrixForm(Y,c)

#Test subtraction, p.15, (1.3-4)
@test X-X == 0X
@test X-Y == -(Y-X)
@test MatrixForm(X-Y,c) == MatrixForm(X,c)-MatrixForm(Y,c)

#Test multiplication, p.15, (1.3-5)
@test X*Y == Y*X
@test MatrixForm(X*X,c) == MatrixForm(X,c)*MatrixForm(X,c)
@test MatrixForm(X*Y,c) == MatrixForm(X,c)*MatrixForm(Y,c)
@test MatrixForm(Y*X,c) == MatrixForm(Y,c)*MatrixForm(X,c)
@test MatrixForm(Y*Y,c) == MatrixForm(Y,c)*MatrixForm(Y,c)

@test X.*Y == Y.*X

#The reciprocal series has associated matrix that is the matrix inverse
#of the original series
#Force reciprocal to exist
X.c[0] = 1
discrepancy = (norm(inv(float(MatrixForm(X,c)))[1, :][:, 1] - tovector(reciprocal(X, c), 0:c-1)))
tol = c*√eps()
if discrepancy > tol
    error(string("Error ", discrepancy, " exceeds tolerance ", tol))
end
#@test norm(inv(float(MatrixForm(X,c)))[1, :]'[:, 1] - tovector(reciprocal(X, c),c)) < c*sqrt(eps())

#Test differentiation
XX = derivative(X)
for (k, v) in XX.c
    k==0 && continue
    @test X.c[k+1] == v/(k+1)
end

#Test product rule [H, Sec.1.4, p.19]
@test derivative(X*Y) == derivative(X)*Y + X*derivative(Y)

#Test right distributive law of composition [H, Sec.1.6, p.38]
@test compose(X,Z)*compose(Y,Z) == compose(X*Y, Z)

#Test chain rule [H, Sec.1.6, p.40]
@test derivative(compose(X,Y)) == compose(derivative(X),Y)*derivative(Y)

end # testset
