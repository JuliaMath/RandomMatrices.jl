using RandomMatrices
using LinearAlgebra: I, tr, Diagonal, QRPackedQ
using Test

@testset "Haar" begin

N=5
A=randn(N,N)
B=randn(N,N)
Q=rand(Haar(1), N)

#Test case where there are no symbolic Haar matrices
@test_broken eval(expectation(:(A*B))) ≈ A*B
#Test case where there is one pair of symbolic Haar matrices
@test_broken tr(eval(expectation(:(A*Q*B*Q')))) ≈ tr(A)*tr(B)/N

println("Case 3")
@test_broken println("E(A*Q*B*Q'*A*Q*B*Q') = ", eval(expectation(N, :Q, :(A*Q*B*Q'*A*Q*B*Q'))))

for T in (Float64, ComplexF64)
    A = Stewart(T, N)
    @test A'A ≈ Matrix{T}(I, N, N)
    A2 = Matrix(A)
    @test A2 ≈ QRPackedQ(A.factors, A.τ)*Diagonal(A.signs)
    C = randn(T, N, N)
    @test A*C ≈ A2*C
    @test C*A ≈ C*A2
    @test A'*C ≈ A2'*C
    @test C*A' ≈ C*A2'
end

end # testset
