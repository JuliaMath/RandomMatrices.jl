using RandomMatrices
using Test

n = 25

for (β, T, N) in [(1, Real, n), (2, Complex, n), (4, Complex, 2n)]
    @testset "Ginibre (β = $(β))" begin
        d = Ginibre(β)
        A = rand(d, n)
        @test eltype(A) <: T
        @test size(A) == (N, N)
    end
end