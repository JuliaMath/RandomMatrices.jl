using RandomMatrices
using Test

@testset "GaussianEnsembles" begin

@test Wigner{3} == GaussianHermite{3}

n = 25

for (β, T, N) in [(1, Real, n), (2, Complex, n), (4, Complex, 2n)]
    d = Wigner(β)
    A = rand(d, n)
    @test eltype(A) <: T
    @test size(A) == (N, N)
    
    At = tridrand(d, n)
    @test eltype(At) <: Real
    @test size(At) == (n, n)

    vals = eigvalrand(d, n)
    @test eltype(vals) <: Real
    @test length(vals) == n

    vd = RandomMatrices.VandermondeDeterminant(vals, β)
    @test isa(vd, Real)

    ed = eigvaljpdf(d, vals)
    @test isa(ed, Real)
end

end # testset
