using RandomMatrices
using Test

@testset "GaussianEnsembles" begin

@test Wigner{3} == GaussianHermite{3}
@test Wishart{3} == GaussianLaguerre{3}

n = 25

for (β, T, N) in [(1, Real, n), (2, Complex, n), (4, Complex, 2n)]
    @testset "Wigner (β = $(β))" begin
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
    @testset "Wishart (β = $(β))" begin
        a = 2(rand(1:5) + β * n)
        d = Wishart(β, a)
        A = rand(d, n)
        @test eltype(A) <: T
        @test size(A) == (N, N)
        
        @test_throws UndefVarError tridrand(d, n) # = At
        #@test eltype(At) <: Real
        #@test size(At) == (n, n)

        @test_throws UndefVarError eigvalrand(d, n) # = vals
        #@test eltype(vals) <: Real
        #@test length(vals) == n

        #vd = RandomMatrices.VandermondeDeterminant(vals, β)
        #@test isa(vd, Real)

        ed = eigvaljpdf(d, vals)
        @test isa(ed, Real)
    end
end

end # testset
