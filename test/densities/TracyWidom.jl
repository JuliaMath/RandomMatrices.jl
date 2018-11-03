using RandomMatrices
using Test

@testset "TracyWidom" begin

# CDF lies in correct range
@test all(i->(0<i<1), cdf(TracyWidom, randn(5)))

#Test far outside support
@test cdf(TracyWidom, -10.1) ≈ 0 atol=1e-14
@test cdf(TracyWidom, 10.1) ≈ 1 atol=1e-14

# Test exact values
# See https://arxiv.org/abs/0904.1581
@test cdf(TracyWidom,0,beta=1) ≈ 0.83190806620295 atol=1e-14
@test cdf(TracyWidom,0,beta=2) ≈ 0.96937282835526 atol=1e-14

@test isfinite(rand(TracyWidom, 10))
@test isfinite(rand(TracyWidom, 100))

end # testset
