using RandomMatrices
using Compat.Test

#Test far outside support
@test cdf(TracyWidom, -10) ≈ 0 atol=1e-14
@test cdf(TracyWidom, 10) ≈ 1 atol=1e-14

# Test exact values
# See https://arxiv.org/abs/0904.1581
@test TWcdf(0,beta=1) ≈ 0.83190806620295 atol=1e-14
@test TWcdf(0,beta=2) ≈ 0.96937282835526 atol=1e-14

@test 0 < cdf(TracyWidom, rand()) < 1

@test isfinite(rand(TracyWidom, 10))
@test isfinite(rand(TracyWidom, 100))
