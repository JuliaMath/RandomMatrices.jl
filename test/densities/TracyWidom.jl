using RandomMatrices
using Compat.Test

#Test far outside support
#Tracy-Widom has support on all x>0, but the integration won't pick it up
@test pdf(TracyWidom, -10) == pdf(TracyWidom, 10) == 0
@test cdf(TracyWidom, -10) == 0
@test cdf(TracyWidom, 10) == 1

if isdefined(:OrdinaryDiffEq) && isa(OrdinaryDiffEq, Module)
    t = rand()
    @test pdf(TracyWidom, t) > 0
    @test 0 < cdf(TracyWidom, t) < 1
end

@test isfinite(rand(TracyWidom, 10))
@test isfinite(rand(TracyWidom, 100))
