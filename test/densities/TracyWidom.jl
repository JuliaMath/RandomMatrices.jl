using RandomMatrices
using Base.Test

#Test far outside support
#Tracy-Widom has support on all x>0, but the integration won't pick it up
@test RandomMatrices.pdf(TracyWidom, -10) == RandomMatrices.pdf(TracyWidom, 10) == 0

if isdefined(:ODE) && isa(ODE, Module)
    t = rand()
    @test RandomMatrices.pdf(TracyWidom, t) > 0
end

@test rand(TracyWidom, 10) > 0
@test rand(TracyWidom, 100) > 0
