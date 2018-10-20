using RandomMatrices
using Test

@testset "StochasticProcess" begin

seed!(1)
dx = 0.001

let
p = WhiteNoiseProcess(dx)
x, S = iterate(p)
@test typeof(x) == typeof(dx)
@test isa(S, Tuple{})
@test iterate(p, S) != nothing
end

let
p = BrownianProcess(dx)
x, S = iterate(p)
@test x == S
@test typeof(x) == typeof(dx)
@test iterate(p, S) != nothing
end

let
p = AiryProcess(dx, 2.0)
x, S = iterate(p)
@test S[1,1] == -2/dx^2
@test next!(p, S) != nothing
end

end # testset
