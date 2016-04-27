using RandomMatrices
using Base.Test

srand(1)
dx = 0.001

p = WhiteNoiseProcess(dx)
@test start(p) == nothing
next(p)
@test done(p) == false

p = BrownianProcess(dx)
@test start(p) == 0
next(p, 0.0)
@test done(p) == false

p = AiryProcess(dx, 2.0)
S = start(p)
@test S == fill(-2/dx^2, 1, 1)
next(p, S)
@test done(p) == false

