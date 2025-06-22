using RandomMatrices
using Test

@testset "throws" begin

@test_throws ArgumentError rand(GaussianLaguerre(10, rand(Float64)), rand(2:5))
@test_throws ErrorException bidrand(GaussianLaguerre(10, rand(Float64)), rand(2:5))
@test_throws ArgumentError rand(Ginibre(10), rand(2:5))
@test_throws ArgumentError rand(Haar(10), rand(2:5))
@test_throws ArgumentError randfast(Haar(10), rand(2:5))

end # testset
