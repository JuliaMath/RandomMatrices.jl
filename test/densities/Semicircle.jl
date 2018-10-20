using RandomMatrices
using Test

@testset "Semicircle" begin

let
    d = Semicircle()
    @test cdf(d, 2.0) == 1.0
    @test cdf(d, -2.0) == 0.0
    @test pdf(d, 2.0) == 0.0
    @test pdf(d, -2.0) == 0.0

    @test entropy(d) == log(2π) - 0.5
    @test mean(d) == median(d) == skewness(d) == 0
    @test modes(d) == [0]
    @test std(d) == var(d) == 1
    @test kurtosis(d) == 2

    @test isfinite(rand(d))
    @test moment(d, 4) == 2
    @test cumulant(d, 4) == 1
    @test moment(d, 5) == cumulant(d, 5) == freecumulant(d, 5) == 0

    @test cf(d, 1.0) == 0.0

    #XXX Throws AmosException - ???
    #@show mgf(d, 1.0)

    x = rand(d)
    @test insupport(d, x)
end

end # testset

