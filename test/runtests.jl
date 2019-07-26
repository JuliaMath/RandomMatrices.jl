using RandomMatrices
using Random: seed!
using Test

@testset "RandomMatrices" begin
	seed!(1)
	include("GaussianEnsembles.jl")
	include("FormalPowerSeries.jl")
	include("Haar.jl")
	include("StochasticProcess.jl")
	
	@testset "densities" begin
		include("densities/Semicircle.jl")
		include("densities/TracyWidom.jl")
	end
end

