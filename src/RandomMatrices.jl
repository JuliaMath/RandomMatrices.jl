require("Distributions")

module RandomMatrices
    using Distributions

    #Classical Gaussian matrix ensembles
    include("GaussianEnsembleSamples.jl")
    #Classical univariate distributions
    include("GaussianDensities.jl")
    #Tracy-Widom distribution
    include("TracyWidom.jl")
    
    #Fast histogrammer for matrix eigenvalues
    include("FastHistogram.jl")
    
    #Formal power series
    include("FormalPowerSeries.jl")
end
