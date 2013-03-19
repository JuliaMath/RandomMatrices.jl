require("Distributions")

module RandomMatrices
    using Distributions

    #Classical Gaussian matrix ensembles
    include("GaussianEnsembleSamples.jl")
    #Classical univariate distributions
    include("GaussianDensities.jl")
    
    #Fast histogrammer for matrix eigenvalues
    include("FastHistogram.jl")
end
