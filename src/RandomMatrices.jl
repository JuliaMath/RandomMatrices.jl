require("Distributions")

module RandomMatrices
    using Distributions

    #Classical Gaussian matrix ensembles
    include("GaussianEnsembleSamples.jl")
    # Classical univariate distributions
    ####################################
    include("densities/Semicircle.jl")
    #Tracy-Widom distribution
    include("densities/TracyWidom.jl")
    
    #Fast histogrammer for matrix eigenvalues
    include("FastHistogram.jl")
    
    #Formal power series
    include("FormalPowerSeries.jl")
end
