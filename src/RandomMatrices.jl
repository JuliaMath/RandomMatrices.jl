require("Distributions")
require("Catalan")
require("GSL")

module RandomMatrices
    importall Distributions
    using Catalan
    
    #Classical Gaussian matrix ensembles
    include("GaussianEnsembleSamples.jl")
    include("GaussianEnsembleDistributions.jl")
    
    # Classical univariate distributions
    ####################################
    include("densities/Semicircle.jl")
    include("densities/TracyWidom.jl")
    
    #Fast histogrammer for matrix eigenvalues
    include("FastHistogram.jl")
    
    #Formal power series
    include("FormalPowerSeries.jl")
    
    #Symbolic integrator of uniform Haar matrices
    include("HaarSymbolic.jl")
    include("HaarMeasure.jl")
end #module
