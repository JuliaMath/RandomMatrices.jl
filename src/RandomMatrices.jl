module RandomMatrices

using Combinatorics
using Compat
using GSL
using ODE

import Base: isinf, rand
import Distributions: ContinuousMatrixDistribution, Chi, pdf

export bidrand,    #Generate random bidiagonal matrix
       tridrand,   #Generate random tridiagonal matrix
       sprand,     #Generate random sparse matrix
       eigvalrand, #Generate random set of eigenvalues
       eigvaljpdf, #Eigenvalue joint probability density
       pdf

typealias Dim2 Tuple{Int, Int} #Dimensions of a rectangular matrix

#Classical Gaussian matrix ensembles
include("GaussianEnsembles.jl")

# Classical univariate distributions
include("densities/Semicircle.jl")
include("densities/TracyWidom.jl")

# Ginibre
include("Ginibre.jl")

# determinantal point processes
include("dpp.jl")

#Generating matrices of Haar measure
include("Haar.jl")
include("HaarMeasure.jl")

#Fast histogrammer for matrix eigenvalues - hist_eig
include("FastHistogram.jl")

#Formal power series
include("FormalPowerSeries.jl")

#Statistical tests based on random matrix theory
include("StatisticalTests.jl")

#Stochastic processes
include("StochasticProcess.jl")

#Invariant ensembles
if Pkg.installed("ApproxFun")!== nothing
    include("InvariantEnsembles.jl")
end

end
