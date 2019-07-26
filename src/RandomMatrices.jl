module RandomMatrices

using Combinatorics
using GSL
using SpecialFunctions, FastGaussQuadrature
using LinearAlgebra

import Base: isinf, rand, convert
import Distributions: ContinuousUnivariateDistribution,
    ContinuousMatrixDistribution,
    Beta, Chi,
    cdf, pdf, entropy, insupport, mean, median, modes, kurtosis, skewness, std, var, moment

export bidrand,    #Generate random bidiagonal matrix
       tridrand,   #Generate random tridiagonal matrix
       sprand,     #Generate random sparse matrix
       eigvalrand, #Generate random set of eigenvalues
       eigvaljpdf, #Eigenvalue joint probability density
       cumulant, freecumulant, cf, mgf,
       cdf, pdf, entropy, insupport, mean, median, modes, kurtosis, skewness, std, var, moment

const Dim2 = Tuple{Int, Int} #Dimensions of a rectangular matrix

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

end
