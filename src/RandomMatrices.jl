require("Distributions")
require("Catalan")
require("GSL")

module RandomMatrices
importall Distributions
using Catalan

#If the GNU Scientific Library is present, turn on additional functionality.
_HAVE_GSL = try
  using GSL
  _HAVE_GSL
end
if _HAVE_GSL
  using GSL
end

export bidrand,    #Generate random bidiagonal matrix
       tridrand,   #Generate random tridiagonal matrix
       sprand,     #Generate random sparse matrix
       eigvalrand, #Generate random set of eigenvalues
       eigvaljpdf  #Eigenvalue joint probability density

typealias Dim2 (Int, Int) #Dimensions of a rectangular matrix

#Classical Gaussian matrix ensembles
include("GaussianEnsembles.jl")

# Classical univariate distributions
include("densities/Semicircle.jl")
include("densities/TracyWidom.jl")

#Fast histogrammer for matrix eigenvalues
include("FastHistogram.jl")

#Formal power series
include("FormalPowerSeries.jl")
end

