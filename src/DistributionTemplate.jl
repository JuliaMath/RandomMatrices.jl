#A template for defining a new distribution
importall Distributions
# Note: these don't exist in Distributions.jl yet
export moment, cumulant, freecumulant, cf, cgf, mgf

#Change type below to an actual type:
# DiscreteUnivariateDistribution    
# ContinuousUnivariateDistribution  
# DiscreteMultivariateDistribution  
# ContinuousMultivariateDistribution
# ContinuousMatrixDistribution      
# DiscreteMatrixDistribution        
type ThisDistribution <: Distribution
    data::Any #Internal fields
    #beta::Real 
    ThisDistribution(data) = new(data) #Constructor
end

#Methods for the distribution
#The methods commented out below have default fallback methods

#Distribution function methods
#ccdf(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))          # complementary cdf i.e. 1 - cdf
cdf(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))           # cumulative distribution function
#cquantile(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))     # complementary quantile (i.e. using prob in right hand tail)
#invlogccdf(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))    # complementary quantile based on log probability
#invlogcdf(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))     # quantile based on log probability
#logccdf(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))       # ccdf returning log-probability
#logcdf(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))        # cdf returning log-probability
#logpdf(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))        # log probability density
#logpmf(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))        # log probability mass
#pmf(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))           # probability mass function (DiscreteDistribution)
pdf(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))           # probability density function (ContinuousDistribution)
quantile(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))      # inverse of cdf (defined for p in (01))
insupport(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))     # predicate is x in the support of the distribution?

#Entropy methods
#binaryentropy(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution")) # entropy of distribution in bits
entropy(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))       # entropy of distribution in nats

#Measures of central tendency methods
mean(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))          # mean of distribution
median(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))        # median of distribution
modes(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))         # mode(s) of distribution as vector
kurtosis(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))      # kurtosis of the distribution
skewness(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))      # skewness of the distribution
#std(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))           # standard deviation of distribution
var(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))           # variance of distribution
moment(X::ThisDistribution, order::Integer)=error(string("Not defined for ", typeof(X), " distribution"))        # moment of distribution
cumulant(X::ThisDistribution, order::Integer)=error(string("Not defined for ", typeof(X), " distribution"))      # cumulant of distribution
freecumulant(X::ThisDistribution, order::Integer)=error(string("Not defined for ", typeof(X), " distribution"))  # free cumulant of distribution

#Generating function methods
cf(X::ThisDistribution, t::Real)=error(string("Not defined for ", typeof(X), " distribution"))   # characteristic function 
cgf(X::ThisDistribution, x::Real)=error(string("Not defined for ", typeof(X), " distribution"))  # cumulant generating function
mgf(X::ThisDistribution, x::Real)=error(string("Not defined for ", typeof(X), " distribution"))  # moment generating function

#Sampling methods
rand(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))          # random sampler

#Regression methods for generalized linear models (GLM)
canonicallink(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution")) # canonical link function for a distribution
#deviance(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))      # deviance of fitted and observed responses
#devresid(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))      # vector of squared deviance residuals
fit(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))           # fit a distribution to data
linkfun(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))       # link function mapping mu to eta the linear predictor
linkinv(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))       # inverse link mapping eta to mu
mueta(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))         # derivative of inverse link function
#mustart(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))       # starting values of mean vector in GLMs
valideta(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))      # validity check on linear predictor
validmu(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))       # validity check on mean vector

