#A template for defining a new distribution
importall Distributions
# Note: these don't exist in Distributions.jl yet
export ThisDistribution, moment, cumulant, freecumulant, cf, cgf, mgf

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

#Standardized distribution - mean 0, variance 1

#Methods for the distribution
#The methods commented out below have default fallback methods

# Distribution function methods
###############################

# complementary cdf i.e. 1 - cdf
#ccdf(X::ThisDistribution, x::Real)=error(string("Not defined for ", typeof(X), " distribution"))

# cumulative distribution function
cdf(X::ThisDistribution, x::Real)=error(string("Not defined for ", typeof(X), " distribution"))

#complementary quantile based on log probability
#invlogccdf(X::ThisDistribution, x::Real)=error(string("Not defined for ", typeof(X), " distribution"))

# quantile based on log probability
#invlogcdf(X::ThisDistribution, x::Real)=error(string("Not defined for ", typeof(X), " distribution"))     

# ccdf returning log-probability
#logccdf(X::ThisDistribution, x::Real)=error(string("Not defined for ", typeof(X), " distribution"))

# cdf returning log-probability
#logcdf(X::ThisDistribution, x::Real)=error(string("Not defined for ", typeof(X), " distribution"))

# log probability density
#logpdf(X::ThisDistribution, x::Real)=error(string("Not defined for ", typeof(X), " distribution"))

# probability density function
pdf(X::ThisDistribution, x::Real)=error(string("Not defined for ", typeof(X), " distribution"))

# inverse of cdf (defined for p in (01))
quantile(X::ThisDistribution, x::Real)=error(string("Not defined for ", typeof(X), " distribution"))

# complementary quantile (i.e. using prob in right hand tail)
#cquantile(X::ThisDistribution, x::Real)=error(string("Not defined for ", typeof(X), " distribution"))

# predicate is x in the support of the distribution?
insupport(X::ThisDistribution, x::Real)=error(string("Not defined for ", typeof(X), " distribution"))     


#Entropy methods
################

# entropy of distribution in bits
#binaryentropy(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# entropy of distribution in nats
entropy(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))       

#Measures of central tendency methods
#####################################

# mean of distribution
mean(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# median of distribution
median(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# mode(s) of distribution as vector
modes(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# standard deviation of distribution
#std(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# variance of distribution
var(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# skewness of the distribution
skewness(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# kurtosis of the distribution
kurtosis(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# excess kurtosis of the distribution
excess_kurtosis(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# moment of distribution
moment(X::ThisDistribution, order::Integer)=error(string("Not defined for ", typeof(X), " distribution"))

# cumulant of distribution
cumulant(X::ThisDistribution, order::Integer)=error(string("Not defined for ", typeof(X), " distribution"))

# free cumulant of distribution
freecumulant(X::ThisDistribution, order::Integer)=error(string("Not defined for ", typeof(X), " distribution"))


#Generating function methods
############################

# characteristic function 
cf(X::ThisDistribution, t::Real)=error(string("Not defined for ", typeof(X), " distribution"))

# cumulant generating function
cgf(X::ThisDistribution, x::Real)=error(string("Not defined for ", typeof(X), " distribution"))

# moment generating function
mgf(X::ThisDistribution, x::Real)=error(string("Not defined for ", typeof(X), " distribution"))

#Sampling methods
#################

# random sampler
rand(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))


#Regression methods for generalized linear models (GLM)
#######################################################

# canonical link function for a distribution
canonicallink(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# deviance of fitted and observed responses
#deviance(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# vector of squared deviance residuals
#devresid(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# fit a distribution to data
fit(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# link function mapping mu to eta the linear predictor
linkfun(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# inverse link mapping eta to mu
linkinv(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# derivative of inverse link function
mueta(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# starting values of mean vector in GLMs
#mustart(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# validity check on linear predictor
valideta(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

# validity check on mean vector
validmu(X::ThisDistribution)=error(string("Not defined for ", typeof(X), " distribution"))

