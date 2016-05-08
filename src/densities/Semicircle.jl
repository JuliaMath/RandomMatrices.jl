#Wigner semicircle distribution

import Distributions
import Distributions: ContinuousUnivariateDistribution

export Semicircle

immutable Semicircle <: ContinuousUnivariateDistribution
    mean::Float64
    radius::Float64
    Semicircle(mu, r) = r > 0 ? new(float64(mu), float64(r)) :
        error("radius r must be positive") #Constructor
end

#Standardized distribution - mean 0, variance 1
Semicircle(mu) = Semicircle(mu, 2.0)
Semicircle() = Semicircle(0.0, 2.0)

# Distribution function methods
###############################

# cumulative distribution function
function cdf(X::Semicircle, x::Real)
    r, a = X.mean, X.radius
    insupport(x) ? 0.5 + (x-a)/(pi*r^2) * sqrt(r^2 - (x-a)^2) + 1/pi * asin((x-a)/r) : (x>a ? 1.0 : 0.0)
end

# probability density function
function pdf(X::Semicircle, x::Real)
    r, a = X.mean, X.radius
    insupport(x) ? 2/(pi*r^2) * sqrt(r^2 - (x-a)^2) : 0.0
end

# predicate is x in the support of the distribution?
insupport(X::Semicircle, x::Real)=abs(x)<=X.radius


#Entropy methods
################

# entropy of distribution in nats
entropy(X::Semicircle)=log(pi*X.radius) - 0.5

#Measures of central tendency methods
#####################################

# mean of distribution
mean(X::Semicircle)=X.mean

# median of distribution
median(X::Semicircle)=X.mean

# mode(s) of distribution as vector
modes(X::Semicircle)=[X.mean]

# kurtosis of the distribution
kurtosis(X::Semicircle)=2

# skewness of the distribution
skewness(X::Semicircle)=0

# standard deviation of distribution
std(X::Semicircle)=X.radius/2

# variance of distribution
var(X::Semicircle)=std(X)^2

# moment of distribution
function moment(X::Semicircle, order::Integer)
    a, r = X.mean, X.radius
    if X.mean != 0
        a^n*hypergeom([(1-n)/2, -n/2], 2, (r/a)^2)
    else
        order%2 ? (0.5*r)^(2n) * catalan(div(order,2)) : 0
    end
end

# cumulant of distribution
function cumulant(X::Semicircle, order::Integer)
  if X.mean != 0 error("not supported") end
  if order%2
    order==0 ? 1 : (0.5*r)^(2n) * lassalle(order/2)
  else
    0
  end
end

# free cumulant of distribution
function freecumulant(X::Semicircle, order::Integer)
  if order == 0
    return 1
  elseif order == 1
    return mean(X)
  elseif order == 2
    return var(X)
  else
    return 0
  end
end


#Generating function methods
############################

# characteristic function
function cf(X::Semicircle, t::Real)
  r = t * X.mean
  2 * besselj(1, r)/r
end

# moment generating function
function mgf(X::Semicircle, t::Real)
  r = t * X.mean
  2 * besseli(1, r)/r
end

#Sampling methods
#################

# random sampler
# Use relationship with beta distribution
function rand(X::Semicircle)
  Y = rand(Beta(1.5, 1.5))
  X.mean + 2 * X.radius * Y - X.radius
end
