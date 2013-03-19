using Distributions

immutable Semicircle <: ContinuousUnivariateDistribution
    mean::Float64
    radius::Float64
    Semicircle(mu, r) = r > 0 ? new(float64(mu), float64(r)) : error("radius r must be positive")
end

#Standard (normalized) - note this is NOT the standard used in to RMT (r = 2)
Semicircle(mu) = Semicircle(mu, 1.0)
Semicircle() = Semicircle(0.0, 1.0)

# Properties
# For convenience, this is ordered like the Wikipedia infobox template ;)

#Probability density function
function pdf(d::Semicircle, x::Real)
    r, a = d.mean, d.radius
    return 2/(pi*r^2) * sqrt(r^2 - (x-a)^2)
end

#Cumulative density function
function cdf(d::Semicircle, x::Real)
    r, a = d.mean, d.radius
    if x < -a
        return 0.0
    elseif x > a
        return 1.0
    else
        return 0.5 + (x-a)/(pi*r^2) * sqrt(r^2 - (x-a)^2) + 1/pi * asin((x-a)/r)
    end
end

function logpdf(d::Laplace, x::Real)
    log(pdf(d, x))
end
function quantile(d::Laplace, p::Real)
    0.0
end

mean(d::Semicircle) = d.mean
median(d::Semicircle) = d.mean
mode(d::Semicircle) = d.mean

std(d::Semicircle) = (d.radius/2)
var(d::Semicircle) = (d.radius/2)^2
skewness(d::Semicircle) = 0.0
kurtosis(d::Semicircle) = 2.0
entropy(d::Semicircle) = log(pi*R) - 0.5

## Moment generating function
function mgf(d::Semicircle, t::Real)
    r = t * d.mean
    2 * besseli(1, r)/r
end

# Characteristic function
function cf(d::Semicircle, t::Real)
    r = t * d.mean
    2 * besselj(1, r)/r
end

## Common methods
#XXX Use the fact that Semicircle(0.5, 0.5) = Beta(1.5, 1.5)
rand(d::Semicircle) = 2*(rand(Beta(1.5, 1.5) - 0.5))
insupport(d::Normal, x::Number) = real_valued(x) & (abs(x) <= d.radius)
