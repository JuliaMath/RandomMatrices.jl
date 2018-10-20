export Semicircle

"""
Wigner semicircle distribution

By default, creates a standardized distribution
with mean 0 and variance 1 (radius 2)
"""
struct Semicircle{T<:Real} <: ContinuousUnivariateDistribution
    mean::T
    radius::T
    Semicircle{T}(μ::T,r::T) where T = new(μ,r)
end
Semicircle(μ::T=0.0, r::T=2.0) where {T<:Real} = r > 0 ? Semicircle{T}(μ, r) :
    throw(ArgumentError("radius r must be positive, got $r"))
    
Semicircle{T}(d::Semicircle) where T = Semicircle{T}(convert(T, d.mean), convert(T, d.radius))
convert(::Type{Semicircle{T}}, d::Semicircle{T}) where T = d
convert(::Type{Semicircle}, d::Semicircle) = d
convert(::Type{Semicircle{T}}, d::Semicircle) where T = Semicircle{T}(convert(T, d.mean), convert(T, d.radius))

# Distribution function methods
###############################

# cumulative distribution function
function cdf(d::Semicircle{T}, x::T) where {T<:Real}
    a, r = d.mean, d.radius
    if insupport(d, x)
        return 0.5 + (x-a)/(π*r^2) * √(r^2 - (x-a)^2) + asin((x-a)/r)/π
    elseif x ≥ a
        return one(T)
    else
	return zero(T)
    end
end

# probability density function
function pdf(d::Semicircle{T}, x::T) where {T<:Real}
    a, r = d.mean, d.radius
    if insupport(d, x)
        return 2/(π*r^2) * √(r^2 - (x-a)^2)
    else
	return zero(T)
    end
end

# predicate is x in the support of the distribution?
insupport(d::Semicircle{T}, x::T) where {T<:Real} = abs(x-d.mean) < d.radius

function cdf(X::Semicircle{T}, x::V) where {T<:Real,V<:Real}
	TV = promote_type(T,V)
	cdf(convert(Semicircle{TV},X), convert(TV,x))
end

# probability density function
function pdf(X::Semicircle{T}, x::V) where {T<:Real,V<:Real}
	TV = promote_type(T,V)
	pdf(convert(Semicircle{TV},X), convert(TV,x))
end

# predicate is x in the support of the distribution?
function insupport(X::Semicircle{T}, x::V) where {T<:Real,V<:Real}
	TV = promote_type(T,V)
	insupport(convert(Semicircle{TV},X), convert(TV,x))
end

#Entropy methods
################

# entropy of distribution in nats
entropy(X::Semicircle)=log(π*X.radius) - 0.5

#Measures of central tendency methods
#####################################

# mean of distribution
mean(X::Semicircle)=X.mean

# median of distribution
median(X::Semicircle)=X.mean

# mode(s) of distribution as vector
modes(X::Semicircle{T}) where {T} = T[X.mean]

# kurtosis of the distribution
kurtosis(X::Semicircle{T}) where {T} = T(2)

# skewness of the distribution
skewness(X::Semicircle{T}) where {T} = T(0)

# standard deviation of distribution
std(X::Semicircle)=X.radius/2

# variance of distribution
var(X::Semicircle)=std(X)^2

# moment of distribution
function moment(X::Semicircle{T}, order::Integer) where {T<:Real}
    a, r = X.mean, X.radius
    if X.mean != 0
        return a^n*hypergeom([(1-n)/2, -n/2], 2, (r/a)^2)
    elseif iseven(order)
        return (0.5*r)^(2order) * T(catalannum(order÷2))
    else
	return zero(T)
    end
end

# cumulant of distribution
function cumulant(d::Semicircle{T}, order::Integer) where {T<:Real}
    if d.mean != 0
        throw(ArgumentError("Non-central Semicircle not supported"))
    end

    if order==0
	return one(T)
    elseif iseven(order)
        return (0.5*d.radius)^(2order) * T(lassallenum(order÷2))
    else
        return zero(T)
    end
end

# free cumulant of distribution
function freecumulant(d::Semicircle{T}, order::Int) where {T<:Real}
    if order == 0
        return one(T)
    elseif order == 1
        return mean(X)
    elseif order == 2
        return var(X)
    else
        return zero(T)
    end
end

#Generating function methods
############################

# characteristic function
function cf(d::Semicircle, t::Real)
  r = t / d.mean
  2 * besselj(1, r)/r
end

# moment generating function
function mgf(d::Semicircle, t::Real)
  r = t * d.mean
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
