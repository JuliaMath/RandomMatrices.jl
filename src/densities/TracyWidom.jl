export TracyWidom

"""
Tracy-Widom distribution

The probability distribution of the normalized largest eigenvalue of a random
matrix with iid Gaussian matrix elements of variance 1/2 and mean 0.

The cdf of Tracy-Widom is given by

``
F_beta (s) = lim_{n→∞} Pr(√2 n^{1/6} (λ_max - √(2n) ≤ s)
``

where beta = 1, 2, or 4 for the orthogonal, unitary, or symplectic ensembles.

References:

1. doi:10.1016/0370-2693(93)91114-3
2. doi:10.1007/BF02100489
3. doi.org/10.1007/BF02099545

Numerical routines adapted from Alan Edelman's course notes for MIT 18.338,
Random Matrix Theory, 2016.
"""
struct TracyWidom <: ContinuousUnivariateDistribution end


"""
Cumulative density function of the Tracy-Widom distribution.

Computes the Tracy-Widom distribution by Bornemann's method of evaluating
a finite dimensional approximation to the Fredholm determinant using quadrature.

doi.org/10.1090/S0025-5718-09-02280-7

# Arguments
* `d::TracyWidom` or `Type{TracyWidom}`: an instance of `TracyWidom` or the type itself
* `s::Real`: The point at which to evaluate the cdf
* `beta::Integer = 2`: The Dyson index defining the distribution. Takes values 1, 2, or 4
* `num_points::Integer = 25`: The number of points in the quadrature
"""
function cdf(d::TracyWidom, s::T; beta::Integer=2, num_points::Integer=25) where {T<:Real}
    beta ∈ (1,2,4) || throw(ArgumentError("Beta must be 1, 2, or 4"))
    quad = gausslegendre(num_points)
    _TWcdf(s, beta, quad)
end

function cdf(d::Type{TracyWidom}, s::T; beta::Integer=2, num_points::Integer=25) where {T<: Real}
    cdf(d(), s, beta=beta, num_points=num_points)
end

function cdf(d::TracyWidom, s_vals::AbstractArray{T}; beta::Integer=2, num_points::Integer=25) where {T<:Real}
    beta ∈ (1,2,4) || throw(ArgumentError("Beta must be 1, 2, or 4"))
    quad = gausslegendre(num_points)
    [_TWcdf(s, beta, quad) for s in s_vals]
end

function cdf(d::Type{TracyWidom}, s_vals::AbstractArray{T}; beta::Integer=2, num_points::Integer=25) where {T<:Real}
    cdf(d(), s_vals, beta=beta, num_points=num_points)
end

function _TWcdf(s::T, beta::Integer, quad::Tuple{Array{Float64,1},Array{Float64,1}}) where {T<:Real}
    if beta == 2
        kernel = ((ξ,η) -> _K2tilde(ξ,η,s))
        return _fredholm_det(kernel, quad)
    elseif beta == 1
        kernel = ((ξ,η) -> _K1tilde(ξ,η,s))
        return _fredholm_det(kernel, quad)
    elseif beta == 4
        kernel2 = ((ξ,η) -> _K2tilde(ξ,η,s*sqrt(2)))
        kernel1 = ((ξ,η) -> _K1tilde(ξ,η,s*sqrt(2)))
        F2 = _fredholm_det(kernel2, quad)
        F1 = _fredholm_det(kernel1, quad)
        return (F1 + F2/F1) / 2
    end
end

function _fredholm_det(kernel::Function, quad::Tuple{Array{T,1},Array{T,1}}) where {T<:Real}
    nodes, weights = quad
    N = length(nodes)
    sqrt_weights = sqrt.(weights)
    weights_matrix = kron(transpose(sqrt_weights),sqrt_weights)
    K_matrix = [kernel(ξ,η) for ξ in nodes, η in nodes]
    det(I - weights_matrix .* K_matrix)
end

_ϕ(ξ, s) =  s + 10*tan(π*(ξ+1)/4)
_ϕprime(ξ) = (5π/2)*(sec(π*(ξ+1)/4))^2

# For beta = 2
function _airy_kernel(x, y)
    if x==y
        return (airyaiprime(x))^2 - x * (airyai(x))^2
    else
        return (airyai(x) * airyaiprime(y) - airyai(y) * airyaiprime(x)) / (x - y)
    end
end

_K2tilde(ξ,η,s) = sqrt(_ϕprime(ξ) * _ϕprime(η)) * _airy_kernel(_ϕ(ξ,s), _ϕ(η,s))

# For beta = 1
_A_kernel(x,y) = airyai((x+y)/2) / 2
_K1tilde(ξ,η,s) = sqrt(_ϕprime(ξ) * _ϕprime(η)) * _A_kernel(_ϕ(ξ,s), _ϕ(η,s))

"""
Samples the largest eigenvalue of the n × n GUE matrix
"""
function rand(d::TracyWidom, n::Int)
    n > 1 || error("Cannot construct $n × $n matrix")
    if n < 100
        k = n
    else #Exploit the fact that the eigenvector is concentrated in the top left corner
	k = round(Int, n-10*n^(1/3)-1)
    end
    a=randn(n-k+1)
    b=[χ(i) for i=(n-1):-1:k]
    v=eigmax(SymTridiagonal(a, b))
end
rand(d::Type{TracyWidom}, t::Integer) = rand(d(), t)
