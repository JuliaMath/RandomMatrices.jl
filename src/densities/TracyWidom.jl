export TracyWidom

include("gradient.jl")

immutable TracyWidom <: ContinuousUnivariateDistribution end


"""
Computes the Tracy-Widom distribution by directly solving the
Painlevé II equation by numerical integration
"""
function pdf{S<:Real}(d::TracyWidom, t::S, t0::S = convert(S, -8.0))
    t≤t0 && return 0.0
    t≥5  && return 0.0

    ts, y = _solve_painleve_ii(t0, t)

    ΔF2=exp(-y[end-1][3]) - exp(-y[end][3]) # the cumulative distribution
    f2=ΔF2/(ts[end]-ts[end-1])              # the density at t
end
pdf(d::Type{TracyWidom}, t::Real) = pdf(d(), t)

"""
Computes the Tracy-Widom distribution by directly solving the
Painlevé II equation by numerical integration
"""
function cdf{S<:Real}(d::TracyWidom, t::S, t0::S = convert(S, -8.0))
    t≤t0 && return 0.0
    t≥5  && return 1.0

    ts, y = _solve_painleve_ii(t0, t)

    F2=exp(-y[end][3])
end
cdf(d::Type{TracyWidom}, t::Real) = cdf(d(), t)

function _solve_painleve_ii{S<:Real}(t0::S, t::S)
    deq(t, y) = [y[2], t*y[1]+2y[1]^3, y[4], y[1]^2]
    a0 = airy(t0)
    T = typeof(a0)
    y0=T[a0, airy(1, t0), 0, airy(t0)^2]    # Initial conditions
    (ts, y)=ode23(deq, y0, [t0, t])         # Solve the ODE
end

"""
Samples the largest eigenvalue of the n x n GUE matrix
"""
function rand(d::TracyWidom, n::Integer)
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
