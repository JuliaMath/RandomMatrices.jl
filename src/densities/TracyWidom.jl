export TracyWidom

"""
Tracy-Widom distribution

The probability distribution of the normalized largest eigenvalue of a random
Hermitian matrix.

The cdf of Tracy-Widom is given by

``
F_2 (s) = lim_{n→∞} Pr(√2 n^{1/6} (λₙ - √(2n) ≤ s)
``

References:

1. doi:10.1016/0370-2693(93)91114-3
2. doi:10.1007/BF02100489

Numerical routines adapted from Alan Edelman's course notes for MIT 18.338,
Random Matrix Theory, 2016.
"""
immutable TracyWidom <: ContinuousUnivariateDistribution end


"""
Probability density function of the Tracy-Widom distribution

Computes the Tracy-Widom distribution by directly solving the
Painlevé II equation using the ode23 numerical integrator

# Arguments
* `d::TracyWidom` or `Type{TracyWidom}`: an instance of `TracyWidom` or the type itself
* `t::Real`: The point at which to evaluate the pdf
* `t0::Real = -8.0`: The point at which to start integrating
"""
function pdf{S<:Real}(d::TracyWidom, t::S, t0::S = convert(S, -8.0))
    t≤t0 && return 0.0
    t≥5  && return 0.0

    sol = _solve_painleve_ii(t0, t)

    ΔF2=exp(-sol[end-1][3]) - exp(-sol[end][3]) # the cumulative distribution
    f2=ΔF2/(sol.t[end]-sol.t[end-1])            # the density at t
end
pdf(d::Type{TracyWidom}, t::Real) = pdf(d(), t)

"""
Cumulative density function of the Tracy-Widom distribution

Computes the Tracy-Widom distribution by directly solving the
Painlevé II equation using the ode23 numerical integrator

See `pdf(::TracyWidom)` for a description of the arguments.
"""
function cdf{S<:Real}(d::TracyWidom, t::S, t0::S = convert(S, -8.0))
    t≤t0 && return 0.0
    t≥5  && return 1.0
    sol = _solve_painleve_ii(t0, t)
    F2=exp(-sol[end][3])
end
cdf(d::Type{TracyWidom}, t::Real) = cdf(d(), t)

# An internal function which sets up the Painleve II differential equation and
# runs it through the ode23 numerical integrator
function _solve_painleve_ii{S<:Real}(t0::S, t::S)
    function deq(t, y, dy)
        dy[1] = y[2]
        dy[2] = t*y[1]+2y[1]^3
        dy[3] = y[4]
        dy[4] = y[1]^2
    end
    a0 = airy(t0)
    T = typeof(big(a0))
    y0=T[a0, airy(1, t0), 0, airy(t0)^2]    # Initial conditions
    prob = ODEProblem(deq,y0,(t0,t))
    solve(prob, Vern8(), abstol=1e-12, reltol=1e-12)         # Solve the ODE
end

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
