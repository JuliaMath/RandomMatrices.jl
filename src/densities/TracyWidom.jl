export TracyWidom

include("gradient.jl")

immutable TracyWidom <: ContinuousUnivariateDistribution end


"""
Computes the Tracy-Widom distribution by directly solving the
Painlevé II equation by numerical integration
"""
function pdf(d::TracyWidom, t::Real)
    t0 = -8.0
    t≤t0 && return 0.0
    t≥5  && return 0.0

    deq(t, y) = [y[2]; t*y[1]+2y[1]^3; y[4]; y[1]^2]

    y0=[airy(t0); airy(1, t0); 0; airy(t0)^2] # Initial conditions
    (ts, y)=ode23(deq, [t0, t], y0)           # Solve the ODE
    F2=exp(-y[:,3])                           # the cumulative distribution
    f2=(F2[end]-F2[end-1])/(ts[end]-ts[end-1])# the density at t
end
pdf(d::Type{TracyWidom}, t::Real) = pdf(d(), t)

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
