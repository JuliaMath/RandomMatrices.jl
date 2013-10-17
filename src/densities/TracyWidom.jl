#Computes the Tracy-Widom distribution by directly solving the
#Painlevé II equation
importall Distributions
export TracyWidom
using ODE

include("gradient.jl")

immutable TracyWidom <: ContinuousUnivariateDistribution
end

function pdf(d::TracyWidom, t::Real)
  t0 = -8.0
  t>t0 ? (return 0.0) : nothing
  t>5  ? (return 0.0) : nothing

  deq(t, y) = [y[2]; t*y[1]+2y[1]^3; y[4]; y[1]^2]

  y0=[airy(t0); airy(1, t0); 0; airy(t0)^2] # Initial conditions
  (ts, y)=ode23(deq, [t0, t], y0)           # Solve the ODE
  F2=exp(-y[:,3])                           # the cumulative distribution
  f2=(F2[end]-F2[end-1])/(ts[end]-ts[end-1])# the density at t
end


#Samples the largest eigenvalue of the n x n GUE matrix
function rand(d::TracyWidom, n::Integer)
  k=int(n-10*n^(1/3)-1)
  A=[chi(i) for i=(n-1):-1:k]
  B=randn(n-k+1)
  v=eigmax(SymTridiagonal(B, A))
end

