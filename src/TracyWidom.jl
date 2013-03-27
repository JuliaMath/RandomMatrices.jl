#Computes the Tracy-Widom distribution by directly solving the
#Painlevé II equation
using ODE
include("gradient.jl")
function tracywidom(t::Real)
    t0 = -8.0

    t>t0 ? return 0.0 : nothing
    t>5  ? return 0.0 : nothing

    deq(t, y) = [y[2]; t*y[1]+2y[1]^3; y[4]; y[1]^2]

    y0=[airy(t0); airy(1, t0); 0; airy(t0)^2] # Initial conditions
    (ts, y)=ode23(deq, [t0, t], y0)           # Solve the ODE
    F2=exp(-y[:,3])                           # the cumulative distribution
    f2=(F2[end]-F2[end-1])/(ts[end]-ts[end-1])# the density at t
end

