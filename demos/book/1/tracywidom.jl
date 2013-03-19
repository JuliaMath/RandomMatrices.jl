#tracywidom.m
#Algorithm 1.3 of Random Eigenvalues by Alan Edelman

#Theory:      Compute and plot the Tracy-Widom distribution

##Parameters
t0=5.         # right endpoint
tn=-8.        # left  endpoint

##Theory: The differential equation solver
#Compute Tracy-Widom distribution
using ODE
include("gradient.jl")
function tracywidom(t0::Float64, tn::Float64)
    function deq(t::Float64, y::Vector{Float64})
        yout = [y[2]; t*y[1]+2y[1]^3; y[4]; y[1]^2]
    end
    
    y0=[airy(t0); airy(1, t0); 0; airy(t0)^2] # Initial conditions
    (t, y)=ode23(deq, [t0, tn], y0)           # Solve the ODE
    F2=exp(-y[:,3])                           # the distribution
    f2=gradient(F2, t)                        # the density
    return (t, f2)
end
(t, f2) = tracywidom(t0, tn)

## Plot
using Winston
p = FramedPlot()
add(p, Curve(t, f2, "linewidth", 2))
file(p, "tracywidom.png")

