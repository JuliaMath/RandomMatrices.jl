##Theory: The differential equation solver
#Computes the Tracy-Widom distribution by directly solving the
#Painlevé II equation
using ODE
include("gradient.jl")
function tracywidom(t0::FloatingPoint, tn::FloatingPoint)
    function deq(t::Float64, y::Vector{Float64})
        yout = [y[2]; t*y[1]+2y[1]^3; y[4]; y[1]^2]
    end
    
    y0=[airy(t0); airy(1, t0); 0; airy(t0)^2] # Initial conditions
    (t, y)=ode23(deq, [t0, tn], y0)           # Solve the ODE
    F2=exp(-y[:,3])                           # the distribution
    f2=gradient(F2, t)                        # the density
    return (t, f2)
end


