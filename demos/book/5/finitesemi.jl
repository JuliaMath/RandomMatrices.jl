#finitesemi.jl
#Algorithm 5.6 of Random Eigenvalues by Alan Edelman

#Experiment:    Eigenvalues of GUE matrices
#Plot:          Histogram of eigenvalues
#Theory:        Semicircle and finite semicircle

## Parameters
n = 3           # size of matrix
s = 10000       # number of samples
d = 0.1         # bin size

## Experiment
function finitesemi_experiment(n,s,d)
    e = Float64[]
    for i = 1:s
        a = randn(n,n)+im*randn(n,n)
        a = (a+a')/(2*sqrt(4n))
        v = eigvals(a)
        append!(e, v)
    end
    return hist(e, -1.5:d:1.5)
end
grid, count = finitesemi_experiment(n,s,d)

## Theory
x = -1:0.01:1
y = sqrt(1-x.^2)

## Plot
using Winston
p = FramedPlot()
h = Histogram(count*pi/(2*d*n*s), step(grid))
h.x0 = first(grid)
add(p, h)
add(p, Curve(x, y, "color", "red", "linewidth", 2))
include("levels.jl")
add(p, Curve(levels(n)..., "color", "blue", "linewidth", 2))
if isinteractive()
    Winston.display(p)
else
    file(p, "spherecomponent.png")
end

