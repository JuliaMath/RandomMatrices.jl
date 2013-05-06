#spherecomponent.jl
#Algorithm 2.3 of Random Eigenvalues by Alan Edelman

#Experiment:    Random components of a vector on a sphere
#Plot:          Histogram of a single component
#Theory:        Beta distribution

## Parameters
t = 100000      # trials
n = 12          # dimensions of sphere
dx = 0.1        # bin size

## Experiment
v = randn(n, t)
v = (v[1,:] ./ sqrt(sum(v.^2,1)))'
grid, count = hist(v, -1:dx:1)

## Theory
x = grid
y = (gamma(n-1)/(2^(n-2)*gamma((n-1)/2)^2))*(1-x.^2).^((n-3)/2)

## Plot
using Winston
p = FramedPlot()
h = Histogram(count/(t*dx), step(grid))
h.x0 = first(grid)
add(p, h)
add(p, Curve(x, y, "color", "blue", "linewidth", 2))
if isinteractive()
    Winston.display(p)
else
    file(p, "spherecomponent.png")
end
