#bellcurve.jl
#Algorithm 1.1 of Random Eigenvalues by Alan Edelman

#Experiment:  Generate random samples from the normal distribution
#Plot:        Histogram of random samples
#Theory:      The normal distribution curve

## Experiment
t = 1000000
dx = .2
v = randn(t)
x = -4:dx:4
grid,count = hist(v,x)

## Theory
y = exp(-x.^2/2)/sqrt(2*pi)

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
    file(p, "bellcurve.png")
end

