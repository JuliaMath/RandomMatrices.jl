#bellcurve.jl
#Algorithm 1.1 of Random Eigenvalues by Alan Edelman

#Experiment:  Generate random samples from the normal distribution
#Plot:        Histogram of random samples
#Theory:      The normal distribution curve

## Experiment
t=1000000
dx=.2
v=randn(t)
grid=[-4:dx:4]
count=hist(v,grid)/(t*dx)

## Theory
x=grid
y=exp(-x.^2/2)/sqrt(2*pi)

## Plot
using Winston
p = FramedPlot()
h = Histogram(count, dx)
h.x0 = grid[1]
add(p, h)
add(p, Curve(x, y, "color", "blue", "linewidth", 2)) 
file(p, "bellcurve.png")
