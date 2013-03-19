#bellcurve.jl
#Algorithm 1.1 of Random Eigenvalues by Alan Edelman

#Experiment:  Generate random samples from the normal distribution
#Plot:        Histogram of random samples
#Theory:      The normal distribution curve

## Experiment
t=1000000
dx=.2
v=randn(t)
count=hist(v,[-4-dx/2:dx:4+dx/2])
count/=sum(count)*dx

## Theory
x=[-4:dx:4]
y=exp(-x.^2/2)/sqrt(2*pi)

## Plot
#XXX Currently Histogram() must start plotting from 0.
#For now, add offset to theoretical curve
#Filed: https://github.com/nolta/Winston.jl/issues/28
#cjh 2013-03-16
using Winston
p = FramedPlot()
add(p, Histogram(count, dx))
add(p, Curve(x+4+dx/2, y, "color", "blue", "linewidth", 2)) 
file(p, "bellcurve.png")
