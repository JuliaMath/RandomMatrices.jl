#largeeig.png
#Algorithm 1.4 of Random Eigenvalues by Alan Edelman

#Experiment:  Largest eigenvalue of random Hermitian matrices
#Plot:        Histogram of the normalized largest eigenvalues
#Theory:      Tracy-Widom as n->infinity

include("tracywidom.jl")

## Parameters
n=100        # matrix size
t=5000       # trials
v=[]         # eigenvalue samples
dx=.2        # binsize
## Experiment
for i=1:t
    a=randn(n,n)+im*randn(n,n) # random nxn complex matrix
    s=(a+a')/2                 # Hermitian matrix
    v=[v; max(eig(s)[1])]         # Largest Eigenvalue
end
v=n^(1/6)*(v-2*sqrt(n))        # normalized eigenvalues

count=hist(v,-5:dx:2)
count/=t*dx

## Theory
(t, f) = tracywidom(5., -8.)

## Plot
#XXX Currently Histogram() must start plotting from 0.
#For now, add offset to theoretical curve
#Filed: https://github.com/nolta/Winston.jl/issues/28
#cjh 2013-03-16
using Winston
p = FramedPlot()
add(p, Histogram(count, dx))
add(p, Curve(t+5, f, "linewidth", 2, "color", "blue")) 
file(p, "largeeig.png")

