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
grid=[-5:dx:2]
## Experiment
for i=1:t
    a=randn(n,n)+im*randn(n,n) # random nxn complex matrix
    s=(a+a')/2                 # Hermitian matrix
    v=[v; max(eig(s)[1])]         # Largest Eigenvalue
end
v=n^(1/6)*(v-2*sqrt(n))        # normalized eigenvalues

count=hist(v,grid)/(t*dx)

## Theory
(t, f) = tracywidom(5., -8.)

## Plot
using Winston
p = FramedPlot()
h = Histogram(count, dx)
h.x0 = grid[1]
add(p, h)
add(p, Curve(t+5, f, "linewidth", 2, "color", "blue")) 
file(p, "largeeig.png")

