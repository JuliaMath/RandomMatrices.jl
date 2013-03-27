#semicircle.m
#Algorithm 1.2 of Random Eigenvalues by Alan Edelman

#Experiment:  Sample random symmetric Gaussian matrices
#Plot:        Histogram of the eigenvalues
#Theory:      Semicircle as n->infinity
## Parameters
n=1000       # matrix size
t=10         # trials
v=[]         # eigenvalue samples
dx=.2        # binsize
grid=[-2.4:dx:2.4]
## Experiment
for i = 1:t
    a=randn(n,n)      # n by n matrix of random Gaussians
    s=(a+a')/sqrt(2*n)# symmetrize and normalize matrix
    v=[v;eig(s)[1]]   # eigenvalues
end
count=hist(grid)/(n*t*dx)
## Theory
x=[-2:dx:2]
y=sqrt(4-x.^2)/(2*pi)

## Plot
using Winston
p = FramedPlot()
h = Histogram(count, dx)
h.x0 = grid[1]
add(p, h)
add(p, Curve(x+2.4+dx/2, y, "color", "blue", "linewidth", 2)) 
file(p, "semicircle.png")
