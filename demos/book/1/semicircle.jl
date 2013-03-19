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
## Experiment
for i = 1:t
    a=randn(n,n)      # n by n matrix of random Gaussians
    s=(a+a')/sqrt(2*n)# symmetrize and normalize matrix
    v=[v;eig(s)[1]]   # eigenvalues
end
count=hist(v,-2.4:dx:2.4)
count/=n*t*dx
## Theory
x=[-2:dx:2]
y=sqrt(4-x.^2)/(2*pi)

## Plot
using Winston
p = FramedPlot()
add(p, Histogram(count, dx))
add(p, Curve(x+2.4+dx/2, y, "color", "blue", "linewidth", 2)) 
file(p, "semicircle.png")
