#semicircle.m
#Algorithm 1.2 of Random Eigenvalues by Alan Edelman

#Experiment:  Sample random symmetric Gaussian matrices
#Plot:        Histogram of the eigenvalues
#Theory:      Semicircle as n->infinity

## Parameters
n=1000       # matrix size
t=10         # trials
dx=.2        # binsize
function semicircle_experiment(n,t,dx)
    ## Experiment
    v=Float64[] # eigenvalue samples
    for i = 1:t
        a=randn(n,n)      # n by n matrix of random Gaussians
        s=(a+a')/sqrt(2*n)# symmetrize and normalize matrix
        append!(v,eig(s)[1])# eigenvalues
    end
    grid=-2.4:dx:2.4
    count=hist(v,grid)[2]/(n*t*dx)
    (grid, count)
end
grid,count = semicircle_experiment(n,t,dx)

## Theory
x=[-2:dx:2]
y=sqrt(4-x.^2)/(2*pi)

## Plot
using Winston
p = FramedPlot()
h = Histogram(count, step(grid))
h.x0 = first(grid)
add(p, h)
add(p, Curve(x, y, "color", "blue", "linewidth", 2)) 
if isinteractive()
    Winston.display(p)
else
    file(p, "semicircle.png")
end
