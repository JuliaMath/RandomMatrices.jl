#mpexperiment.jl
#Algorithm 4.2 of Random Eigenvalues by Alan Edelman

#Experiment:    Marcenko-Pastur for applications
#Plot:          Histogram eigenvalues of the covariance matrix
#Theory:        Marcenko-Pastur as n->Infinity

## Parameters
t = 100         # trials
r = 1           # aspect ratio
n = 100         # matrix column size
dx = 0.1        # binsize

function mp_experiment(n,r,t,dx)
    m = iround(n/r)
    ## Experiment
    v = Float64[]          # eigenvalue samples
    for i = 1:t
        A = randn(m,n) + 4*sqrt(n)*diagm(((1:n).<10))
        A = A + sqrt(n) * diagm((1:n).>(n-1)) * 3 #3+0.1*randn(n,1) #3/sqrt(n)
        append!(v, svd(A)[2]) # eigenvalues
    end
    v = v / sqrt(m) # normalized eigenvalues
    ## Theory
    a = 1-sqrt(r)
    b = 10
    x = a-dx/2:dx:b
    y = real(sqrt((x.^2-a^2).*(2^2-x.^2) + 0im)./(pi*x*r))
    return (hist(v, x), (x, y))
end
((grid, count), (x,y)) = mp_experiment(n,r,t,dx)

## Plot
using Winston
p = FramedPlot()
h = Histogram(count/(t*n*dx), step(grid))
h.x0 = first(grid)
add(p, h)
add(p, Curve(x, y, "linewidth", 2, "color", "blue"))
last = length(count)-2
while count[last] == 0
    last -= 1
end
setattr(p, "xrange", (0, grid[last+2]))
setattr(p, "yrange", (-.1, ceil(max(count/(t*n*dx)))))
if isinteractive()
    Winston.display(p)
else
    file(p, "mpexperiment.png")
end
