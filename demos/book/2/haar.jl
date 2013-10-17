#haar.jl
#Algorithm 2.4 of Random Eigenvalues by Alan Edelman

#Experiment:    Generate random orthogonal/unitary matrices
#Plot:          Histogram eigenvalues
#Theory:        Eigenvalues are on unit circle

## Parameters
t = 5000        # trials
dx = 0.05       # binsize
n = 10          # matrix size

function haar_experiment(t,dx,n)
    v = Complex{Float64}[] # eigenvalue samples
    for i = 1:t
        # Sample random orthogonal matrix
        # X = QR(randn(n,n))[:Q]
        #  If you have non-uniformly sampled eigenvalues, you may need this fix:
        #X = X*diagm(sign(randn(n,1)))
        #
        # Sample random unitary matrix
        X = QR(randn(n,n)+im*randn(n,n))[:Q]
        #  If you have non-uniformly sampled eigenvalues, you may need this fix:
        #X = X*diagm(sign(randn(n)+im*randn(n)))
        append!(v, eigvals(full(X)))
    end
    x = ((-1+dx/2):dx:(1+dx/2))*pi
    x, v
end
x, v = haar_experiment(t,dx,n)
grid, count = hist(angle(v), x)

## Theory
h2 = t*n*dx/2*x.^0

## Plot
using Winston
p = FramedPlot()
h = Histogram(count, step(grid))
h.x0 = first(grid)
add(p, h)
add(p, Curve(x, h2, "color", "blue", "linewidth", 2)) 
if isinteractive()
    Winston.display(p)
else
    file(p, "haar.png")
end
