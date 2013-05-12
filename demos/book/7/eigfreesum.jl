#eigfreesum.jl
#Algorithm 7.1 of Random Eigenvalues by Alan Edelman

#Experiment:    n by n diagonal matrices with random entries +1 or -1
#Plot:          classical and free sums
#Theory:        arcsine distribution

## Parameters
t = 1           # number of trials
n = 1000        # size of matrix
dx = .1

## Experiment
function eigfreesum_experiment(t, n, dx)
    v = Float64[]
    w = Float64[]
    for i = 1:t
        a = sign(randn(n))
        b = sign(randn(n))
        Q = full(QR(randn(n,n))[:Q])                    # random orthogonal matrix
        append!(v, a+b)                                 # classical sum
        append!(w, eigvals(diagm(a) + Q'*diagm(b)*Q))   # free sum
    end
    x = -2-dx/2:dx:2+dx/2
    _,countw = hist(w, x)
    _,countv = hist(v, x)
    x, countw, countv
end
grid, countw, countv = eigfreesum_experiment(t, n, dx)

## Theory
x = -2+dx/20:dx/10:2-dx/20
y = (1/pi) ./ sqrt(4 - x.^2)

## Plot
using Winston
p = FramedPlot()
hw = Histogram(countw*pi/(n*t*dx), step(grid), "color", "blue")
hw.x0 = first(grid)
add(p, hw)
hv = Histogram(countv*pi/(n*t), step(grid), "color", "green")
hv.x0 = first(grid)
add(p, hv)
add(p, Curve(x, y, "color", "red", "linewidth", 2))
if isinteractive()
    Winston.display(p)
else
    file(p, "eigfreesum.png")
end
