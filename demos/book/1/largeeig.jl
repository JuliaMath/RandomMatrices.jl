#largeeig.png
#Algorithm 1.4 of Random Eigenvalues by Alan Edelman

#Experiment:  Largest eigenvalue of random Hermitian matrices
#Plot:        Histogram of the normalized largest eigenvalues
#Theory:      Tracy-Widom as n->infinity

include("tracywidom.jl")

function largeeig_experiment(
    ## Parameters
    n=100,       # matrix size
    t=5000,      # trials
    dx=.2,       # binsize
)
    ## Experiment
    v=Float64[] # eigenvalue samples
    for i=1:t
        a=randn(n,n)+im*randn(n,n) # random nxn complex matrix
        s=(a+a')/2                 # Hermitian matrix
        push!(v,max(eig(s)[1]))    # Largest Eigenvalue
    end
    v=n^(1/6)*(v-2*sqrt(n))        # normalized eigenvalues

    grid=-5:dx:2
    count=hist(v,grid)[2]/(t*dx)
    (grid, count)
end
(grid, count) = largeeig_experiment()

## Theory
(t, f) = tracywidom(5., -8.)

## Plot
using Winston
p = FramedPlot()
h = Histogram(count, step(grid))
h.x0 = first(grid)
add(p, h)
add(p, Curve(t, f, "linewidth", 2, "color", "blue"))
if isinteractive()
    Winston.display(p)
else
    file(p, "largeeig.png")
end

