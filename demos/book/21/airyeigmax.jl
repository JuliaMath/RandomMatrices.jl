#airyeigmax.jl
#Algorithm <21.1> of Random Eigenvalues by Alan Edelman

#Experiment:    Largets eigenvalue of a Stochastic Airy Operator
#Plot:          Histogram of the largest eigenvalues
#Theory:        The Tracy-Widom Law

## Parameters   # Include the most interesting configurable values at the top, with brief descriptions
t = 10000       # trials
n = 1e9         # level of discretization
beta = 2
h = n^(-1/3)    # h serves as dx

function airyeigmax_experiment(t,n,beta,h) # wrap the experimental logic in a function to enable faster JIT
    ## Experiment
    v = zeros(t)    # samples
    x = 0:h:10      # discretization of x
    N = length(x)

    #Parallel experiment
    b = (1 / h^2) * ones(N-1)
    v = pmap(1:t) do i
        ## discretize stochastic airy operator
        # discretize airy operator
        a = -(2 / h ^ 2) * ones(N);  # differential operator: d^2 / dx^2
        a -= x                       # d^2 / dx^2 - x
        # add the stochastic part
        dW = randn(N) * sqrt(h)
        a += (2 / sqrt(beta)) * dW / h
        ## calculate the largest eigenvalue of tridiagonal matrix T
        # diagonal of T: a
        # subdiagonal of T: b
        return eigmax(SymTridiagonal(a,b))
    end

    binsize = 1/6
    return hist(v, -6:binsize:6)
end
grid, count = airyeigmax_experiment(t,n,beta,h) #run the experiment, making the global variables local for speed

## Plot
include("../1/tracywidom.jl")
h = Histogram(count/(sum(count)*step(grid)), step(grid))    # add the histogram of data
h.x0 = first(grid)
add(p, h)

if isinteractive()
    Winston.display(p)
else
    file(p, "<filename>.png")
end
