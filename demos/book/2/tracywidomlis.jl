#tracywidomlis.jl
#Algorithm 2.7 of Random Eigenvalues by Alan Edelman

#Experiment:    Sample random permutations
#Plot:          Histogram of lengths of longest increasing subsequences
#Theory:        The Tracy-Widom law

## Parameters
t = 10000       # number of trials
n = 6^6         # length of permutations
dx = 1/6        # bin size

## Experiment
include("patiencesort.jl")
function tracywidomlis(t, n, dx)
    v = zeros(t) # samples
    for i=1:t
        v[i] = patiencesort(randperm(n))
    end
    w = (v-2sqrt(n))/n^(1/6)
    return hist(w, -5:dx:2)
end
grid, count = tracywidomlis(t, n, dx)

## Plot
h = Histogram(count/(t*dx), step(grid))
h.x0 = first(grid)
include("../1/tracywidom.jl")
add(p, h)
if isinteractive()
    Winston.display(p)
else
    file(p, "tracywidomlis.png")
end
