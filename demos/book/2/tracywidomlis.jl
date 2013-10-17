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
require("patiencesort.jl")
function tracywidomlis(t, n, dx)
    ## single processor: 1x speed
    #v = [patiencesort(randperm(n)) for i = 1:t]

    ## simple parallelism: Nx speed * 130%
    v = pmap((i)->patiencesort(randperm(n)), 1:t)

    ## maximum parallelism: Nx speed
    #grouped = floor(t/(nprocs()-2))
    #println(grouped)
    #v = vcat(pmap((i)->[patiencesort(randperm(n)) for j = 1:grouped], 1:t/grouped)...)

    w = (v-2sqrt(n))/n^(1/6)
    return hist(w, -5:dx:2)
end
@time grid, count = tracywidomlis(t, n, dx)

## Plot
using Winston
p = FramedPlot()
h = Histogram(count/(t*dx), step(grid))
h.x0 = first(grid)
include("../1/tracywidom.jl")
add(p, h)
if isinteractive()
    Winston.display(p)
else
    file(p, "tracywidomlis.png")
end
