#unitarylis.jl
#Algorithm 2.5 of Random Eigenvalues by Alan Edelman

#Experiment:    Generate random orthogonal/unitary matrices
#Theory:        Counts longest increasing subsequence statistics

## Parameters
t = 200000      # Number of trials
n = 4           # permutation size
k = 2           # length of longest increasing subsequence

include("patiencesort.jl")

function unitarylis(t, n, k)
    v = zeros(t) # samples
    ## Experiment
    for i = 1:t
        X = QR(complex(randn(k,k), randn(k,k)))[:Q]
        X = X*diagm(sign(complex(randn(k),randn(k))))
        v[i] = abs(trace(X)) ^ (2n)
    end
    z = mean(v)
    c = 0
    for i=1:factorial(n)
        c = c + int(patiencesort(nthperm([1:n],i))<=k)
    end
    return (z, c)
end

println(unitarylis(t, n, k))
