# Random samples from determinantal point processes

"""
Computes a random sample from the determinantal point process defined by the
spectral factorization object `L`.

Inputs:

    `L`: `Eigen` factorization object of an N x N matrix

Output:

    `Y`: A `Vector{Int}` with entries in [1:N].

References:

    Algorithm 18 of \\cite{HKPV05}, as described in Algorithm 1 of \\cite{KT12}.

    @article{HKPV05,
        author = {Hough, J Ben and Krishnapur, Manjunath and Peres, Yuval and Vir\'{a}g, B\'{a}lint},
        doi = {10.1214/154957806000000078},
        journal = {Probability Surveys},
        pages = {206--229},
        title = {Determinantal Processes and Independence},
        volume = {3},
        year = {2005}
        archivePrefix = {arXiv},
        eprint = {0503110},
    }

    @article{KT12,
        author = {Kulesza, Alex and Taskar, Ben},
        doi = {10.1561/2200000044},
        journal = {Foundations and Trends in Machine Learning},
        number = {2-3},
        pages = {123--286},
        title = {Determinantal Point Processes for Machine Learning},
        volume = {5},
        year = {2012},
        archivePrefix = {arXiv},
        eprint = {1207.6083},
    }

    TODO Check loss of orthogonality - a tip from Zelda Mariet
"""
function rand(L::LinearAlgebra.Eigen{S,T}) where {S<:Real,T}
    N = length(L.values)
    J = Int[]
    for n=1:N
        λ = L.values[n]
        rand() < λ/(λ+1) && push!(J, n)
    end

    V = L.vectors[:, J]
    Y = Int[]
    nV = size(V, 2)
    while true
        # Select i from 𝒴=[1:N] (ground set) with probabilities
        # Pr(i) = 1/|V| Σ_{v∈V} (v⋅eᵢ)²

        #Compute selection probabilities
        Pr = zeros(N)
        for i=1:N
            for j=1:nV #TODO this loop is a bottleneck - why?
                Pr[i] += (V[i,j])^2 #ith entry of jth eigenvector
            end
            Pr[i] /= nV
        end
        @assert abs(1-sum(Pr)) < N*eps() #Check normalization

        #Simple discrete sampler
        i, ρ = N, rand()
        for j=1:N
            if ρ < Pr[j]
                i = j
                break
            else
                ρ -= Pr[j]
            end
        end
        push!(Y, i)
        nV == 1 && break #Done
        #V = V⊥ #an orthonormal basis for the subspace of V ⊥ eᵢ
        V[i, :] = 0 #Project out eᵢ
        V = full(qrfact!(V)[:Q])[:, 1:nV-1]
        nV = size(V, 2)
    end
    return Y
end
