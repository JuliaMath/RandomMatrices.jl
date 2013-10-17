#patiencesort.jl
#Algorithm 2.5 of Random Eigenvalues by Alan Edelman

function patiencesort(p)
    piles = similar(p, 0)
    for pi in p
        idx = 1
        for pp in piles
            if pi > pp
                idx += 1
            end
        end
        if idx > length(piles)
            d = length(piles)
            resize!(piles, idx)
            piles[d+1:end] = 0
        end
        piles[idx] = pi
    end
    return length(piles)
end
