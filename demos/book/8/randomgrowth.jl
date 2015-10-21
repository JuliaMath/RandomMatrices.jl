#randomgrowth.jl
using Winston

function random_growth(M, N, q)
    G = zeros(M, N)
    T = 1000

    G[1,1] = true
    # update the possible sets
    display(imagesc((0,N), (M,0), G))
    for t = 1:T
        sets = next_possible_squares(G)
        ## Highlights all the possible squares
        for i = 1:length(sets)
            idx = sets[i]::(Int,Int)
            G[idx[1], idx[2]] = 0.25
        end
        display(imagesc((0,N), (M,0), G))
        sleep(.01)

        ## Actual growth
        for i = 1:length(sets)
            ison = 0.5 * (rand() > (1-q))
            if ison > 0
                idx = sets[i]::(Int,Int)
                G[idx[1], idx[2]] = ison
            end
        end
        display(imagesc((0,N), (M,0), G))
        G[G .== 0.5] = 1
        G[G .== 0.25] = 0
        sleep(.01)
    end
    return G
end

function next_possible_squares(G)
    M, N = size(G)::(Int,Int)
    sets = Array((Int,Int),0)
    for ii = 1:M
        for jj = 1:N
            if G[ii, jj] == 0
                if  (ii == 1 && jj > 1 && G[ii, jj-1] == 1) ||
                    (jj == 1 && ii > 1 && G[ii-1, jj] == 1) ||
                    (ii > 1 && jj > 1 && G[ii-1, jj] == 1 && G[ii, jj-1] == 1)
                    push!(sets, (ii,jj))
                end
            end
        end
    end
    sets
end
