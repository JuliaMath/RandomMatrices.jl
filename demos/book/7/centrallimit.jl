#centrallimit.jl
#Algorithm 7.2 of Random Eigenvalues by Alan Edelman

using Winston

function centrallimit(m, t)
    # sums m random matrices and histograms their eigenvalues
    # normalized by 1/sqrt(m)
    # Input :
    #  m : number of matrices to sum
    #  t : number of trials to perform
    n = 100     # matrix size
    dx = .1     # bin size
    v = Float64[]
    # we choose a random diagonal {-1,1} matrx as our prototype
    # A = diagm(sign(randn(1,n)))
    # proto = kron([1 0; 0 -1], eye(n/2))
    for i = 1:t                         # for each trial
        B = zeros(n,n)                  # initialize the accumulator
        for z = 1:m                     # sum up the random matrices
            Q = full(QR(randn(n,n) + im*randn(n,n))[:Q])  # (piecewise) Haar orthogonal
            A = diagm(sign(randn(n,1)))
            B = B + Q'*A*Q              # A and Q'AQ are asymptotically free
        end
        append!(v, real(eigvals(B)/sqrt(m)))  # normalize eigenvalues
    end
    grid, count = hist(v, -2-dx/2:dx:2)

    p = FramedPlot()
    h = Histogram(count*pi/2/(dx*n*t), step(grid))
    h.x0 = first(grid)
    add(p, h)
    if isinteractive()
        Winston.display(p)
    else
        file(p, "centrallimit_m=$(m)_t=$(t).png")
    end
end
