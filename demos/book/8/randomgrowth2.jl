#randomgrowth2.jl
using Winston
using ODE
include("../1/gradient.jl")

function randomgrowth2()
    num_trials = 2000
    g = 1
    q = 0.7
    Ns = 80 # [50, 100, 200]
    # [-1.771086807411  08131947928329]
    # Gap = c d^N
    for jj = 1:length(Ns)
        N = Ns[jj]
        M = round(N*g)
        B = zeros(num_trials)
        @time begin
            for ii = 1:num_trials
                # G = G[M, N, q]
                # B[ii] = G[M,N]
                B[ii] = G(M, N, q)
            end
        end
        C = (B - N*om(g,q)) / (sigma_growth(g,q)*N^(1/3))
        d = 0.2
        x, n = hist(C - exp(-N*0.0025 - 3), -6:d:4)
        p = FramedPlot()
        h = Histogram(n/(d*num_trials), step(x))
        h.x0 = first(x)
        add(p, h)

        ## Theory
        t0 = 4
        tn = -6
        dx = 0.005
        deq = (t, y) -> [y[2]; t*y[1]+2*y[1]^3; y[4]; y[1]^2]
        y0 = [airy(t0); airy(1,t0); 0; airy(t0)^2]      # boundary conditions
        t, y = ode23(deq, t0:-dx:tn, y0)                # solve
        F2 = exp(-y[:,3])                               # the distribution
        f2 = gradient(F2, t)                            # the density
       
        add(p, Curve(t, f2, "color", "red", "linewidth", 3))
        Winston.display(p)
        println(mean(C))
    end
end

function G(N, M, q)
    # computes matrix G[N,M] for a given q
    GG = floor(log(rand(N,M))/log(q))
    #float(rand(N,M) < q)
    
    # Compute the edges: the boundary
    for ii = 2:N
        GG[ii, 1] = GG[ii, 1] + GG[ii-1, 1]
        GG[1, ii] = GG[1, ii] + GG[1, ii-1]
    end
    # Compute the inside
    for ii = 2:N
        for jj = 2:M
            GG[ii, jj] = GG[ii, jj] + max(GG[ii, jj-1], GG[ii-1, jj])
        end
    end
    ## Plot
    # imagesc((0,N),(M,0),GG)
    # sleep(.01)
    
    return  GG[N, M]
end

# helper functions sigma.m and om, which describe mean and
# standard deviation of G[N, M]

# Computes G[qN,N]/N for N->inf
om(g,q) = ((1 + sqrt(g * q)) ^ 2) / (1 - q) - 1

sigma_growth(r, q) =
    q^(1/6) * r^(1/6) / (1 - q) *
    (sqrt(r) + sqrt(q))^(2/3) *
    (1 + sqrt(q * r))^(2/3)


randomgrowth2()
