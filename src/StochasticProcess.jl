#Matrices related to stochastic processes

abstract StochasticProcess
export evolve

immutable BrownianProcess <: StochasticProcess
end

# Generates Brownian motion
evolve(p::BrownianProcess, dx::Real, xend::Real) = cumsum(WhiteNoiseProcess(dx, xend))

# Generates a white noise process
function WhiteNoiseProcess(dx::Real, xend::Real)
  x=[0:dx:xend]
  randn(length(x),1)*sqrt(dx)
end

################
# Airy process #
################

immutable AiryProcess <: StochasticProcess
  beta::Real
end

# Calculates the largest eigenvalue of a stochastic Airy process
# with Brownian noise
function evolve(p::AiryProcess, dx::Real, xend::Real)
    x=[0:dx:xend]
    N=length(x)

    #Discretized Airy operator
    a=-(2/dx^2)*ones(N) - x
    b=+(1/dx^2)*ones(N-1)
    #Plus noise
    dW=WhiteNoiseProcess(dx, xend)
    a+=(2/sqrt(p.beta))*dW/dx

    maxeig(SymTridiagonal(a,b))
end

#Sample program
#t=10000 #number of trials
#v=[AiryProcess(0.001, 10, 2) for i=1:t]
#binsize=.2
#grid=[-5:binsize:2]
#x=hist(v, grid)
#for i=1:length(grid)
#    @printf("%10.5f  %8d  %s\n",grid[i],x[i],repeat("*",int(x[i]/(t*binsize)*200)))
#end
 
