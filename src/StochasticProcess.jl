#Matrices related to stochastic processes

import Base: start, next, done
export AiryProcess, BrownianProcess, WhiteNoiseProcess

abstract StochasticProcess{T<:Real}

immutable WhiteNoiseProcess{T<:Real} <: StochasticProcess{T}
    dt::T
end

immutable BrownianProcess{T<:Real} <: StochasticProcess{T}
    dt::T
end

type AiryProcess{S<:Real, T<:Real} <: StochasticProcess{T}
    dt::T
    beta::S
end

done{T}(p::StochasticProcess{T}) = false #Processes can always go on forever

#######################
# White noise process #
#######################

start{T}(p::WhiteNoiseProcess{T}) = nothing
next{T}(p::WhiteNoiseProcess{T}, x::Void=nothing) = (randn()*sqrt(p.dt), nothing)

####################
# Brownian process #
####################

start{T}(p::BrownianProcess{T}) = zero(T)
function next{T}(p::BrownianProcess{T}, x::T)
    newx = x + randn()*sqrt(p.dt)
    newx, newx
end

################
# Airy process #
################
start{T}(p::AiryProcess{T}) = SymTridiagonal(T[-(2/p.dt^2)], T[])
function next{T}(p::AiryProcess{T}, S::SymTridiagonal{T})
    t = (size(S, 1)-1)*p.dt

    #Discretized Airy operator plus diagonal noise
    x = inv(p.dt^2)
    push!(S.dv, -2x - t + 2/sqrt(p.beta/p.dt)*randn())
    push!(S.ev, x)

    (eigmax(S)*p.dt^2, S)
end

