#Matrices related to stochastic processes

import Base: iterate
export AiryProcess, BrownianProcess, WhiteNoiseProcess, next!

abstract type StochasticProcess{T<:Real} end

struct WhiteNoiseProcess{T<:Real} <: StochasticProcess{T}
    dt::T
end

struct BrownianProcess{T<:Real} <: StochasticProcess{T}
    dt::T
end

struct AiryProcess{S<:Real, T<:Real} <: StochasticProcess{T}
    dt::T
    beta::S
end

done(p::StochasticProcess{T}, x...) where {T} = false #Processes can always go on forever

#######################
# White noise process #
#######################

iterate(p::WhiteNoiseProcess, state=()) = (randn()*sqrt(p.dt), ())

####################
# Brownian process #
####################

function iterate(p::BrownianProcess{T}, state=zero(T)) where {T}
    newx = state + randn()*sqrt(p.dt)
    return newx, newx
end

################
# Airy process #
################

"""
Like next, but update only the state of the AiryProcess

Skip the eigenvalue computation, which gets expensive
"""
function next!(p::AiryProcess{T}, S::SymTridiagonal{T}=SymTridiagonal(T[-(2/p.dt^2)], T[])) where {T} # TODO maybe change this to iterate!
    t = (size(S, 1)-1)*p.dt

    #Discredited Airy operator plus diagonal noise
    x = inv(p.dt^2)
    push!(S.dv, -2x - t + 2/sqrt(p.dt*p.beta)*randn())
    push!(S.ev, x)

    return S
end
function iterate(p::AiryProcess{T}, S::SymTridiagonal{T}= SymTridiagonal(T[-(2/p.dt^2)], T[])) where {T}
    t = (size(S, 1)-1)*p.dt

    #Discredited Airy operator plus diagonal noise
    x = inv(p.dt^2)
    push!(S.dv, -2x - t + 2/sqrt(p.dt*p.beta)*randn())
    push!(S.ev, x)

    return (eigmax(S), S)
end
