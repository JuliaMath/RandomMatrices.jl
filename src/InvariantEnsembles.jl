# Computes samples of unitary invariant ensembles matrices of size nxn
#
# InvariantEnsemble(str,n) represents an n x n unitary invariant ensemble
# with distribution
#
#   exp(- Tr Q(M)) dM
#
# str specifies an ensemble with precomputed recurrence coefficients.
#  The currently include ensembles are
#
#       Quartic:            Q(M) = n M^4
#       Eight:              Q(M) = n M^8
#       HODecay:            Q(M) = n (M^4/20 - 4/15M^3 +M^2/5 + 8/5M)
#       CoshUnscaled:       Q(M) = cosh(M)
#       QuarticUnscaled:    Q(M) = M^4
#       EightUnscaled:      Q(M) = M^8
#
# References:
#    Olver, Rao & Trogdon 2014 arXiv:1404.007



module InvariantEnsembles
    using RandomMatrices, RandomMatrices.ApproxFun

export InvariantEnsemble





## Construct orthogonal polynomials as Funs from first moment and recurrence relationship
function orthonormalpolynomials(μ0,α::Vector,β::Vector,d)
    n=length(α)+1

    p=Array(Fun{Chebyshev,Float64},n)
    p[1] = Fun(d,[μ0])
    p[2] = multiplybyx(p[1])./β[1] - p[1].*α[1]./β[1]
    for k = 3:n
       p[k] =  multiplybyx(p[k-1])./β[k-1] - p[k-1].*α[k-1]./β[k-1] - p[k-2].*β[k-2]./β[k-1]
    end
    p
end

function orthonormalpolynomials(μ0,α::Function,β::Function,d,n)
    orthonormalpolynomials(μ0,map(α,[1:n-1]),map(β,[1:n-1]),d)
end

function legendrepolynomials(n)
  orthonormalpolynomials(1 ./ sqrt(2), k->0, k->sqrt(k^2 ./ (4k^2-1)), [-1.,1.], n)
end
function scaledhermitepolynomials(n)
  orthonormalpolynomials(n^(1/4)/(1.0π)^(1/4),k->0,k->sqrt((k)/(2n)),[-3.,3.],n)
end


## Evaluate orthogonal polynomials at points x from first moment and recurrence relationship
function orthonormalpolynomialsvalues(μ0,α::Vector,β::Vector,x)
    n=length(α)+1
    m=length(x)

    p=Array(Float64,m,n)

    p[:,1] = μ0*ones(m)
    p[:,2] = x.*p[:,1]./β[1] - p[:,1].*α[1]./β[1]
    for k = 3:n
       p[:,k] =  x.*p[:,k-1]./β[k-1] - p[:,k-1].*α[k-1]./β[k-1] - p[:,k-2].*β[k-2]./β[k-1]
    end
    p
end

## Ensembles

mutable struct InvariantEnsemble{D}
    basis::Array{Float64,2}         # the m x n array of the first n weighted OPs evaluated at m Chebyshev points
    domain::D                       # the sub domain  of the real line
    InvariantEnsemble{D}(b,d) where D = new(b,d)
end

InvariantEnsemble(basis::Matrix,d::Domain) =
    InvariantEnsemble{typeof(d)}(basis,d)

InvariantEnsemble(basis::Matrix,d::Vector) =
    InvariantEnsemble(basis,Segment(d[1],d[2]))

InvariantEnsemble(basis,d::Vector) =
    InvariantEnsemble(basis,Segment(d[1],d[2]))


## n, where the invariant ensemble is n x n
Base.size(ie::InvariantEnsemble)=size(ie.basis,2),size(ie.basis,2)
Base.size(ie::InvariantEnsemble,n)=size(ie.basis,2)



#Takes in list of OPs, constructs phis
# can be unstable for large n
function InvariantEnsemble(p::Array{F},V::Function,d,n::Integer) where {F<:Fun}
    error("Reimplement with matrix")

#     wsq=IFun(x->exp(-n/2.*V(x)),d)
#
#     #We now extend the lengths
#     m=length(wsq)+2n
#     p=map(f->pad(f,m),p)
#
#     # We want to multiply, but in value space to
#     # take advantage of wsq being very small
#
#     wsqv=exp(-n/2.*V(points(wsq.domain,m)))
#     q=map(f->IFun(chebyshevtransform(wsqv.*values(f)),f.domain),p)
#
#     InvariantEnsemble(q)

end






##Adaptively construct invariant ensemble from weight, first moment and recurrence relationship
function adaptiveie(w,μ0,α,β,d)
  for logn = 4:20
    m=2^logn + 1
    pts=points(Segment(d[1],d[2]),m)
    pv=orthonormalpolynomialsvalues(μ0,α,β,pts)

    wv=w(pts)

    cfs=chebyshevtransform(pv[:,end] .^2 .* wv)

    if maximum(abs(cfs[end-8:end])) < 200*eps()
        ## chop to minimize oversample

        m=length(chop!(cfs,200*eps()))

        pts=points(Segment(d[1],d[2]),m)
        pv=orthonormalpolynomialsvalues(μ0,α,β,pts)
        wv=w(pts)

        return InvariantEnsemble(diagm(sqrt(wv))*pv,d)
    end
  end

  error("Did not converge")
end

#Reads in recurrence relationship and constructs OPs
function InvariantEnsemble(str::AbstractString,V::Function,d,n::Integer)
    file = joinpath(dirname(@__FILE__), "../data/CoefficientDatabase/" * str * "/" * string(n))
    μ0=readdlm(file * "norm.csv")[1]
    A=readdlm(file * "rc.csv",',')
    a=reshape((A[1,1:end-1]),size(A)[2]-1)
    b=reshape(A[2,2:end],size(A)[2]-1)

    adaptiveie(x->exp(-n.*V(x)),μ0,a,sqrt(b),d)
end


# For constructing InvariantEnsembles that do not scale with n
function InvariantEnsembleUnscaled(str::AbstractString,V::Function,d,n::Integer)
    file = joinpath(dirname(@__FILE__), "../data/CoefficientDatabaseUnscaled/" * str * "/")
    μ0=readdlm(file * "norm.csv")[1]
    A=readdlm(file * "rc.csv",',')
    a=reshape((A[1,1:end-1]),size(A)[2]-1)
    b=reshape(A[2,2:end],size(A)[2]-1)


    adaptiveie(x->exp(-V(x)),μ0,a[1:n-1],sqrt(b[1:n-1]),d)
end



#Decides whether to use built in recurrence or read it in
# Also contains ensemble data
function InvariantEnsemble(str::AbstractString,n::Integer)
    if(str == "GUE")
        InvariantEnsemble(str,x->x.^2,[-3.,3.],n)
    elseif(str == "Quartic")
        InvariantEnsemble(str,x->x.^4,[-3.,3.],n)
    elseif(str == "Eight")
        InvariantEnsemble(str,x->x.^8,[-3.,3.],n)
    elseif(str == "HODecay")
        InvariantEnsemble(str,x->x.^4/20 -4/15*x.^3 + x.^2/5+8/5*x,[-4.,6.],n)
    elseif(str == "Legendre")
        error("Legendre not implemented")
        #InvariantEnsembles(map(f->pad(f,2n),legendrepolynomials(n)))
    elseif(str == "CoshUnscaled")
        InvariantEnsembleUnscaled("Cosh",x->cosh(x),[-2acosh(1.0n),2acosh(1.0n)],n)
    elseif(str == "EightUnscaled")
        InvariantEnsembleUnscaled("Eight",x->x.^8,[-3.0n^(1/8),3.0n^(1/8)],n)
    elseif(str == "QuarticUnscaled")
        InvariantEnsembleUnscaled("Quartic",x->x.^4,[-3.0n^(1/4),3.0n^(1/4)],n)
    end
end


##  Construct the invariant ensemble kernel K_n(x,x)
iekernel(q::Array{Float64,2},d)=iekernel(q,d,plan_chebyshevtransform(q[:,1]))
function iekernel(q::Array{Float64,2},d,plan)
    n=size(q)[1]
    m=size(q)[2]
    ret=zeros(n)
    for i = 1:m
        for k = 1:n
            ret[k] += q[k,i].*q[k,i]
        end
    end

  Fun(d,plan*ret)
end


samplespectra(str::AbstractString,n::Integer,m::Integer) = samplespectra(InvariantEnsemble(str,n),m)


# Sample eigenvalues of invariant ensemble, m times
function RandomMatrices.eigvalrand(p::InvariantEnsemble,m::Integer)
    q = p.basis
    plan = plan_chebyshevtransform(q[:,1])
    pts=chebyshevpoints(size(q,1))

    hcat([samplespectra(p.basis,p.domain,plan,pts) for i=1:m]...)'
end

RandomMatrices.eigvalrand(p::InvariantEnsemble)=[eigvalrand(p,1)[1,:]...]


# Sample invariant ensemble
function Base.rand(p::InvariantEnsemble)
    Q=rand(Haar(2),size(p,1))
    Q*diagm(eigvalrand(p))*Q'
end

# Sample invariant ensemble, m times
function Base.rand(p::InvariantEnsemble,m::Integer)
    ei = eigvalrand(p,m)
    Array{Complex{Float64},2}[(    Q=rand(Haar(2),size(p,1));
    Q*diagm([ei[k,:]...])*Q') for k=1:size(ei,1)]
end



## samplespectra is back end for eigvalrand
samplespectra(q::Array{Float64,2},d)=samplespectra(q,d,plan_chebyshevtransform(q[:,1]),chebyshevpoints(size(q,1)))

function samplespectra(q::Array{Float64,2},d,plan,pts)
    n = size(q,2)
    r=Array(Float64,n)


    for k=1:n-1
        m = n - k + 1
        r[k] = samplecdf(normalizedcumsum!(iekernel(q,d,plan).coefficients/m))
        f=Float64[bary(q[:,j],pts,r[k]) for j = 1:m]
        r[k] = fromcanonical(d,r[k])
        Q=nullspace(f')
        q=q*Q
    end

    r[n] = sample(iekernel(q,d,plan))

    r
end


function iekernel(p::Array{F}) where {F<:Fun}
    ret = 0 .* p[1]
    for i = 1:length(p)
        ret += fasttimes(p[i],p[i])
    end

    ret
end

function samplespectra(p::Array{F}) where {F<:Fun}
    n = length(p)
    r=sample(iekernel(p)/n)
    if n==1
        r
    else
        f=map(q->q[r],p)
        Q=nullspace(f')'
        [r ;samplespectra(Q*p)]
    end
end

iekernel(p::InvariantEnsemble)=iekernel(p.basis,p.domain)
samplespectra(p::InvariantEnsemble)=samplespectra(p.basis,p.domain)





## Spectra database supports building up a database of precomputed samples


function spectradatabase(str,n::Integer,m::Colon)
  file = joinpath(dirname(@__FILE__), "../data/SpectraDatabase/" * str * "/" * string(n) * ".csv");

  if(filesize(file) == 0)
    []
  else
    readdlm(file,',')
  end
end

function spectradatabase(str,n::Integer,m::Integer)
    if  filesize(joinpath(dirname(@__FILE__), "../data/CoefficientDatabase")) == 0 ||
        filesize(joinpath(dirname(@__FILE__), "../data/CoefficientDatabaseUnscaled")) == 0
        error("No coefficient database found.  InvariantEnsembles package corrupted.")
    end


    dir = joinpath(dirname(@__FILE__), "../data/SpectraDatabase")
    if filesize(dir) == 0
        warn("Creating SpectraDatabase folder " * dir)

        mkdir(dir)
    end


    dir = joinpath(dirname(@__FILE__), "../data/SpectraDatabase/" * str)

    if filesize(dir) == 0
        warn("Creating folder " * dir)

        mkdir(dir)
    end

  file = dir * "/" * string(n) * ".csv"


  if filesize(file) == 0
    warn("No existing samples found.  Creating sample database.");

    ie = InvariantEnsemble(str,n)

    A=eigvalrand(ie,m+1)
    writedlm(file,A,',')
  else
    A = readdlm(file,',')

    if ( m > size(A)[1])
      warn(string(m) * " is greater than current samples " * string(size(A)[1]) * ".  Generating more samples");

      ie = InvariantEnsemble(str,n)

      An=eigvalrand(ie,m-size(A)[1]+1)

      A = vcat(A,An)

      writedlm(file,A,',')
    end
  end

    return A[1:m,:]
end



end #module
