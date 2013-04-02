################################################
# Hypergeometric Function of a Matrix Argument #
################################################
#
# Translated from v1.3 of the MATLAB code by Plamen Koev, SJSU
#
# Reference:
# Plamen Koev and Alan Edelman, The Efficient Evaluation of the 
# Hypergeometric Function of a Matrix Argument, Mathematics of Computation
# 75 (2006), 833-846.
#
# MEX function [s,coef]=mhg([MAX,K],alpha,p,q,x,y)
# computes the truncated hypergeometric function pFq ^alpha(pqxy)
# of one or two matrix arguments.
# 
# The sum is only over those partitions kappa such that 
# a) |kappa|<=MAX 
# b) kappa has no more than n=length(x) parts
# c) kappa_1<=K (only if K is specified, if K is omitted, this restriction 
#    is not applied)
# 
# p and q are arrays, so mhg(30,9,[3 4],[5 6 7],[0.5 0.6],[0.8,0.9]) is 
# 2F3^9([3 4],[5 6 7][0.5 0.6], [0.8 0.9]) summed over all kappa with 
# |kappa|<=30
#
# K and y may be omitted. If y is omitted, the hypergeometric function of one 
# matrix argument is computed.

mhg{T1<:Real, T2<:Real, T3<:Real}(MAX::Integer, alpha::Real, p::Vector{T1}, q::Vector{T2}, x::Vector{T3})=mhg(MAX, MAX, alpha, p, q, x, [])
mhg{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(MAX::Integer, alpha::Real, p::Vector{T1}, q::Vector{T2}, x::Vector{T3}, y::Vector{T4})=mhg(MAX, MAX, alpha, p, q, x, y)
mhg{T1<:Real, T2<:Real, T3<:Real}(MAX::Integer, k::Integer, alpha::Real, p::Vector{T1}, q::Vector{T2}, x::Vector{T3})=mhg(MAX, K, alpha, p, q, x, [])
function mhg{T1<:Real, T2<:Real, T3<:Real, T4<:Real}(MAX::Integer, k::Integer, alpha::Real, p::Vector{T1}, q::Vector{T2}, x::Vector{T3}, y::Vector{T4})
    n, np, nq=length(x), length(p), length(q)
    
    coef=zeros(MAX+1) # set to zero, these are the coefficients of the polynomial  
    coef[1]=1 # free term equals one 
    
    # figure out the number of partitions |kappa|<= MAX with at most n parts 
    f=Int64[0:MAX+1]
    for i=2:n-1
        for j=i+1:MAX+1
            f[j+1]+=f[j-i+1]
        end
    end
    w=f[MAX+2]
    
    D     = zeros(Int64, w+1)
    Sx    = zeros(n*(w+1))
    xn    = zeros((n+1)*(MAX+2))
    prodx = zeros(n+1)

    prodx[2]=x[1]
    for i=2:n prodx[i+1]=prodx[i]*x[i] end
    for i=1:n
        Sx[n+i]=1
        xn[(MAX+2)*i+2]=1
        for j=2:MAX+1
            xn[(MAX+2)*i+j+1]=xn[(MAX+2)*i+j]*x[i]
        end
    end

    if y != []
       Sy    = zeros(n*(w+1))
       yn    = zeros((n+1)*(MAX+2))
       prody = zeros(n+1)
       prody[2]=y[1]
       for i=2:n prody[i+1]=prody[i]*y[i] end
       for i=1:n
           Sy[n+i]=1
           yn[(MAX+2)*i+2]=1
           for j=2:MAX+1
               yn[(MAX+2)*i+j+1]=yn[(MAX+2)*i+j]*y[i]
           end
       end
    end
    
    l = zeros(Int64, n+1)
    l[1]=k
    
    # this is what limits l[1] by the second element of MAX if needed and 
    # allows for the check l[i]<l[i-1] to be OK even for i=1 
    
    z     = ones(n+1)
    mu    = zeros(Int64, n+1)
    kt    = Float64[0:-1:-n]
    ww    = ones(Int64, n+1)
    d     = zeros(Int64, n+1)
    g     = zeros(Int64, n+1)
    mt    = zeros(n+1)
    blm   = zeros(n+1)
    lmd   = zeros(Int64, n+1)
    
    heap  = MAX+2
    cc, h, sl=1, 1, 1 # sl= sum(l) 
    while (h>0) 
        if ((l[h+1]<l[h]) && (MAX>=sl) && (z[h+1]!=0)) 
            l[h+1]+=1
            if ((l[h+1]==1) && (h>1) && (h<n)) 
                D[ww[h+1]]=heap
                ww[h+1]=heap
                k=MAX-sl+l[h+1]
                if (k>l[h]) k=l[h] end
                heap+=k
            else
                ww[h+1]+=1
            end
            w=ww[h+1]
                
            # Update Q 
            c=(1-h)/alpha+l[h+1]-1
            zn=alpha
            dn=kt[h+1]+h+1
            for j=1:np  zn*=p[j]+c end
            for j=1:nq  dn*=q[j]+c end
            if y!=[]
                zn*=alpha*l[h+1]
                dn*=n+alpha*c
                for j=2:h
                    t=kt[j]-kt[h+1]
                    zn*=t
                    dn*=t-1
                end
                zn/=dn
                dn=1 # trying to prevent overflow 
            end
            kt[h+1]+=alpha
            for j=2:h
                t=kt[j]-kt[h+1]
                zn*=t
                dn*=t+1
            end
            z[h+1]*=zn/dn
            # Working hard only when l has less than n parts 
            if (h<n) 
                t=h+1-alpha
                cc=1
                for j=2:h+1 cc*=(t+kt[j])/(h+kt[j]) end
                # computing the index of l-ones(1,h) 
                nmu=l[2]
                k=2
                while ((k<=h)&&(l[k+1]>1))
                    k+=1
                    nmu=D[nmu+1]+l[k]-2
                end
                Sx[w*n+h]=cc*prodx[h+1]*Sx[nmu*n+h]

                if (y!=[]) Sy[w*n+h]=cc*prody[h+1]*Sy[nmu*n+h] end
                cc=1 # this way we can just update from 1 in the h=n case

                d[h]-=1# technically this has to execute only when h>1 
                         # but is OK if it is always executed d[0] will 
                         # end up being -MAX at the end of the code 
                d[h+1]=l[h+1]# for (k=1k<hk+=1) d[k]=l[k]-l[k+1] 
                         # this happens automatically now via updates 
                
                lg=0
                for k=2:h+1
                    if (d[k]>0)
                        lg+=1
                        g[lg+1]=k-1
                    end
                end
                slm=1 # this is sum(l-mu) 
                nhstrip=1
                for k=2:lg+1
                    nhstrip*=d[g[k]+1]+1
                    nhstrip-=1
                end

                mu[2:h+1]=l[2:h+1]
                mt[2:h+1]=kt[2:h+1]
                for k=2:lg+1
                    blm[k]=1
                    lmd[k]=l[g[k]+1]-d[g[k]+1]
                end
                
                for i=1:nhstrip
                    j=lg
                    gz=g[lg+1]
                    while (mu[gz+1]==lmd[j+1]) 
                        mu[gz+1]=l[gz+1]
                        mt[gz+1]=kt[gz+1]
                        slm-=d[gz+1]
                        gz=g[j]
                        j-=1
                    end
                    t=kt[gz+1]-mt[gz+1]
                    
                    zn=1+t
                    dn=t+alpha
                    for k=2:gz
                        q1=mt[k]-mt[gz+1]
                        q2=kt[k]-mt[gz+1]
                        zn*=(alpha-1+q1)*(1+q2)
                        dn*=q1*(alpha+q2) 
                    end
                    blm[j+1]*=zn/dn
                    mu[gz+1]-=1
                    mt[gz+1]-=alpha
                    slm+=1
                    for k=j+1:lg blm[k+1]=blm[j+1] end 
                    
                    # next, find the index of mu 
                    nmu=mu[2]+1
                    for k=2:h-((mu[h+1]==0)?1:0)
                        nmu=D[nmu+1]+mu[k+1]-1
                    end
                    for k=h+1:n
                        Sx[w*n+k]+=blm[j+1]*Sx[nmu*n+k-1]*xn[k*(MAX+2)+slm+1]
                    end
                    if y!=[]
                        for k=h+1:n
                            Sy[w*n+k]+=blm[j+1]*Sy[nmu*n+k-1]*yn[k*(MAX+2)+slm+1]
                        end
                    end
                end
                
                for k=h:n-1 Sx[w*n+k+1]+=Sx[w*n+k] end
                if y!=[]
                    for k=h:n-1 Sy[w*n+k+1]+=Sy[w*n+k] end
                    coef[sl+1]+=z[h+1]*Sx[w*n+n]*Sy[w*n+n]
                else
                    coef[sl+1]+=z[h+1]*Sx[w*n+n]
                end
            else #h>=n 
                # computing the index of the partition l-l[n]*ones(1,n) 
                nmu=l[2]-l[n+1]+1
                k=2
                while ((k<n)&&(l[k+1]>l[n+1]))
                    k+=1
                    nmu=D[nmu+1]+l[k]-1-l[n+1]
                end
                # cc is 1 if l[n+1]==1, (guaranteed by the h<n case) 
                #   we then update from the previous  
              
                if y!=[]
                    t=(1/alpha+l[n+1]-1)/l[n+1]
                    for k=1:n-1 t*=(1+kt[k+1]-kt[n+1])/(alpha+kt[k+1]-kt[n+1]) end
                    cc*=t*t*prodx[n+1]*prody[n+1]
                    coef[sl+1]+=z[n+1]*cc*Sx[nmu*n+n]*Sy[nmu*n+n]
                else 
                    cc*=(1/alpha+l[n+1]-1)*prodx[n+1]/l[n+1]
                    for k=2:n cc*=(1+kt[k]-kt[n+1])/(alpha+kt[k]-kt[n+1]) end
                    coef[sl+1]+=z[n+1]*cc*Sx[nmu*n+n]
                  @printf("%d %f\n", sl+1, coef[sl+1])
	          @printf("%d %f %f %d %d\n", n, z[n+1], cc, nmu*n+n, Sx[nmu*n+n])
                 end
            end
            if (h<n) 
                z[h+2]=z[h+1]
                h+=1
                ww[h+1]=w
            end
            sl+=1
        else  
            sl-=l[h+1]
            l[h+1]=0
            kt[h+1]=-h
            h-=1
        end
    end # of while h>0 
    s=sum(coef)
    s, coef
end

#println(mhg(30,9,[3, 4],[5, 6, 7],[0.5, 0.6],[0.8,0.9]))
#println(1.054910425171565)
println(mhg(20,2,[0.1,0.2],[0.3,0.4],[0.5,0.6,0.7]))
println(3.134912003261897)
println(mhg(20,1,[],[0.1],[0.2,0.3,0.4],[0.5,0.6,0.7]))
println(6.626503739912625)
println(mhg(10,2.0,[],[],[1 0;0 1]))
