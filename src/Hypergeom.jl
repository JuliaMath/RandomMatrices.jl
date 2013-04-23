#Lemma 3.1 - Update the Q
function hypergeom(m, alpha, a, b, x)

    #Compute P_{mn} using (4.2)
    n = length(x)
    H = m+1
    s = 1
    D = zeros(Pmn)
    J = zeros(Pmn, n) #J[N_kappa, i) = J^(alpha)_kappa (x_[1:i])
    J[1,:]=1 
    summation(1,1,m)
end

function summation(i, z, j)
    r = min(kappa[i-1], j)
    for kappa[i]=1:r
        if kappa[i]==1 && i>1
            N_kappa = D[N_kappa] = H
            H += r
        else
            N_kappa += 1
        end
        z *= T #T is rhs of (3.5)
        if kappa_prime[1]==1
            J[N_kappa, 1] = x[1]*(1+alpha*(kappa[1]-1))*J[N_kappa-1,1]
        end
        for t=2:n
            jack(0,1,0,t) #Computes J^(alpha)_kappa(x_[1:t])
        end
        s += z*J[N_kappa, 1]
        if j>kappa[i] && i<n summation(i+1,z,j-kappa[i]) end
    end
end

function jack(k, beta_kappamu, c, t)
    for i=k:mu_prime[1]
        if (k>0) && (mu[i]>mu[i+1])
            d=N_mu #Compute rhs using 4.1
            mu[i] -= 1
            #Update beta_kappamu with (3.10)
            if mu[i]>0
                jack(i, beta_kappamu, c+1, t)
            else
                J[N_kappa, t] += beta_kappamu * J[N_mu, t-1] * x[t]^(c+1)
            end
            mu[i] += 1
            N_mu = d
        end
    end

    if k==0
        J[N_kappa, t] += J[N_kappa, t-1]
    else
        J[N_kappa, t] += beta_kappamu * x[t]^c
    end
end

#    c = -(i-1)/alpha + kappa[i] - 1
#    d = kappa[i]*alpha - i
#    e[j] = d - j*alpha + kappaprime[j]
#    g[j] = e[j]+1
#    f[j] = kappa[j]*alpha - j - d
#    h[j] = f[j] + alpha
#    l[j] = h[j]*f[j]
#    
#    proda = prod([a[j]+c for j=1:p])
#    prodb = prod([b[j]+c for j=1:q])
#    prodg = prod([(g[j]-alpha)*e[j]/(g[j]*(e[j]+alpha)) for j=1:kappa[i]-1])
#    prodl = prod([(l[j]-f[j])/(l[j]+h[j]) for j=1:i-1])
#    
#    Q_kappa = Q_kappa[i] * proda/prodb * prodg * prodl
#    
#    #Lemma 3.2 - Update the Jack function
#    
#    #SPecial case of jack update, matrix = x I, prop. to identity 
#    J = J_kappa[i] * x * (n-i+1+alpha*(kappa[i]-1))
