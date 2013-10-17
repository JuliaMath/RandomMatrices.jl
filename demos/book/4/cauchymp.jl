#cauchymp.jl
#Algorithm 4.1 of Random Eigenvalues by Alan Edelman

#Experiment:    Calculate Cauchy transform as trace of resolvent of Wishart matrix
#Theory:        Marcenko-Pastur law

m = 2000            # larger dimension of matrix
n = 1000            # smaller dimension of matrix
r = n/m             # aspect ratio
a = (1-sqrt(r))^2
b = (1+sqrt(r))^2
z = 3               # should be outside [a,b]
X = randn(m,n)
W = X'*X/m          # Generate Wishart matrix
println(( trace(inv(z*eye(n)-W))/n,
          (z-1+r-sqrt((z-a)*(z-b)))/(2*r*z) ))
