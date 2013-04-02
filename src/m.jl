using RandomMatrices

n=100
m=300
v = eigvals(GaussianLaguerreMatrix(n,m,1))
println(length(v))

l=n/m
vv = GaussianLaguerreSamples(n, m, 1)
#MarcenkoPastur says the support is in

println(min(v),'\t',(1-sqrt(l))^2)
println(max(v),'\t',(1+sqrt(l))^2)

#How many eigenvalues at 0?

