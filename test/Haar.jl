if isdefined(:_HAVE_GSL)
N=5
A=randn(N,N)
B=randn(N,N)
Q=UniformHaar(2, N)

#Test case where there are no symbolic Haar matrices
@test_approx_eq eval(expectation(N, :Q, :(A*B))) A*B
#Test case where there is one pair of symbolic Haar matrices
@test_approx_eq eval(expectedtrace(N, :Q, :(A*Q*B*Q'))) trace(A)*trace(B)/N 

println("Case 3")
println("E(A*Q*B*Q'*A*Q*B*Q') = ", eval(expectation(N, :Q, :(A*Q*B*Q'*A*Q*B*Q'))))

end #_HAVE_GSL
