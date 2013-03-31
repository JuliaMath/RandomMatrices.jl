
#Samples the largest eigenvalue of the GUE matrix
function SampleTracyWidom(n::Integer)
    k=int(n-10*n^(1/3)-1)
    A=[chi(i) for i=(n-1):-1:k]
    B=randn(n-k+1)
    v=max(eigvals(SymTridiagonal(B, A)))
end

