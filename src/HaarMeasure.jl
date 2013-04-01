
# Computes samplex of real or complex Haar matrices of size nxn
function HaarMatrix(n::Integer, beta::Integer)
    if beta==1
        M=randn(n,n)
    elseif beta==2
        M=randn(n,n)+im*randn(n,n)
    else
        error(string("beta = ",beta, " not implemented."))
    end
    qr(M)[1]
end



