export hist_eig

#Uses the method of Sturm sequences to compute the histogram of eigenvalues of a matrix
#
# Reference
#   Cy P. Chan, "Sturm Sequences and the Eigenvalue Distribution of
#   the Beta-Hermite Random Matrix Ensemble", M.Sc. thesis (MIT), 2007
#
# For better numerical stability, this algorithm is modified so that only the sign
# of the determinant is stored.
#
#For general matrices, this is slower than hist(eigvals(M), bins) and is NOT recommended
function hist_eig(M::AbstractMatrix, bins::Vector{GridPoint}) where {GridPoint <: Number}
  n = size(M)[1]
  NumBins = length(bins)
  histogram = zeros(NumBins)
  SturmSequence = zeros(n)
  for BinId=1:NumBins
    K = M - bins[BinId]*eye(n)
    #Progression of determinants of lower-right submatrices
    #SturmSequence = [det(K[end-j+1:end,end-j+1:end]) for j=1:n]
    #SturmRatioSequence = [SturmSequence[j+1]/SturmSequence[j] for j=1:n]

    SturmSequenceSigns = [int(sign(det(K[end-j+1:end,end-j+1:end]))) for j=1:n]
    SturmSequenceSigns = [1; SturmSequenceSigns]
    SturmRatioSequence = [SturmSequenceSigns[j+1]/SturmSequenceSigns[j] for j=1:n]
    histogram[BinId] = sum([r < 0 for r in SturmRatioSequence])
  end
  histogram = int(diff([histogram; n]))
end



#Uses the method of Sturm sequences to compute the histogram of eigenvalues of
#a symmetric tridiagonal mmatrix
#
# Reference
#   Cy P. Chan, "Sturm Sequences and the Eigenvalue Distribution of
#   the Beta-Hermite Random Matrix Ensemble", M.Sc. thesis (MIT), 2007
#
function hist_eig(M::SymTridiagonal, bins::Vector{GridPoint}) where {GridPoint <: Number}
  n = size(M)[1]
  NumBins = length(bins)
  histogram = zeros(NumBins)
  r = zeros(n)
  for BinId=1:NumBins
    a, b = M.dv - bins[BinId], M.ev
    #Formula (2.3)
    r[1] = a[1]
    for i=2:n
      r[i] = a[i] - b[i-1]^2/r[i-1]
    end
    histogram[BinId] = sum([SturmRatio < 0 for SturmRatio in r])
  end
  histogram = int(diff([histogram; n]))
end

