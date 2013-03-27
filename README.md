RandomMatrices.jl
=================

Random matrix package for [Julia](http://julialang.org).

This extends the [Distributions](https://github.com/JuliaStats/Distributions.jl)
package to provide methods for working with matrix-valued random variables,
a.k.a. random matrices. State of the art methods for computing random matrix
samples and their associated distributions are provided.

## License
Copyright (c) 2013 [Jiahao Chen](https://github.com/jiahao) <jiahao@mit.edu>

Distributed under the [MIT License](http://opensource.org/licenses/MIT).

# Gaussian matrix ensembles

## Joint probability density functions (jpdfs)

Given eigenvalues `lambda` and the `beta` parameter of the random matrix distribution:

- `VandermondeDeterminant(lambda, beta)` computes the Vandermonde determinant
- `HermiteJPDF(lambda, beta)` computes the jpdf for the Hermite ensemble
- `LaguerreJPDF(lambda, n, beta)` computes the jpdf for the Laguerre(n) ensemble
- `JacobiJPDF(lambda, n1, n2, beta)` computes the jpdf for the Jacobi(n1, n2) ensemble

## Matrix samples

Constructs samples of random matrices corresponding to the classical Gaussian
Hermite, Laguerre(m) and Jacobi(m1, m2) ensembles.

- `GaussianHermiteMatrix(n, beta)`, `GaussianLaguerreMatrix(n, m, beta)`,
   `GaussianJacobiMatrix(n, m1, m2, beta)`
   each construct a sample dense `n`x`n` matrix for the corresponding matrix
   ensemble with `beta=1,2,4`

- `GaussianHermiteTridiagonalMatrix(n, beta)`,
  `GaussianLaguerreTridiagonalMatrix(n, m, beta)`,
  `GaussianJacobiSparseMatrix(n, m1, m2, beta)`
  each construct a sparse `n`x`n` matrix for the corresponding matrix ensemble
  for arbitrary positive finite `beta`. 
  `GaussianHermiteTridiagonalMatrix(n, Inf)` is also allowed.
  These sampled matrices have the same eigenvalues as above but are much faster
  to diagonalize oweing to their sparsity. They also extend Dyson's threefold
  way to arbitrary `beta`.
- `GaussianHermiteSamples(n, beta)`,
  `GaussianLaguerreSamples(n, m, beta)`,
  `GaussianJacobiSamples(n, m1, m2, beta)`
  return a set of `n` eigenvalues from the sparse random matrix samples

(Note the parameters of the Laguerre and Jacobi ensembles are not yet defined consistently.
For the first set of methods they are integers but for the rest they are reals.)

# Formal power series

Allows for manipulations of formal power series (fps) and formal Laurent series
(fLs), which come in handy for the computation of free cumulants.

## Types
- `FormalPowerSeries`: power series with coefficients allowed only for
  non-negative integer powers
- `FormalLaurentSeries`: power series with coefficients allowed for all
  integer powers

## FormalPowerSeries methods

### Elementary operations
- basic arithmetic operations `==`, `+`, `-`, `^`
- `*` computes the Cauchy product (discrete convolution)
- `.*` computes the Hadamard product (elementwise multiplication)
- `compose(P,Q)` computes the series composition P.Q
- `derivative` computes the series derivative
- `reciprocal` computes the series reciprocal

### Utility methods
- `trim(P)` removes extraneous zeroes in the internal representation of `P`
- `isalmostunit(P)` determines if `P` is an almost unit series
- `isconstant(P)` determines if `P` is a constant series
- `isnonunit(P)` determines if `P` is a non-unit series
- `isunit(P)` determines if `P` is a unit series
- `MatrixForm(P)` returns a matrix representation of `P` as an upper triangular
  Toeplitz matrix
- `tovector` returns the series coefficients

# Densities

Famous distributions in random matrix theory

- `Semicircle` provides the semicircle distribution
- `TracyWidom` computes the Tracy-Widom density distribution
  by brute-force integration of the Painlevé II equation

# Utility functions

- `hist_eig` computes the histogram of eigenvalues of a matrix using the
  method of Sturm sequences.
  This is recommended for `SymTridiagonal` matrices as it is significantly
  faster than `hist(eigvals())`
  This is also implemented for dense matrices, but it is pretty slow and
  not really practical.

# References
- James Albrecht, Cy Chan, and Alan Edelman,
    "Sturm Sequences and Random Eigenvalue Distributions",
    *Foundations of Computational Mathematics*,
    vol. 9 iss. 4 (2009), pp 461-483.
  [[pdf]](www-math.mit.edu/~edelman/homepage/papers/sturm.pdf)
  [[doi]](http://dx.doi.org/10.1007/s10208-008-9037-x)

- Ioana Dumitriu and Alan Edelman,
    "Matrix Models for Beta Ensembles",
    *Journal of Mathematical Physics*,
    vol. 43 no. 11 (2002), pp. 5830-5547
  [[doi]](http://dx.doi.org/doi: 10.1063/1.1507823)
  [arXiv:math-ph/0206043](http://arxiv.org/abs/math-ph/0206043)

- Alan Edelman, Per-Olof Persson and Brian D Sutton,
    "The fourfold way",
    *Journal of Mathematical Physics*,
    submitted (2013).
  [[pdf]](http://www-math.mit.edu/~edelman/homepage/papers/ffw.pdf)

- Alan Edelman and Brian D. Sutton,
    "The beta-Jacobi matrix model, the CS decomposition,
     and generalized singular value problems",
    *Foundations of Computational Mathematics*,
    vol. 8 iss. 2 (2008), pp 259-285.
  [[pdf]](http://www-math.mit.edu/~edelman/homepage/papers/betajacobi.pdf)
  [[doi]](http://dx.doi.org/10.1007/s10208-006-0215-9)

- Peter Henrici,
    *Applied and Computational Complex Analysis,
    Volume I: Power Series---Integration---Conformal Mapping---Location of Zeros*,
    Wiley-Interscience: New York, 1974
  [[worldcat]](http://www.worldcat.org/oclc/746035)

