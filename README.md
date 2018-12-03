# f90_quadratures
Modern Fortran modules to generate quadrature rules 

# Functionality

* Gauss-Legendre - this is via Golub-Welsh and Townsend-Hale. Golub-Welsch ok upto n ~ 10^3 [O(n^2) scaling] and Townsend-Hale ok for 10^3 < n <~ 10^9 [O(n) scaling]. Use cudabigquad for n ~> 10^6 as much faster.

* Gauss-Jacobi - this is via Golub-Welsch with O(n^2) complexity which is quick upto n ~ 10^3

* Gauss-Christofel - currently only polynomial weight functions are supported scaling is currently O(n^2), runtime is ~double that of Gauss-Legendre.

* Clenshaw-Curtis - complexity O(n^2) but slower than Gauss-Legendre 

* Gauss-Kronrod - this is through Golub-Welsch and is agian O(n^2)

* Others - Chebyshev,Gegenbauer, Gen' Laguerre, Gen' Hermite, Exponential and Rational (beta distribution weights) can also be calculated. These are accessable through the subroutine `quad`. (this is using code by S. Elhay, J. Kautsky and J. Burkardt which uses Golub-Welsch under the hood).

# TODO
* Test LAPACK Divide and Conquer method for symetric-tridiagonal eigenvaules to get O(nlog(n)) Golub-Welsch.
* Extend Large N Gauss-Legendre quad to Gauss-Jacobi