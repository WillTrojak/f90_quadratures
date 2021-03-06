# f90_quadratures
Modern Fortran modules to generate quadrature rules.

# Functionality

* Gauss-Legendre - this is via Golub-Welsh and Townsend-Hale. Golub-Welsch ok upto n ~ 10^3 [O(n^2) scaling] and Townsend-Hale ok for 10^3 < n <~ 10^9 [O(n) scaling]. Use cudabigquad for n ~> 10^6 as much faster.

* Gauss-Jacobi - this is via Golub-Welsch with O(n^2) complexity which is quick upto n ~ 10^3

* Gauss-Christofel - currently only polynomial weight functions are supported scaling is currently O(n^2), runtime is ~double that of Gauss-Legendre. Uses LAPACK Divide-and-conquer for n > 500, giving O(nlog(n)) scaling. 

* Clenshaw-Curtis - complexity O(n^2) but slower than Gauss-Legendre 

* Gauss-Kronrod - this is through Golub-Welsch and is agian O(n^2)

* Others - Chebyshev,Gegenbauer, Gen' Laguerre, Gen' Hermite, Exponential and Rational (beta distribution weights) can also be calculated. These are accessable through the subroutine `quad`. (this is using code by S. Elhay, J. Kautsky and J. Burkardt which uses Golub-Welsch under the hood).

# Requirements for `tests.F90`
* pgfortran
* CUDA 9.0
* LAPACK and BLAS

# TODO
* Extend Large N Gauss-Legendre quad to Gauss-Jacobi

