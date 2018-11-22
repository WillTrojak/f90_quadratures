# f90_quadratures
Modern Fortran modules to generate quadrature rules 

# Functionality

* Gauss-Legendre - this is via Golub-Welsh and Townsend-Hale. Golub-Welsch ok upto n ~ 10^3 [O(n^2) scaling] and Townsend-Hale ok for 10^3 < n <~ 10^9 [O(n) scaling] 

* Gauss-Jacobi - this is via Golub-Welsch with O(n^2) complexity which is quick upto n ~ 10^3

* Gauss-Christofel - currently only polynomial weight functions are supported scaling is currently O(n^3) therefore ok for n < ~500. This uses a niave approach but is suitable ofor most applications

* Clenshaw-Curtis - complexity O(n^2) but slower than Gauss-Legendre 

* Gauss-Kronrod - this is through Golub-Welsch and is agian O(n^2)

* Others - Chebyshev,Gegenbauer, Gen' Laguerre, Gen' Hermite, Exponential and Rational (beta distribution weights) can also be calculated. These are accessable through the subroutine `quad`. (this is using code by S. Elhay, J. Kautsky and J. Burkardt which uses Golub-Welsch under the hood).

# TODO
* Divide and Conquer method for symetric-tridiagonal eigenvaules to get O(nlog(n)) Golub-Welsch.
* Extend Large N Gauss-Legendre quad to Gauss-Jacobi