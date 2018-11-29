pgfortran -llapack -lblas -Mcuda:cc60,cuda9.0,fastmath -Mvect=simd -Mlarge_arrays -Minline -ta=tesla:cc60 mod_precision.F90 mod_maths.F90 mod_polynomial.F90 mod_quadrature.F90 mod_bigquad.F90 mod_cudabigquad.F90 tests.F90 -O1 -o test #-Minfo


