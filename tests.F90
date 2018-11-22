program main
   use precision
   use bigquad
   implicit none

   integer(kind=int1), parameter :: n = 1000000

   integer(kind=int1) :: i

   real(kind=real2) :: sumx,sumw
   real(kind=real2) :: x(n),w(n)


   print *,'TESTING BIGQ, N = ',n
   call gauss_Legendre(n,x,w)
   print *,'COMPLETE'
   
   sumx = 0d0; sumw = 0d0
   do i=1,n
      sumx = sumx + x(i)
      sumw = sumw + w(i)
   enddo

   print *,'NODE SUM ERROR = ',sumx
   print *,'WEIGHT SUM ERROR = ',sumw-2d0

   
   return
end program main
