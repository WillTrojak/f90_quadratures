program main
   use precision
   use bigquad
   implicit none

   integer(kind=int1), parameter :: n = 1000000

   integer(kind=int1) :: i

   real(kind=real2) :: sumx
   real(kind=real2) :: x(n),t(n)


   print *,'TESTING BIGQ, N = ',n
   call gaussl_nodes(n,x,t)
   print *,'COMPLETE'
   
   sumx = 0d0
   do i=1,n
      sumx = sumx + x(i)
   enddo

   print *,'NODE SUM ERROR = ',sumx

   
   return
end program main
