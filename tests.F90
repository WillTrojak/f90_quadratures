program main
   use precision
   use bigquad, only : bigL => gauss_Legendre
   use quadrature,only : smallL => gauss_Legendre
   implicit none

   integer(kind=int1), parameter :: n1 = 1000000, n2 = 5

   integer(kind=int1) :: i

   real(kind=real2) :: sumx,sumw
   real(kind=real2) :: x(n1),w(n1)

   real(kind=real2) :: xs(n2),ws(n2)

   print *,'TESTING BIGQ, N = ',n1
   call bigL(n1,x,w)
   print *,'COMPLETE'
   
   sumx = 0d0; sumw = 0d0
   do i=1,n1
      sumx = sumx + x(i)
      sumw = sumw + w(i)
   enddo

   print *,'NODE SUM ERROR = ',sumx
   print *,'WEIGHT SUM ERROR = ',sumw-2d0

   print *
   print *,'TESTING SMALLQ, N = ',n2
   call smallL(n2,xs,ws)
   print *,'COMPLETE'
   print *
   
   do i=1,n2
      print *,'xs',i,xs(i)
   enddo
   print *
   
   do i=1,n2
      print *,'ws',i,ws(i)
   enddo
   
   
   return
end program main
