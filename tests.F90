program main
   use precision
   use cudabigquad, only : bigLc => gauss_Legendre
   use bigquad, only : bigL => gauss_Legendre
   use quadrature,only : smallL => gauss_Legendre
   implicit none

   integer(kind=int1), parameter :: n1 = 100000000, n2 = 5

   integer(kind=int1) :: i

   real(kind=real2) :: time0,time1
   real(kind=real2) :: sumx,sumw
   real(kind=real2) :: x(n1),w(n1),t(n1)

   real(kind=real2) :: xs(n2),ws(n2)

   print *,'TESTING CUDA BIGQ, N = ',n1
   call cpu_time(time0)
   call bigLc(n1,x,w)
   call cpu_time(time1)
   print *,'runtime: ',time1 - time0
   print *,'COMPLETE'
   
   sumx = 0d0; sumw = 0d0
   do i=1,n1
      sumx = sumx + x(i)
      sumw = sumw + w(i)
   enddo

!!$   print *,'NODE SUM ERROR = ',sumx
!!$   print *,'WEIGHT SUM ERROR = ',sumw-2d0
!!$   
!!$   print *,'TESTING BIGQ, N = ',n1
!!$   call cpu_time(time0)
!!$   call bigL(n1,x,w)
!!$   call cpu_time(time1)
!!$   print *,'runtime: ',time1 - time0
!!$   print *,'COMPLETE'
!!$   
!!$   sumx = 0d0; sumw = 0d0
!!$   do i=1,n1
!!$      sumx = sumx + x(i)
!!$      sumw = sumw + w(i)
!!$   enddo

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
