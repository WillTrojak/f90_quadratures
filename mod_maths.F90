module maths  
   use precision
   implicit none
   
   real(kind=real2), parameter :: pi = 4d0*atan(1d0)
   real(kind=real2), parameter :: ex = exp(1d0)

   private
      
   public :: choose, determinant, enorm, ex, frobeniusNorm, factorial, factorial2, &
             facBig, factorialArrayReal, factorialReal, genHypGeo, heapsort, &
             integrate, invert, inverse, invert_diag, isnan, isinf, linspace, &
             Lpnorm_vector, new_random_seed, pi, pochhammerR, pseudo_inverse, &
             pseudo_inverse_R, random_int, random_int_range, sortrows, &
             sortrows_real
   
   interface heapsort
      module procedure heapsort_1array
      module procedure heapsort_2array
   end interface heapsort

contains
   !**********************************************************************
   function choose(n,k) result(c)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n,k

      integer(kind=int1) :: c

      c = nint(gamma(n+1d0)/gamma(k+1d0)/gamma(n-k+1d0))

      return
   end function choose
   !**********************************************************************
   !> @brief calculates the determinant or perminant of a matrix
   !> @par to calculate the determinant p=-1, for the perminant p=1, this is
   !! a very expensive method with complexity O(n!) so only use for small
   !! matrices. A better method would use LU decompostion. 
   !> @param[in] a matrix
   !> @param[in] n size
   !> @param[in] p det/per selector
   !> @param[out] sum the cumulating sum giving the det/per
   !**********************************************************************
   recursive function determinant(a,n,p) result(sum)
      ! p =  1 => permanent
      ! p = -1 => determinant
      use precision
      implicit none
      
      integer(kind=int1), intent(in) :: n,p
      real(kind=real2), intent(in) :: a(n,n)
      
      real(kind=real2) :: sum
      
      integer(kind=int1) :: i, sgn
      real(kind=real2) :: b(n-1,n-1)
      
      if(n .eq. 1) then
         sum = a(1,1)
      else
         sum = 0
         sgn = 1
         
         do i=1,n
            b(:,:(i-1)) = a(2:,:i-1)
            b(:,i:) = a(2:,i+1:)
            
            sum = sum + sgn*a(1, i)*determinant(b,n-1,p)
            sgn = sgn*p
         enddo
      endif

      return
   end function determinant
   !**********************************************************************
   !> @brief calculates the euclidian or l_2 norm
   !> @par as this norm is more commonly used, added a more efficient implementation
   !> @param[in] x vector
   !> @param[out] n l_2 norm
   !**********************************************************************
   function enorm(x) result(n)
      use precision
      implicit none

      real(kind=real2), intent(in) :: x(:)

      real(kind=real2) :: n
      
      n = sum(x(:)*x(:))
      n = sqrt(n)

      return
   end function enorm
   !**********************************************************************
   function facBig(n) result(fac)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n
      
      integer*8 :: fac,i

      fac = 1
      if(n .gt. 0)then
         do i=1,n
            fac = fac*i
         enddo
      endif

      return
   end function facBig
   !**********************************************************************
   pure recursive double precision function factorial(n) result(fac)  
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n

      if(n .eq. 0) then
         fac = 1d0
      else
         fac = dble(n)*factorial(n-1)
      endif

      return
   end function factorial
   !**********************************************************************
   pure recursive double precision function factorial2(n) result(fac)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n
      
      if(n .ge. 2)then
         fac = dble(n)*factorial2(n-2)
      else
         fac = 1d0
      endif
     
     return
   end function factorial2
   !**********************************************************************
   function factorialArrayReal(n) result(facA)
      use precision
      implicit none

      integer(kind=int1) :: n

      real(kind=real2) :: facA(n+1)

      integer(kind=int1) :: i

      facA(1) = 1
      do i=1,n
         facA(i+1) = facA(i)*dble(i)
      enddo

      return
   end function factorialArrayReal
   !**********************************************************************
   pure function factorialReal(n) result(fac)  
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n

      real(kind=real2) :: fac

      integer(kind=int1) :: in
      
      fac = 1d0
      if(n .eq. 0) then
         continue
      else
         do in=1,n
            fac = fac*dble(in)
         enddo
      endif

      return
   end function factorialReal
   !**********************************************************************
   !> @breif calculates the frobenius norm of a matrix
   !> @param[in] a matrix
   !> @param[out] norm the norm
   !**********************************************************************
   function frobeniusNorm(a) result(norm)
      use precision
      implicit none

      real(kind=real2), intent(in) :: a(:,:)

      real(kind=real2) :: norm

      integer(kind=int1) :: n,m,in,im
      
      n = size(a,1); m = size(a,2)
      
      norm = 0d0

      do im=1,m
         do in=1,n
            norm = norm + abs(a(in,im))**2
         enddo
      enddo

      norm = sqrt(norm)
      
      return
   end function frobeniusNorm
   !**********************************************************************
   function genHypGeo(p,q,a,b,z) result(f)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: p,q

      real(kind=real2), intent(in) :: a(p),b(q),z

      real(kind=real2) :: f

      integer(kind=int1), parameter :: nmax = 50
      integer(kind=int1) :: ip,iq,n

      real(kind=real2), parameter :: tol = 1d-8
      real(kind=real2) :: delf,an,bn,fac(nmax+1)

      logical :: conv
      
      conv = .false.
      
      fac = factorialArrayReal(nmax)

      f = 0d0
      do n=0,nmax
         an = 1d0
         do ip=1,p
            an = an*pochhammerR(a(ip),n)
         enddo
         bn = 1d0
         do iq=1,q
            bn = bn*pochhammerR(b(iq),n)
         enddo
         delf = (an/bn)*z**n/fac(n+1)
         f = f + delf
         
         if(abs(delf) .lt. abs(f)*tol) exit
      enddo

      return
   end function genHypGeo
   !**********************************************************************
   subroutine heapsort_1array(a)
      use precision
      implicit none
      
      real(kind=real2), intent(inout) :: a(0:)

      integer(kind=int1) :: start,n,bottom
      
      n = size(a)
      do start=(n - 2)/2,0,-1
         call siftdown(start,n,a)
      enddo
      
      do bottom=n-1,1,-1
         call swap_entries(a(0),a(bottom))
         call siftdown(0,bottom,a)
      enddo

      return
   end subroutine heapsort_1array
   !**********************************************************************
   subroutine heapsort_2array(x,a)
      use precision
      implicit none
      
      real(kind=real2), intent(inout) :: x(0:)
      real(kind=real2), intent(inout) :: a(0:)

      integer(kind=int1) :: start,n,bottom
      
      n = size(a)
      do start=(n - 2)/2,0,-1
         call siftdown_array(start,n,x,a)
      enddo
      
      do bottom=n-1,1,-1
         call swap_entries(a(0),a(bottom))
         call swap_entries(x(0),x(bottom))
         call siftdown_array(0,bottom,x,a)
      enddo

      return
   end subroutine heapsort_2array
   !**********************************************************************
   function integrate(f1,g1,w) result(v)
      use precision
      implicit none

      real(kind=real2), intent(in) :: f1(:,:), g1(:,:), w(:)
    
      real(kind=real2) :: v(size(f1,1),size(g1,1))

      integer(kind=int1) :: in, jn, kn

      ! Loop over the points in v
      do in=1,size(f1,1)
      do jn=1,size(g1,1)

         ! Accumulate the gaussian integrator
         v(in,jn) = 0d0
         do kn=1,size(w,1)
            v(in,jn) = v(in,jn) + f1(in,kn)*g1(jn,kn)*w(kn)
         enddo
      enddo
      enddo

      return
   end function integrate
   !**********************************************************************
   !> @brief matrix inversion
   !> @param[in] n matrix size
   !> @param[in] a matrix to invert
   !> @param[out] b inversion
   !> @param[out] safe optional error check
   !**********************************************************************
   subroutine invert(n,a,b)  
      use precision
      use logging, only : log
      implicit none
       
      integer(kind=int1), intent(in) :: n

      real(kind=real2), intent(inout) :: a(n,n),b(n,n)
      
      integer(kind=int1) :: ipiv(n),info
      real(kind=real2) :: c(n,n),work(n)
      
      b = a
      call dgetrf(n,n,a,n,ipiv,info)
      call dgetri(n,a,n,ipiv,work,n,info)
      c = b
      b = a
      a = c
         
      return
   end subroutine invert
   !*******************************************************************
   function invert_diag(n,a) result(ai)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n

      real(kind=real2), intent(in) :: a(n,n)

      real(kind=real2) :: ai(n,n)

      integer(kind=real2) :: i
      
      ai = 0d0
      do i=1,n
         ai(i,i) = 1d0/a(i,i)
      enddo
      
      return
   end function invert_diag
   !*******************************************************************
   !> @breif functional implementation of matrix inversion
   !> @param[in] a matrix to invert
   !> @param[out] ai inversion
   !*******************************************************************
   function inverse(a) result(ai)
      use precision
      implicit none

      real(kind=real2), intent(in) :: a(:,:)
      
      real(kind=real2) :: ai(size(a,1),size(a,2))

      integer(kind=int1) :: ipiv(size(a,1)), info, n
      real(kind=real2) :: c(size(a,1), size(a,2)),work(size(a,1))

      n = size(a,1)
      
      c = a
      call dgetrf(n,n,c,n,ipiv,info)
      call dgetri(n,c,n,ipiv,work,n,info)
      ai = c
      
      return
   end function inverse
   !**********************************************************************
   !> @breif calculates normed errro of inversion
   !> @param[in] a matrix
   !> @param[in] ai inversion
   !> @param[out] error inversion error
   !**********************************************************************
   function inversionError(a,ai) result(error)
      use precision
      implicit none

      real(kind=real2), intent(in) :: a(:,:),ai(:,:)

      real(kind=real2) :: error

      integer(kind=int1) :: n,in
      real(kind=real2), allocatable :: eye(:,:)

      n = size(a,1)
      allocate(eye(n,n))

      eye = 0
      do in=1,n
         eye(in,in) = 1d0
      enddo

      eye = eye - matmul(a,ai)
      error = frobeniusNorm(eye)
      
      return
   end function inversionError
   !**********************************************************************
   logical function isinf(y)
      use, intrinsic :: iso_fortran_env
      use precision
      implicit none

      class(*), intent(in) :: y
      
      real(kind=4) :: fpinf,fminf
      real(kind=8) :: dpinf,dminf
      
      fpinf = b'01111111100000000000000000000000'   ! +infinity
      fminf = b'11111111100000000000000000000000'   ! -infinity
      dpinf = b'0111111111110000000000000000000000000000000000000000000000000000' ! +infinity
      dminf = b'1111111111110000000000000000000000000000000000000000000000000000' ! -infinity
      
      isinf = .false.

      select type(y)
      type is (real(real32))
         if(y .eq. fpinf) isinf = .true.
         if(y .eq. fminf) isinf = .true.
      type is (real(real64))
         if(y .eq. dpinf) isinf = .true.
         if(y .eq. dminf) isinf = .true.
      end select
      
      return
   end function isinf
   !**********************************************************************
   logical function isnan(y)
      use, intrinsic :: iso_fortran_env
      use precision
      implicit none

      class(*), intent(in) :: y
      
      isnan = .false.

      select type(y)
      type is (real(real32))
         if((y .gt. 0e0) .eqv. (y .le. 0e0)) isnan = .true.
         if(y .ne. y) isnan = .true.
      type is (real(real64))
         if((y .gt. 0d0) .eqv. (y .le. 0d0)) isnan = .true.
         if(y .ne. y) isnan = .true.
      end select
      
      return
   end function isnan
   !**********************************************************************
   function linspace(n,x0,x1) result(x)
      implicit none

      integer*4, intent(in) :: n
      real*8, intent(in) :: x0,x1

      integer*4 :: i
      real*8 :: x(n),nd,dx

      nd = dble(n) - 1d0
      dx = (x1-x0)/nd

      do i=1,n
         x(i) = x0 + dx*(i-1d0)
      enddo

      return
   end function linspace
   !**********************************************************************
   !> @brief calculates the L_p norm of a vector
   !> @par calculates (x(1)^p + x(2)^p + ...)^(1/p)
   !> @param[in] x vector
   !> @param[in] p norm type
   !> @param[out] lp norm value
   !**********************************************************************
   function Lpnorm_vector(x,p) result(lp)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: p
      real(kind=real2), intent(in) :: x(:) 

      real(kind=real2) :: lp

      real(kind=real2), allocatable :: xp(:)

      allocate(xp,mold=x)
      xp(:) = x(:)**p
      lp = sum(xp)
      lp = lp**(1d0/dble(p))
      
      return
   end function Lpnorm_vector
   !**********************************************************************
   function monteCarloIntegration(npnts,stride,v,u) result(int)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: npnts,stride

      real(kind=real2), intent(in) :: v,u(npnts)

      real(kind=real2) :: int

      integer(kind=int1) :: ic,in
      real(kind=real2) :: temp

      int = 0d0
      do ic=1,npnts,stride
         temp = 0d0
         do in=1,stride
            temp = temp + u(ic-1+in)
         enddo
         int = int + temp
      enddo
      int = v*int/dble(npnts)
      
      return
   end function monteCarloIntegration
   !**********************************************************************
   recursive function multiindex(ns,ni,sum) result(m)
      ! ns - number of state, ni - number of indexes, sum = sum of multiindexs
      use precision
      implicit none

      integer(kind=int1), intent(in) :: ns,ni,sum

      integer(kind=int1) :: m(ni,ns)

      integer(kind=int1) :: ma,s,is
      integer(kind=int1), allocatable :: temp(:,:)

      if(ni .eq. 1)then
         m = sum
      else
         ma = 0
         do is=0,sum

            s = simplicy_number(sum+1-is,ni-2)
            
            allocate(temp(ni-1,s))
            temp = multiindex(s,ni-1,sum-is)
            
            m(1,ma+1:ma+s) = is
            m(1+1:ni,ma+1:ma+s) = temp(1:ni-1,1:s)
            ma = ma + s

            deallocate(temp)
         enddo
      endif
         
      return
   end function multiindex
   !**********************************************************************
   subroutine new_random_seed
      use precision
      implicit none
            
      integer(kind=int1) :: tvalues(8),n_seed
      integer(kind=int1), allocatable :: seed(:)

      ! Set new random number seed
      call date_and_time(values=tvalues)
      call random_seed(size=n_seed)
      allocate(seed(n_seed))
      seed(1) = tvalues(8)*tvalues(7)*tvalues(6)
      if(n_seed .gt. 2)then
         seed(2) = (tvalues(1)+tvalues(2))*tvalues(8)*tvalues(3)
         seed(n_seed) = (tvalues(3)+tvalues(4))*tvalues(7)*tvalues(5)
      endif
      call random_seed(put=seed)
      deallocate(seed)

      return
   end subroutine new_random_seed
   !**********************************************************************
   pure function pochhammerR(x,n) result(p)
      use precision
      implicit none
      
      integer(kind=int1), intent(in) :: n
      real(kind=real2), intent(in) :: x

      real(kind=real2) :: p

      integer(kind=int1) :: i

      p = 1d0
      if(n .eq. 0) return
      do i=0,n-1
         p = p*(x+dble(i))
      enddo

      return
   end function pochhammerR
   !**********************************************************************
   function pseudo_inverse_R(b) result(ai)
      use precision
      use logging, only : log
      implicit none
      
      real(kind=real2), intent(in) :: b(:,:)
      
      real(kind=real2) :: ai(size(b,2),size(b,1))

      integer(kind=int1) :: m,n,i,j,lwork,lwmax
      
      real(kind=real2) :: a(size(b,1),size(b,2)),s(size(b,2))
      real(kind=real2) :: sigma(size(b,1),size(b,2)),work(10000),sigmi(size(b,2),size(b,1))
      real(kind=real2) :: u(size(b,1),size(b,1)),ut(size(b,1),size(b,1))
      real(kind=real2) :: v(size(b,2),size(b,2)),vt(size(b,2),size(b,2))
      real(kind=real2), parameter :: eps = 1d2*epsilon(b(1,1)) 
      
      m = size(b,1)
      n = size(b,2)

      !if(n .gt. m) call log(event_type=5,event='PSEUDO INVERSE DIMENSION INVALID')
      
      a = b
      
      lwork = -1
      lwmax = 10000
      
      call dgesvd('A','A',m,n,a,m,s,u,m,vt,n,work,lwork,i)
      lwork = min(lwmax,int(work(1)))
      
      call dgesvd('A','A',m,n,a,m,s,u,m,vt,n,work,lwork,i)
      
      sigma = 0d0
      do i=1,n
         sigma(i,i) = s(i)
      enddo
      
      ! find pseudo inverse of sigma
      sigmi = transpose(sigma)
      do i=1,n
         if(sigmi(i,i) .lt. eps) then
            sigmi(i,i) = 0d0
         else
            sigmi(i,i) = 1d0/sigmi(i,i)
         endif
      enddo
      
      v = transpose(vt)
      ut = transpose(u)
      
      ai = matmul(v,matmul(sigmi,ut))
      
      return  
   end function pseudo_inverse_R
   !**********************************************************************
   function pseudo_inverse(b) result(ai)
      use precision
      use logging, only : log
      implicit none
      
      real(kind=real2), intent(in) :: b(:,:)
      
      real(kind=real2) :: ai(size(b,2),size(b,1))

      integer(kind=int1) :: m,n,i,j,lwork,lwmax
      
      real(kind=real2) :: a(size(b,1),size(b,2)),s(size(b,2))
      real(kind=real2) :: sigma(size(b,1),size(b,2)),work(10000),sigmi(size(b,2),size(b,1))
      real(kind=real2) :: u(size(b,1),size(b,1)),ut(size(b,1),size(b,1))
      real(kind=real2) :: v(size(b,2),size(b,2)),vt(size(b,2),size(b,2))
      real(kind=real2), parameter :: eps = 1d2*epsilon(b(1,1)) 
      
      m = size(b,1)
      n = size(b,2)

      if(n .gt. m) call log(event_type=5,event='PSEUDO INVERSE DIMENSION INVALID')
      
      a = b
      
      lwork = -1
      lwmax = 10000
      
      call dgesvd('A','A',m,n,a,m,s,u,m,vt,n,work,lwork,i)
      lwork = min(lwmax,int(work(1)))
      
      call dgesvd('A','A',m,n,a,m,s,u,m,vt,n,work,lwork,i)
      
      sigma = 0d0
      do i=1,n
         sigma(i,i) = s(i)
      enddo
      
      ! find pseudo inverse of sigma
      sigmi = transpose(sigma)
      do i=1,n
         if(sigmi(i,i) .lt. eps) then
            sigmi(i,i) = 0d0
         else
            sigmi(i,i) = 1d0/sigmi(i,i)
         endif
      enddo
      
      v = transpose(vt)
      ut = transpose(u)
      
      ai = matmul(v,matmul(sigmi,ut))
      
      return  
   end function pseudo_inverse
   !**********************************************************************
   subroutine random_int(harvest)
      ! Use random_number to creat a random positive integer  
      use precision
      implicit none

      integer(kind=int1), intent(out) :: harvest

      real(kind=real2) :: a

      real(kind=real2) :: b,c

      b = 2d0**32d0; c = b/2d0

      call random_number(a)
      a = a*b
      a = a - c
      harvest = nint(a)

      return
   end subroutine random_int
   !**********************************************************************
   subroutine random_int_range(i1,i2,harvest)
      ! Use random_number to creat a random positive integer  
      use precision
      implicit none

      integer(kind=int1), intent(in) :: i1,i2
      integer(kind=int1), intent(out) :: harvest

      real(kind=real2) :: a

      real(kind=real2) :: di

      di = dble(i2-i1)

      call random_number(a)
      a = a*di
      a = a + dble(i1)
      harvest = nint(a)

      return
   end subroutine random_int_range
   !**********************************************************************
   subroutine siftdown(start,bottom,a)
      use precision
      implicit none
      
      integer(kind=int1), intent(in) :: start,bottom

      real(kind=real2), intent(inout) :: a(0:)

      integer(kind=int1) :: child,root
      
      root = start
      do while(2*root + 1 .lt. bottom)
         child = 2*root + 1
         
         if(child + 1 .lt. bottom) then
            if (a(child) .lt. a(child+1)) child = child + 1
         endif
         
         if(a(root) .lt. a(child)) then
            call swap_entries(a(child),a(root))
            root = child
         else
            return
         endif
      enddo

      return
   end subroutine siftdown
   !**********************************************************************
   subroutine siftdown_array(start,bottom,x,a)
      use precision
      implicit none
      
      integer(kind=int1), intent(in) :: start,bottom

      real(kind=real2), intent(inout) :: x(0:)
      real(kind=real2), intent(inout) :: a(0:)

      integer(kind=int1) :: child,root
      
      root = start
      do while(2*root + 1 .lt. bottom)
         child = 2*root + 1
         
         if(child + 1 .lt. bottom) then
            if (a(child) .lt. a(child+1)) child = child + 1
         endif
         
         if(a(root) .lt. a(child)) then
            call swap_entries(a(child),a(root))
            call swap_entries(x(child),x(root))            
            root = child
         else
            return
         endif
      enddo

      return
   end subroutine siftdown_array
   !**********************************************************************
   pure function simplicy_number(n,d) result(s)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n,d

      integer(kind=int1) :: s

      s = nint(pochhammerR(1d0*n,d))/factorial(d)
      
      return
   end function simplicy_number
   !**********************************************************************
   recursive subroutine sortrows(n,m,a) ! with bubble based sort
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n,m 

      integer(kind=int1), intent(inout) :: a(n,m)

      integer(kind=int1) :: temp(m)
      integer(kind=int1) :: i,j,k
      logical :: swapped
  
      do j=n-1,1,-1

         swapped = .false.
         do i=1,j
            if(a(i,1) .gt. a(i+1,1))then
               temp(:) = a(i,:)
               a(i,:) = a(i+1,:)
               a(i+1,:) = temp(:)
               swapped = .true.
            endif
            
            if((a(i,1) .eq. a(i+1,1)) .and. (m .gt. 1))then
               do k=0,i-1
                  if((a(i-k,1) .ne. a(i-k+1,1))) exit
               enddo
               k = k - 1
               call sortrows(2+k,m-1,a(i-k:i+1,2:))
            endif
         enddo
         if(.not. swapped) exit
      enddo

      return
   end subroutine sortrows
   !**********************************************************************
   subroutine sortrows_flip(n,m,a) ! with bubble based sort
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n,m 

      integer(kind=int1), intent(inout) :: a(n,m)

      integer(kind=int1) :: b(n,m)

      b = a(:,m:1:-1)
      
      call sortrows(n,m,b)
      a = b(:,m:1:-1)

      return
   end subroutine sortrows_flip
   !**********************************************************************
   recursive subroutine sortrows_real(n,m,a) ! with bubble based sort
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n,m 

      real(kind=real2), intent(inout) :: a(n,m)

      integer(kind=int1) :: i,j,k

      real(kind=real2) :: temp(m)
      real(kind=real2), parameter :: tol = 100d0*epsilon(a(1,1))

      logical :: swapped
  
      do j=n-1,1,-1
         swapped = .false.
         do i=1,j
            if((a(i,1) .gt. a(i+1,1)) .and. (abs(a(i,1)-a(i+1,1)) .gt. tol))then
               temp(:) = a(i,:)
               a(i,:) = a(i+1,:)
               a(i+1,:) = temp(:)
               swapped = .true.
            endif
            
            if((abs(a(i,1)-a(i+1,1)) .lt. tol)  .and. (m .gt. 1))then
               do k=0,i-1
                  if((a(i-k,1) .ne. a(i-k+1,1))) exit
               enddo
               k = k - 1
               call sortrows_real(2+k,m-1,a(i-k:i+1,2:))
            endif
         enddo
         if(.not. swapped) exit
      enddo

      return
   end subroutine sortrows_real
   !**********************************************************************
   subroutine sortrows_real_flip(n,m,a) ! with bubble based sort
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n,m 

      real(kind=real2), intent(inout) :: a(n,m)

      real(kind=real2) :: b(n,m)

      b = a(:,m:1:-1)
      
      call sortrows_real(n,m,b)
      a = b(:,m:1:-1)

      return
   end subroutine sortrows_real_flip
   !**********************************************************************
   pure subroutine swap_entries(a,b)
      use precision
      implicit none

      real(kind=real2), intent(inout) :: a,b

      real(kind=real2) :: temp

      temp = a
      a = b
      b = temp
      
      return
   end subroutine swap_entries
   !********************************************************************** 
end module maths
!*************************************************************************
