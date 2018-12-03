module quadrature
   use precision
   implicit none

   private
   
   real(kind=real2), parameter :: pi = 4d0*atan(1d0)

   public :: quad
   public :: gaussj_nodes,gaussj_weights,gauss_Jacobi
   public :: gaussl_nodes,gaussl_weights,gauss_Legendre
   public :: lobatto_nodes,lobatto_weights
   public :: kronrod_nodes,kronrod_weights
   public :: christoffel_nodes,christoffel_weights,gauss_Christoffel
   public :: clenshaw_curtis
   
contains
   !**********************************************************************
   !> @brief Calculate Gauss-Legendre Nodes
   !> @par Calculates the Gauss-Legendre nodes using Gauss-Jacobi
   !> @param[in] n order (number of points)
   !> @param[out] x nodes
   !**********************************************************************
   function gaussl_nodes(n) result(x)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n
      
      real(kind=real2) :: x(n)

      real(kind=real2) :: w(n)
      
      call quad('legendre',n,x,w,-1d0,1d0)
      
      return
   end function gaussl_nodes
   !**********************************************************************
   !> @brief Calculate Gauss-Legendre Weights
   !> @par Calculates the Gauss-Legendre weights using Gauss-Jacobi
   !> @param[in] n order (number of points)
   !> @param[out] w weights
   !**********************************************************************
   function gaussl_weights(n) result(w)
      use precision
      implicit none
     
      integer(kind=int1), intent(in) :: n
      
      real(kind=real2) :: w(n)

      real(kind=real2) :: x(n)
      
      call quad('legendre',n,x,w,-1d0,1d0)
      
      return
   end function gaussl_weights
   !**********************************************************************
   !> @brief Calculate Gauss-Legendre Quadrature
   !> @par Calculates the Gauss-Legendre nodes and weights using Gauss-Jacobi
   !> @param[in] n order (number of points)
   !> @param[out] x nodes
   !> @param[out] w weights
   !**********************************************************************
   subroutine gauss_Legendre(n,x,w)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n

      real(kind=real2), intent(out) :: x(n),w(n)

      call quad('legendre',n,x,w,-1d0,1d0)

      return
   end subroutine gauss_Legendre
   !**********************************************************************
   !> @brief Calculate Gauss-Jacobi Nodes
   !> @par Calculates the Gauss-Jacobi nodes using Golub-Welsch 
   !> @param[in] n order (number of points)
   !> @param[in] alpha Jacobi weight function alpha parameter 
   !> @param[in] beta Jacobi weight function beta parameter
   !> @param[out] x nodes
   !**********************************************************************
   function gaussj_nodes(n,alpha,beta) result(x)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n
      real(kind=real2), intent(in) :: alpha,beta

      real(kind=real2) :: x(n)

      real(kind=real2) :: w(n)

      call gauss_Jacobi(n,alpha,beta,x,w)

      return
   end function gaussj_nodes
   !**********************************************************************
   !> @brief Calculate Gauss-Jacobi Weights
   !> @par Calculates the Gauss-Jacobi weights using Golub-Welsch
   !> @param[in] n order (number of points)
   !> @param[in] alpha Jacobi weight function alpha parameter 
   !> @param[in] beta Jacobi weight function beta parameter
   !> @param[out] w weights
   !**********************************************************************
   function gaussj_weights(n,alpha,beta) result(w)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n
      real(kind=real2), intent(in) :: alpha,beta

      real(kind=real2) :: w(n)

      real(kind=real2) :: x(n)

      call gauss_Jacobi(n,alpha,beta,x,w)

      return
   end function gaussj_weights
   !**********************************************************************
   !> @brief Calculate Gauss-Lobatto Nodes
   !> @par Calculates the Gauss-Lobatto nodes using Gauss-Legendre
   !> @param[in] n order (number of points)
   !> @param[out] x nodes
   !**********************************************************************
   function lobatto_nodes(n) result(x)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n

      real(kind=real2) :: x(n)

      if(n .ge. 3)then
         x(1) = -1d0
         x(2:n-1) = gaussl_nodes(n-2)
         x(n) = 1d0
      elseif(n .eq. 2)then
         x(1) = -1d0
         x(n) = 1d0
      else
         print *,'ERROR :: INVALID LOBATTO QUAD ORDER'
         stop
      endif
      
      return
   end function lobatto_nodes
   !**********************************************************************
   !> @brief Calculate Gauss-Lobatto Weights
   !> @par Calculates the Gauss-Lobatto nodes using Gauss-Legendre
   !> @param[in] n order (number of points)
   !> @param[out] w weight
   !**********************************************************************
   function lobatto_weights(n) result(w)
      use precision
      use polynomial, only : legendre
      implicit none

      integer(kind=int1), intent(in) :: n

      real(kind=real2) :: w(n)

      integer(kind=int2) :: in
      real(kind=real2) :: x(n),rnnm1

      x = lobatto_nodes(n)

      rnnm1 = 2d0/dble(n*(n-1d0))
      do in=1,n
         w(in) = 1d0/(legendre(x(in),n-1))
         w(in) = rnnm1*w(in)*w(in)
      enddo
      
      return
   end function lobatto_weights
   !**********************************************************************
   function kronrod_nodes(n) result(x)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n

      real(kind=real2) :: x(n)

      integer(kind=int1) :: m,in,i
      real(kind=real2), allocatable :: x_in(:),w1_in(:),w2_in(:)

      m = n/2
      allocate(x_in(m+1)); allocate(w1_in(m+1)); allocate(w2_in(m+1))

      call kronrod(m,x_in,w1_in,w2_in)

      i = 0
      do in=1,m+mod(m,2)
         i = i + 1
         x(i) = x_in(in)
      enddo
      do in=m,1,-1
         i = i + 1
         x(i) = -x_in(in)
      enddo
      
      return
   end function kronrod_nodes
   !**********************************************************************
   function kronrod_weights(n) result(w)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n

      real(kind=real2) :: w(n)

      integer(kind=int1) :: m,in,i
      real(kind=real2), allocatable :: x_in(:),w1_in(:),w2_in(:)

      m = n/2
      allocate(x_in(m+1)); allocate(w1_in(m+1)); allocate(w2_in(m+1))

      call kronrod(m,x_in,w1_in,w2_in)

      i = 0
      do in=1,m+mod(m,2)
         i = i + 1
         w(i) = w1_in(in)
      enddo
      do in=m,1,-1
         i = i + 1
         w(i) = w1_in(in)
      enddo
      
      return
   end function kronrod_weights
   !**********************************************************************
   function christoffel_nodes(n,r) result(x)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n

      real(kind=real2), intent(in) :: r(:)

      real(kind=real2) :: x(n)

      real(kind=real2) :: w(n)

      call gauss_Christoffel(n,r,x,w)

      return
   end function christoffel_nodes
   !**********************************************************************
   function christoffel_weights(n,r) result(w)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n

      real(kind=real2), intent(in) :: r(:)

      real(kind=real2) :: w(n)

      real(kind=real2) :: x(n)

      call gauss_Christoffel(n,r,x,w)

      return
   end function christoffel_weights
   !**********************************************************************
   !> @breif computes Gauss-Christoffel quadrature for measure desrcibed in r
   !
   !> @par this computes the nodes and weights according to gauss-christoffel
   !! for a polynomial measure set out in r. The form or r are the monomial
   !! coefficents of the measure. Such that:
   !! measure = r(1)*x^0 + r(2)*x^1 + r(3)*x^2 + ...
   !! This assumes the range of integration is [-1,1] 
   !! Complexity is O(n^3). Could be better, but ok for n<~500.
   !
   !> @param[in] n order
   !> @param[in] r array of monomial coefficents of measure
   !> @param[out] x nodes
   !> @param[out] w weights
   !**********************************************************************
   subroutine gauss_Christoffel(n,r,x,w)
      use precision
      use maths, only : eigen_std
      implicit none

      integer(kind=int1), intent(in) :: n
      real(kind=real2), intent(in) :: r(:)

      real(kind=real2), intent(out) :: x(n),w(n)

      integer(kind=int1) :: nq,i
      real(kind=real2), allocatable :: xg(:),wg(:)
      real(kind=real2) :: a(n),b(n),mu0,v(n,n)
      
      nq = ceiling(dble(2*(n-1)+size(r))/2d0) + 1

      allocate(wg(nq)); allocate(xg(nq))
      call gauss_Legendre(nq,xg,wg)
      
      call christoffel_ab(n,r,xg,wg,a,b,mu0)

      call sgqf(n,a,b,mu0,x,w)
              
      return
   end subroutine gauss_Christoffel
   !**********************************************************************
   !> @brief calculates the normalised polynomial coefficents of polynomial
   !! recurrsion
   !
   !> @param[in] n order
   !> @param[in] w quadrature measure
   !> @param[in] xg quadrature points used to integrate
   !> @param[in] wg quadrature weights used to integrate
   !> @param[out] a first recurrsion coefficient
   !> @param[out] b second recurrsion coefficient
   !> @param[out] mu0 zeroth order moment of measure
   !**********************************************************************
   subroutine christoffel_ab(n,w,xg,wg,a,b,mu0)
      use precision
      use polynomial, only : inproduct,polyval,poly_recursion_single
      implicit none

      integer(kind=int1), intent(in) :: n
      real(kind=real2), intent(in) :: w(:),xg(:),wg(:)

      real(kind=real2), intent(out) :: a(n),b(n),mu0

      integer(kind=int1) :: i,nx
      
      real(kind=real2) :: at(n+1),bt(n+1)

      real(kind=real2) :: int,xint,intm1
      real(kind=real2) :: wx(size(xg)),ww(size(w)+1),wxx(size(xg))
      real(kind=real2) :: pm1(size(xg)),p(size(xg)),pp1(size(xg))
      
      real(kind=real2), parameter :: tol = sqrt(epsilon(int))
      
      nx = size(xg)
      ww(1) = 0d0
      ww(2:) = w(1:)
      wx  = polyval(w,xg)
      wxx = polyval(ww,xg)

      pm1(:) = 0d0
      p(:) = 1d0      
      
      int = inproduct(nx,p,wx,xg,wg); xint = inproduct(nx,p,wxx,xg,wg)
      
      mu0 = int
      at(1) = xint/int
      bt(1) = 1d0
      
      do i=1,n
         intm1 = int

         pp1 = poly_recursion_single(at(i),bt(i),xg,p,pm1)        
         
         pm1 = p
         p = pp1

         int  = inproduct(nx,p,wx ,xg,wg)
         xint = inproduct(nx,p,wxx,xg,wg)
         
         at(i+1) = xint/int
         bt(i+1) = int/intm1
      enddo

      a(1:n) = at(1:n)
      do i=1,n         
         if(bt(i+1) .gt. tol) then
            b(i) = sqrt(abs(bt(i+1)))
         else
            b(i) = 0d0
         endif
      enddo
            
      return
   end subroutine christoffel_ab
   !**********************************************************************
   !> @breif computes a Clenshaw Curtis quadrature rule.
   !
   !> @par  Computes a Clenshaw Curtis quadrature rule with the convention
   !! that the abscissas are numbered from left to right
   !
   !> @author John Burkardt, Will Trojak
   !
   !> @param[in] order the order of the 1D quadrature
   !> @param[out] x the abscissas
   !> @param[out] w the weights
   !**********************************************************************
   subroutine clenshaw_curtis(order,x,w)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: order
      
      real(kind=real2), intent(out) :: x(order),w(order)
      
      real(kind=real2), parameter :: pi = 4d0*atan(1d0)

      integer(kind=int1) :: i,j
      
      real(kind=real2) :: theta,b,rom1
      
      if(order .lt. 1) then
         print *,'ERROR: CLENSHAW QUAD ORDER TOO SMALL'
         stop
      endif
      
      if(order .eq. 1) then
         x(1) = 0d0
         w(1) = 2d0
         return
      endif

      rom1 = 1d0/real(order - 1)
      
      do i=1,order
         x(i) = cos(real(order - i)*pi*rom1)
         theta = real(i - 1)*pi*rom1
         
         w(i) = 1d0

         b = 2d0
         do j=1,(order - 1)/2-1            
            w(i) = w(i) - b*cos(2d0*real(j)*theta)/real(4*j*j - 1)
         enddo
         b = 1d0
         j = (order - 1)/2
         w(i) = w(i) - b*cos(2d0*real(j)*theta)/real(4*j*j - 1)
         
         w(i) = 2d0*w(i)*rom1
      enddo

      x(1) = -1d0
      if(mod(order,2) .eq. 1)then
         x((order+1)/2) = 0d0
      endif
      x(order) = 1d0
      
      w(1) = 0.5d0*w(1)
      
      !w(2:order-1) = 2d0*w(2:order-1)/real(order - 1)
      
      w(order) = 0.5d0*w(order)    
      
      return
   end subroutine clenshaw_curtis
   !**********************************************************************
   subroutine gauss_Jacobi(n,alpha,beta,x,w)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n
      real(kind=real2), intent(in) :: alpha,beta

      real(kind=real2), intent(out) :: x(n),w(n)

      integer(kind=int1) :: qtype
      real(kind=real2) :: a,b

      a = -1d0; b = 1d0

      qtype = 4

      if(alpha .le. -1d0) then
         print *,'ERROR: JACOBI POLY \alpha < -1'
         stop
      endif
      if(beta  .le. -1d0) then
         print *,'ERROR: JACOBI POLY \beta < -1'
         stop
      endif
      call cgqf(n,qtype,alpha,beta,a,b,x,w)
      
      return
   end subroutine gauss_Jacobi
   !**********************************************************************
   subroutine quad(qtype_str,n,x,w,a0,b0,c0,d0)
      use precision
      implicit none

      character(*), intent(in) :: qtype_str

      integer(kind=int1), intent(in) :: n

      real(kind=real2), intent(out) :: x(n),w(n) 

      real(kind=real2), optional, intent(in) :: a0,b0,c0,d0

      integer(kind=int1) :: qtype
      real(kind=real2) :: a,b,c,d

      a = 0d0
      b = 0d0
      if(present(a0)) a = a0
      if(present(b0)) b = b0
      
      select case(qtype_str)
      case('legendre') ! x = (c,d), w = 1.0
         qtype = 1

         c = -1d0; d = 1d0
         if(present(c0)) c = c0
         if(present(d0)) d = d0
      case('chebyshev') ! x = (c,d), w = ((d-x)*(x-c))^(-0.5)
         qtype = 2
         
         c = -1d0; d = 1d0
         if(present(c0)) c = c0
         if(present(d0)) d = d0
      case('gegenbauer') ! x = (c,d), w = ((d-x)*(x-c))^a
         qtype = 3
         
         c = -1d0; d = 1d0
         if(present(c0)) c = c0
         if(present(d0)) d = d0
      case('jacobi') ! x = (c,d), w = (d-x)^a*(x-c)^b
         qtype = 4
         
         c = -1d0; d = 1d0
         if(present(c0)) c = c0
         if(present(d0)) d = d0
      case('laguerre') ! x = (c,inf), w = (x-c)^a*exp(-d*(x-c))
         qtype = 5
         
         c = 0d0; d = 0d0
         if(present(c0)) c = c0
         if(present(d0)) d = d0
      case('hermite') ! x = (-inf,inf), w = |x-c|^a*exp(-d*(x-c)^2)
         qtype = 6
         
         c = 0d0; d = 0d0
         if(present(c0)) c = c0
         if(present(d0)) d = d0
      case('exponential') ! x = (c,d), w =  |x-(c+d)/2.0|^a
         qtype = 7
         
         c = -1d0; d = 1d0
         if(present(c0)) c = c0
         if(present(d0)) d = d0
      case('rational') ! x = (c,inf), w =  (x-c)^a*(x+d)^b         
         qtype = 8

         c = 0d0; d = 0d0
         if(present(c0)) c = c0
         if(present(d0)) d = d0
      case default
         qtype = 1

         !print *,'defualting to Guass-Legendre'
         
         c = -1d0; d = 1d0
         if(present(c0)) c = c0
         if(present(d0)) d = d0
      end select
      
      call cgqf(n,qtype,a,b,c,d,x,w)

      return
   end subroutine quad
   !**********************************************************************
   !! CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
   !
   !  Discussion:
   !
   !    This routine computes all the knots and weights of a Gauss quadrature
   !    formula with a classical weight function with default values for A and B,
   !    and only simple knots.
   !
   !    There are no moments checks and no printing is done.
   !
   !    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
   !  Author: Sylvan Elhay, Jaroslav Kautsky and, John Burkardt (F90)
   !
   !  Parameters:
   !
   !    Input, integer(kind=int1) :: NT, the number of knots.
   !    Input, integer(kind=int1) :: KIND, the rule.
   !    1, Legendre,             (a,b)       1.0
   !    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
   !    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
   !    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
   !    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
   !    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
   !    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
   !    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
   !
   !    Input, real(kind=real2) :: ALPHA, the value of Alpha, if needed.
   !    Input, real(kind=real2) :: BETA, the value of Beta, if needed.
   !    Output, real(kind=real2) :: T(NT), the knots.
   !    Output, real(kind=real2) :: WTS(NT), the weights.
   !**********************************************************************
   subroutine cdgqf(nt,qtype,alpha,beta,t,wts)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: nt,qtype

      real(kind=real2), intent(in) :: alpha,beta

      real(kind=real2), intent(out) :: t(nt),wts(nt)

      real(kind=real2) :: zemu,aj(nt),bj(nt)

      ! Get the Jacobi matrix and zero-th moment.
      call class_matrix(qtype,nt,alpha,beta,aj,bj,zemu)

      ! Compute the knots and weights.
      call sgqf(nt,aj,bj,zemu,t,wts)

      return
   end subroutine cdgqf
   !**********************************************************************
   !! CGQF computes knots and weights of a Gauss quadrature formula.
   !
   !  Licensing:
   !    This code is distributed under the GNU LGPL license. 
   !
   !  Author: Sylvan Elhay, Jaroslav Kautsky and, John Burkardt (F90)
   !
   !  Reference:
   !    Sylvan Elhay, Jaroslav Kautsky,
   !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
   !    Interpolatory Quadrature,
   !    ACM Transactions on Mathematical Software,
   !    Volume 13, Number 4, December 1987, pages 399-415.
   !
   !  Parameters:
   !
   !    Input, integer(kind=int1) :: NT, the number of knots.
   !    Input, real(kind=real2) :: ALPHA, the value of Alpha, if needed.
   !    Input, real(kind=real2) :: BETA, the value of Beta, if needed.
   !    Input, real(kind=real2) :: A, B, the interval endpoints, or
   !    other parameters.
   !    Output, real(kind=real2) :: T(NT), the knots.
   !    Output, real(kind=real2) :: WTS(NT), the weights.
   !**********************************************************************
   subroutine cgqf(nt,qtype,alpha,beta,a,b,t,wts)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: nt,qtype

      real(kind=real2), intent(in) :: alpha,beta,a,b

      real(kind=real2), intent(out) :: t(nt),wts(nt)

      integer(kind=int1) :: i
      integer(kind=int1), allocatable :: mlt(:),ndx(:)

      ! Compute the Gauss quadrature for default a & b
      call cdgqf(nt,qtype,alpha,beta,t,wts)

      !  Prepare to scale the quadrature to other weight function  a & b
      allocate(mlt(nt))
      mlt = 1

      allocate(ndx(nt))
      do i=1,nt 
         ndx(i) = i
      enddo

      call scqf(nt,t,mlt,wts,nt,ndx,wts,t,qtype,alpha,beta,a,b)

      deallocate(mlt)
      deallocate(ndx)

      return
   end subroutine cgqf
   !**********************************************************************
   !! CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
   !
   !  Discussion:
   !
   !    This routine computes the diagonal AJ and sub-diagonal BJ
   !    elements of the order M tridiagonal symmetric Jacobi matrix
   !    associated with the polynomials orthogonal with respect to
   !    the weight function specified by KIND.
   !
   !    For weight functions 1-7, M elements are defined in BJ even
   !    though only M-1 are needed.  For weight function 8, BJ(M) is
   !    set to zero.
   !
   !    The zero-th moment of the weight function is returned in ZEMU.
   !
   !  Author: Sylvan Elhay, Jaroslav Kautsky and, John Burkardt (F90)
   !
   !  Reference:
   !    Sylvan Elhay, Jaroslav Kautsky,
   !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
   !    Interpolatory Quadrature,
   !    ACM Transactions on Mathematical Software,
   !    Volume 13, Number 4, December 1987, pages 399-415.
   !
   !  Parameters:
   !    Input, integer(kind=int1) :: M, the order of the Jacobi matrix.
   !    Input, real(kind=real2) :: ALPHA, the value of Alpha, if needed.
   !    Input, real(kind=real2) :: BETA, the value of Beta, if needed.
   !    Output, real(kind=real2) :: AJ(M), BJ(M), the diagonal and subdiagonal
   !    of the Jacobi matrix.
   !    Output, real(kind=real2) :: ZEMU, the zero-th moment.
   !**********************************************************************
   subroutine class_matrix(qtype,m,alpha,beta,aj,bj,zemu)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: m,qtype

      real(kind=real2), intent(in) :: alpha,beta

      real(kind=real2), intent(out) :: aj(m),bj(m),zemu

      integer(kind=int1) :: i
      real(kind=real2) :: a2b2,ab,aba,abi,abj,abti,apone
      real(kind=real2) :: temp,temp2
      real(kind=real2), parameter :: pi = 4d0*atan(1d0)

      temp = epsilon(temp)
      temp2 = 0.5d0

      if(500d0*temp .lt. abs(gamma(temp2)**2 - pi)) then
         print *,'ERROR: GAMMA FUNCTION INACCURATE'
         stop
      endif

      select case(qtype)
      case(1) ! Legendre
         ab = 0d0
         zemu = 2d0/(ab + 1d0)

         aj = 0d0
         do i=1,m
            abi = i + ab*mod(i,2)
            abj = 2*i + ab
            bj(i) = abi*abi/(abj*abj - 1d0)
         enddo
         bj(:) = sqrt(bj(:))

      case(2) ! Chebyshev
         zemu = pi

         aj = 0d0

         bj(1) = sqrt(0.5d0)
         bj(2:m) = 0.5d0
      case(3) ! Gegenbauer
         ab = alpha*2d0

         zemu = 2d0**(ab + 1d0)*gamma(alpha + 1d0)**2/gamma(ab + 2d0)

         aj = 0d0
         bj(1) = 1d0/(2d0*alpha + 3d0)
         do i=2,m
            bj(i) = i*(i + ab)/(4d0*(i + alpha)**2 - 1d0)
         enddo
         bj(:) = sqrt(bj(:))
      case(4) ! Jacobi
         ab = alpha + beta
         abi = 2d0 + ab
         zemu = 2d0**(ab + 1d0)*gamma(alpha + 1d0)*gamma(beta + 1d0)/gamma(abi)

         aj(1) = (beta - alpha)/abi
         bj(1) = 4d0*(1d0 + alpha)*(1d0 + beta)/((abi + 1d0)*abi*abi)
         a2b2 = beta*beta - alpha*alpha

         do i=2,m
            abi = 2d0*i + ab
            aj(i) = a2b2/((abi - 2d0)*abi)
            abi = abi**2
            bj(i) = 4d0*i*(i + alpha)*(i + beta)*(i + ab)/((abi - 1d0)*abi)
         enddo
         bj(:) = sqrt(bj(:))

      case(5) ! Generalized Laguerre
         zemu = gamma(alpha + 1d0)

         do i=1,m
            aj(i) = 2d0*i - 1d0 + alpha
            bj(i) = i*(i + alpha)
         enddo
         bj(:) = sqrt(bj(:))

      case(6) ! Generalized Hermite
         zemu = gamma((alpha + 1d0)/2d0)
         aj = 0d0
         
         do i=1,m
            bj(i) = 0.5d0*(i + alpha*mod(i,2))
         enddo
         bj(:) = sqrt(bj(:))

      case(7) ! Exponential
         ab = alpha
         zemu = 2d0/(ab + 1d0)
         aj = 0d0

         do i=1,m
            abi = i + ab*mod(i,2)
            abj = 2*i + ab
            bj(i) = abi*abi/(abj*abj - 1d0)
         enddo
         bj(:) = sqrt(bj(:))

      case(8) ! Rational
         ab = alpha + beta

         zemu = gamma(alpha + 1d0)*gamma(-(ab + 1d0))/gamma(-beta)
         apone = alpha + 1d0
         aba = ab*apone
         aj(1) = -apone/(ab + 2d0)
         bj(1) = -aj(1)*(beta + 1d0)/(ab + 2d0)/(ab + 3d0)
         do i=2,m
            abti = ab + 2d0*i
            aj(i) = aba + 2d0*(ab + i)*(i - 1)
            aj(i) = -aj(i)/abti/(abti - 2d0)
         enddo

        do i=2,m-1
           abti = ab + 2d0*i
           bj(i) = i*(alpha + i)/(abti - 1d0)*(beta + i)/(abti**2)*(ab + i)/(abti + 1d0)
        enddo

        bj(m) = 0d0
        bj(:) = sqrt(bj(:))

      end select

      return
   end subroutine class_matrix
   !**********************************************************************
   !! IMTQLX diagonalizes a symmetric tridiagonal matrix.
   !
   !  Discussion:
   !
   !    This routine is a  modified version of EISPACK implicit QL algorithm for
   !    symmetric tridiagonal matrix. 
   !
   !    It is modified to produce the product Q' * Z, where Z is an input 
   !    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
   !    The changes consist (essentially) of applying the orthogonal 
   !    transformations directly to Z as they are generated.
   !
   !  Author: Sylvan Elhay, Jaroslav Kautsky and, John Burkardt (F90)
   !
   !  Reference:
   !
   !    Sylvan Elhay, Jaroslav Kautsky,
   !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
   !    Interpolatory Quadrature,
   !    ACM Transactions on Mathematical Software,
   !    Volume 13, Number 4, December 1987, pages 399-415.
   !
   !    Roger Martin, James Wilkinson,
   !    The Implicit QL Algorithm,
   !    Numerische Mathematik,
   !    Volume 12, Number 5, December 1968, pages 377-383.
   !
   !  Parameters:
   !    Input, integer(kind=int1) :: N, the order of the matrix.
   !    Input/output, real(kind=real2) :: D(N), the diagonal entries of the matrix.
   !    On output, the information in D has been overwritten.
   !
   !    Input/output, real(kind=real2) :: E(N), the subdiagonal entries of the 
   !    matrix, in entries E(1) through E(N-1).  On output, the information in
   !    E has been overwritten.
   !
   !    Input/output, real(kind=real2) :: Z(N).  On input, a vector.  On output,
   !    the value of Q' * Z, where Q is the matrix that diagonalizes the
   !    input symmetric tridiagonal matrix.
   !**********************************************************************
   subroutine imtqlx(n,d,e,z)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n

      real(kind=real2), intent(inout) :: d(n),e(n),z(n)

      integer(kind=int1) :: i,ii,j,k,l,m,mml
      integer(kind=int1), parameter :: itn = 30
      real(kind=real2) :: b,c,f,g,p,prec,r,s

      prec = epsilon(prec)

      if(n .eq. 1) return

      e(n) = 0d0
      do l=1,n
         j = 0
         
         do
            do m=l,n-1
               !if(m .eq. n) exit
               if(abs(e(m)) .le. prec*(abs(d(m)) + abs(d(m+1)))) exit
            enddo
            
            p = d(l)
            if(m .eq. l) exit
            
            if(itn .le. j) then
               print *,'ERROR: IMTQLX CONV FAIL'
               stop
            endif
            
            j = j + 1
            g = (d(l+1) - p)/(2d0*e(l))
            r =  sqrt(g*g + 1d0)
            g = d(m) - p + e(l)/(g + sign(r,g))
            s = 1d0; c = 1d0; p = 0d0
            mml = m - l
            
            do ii=1,mml
               i = m - ii
               f = s*e(i); b = c*e(i)
               
               if(abs(g) .le. abs(f)) then
                  c = g/f
                  r = sqrt(c*c + 1d0)
                  e(i+1) = f*r
                  s = 1d0/r
                  c = c*s
               else
                  s = f/g
                  r = sqrt(s*s + 1d0)
                  e(i+1) = g*r
                  c = 1d0/r
                  s = s*c
               end if
               
               g = d(i+1) - p
               r = (d(i) - g)*s + 2d0*c*b
               p = s*r
               d(i+1) = g + p
               g = c*r - b
               f = z(i+1)
               
               z(i+1) = s*z(i) + c*f
               z(i) = c*z(i) - s*f
            enddo
            
            d(l) = d(l) - p; e(l) = g; e(m) = 0d0
         enddo
      enddo

      ! Sorting.
      do ii=2,n
         i = ii - 1; k = i
         p = d(i)
         
         do j=ii,n
            if(d(j) .lt. p) then
               k = j; p = d(j)
            endif
         enddo
         
         if(k .ne. i) then
            d(k) = d(i); d(i) = p
            p = z(i)
            z(i) = z(k); z(k) = p
         endif

      enddo

      return
   end subroutine imtqlx
   !**********************************************************************
   !! SCQF scales a quadrature formula to a nonstandard interval.
   !
   !  Discussion: The arrays WTS and SWTS may coincide, the arrays T and ST may coincide.
   !
   !  Author: Sylvan Elhay, Jaroslav Kautsky and, John Burkardt (F90)
   !
   !  Reference:
   !    Sylvan Elhay, Jaroslav Kautsky,
   !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
   !    Interpolatory Quadrature,
   !    ACM Transactions on Mathematical Software,
   !    Volume 13, Number 4, December 1987, pages 399-415.
   !
   !  Parameters:
   !    Input, integer(kind=int1) :: NT, the number of knots.
   !    Input, real(kind=real2) :: T(NT), the original knots.
   !    Input, integer(kind=int1) :: MLT(NT), the multiplicity of the knots.
   !    Input, real(kind=real2) :: WTS(NWTS), the weights.
   !    Input, integer(kind=int1) :: NWTS, the number of weights.
   !    Input, integer(kind=int1) :: NDX(NT), used to index the array WTS.  
   !    For more details see the comments in CAWIQ.
   !    Output, real(kind=real2) :: SWTS(NWTS), the scaled weights.
   !    Output, real(kind=real2) :: ST(NT), the scaled knots.
   !    Input, integer(kind=int1) :: KIND, the rule.
   !    Input, real(kind=real2) :: ALPHA, the value of Alpha, if needed.
   !    Input, real(kind=real2) :: BETA, the value of Beta, if needed.
   !    Input, real(kind=real2) :: A, B, the interval endpoints.
   !**********************************************************************
   subroutine scqf(nt,t,mlt,wts,nwts,ndx,swts,st,qtype,alpha,beta,a,b)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: nt,nwts,qtype
      integer(kind=int1), intent(in) :: mlt(nt),ndx(nt)

      real(kind=real2), intent(in) :: alpha,beta,a,b
      real(kind=real2), intent(in) :: t(nt),wts(nwts)

      real(kind=real2), intent(out) :: swts(nwts),st(nt)

      integer(kind=int1) :: i,k,l
      real(kind=real2) :: al,be,p,shft,slp,temp,tmp
      logical :: quality

      temp = epsilon(temp)

      quality = .true.
      select case(qtype)
      case(1,2,3,4,7,9)
         if(qtype .eq. 1) al = 0d0; be = 0d0
         if(qtype .eq. 2) al = -0.5d0; be = -0.5d0
         if(qtype .eq. 3) al = alpha; be = alpha
         if(qtype .eq. 4) al = alpha; be = beta
         if(qtype .eq. 7) al = alpha; be = 0d0
         if(qtype .eq. 9) al = 0.5d0; be = 0.5d0

         shft = 0.5d0*(a + b); slp = 0.5d0*(b - a)
         if(abs(slp) .le. temp) quality = .false.
      case(5)
         al = alpha; be = 0d0

         shft = a; slp = 1d0/b
         if(b .le. 0d0) quality = .false.
      case(6)
         al = alpha; be = 0d0

         shft = a; slp = 1d0/sqrt(b)
         if(b .le. 0d0) quality = .false.
      case(8)
         al = alpha; be = beta

         shft = a; slp = a + b
         if(slp .le. 0d0) quality = .false.
      end select

      if(.not. quality) then
         print *,'ERROR: SCQF QUAL ERROR'
         stop
      endif

       p = slp**(al + be + 1d0)

      do k=1,nt
         st(k) = shft + slp*t(k)
         l = abs(ndx(k))
         
         if(l .ne. 0) then
            tmp = p
            do i=l,l+mlt(k)-1
               swts(i) = wts(i)*tmp
               tmp = tmp*slp
            enddo
         endif
      enddo

      return
   end subroutine scqf
   !**********************************************************************
   !! SGQF computes knots and weights of a Gauss Quadrature formula.
   !  Discussion:
   !    This routine computes all the knots and weights of a Gauss quadrature
   !    formula with simple knots from the Jacobi matrix and the zero-th
   !    moment of the weight function, using the Golub-Welsch technique.
   !
   !  Author: Sylvan Elhay, Jaroslav Kautsky and, John Burkardt (F90)
   !
   !  Reference:
   !    Sylvan Elhay, Jaroslav Kautsky,
   !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
   !    Interpolatory Quadrature,
   !    ACM Transactions on Mathematical Software,
   !    Volume 13, Number 4, December 1987, pages 399-415.
   !
   !  Parameters:
   !    Input, integer(kind=int1) :: NT, the number of knots.
   !    Input, real(kind=real2) :: AJ(NT), the diagonal of the Jacobi matrix.
   !    Input/output, real(kind=real2) :: BJ(NT), the subdiagonal of the Jacobi 
   !    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
   !    Input, real(kind=real2) :: ZEMU, the zero-th moment of the weight function.
   !    Output, real(kind=real2) :: T(NT), the knots.
   !    Output, real(kind=real2) :: WTS(NT), the weights.
   !**********************************************************************
   subroutine sgqf(nt,aj,bj,zemu,t,wts)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: nt
      real(kind=real2), intent(in) :: aj(nt),zemu

      real(kind=real2), intent(out) :: t(nt),wts(nt)

      real(kind=real2), intent(inout) :: bj(nt)

      integer(kind=int1) :: i

      ! Exit if the zero-th moment is not positive.
      if(zemu .le. 0d0) then
         print *,'ERROR: ZEROTH MOMENT <= 0'
         stop
      endif
      
      ! Set up vectors for IMTQLX.
      t(:) = aj(:)

      wts(1) = sqrt(zemu)
      wts(2:nt) = 0d0

      ! Diagonalize the Jacobi matrix.
      call imtqlx(nt,t,bj,wts)
      wts(:) = wts(:)**2

      return
   end subroutine sgqf
   !**********************************************************************
   !> @author Robert Piessens (F77), Maria Branders(F77) and John Burkardt (F90)
   ! 
   !> @brief Calculate the n^th order Gauss-Kronrod quadrature nodes and weights.
   !
   !> @param[in] n the order of the Gauss rule
   !> @param[in] eps the requested absolute accuracy of the abscissas
   !> @param[out] x array of abscissas
   !> @param[out] w1 array of Gauss-Kronrod weights
   !> @param[out] w2 array Gauss weights 
   !**********************************************************************
   subroutine kronrod(n,x,w1,w2)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n

      real(kind=real2), intent(out) :: x(n+1),w1(n+1),w2(n+1)

      integer(kind=int1), parameter :: niter = 100
      integer(kind=int1) :: i,k,l,ll,m

      real(kind=real2), parameter :: half_pi = 2d0*atan(1d0) 
      real(kind=real2) :: ak,an,d
      real(kind=real2) :: b(((n+1)/2)+1),bb
      real(kind=real2) :: coef,coef2,c
      real(kind=real2) :: tau((n + 1)/2)
      real(kind=real2) :: x1,xx,s,y

      logical :: even      
      
      m = (n + 1)/2
      even = .false.
      if(mod(n,2) .eq. 0) even = .true.
      
      d = 2d0
      an = 0d0
      do k=1,n
         an = an + 1d0
         d = d*an/(an + 0.5d0)
      end do

      ! Calculate Chebyshev coefficients of the orthogonal polynomial.
      tau(1) = (an + 2d0)/(2d0*an + 3d0)
      b(m) = tau(1) - 1d0
      ak = an
      
      do l=1,m-1         
         ak = ak + 2d0
         tau(l + 1) = ((ak - 1d0)*ak - an*(an + 1d0))*(ak + 2d0)*tau(l)/ &
            (ak*((ak + 3d0)*(ak + 2d0) - an*(an + 1d0)))
         b(m - l) = tau(l + 1)
         
         do ll=1,l
            b(m - l) = b(m - l) + tau(ll) * b(m - l + ll)
         enddo         
      enddo      
      b(m+1) = 1d0

      ! Calculate approximate values for the abscissas.
      bb = sin(half_pi/(2d0*an + 1d0))
      x1 = sqrt(1d0 - bb*bb)
      s = 2d0*bb*x1
      c = sqrt(1d0 - s*s)
      coef = 1d0 - (1d0 - 1d0/an)/(8d0*an*an)
      xx = coef*x1

      ! Coefficient needed for weights. COEF2 = 2^(2*n+1) * n! * n! / (2n+1)!
      coef2 = 2d0/dble(2*n + 1)
      do i=1,n
         coef2 = coef2*4d0*dble(i)/dble(n + i)
      enddo
      
      ! Calculate K-th abscissa (a Kronrod abscissa) and weight.
      do k=1,n,2
         call abwe1(n,m,niter,coef2,even,b,xx,w1(k))
         w2(k) = 0d0
         
         x(k) = xx
         y = x1
         x1 = y*c - bb*s
         bb = y*s + bb*c
         
         if(k .eq. n)then
            xx = 0d0
         else
            xx = coef*x1
         endif

         ! Calculate K+1 abscissa (a Gaussian abscissa) and  weights.
         call abwe2(n,m,niter,coef2,even,b,xx,w1(k+1),w2(k+1))
         
         x(k+1) = xx
         y = x1
         x1 = y*c - bb*s
         bb = y*s + bb*c
         xx = coef*x1
      enddo

      !  If N is even, compute extra Kronrod abscissa 
      if(even) then
         xx = 0d0
         call abwe1(n,m,niter,coef2,even,b,xx,w1(n+1))
         w2(n+1) = 0d0
         x(n+1) = xx
      endif
      
      return
   end subroutine kronrod
   !**********************************************************************
   !> @author John Burkardt
   !
   !> @brief adjusts Gauss-Kronrod rule from (-1,+1) to (a,b)
   !
   !> @param[in] a interval start
   !> @param[in] b interval end
   !> @param[in] n order of quadraute
   !> @param[in,out] x abscissas
   !> @param[in,out] w1 weight
   !> @param[in,out] w2 weight
   !**********************************************************************
   subroutine kronrod_adjust(a,b,n,x,w1,w2)
      use precision
      implicit none
      
      integer(kind=int1), intent(in) :: n
      real(kind=real2), intent(in) :: a,b

      real(kind=real2), intent(inout) :: x(n+1),w1(n+1),w2(n+1)
      
      x(:) = ((1d0 - x(:))*a + (1d0 + x(:))*b)/2d0
      
      w1(:) = ((b - a)/2d0)*w1(:)
      w2(:) = ((b - a)/2d0)*w2(:)

      return
   end subroutine kronrod_adjust
   !**********************************************************************
   !> @author Robert Piessens (F77), Maria Branders(F77) and John Burkardt (F90)
   !
   !> @brief Calculate Gauss-Kronrod abscissa and weight
   !
   !> @param[in] n order of the Gauss rule.
   !> @param[in] m (n + 1)/2
   !> @param[in] eps  accuracy of the abscissas
   !> @param[in] coef2 value needed to compute weights
   !> @param[in] even TRUE if n is even
   !> @param[in] b array of Chebyshev coefficients
   !> @param[in,out] x input is estimate of abscissa, output is computed abscissa
   !> @param[out] w weight
   !**********************************************************************
   subroutine abwe1(n,m,niter,coef2,even,b,x,w)
      use precision
      implicit none
      
      integer(kind=int1), intent(in) :: m,n,niter
      real(kind=real2), intent(in) :: coef2
      real(kind=real2), intent(in) :: b(m+1)
      logical, intent(in) :: even

      real(kind=real2), intent(out) :: w
      
      real(kind=real2), intent(inout) :: x

      real(kind=real2), parameter :: eps = 10d0*epsilon(b(1))
      
      integer(kind=int1) :: i,k,ka,iter

      real(kind=real2) :: ai
      real(kind=real2) :: c(3),d(3)
      real(kind=real2) :: delta,dif
      real(kind=real2) :: f,fd
      real(kind=real2) :: yy
      
      ka = 0
      if(x .eq. 0d0) ka = 1

      ! Iterative process for the computation of a Kronrod abscissa.
      do iter=1,niter

         c(2) = 0d0; c(3) = b(m+1)
         yy = 4d0*x*x - 2d0
         d(2) = 0d0
         
         if(even)then
            ai  = 2d0*m + 1d0
            d(3)  = ai*b(m+1)
            dif = 2d0
         else
            ai = m + 1d0
            d(3) = 0d0
            dif = 1d0
         endif
         
         do k=1,m
            ai = ai - dif
            i  = m - k + 1
            
            c = eoshift(c,1); d = eoshift(d,1)
            c(3) = yy*c(2) - c(1) + b(i)

            if(.not. even)then
               i = i + 1
            endif
            d(3) = yy*d(2) - d(1) + ai*b(i)
         enddo
         
         if(even)then
            f  = x*(c(3) - c(2))
            fd = d(3) + d(2)
         else
            f = 0.5d0*(c(3) - c(1))
            fd = 4d0*x*d(3)
         endif
         
         ! Newton correction.
         delta = f/fd
         x = x - delta  

         if(ka .eq. 1 .or. abs(delta) .le. eps)then
            ka = 1
            exit
         endif
      end do

      ! Catch non-convergence.
      if(ka .ne. 1)then
         print *,'ERROR: ABWE1 CONV FAIL'
         stop
      endif

      ! Computation of the weight.
      d(1) = 1d0; d(2) = x
      ai = 0d0
      do k=2,n
         ai = ai + 1d0
         d(3) = ((2d0*ai + 1d0)*x*d(2) - ai*d(1))/(ai + 1d0)
         d = eoshift(d,1)
      end do
      
      w = coef2/(fd*d(2))
      
      return
   end subroutine abwe1
   !**********************************************************************
   !> @author Robert Piessens (F77), Maria Branders(F77) and John Burkardt (F90)
   ! 
   !> @brief  Calculate Gaussian abscissa and two weights
   !
   !> @param[in] n order of the Gauss rule
   !> @param[in] m (N + 1)/2
   !> @param[in] eps accuracy of the abscissas
   !> @param[in] coef2 value needed to compute weights
   !> @param[in] even TRUE if N is even.
   !> @param[in] b(m+1) Chebyshev coefficients.
   !> @param[in,out] x input is estimate of abscissa, output is the computed abscissa
   !> @param[out] w1 Gauss-Kronrod weight
   !> @param[out] w2 Gauss weight.
   !**********************************************************************
   subroutine abwe2(n,m,niter,coef2,even,b,x,w1,w2)
      use precision
      implicit none
      
      integer(kind=int1), intent(in) :: n,m,niter
      real(kind=real2), intent(in) :: coef2
      real(kind=real2), intent(in) :: b(m+1)
      logical, intent(in) :: even
      
      real(kind=real2), intent(out) :: w1,w2

      real(kind=real2), intent(inout) :: x

      integer(kind=int1) :: i,k,ka,iter

      real(kind=real2), parameter :: eps = 10d0*epsilon(b(1))
      real(kind=real2) :: p(3),pd(3)
      real(kind=real2) :: ai,yy,delta
      
      k = 0d0
      if(x .eq. 0d0) ka = 1

      ! Iterative process for the computation of a Gaussian abscissa.
      do iter=1,niter      
         p(1) = 1d0
         p(2) = x
         pd(1) = 0d0
         pd(2) = 1d0

         ! When N is 1, we need to initialize P2 and PD2 to avoid problems with DELTA.
         if(n .le. 1)then
            if(eps .lt. abs (x))then
               p(3)  = (3d0*x*x - 1d0)/2d0
               pd(3) = 3d0*x
            else
               p(3)  = 3d0*x
               pd(3) = 3d0
            endif
         endif
         
         ai = 0d0
         do k=2,n
            ai  = ai + 1d0
            p(3)  = ((2d0*ai + 1d0)*x*p(2) - ai*p(1))/(ai + 1d0)
            pd(3) = ((2d0*ai + 1d0)*(p(2) + x*pd(2)) - ai*pd(1))/(ai + 1d0)

            p = eoshift(p,1)
            pd = eoshift(pd,1)
         end do

         ! Newton correction.
         delta = p(2)/pd(2)
         x = x - delta

         if(ka .eq. 1 .or. abs(delta) .le. eps)then
            ka = 1
            exit
         endif
      enddo

      ! Catch non-convergence.
      if(ka .ne. 1)then
         print *, 'ERROR: ABWE2 CONV FAIL'
         stop
      endif
      
      ! Computation of the weight.      
      w2 = 2d0/(n*pd(2)*p(1))
      
      p(2) = 0d0
      p(3) = b(m+1)
      yy = 4d0*x*x - 2d0
      do k=1,m
         i = m - k + 1
         
         p = eoshift(p,1)
         p(3) = yy*p(2) - p(1) + b(i)
      enddo
      
      if(even)then
         w1 = w2 + coef2/(pd(2)*x*(p(3) - p(2)))
      else
         w1 = w2 + 2d0*coef2/(pd(2)*(p(3) - p(1)))
      end if
      
      return
   end subroutine abwe2
   !********************************************************************** 
end module quadrature
