module polynomial
   use precision
   implicit none

   private

   public :: dlegendre
   public :: jacobi,jacobiD,jacobiDC,jacobi_int,jacobiProduct,jacobi_PDP
   public :: jacobi2poly,jacobi_swap,jacobi_swap_c
   public :: lagrange
   public :: legendre,legendreD,legendreDArray,legendreLambda,legendreProduct
   public :: legendreProductArray,legendreProductTrunc,legendre2poly
   public :: legendre_poly_weigth,legendreDPI,legendreDPI
   public :: polyval,poly2jacobi,poly2legendre,swap_indices
   public :: vandermonde_mono_row  ,vandermonde_jac
   public :: vandermonde_mono_row_d,vandermonde_jac_d
   public :: vandermonde_orth,vandermonde_orth_d
   
   interface jacobi
      module procedure jacobi_array
      module procedure jacobi_single
   end interface jacobi
   
contains
   !**********************************************************************
   double precision function dlegendre(x,dx,n) result(dP)  
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n
      real(kind=real2), intent(in) :: x,dx

      real(kind=real2) :: x1,x2

      x1 = x - 0.5d0*dx
      x2 = x + 0.5d0*dx
      dP = (legendre(x2,n) - legendre(x1,n))/dx

      return
   end function dlegendre
   !**********************************************************************
   function jacobi_single(a,b,n,x) result(J)
      use precision
      use maths, only : choose,fac => factorialReal
      implicit none

      integer(kind=int1), intent(in) :: n
      real(kind=real2), intent(in) :: a,b,x

      real(kind=real2) :: J
      
      integer(kind=int1) :: m
      real(kind=real2) :: c
      
      J = 0d0
      do m=0,n
         c = choose(n,m)*gamma(a + b + n + m + 1d0)/gamma(a + m + 1d0)
         J = J + c*((0.5d0*(x-1))**m)
      enddo
      J = J*gamma(a+n+1)/gamma(a+b+n+1)/fac(n)

      return
   end function jacobi_single
   !**********************************************************************
   function jacobi_array(a,b,n,x) result(J)
      use precision
      use maths, only : choose,fac => factorialReal
      implicit none

      integer(kind=int1), intent(in) :: n
      real(kind=real2), intent(in) :: a,b,x(:)

      real(kind=real2) :: J(size(x))
      
      integer(kind=int1) :: m
      real(kind=real2) :: c
      
      J = 0d0
      do m=0,n
         c = choose(n,m)*gamma(a + b + n + m + 1d0)/gamma(a + m + 1d0)
         J(:) = J(:) + c*((0.5d0*(x(:)-1))**m)
      enddo
      J = J*gamma(a+n+1)/gamma(a+b+n+1)/fac(n)

      return
   end function jacobi_array
   !**********************************************************************
   function jacobiD(q,n,a,b) result(d)
      use precision
      use maths, only : pochhammerR
      implicit none
      
      integer(kind=int1), intent(in) :: q,n
      real(kind=real2), intent(in) :: a,b

      real(kind=real2) :: d(n-q+1)

      integer(kind=int1) :: i
      real(kind=real2) :: g

      do i=0,n-q
         d(i+1) = jacobiDC(n,q,i,a,b)
      enddo
      d = d*(2d0**(-q))*pochhammerR(n+a+b+1,q)

      return
   end function jacobiD
   !**********************************************************************
   function jacobiDC(n,q,i,a,b) result(c)
      use precision
      use maths, only : pochhammerR,genHypGeo,factorial
      implicit none

      integer(kind=int1), intent(in) :: n,q,i

      real(kind=real2), intent(in) :: a,b

      real(kind=real2) :: c,a3(3),b2(2)

      c = pochhammerR(n+q+a+b+1d0,i)*pochhammerR(i+q+a+1d0,n-i-q)/&
          dble(factorial(n-i-q))/pochhammerR(i+a+b+1d0,i)

      a3 = (/dble(-n+q+i),dble(n+i+q+a+b+1),dble(i+a+1)/)
      b2 = (/dble(i+q+a+1),dble(2*i+a+b+2)/)
      c = c*genHypGeo(3,2,a3,b2,1d0)

      return
   end function jacobiDC
   !**********************************************************************
   function jacobi_int(m,n,a,b) result(qp)
      use precision
      use maths,only : facBig
      implicit none

      integer(kind=int1), intent(in) :: m,n
      real(kind=real2), intent(in) :: a,b

      real(kind=real2) :: qp

      if(m .ne. n) then
         qp = 0d0
         return
      else
         qp = (2d0**(a+b+1d0))/(2d0*n+a+b+1d0)
         qp = qp*gamma(n+a+1d0)*gamma(n+b+1d0)/facBig(n)/gamma(n+a+b+1d0)
      endif
      
      return
   end function jacobi_int
   !**********************************************************************
   function jacobiProduct(p,q,a,b) result(j)
      use precision
      use maths, only : invert
      implicit none

      integer(kind=int1), intent(in) :: p,q

      real(kind=real2), intent(in) :: a,b

      real(kind=real2) :: j(p+q+1)

      integer(kind=int1) :: i,io
      real(kind=real2), parameter :: pi = 4d0*atan(1d0)
      real(kind=real2) :: xc,Jp(1,p+q+1),X(p+q+1,p+q+1),temp(1,p+q+1),L(p+q+1,p+q+1)

      Jp = 0d0
      do i=1,p+q+1
         xc = cos(dble(2*i-1)*pi/dble(2*(p+q+1)))
         Jp(1,i) = jacobi(a,b,p,xc)*jacobi(a,b,q,xc)

         do io=0,p+q
            X(io+1,i) = jacobi(a,b,io,xc)
         enddo
      enddo

      call invert(p+q+1,X,L)

      temp = matmul(Jp,L)

      j(:) = temp(1,:)

      return
   end function jacobiProduct
   !**********************************************************************
   function jacobi_PDP(p,a,b) result(ap)
      use precision
      use maths, only : pochhammerR
      implicit none

      integer(kind=int1), intent(in) :: p
      real(kind=real2), intent(in) :: a,b

      real(kind=real2) :: ap

      ap = (2d0**(-p))*pochhammerR(p+a+b+1d0,p)
      
      return
   end function jacobi_PDP
   !**********************************************************************
   function jacobi2poly(order,a,b,J) result(Y)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: order
      real(kind=real2), intent(in) :: a,b
      real(kind=real2), intent(in) :: J(order+1)

      real(kind=real2) :: Y(order+1)

      integer(kind=int1) :: io
      real(kind=real2) :: L(order+1),S(order+1)

      L = 0d0
      do io=0,order
         S = 0d0
         S(1:io+1) = jacobi_swap(io,a,b,0d0,0d0)
         L(:) = L(:) + J(io+1)*S(:)
      enddo

      Y = legendre2poly(order,L)
      
      return
   end function jacobi2poly
   !**********************************************************************
   function jacobi_swap(n,a,b,g,d) result(J)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n
      real(kind=real2), intent(in) :: a,b,g,d

      real(kind=real2) :: J(n+1)

      integer(kind=int1) :: ik
      
      J = 0d0
      do ik=0,n
         J(ik+1) = jacobi_swap_c(n,ik,a,b,g,d)
      enddo
      
      return
   end function jacobi_swap
   !**********************************************************************
   function jacobi_swap_c(n,k,a,b,g,d) result(c)
      use precision
      use maths, only : pochhammerR,genHypGeo,facBig
      implicit none

      integer(kind=int1), intent(in) :: n,k
      real(kind=real2), intent(in) :: a,b,g,d

      real(kind=real2) :: c

      c = pochhammerR(k+1d0+a,n-k)*pochhammerR(n+1d0+a+b,k)*gamma(g+d+k+1d0)
      c = c/(facBig(n-k)*gamma(g+d+2d0*k+1d0))

      c = c*genHypGeo(3,2,[1d0*(-n+k),n+k+a+b+1d0,g+k+1d0],[a+k+1d0,g+d+2d0*k+2d0],1d0)
      
      return
   end function jacobi_swap_c
   !**********************************************************************
   function lagrange(nd,nb,np,ind,pb,w,x) result(l)
      use precision
      use maths, only : pseudo_inverse_R
      implicit none

      integer(kind=int1), intent(in) :: nd,nb,np
      integer(kind=int1), intent(in) :: ind(nd,nb)
      real(kind=real2), intent(in) :: pb(2)
      real(kind=real2), intent(in) :: w(nd,np),x(nd,np)

      real(kind=real2) :: l(np,nb)

      real(kind=real2) :: v(nb,np)
      
      v = vandermonde_jac(nd,nb,np,ind,pb(1),pb(2),x)
      l = pseudo_inverse_R(v)
            
      return
   end function lagrange
   !**********************************************************************
   function lagrange_orth(nd,nb,np,ind,pb,w,x) result(l)
      use precision
      use maths, only : pseudo_inverse_R
      implicit none

      integer(kind=int1), intent(in) :: nd,nb,np
      integer(kind=int1), intent(in) :: ind(nd,nb)
      real(kind=real2), intent(in) :: pb(2)
      real(kind=real2), intent(in) :: w(nd,np),x(nd,np)

      real(kind=real2) :: l(np,nb)

      integer(kind=int1) :: in,id
      real(kind=real2) :: v(nb,np),b(np,np),p(np,nb)
      
      v = vandermonde_orth(nd,nb,np,ind,pb(1),pb(2),w,x)
      p = pseudo_inverse_R(v)
      
      b = 0d0
      do in=1,np
         b(in,in) = 1d0
         do id=1,nd
            b(in,in) = b(in,in)*sqrt(w(id,in))
         enddo
      enddo      
      
      l = matmul(b,p)
      
      return
   end function lagrange_orth
   !**********************************************************************
   recursive double precision function legendre(x,n) result(P)  
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n
      real(kind=real2), intent(in) :: x
      
      if(n .eq. 0)then
         P = 1d0
      elseif(n .eq. 1)then
         P = x
      else
         P = ((2d0*dble(n) - 1d0)*x*legendre(x,n-1) - &
              (dble(n)-1d0)*legendre(x,n-2))/dble(n)
      endif 
      
      return
   end function legendre
   !**********************************************************************
   function legendreD(n,j) result(q)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n,j
      
      real(kind=real2) :: q
      
      if(mod(n,2) .eq. mod(j,2)) then
         q = 0d0
      elseif(j .gt. n)then
         q = 0d0
      else
         q = 2d0*j + 1d0
      endif

      return
   end function legendreD
   !**********************************************************************
   function legendreDArray(n) result(q)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n

      real(kind=real2) :: q(n+1)

      integer(kind=int1) :: j

      do j=1,n+1
         q(j) = legendreD(n,j-1)
      enddo

      return
   end function legendreDArray
   !**********************************************************************
   function legendreLambda(n) result(l)
      use precision
      use maths, only : factorial
      implicit none

      integer(kind=int1), intent(in) :: n

      real(kind=real2) :: l

      l = dble(factorial(2*n))/dble(((2**n)*(factorial(n)**2)))

      return
   end function legendreLambda
   !**********************************************************************
   function legendreProduct(p,q) result(b)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: p,q

      real(kind=real2) :: b(p+q+1)

      integer(kind=int1) :: kmax,i,k

      b = 0d0
      kmax = floor(dble(p+q)/2d0)
      do k=0,kmax
         i = p + q - 2*k + 1
         b(i) = legendreProductA(p,q,k)
      enddo

      return
   end function legendreProduct
   !**********************************************************************
   function legendreProductA(p,q,k) result(a)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: p,q,k

      real(kind=real2) :: a

      a = dble(2*p + 2*q - 4*k + 1)/dble(2*p + 2*q - 2*k + 1)
      a = a*legendreLambda(k)*legendreLambda(p-k)*legendreLambda(q-k)/ &
          legendreLambda(p+q-k)

      return
   end function legendreProductA
   !**********************************************************************
   !> @brief Coeffiectent Matrix for U(1:p)^2 truncated to order n
   !**********************************************************************
   function legendreProductArray(p,n) result(b)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: p,n

      real(kind=real2) :: b(p*p,n)
      
      integer(kind=int1) :: mb,i,j

      do j=1,p
         do i=1,p
            mb = (j-1)*p + i
            b(mb,:) = legendreProductTrunc(n,i-1,j-1)
         enddo
      enddo

      return
   end function legendreProductArray
   !**********************************************************************
   function legendreProductTrunc(n,p,q) result(bn)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n,p,q

      real(kind=real2) :: bn(n)

      integer(kind=int1) :: kmax,i,k

      bn = 0d0
      kmax = floor(dble(p+q)/2d0)
      do k=0,kmax
         i = p + q - 2*k + 1
         if(i .gt. n) cycle

         bn(i) = legendreProductA(p,q,k)
      enddo

      return
   end function legendreProductTrunc
   !**********************************************************************
   function legendre2poly(p,h) result(l)
      use precision
      use maths, only : invert
      implicit none

      integer(kind=int1), intent(in) :: p

      real(kind=real2), intent(in) :: h(p+1)

      real(kind=real2) :: l(p+1)

      integer(kind=int1) :: io,il
      real(kind=real2) :: Lx(p+1,p+1),Lp(p+1,p+1)
      
      Lx = 0d0
      
      do io=0,p
         do il=io,0,-2
            Lx(io+1,il+1) = legendre_poly_weigth(io,il)
         enddo
      enddo
      call invert(p+1,Lx,Lp)

      do io=1,p+1
         Lp(io,:) = h(io)*Lp(io,:)
      enddo

      l = sum(Lp,1)
      
      return
   end function legendre2poly
   !**********************************************************************
   double precision function legendre_poly_weigth(n,l) result(p)
      use precision
      use maths, only : factorial,factorial2
      implicit none

      integer(kind=int1), intent(in) :: n,l

      integer(kind=int1) :: nl05
     
      nl05 = (n-l)/2
     
      p = (2d0*l + 1d0)*factorial(n)/(1d0*(2**nl05)*factorial(nl05)*factorial2(l + n + 1))
     
      return
   end function legendre_poly_weigth
   !**********************************************************************
   ! calculates \int^1_{-1} d^n1(\psi_m1)/dx^n1 * d^n2(\psi_m2)/dx^n2 dx
   ! where \psi_i is the i^th Legendre polynomial of the first kind
   !**********************************************************************
   function legendreDPI(n1,m1,n2,m2) result(i)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n1,m1,n2,m2

      real(kind=real2) :: i

      integer(kind=int1) :: n,m,in,jn,imax,jmax
      
      n = n1 + n2
      m = m1 + m2
    
      imax = floor(0.5d0*(n1-m1))
      jmax = floor(0.5d0*(n2-m2))
    
      i = 0d0
      do in=0,imax
         do jn=0,jmax
            i = i + legendreDPIc(in,m1,n1)*legendreDPIc(jn,m2,n2)* &
                 (1d0-(-1d0)**(N-M-2*(in+jn)+1))/dble(N-M-2*(in+jn)+1)
         enddo
      enddo
    
      return
   end function legendreDPI
   !**********************************************************************
   function legendreDPIc(i,m,n) result(c)
      use precision
      use maths, only : factorial2,factorial
      implicit none

      integer(kind=int1), intent(in) :: i,m,n
      
      real(kind=real2) :: c
      
      c = ((-1d0)**i)*factorial(2*(n-i))/ &
           dble((2d0**n)*factorial(n-m-2*i)*factorial(n-i)*factorial(i))

      return
   end function legendreDPIc
   !**********************************************************************
   !> @breif Evaluated the polynomial defined by p at x
   !> @par this differ slightly from matlab in that it calculates:
   !! y(x) = p(1)*x^0 + p(2)*x^1 + ....
   !**********************************************************************
   function polyval(p,x) result(y)
      use precision
      implicit none

      real(kind=real2), intent(in) :: p(:),x(:)

      real(kind=real2) :: y(size(x))

      integer(kind=int1) :: i

      y = 0d0
      do i=1,size(p)
         y(:) = y(:) + p(i)*(x(:)**(i-1))
      enddo

      return
   end function polyval
   !**********************************************************************
   function poly2jacobi(order,a,b,Y) result(J)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: order
      real(kind=real2), intent(in) :: a,b
      real(kind=real2), intent(in) :: Y(order+1)

      real(kind=real2) :: J(order+1)

      integer(kind=int1) :: io
      real(kind=real2) :: L(order+1),S(order+1)

      J = 0d0
      
      L = poly2legendre(order,Y)

      do io=0,order
         S = 0d0
         S(1:io+1) = jacobi_swap(io,0d0,0d0,a,b)
         
         J(:) = J(:) + L(io+1)*S(:)
      enddo

      
      return
   end function poly2jacobi
   !**********************************************************************
   function poly2legendre(n,p) result(l)
      use precision
      use maths, only : invert
      implicit none

      integer(kind=int1), intent(in) :: n

      real(kind=real2), intent(in) :: p(n+1)

      real(kind=real2) :: l(n+1)

      integer(kind=int1) :: io,il
      real(kind=real2) :: Lx(n+1,n+1),Lp(n+1,n+1)
      
      Lx = 0d0
      
      do io=0,n
         do il=io,0,-2
            Lx(io+1,il+1) = legendre_poly_weigth(io,il)
         enddo
      enddo

      do io=1,n+1
         Lp(io,:) = p(io)*Lx(io,:)
      enddo

      l = sum(Lp,1)
      
      return
   end function poly2legendre
   !**********************************************************************
   function swap_indices(nd,nb1,ind1,nb2,ind2,q1) result(q2)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: nd,nb1,nb2
      integer(kind=int1), intent(in) :: ind1(nd,nb1),ind2(nd,nb2)

      real(kind=real2), intent(in) :: q1(nb1)

      real(kind=real2) :: q2(nb2)

      integer(kind=int1) :: i1,i2

      q2 = 0d0
      do i1=1,nb1
         do i2=1,nb2
            if(all(ind1(:,i1) .eq. ind2(:,i2))) q2(i2) = q1(i1)
         enddo
      enddo

   end function swap_indices
   !**********************************************************************
   function vandermonde_jac(ndims,nb,np,ind,alph,beta,x) result(v)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: ndims,nb,np
      integer(kind=int1), intent(in) :: ind(ndims,nb)
      
      real(kind=real2), intent(in) :: alph,beta 
      real(kind=real2), intent(in) :: x(ndims,np)

      real(kind=real2) :: v(nb,np)

      integer(kind=int1) :: ib,id
      
      v = 1d0
      do id=1,ndims
         do ib=1,nb
            v(ib,1:np) = v(ib,1:np)*jacobi(alph,beta,ind(id,ib),x(id,1:np))
         enddo
      enddo
           
      return
   end function vandermonde_jac
   !**********************************************************************
   function vandermonde_orth(ndims,nb,np,ind,alph,beta,w,x) result(v)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: ndims,nb,np
      integer(kind=int1), intent(in) :: ind(ndims,nb)
      
      real(kind=real2), optional, intent(in) :: alph,beta 
      real(kind=real2), intent(in) :: w(ndims,np)
      real(kind=real2), intent(in) :: x(ndims,np)

      real(kind=real2) :: v(nb,np)

      integer(kind=int1) :: ib,id,in
      real(kind=real2) :: a,b,j(np)

      a = 0d0
      b = 0d0
      if(present(alph)) a = alph
      if(present(beta)) b = beta
      
      v = 1d0
      
      do id=1,ndims
         do ib=1,nb
            j = jacobi(a,b,ind(id,ib),x(id,1:np))
            do in=1,np
               v(ib,in) = v(ib,in)*j(in)*sqrt(w(id,in))
            enddo
         enddo
      enddo
           
      return
   end function vandermonde_orth
   !**********************************************************************
   function vandermonde_mono_row(ndims,nb,ind,x) result(vr)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: ndims,nb
      integer(kind=int1), intent(in) :: ind(ndims,nb)

      real(kind=real2), intent(in) :: x(ndims)

      real(kind=real2) :: vr(nb)

      integer(kind=int1) :: ib,id
      
      do ib=1,nb
         vr(ib) = 1d0
         do id=1,ndims
            vr(ib) = vr(ib)*(x(id)**dble(ind(id,ib)))
         enddo
      enddo
      
      return
   end function vandermonde_mono_row
   !**********************************************************************
   function vandermonde_mono_row_d(ndims,d,nb,ind,x) result(vr)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: ndims,d,nb
      integer(kind=int1), intent(in) :: ind(ndims,nb)

      real(kind=real2), intent(in) :: x(ndims)

      real(kind=real2) :: vr(nb)

      integer(kind=int1) :: ib,id,fac
      
      do ib=1,nb
         vr(ib) = 1d0
         do id=1,ndims
            if(id .eq. d) then
               fac = max(ind(id,ib)-1,0)
               vr(ib) = vr(ib)*ind(id,ib)*(x(id)**dble(fac))
            else
               vr(ib) = vr(ib)*(x(id)**dble(ind(id,ib)))
            endif
         enddo
      enddo
      
      return
   end function vandermonde_mono_row_d
   !**********************************************************************
   function vandermonde_jac_d(ndims,d,nb,np,order,ind,alph,beta,x) result(v)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: ndims,d,nb,np,order
      integer(kind=int1), intent(in) :: ind(ndims,nb)

      real(kind=real2), intent(in) :: alph,beta
      real(kind=real2), intent(in) :: x(ndims,np)

      real(kind=real2) :: v(nb,np)

      integer(kind=int1) :: ib,id,ij,p
      real(kind=real2) :: dj(order),djs(np)
      
      v = 1d0
      do id=1,ndims
         do ib=1,nb
            p = ind(id,ib)
            
            if(id .eq. d) then
               dj(1:p) = jacobiD(1,p,alph,beta)
               djs = 0d0
               do ij=1,p
                  djs(1:np) = djs(1:np) + dj(ij)*jacobi(alph,beta,ij-1,x(id,1:np))
               enddo
               v(ib,1:np) = v(ib,1:np)*djs(1:np)
            else
               v(ib,1:np) = v(ib,1:np)*jacobi(alph,beta,p,x(id,1:np))
            endif
            
         enddo
      enddo
      
      return
   end function vandermonde_jac_d
   !**********************************************************************
   function vandermonde_orth_d(ndims,d,nb,np,order,ind,alph,beta,w,x) result(v)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: ndims,d,nb,np,order
      integer(kind=int1), intent(in) :: ind(ndims,nb)

      real(kind=real2), optional, intent(in) :: alph,beta
      real(kind=real2), intent(in) :: w(ndims,np)
      real(kind=real2), intent(in) :: x(ndims,np)

      real(kind=real2) :: v(ib,np)

      integer(kind=int1) :: ib,id,ij,p
      real(kind=real2) :: dj(order),djs(np),a,b

      a = 0d0
      b = 0d0
      if(present(alph)) a = alph
      if(present(beta)) b = beta
      
      v = 1d0
      do id=1,ndims
         do ib=1,nb
            p = ind(id,ib)
            
            if(id .eq. d) then
               dj(1:p) = jacobiD(1,p,a,b)
               djs = 0d0
               do ij=1,p
                  djs(1:np) = djs(1:np) + dj(ij)*jacobi(a,b,ij-1,x(id,1:np))
               enddo
               v(ib,1:np) = v(ib,1:np)*djs(1:np)*sqrt(w(id,1:np))
            else
               v(ib,1:np) = v(ib,1:np)*jacobi(a,b,p,x(id,1:np))*sqrt(w(id,1:np))
            endif
            
         enddo
      enddo


   end function vandermonde_orth_d
   !**********************************************************************
end module polynomial
