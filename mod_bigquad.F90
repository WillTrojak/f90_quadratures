module bigquad
   use precision
   implicit none

   private

   real(kind=real2), save :: faca(31)
   real(kind=real2), parameter :: tol = epsilon(faca(1))
   
   public :: gauss_Legendre,gaussl_nodes,gaussl_weights,stieltjes_c
   public :: baratella_dleg_approx
  
contains
   !**********************************************************************
   !> @breif Calculates Gauss-Legendre nodes for Large N
   !> @par This calculates the Gauss-Legendre nodes for large n based on
   !! Townsend & Hale (2013) and Bogeart et al. (2013). Recomended for N < 100.
   !! Although Golub-Welsch is still acceptable until N~=1000
   !> @param[in] n order (number of points)
   !> @param[out] x nodes, preallocated for speed
   !> @param[out] t theta of nodes such that x = cos(t), preallocated for speed
   !**********************************************************************
   subroutine gaussl_nodes(n,x,t)
      use precision
      use maths, only : pi,factorialArrayReal
      implicit none

      integer(kind=int1), intent(in) :: n 

      real(kind=real2), intent(out) :: x(n),t(n)

      integer(kind=int1) :: n2,k
      
      n2 = floor(real(n,kind=real2)/2d0)

      faca = factorialArrayReal(30)
      
      do k=1,n2
         t(k) = theta0(n,k)
      enddo

      call gaussl_node_adjust(n,n2,t)

      do k=1,n2
         t(n-k+1) = pi - t(k)
         
         x(k) = cos(t(k))
         x(n-k+1) = -cos(t(k))
      enddo

      if(mod(n,2) .ne. 0) x(n2 + 1) = 0d0
      
      return
   end subroutine gaussl_nodes
   !**********************************************************************
   !> @breif Calculates Gauss-Legendre nodes for Large N
   !> @par This calculates the Gauss-Legendre nodes for large n based on
   !! Townsend & Hale (2013) and Bogeart et al. (2013). Recomended for N < 100.
   !! Although Golub-Welsch is still acceptable until N~=1000
   !> @param[in] n order (number of points)
   !> @param[out] x nodes, preallocated for speed
   !> @param[out] t theta of nodes such that x = cos(t), preallocated for speed
   !**********************************************************************
   subroutine gaussl_weights(n,w)
      use precision
      use maths, only : pi,factorialArrayReal
      implicit none

      integer(kind=int1), intent(in) :: n 

      real(kind=real2), intent(out) :: w(n)

      integer(kind=int1) :: n2,k

      real(kind=real2) :: t(n),dpn
      
      n2 = floor(real(n,kind=real2)/2d0)

      faca = factorialArrayReal(20)
      
      do k=1,n2
         t(k) = theta0(n,k)
      enddo

      call gaussl_node_adjust(n,n2,t,w)

      do k=1,n2
         w(n-k+1) = w(k)
      enddo

      if(mod(n,2) .ne. 0) then
         dpn = baratella_dleg_approx(n,pi/2d0)
         w(n2 + 1) = 2d0/(dpn*dpn)
      endif
         
      return
   end subroutine gaussl_weights
   !**********************************************************************
   !> @breif Calculates Gauss-Legendre quadrature for Large N
   !> @par This calculates the Gauss-Legendre nodes adn weights for large n based on
   !! Townsend & Hale (2013) and Bogeart et al. (2013). Recomended for N < 100.
   !! Although Golub-Welsch is still acceptable until N~=1000
   !> @param[in] n order (number of points)
   !> @param[out] x nodes
   !> @param[out] w weigths
   !**********************************************************************
   subroutine gauss_Legendre(n,x,w)
      use precision
      use maths, only : pi,factorialArrayReal
      implicit none

      integer(kind=int1), intent(in) :: n 

      real(kind=real2), intent(out) :: x(n),w(n)

      integer(kind=int1) :: n2,k

      real(kind=real2) :: t(n),dpn
      
      n2 = floor(real(n,kind=real2)/2d0)

      faca = factorialArrayReal(20)
      
      do k=1,n2
         t(k) = theta0(n,k)
      enddo

      call gaussl_node_adjust(n,n2,t,w)

      do k=1,n2
         x(k) = cos(t(k))
         x(n-k+1) = -cos(t(k))
         
         w(n-k+1) = w(k)
      enddo

      if(mod(n,2) .ne. 0) then
         x(n2 + 1) = 0d0
         
         dpn = baratella_dleg_approx(n,pi/2d0)
         w(n2 + 1) = 2d0/(dpn*dpn)
      endif
         
      return
   end subroutine gauss_Legendre
   !**********************************************************************
   !> @brief Adjust inital guess of node
   !> @par Adjust the initial guess of thenodes location using Newtons method
   !! together with suitable methods for estimating the value of the legendre
   !! polynomial and its gradient.
   !> @param[in] n order
   !> @param[in] n2 floor(order/2)
   !> @param[in] t0 initial guess of node theta
   !> @param[out] t refined vaule of theta
   !**********************************************************************
   subroutine gaussl_node_adjust(n,n2,t,w)
      use precision
      use maths, only : pi
      implicit none

      integer(kind=int1), intent(in) :: n,n2

      real(kind=real2), intent(inout) :: t(n)

      real(kind=real2), optional, intent(out) :: w(n) 

      integer(kind=int1) :: k,lim
      
      real(kind=real2) :: cn
      
      lim = max(10,floor(n/2000d0))
      lim = 10
      
      cn = stieltjes_c(n)

      if(present(w))then
         do k=1,n2
            if(k .gt. lim) then
               call gaussl_newton_central(n,k,cn,t(k),w(k))
            else
               call gaussl_newton_boundary(n,k,t(k),w(k))
            endif
         enddo
      else
         do k=1,n2
            if(k .gt. lim) then
               call gaussl_newton_central(n,k,cn,t(k))
            else
               call gaussl_newton_boundary(n,k,t(k))
            endif
         enddo
      endif
      
      return
   end subroutine gaussl_node_adjust
   !**********************************************************************
   subroutine gaussl_newton_central(n,k,cn,t,w)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n,k

      real(kind=real2), intent(in) :: cn
      
      real(kind=real2), intent(inout) :: t

      real(kind=real2), optional, intent(out) :: w

      integer(kind=int1) :: j
      integer(kind=int1), parameter :: jmax = 20

      real(kind=real2) :: delt,pn,dpn,ttol
      
      delt = 1d0

      j = 0
      ttol = t*tol
      do while(abs(delt) .gt. ttol)
         j = j + 1
         call stieltjes_leg_approx(n,k,cn,t,pn,dpn)
         delt = pn/(dpn)

         t = t + delt
         
         if(j .gt. jmax)then
            print *,'ERROR: BIGQ CONV FAIL, CENTRAL',pn,k
            exit
         endif
      enddo

      ! For more accurate weights uncomment this line
      ! WARNING, this will be much slower
      !dpn = baratella_dleg_approx(n,t)
      if(present(w)) w = 2d0/(dpn*dpn)
      
      return
   end subroutine gaussl_newton_central
   !**********************************************************************
   subroutine gaussl_newton_boundary(n,k,t,w)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n,k

      real(kind=real2), intent(inout) :: t

      real(kind=real2), optional, intent(out) :: w

      integer(kind=int1), parameter :: jmax = 20
      integer(kind=int1) :: j
      
      real(kind=real2) :: delt,pn,dpn,ttol

      delt = 1d0
      j = 0
      ttol = t*tol
      do while(abs(delt) .gt. ttol)
         j = j + 1
         pn  = bogaert_leg_approx(n,t)
         dpn = baratella_dleg_approx(n,t)
         delt = pn/(dpn)
         
         t = t + delt
         
         if(j .gt. jmax)then
            print *,'ERROR: BIGQ CONV FAIL, BOUNDARY',pn,k
            exit
         endif
      enddo
      
      ! For more accurate weights uncomment this line
      ! WARNING, this will be much slower
      !dpn = baratella_dleg_approx(n,t)
      if(present(w)) w = 2d0/(dpn*dpn)
      
      return
   end subroutine gaussl_newton_boundary
   !**********************************************************************
   subroutine stieltjes_leg_approx(n,k,cn,t,pn,dpn) ! best for |x|< 0.866 (=cos(pi/6))
      use precision
      use maths, only : pi
      implicit none

      integer(kind=int1), intent(in) :: n,k
      real(kind=real2), intent(in) :: cn,t

      real(kind=real2), intent(out) :: pn,dpn

      integer(kind=int1) :: nm,m

      real(kind=real2), parameter :: pisixth = pi/6d0
      real(kind=real2) :: rn,rm,beta,del,hnm,z
      real(kind=real2) :: sina,cosa,sinz,cosz,cott,r2sin,r2sinm
      
      rn = real(n,kind=real2)

      z = (rn + 0.5d0)*t - 0.25d0*pi
      del = z - real(k,kind=real2)*pi
      
      sinz = taylor_sin_kpi(k,del)
      cosz = taylor_cos_kpi(k,del)
      
      cott = 1d0/tan(t)

      r2sin  = 1d0/(2d0*sin(t))
      r2sinm = sqrt(r2sin)
      
      pn  = 0d0
      dpn = 0d0

      nm = 15
      if(t .gt. pisixth) nm = 7
      
      do m=0,nm-1
         rm = real(m,kind=real2)

         beta = rm*(t - 0.5d0*pi)

         sina = sinz*cos(beta) + cosz*sin(beta)
         cosa = cosz*cos(beta) - sinz*sin(beta)
         
         hnm = stieltjes_h(n,m)
         pn  = pn  + hnm*cosa*r2sinm
         dpn = dpn + hnm*((rm - 0.5d0)*cosa*cott + (rn + rm - 0.5d0)*sina)*r2sinm
         
         r2sinm = r2sinm*r2sin
      enddo

      pn  = cn*pn
      dpn = cn*dpn
      
      return
   end subroutine stieltjes_leg_approx
   !**********************************************************************
   function taylor_sin_kpi(k,del) result(z)
      use precision
      use maths, only : pi
      implicit none

      integer(kind=int1), intent(in) :: k
      real(kind=real2), intent(in) :: del

      integer(kind=int1), parameter :: nm = 10
      integer(kind=int1) :: i
      real(kind=real2) :: z,c,m1,delz

      m1 = 1d0
      if(mod(k,2) .eq. 0) m1 = -1d0
            
      c = 1d0/del
      z = sin(real(k)*pi)
      do i=1,nm
         m1 = -m1
         c = c*del*del
         
         delz = m1*c/faca(2*i)

         z = z + delz
      enddo      
      
      return
   end function taylor_sin_kpi
   !**********************************************************************
   function taylor_cos_kpi(k,del) result(z)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: k
      real(kind=real2), intent(in) :: del

      integer(kind=int1), parameter :: nm = 10
      integer(kind=int1) :: i
      real(kind=real2) :: z,c,m1,delz

      m1 = -1d0
      if(mod(k,2) .eq. 0) m1 = 1d0

      c = 1d0
      z = m1
      do i=1,nm
         m1 = -m1
         c = c*del*del
         
         delz = m1*c/faca(2*i+1)
         
         z = z + delz
      enddo      
      
      return
   end function taylor_cos_kpi
   !**********************************************************************
   function stieltjes_c(n) result(cn)
      use precision
      use maths, only : ex,pi
      implicit none

      integer(kind=int1), intent(in) :: n

      real(kind=real2) :: cn

      real(kind=real2) :: rn
      
      rn = real(n,kind=real2) 

      cn = quotentent_large_indice(n,0.5d0)
      cn = cn*sqrt(rn/(rn + 0.5d0))*sqrt(ex/rn)
      cn = cn*stirling_series(rn)/stirling_series(rn + 0.5d0)      

      cn = cn*sqrt(4d0/pi)
      
      return
   end function stieltjes_c
   !**********************************************************************
   function stirling_series(x) result(s)
      use precision
      implicit none

      real(kind=real2), intent(in) :: x

      real(kind=real2) :: s

      ! First 10 terms in taylor expansion for stirlings series
      ! cn is normally only calculated onceper quad so 10 is not expensive
      real(kind=real2), parameter :: si(10) = [ &
            8.3333333333333329E-002, &
            3.4722222222222220E-003, &
           -2.6813271604938273E-003, &
           -2.2947209362139917E-004, &
            7.8403922172006662E-004, &
            6.9728137583658571E-005, &
           -5.9216643735369393E-004, &
           -5.1717909082605919E-005, &
            8.3949872067208726E-004, &
            7.2048954160200109E-005]
      
      integer(kind=int1) :: i
      real(kind=real2) :: rx,rxn

      rx = 1d0/x
      rxn = rx
      
      s = 1d0
      do i=1,10
         s = s + si(i)*rxn
         rxn = rx*rxn
      enddo
      
      return
   end function stirling_series
   !**********************************************************************
   function quotentent_large_indice(n,a) result(z)
      ! calculates (n/(n+a))**(n+a) for large values of n
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n
      real(kind=real2), intent(in) :: a

      real(kind=real2) :: z

      integer(kind=int1) :: j
      real(kind=real2) :: s,rn,ds

      s = 0d0
      rn = real(n,kind=real2)
      ds = 1d0
      
      j = 0
      do while(abs(ds) .gt. tol)
         j = j + 1
         ds = ((-a/rn)**real(j,kind=real2))/real(j,kind=real2)
         s = s + ds
      enddo

      z = (rn + a)*s
      z = exp(z)
      
      return
   end function quotentent_large_indice
   !**********************************************************************
   function stieltjes_h(n,m) result(h)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n,m

      real(kind=real2) :: h

      integer(kind=int1) :: j
      real(kind=real2) :: rn,rj
      
      h = 1d0
      
      if(m .gt. 0)then
         rn = real(n,kind=real2)
         do j=1,m
            rj = real(j,kind=real2)
            
            h = h*(rj - 0.5d0)*(rj - 0.5d0)
            h = h/(rj*(rn + rj + 0.5d0))
         enddo
      endif
      
      return
   end function stieltjes_h
   !**********************************************************************
   function bogaert_leg_approx(n,t) result(pn) ! best for cos(pi/6) <= |x| <= 1 as expensive
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n
      real(kind=real2), intent(in) :: t

      real(kind=real2) :: pn

      integer(kind=int1) :: i
      real(kind=real2) :: v,y,fn(7),c,rv2

      v = real(n,kind=real2) + 0.5d0
      rv2 = 1/v/v
      y = t*v

      fn = bogaert_f(y)

      c = 1d0
      pn = 0d0
      do i=0,6
         pn = pn + fn(i+1)*c
         c = c*rv2
      enddo
      
      return
   end function bogaert_leg_approx
   !**********************************************************************
   function bogaert_f(y) result(fn)
      use precision
      implicit none

      real(kind=real2), intent(in) :: y

      real(kind=real2) :: fn(7)

      ! Bessel function approxiamtion parameters from Bogaert, Michiels and Fostier
      real(kind=real2), parameter :: c21 =  1d0/8d0
      real(kind=real2), parameter :: c22 = -1d0/12d0
      real(kind=real2), parameter :: c42 =  11d0/384d0 
      real(kind=real2), parameter :: c43 = -7d0/160d0
      real(kind=real2), parameter :: c44 =  1d0/160d0
      real(kind=real2), parameter :: c63 =  173d0/15360d0
      real(kind=real2), parameter :: c64 = -101d0/3584d0 
      real(kind=real2), parameter :: c65 =  671d0/80640d0
      real(kind=real2), parameter :: c66 = -61d0/120960d0
      real(kind=real2), parameter :: c84 =  22931d0/3440640d0
      real(kind=real2), parameter :: c85 = -90497d0/3870720d0
      real(kind=real2), parameter :: c86 =  217d0/20480d0
      real(kind=real2), parameter :: c87 = -1261d0/967680d0
      real(kind=real2), parameter :: c88 =  1261d0/29030400d0
      real(kind=real2), parameter :: cX5 =  1319183d0/247726080d0
      real(kind=real2), parameter :: cX6 = -10918993/454164480d0
      real(kind=real2), parameter :: cX7 =  1676287d0/113541120d0
      real(kind=real2), parameter :: cX8 = -7034857d0/2554675200d0
      real(kind=real2), parameter :: cX9 =  1501d0/8110080d0
      real(kind=real2), parameter :: cXX = -79d0/20275200d0
      real(kind=real2), parameter :: cZ6 =  233526463d0/43599790080d0
      real(kind=real2), parameter :: cZ7 = -1396004969d0/47233105920d0
      real(kind=real2), parameter :: cZ8 =  2323237523d0/101213798400d0
      real(kind=real2), parameter :: cZ9 = -72836747d0/12651724800d0
      real(kind=real2), parameter :: cZX =  3135577d0/5367398400d0
      real(kind=real2), parameter :: cZY = -1532789d0/61993451520d0
      real(kind=real2), parameter :: cZZ =  66643d0/185980354560d0

      integer(kind=int1) :: i
      real(kind=real2) :: h(13),c

      c = 1d0
      do i=1,13
         h(i) = bessel_jn(i-1,y)*c
         c = c*y
      enddo
      
      fn(1) = h(1)
      fn(2) = c21*h(2) + c22*h(3)
      fn(3) = c42*h(3) + c43*h(4) + c44*h(5)
      fn(4) = c63*h(4) + c64*h(5) + c65*h(6) + c66*h(7)
      fn(5) = c84*h(5) + c85*h(6) + c86*h(7) + c87*h(8) + c88*h(9)
      fn(6) = cX5*h(6) + cX6*h(7) + cX7*h(8) + cX8*h(9) + cX9*h(10) + cXX*h(11)
      fn(7) = cZ6*h(7) + cZ7*h(8) + cZ8*h(9) + cZ9*h(10) + cZX*h(11) + &
              cZY*h(12) + cZZ*h(13)
      
      return
   end function bogaert_f
   !**********************************************************************
   function baratella_dleg_approx(n,t) result(dpn) ! best for cos(pi/6) <= |x| <=1
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n

      real(kind=real2), intent(in) :: t

      real(kind=real2) :: dpn

      integer(kind=int1), parameter :: lmax = 15
      integer(kind=int1) :: l
      
      real(kind=real2) :: delj,rt,mt,mtl,rfacl

      rt = (real(n,kind=real2)+ 0.5d0)*t
      rfacl = 1d0
      
      mt  = -t
      mtl = 1d0
      
      dpn = 0d0
      l = 0
      delj = 1d0

      do l=1,lmax
         rfacl = rfacl/real(l,kind=real2)

         delj = bessel0d(l,rt)*mtl*rfacl
        
         dpn = dpn + delj
         
         mtl = mtl*mt
      enddo

      dpn = (t*dpn + (cos(t) - 1d0)*bessel_jn(0,rt))/sin(t)
      dpn = dpn*sqrt(t/sin(t))
      dpn = -real(n,kind=real2)*dpn

      
      return
   end function baratella_dleg_approx
   !**********************************************************************
   function bessel0d(l,z) result(jl)
      use precision 
      implicit none

      integer(kind=int1), intent(in) :: l

      real(kind=real2), intent(in) :: z

      real(kind=real2) :: jl

      integer(kind=int1) :: j
      real(kind=real2) :: m1,c

      jl =  0d0
      m1 = -1d0
      c = 1d0
      do j=0,l
         m1 = -m1
         
         jl = jl + bessel_jn((-l + 2*j),z)*m1*c
         c = c*real(l-j,kind=real2)/real(j+1,kind=real2)
      enddo
      jl = jl*(2d0**(-l))

      return
   end function bessel0d
   !**********************************************************************
   function theta0(n,k) result(tk)
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n,k

      real(kind=real2) :: tk
      
      real(kind=real2) :: xk

      xk = tricomi_approx(n,k)
      if(abs(xk) .gt. 0.5d0) xk = olver_approx(n,k)
         
      tk = acos(xk)      
      
      return
   end function theta0
   !**********************************************************************
   function tricomi_approx(n,k) result(xk) ! best for |xk| < 0.5
      use precision
      use maths, only : pi
      implicit none

      integer(kind=int1), intent(in) :: n,k

      real(kind=real2) :: xk

      real(kind=real2) :: phi,rn

      rn = real(n,kind=real2)
      phi = (real(k,kind=real2) - 0.25d0)*pi/(rn + 0.5d0)

      xk = 39d0 - 28d0/(sin(phi)**2d0)
      xk = -xk/(384d0*rn**4d0)
      xk = xk - (rn - 1d0)/(8d0*rn**3d0)
      xk = 1d0 + xk
      xk = xk*cos(phi)
      
      return
   end function tricomi_approx
   !**********************************************************************
   function olver_approx(n,k) result(xk) ! best for 0.5 <= |xk| <=1
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n,k

      real(kind=real2) :: xk

      real(kind=real2) :: psi,nph

      nph = real(n,kind=real2) + 0.5d0
      psi = bessel0_root(k)/nph

      xk = psi/tan(psi) - 1d0
      xk = xk/(8d0*psi*nph*nph)
      xk = xk + psi
      xk = cos(xk)
      
      return
   end function olver_approx
   !**********************************************************************
   function bessel0_root(k) result(j0)
      use precision
      use maths, only : pi
      implicit none

      integer(kind=int1), intent(in) :: k

      real(kind=real2) :: j0

      ! Approximation parameters from JCP 42, pp. 403-405 (1981) 
      real(kind=real2), parameter :: a0 = 0.682894897349453d-1
      real(kind=real2), parameter :: a1 = 0.131420807470708d+0
      real(kind=real2), parameter :: a2 = 0.245988241803681d-1
      real(kind=real2), parameter :: a3 = 0.813005721543268d-3

      real(kind=real2), parameter :: b0 = 1d0
      real(kind=real2), parameter :: b1 = 0.116837242570470d+1
      real(kind=real2), parameter :: b2 = 0.200991122197811d+0
      real(kind=real2), parameter :: b3 = 0.650404577261471d-2
      
      real(kind=real2) :: beta
      
      beta = (real(k,kind=real2) - 0.25d0)*pi

      j0 = a0 + a1*beta**2 + a2*beta**4 + a3*beta**6
      j0 = j0/(beta*(b0 + b1*beta**2 + b2*beta**4 + b3*beta**6))
      j0 = beta + j0      
      
      return
   end function bessel0_root
   !**********************************************************************
end module bigquad
