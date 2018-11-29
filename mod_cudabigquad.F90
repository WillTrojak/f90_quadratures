module cudabigquad
   use precision
   use cudafor
   implicit none

   real(kind=real2), device, allocatable :: xd(:),td(:),wd(:)
   real(kind=real2), device, save :: faca(21)
   real(kind=real2),  parameter :: tol = 3d-16
   
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
      use bigquad, only : stieltjes_c
      use cudafor
      implicit none

      integer(kind=int1), intent(in) :: n 

      real(kind=real2), intent(out) :: x(n),t(n)

      integer(kind=int1) :: n2,k,lim
      real(kind=real2) :: cn

      integer(kind=int1), device :: limd,nd,n2d
      real(kind=real2), device :: cnd

      type(dim3) :: blockt,grid
      
      allocate(xd(n))
      allocate(td(n))
      allocate(wd(n))
      
      n2 = floor(real(n,kind=real2)/2d0)
      blockt = dim3(256,1,1)
      grid = dim3(ceiling(real(n2)/blockt%x),1,1)

      faca = factorialArrayReal(20)
      
      cn = stieltjes_c(n)
      cnd = cn
      lim = max(10,floor(n/2000d0))
      limd = lim
      nd = n; n2d = n2
      
      call theta0<<<grid,blockt>>>(nd,n2d,td)      
      call gaussl_node_adjust<<<grid,blockt>>>(nd,n2d,limd,cnd,td,xd)

      t = td
      x = xd
     
      do k=1,n2
         t(n-k+1) = pi - t(k)
         x(n-k+1) = -x(k)
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
      use bigquad, only : stieltjes_c,baratella_dleg_approx
      use cudafor
      implicit none

      integer(kind=int1), intent(in) :: n 

      real(kind=real2), intent(out) :: w(n)

      integer(kind=int1) :: n2,k,lim
      real(kind=real2) :: cn,dpn

      integer(kind=int1), device :: limd,nd,n2d
      real(kind=real2), device :: cnd

      type(dim3) :: blockt,grid
      
      allocate(xd(n))
      allocate(td(n))
      allocate(wd(n))
      
      n2 = floor(real(n,kind=real2)/2d0)
      blockt = dim3(256,1,1)
      grid = dim3(ceiling(real(n2)/blockt%x),1,1)

      faca = factorialArrayReal(20)
      
      cn = stieltjes_c(n)
      cnd = cn
      lim = max(10,floor(n/2000d0))
      limd = lim
      nd = n; n2d = n2
      
      call theta0<<<grid,blockt>>>(nd,n2d,td)      
      call gaussl_node_adjust<<<grid,blockt>>>(nd,n2d,limd,cnd,td,xd,wd)

      w = wd
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
      use bigquad, only : stieltjes_c,baratella_dleg_approx
      use cudafor
      implicit none

      integer(kind=int1), intent(in) :: n 

      real(kind=real2), intent(out) :: x(n),w(n)

      integer(kind=int1) :: n2,k,lim,ierrc
      real(kind=real2) :: cn,dpn

      integer(kind=int1), device :: limd,nd,n2d
      real(kind=real2), device :: cnd

      type(dim3) :: blockt,grid
      
      allocate(xd(n))
      allocate(td(n))
      allocate(wd(n))
      
      n2 = floor(real(n,kind=real2)/2d0)
      blockt = dim3(256,1,1)
      grid = dim3(ceiling(real(n2)/blockt%x),1,1)

      faca = factorialArrayReal(20)
      
      cn = stieltjes_c(n)
      cnd = cn
      lim = max(10,floor(n/2000d0))
      limd = lim
      nd = n; n2d = n2
      
      call theta0<<<grid,blockt>>>(nd,n2d,td)
      ierrc = cudaDeviceSynchronize()
      call gaussl_node_adjust<<<grid,blockt>>>(nd,n2d,limd,cnd,td,xd,wd)

      x = xd
      w = wd
      
      do k=1,n2
         x(n-k+1) = -x(k)
         w(n-k+1) = w(k)
      enddo

      !print *,w

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
   attributes(global) subroutine gaussl_node_adjust(n,n2,lim,cn,t,x,w)
      use precision
      use cudafor
      use maths, only : pi      
      implicit none

      integer(kind=int1), device, intent(in) :: n,n2,lim
      
      real(kind=real2), device, intent(in) :: cn
      real(kind=real2), device, intent(inout) :: t(n)

      real(kind=real2), device, intent(out) :: x(n)
      real(kind=real2), device, optional, intent(out) :: w(n) 

      integer(kind=int1), device :: k

      k = blockDim%x*(blockIdx%x - 1) + threadIdx%x

      if(k .le. n2)then
         if(present(w))then
            call gaussl_newton(n,k,lim,cn,t(k),x(k),w(k))
         else
            call gaussl_newton(n,k,lim,cn,t(k),x(k))
         endif
      endif
      
      return
   end subroutine gaussl_node_adjust
   !**********************************************************************
   attributes(device) subroutine gaussl_newton(n,k,lim,cn,t,x,w)
      use precision
      use cudafor
      implicit none

      integer(kind=int1), intent(in) :: n,k,lim

      real(kind=real2), intent(in) :: cn
      
      real(kind=real2), intent(inout) :: t

      real(kind=real2), intent(out) :: x
      real(kind=real2), optional, intent(out) :: w

      integer(kind=int1) :: j
      integer(kind=int1), parameter :: jmax = 20

      real(kind=real2) :: delt,pn,dpn,ttol
      
      delt = 1d0

      j = 0
      ttol = t*tol
      if(k .gt. lim) then
         do while(abs(delt) .gt. ttol)
            j = j + 1
            call stieltjes_leg_approx(n,k,cn,t,pn,dpn)
            delt = pn/(dpn)
            
            t = t + delt
         enddo
         call stieltjes_leg_approx(n,k,cn,t,pn,dpn)
      else
         do while(abs(delt) .gt. ttol)
            j = j + 1
            pn  = bogaert_leg_approx(n,t)
            dpn = baratella_dleg_approx(n,t)
            delt = pn/(dpn)
            
            t = t + delt
         enddo
         dpn = baratella_dleg_approx(n,t)
      endif
      
      x = cos(t)
      ! For more accurate weights uncomment this line
      ! WARNING, this will be much slower
      !dpn = baratella_dleg_approx(n,t)
      if(present(w)) w = 2d0/(dpn*dpn)
      
      return
   end subroutine gaussl_newton
   !**********************************************************************
   !> @brief Approximate legnedre polynomial using Stieltjes
   !> @par This approximates a legendre polynomial using Stieltjes approximation,
   !! it is only accurate for large n and is note valid near the boundaries.
   !! Best performance is when  |x| < 0.866 (=cos(pi/6))
   !**********************************************************************
   attributes(device) subroutine stieltjes_leg_approx(n,k,cn,t,pn,dpn)
      use precision
      use cudafor
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
   attributes(device) function taylor_sin_kpi(k,del) result(z)
      use precision
      use cudafor
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
   attributes(device) function taylor_cos_kpi(k,del) result(z)
      use precision
      use cudafor
      implicit none

      integer(kind=int1), intent(in) :: k
      real(kind=real2), intent(in) :: del

      integer(kind=int1), parameter :: nm = 9
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
   attributes(device) function stieltjes_h(n,m) result(h)
      use precision
      use cudafor
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
   attributes(device) function bogaert_leg_approx(n,t) result(pn) ! best for cos(pi/6) <= |x| <= 1 as expensive
      use precision
      use cudafor
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
   attributes(device) function bogaert_f(y) result(fn)
      use precision
      use cudafor
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
   attributes(device) function baratella_dleg_approx(n,t) result(dpn) ! best for cos(pi/6) <= |x| <=1
      use precision
      use cudafor
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
   attributes(device) function bessel0d(l,z) result(jl)
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
         
         jl = jl + besseljn((-l + 2*j),z)*m1*c
         c = c*real(l-j,kind=real2)/real(j+1,kind=real2)
      enddo
      jl = jl*(2d0**(-l))

      return
   end function bessel0d
   !**********************************************************************
   attributes(global) subroutine theta0(n,n2,t)
      use precision
      use cudafor
      implicit none

      integer(kind=int1), intent(in) :: n,n2

      real(kind=real2), device, intent(out) :: t(:)

      integer(kind=int1) :: k
      real(kind=real2) :: xk

      k = blockDim%x*(blockIdx%x - 1) + threadIdx%x

      if(k .le. n2) then
         xk = tricomi_approx(n,k)
         if(abs(xk) .gt. 0.5d0) xk = olver_approx(n,k)
         
         t(k) = acos(xk)      
      endif
      
      return
   end subroutine theta0
   !**********************************************************************
   attributes(device) function tricomi_approx(n,k) result(xk) ! best for |xk| < 0.5
      use precision
      use cudafor
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
   attributes(device) function olver_approx(n,k) result(xk) ! best for 0.5 <= |xk| <=1
      use precision
      use cudafor
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
   attributes(device) function bessel0_root(k) result(j0)
      use precision
      use cudafor
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
   attributes(device) function besseljn(n,x) result(jn)
      ! this is a temporary fix to a bug in CUDA that is under review
      use precision
      implicit none

      integer(kind=int1), intent(in) :: n
      real(kind=real2), intent(in) :: x

      real(kind=real2) :: jn

      integer(kind=int1) :: m
      
      if(n .ge. 0)then
         jn = bessel_jn(n,x)
      else
         m = -n
         jn = bessel_jn(m,x)
         jn = (1d0 - 2d0*mod(m,2))*jn
      endif
      
      return
   end function besseljn
   !**********************************************************************
end module  cudabigquad
