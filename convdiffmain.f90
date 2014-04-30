program convdiffmain

      use solvergmres
      use convdiff
      implicit none

      integer :: i
      integer :: np
      integer :: m,n,it
      real*8, allocatable :: x0(:),x(:),b(:)
      integer :: precondoption

      print *,"**********************************************"
      print *,"Solving convection diffusion boundary layer equations &
		      using upwind discretization"
      print *,"**********************************************"
      print *,"preconditioning options:"
      print *,"1 no preconditioning"
      print *,"2 mgrid preconditioning"
      print *,"3 pjacobi preconditioning"
      print *,"4 sor preconditioning"
      
      read(*,*)precondoption
      
      open(unit=5,file='soln.dat')

      np = 4097
      h = 1.d0/(np-1)
      k = 0.1
      c = 1.d0
      phi0 = 0.d0
      phiL = 1.d0

      n = np
      m = 10
      it = 40

      allocate(x0(n))
      allocate(x(n))
      allocate(b(n))

      x0 = 0.d0
      x  = 0.d0
      b  = 0.d0

      call findBupwind(b,n)

      if(precondoption .eq. 1) then
      	call performgmres(b,x0,x,m,n,it,findAXupwind,noprecond)
      else if(precondoption .eq. 2) then
      	call performgmres(b,x0,x,m,n,it,findAXupwind,mgridprecond)
      else if(precondoption .eq. 3) then
      	call performgmres(b,x0,x,m,n,it,findAXupwind,pjacobiprecond)
      else if(precondoption .eq. 4) then
      	call performgmres(b,x0,x,m,n,it,findAXupwind,sorprecond)
      endif	
      !call performgmres(b,x0,x,m,n,it,findAXcentral,noprecond)

      do i=1,n
      	write(5,'(F10.5 F10.5)'),(i-1)*h,x(i)
      enddo
      
      close(5)	

end program convdiffmain 
