program convdiffmain

      use solvergmres
      use convdiff
      implicit none

      integer :: i
      integer :: np
      integer :: m,n,it
      real*8, allocatable :: x0(:),x(:),b(:)
      !external :: findAXupwind
      !external :: noprecond

      open(unit=5,file='soln.dat')

      np = 101
      h = 1.d0/(np-1)
      k = 0.1
      c = 1.0
      phi0 = 0.0
      phiL = 1.0

      n = np-2
      m = n/2
      it = n

      allocate(x0(n))
      allocate(x(n))
      allocate(b(n))

      x0 = 0.d0
      x  = 0.d0
      b  = 0.d0

      !call findBupwind(b,n)
      call findBcentral(b,n)

      !call performgmres(b,x0,x,m,n,it,findAXupwind,noprecond)
      call performgmres(b,x0,x,m,n,it,findAXcentral,noprecond)

      write(5,'(F10.5 F10.5)'),0.d0,phi0
      do i=1,n
      	write(5,'(F10.5 F10.5)'),i*h,x(i)
      enddo
      write(5,'(F10.5 F10.5)'),1.d0,phiL
      
      close(5)	

end program convdiffmain 
