!=========================================================================
subroutine noprecond(MinvX,X,n)

      implicit none

      integer,intent(in) :: n
      real*8, intent(in) :: X(n)
      real*8,intent(out) :: MinvX(n)

      MinvX=X

end subroutine noprecond
!=========================================================================
subroutine findAX(AX,X,n)
      
	implicit none

	integer,intent(in) :: n
        real*8, intent(in) :: X(n)
        real*8,intent(out) :: AX(n)

	real*8 :: sup,sub,diag
	integer :: i

	sup  =  1.0
	sub  =  1.0
	diag = -2.0

	AX = 0.d0

	do i=1,n

		if(i .eq. 1) then
			AX(i)=diag*x(i)+sup*x(i+1)
		else if(i .eq. n) then
			AX(i)=diag*x(i)+sub*x(i-1)
		else
			AX(i)=sub*x(i-1)+diag*x(i)+sup*x(i+1)
		endif
	enddo

end subroutine findAX
!=========================================================================
subroutine findbvec(b,n)
	
	integer, intent(in) :: n
	real*8, intent(out)  :: b(n)

	b=0.d0

	b(1) = 0.d0
	b(n) =-1.d0

end subroutine findbvec
!=========================================================================
program testgmres
      use solvergmres
      implicit none

      integer :: m,n,it,i
      real*8,allocatable :: x0(:),x(:),b(:)
      external :: findAX,noprecond

      open(unit=5,file="soln.dat")

      m=40
      n=50
      it=20

      allocate(x0(n))
      allocate(x(n))
      allocate(b(n))

      x0 = 0.d0
      x  = 0.d0

      call findbvec(b,n)

      call performgmres(b,x0,x,m,n,it,findAX,noprecond)

      do i=1,n
      	write(5,'(I5 F10.5)'),i,x(i)
      enddo

      close(5)

end program testgmres
!=========================================================================
