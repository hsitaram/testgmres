module convdiff

implicit none

integer,parameter :: minsize=5
real*8 :: k
real*8 :: c
real*8 :: h
real*8 :: phi0
real*8 :: phiL

contains

!===============================================================
subroutine findAXupwind(AX,X,n)

      integer, intent(in) :: n
      real*8, intent(out) :: AX(n)
      real*8, intent(in)  :: X(n)

      real*8 :: sup,sub,diag
      integer :: i

      sub  = (k + c*h)/h/h
      diag = (-2.d0*k - c*h)/h/h
      sup  = k/h/h

      AX(1)=X(1)
      do i=2,n-1
	AX(i)=sub*X(i-1)+diag*X(i)+sup*X(i+1)
      enddo
      AX(n)=X(n)

end subroutine findAXupwind
!===============================================================
subroutine findAXcentral(AX,X,n)

	integer, intent(in) :: n
	real*8, intent(out) :: AX(n)
	real*8, intent(in)  :: X(n)

	real*8 :: sup,sub,diag
	integer :: i

	sub  = (k + 0.5*c*h)/h*h
	diag = (-2.d0*k)/h/h
	sup  = (k - 0.5*c*h)/h/h

	AX(1) = X(1)
	do i=2,n-1
  	   AX(i)=sub*X(i-1)+diag*X(i)+sup*X(i+1)
        enddo
	AX(n) = X(n)

end subroutine findAXcentral
!===============================================================
subroutine findBupwind(b,n)

      integer, intent(in) :: n
      real*8, intent(out) :: b(n)
      integer :: i

      b(1) = phi0
      b(n) = phiL
      !do i=1,n
      !	b(i)=phi0+(i-1)*(phiL-phi0)/(n-1)
      !enddo	

end subroutine findBupwind
!===============================================================
subroutine findBcentral(b,n)

	integer,intent(in) :: n
	real*8, intent(out) :: b(n)

	b(1) = phi0
	b(n) = phiL

end subroutine findBcentral
!===============================================================
subroutine pjacobiprecond(MinvX,X,n)

      integer, intent(in) :: n
      real*8, intent(out) :: MinvX(n)
      real*8, intent(in)  :: X(n)

      real*8 :: diag

      diag=-(2.d0*k + c*h)/h/h
      !diag=-2.d0*k/h/h
      !diag = 1.d0

      MinvX = X
      MinvX(2:n-1) = X(2:n-1)/diag

end subroutine pjacobiprecond
!===============================================================
subroutine noprecond(MinvX,X,n)

      integer, intent(in) :: n
      real*8, intent(out) :: MinvX(n)
      real*8, intent(in)  :: X(n)

      MinvX = X

end subroutine noprecond
!===============================================================
subroutine sorprecond(MinvX,X,n)

	integer, intent(in) :: n
	real*8, intent(out) :: MinvX(n)
	real*8, intent(in)  :: X(n)

	real*8 :: sub,diag,w
	integer :: i

	w=1.d0

	sub = (k+c*h)/h/h
        diag=-(2.d0*k + c*h)/h/h

	MinvX(1) = X(1)
	MinvX(n) = X(n)

	do i=2,n-1
		MinvX(i)=w*(X(i)-sub*MinvX(i-1))/diag
	enddo

end subroutine sorprecond
!===============================================================
subroutine gauss_seidel_smoothing(res,b,X,h,n)

	integer, intent(in)   :: n
	real*8, intent(in)    :: h
	real*8, intent(in)    :: b(n)
	real*8, intent(inout) :: X(n)
	real*8, intent(out)   :: res(n)

	integer :: i

	real*8 :: sup,sub,diag

	diag = (-2.0*k-c*h)/h**2
	sub  = (k+c*h)/h**2
	sup  = k/h**2

	X(1) = b(1)
	X(n) = b(n)
	do i=2,n-1
	   X(i) = ( b(i) - sub*X(i-1) - sup*X(i+1) )/diag
	enddo

	res = 0.d0

	do i=2,n-1
	   res(i) = b(i) - diag*X(i) - sub*X(i-1) - sup*X(i+1)
	enddo
		

end subroutine gauss_seidel_smoothing
!===============================================================
subroutine restriction(Xh,X2h,n2h) 

	integer, intent(in) :: n2h
	real*8, intent(inout) :: Xh(2*n2h-1)
	real*8, intent(inout) :: X2h(n2h)

	integer :: i

	X2h(1)     = Xh(1)
	X2h(n2h) = Xh(2*n2h-1)

	do i=2,n2h-1
	   X2h(i)=0.25*(Xh(2*i-2) + 2.0*Xh(2*i-1) + Xh(2*i))
	enddo


end subroutine restriction
!===============================================================
subroutine prolong(Xh,X2h,n2h)

	integer, intent(in) :: n2h
	real*8, intent(inout) :: Xh(2*n2h-1)
	real*8, intent(inout) :: X2h(n2h)

	integer :: i

	do i=1,n2h-1
	   Xh(2*i-1) = X2h(i)
	   Xh(2*i)   = 0.5*(X2h(i) + X2h(i+1))
	enddo

	Xh(2*n2h-1) = X2h(n2h)
	
end subroutine prolong
!===============================================================
recursive subroutine dovcycle(X,b,dx,n)

	integer, intent(in) :: n
	real*8, intent(inout) :: X(n)
	real*8, intent(inout) :: b(n)
	real*8, intent(in) :: dx

	real*8  :: resh(n)
	real*8 :: res2h(n/2+1)
	real*8 :: e2h(n/2+1)
	real*8 :: eh(n)

	if(n .le. minsize) then
		call exactgaussjordansolve(b,X,dx,n)
	else

		!initial smoothing
		call gauss_seidel_smoothing(resh,b,X,dx,n)

		!restriction of residual from fine to coarse grid
        	call restriction(resh,res2h,n/2+1)

		e2h = 0.d0
		call dovcycle(e2h,res2h,2*dx,n/2+1)

		!prolong error from coarse to fine grid
		call prolong(eh,e2h,n/2+1)

		!update
		X(:) = X(:) + eh(:)	
		!post smooth
		call gauss_seidel_smoothing(resh,b,X,dx,n)
	endif

end subroutine dovcycle
!===============================================================
subroutine mgridprecond(MinvX,X,n)
	
	integer, intent(in)    :: n
	real*8, intent(out)    :: MinvX(n)
	real*8, intent(inout)  :: X(n)
	
	integer :: nvcycles,i

	nvcycles=4

	do i=1,nvcycles
		call dovcycle(MinvX,X,h,n)
	enddo

end subroutine mgridprecond
!===============================================================
subroutine exactgaussjordansolve(b,X,h,n)

        implicit none

	integer, intent(in)    :: n
	real*8, intent(in)     :: b(n)
	real*8, intent(inout)  :: X(n)
	real*8, intent(in)     :: h

	real*8 :: diag(n),sub(2:n),sup(1:n-1)
	integer :: i
	real*8  :: factor
	real*8  :: new_b(n)

	diag = (-2.0*k - c*h)/h**2
	sub  = (k + c*h)/h**2
	sup  = k/h**2
	!diag = -2.0*k
	!sub  = (k+0.5*c*h)
	!sup  = (k-0.5*c*h)


	diag(1) = 1.d0
	diag(n) = 1.d0

	sup(1)  = 0.d0
	sub(n)  = 0.d0

	new_b   = b

	!make sub diagonal 0
	do i=2,n

		if(sub(i) .ne. 0) then	

			factor   = diag(i-1)/sub(i)
			sub(i)   = factor*sub(i)   - diag(i-1)
			diag(i)  = factor*diag(i)  - sup(i-1)
			sup(i)   = factor*sup(i)
			new_b(i) = factor*new_b(i) - new_b(i-1)

		endif

	enddo

	do i=n,1,-1
		
		X(i) = new_b(i)
		
		if(i .le. (n-1)) then
			X(i)=X(i)-sup(i)*X(i+1)
		endif

		X(i) = X(i)/diag(i)
	enddo

end subroutine exactgaussjordansolve
!===============================================================

end module convdiff
