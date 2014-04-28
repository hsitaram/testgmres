module convdiff

implicit none

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

      sub  = k + c*h
      diag = -2.d0*k - c*h
      sup  = k

      do i=1,n
	
	if(i .eq. 1) then
	  	AX(i)=diag*X(i)+sup*X(i+1)
	else if(i .eq. n) then
		AX(i)=sub*X(i-1)+diag*X(i)
	else
		AX(i)=sub*X(i-1)+diag*X(i)+sup*X(i+1)
	endif
      enddo

end subroutine findAXupwind
!===============================================================
subroutine findAXcentral(AX,X,n)

	integer, intent(in) :: n
	real*8, intent(out) :: AX(n)
	real*8, intent(in)  :: X(n)

	real*8 :: sup,sub,diag
	integer :: i

	sub  = k + 0.5*c*h
	diag = -2.d0*k
	sup  = k - 0.5*c*h

	do i=1,n
	 
          if(i .eq. 1) then
	  	AX(i)=diag*X(i)+sup*X(i+1)
      	  else if(i .eq. n) then
  		AX(i)=sub*X(i-1)+diag*X(i)
	  else
  		AX(i)=sub*X(i-1)+diag*X(i)+sup*X(i+1)
  	  endif
       enddo

end subroutine findAXcentral
!===============================================================
subroutine findBupwind(b,n)

      integer, intent(in) :: n
      real*8, intent(out) :: b(n)

      b(1) = -(k + c*h)*phi0
      b(n) =  -k*phiL

end subroutine findBupwind
!===============================================================
subroutine findBcentral(b,n)

	integer,intent(in) :: n
	real*8, intent(out) :: b(n)

	b(1) = -(k+0.5*c*h)*phi0
	b(n) = -(k-0.5*c*h)*phiL

end subroutine findBcentral
!===============================================================
subroutine noprecond(MinvX,X,n)

      integer, intent(in) :: n
      real*8, intent(out) :: MinvX(n)
      real*8, intent(in)  :: X(n)

      MinvX = X

end subroutine noprecond
!===============================================================

end module convdiff
