FUNCTION LINTERPARR(xin,yin,xout)

  !routine to linearly interpolate a function yin(xin) at xout

  USE sps_vars; USE sps_utils, ONLY: locate
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(in) :: xin,yin
  REAL(SP), INTENT(in), DIMENSION(:) :: xout
  REAL(SP), DIMENSION(SIZE(xout)) :: linterparr
  INTEGER :: klo,n,n2,i

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  n   = SIZE(xin)
  n2  = SIZE(xout)

  DO i=1,n2
     klo = MAX(MIN(locate(xin,xout(i)),n-1),1)
     linterparr(i) = yin(klo) + (yin(klo+1)-yin(klo))*&
          (xout(i)-xin(klo))/(xin(klo+1)-xin(klo))
  ENDDO

END FUNCTION LINTERPARR
