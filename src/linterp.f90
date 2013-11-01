FUNCTION LINTERP(xin,yin,xout)

  !routine to linearly interpolate a function yin(xin) at xout

  USE sps_vars; USE sps_utils, ONLY: locate
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(in) :: xin,yin
  REAL(SP), INTENT(in)  :: xout
  REAL(SP) :: linterp
  INTEGER :: klo,n

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  n   = SIZE(xin)
  klo = MAX(MIN(locate(xin,xout),n-1),1)

  linterp = yin(klo) + (yin(klo+1)-yin(klo))*&
       (xout-xin(klo))/(xin(klo+1)-xin(klo))

END FUNCTION LINTERP
