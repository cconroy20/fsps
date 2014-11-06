FUNCTION TSUM(xin,yin)

  !simple trapezoidal integration of tabulated function (xin,yin)
  !NB: notice the abs(dx) in the equation below.  The x-axis can
  !therefore be in reverse order

  USE sps_vars
  IMPLICIT NONE

  REAL(SP), DIMENSION(:), INTENT(in) :: xin,yin
  REAL(SP) :: tsum
  INTEGER  :: nn

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  nn = SIZE(xin)

  tsum = SUM( ABS((xin(2:nn)-xin(1:nn-1))) * &
       (yin(2:nn)+yin(1:nn-1))/2. )


END FUNCTION TSUM
