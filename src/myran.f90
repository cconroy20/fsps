FUNCTION MYRAN()

  !turns ran1 into a function, rather than a subroutine

  USE sps_vars
  USE nr, ONLY : ran1
  IMPLICIT NONE

  REAL(SP) :: myran

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  CALL ran1(myran)

END FUNCTION MYRAN
