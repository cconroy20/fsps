FUNCTION GET_LUMDIST(z)

  !compute luminosity distance to redshift z
  !assumes flat universe w/ only matter and lambda
  !assumes om0,ol0,H0 set in sps_vars.f90
  
  USE sps_vars; USE sps_utils, ONLY : tsum
  IMPLICIT NONE
  INTEGER :: i
  INTEGER, PARAMETER :: ii=10000
  REAL(SP), INTENT(in) :: z
  REAL(SP) :: get_lumdist, dhub
  REAL(SP), DIMENSION(ii) :: zz, hub

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  get_lumdist = 0.0

  !Hubble distance in pc
  dhub = clight/1E13/H0*1E6

  DO i=1,ii
     zz(i) = REAL(i)/ii*z
  ENDDO
  
  hub = SQRT( om0*(1+zz)**3 + ol0 )

  get_lumdist = TSUM(zz,1/hub) * (1+z) * dhub


END FUNCTION GET_LUMDIST
