FUNCTION GET_TUNIV(z)

  !compute age of Universe in Gyr at redshift z
  !assumes flat universe w/ only matter and lambda
  !assumes om0,ol0,H0 set in sps_vars.f90
  
  USE sps_vars
  IMPLICIT NONE
  INTEGER :: i
  INTEGER, PARAMETER :: ii=10000
  REAL(SP), INTENT(in) :: z
  REAL(SP) :: get_tuniv, thub
  REAL(SP), DIMENSION(ii) :: lnstig, hub

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  get_tuniv = 0.0

  !Hubble time in Gyr
  thub = 0.978E3 / H0

  DO i=1,ii
     lnstig(i) = REAL(i)/ii*(LOG(1E4)-LOG(1+z))+LOG(1+z)
  ENDDO
  
  hub = SQRT( om0*EXP(lnstig)**3 + ol0 )

  DO i=1,ii-1
     get_tuniv = get_tuniv + 0.5*(1/hub(i)+1/hub(i+1))
  ENDDO
  get_tuniv = get_tuniv * thub * (lnstig(2)-lnstig(1))


END FUNCTION GET_TUNIV
