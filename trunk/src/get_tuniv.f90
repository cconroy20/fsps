FUNCTION GET_TUNIV(z)

  !compute age of Universe at redshift z
  !assumes flat universe w/ only matter and lambda (from WMAP5).
  !assumes t_H set in sps_vars.f90
  
  USE nrtype
  IMPLICIT NONE
  INTEGER :: i
  REAL(SP), INTENT(in) :: z
  REAL(SP) :: get_tuniv
  REAL(SP), DIMENSION(10000) :: lnstig, hub, weights

  !set up integration weights
  weights                = 1.0
  weights(1:3)           = (/3./8,7./6,23./24/)
  weights(10000-2:10000) = (/23./24,7./6,3./8/)
  
  DO i=1,10000
     lnstig(i) = REAL(i)/1E4*(LOG(1E8)-LOG(1+z))+LOG(1+z)
  ENDDO
  
  hub       = SQRT( 0.258*EXP(lnstig)**3 + 0.742 )
  get_tuniv = (lnstig(2)-lnstig(1)) * SUM(weights/hub) * 13.7
  
END FUNCTION GET_TUNIV
