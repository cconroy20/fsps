FUNCTION COMPUTE_TAU1(mact,logt,logl,logg)

  !routine to compute tau at 1um from input
  !stellar parameters.  See Piovan et al. (2003)
  !for the origins of these equations
  
  USE sps_vars
  IMPLICIT NONE
  REAL(SP), INTENT(in)  :: mact,logt,logl,logg
  REAL(SP) :: compute_tau1

 
  compute_tau1 = 1.0

END FUNCTION COMPUTE_TAU1

!------------------------------------------------------------!
!------------------------------------------------------------!
!------------------------------------------------------------!

SUBROUTINE ADD_AGB_DUST(tspec,mact,logt,logl,logg,tco)

  USE sps_vars
  IMPLICIT NONE

  REAL(SP), DIMENSION(nspec), INTENT(out) :: tspec
  REAL(SP), INTENT(in)  :: mact,logt,logl,logg,tco
  INTEGER :: cstar,j,k
  REAL(SP) :: tau1,dj,dk, compute_tau1
  REAL(SP), DIMENSION(nspec) :: dusty
  
  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  !compute tau1 based on input stellar parameters
  tau1 = compute_tau1(mact,logt,logl,logg)

  IF (tco.GT.1) THEN 
     cstar=1 
  ELSE 
     cstar=0
  ENDIF

  !find dusty model given tau1,tco,Teff

  !need a dusty_tau1 grid and a dusty_teff grid

  !modify input stellar spectrum
  tspec = tspec * dusty


END SUBROUTINE ADD_AGB_DUST
