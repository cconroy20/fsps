FUNCTION ADD_AGN_DUST(lam,spec,pset)

  USE sps_vars
  IMPLICIT NONE

  REAL(SP), DIMENSION(nspec), INTENT(in) :: lam,spec
  TYPE(PARAMS), INTENT(in) :: pset
  REAL(SP), DIMENSION(nspec), INTENT(out) :: add_agn_dust
 
  !--------------------------------------------------------------!




END FUNCTION ADD_AGN_DUST
  
