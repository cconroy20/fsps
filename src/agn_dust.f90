FUNCTION AGN_DUST(lam,spec,pset,lbol_csp)

  USE sps_vars; USE sps_utils, ONLY: locate,attn_curve
  IMPLICIT NONE

  REAL(SP), DIMENSION(nspec), INTENT(in) :: lam,spec
  REAL(SP), INTENT(in)       :: lbol_csp
  TYPE(PARAMS), INTENT(in)   :: pset
  REAL(SP), DIMENSION(nspec) :: agn_dust,agnspeci
  INTEGER  :: jlo
  REAL(SP) :: dj
 
  !--------------------------------------------------------------!

  !interpolate in tau_agn
  jlo = MIN(MAX(locate(agndust_tau,pset%agn_tau),1),&
       nagndust-1)
  dj  = (pset%agn_tau-agndust_tau(jlo)) / &
       (agndust_tau(jlo+1)-agndust_tau(jlo))
  dj  = MAX(MIN(dj,1.0),0.0) !no extrapolation

  agnspeci  = (1-dj)*agndust_spec(:,jlo) + dj*agndust_spec(:,jlo+1)

  !attenuate the AGN emission by the diffuse dust
  agnspeci = agnspeci*EXP(-attn_curve(spec_lambda,dust_type,pset))

  agn_dust = spec + 10**lbol_csp*pset%fagn*agnspeci


END FUNCTION AGN_DUST
  
