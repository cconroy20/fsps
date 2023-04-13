SUBROUTINE ADD_XRB(pset,sspi,sspo)

  ! Routine to add emission from X-ray binaries

  USE sps_vars
  USE sps_utils, ONLY : locate,tsum
  IMPLICIT NONE

  INTEGER :: t,a1,z1
  REAL(SP) :: da,dz,tmpz
  TYPE(PARAMS), INTENT(in) :: pset
  REAL(SP), INTENT(in), DIMENSION(nspec,ntfull)    :: sspi
  REAL(SP), INTENT(inout), DIMENSION(nspec,ntfull) :: sspo
  REAL(SP), DIMENSION(nspec) :: tmpspec

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  !set up the interpolation variables for logZ
  tmpz = log10(zlegend(pset%zmet)/zsol)
  z1   = MAX(MIN(locate(zmet_xrb,tmpz),nz_xrb-1),1)
  dz   = (tmpz-zmet_xrb(z1))/(zmet_xrb(z1+1)-zmet_xrb(z1))
  dz   = MAX(MIN(dz,1.0),0.0) !no extrapolation

  sspo = sspi
 
  DO t=1,nt

     !set up age interpolant
     a1 = MAX(MIN(locate(ages_xrb,time_full(t)),nt_xrb-1),1)
     da = (time_full(t)-ages_xrb(a1))/(ages_xrb(a1+1)-ages_xrb(a1))
     
     IF (da.LT.0.0.OR.da.GT.1.0) CYCLE

     tmpspec = &   !interpolate in logZ and time
          (1-da)*(1-dz)* spec_xrb(:,a1,z1)+&
          da*(1-dz)* spec_xrb(:,a1+1,z1)+&
          (1-da)*dz* spec_xrb(:,a1,z1+1)+&
          da*dz* spec_xrb(:,a1+1,z1+1)

     sspo(:,t) = sspo(:,t) + pset%frac_xrb * tmpspec

  ENDDO


END SUBROUTINE ADD_XRB
