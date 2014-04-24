SUBROUTINE ZTINTERP(zpos,spec,lbol,mass)

  !linearly interpolate the SSP(Z) to the correct metallicity

  USE sps_vars; USE sps_utils, ONLY : locate
  IMPLICIT NONE

  REAL(SP),INTENT(in) :: zpos
  REAL(SP),INTENT(inout),DIMENSION(ntfull) :: mass, lbol
  REAL(SP),INTENT(inout),DIMENSION(nspec,ntfull) :: spec
  INTEGER  :: zlo,klo
  REAL(SP) :: dz

  !------------------------------------------------------------!

  zlo = MAX(MIN(locate(LOG10(zlegend/0.0190),zpos),nz-1),1)
  dz  = (zpos-LOG10(zlegend(zlo)/0.0190)) / &
       ( LOG10(zlegend(zlo+1)/0.0190) - LOG10(zlegend(zlo)/0.0190) )

  mass = (1-dz)*mass_ssp_zz(:,zlo)   + dz*mass_ssp_zz(:,zlo+1)
  lbol = (1-dz)*lbol_ssp_zz(:,zlo)   + dz*lbol_ssp_zz(:,zlo+1)
  spec = (1-dz)*spec_ssp_zz(:,:,zlo) + dz*spec_ssp_zz(:,:,zlo+1)

  spec = 10**spec-tiny_number

END SUBROUTINE ZTINTERP
