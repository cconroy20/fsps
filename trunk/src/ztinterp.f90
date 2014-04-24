SUBROUTINE ZTINTERP(zpos,spec,lbol,mass,tpos)

  !linearly interpolate the SSP(Z) to the correct metallicity

  USE sps_vars; USE sps_utils, ONLY : locate
  IMPLICIT NONE

  REAL(SP),INTENT(in) :: zpos
  REAL(SP),INTENT(in), OPTIONAL :: tpos
  REAL(SP),INTENT(inout),DIMENSION(:) :: mass, lbol
  REAL(SP),INTENT(inout),DIMENSION(:,:) :: spec
  INTEGER  :: zlo,tlo
  REAL(SP) :: dz,dt

  !------------------------------------------------------------!

  zlo = MAX(MIN(locate(LOG10(zlegend/zsol),zpos),nz-1),1)
  dz  = (zpos-LOG10(zlegend(zlo)/zsol)) / &
       ( LOG10(zlegend(zlo+1)/zsol) - LOG10(zlegend(zlo)/zsol) )

  IF (PRESENT(tpos)) THEN
     
     IF (SIZE(mass).GT.1.OR.SIZE(lbol).GT.1.OR.SIZE(spec(:,1)).GT.nspec) THEN
        WRITE(*,*) 'ZTINERP ERROR: you specified an age but are '//&
             'asking for the full age array as output!'
        STOP
     ENDIF
     
     tlo  = MAX(MIN(locate(time_full,tpos),ntfull-1),1)
     dt   = (tpos - time_full(tlo)) / (time_full(tlo+1) - time_full(tlo))
     
     mass = (1-dz)*(1-dt)*mass_ssp_zz(tlo,zlo) + &
          dz*(1-dt)*mass_ssp_zz(tlo,zlo+1) + &
          (1-dz)*dt*mass_ssp_zz(tlo+1,zlo) + &
          dz*dt*mass_ssp_zz(tlo+1,zlo+1)
     
     lbol = (1-dz)*(1-dt)*lbol_ssp_zz(tlo,zlo) + &
          dz*(1-dt)*lbol_ssp_zz(tlo,zlo+1) + &
          (1-dz)*dt*lbol_ssp_zz(tlo+1,zlo) + &
          dz*dt*lbol_ssp_zz(tlo+1,zlo+1)
     
     spec(:,1) = (1-dz)*(1-dt)*spec_ssp_zz(:,tlo,zlo) + &
          dz*(1-dt)*spec_ssp_zz(:,tlo,zlo+1) + &
          (1-dz)*dt*spec_ssp_zz(:,tlo+1,zlo) + &
          dz*dt*spec_ssp_zz(:,tlo+1,zlo+1)
     
  ELSE
     
     mass = (1-dz)*mass_ssp_zz(:,zlo)   + dz*mass_ssp_zz(:,zlo+1)
     lbol = (1-dz)*lbol_ssp_zz(:,zlo)   + dz*lbol_ssp_zz(:,zlo+1)
     spec(:,:) = (1-dz)*spec_ssp_zz(:,:,zlo) + dz*spec_ssp_zz(:,:,zlo+1)
     
  ENDIF


END SUBROUTINE ZTINTERP
