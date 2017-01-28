SUBROUTINE ZTINTERP(zpos,spec,lbol,mass,tpos,zpow)

  !Linearly interpolate a grid of SSPs with the following options:
  !1) single metallicity (zpos) for a single age (tpos)
  !2) integrate over an MDF (zpos,zpow) for a grid of ages
  !3) single metallicity (zpos) for a grid of ages

  USE sps_vars; USE sps_utils, ONLY : locate, tsum
  IMPLICIT NONE

  REAL(SP),INTENT(in) :: zpos
  REAL(SP),INTENT(in), OPTIONAL :: tpos,zpow
  REAL(SP),INTENT(inout),DIMENSION(:) :: mass, lbol
  REAL(SP),INTENT(inout),DIMENSION(:,:) :: spec
  INTEGER  :: zlo,zhi,tlo,i
  REAL(SP) :: dz,dt,z0,imdf,w1=0.25,w2=0.5,w3=0.25
  REAL(SP), DIMENSION(nz) :: mdf

  !------------------------------------------------------------!


  !interpolate to a single metallicity and a single time
  IF (PRESENT(tpos)) THEN
     
     IF (SIZE(mass).GT.1.OR.SIZE(lbol).GT.1.OR.SIZE(spec(:,1)).GT.nspec) THEN
        WRITE(*,*) 'ZTINERP ERROR: you specified an age but are '//&
             'asking for the full age array as output!'
        STOP
     ENDIF
     
     IF (PRESENT(zpow)) THEN
        WRITE(*,*) 'ZTINERP ERROR: you cannot specify both an age and'//&
             'an MDF at the same time!'
        STOP
     ENDIF

     tlo = MAX(MIN(locate(time_full,tpos),ntfull-1),1)
     dt  = (tpos - time_full(tlo)) / (time_full(tlo+1) - time_full(tlo))
     zlo = MAX(MIN(locate(LOG10(zlegend/zsol),zpos),nz-1),1)
     dz  = (zpos-LOG10(zlegend(zlo)/zsol)) / &
          ( LOG10(zlegend(zlo+1)/zsol) - LOG10(zlegend(zlo)/zsol) )

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

     
     !integrate over an MDF and return a grid of ages
     !Note: this is a very crude implementation of an integral
     !over an MDF.  Really the MDF is just specifying a set of 
     !weights over which we sum the Z-dependent SSPs.
     IF (PRESENT(zpow)) THEN

        IF (zpow.LT.0) THEN
           ! This smooths the SSPs in metallicity using a simple
           ! triangular kernel given by w1,w2,w3.  The smoothed SSPs
           ! are then interpolated to the target metallicity.  
           mdf = 0.
           zlo = MAX(MIN(locate(LOG10(zlegend/zsol),zpos),nz-1),1)
           dz  = (zpos-LOG10(zlegend(zlo)/zsol)) / &
                ( LOG10(zlegend(zlo+1)/zsol) - LOG10(zlegend(zlo)/zsol) )
           ! Here the weights for the mdf are a combination of the
           ! triangular kernel weights and the the linear interpolation
           ! weights.  If the kernel extends off the metallicity grid
           ! then the out-of-bounds weights are accumulated at the edges.
           mdf(MAX(zlo-1,1)) = w1*(1-dz)
           mdf(MIN(zlo+2,nz)) = w3*dz
           mdf(zlo)   = mdf(zlo) + w2*(1-dz) + w1*dz
           mdf(zlo+1) = mdf(zlo+1) + w3*(1-dz) + w2*dz
           ! order matters here
           zhi = min(zlo+2, nz)
           zlo = max(zlo-1, 1)
        ELSE
           z0  = 10**zpos*zsol
           !mdf = ds/dlogZ
           mdf = ( zlegend/z0 * EXP(-zlegend/z0) )**zpow
           zlo = 1
           zhi = nz
        ENDIF
        mdf = mdf / SUM(mdf)

        mass = 0.
        lbol = 0.
        spec = 0.
        DO i=zlo,zhi
           mass = mass + mdf(i)*mass_ssp_zz(:,i)
           lbol = lbol + mdf(i)*lbol_ssp_zz(:,i)
           spec = spec + mdf(i)*spec_ssp_zz(:,:,i)
        ENDDO
 
     !interpolate to a single metallicity and return a grid of ages
     ELSE

        zlo = MAX(MIN(locate(LOG10(zlegend/zsol),zpos),nz-1),1)
        dz  = (zpos-LOG10(zlegend(zlo)/zsol)) / &
             ( LOG10(zlegend(zlo+1)/zsol) - LOG10(zlegend(zlo)/zsol) )
        
        mass = (1-dz)*mass_ssp_zz(:,zlo)   + dz*mass_ssp_zz(:,zlo+1)
        lbol = (1-dz)*lbol_ssp_zz(:,zlo)   + dz*lbol_ssp_zz(:,zlo+1)
        spec(:,:) = (1-dz)*spec_ssp_zz(:,:,zlo) + dz*spec_ssp_zz(:,:,zlo+1)
     
     ENDIF

  ENDIF


END SUBROUTINE ZTINTERP
