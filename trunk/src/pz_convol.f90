SUBROUTINE PZ_CONVOL(yield,zave,spec_pz,lbol_pz,mass_pz)

  !routine to weight SSP(Z) by the MDF: P(Z)=Z**zpow*EXP(-Z/p).  
  !The yield is the only input. 
  !zpow is set in sps_vars.f90; zpow=1 yields a closed-box
  !requires that all the SSPs are set up in the common block
  !variables spec_ssp_zz, mass_ssp_zz, and lbol_ssp_zz
  !The average metallicity is returned as zave

  USE sps_vars; USE sps_utils, ONLY : linterp
  IMPLICIT NONE
  
  INTEGER  :: i,t,z
  REAL(SP) :: norm
  REAL(SP), INTENT(out), DIMENSION(nspec,ntfull) :: spec_pz
  REAL(SP), INTENT(out), DIMENSION(ntfull) :: mass_pz, lbol_pz
  REAL(SP), INTENT(out)    :: zave
  REAL(SP), DIMENSION(nz)  ::  pzz1
  REAL(SP), DIMENSION(100) :: pzz2,zz2,zzspec
  REAL(SP), INTENT(in) :: yield

  !-----------------------------------------------------------!

  spec_pz = 0.0  
  lbol_pz = 0.0
  mass_pz = 0.0
  norm    = 0.0
  zave    = 0.0

  !if using the Padova+BaSeL model, just use the native Z grid
  IF (nz.EQ.22) THEN

     !define P(Z)
     pzz1 = zlegend**zpow2 * EXP(-zlegend/yield)
     
     !integrate over P(Z)
     DO z=1,nz-1
        spec_pz = spec_pz + (LOG(zlegend(z+1))-LOG(zlegend(z))) * &
             (pzz1(z+1)*spec_ssp_zz(:,:,z+1) + pzz1(z)*spec_ssp_zz(:,:,z))/2.
        lbol_pz = lbol_pz + (LOG(zlegend(z+1))-LOG(zlegend(z))) * &
             (pzz1(z+1)*lbol_ssp_zz(:,z+1) + pzz1(z)*lbol_ssp_zz(:,z))/2.
        mass_pz = mass_pz + (LOG(zlegend(z+1))-LOG(zlegend(z))) * &
             (pzz1(z+1)*mass_ssp_zz(:,z+1) + pzz1(z)*mass_ssp_zz(:,z))/2.
        norm    = norm + (LOG(zlegend(z+1))-LOG(zlegend(z))) * &
          (pzz1(z+1)+pzz1(z))/2.
        zave    = zave + (LOG(zlegend(z+1))-LOG(zlegend(z))) * &
          (pzz1(z+1)*zlegend(z+1)+pzz1(z)*zlegend(z))/2.
     ENDDO
     
  ELSE

     !set up Z grid
     DO i=1,100
        zz2(i) = 10**(REAL(i)/100.*(-1+4.0)-4)
     ENDDO

     !define P(Z)
     pzz2 = zz2**zpow2 * EXP(-zz2/yield)
     !pzz2 = zz2/ (1/zz2 + 1E4*zz2**zpow)

     DO t=nt,nt
        DO i=1,nspec
           !interpolate
           DO z=1,100
              zzspec(z) = 10**linterp(log10(zlegend),&
                   log10(spec_ssp_zz(i,t,:)),log10(zz2(z)))
           ENDDO
           !integrate
           DO z=1,100-1
              spec_pz(i,t) = spec_pz(i,t) + (LOG(zz2(z+1))-LOG(zz2(z))) * &
                   (pzz2(z+1)*zzspec(z+1) + pzz2(z)*zzspec(z))/2.
           ENDDO
        ENDDO
     ENDDO

     DO z=1,100-1
        norm = norm + (LOG(zz2(z+1))-LOG(zz2(z))) * &
             (pzz2(z+1)+pzz2(z))/2.
        zave = zave + (LOG(zz2(z+1))-LOG(zz2(z))) * &
             (pzz2(z+1)*zz2(z+1)+pzz2(z)*zz2(z))/2.
     ENDDO

  ENDIF

 
  spec_pz = spec_pz / norm
  lbol_pz = lbol_pz / norm
  mass_pz = mass_pz / norm
  zave    = zave    / norm

END SUBROUTINE PZ_CONVOL
