SUBROUTINE COMPSP_GRID(pset,nti,specout)

  !This is a highly specialized routine to pre-compute and then
  !use sfh=4 tau-only models.  A grid of CSPs as a function
  !of age, tau, and metallicity are pre-computed and in subsequent
  !calls CSPs are returned by interpolating within the precomputed grid

  !Note that to get these large arrays to fit I switched to single precision

  USE sps_vars
  USE sps_utils, ONLY : intspec,locate
  IMPLICIT NONE

  TYPE(PARAMS), INTENT(in) :: pset
  INTEGER, INTENT(in) :: nti
  REAL(SP), DIMENSION(nspec), INTENT(inout) :: specout 
  INTEGER :: i,j,t,tlo,zlo
  TYPE(PARAMS) :: tpset
  REAL(SP), DIMENSION(nspec,ntfull) :: spec_ssp, tmp
  REAL(SP), DIMENSION(ntfull)   :: mass_ssp,lbol_ssp
  REAL(SP) :: mass_csp,lbol_csp,delt_burst,tau,const,mass_burst
  REAL(SP) :: lbol_burst,mdust,dt,sfstart,sftrunc,dz
  REAL(SP), DIMENSION(nspec) :: spec_burst,spec_csp
  REAL, DIMENSION(nspec,ntfull,ntaugrid,nz) :: tcsp_grid=0.0

  !-------------------------------------------------------------!
  !-------------------------------------------------------------!

  !set up the grid the first time this routine is called
  IF (csp_grid_flag.EQ.0) THEN 

     delt_burst = 0.0
     const      = 0.0
     mass_burst = 0.0
     spec_burst = 0.0
     lbol_burst = 0.0
     sfstart    = 0.0
     sftrunc    = 10**time_full(ntfull)

     DO j=1,nz

        tpset%zmet=j
        CALL SSP_GEN(tpset,mass_ssp,lbol_ssp,spec_ssp)
        IF (add_neb_emission.EQ.1) THEN
           tmp = spec_ssp
           CALL ADD_NEBULAR(pset,tmp,spec_ssp)
        ENDIF

        tpset     = pset
        tpset%sfh = 4

        DO t=1,ntaugrid

           !define the tau grid
           taugrid(t) = 10**(REAL(t)/ntaugrid*3-1)

           DO i=1,ntfull

              !compute the CSP
              CALL INTSPEC(tpset,i,spec_ssp,spec_csp,mass_ssp,lbol_ssp,&
                   mass_csp,lbol_csp,spec_burst,mass_burst,lbol_burst,&
                   delt_burst,sfstart,taugrid(t),const,sftrunc,mdust)
              tcsp_grid(:,i,t,j) = spec_csp

           ENDDO

        ENDDO

     ENDDO
     
     !now smooth in metallicity
     csp_grid = tcsp_grid
     DO j=2,nz-1
        csp_grid(:,:,:,j) = 0.2*tcsp_grid(:,:,:,j-1)+&
             0.6*tcsp_grid(:,:,:,j)+&
             0.2*tcsp_grid(:,:,:,j+1)
     ENDDO

     csp_grid_flag=1

  ELSE

     !the grid has already been defined, now interpolate in tau and logZ

     !set up interpolator in tau and Z
     tlo = MAX(MIN(locate(taugrid,pset%tau),ntaugrid-1),1)
     dt  = (pset%tau-taugrid(tlo))/(taugrid(tlo+1)-taugrid(tlo))
     zlo = MAX(MIN(locate(LOG10(zlegend/zsol),pset%logzsol),nz-1),1)
     dz  = (pset%logzsol-LOG10(zlegend(zlo)/zsol)) / &
          ( LOG10(zlegend(zlo+1)/zsol) - LOG10(zlegend(zlo)/zsol) )

     !interpolate in tau and Z
     specout = (1-dt)*(1-dz)*csp_grid(:,nti,tlo,zlo) + &
          (1-dt)*dz*csp_grid(:,nti,tlo,zlo+1) + &
          dt*(1-dz)*csp_grid(:,nti,tlo+1,zlo) + &
          dt*dz*csp_grid(:,nti,tlo+1,zlo+1)

  ENDIF


END SUBROUTINE COMPSP_GRID

