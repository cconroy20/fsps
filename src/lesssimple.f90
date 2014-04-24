PROGRAM LESSSIMPLE

  !set up modules
  USE sps_vars; USE sps_utils
  
  IMPLICIT NONE

  !NB: the various structure types are defined in sps_vars.f90
  !    variables not explicitly defined here are defined in sps_vars.f90
  INTEGER :: i
  !define variable for SSP spectrum
  REAL(SP), DIMENSION(ntfull,nspec)  :: spec_pz
  !define variables for Mass and Lbol info
  REAL(SP), DIMENSION(ntfull)    :: mass_pz,lbol_pz
  CHARACTER(100) :: file2=''
  !structure containing all necessary parameters
  TYPE(PARAMS) :: pset
  !define structure for CSP spectrum
  TYPE(COMPSPOUT), DIMENSION(ntfull) :: ocompsp
  REAL(SP) :: zave

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
  
  ! Now we're going to show you how to use full  
  ! metallicity-dependent info                   
  
  imf_type = 0    ! Salpeter IMF
  pset%sfh = 0    ! compute SSP

  !here we have to read in all the librarries
  CALL SPS_SETUP(-1)

  !compute all SSPs (i.e. at all Zs)
  !nz and the various *ssp_zz arrays are stored 
  !in the common block set up in sps_vars.f90
  DO i=1,nz
     pset%zmet = i
     CALL SSP_GEN(pset,mass_ssp_zz(:,i),&
          lbol_ssp_zz(:,i),spec_ssp_zz(:,:,i))
  ENDDO

  !define the yield for a closed box distribution
  pset%pmetals = 0.02
  !compute SSP convolved with a closed box  
  CALL PZ_CONVOL(pset,zave,spec_pz,lbol_pz,mass_pz)
  file2    = 'SSP_pz.out'
  !now compute magnitudes for this SSP
  CALL COMPSP(1,1,file2,mass_pz,lbol_pz,spec_pz,pset,ocompsp)

  !run compsp for a tabulated sfh with a metallicity history
  !NB: one must have setup all the SSPs, as was done in the DO-loop above
  pset%sfh = 2
  file2    = 'CSP_tabsfh.out'
  CALL COMPSP(1,nz,file2,mass_ssp_zz,lbol_ssp_zz,&
          spec_ssp_zz,pset,ocompsp)


END PROGRAM LESSSIMPLE
