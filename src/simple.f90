PROGRAM SIMPLE

    !set up modules
    USE sps_vars; USE sps_utils
  
    IMPLICIT NONE

    INTEGER :: i
    !define variable for SSP spectrum
    REAL(SP), DIMENSION(ntfull,nspec)  :: spec_ssp
    !define variables for Mass and Lbol info
    REAL(SP), DIMENSION(ntfull)    :: mass_ssp,lbol_ssp
    CHARACTER(100) :: file1='', file2=''
    !structure containing all necessary parameters
    TYPE(PARAMS) :: pset
    !define structure for CSP spectrum
    TYPE(COMPSPOUT), DIMENSION(ntfull) :: ocompsp
    REAL(SP) :: ssfr6,ssfr7,ssfr8,ave_age

    REAL, DIMENSION(16) :: dust
    CHARACTER, DIMENSION(10) :: files
    CHARACTER(15) :: tmp
    !CHARACTER(12) :: file1
    !---------------------------------------------------------------!
    !---------------------------------------------------------------!
  
    dust = (/0.000, 0.001, 0.003, 0.004, 0.01, &
             0.013, 0.023, 0.039, 0.067, 0.116, 0.20, 0.342, 0.6, & 
             1.014, 1.744, 3.0/)

    imf_type = 2    ! Salpeter IMF
    pset%zmet = 20    ! compute SSP

    !here we have to read in all the librarries
    CALL SPS_SETUP(pset%zmet)
    pset%sfh   = 0     !set SFH to "SSP"
    pset%zred  = 0.0  !redshift  
    pset%dust1 = 0.0   !dust parameter 1
    pset%dust2 = 0.0   !dust parameter 2

    pset%dell  = 0.0   !shift in log(L) for TP-nod stars
    pset%delt  = 0.0   !shift in log(Teff) for TP-nod stars
    pset%fbhb  = 0.0   !fraction of blue HB stars
    pset%sbss  = 0.0   !specific frequency of BS stars

    DO i=1,size(dust, 1)
        pset%dust2 = dust(i)
        pset%dust1 = dust(i)*3
        write(tmp, '(F10.5)'), dust(i)
        file1 = 'CompareDust/norm2/SSP_'//tmp(4:9)//'_nod.out'
        write(*,*) file1
        CALL SSP_GEN(pset, mass_ssp, lbol_ssp, spec_ssp)
        CALL COMPSP(3, 1, file1, mass_ssp, lbol_ssp, spec_ssp, pset, ocompsp)
    ENDDO 
    
    imf_type = 2    ! Salpeter IMF
    pset%zmet = 20   ! compute SSP
    pset%sfh = 1
    pset%tau = 1.0
    !compute all SSPs for different values of dust2
    !in the common block set up in sps_vars.f90
    DO i=1,size(dust, 1)
        pset%dust2 = dust(i)
        pset%dust1 = dust(i)*3
        write(tmp, '(F10.5)'), dust(i)
        file2 = 'CompareDust/norm2/CSP_'//tmp(4:9)//'_nod.out'
        write(*,*) file2
        CALL SSP_GEN(pset, mass_ssp, lbol_ssp, spec_ssp)
        CALL COMPSP(1, 1, file2, mass_ssp, lbol_ssp, spec_ssp, pset, ocompsp)
    ENDDO

END PROGRAM SIMPLE
