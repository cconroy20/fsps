PROGRAM GENERATE_TEST_DATA

  ! Generates reference data for FSPS regression testing.
  ! Uses allocatable arrays to support multiple compile-time configurations.

  USE sps_vars
  USE sps_utils
  IMPLICIT NONE

  ! Variables for SSP generation (allocatable)
  REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: spec_ssp
  REAL(SP), ALLOCATABLE, DIMENSION(:)   :: mass_ssp, lbol_ssp
  
  ! Variables for CSP generation (allocatable)
  TYPE(COMPSPOUT), ALLOCATABLE, DIMENSION(:) :: ocompsp
  
  ! Control variables
  TYPE(PARAMS) :: pset
  INTEGER :: i, unit_out
  CHARACTER(LEN=50) :: filename_out
  CHARACTER(LEN=100) :: csp_dummy_file

  ! Setup and initialization
  
  ! Name of the reference output file
  filename_out = 'sps_test_output.bin'
  unit_out = 40

  csp_dummy_file = 'dummy_csp.out'

  WRITE(*,*) '--------------------------------------------------'
  WRITE(*,*) 'FSPS REFERENCE DATA GENERATOR'
  
  ! Initialize FSPS parameters (MIST/MILES defaults)
  ! We call this FIRST so we can be sure parameters like ntfull/nspec are set
  ! before we allocate (though they are static in the current codebase).
  imf_type = 1          
  pset%zmet = 10        
  CALL SPS_SETUP(pset%zmet)

  ! Memory allocation
  WRITE(*,*) 'Allocating memory (nspec:', nspec, ' ntfull:', ntfull, ')...'
  ALLOCATE(spec_ssp(ntfull,nspec))
  ALLOCATE(mass_ssp(ntfull))
  ALLOCATE(lbol_ssp(ntfull))
  ALLOCATE(ocompsp(ntfull))

  ! Write header
  WRITE(*,*) 'Writing to file: ', TRIM(filename_out)
  OPEN(UNIT=unit_out, FILE=TRIM(filename_out), STATUS='REPLACE', &
       FORM='UNFORMATTED', ACCESS='STREAM')

  WRITE(*,*) 'Writing Header...'
  WRITE(unit_out) nspec
  WRITE(unit_out) ntfull
  WRITE(unit_out) nbands

  ! Test Case 1: Simple SSP (Solar, Chabrier)
  WRITE(*,*) 'Running Test Case 1: SSP (Solar, Chabrier)...'

  ! Define SSP Parameters
  pset%sfh   = 0     ! SSP
  pset%const = 0.0
  pset%zred  = 0.0
  pset%dust1 = 0.0
  pset%dust2 = 0.0
  add_neb_emission = 1 

  ! Compute SSP
  CALL SSP_GEN(pset, mass_ssp, lbol_ssp, spec_ssp)

  ! Write SSP Data
  WRITE(*,*) 'Saving SSP results...'
  WRITE(unit_out) mass_ssp
  WRITE(unit_out) lbol_ssp
  WRITE(unit_out) spec_ssp

  ! Test Case 2: CSP (Tau Model, Dusty)
  WRITE(*,*) 'Running Test Case 2: CSP (Tau=2.0, Dust=1.0)...'

  pset%sfh   = 1     ! Tau model
  pset%tau   = 2.0   
  pset%dust1 = 1.0   
  pset%dust2 = 0.3
  
  ! Re-run SSP_GEN
  CALL SSP_GEN(pset, mass_ssp, lbol_ssp, spec_ssp)

  ! Compute CSP
  CALL COMPSP(3, 1, csp_dummy_file, mass_ssp, lbol_ssp, spec_ssp, pset, ocompsp)

  ! Write CSP Data
  WRITE(*,*) 'Saving CSP results...'
  DO i = 1, ntfull
     WRITE(unit_out) ocompsp(i)%age
     WRITE(unit_out) ocompsp(i)%mass_csp
     WRITE(unit_out) ocompsp(i)%lbol_csp
     WRITE(unit_out) ocompsp(i)%sfr
     WRITE(unit_out) ocompsp(i)%mdust
     WRITE(unit_out) ocompsp(i)%mformed
     WRITE(unit_out) ocompsp(i)%mags      
     WRITE(unit_out) ocompsp(i)%spec      
     WRITE(unit_out) ocompsp(i)%indx      
     WRITE(unit_out) ocompsp(i)%emlines   
  END DO

  ! Cleanup
  CLOSE(unit_out)
  DEALLOCATE(spec_ssp, mass_ssp, lbol_ssp, ocompsp)

  WRITE(*,*) 'Complete. Data saved.'

END PROGRAM GENERATE_TEST_DATA