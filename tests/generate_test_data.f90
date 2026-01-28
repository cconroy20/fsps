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
  INTEGER :: i, unit_out, status, arg_count
  CHARACTER(LEN=50) :: filename_out
  CHARACTER(LEN=100) :: csp_dummy_file
  CHARACTER(LEN=255) :: arg_val
  CHARACTER(LEN=20) :: isoc_arg, spec_arg, dust_arg
  LOGICAL :: isoc_set, spec_set, dust_set

  ! Exit codes
  INTEGER, PARAMETER :: EXIT_FAILURE = 1

  ! Setup and initialization
  
  ! Name of the reference output file
  filename_out = 'sps_test_output.bin'
  unit_out = 40

  csp_dummy_file = 'dummy_csp.out'
  isoc_set = .FALSE.
  spec_set = .FALSE.
  dust_set = .FALSE.

  WRITE(*,*) '--------------------------------------------------'
  WRITE(*,*) 'FSPS REFERENCE DATA GENERATOR'

  ! Parse command line arguments
  arg_count = COMMAND_ARGUMENT_COUNT()
  i = 1
  DO WHILE (i <= arg_count)
     CALL GET_COMMAND_ARGUMENT(i, arg_val, STATUS=status)
     IF (status /= 0) EXIT
     
     IF (TRIM(arg_val) == '--isoc') THEN
        i = i + 1
        IF (i <= arg_count) THEN
           CALL GET_COMMAND_ARGUMENT(i, isoc_arg, STATUS=status)
           isoc_set = .TRUE.
        ELSE
           WRITE(*,*) 'ERROR: --isoc requires an argument'
           STOP EXIT_FAILURE
        END IF
     ELSE IF (TRIM(arg_val) == '--spec') THEN
        i = i + 1
        IF (i <= arg_count) THEN
           CALL GET_COMMAND_ARGUMENT(i, spec_arg, STATUS=status)
           spec_set = .TRUE.
        ELSE
           WRITE(*,*) 'ERROR: --spec requires an argument'
           STOP EXIT_FAILURE
        END IF
     ELSE IF (TRIM(arg_val) == '--dust') THEN
        i = i + 1
        IF (i <= arg_count) THEN
           CALL GET_COMMAND_ARGUMENT(i, dust_arg, STATUS=status)
           dust_set = .TRUE.
        ELSE
           WRITE(*,*) 'ERROR: --dust requires an argument'
           STOP EXIT_FAILURE
        END IF
     END IF
     i = i + 1
  END DO
  
  ! Initialize FSPS parameters (MIST/MILES defaults)
  ! We call this FIRST so we can be sure parameters like ntfull/nspec are set
  ! before we allocate (though they are static in the current codebase).
  imf_type = 1          
  pset%zmet = 10        
  
  ! Use optional arguments if set
  IF (isoc_set .AND. spec_set .AND. dust_set) THEN
      CALL SPS_SETUP(pset%zmet, isoc_type_in=TRIM(isoc_arg), spec_type_in=TRIM(spec_arg), dust_type_in=TRIM(dust_arg))
  ELSE IF (isoc_set .AND. spec_set) THEN
      CALL SPS_SETUP(pset%zmet, isoc_type_in=TRIM(isoc_arg), spec_type_in=TRIM(spec_arg))
  ELSE IF (isoc_set .AND. dust_set) THEN
      CALL SPS_SETUP(pset%zmet, isoc_type_in=TRIM(isoc_arg), dust_type_in=TRIM(dust_arg))
  ELSE IF (spec_set .AND. dust_set) THEN
      CALL SPS_SETUP(pset%zmet, spec_type_in=TRIM(spec_arg), dust_type_in=TRIM(dust_arg))
  ELSE IF (isoc_set) THEN
      CALL SPS_SETUP(pset%zmet, isoc_type_in=TRIM(isoc_arg))
  ELSE IF (spec_set) THEN
      CALL SPS_SETUP(pset%zmet, spec_type_in=TRIM(spec_arg))
  ELSE IF (dust_set) THEN
      CALL SPS_SETUP(pset%zmet, dust_type_in=TRIM(dust_arg))
  ELSE
      CALL SPS_SETUP(pset%zmet)
  END IF

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
  
  ! Allocate pset allocatable components
  IF (ALLOCATED(pset%mag_compute)) DEALLOCATE(pset%mag_compute)
  ALLOCATE(pset%mag_compute(nbands))
  pset%mag_compute = 1
  
  IF (ALLOCATED(pset%ssp_gen_age)) DEALLOCATE(pset%ssp_gen_age)
  ALLOCATE(pset%ssp_gen_age(nt))
  pset%ssp_gen_age = 1

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

  ! Manually allocate components of ocompsp array elements
  DO i = 1, ntfull
     IF (.NOT. ALLOCATED(ocompsp(i)%mags)) ALLOCATE(ocompsp(i)%mags(nbands))
     IF (.NOT. ALLOCATED(ocompsp(i)%spec)) ALLOCATE(ocompsp(i)%spec(nspec))
     IF (.NOT. ALLOCATED(ocompsp(i)%indx)) ALLOCATE(ocompsp(i)%indx(nindx))
     IF (.NOT. ALLOCATED(ocompsp(i)%emlines)) ALLOCATE(ocompsp(i)%emlines(nemline))
  END DO

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
  IF (ALLOCATED(pset%mag_compute)) DEALLOCATE(pset%mag_compute)
  IF (ALLOCATED(pset%ssp_gen_age)) DEALLOCATE(pset%ssp_gen_age)
  
  CALL SPS_TAKEDOWN()

  WRITE(*,*) 'Complete. Data saved.'

END PROGRAM GENERATE_TEST_DATA
