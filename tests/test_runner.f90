PROGRAM TEST_RUNNER

  ! FSPS Regression Test Runner
  ! Reads reference data from disk, runs current FSPS code with identical
  ! parameters, and compares outputs with a relative tolerance.

  USE sps_vars
  USE sps_utils
  IMPLICIT NONE

  ! Exit codes
  INTEGER, PARAMETER :: EXIT_SUCCESS = 0
  INTEGER, PARAMETER :: EXIT_FAILURE = 1
  
  ! Default tolerance (can be overridden by env var FSPS_TEST_RTOL)
  REAL(SP), PARAMETER :: DEFAULT_RTOL = 1.0E-5

  ! Test arrays (allocatable)
  ! Reference data (read from disk)
  REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: ref_spec_ssp
  REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ref_mass_ssp, ref_lbol_ssp
  TYPE(COMPSPOUT), ALLOCATABLE, DIMENSION(:) :: ref_ocompsp

  ! New data (computed on the fly)
  REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: new_spec_ssp
  REAL(SP), ALLOCATABLE, DIMENSION(:)   :: new_mass_ssp, new_lbol_ssp
  TYPE(COMPSPOUT), ALLOCATABLE, DIMENSION(:) :: new_ocompsp

  ! Control variables
  TYPE(PARAMS) :: pset
  INTEGER :: i, unit_in, status
  CHARACTER(LEN=255) :: filename_in, env_buffer
  CHARACTER(LEN=100) :: csp_dummy_file
  
  ! Dimensions read from file
  INTEGER :: file_nspec, file_ntfull, file_nbands

  ! Comparison stats
  REAL(SP) :: rtol
  LOGICAL :: test_passed

  ! Configuration
  unit_in = 40
  csp_dummy_file = 'dummy_csp.out'
  test_passed = .TRUE.

  WRITE(*,*) '========================================='
  WRITE(*,*) 'FSPS TEST RUNNER'
  WRITE(*,*) '========================================='

  ! Read filename from command line
  CALL GET_COMMAND_ARGUMENT(1, filename_in, STATUS=status)
  IF (status /= 0 .OR. LEN_TRIM(filename_in) == 0) THEN
     WRITE(*,*) 'ERROR: Must provide reference filename as argument.'
     WRITE(*,*) 'Usage: ./test_runner tests/data/sps_ref_XXX.bin'
     STOP EXIT_FAILURE
  END IF

  ! Read tolerance from environment variable
  CALL GET_ENVIRONMENT_VARIABLE("FSPS_TEST_RTOL", VALUE=env_buffer, STATUS=status)
  IF (status == 0) THEN
     READ(env_buffer, *) rtol
     WRITE(*,*) 'Using RTOL from environment: ', rtol
  ELSE
     rtol = DEFAULT_RTOL
     WRITE(*,*) 'Using default RTOL: ', rtol
  END IF

  ! Initialize FSPS and check dimensions
  ! Note: We must initialize FSPS before allocating, but we must read the 
  ! file header before we know if dimensions match.
  
  imf_type = 1          
  pset%zmet = 10        
  CALL SPS_SETUP(pset%zmet)

  ! Open the reference file
  OPEN(UNIT=unit_in, FILE=TRIM(filename_in), STATUS='OLD', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
  IF (status /= 0) THEN
     WRITE(*,*) 'ERROR: Could not open file: ', TRIM(filename_in)
     STOP EXIT_FAILURE
  END IF

  ! Read Header
  READ(unit_in) file_nspec
  READ(unit_in) file_ntfull
  READ(unit_in) file_nbands

  WRITE(*,*) 'Reference Dimensions: nspec=', file_nspec, ' nt=', file_ntfull
  WRITE(*,*) 'Compiled Dimensions:  nspec=', nspec,      ' nt=', ntfull

  ! Strict dimension check
  IF (file_nspec /= nspec .OR. file_ntfull /= ntfull) THEN
     WRITE(*,*) 'FATAL: Binary dimensions do not match compiled FSPS dimensions.'
     WRITE(*,*) 'Ensure you are running the test with the same flags used to generate the data.'
     STOP EXIT_FAILURE
  END IF

  ! Allocate and read reference data
  ALLOCATE(ref_spec_ssp(ntfull,nspec))
  ALLOCATE(ref_mass_ssp(ntfull))
  ALLOCATE(ref_lbol_ssp(ntfull))
  ALLOCATE(ref_ocompsp(ntfull))

  ALLOCATE(new_spec_ssp(ntfull,nspec))
  ALLOCATE(new_mass_ssp(ntfull))
  ALLOCATE(new_lbol_ssp(ntfull))
  ALLOCATE(new_ocompsp(ntfull))

  WRITE(*,*) 'Reading SSP reference data...'
  READ(unit_in) ref_mass_ssp
  READ(unit_in) ref_lbol_ssp
  READ(unit_in) ref_spec_ssp

  WRITE(*,*) 'Reading CSP reference data...'
  DO i = 1, ntfull
     READ(unit_in) ref_ocompsp(i)%age
     READ(unit_in) ref_ocompsp(i)%mass_csp
     READ(unit_in) ref_ocompsp(i)%lbol_csp
     READ(unit_in) ref_ocompsp(i)%sfr
     READ(unit_in) ref_ocompsp(i)%mdust
     READ(unit_in) ref_ocompsp(i)%mformed
     READ(unit_in) ref_ocompsp(i)%mags      
     READ(unit_in) ref_ocompsp(i)%spec      
     READ(unit_in) ref_ocompsp(i)%indx      
     READ(unit_in) ref_ocompsp(i)%emlines   
  END DO
  CLOSE(unit_in)

  ! Generate new data
  WRITE(*,*) 'Generating new SSP data...'
  pset%sfh   = 0     
  pset%const = 0.0
  pset%zred  = 0.0
  pset%dust1 = 0.0
  pset%dust2 = 0.0
  add_neb_emission = 1 
  CALL SSP_GEN(pset, new_mass_ssp, new_lbol_ssp, new_spec_ssp)

  WRITE(*,*) 'Generating new CSP data...'
  pset%sfh   = 1    
  pset%tau   = 2.0   
  pset%dust1 = 1.0   
  pset%dust2 = 0.3
  CALL SSP_GEN(pset, new_mass_ssp, new_lbol_ssp, new_spec_ssp)
  CALL COMPSP(3, 1, csp_dummy_file, new_mass_ssp, new_lbol_ssp, new_spec_ssp, pset, new_ocompsp)

  ! Compare results
  WRITE(*,*) 'Verifying results (RTOL = ', rtol, ')...'
  CALL CHECK_ARRAY_2D("SSP Spectra", ref_spec_ssp, new_spec_ssp, ntfull, nspec)
  CALL CHECK_ARRAY_1D("SSP Mass", ref_mass_ssp, new_mass_ssp, ntfull)
  CALL CHECK_ARRAY_1D("SSP Lbol", ref_lbol_ssp, new_lbol_ssp, ntfull)

  ! Check CSP structure components
  DO i = 1, ntfull
     CALL CHECK_VAL("CSP Lbol", i, ref_ocompsp(i)%lbol_csp, new_ocompsp(i)%lbol_csp)
     CALL CHECK_VAL("CSP Mass", i, ref_ocompsp(i)%mass_csp, new_ocompsp(i)%mass_csp)
     
     ! Check spectra for the last time step as a proxy
     IF (i == ntfull) THEN
        CALL CHECK_ARRAY_1D("CSP Spec (Final)", ref_ocompsp(i)%spec, new_ocompsp(i)%spec, nspec)
     END IF
  END DO

  ! Report results
  WRITE(*,*) '--------------------------------------------------'
  IF (test_passed) THEN
     WRITE(*,*) 'TEST RESULT: PASS'
     STOP EXIT_SUCCESS
  ELSE
     WRITE(*,*) 'TEST RESULT: FAIL'
     STOP EXIT_FAILURE
  END IF


CONTAINS

  SUBROUTINE CHECK_ARRAY_2D(label, ref, new, d1, d2)
    CHARACTER(*), INTENT(IN) :: label
    INTEGER, INTENT(IN) :: d1, d2
    REAL(SP), DIMENSION(d1,d2), INTENT(IN) :: ref, new
    REAL(SP) :: delta, threshold
    INTEGER :: j, k

    DO k = 1, d2
       DO j = 1, d1
          delta = ABS(ref(j,k) - new(j,k))
          ! If ref is close to zero, use absolute tolerance, else relative
          threshold = MAX(ABS(ref(j,k)) * rtol, 1.0E-30) 
          
          IF (delta > threshold) THEN
             WRITE(*,*) 'FAIL: ', label, ' mismatch at index (',j,',',k,')'
             WRITE(*,*) '  Ref:', ref(j,k), ' New:', new(j,k), ' Diff:', delta
             test_passed = .FALSE.
             RETURN ! Return early on failure to avoid log spam
          END IF
       END DO
    END DO
  END SUBROUTINE CHECK_ARRAY_2D

  SUBROUTINE CHECK_ARRAY_1D(label, ref, new, d1)
    CHARACTER(*), INTENT(IN) :: label
    INTEGER, INTENT(IN) :: d1
    REAL(SP), DIMENSION(d1), INTENT(IN) :: ref, new
    REAL(SP) :: delta, threshold
    INTEGER :: j

    DO j = 1, d1
       delta = ABS(ref(j) - new(j))
       threshold = MAX(ABS(ref(j)) * rtol, 1.0E-30)
       
       IF (delta > threshold) THEN
          WRITE(*,*) 'FAIL: ', label, ' mismatch at index (',j,')'
          WRITE(*,*) '  Ref:', ref(j), ' New:', new(j), ' Diff:', delta
          test_passed = .FALSE.
          RETURN
       END IF
    END DO
  END SUBROUTINE CHECK_ARRAY_1D

  SUBROUTINE CHECK_VAL(label, idx, r, n)
    CHARACTER(*), INTENT(IN) :: label
    INTEGER, INTENT(IN) :: idx
    REAL(SP), INTENT(IN) :: r, n
    REAL(SP) :: delta, threshold

    delta = ABS(r - n)
    threshold = MAX(ABS(r) * rtol, 1.0E-30)

    IF (delta > threshold) THEN
       WRITE(*,*) 'FAIL: ', label, ' mismatch at step ', idx
       WRITE(*,*) '  Ref:', r, ' New:', n, ' Diff:', delta
       test_passed = .FALSE.
    END IF
  END SUBROUTINE CHECK_VAL

END PROGRAM TEST_RUNNER