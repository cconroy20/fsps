SUBROUTINE init_random_seed()

  USE ran_state, ONLY : ran_seed
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed

  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  CALL SYSTEM_CLOCK(COUNT=clock)
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RAN_SEED(sequence=seed(1))
  
  DEALLOCATE(seed)

END SUBROUTINE init_random_seed
