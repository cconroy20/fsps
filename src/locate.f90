FUNCTION LOCATE(xx,x)

  USE sps_vars
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: xx
  REAL(SP), INTENT(IN) :: x
  INTEGER :: locate
  INTEGER :: n,jl,jm,ju
  LOGICAL :: ascnd
  
  n = SIZE(xx)
  ascnd = (xx(n) >= xx(1))
  jl=0
  ju=n+1
  
  DO
     IF (ju-jl <= 1) EXIT
     jm=(ju+jl)/2
     IF (ascnd.EQV.(x >= xx(jm))) THEN
        jl=jm
     ELSE
        ju=jm
     ENDIF
  ENDDO

  IF (x == xx(1)) THEN
     locate=1
  ELSE IF (x == xx(n)) THEN
     locate=n-1
  ELSE
     locate=jl
  ENDIF

END FUNCTION LOCATE
