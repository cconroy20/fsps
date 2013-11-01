FUNCTION FUNCINT(func,a,b)

  USE sps_vars; USE sps_utils, ONLY : tsum
  IMPLICIT NONE
  INTEGER, PARAMETER :: jmax=20, init=100,imax=10000000
  INTEGER :: nn,i,j
  REAL(SP), INTENT(IN) :: a,b
  REAL(SP) :: funcint, last, itol=1E-3
  REAL(SP), DIMENSION(imax) :: x,y

  INTERFACE
     FUNCTION func(x)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(SIZE(x)) :: func
     END FUNCTION func
  END INTERFACE

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  last = -99999.

  IF (a.EQ.b) THEN
     funcint=0.0
     RETURN
  ENDIF

  DO j=1,jmax
  
     nn = init*INT(2**j)
     IF (nn.GT.imax) THEN 
        WRITE(*,*) 'FUNCTINT ERROR: nn>imax'
        STOP
     ENDIF

     DO i=1,nn
        x(i) = i/REAL(nn)*(b-a)+a
     ENDDO
     funcint = tsum(x(1:nn),func(x(1:nn)))

     !write(*,*) nn,funcint,abs(last-funcint)/funcint

     IF ( abs(last-funcint)/funcint.LT.itol) RETURN
     last = funcint
     
  ENDDO

  WRITE(*,*) 'FUNCINT ERROR: too many steps'
  STOP


END FUNCTION FUNCINT
