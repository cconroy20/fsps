FUNCTION MYARTH(first,increment,n)

  USE sps_vars
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: first,increment
  INTEGER, INTENT(IN) :: n
  INTEGER, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
  REAL(SP), DIMENSION(n) :: myarth
  INTEGER :: k,k2
  REAL(SP) :: temp

  IF (n > 0) myarth(1)=first
  IF (n <= NPAR_ARTH) THEN
     DO k=2,n
        myarth(k)=myarth(k-1)+increment
     ENDDO
  ELSE
     DO k=2,NPAR2_ARTH
        myarth(k)=myarth(k-1)+increment
     ENDDO
     temp=increment*NPAR2_ARTH
     k=NPAR2_ARTH
     DO
        IF (k >= n) EXIT
        k2=k+k
        myarth(k+1:MIN(k2,n))=temp+myarth(1:MIN(k,n-k))
        temp=temp+temp
        k=k2
     ENDDO
  ENDIF

END FUNCTION MYARTH

!---------------------------------------------------------------!
!---------------------------------------------------------------!

SUBROUTINE MYTRAPZD(func,a,b,s,n)

  USE sps_vars
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: a,b
  REAL(SP), INTENT(INOUT) :: s
  INTEGER, INTENT(IN) :: n
  INTERFACE
     FUNCTION func(x)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(SIZE(x)) :: func
     END FUNCTION func
  END INTERFACE
  INTERFACE
     FUNCTION myarth(first,increment,n)
       USE sps_vars
       REAL(SP), INTENT(IN) :: first,increment
       INTEGER, INTENT(IN) :: n
       REAL(SP), DIMENSION(n) :: myarth
     END FUNCTION myarth
  END INTERFACE
  REAL(SP) :: del,fsum
  INTEGER :: it

  IF (n == 1) THEN
     s=0.5*(b-a)*SUM(func( (/ a,b /) ))
  ELSE
     it=2**(n-2)
     del=(b-a)/it
     fsum=SUM(func(myarth(a+0.5*del,del,it)))
     s=0.5*(s+del*fsum)
  ENDIF

END SUBROUTINE MYTRAPZD

!---------------------------------------------------------------!
!---------------------------------------------------------------!

SUBROUTINE MYPOLINT(xa,ya,x,y,dy)

  USE sps_vars
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
  REAL(SP), INTENT(IN) :: x
  REAL(SP), INTENT(OUT) :: y,dy
  INTEGER :: m,n,ns
  INTEGER, DIMENSION(1) :: imin
  REAL(SP), DIMENSION(SIZE(xa)) :: c,d,den,ho

  n  = SIZE(xa)
  c  = ya
  d  = ya
  ho = xa-x
  imin = MINLOC(ABS(x-xa))
  ns = imin(1)
  y  = ya(ns)
  ns = ns-1

  DO m=1,n-1
     den(1:n-m)=ho(1:n-m)-ho(1+m:n)
     IF (ANY(den(1:n-m) == 0.0)) &
          WRITE(*,*) 'POLINT ERROR'
     den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
     d(1:n-m)=ho(1+m:n)*den(1:n-m)
     c(1:n-m)=ho(1:n-m)*den(1:n-m)
     IF (2*ns < n-m) THEN
        dy=c(ns+1)
     ELSE
        dy=d(ns)
        ns=ns-1
     ENDIF
     y=y+dy
  ENDDO
  
END SUBROUTINE MYPOLINT

!---------------------------------------------------------------!
!---------------------------------------------------------------!

FUNCTION FUNCINT(func,a,b)

  USE sps_vars
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: a,b
  REAL(SP) :: funcint
  INTERFACE
     FUNCTION func(x)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP), DIMENSION(size(x)) :: func
     END FUNCTION func
  END INTERFACE
  INTERFACE 
     SUBROUTINE MYPOLINT(xa,ya,x,y,dy)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
       REAL(SP), INTENT(IN) :: x
       REAL(SP), INTENT(OUT) :: y,dy
     END SUBROUTINE MYPOLINT
  END INTERFACE
  INTEGER, PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
  REAL(SP), PARAMETER :: EPS=1.0e-7
  REAL(SP), DIMENSION(JMAXP) :: h,s
  REAL(SP) :: dqromb,zero=0.0
  INTEGER :: j

  h(1)=1.0
  DO j=1,JMAX
     CALL mytrapzd(func,a,b,s(j),j)
     IF (j >= K) THEN
        CALL mypolint(h(j-KM:j),s(j-KM:j),zero,funcint,dqromb)
        IF (abs(dqromb) <= EPS*ABS(funcint)) RETURN
     ENDIF
     s(j+1)=s(j)
     h(j+1)=0.25*h(j)
  ENDDO

  WRITE(*,*) 'FUNCINT ERROR:',a,b

END FUNCTION FUNCINT


