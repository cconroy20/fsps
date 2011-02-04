FUNCTION splint(xa,ya,y2a,x)
  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
  USE nr, ONLY: locate
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
  REAL(SP), INTENT(IN) :: x
  REAL(SP) :: splint
  INTEGER(I4B) :: khi,klo,n
  REAL(SP) :: a,b,h

  n=assert_eq(size(xa),size(ya),size(y2a),'splint')
  klo=max(min(locate(xa,x),n-1),1)
  khi=klo+1
  h=xa(khi)-xa(klo)
  if (h == 0.0) call nrerror('bad xa input in splint')
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+&
       (b**3-b)*y2a(khi))*(h**2)/6.0_sp

END FUNCTION splint
