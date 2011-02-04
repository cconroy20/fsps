SUBROUTINE spline(x,y,yp1,ypn,y2)
  USE nrtype; USE nrutil, ONLY : assert_eq
  USE nr, ONLY : tridag
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
  REAL(SP), INTENT(IN) :: yp1,ypn
  REAL(SP), DIMENSION(:), INTENT(OUT) :: y2
  INTEGER(I4B) :: n
  REAL(SP), DIMENSION(size(x)) :: a,b,c,r

  n=assert_eq(size(x),size(y),size(y2),'spline')
  c(1:n-1)=x(2:n)-x(1:n-1)
  r(1:n-1)=6.0_sp*((y(2:n)-y(1:n-1))/c(1:n-1))
  r(2:n-1)=r(2:n-1)-r(1:n-2)
  a(2:n-1)=c(1:n-2)
  b(2:n-1)=2.0_sp*(c(2:n-1)+a(2:n-1))
  b(1)=1.0
  b(n)=1.0
  if (yp1 > 0.99e30_sp) then
     r(1)=0.0
     c(1)=0.0
  else
     r(1)=(3.0_sp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
     c(1)=0.5
  end if
  if (ypn > 0.99e30_sp) then
     r(n)=0.0
     a(n)=0.0
  else
     r(n)=(-3.0_sp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
     a(n)=0.5
  end if
  call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
END SUBROUTINE spline
