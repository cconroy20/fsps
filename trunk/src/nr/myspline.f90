SUBROUTINE myspline(x,y,y2)

  USE nrtype; USE nr, ONLY : tridag
  USE sps_vars
  IMPLICIT NONE
  REAL(SP), DIMENSION(nz), INTENT(IN) :: x,y
  REAL(SP), DIMENSION(nz), INTENT(OUT) :: y2
  REAL(SP), DIMENSION(nz) :: a,b,c,r

  c(1:nz-1)=x(2:nz)-x(1:nz-1)
  r(1:nz-1)=6.0*((y(2:nz)-y(1:nz-1))/c(1:nz-1))
  r(2:nz-1)=r(2:nz-1)-r(1:nz-2)
  a(2:nz-1)=c(1:nz-2)
  b(2:nz-1)=2.0*(c(2:nz-1)+a(2:nz-1))

  b(1) =1.0
  b(nz)=1.0
  r(1) =0.0
  c(1) =0.0
  r(nz)=0.0
  a(nz)=0.0

  call tridag(a(2:nz),b(1:nz),c(1:nz-1),r(1:nz),y2(1:nz))

END SUBROUTINE myspline 
