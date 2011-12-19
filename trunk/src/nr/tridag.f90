SUBROUTINE tridag(a,b,c,r,u)
  
  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
  REAL(SP), DIMENSION(:), INTENT(OUT) :: u
  REAL(SP), DIMENSION(size(b)) :: gam
  INTEGER(I4B) :: n,j
  REAL(SP) :: bet
  n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser')
  bet=b(1)
  if (bet == 0.0) call nrerror('tridag_ser: Error at code stage 1')
  u(1)=r(1)/bet
  do j=2,n
     gam(j)=c(j-1)/bet
     bet=b(j)-a(j-1)*gam(j)
     if (bet == 0.0) &
          call nrerror('tridag_ser: Error at code stage 2')
     u(j)=(r(j)-a(j-1)*u(j-1))/bet
  end do
  do j=n-1,1,-1
     u(j)=u(j)-gam(j+1)*u(j+1)
  end do

END SUBROUTINE tridag
