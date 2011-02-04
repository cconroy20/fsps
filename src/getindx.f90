REAL FUNCTION INTIND(lam,func,lo,hi)

  !perform integral over spectrum, for index computation
  !uses trapezoidal rule

  USE sps_vars; USE nrtype; USE nr, ONLY : locate
  IMPLICIT NONE
  INTEGER :: l1,l2,i,n
  REAL(SP), INTENT(in), DIMENSION(:) :: lam,func
  REAL(SP), INTENT(in) :: lo,hi
  REAL(SP) :: f1,f2

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  l1 = locate(lam,lo)
  l2 = locate(lam,hi)
  f1 = (func(l1+1)-func(l1))/(lam(l1+1)-lam(l1))*&
       (lo-lam(l1)) + func(l1)
  f2 = (func(l2+1)-func(l2))/(lam(l2+1)-lam(l2))*&
       (hi-lam(l2)) + func(l2)

  IF (l1.EQ.l2) THEN
     intind = (f2+f1)/2.*(hi-lo)
  ELSE
     intind = 0.0
     DO i=l1+1,l2-1
        intind = intind + (lam(i+1)-lam(i)) * &
             (func(i+1) + func(i))/2.
     ENDDO
     intind = intind + (lam(l1+1)-lo)*(f1+func(l1+1))/2.
     intind = intind + (hi-lam(l2))*(f2+func(l2))/2.
  ENDIF

END FUNCTION INTIND

!------------------------------------------------------------!
!------------------------------------------------------------!

SUBROUTINE GETINDX(lambda,spec,indices)

  !routine to calculate indices from an input spectrum
  !bandpasses and units are taken from the SDSS MPA reductions website:
  !http://www.mpa-garching.mpg.de/SDSS/DR4/SDSS_fixedidx.html
  !the input spectrum is assumed to be in Fnu

  USE nrtype; USE nrutil, ONLY : assert_eq; USE sps_vars
  USE sps_utils
  USE nr, ONLY : spline, splint,locate
  IMPLICIT NONE

  INTEGER :: i,j
  REAL(SP), INTENT(in), DIMENSION(:) :: spec,lambda
  REAL(SP), INTENT(inout), DIMENSION(nindsps) :: indices
  REAL(SP) :: intfifc,cb,cr,lr,lb

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  indices = -99.

  IF (spec_type.NE.'miles') & 
       WRITE(*,*) 'GETINDX WARNING: you are attempting to '//&
       'compute indices on the BaSeL library; are '//&
       'you sure you want to do that??'

  DO j=1,nindsps

     !red continuum 
     cr = intind(lambda,spec,indexdefined(5,j),indexdefined(6,j))
     cr = cr / (indexdefined(6,j)-indexdefined(5,j))
     lr = (indexdefined(5,j)+indexdefined(6,j))/2.
     
     !blue continuum
     cb = intind(lambda,spec,indexdefined(3,j),indexdefined(4,j))
     cb = cb / (indexdefined(4,j)-indexdefined(3,j))     
     lb = (indexdefined(3,j)+indexdefined(4,j))/2.
     
     !compute integral(fi/fc)
     !NB: fc here is a linear interpolation between the red and blue.
     intfifc = intind(lambda,spec/((cr-cb)/(lr-lb)*(lambda-lb) + cb),&
          indexdefined(1,j),indexdefined(2,j))
     
     IF (indexdefined(7,j).EQ.1.) THEN
        !compute magnitude
        indices(j) = -2.5* &
             log10(1/(indexdefined(2,j)-indexdefined(1,j))*intfifc)
     ELSE IF (indexdefined(7,j).EQ.2.) THEN
        !compute EW (in Ang)
        indices(j) = (indexdefined(2,j)-indexdefined(1,j)) - intfifc
     ELSE IF (indexdefined(7,j).EQ.3.) THEN
        !compute Dn4000 and CO indices
        !NB: this only works with cr and cb computed above b/c 
        !the wavelength intervals for blue and red are the 
        !same for these indices, otherwise a slightly 
        !different cr and cb would have to be computed
        indices(j) = cr/cb
     ENDIF

  ENDDO

END SUBROUTINE GETINDX
