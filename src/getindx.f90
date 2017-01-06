FUNCTION INTIND(lam,func,lo,hi)

  !perform integral over spectrum for index computation

  USE sps_vars
  USE sps_utils, ONLY : tsum, locate
  IMPLICIT NONE

  INTEGER :: l1,l2,i
  REAL(SP), INTENT(in), DIMENSION(nspec) :: lam, func
  REAL(SP), INTENT(in) :: lo,hi
  REAL(SP) :: f1,f2,intind

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !take care of the ends
  l1 = MAX(MIN(locate(lam,lo),nspec-1),1)
  l2 = MAX(MIN(locate(lam,hi),nspec-1),2)
  f1 = (func(l1+1)-func(l1))/(lam(l1+1)-lam(l1))*&
       (lo-lam(l1))+func(l1)
  f2 = (func(l2+1)-func(l2))/(lam(l2+1)-lam(l2))*&
       (hi-lam(l2))+func(l2)

  IF (l1.EQ.l2) THEN
     intind = (f2+f1)/2.*(hi-lo)
  ELSE
     intind = TSUM(lam(l1+1:l2),func(l1+1:l2))
     intind = intind + (lam(l1+1)-lo)*(f1+func(l1+1))/2.
     intind = intind + (hi-lam(l2))*(f2+func(l2))/2.
  ENDIF

END FUNCTION INTIND

!------------------------------------------------------------!
!------------------------------------------------------------!
!------------------------------------------------------------!

SUBROUTINE GETINDX(lambda,spec,indices)

  !routine to calculate indices from an input spectrum
  !indices are defined in fsps/data/allindices.dat

  USE sps_vars
  USE sps_utils, ONLY : intind, locate
  IMPLICIT NONE

  INTEGER :: j
  REAL(SP), INTENT(in), DIMENSION(nspec) :: spec,lambda
  REAL(SP), INTENT(inout), DIMENSION(nindx) :: indices
  REAL(SP) :: intfifc,cb,cr,lr,lb

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  indices = 999.

  DO j=1,nindx

     !blue continuum
     cb = intind(lambda,spec,indexdefined(3,j),indexdefined(4,j))
     cb = cb / (indexdefined(4,j)-indexdefined(3,j))     
     lb = (indexdefined(3,j)+indexdefined(4,j))/2.

     IF (indexdefined(7,j).NE.4) THEN 

        !red continuum 
        cr = intind(lambda,spec,indexdefined(5,j),indexdefined(6,j))
        cr = cr / (indexdefined(6,j)-indexdefined(5,j))
        lr = (indexdefined(5,j)+indexdefined(6,j))/2.
        
        !compute integral(fi/fc)
        !NB: fc here is a linear interpolation between the red and blue.
        intfifc = intind(lambda,spec/((cr-cb)/(lr-lb)*(lambda-lb) + cb),&
             indexdefined(1,j),indexdefined(2,j))
     
     ELSE
        intfifc = intind(lambda,spec,indexdefined(1,j),indexdefined(2,j))
        intfifc = intfifc/(indexdefined(2,j)-indexdefined(1,j))
     ENDIF

     IF (indexdefined(7,j).EQ.1.) THEN
        !compute magnitude
        indices(j) = &
             -2.5*LOG10(intfifc/(indexdefined(2,j)-indexdefined(1,j)))
     ELSE IF (indexdefined(7,j).EQ.2.) THEN
        !compute EW (in Ang)
        indices(j) = (indexdefined(2,j)-indexdefined(1,j)) - intfifc
     ELSE IF (indexdefined(7,j).EQ.3.) THEN
        !compute Dn4000
        !NB: this only works with cr and cb computed above b/c 
        !the wavelength intervals for blue and red are the 
        !same for these indices, otherwise a slightly 
        !different cr and cb would have to be computed
        indices(j) = cr/cb
     ELSE IF (indexdefined(7,j).EQ.4.) THEN
        !compute magnitude from a flux ratio
        indices(j) = -2.5*LOG10(intfifc/cb)
     ENDIF

     !set dummy values for indices defined off of the wavelength grid
     IF (indexdefined(6,j).GT.lambda(nspec)) indices(j) = 999.0
     IF (indexdefined(3,j).LT.lambda(1)) indices(j) = 999.0

  ENDDO

END SUBROUTINE GETINDX
