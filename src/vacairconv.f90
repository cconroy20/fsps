FUNCTION AIRTOVAC(lam)
  
  !convert wavelengths from air to vac
  !see Morton (1991 Ap.J. Suppl. 77, 119)
  !this code was adapted from the IDL routine airtovac.pro

  USE sps_vars; USE sps_utils, ONLY : locate
  IMPLICIT NONE

  INTEGER :: vv,nn
  REAL(SP), DIMENSION(:), INTENT(in) :: lam
  REAL(SP), DIMENSION(SIZE(lam)) :: airtovac,fact,sigma2

  !------------------------------------------------------!

  nn = SIZE(lam)

  !Convert to wavenumber squared
  sigma2(1:nn) = (1E4/lam(1:nn) )**2.   

  !Compute conversion factor
  fact(1:nn) = 1.D0 + 6.4328D-5 + 2.94981D-2/(146.d0 - sigma2(1:nn)) + &
       2.5540D-4/( 41.d0 - sigma2(1:nn))
  
  !no conversion for wavelengths <2000A
  IF (lam(1).LT.2000.D0) THEN
     vv = MIN(MAX(locate(lam(1:nn),2000.D0),1),nn)
     fact(1:vv) = 1.0
  ENDIF

  !Convert wavelength
  airtovac(1:nn) = lam(1:nn)*fact(1:nn)

END FUNCTION AIRTOVAC

!---------------------------------------------------------------!
!---------------------------------------------------------------!
!---------------------------------------------------------------!

FUNCTION VACTOAIR(lam)

  !convert wavelengths from air to vac
  !see Morton (1991 Ap.J. Suppl. 77, 119)
  !this code was adapted from the IDL routine vactoair.pro

  USE sps_vars; USE sps_utils, ONLY : locate
  IMPLICIT NONE

  INTEGER :: vv,nn
  REAL(SP), DIMENSION(:), INTENT(in) :: lam
  REAL(SP), DIMENSION(SIZE(lam)) :: vactoair,fact

  !------------------------------------------------------!

  nn = SIZE(lam)

  fact(1:nn) = 1.D0 + 2.735182D-4 + 131.4182/lam(1:nn)**2 + &
       2.76249D8/lam(1:nn)**4

  !no conversion for wavelengths <2000A
   IF (lam(1).LT.2000.D0) THEN
     vv = MIN(MAX(locate(lam(1:nn),2000.D0),1),nn)
     fact(1:vv) = 1.0
  ENDIF

  !Convert wavelengths
  vactoair(1:nn) = lam(1:nn)/fact(1:nn)

END FUNCTION VACTOAIR
