SUBROUTINE IMF_WEIGHT(mini,wght,nmass)

  !weight each star by the initial mass function (IMF)
  !such that the total initial population consists of 
  !one solar mass of stars.

  !This weighting scheme assumes that the luminosity, mass, etc.
  !does not vary within the mass bin.  The point is that we
  !want each element to represent the whole bin, from 
  !mass+/-0.5dm, rather than just the values at point i.
  !Then every intergral over mass is just a sum.

  USE sps_vars; USE sps_utils, ONLY : imf
  USE nr, ONLY : qromb
  IMPLICIT NONE

  REAL(SP), INTENT(inout), DIMENSION(nm) :: wght
  REAL(SP), INTENT(in), DIMENSION(nm)    :: mini
  INTEGER, INTENT(in) :: nmass
  INTEGER  :: i
  REAL(SP) :: m1,m2, i1,i2
  REAL(SP), DIMENSION(1) :: tm1,tm2,ti1,ti2

  !--------------------------------------------------------!
  !--------------------------------------------------------!

  wght = 0.0

  DO i=1,nmass

     IF (mini(i).LT.imf_lower_limit.OR.&
          mini(i).GT.imf_upper_limit) CYCLE

     IF (i.EQ.1) THEN
        m1 = imf_lower_limit
     ELSE
        m1 = mini(i) - 0.5*(mini(i)-mini(i-1))
     ENDIF
     IF (i.EQ.nmass) THEN
        m2 = mini(i)
     ELSE
        m2 = mini(i) + 0.5*(mini(i+1)-mini(i))
     ENDIF

     wght(i) = qromb(imf,m1,m2)

     !tm1 = m1
     !tm2 = m2
     !ti1 = imf(tm1)
     !ti2 = imf(tm2)
     !write(*,*) m1,m2,wght(i)/ ((m2-m1)*(ti1+ti2)/2.)-1

  ENDDO

  !normalize the weights
  imf_type = imf_type + 10
  wght = wght / qromb(imf,imf_lower_limit,imf_upper_limit)
  imf_type = imf_type - 10

  RETURN
END SUBROUTINE IMF_WEIGHT

