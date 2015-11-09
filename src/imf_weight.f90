SUBROUTINE IMF_WEIGHT(mini,wght,nmass)

  !weight each star by the initial mass function (IMF)
  !such that the total initial population consists of 
  !one solar mass of stars.

  !This weighting scheme assumes that the luminosity, mass, etc.
  !does not vary within the mass bin.  The point is that we
  !want each element to represent the whole bin, from 
  !mass+/-0.5dm, rather than just the values at point i.
  !Then every intergral over mass is just a sum.

  USE sps_vars; USE sps_utils, ONLY : imf, funcint
  IMPLICIT NONE

  REAL(SP), INTENT(inout), DIMENSION(nm) :: wght
  REAL(SP), INTENT(in), DIMENSION(nm)    :: mini
  INTEGER, INTENT(in) :: nmass
  INTEGER  :: i
  REAL(SP) :: m1,m2

  !--------------------------------------------------------!
  !--------------------------------------------------------!

  wght = 0.0

  DO i=1,nmass

     IF (mini(i).LT.imf_lower_limit.OR.&
          mini(i).GT.imf_upper_limit) CYCLE

     IF (i.EQ.1) THEN
        !note that this is not equal to imf_lower_limit
        !only for the Geneva models, which do not extend below 1.0 Msun
        m1 = imf_lower_bound
     ELSE
        m1 = mini(i) - 0.5*(mini(i)-mini(i-1))
     ENDIF
     IF (i.EQ.nmass) THEN
        m2 = mini(i)
     ELSE
        m2 = mini(i) + 0.5*(mini(i+1)-mini(i))
     ENDIF

     IF (m2.LT.m1) THEN
        WRITE(*,*) 'IMF_WEIGHT WARNING: non-monotonic mass!',m1,m2,m2-m1
        CYCLE
     ENDIF

     IF (m2.EQ.m1) CYCLE

     wght(i) = funcint(imf,m1,m2)

  ENDDO

  !normalize the weights as an integral from lower to upper limits
  imf_type = imf_type + 10
  wght = wght / funcint(imf,imf_lower_limit,imf_upper_limit)
  imf_type = imf_type - 10

  RETURN

END SUBROUTINE IMF_WEIGHT

