SUBROUTINE ADD_REMNANTS(mass,maxmass)

  !add remnant (WD, NS, BH) masses back into the total mass
  !of the SSP.  These initial-mass-dependent remnant 
  !formulae are taken from Renzini & Ciotti 1993.
  
  USE sps_vars; USE nrtype
  USE nr, ONLY : qromb
  USE sps_utils, ONLY : imf
  IMPLICIT NONE

  REAL(SP), INTENT(inout) :: mass
  !maximum mass still in the isochrone
  REAL(SP), INTENT(in) :: maxmass
  REAL(SP) :: minmass, imfnorm

  !normalize the weights
  imf_type = imf_type + 10
  imfnorm  = qromb(imf,mlo,mup)
  imf_type = imf_type - 10

  !BH remnants
  minmass = MAXVAL((/40.,maxmass/))
  mass = mass + 0.48*qromb(imf,minmass,100.)/imfnorm
  imf_type = imf_type+10
  mass = mass + 0.077*qromb(imf,minmass,100.)/imfnorm
  imf_type = imf_type-10

  !NS remnants
  IF (maxmass.LE.40) THEN
     minmass = MAXVAL((/8.5,maxmass/))
     mass = mass + 1.4*qromb(imf,minmass,40.)/imfnorm
  ENDIF

  !WD remnants
  IF (maxmass.LE.8.5) THEN
     mass = mass + 0.5*qromb(imf,maxmass,8.5)/imfnorm
  ENDIF

  RETURN
END SUBROUTINE ADD_REMNANTS
