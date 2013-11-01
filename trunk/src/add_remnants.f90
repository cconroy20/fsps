SUBROUTINE ADD_REMNANTS(mass,maxmass)

  !add remnant (WD, NS, BH) masses back into the total mass
  !of the SSP.  These initial-mass-dependent remnant 
  !formulae are taken from Renzini & Ciotti 1993.
  
  USE sps_vars
  USE sps_utils, ONLY : imf, funcint
  IMPLICIT NONE

  REAL(SP), INTENT(inout) :: mass
  !maximum mass still alive
  REAL(SP), INTENT(in) :: maxmass
  REAL(SP) :: minmass, imfnorm

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  !normalize the weights
  imf_type = imf_type+10
  imfnorm  = funcint(imf,imf_lower_limit,imf_upper_limit)
  imf_type = imf_type-10

  !BH remnants
  !40<M<imf_up leave behind a 0.5*M BH
  minmass = MAXVAL((/mlim_bh,maxmass/))
  imf_type = imf_type+10
  mass = mass + 0.5*funcint(imf,minmass,imf_upper_limit)/imfnorm
  imf_type = imf_type-10

  !NS remnants
  !8.5<M<40 leave behind 1.4 Msun NS
  IF (maxmass.LE.mlim_bh) THEN
     minmass = MAXVAL((/mlim_ns,maxmass/))
     mass = mass + 1.4*funcint(imf,minmass,mlim_bh)/imfnorm
  ENDIF

  !WD remnants
  !M<8.5 leave behind 0.077*M+0.48 WD
  IF (maxmass.LE.8.5) THEN
     mass = mass + 0.48*funcint(imf,maxmass,mlim_ns)/imfnorm
     imf_type = imf_type+10
     mass = mass + 0.077*funcint(imf,maxmass,mlim_ns)/imfnorm
     imf_type = imf_type-10

  ENDIF

  RETURN

END SUBROUTINE ADD_REMNANTS
