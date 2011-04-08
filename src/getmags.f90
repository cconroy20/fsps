SUBROUTINE GETMAGS(zred,spec,mags)

  !routine to calculate magnitudes in the Vega or AB systems,
  !given an input spectrum and redshift.
  !see parameter compute_vega_mags in sps_vars.f90
  !magnitudes defined in accordance with Fukugita et al. 1996, Eqn 7

  USE sps_vars 
  USE nrtype; USE nrutil, ONLY : assert_eq
  USE nr, ONLY : spline, splint
  IMPLICIT NONE

  INTEGER :: i
  REAL(SP), INTENT(in), DIMENSION(nspec) :: spec
  REAL(SP), INTENT(in) :: zred
  REAL(SP), INTENT(inout), DIMENSION(nbands) :: mags
  REAL(SP), DIMENSION(nspec)  :: tmpspl, tspec

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  !redshift the spectra
  IF (zred.NE.0.0) THEN
     tmpspl = 0.0
     CALL SPLINE(spec_lambda*(1+zred),spec,1.0e30_sp,1.0e30_sp,tmpspl)
     DO i=1,nspec
        tspec(i) = splint(spec_lambda*(1+zred),spec,tmpspl,spec_lambda(i))
        tspec(i) = MAX(tspec(i),0.0)
     ENDDO
  ELSE
     tspec = spec  
  ENDIF

  !the units of the spectra are Lsun/Hz; convert to
  !erg/s/cm^2/Hz, at 10pc for absolute mags
  tspec = tspec*lsun/4.0 /mypi/(pc2cm*pc2cm)/100.0

  !integrate over each filter
  DO i=1,nbands
     mags(i) = SUM( (spec_lambda(2:nspec)-spec_lambda(1:nspec-1)) * &
          (tspec(2:nspec)*bands(i,2:nspec)/spec_lambda(2:nspec)+&
          tspec(1:nspec-1)*bands(i,1:nspec-1)/spec_lambda(1:nspec-1))/2. )
     mags(i)  = MAX(mags(i),tiny_number)     
  ENDDO

  mags  = -2.5*LOG10(mags)-48.60

  !put magnitudes in the Vega system
  !(V-band is the first element in the array)
  IF (compute_vega_mags.EQ.1) &
       mags(2:nbands) = (mags(2:nbands)-mags(1)) - &
       (magvega(2:nbands)-magvega(1)) + mags(1)
     
END SUBROUTINE GETMAGS
