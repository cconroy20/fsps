SUBROUTINE GETMAGS(zred,spec,mags)

  !routine to calculate magnitudes in the Vega or AB systems,
  !given an input spectrum and redshift.
  !see parameter compute_vega_mags in sps_vars.f90

  USE sps_vars 
  USE nrtype; USE nrutil, ONLY : assert_eq
  USE nr, ONLY : spline, splint
  IMPLICIT NONE

  INTEGER :: i,j
  REAL(SP), INTENT(in), DIMENSION(nspec) :: spec
  REAL(SP), INTENT(in) :: zred
  REAL(SP), DIMENSION(nspec) :: tspec
  REAL(SP), DIMENSION(nbands) :: mags, vmags
  REAL(SP), DIMENSION(nspec)  :: tmpspl

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  !redshift the spectra
  IF (zred.NE.0.0) THEN
     tmpspl = 0.0
     CALL spline(spec_lambda*(1+zred),spec,1.0e30_sp,1.0e30_sp,tmpspl)
     DO i=1,nspec
        tspec(i) = splint(spec_lambda*(1+zred),spec,&
             tmpspl,spec_lambda(i))
        IF (tspec(i).LT.0.0) tspec(i)=0.
     ENDDO
     !the units of the spectra are Lsun/Hz; convert to
     !erg/s/cm^2/Hz, at 10pc for absolute mags
     tspec = tspec * lsun / 4.0 / mypi / (pc2cm*pc2cm) / 100.0
  ELSE
     tspec = spec  * lsun / 4.0 / mypi / (pc2cm*pc2cm) / 100.0
  ENDIF

  !integrate over each band-pass using trapezoidal rule.  
  !more sophisticated integration using splines yields
  !similar results to ~0.001 mags in FUV-K bands.
  mags  = 0.0
  vmags = 0.0
  DO i=1,nbands
     DO j=1,nspec-1
        mags(i) = mags(i) + (spec_lambda(j+1)-spec_lambda(j)) * &
             (tspec(j+1)*bands(i,j+1)/spec_lambda(j+1)+&
             tspec(j)*bands(i,j)/spec_lambda(j))/2
        IF (compute_vega_mags.EQ.1) &
             vmags(i) = vmags(i) + (spec_lambda(j+1)-spec_lambda(j)) * &
             (vega_spec(j+1)*bands(i,j+1)/spec_lambda(j+1)+&
             vega_spec(j)*bands(i,j)/spec_lambda(j))/2        
     ENDDO
     IF (mags(i).LE.1E-35)  mags(i)  = 1E-35
     IF (vmags(i).LE.1E-35) vmags(i) = 1E-35
  ENDDO

  !define AB magnitudes; see Fukugita et al. 1996, Oke & Gunn
  mags  = -2.5 * LOG10(mags)  - 48.60
  vmags = -2.5 * LOG10(vmags) - 48.60

  IF (compute_vega_mags.EQ.1) THEN
     !set zero points for Vega magnitudes
     !V-band is the first element in the array
     DO i=2,nbands
        mags(i) = (mags(i)-mags(1)) - (vmags(i)-vmags(1)) + &
             mags(1)
     ENDDO
  ENDIF
     
END SUBROUTINE GETMAGS
