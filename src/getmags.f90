SUBROUTINE GETMAGS(zred,spec,mags)

  !routine to calculate magnitudes in the Vega or AB systems,
  !given an input spectrum and redshift.
  !see parameter compute_vega_mags in sps_vars.f90
  !magnitudes defined in accordance with Fukugita et al. 1996, Eqn 7

  USE sps_vars; USE sps_utils, ONLY : linterp
  IMPLICIT NONE

  INTEGER  :: i
  REAL(SP) :: d
  REAL(SP), INTENT(in), DIMENSION(nspec) :: spec
  REAL(SP), INTENT(in) :: zred
  REAL(SP), INTENT(inout), DIMENSION(nbands) :: mags
  REAL(SP), DIMENSION(nspec)  :: tmpspl, tspec
  REAL(SP), DIMENSION(14) :: lami
  INTEGER,  DIMENSION(14) :: ind

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  !redshift the spectra
  IF (zred.NE.0.0) THEN
     tmpspl = 0.0
     DO i=1,nspec
        tspec(i) = MAX(linterp(spec_lambda*(1+zred),spec,&
        spec_lambda(i)),0.0)
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

  !normalize the IRAC, PACS, SPIRE, and IRAS photometry to nu*fnu=const
  lami = (/3.550,4.493,5.731,7.872,70.0,100.0,160.0,250.0,350.0,500.0,&
       12.0,25.0,60.0,100.0/)*1E4
  DO i=1,14
     ind=(/19,20,21,22,71,72,73,74,75,76,77,78,79,80/)
     d = SUM( (spec_lambda(2:nspec)-spec_lambda(1:nspec-1)) * &
          ((spec_lambda(2:nspec)/lami(i))**(-1.0)*bands(ind(i),2:nspec)/&
          spec_lambda(2:nspec)+(spec_lambda(1:nspec-1)/lami(i))**(-1.0)*&
          bands(ind(i),1:nspec-1)/spec_lambda(1:nspec-1))/2. )
     mags(ind(i)) = mags(ind(i)) / d
  ENDDO

  !normalize the MIPS photometry to a BB (beta=2)
  lami(1:3) = (/23.68,71.42,155.9/)*1E4
  DO i=1,3
     ind(1:3) = (/66,67,68/)
     d = SUM( (spec_lambda(2:nspec)-spec_lambda(1:nspec-1)) * &
          ((spec_lambda(2:nspec)/lami(i))**(-2.0)*bands(ind(i),2:nspec)/&
          spec_lambda(2:nspec)+(spec_lambda(1:nspec-1)/lami(i))**(-1.0)*&
          bands(ind(i),1:nspec-1)/spec_lambda(1:nspec-1))/2. )
     mags(ind(i)) = mags(ind(i)) / d
  ENDDO

  !convert to magnitudes in AB system
  mags  = -2.5*LOG10(mags)-48.60

  !put magnitudes in the Vega system if keyword is set
  !(V-band is the first element in the array)
  IF (compute_vega_mags.EQ.1) &
       mags(2:nbands) = (mags(2:nbands)-mags(1)) - &
       (magvega(2:nbands)-magvega(1)) + mags(1)
     
END SUBROUTINE GETMAGS
