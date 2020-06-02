SUBROUTINE GETMAGS(zred,spec,mags,mag_compute)

  !routine to calculate magnitudes in the Vega or AB systems,
  !given an input spectrum and redshift.
  !see parameter compute_vega_mags in sps_vars.f90
  !magnitudes defined in accordance with Fukugita et al. 1996, Eqn 7
  !This routine also redshifts the spectrum, if necessary.

  USE sps_vars; USE sps_utils, ONLY : linterp, tsum
  IMPLICIT NONE

  INTEGER  :: i
  REAL(SP), INTENT(in) :: zred
  REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec
  REAL(SP), INTENT(inout), DIMENSION(nbands) :: mags
  INTEGER, DIMENSION(nbands), INTENT(in), OPTIONAL  :: mag_compute
  INTEGER, DIMENSION(nbands) :: magflag
  REAL(SP), DIMENSION(nspec)  :: tspec
  REAL(SP) :: const, dm

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  mags = 99.
  const= 0.0

  !set up the flags determining which mags are computed
  IF (PRESENT(mag_compute)) THEN
     magflag = mag_compute
     IF (compute_vega_mags.EQ.1) &
          magflag(1) = 1  !force V band to be computed
     IF (MAXVAL(magflag).EQ.0) & !no mags being computed so exit
          RETURN
  ELSE
     magflag = 1
  ENDIF

  !redshift the spectrum
  IF (ABS(zred).GT.tiny_number) THEN
     !write(*,*) "getmags: interpolating"
     DO i=1,nspec
        tspec(i) = MAX(linterp(spec_lambda*(1+zred),spec,&
                       spec_lambda(i)),0.0)
     ENDDO

     !compute additional terms for cosmological mags
     dm    = 5*LOG10(linterp(cosmospl(:,1),cosmospl(:,3),zred)/10.)
     const = dm - 2.5*LOG10(1+zred)

  ELSE

     tspec = spec

  ENDIF

  !integrate over each filter
  DO i=1,nbands
     IF (magflag(i).EQ.0) CYCLE
     mags(i) = TSUM(spec_lambda,tspec*bands(:,i)/spec_lambda)
     IF (mags(i).LE.tiny_number) THEN
        mags(i) = 99.0
     ELSE
        IF (compute_light_ages.EQ.0) THEN
           !the mag2cgs var converts from Lsun/Hz to cgs at 10pc
           mags(i) = -2.5*LOG10(mags(i)) - 48.60 - 2.5*mag2cgs + const
        ENDIF
     ENDIF
  ENDDO

  !put magnitudes in the Vega system if keyword is set
  !(V-band is the first element in the array)
  IF (compute_vega_mags.EQ.1.AND.compute_light_ages.EQ.0) &
       mags(2:nbands) = (mags(2:nbands)-mags(1)) - &
       (magvega(2:nbands)-magvega(1)) + mags(1)


END SUBROUTINE GETMAGS
