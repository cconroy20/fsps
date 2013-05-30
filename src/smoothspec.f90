SUBROUTINE SMOOTHSPEC(lambda,spec,sigma)

  !routine to compute velocity broadening of an input spectrum
  !the PSF kernel has a width of m*sigma, where m=6

  USE sps_vars; USE nr, ONLY : locate
  USE sps_utils, ONLY : linterp
  IMPLICIT NONE
  
  REAL(SP), INTENT(in), DIMENSION(nspec) :: lambda
  REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec
  REAL(SP), INTENT(in) :: sigma
  REAL(SP), DIMENSION(nspec) :: tspec,tnspec,vel,func,gauss,psf
  REAL(SP) :: cg,xmax,xmin,fwhm,psig
  INTEGER :: i,j,il,ih,m=6,grange

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  IF (sigma.GT.0) THEN
  
     DO i=1,nspec
        tspec(i) = linterp(LOG(lambda(1:nspec)),spec(1:nspec),lnlam(i))
     ENDDO
     
     fwhm   = sigma*2.35482/clight*1E5/dlstep
     psig   = fwhm/2.d0/SQRT(-2.d0*LOG(0.5d0)) ! equivalent sigma for kernel
     grange = FLOOR(m*psig)	               ! range for kernel (-range:range)
     DO i=1,2*grange+1
        psf(i) = 1.d0/SQRT(2.d0*mypi)/psig*EXP(-((i-grange-1)/psig)**2/2.d0)
     ENDDO
     psf(1:2*grange+1) = psf(1:2*grange+1) / SUM(psf(1:2*grange+1))
     
     DO i=grange+1,nspec-grange
        tnspec(i) = SUM( psf(1:2*grange+1)*tspec(i-grange:i+grange) )
     ENDDO
     
     !interpolate back to the main array
     DO i=1,nspec
        spec(i) = linterp(EXP(lnlam),tnspec,lambda(i))
     ENDDO

  ENDIF

END SUBROUTINE SMOOTHSPEC
