SUBROUTINE SMOOTHSPEC(lambda,spec,sigma,minl,maxl,ires)

  !routine to compute velocity broadening of an input spectrum
  !this is only approximate b/c we are ignoring the variation
  !in dlambda with lambda within each integration
  !integration is truncated at +/-4*sigma.
  !If optional input ires is present, then the spectrum will be
  !smoothed by a wavelength dependent velocity dispersion.

  USE sps_vars; USE sps_utils, ONLY : locate,linterp,tsum,linterparr
  IMPLICIT NONE
  
  REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec
  REAL(SP), INTENT(in), DIMENSION(nspec)    :: lambda
  REAL(SP), INTENT(in), DIMENSION(nspec), OPTIONAL :: ires
  REAL(SP), INTENT(in) :: sigma,minl,maxl
  REAL(SP), DIMENSION(nspec) :: tspec,tnspec,vel,func,gauss,psf,lnlam
  REAL(SP) :: ckms,cg,xmax,xmin,fwhm,psig,dlstep,sigmal
  INTEGER :: i,j,il,ih,m=4,grange

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
 
  IF (sigma.LE.tiny_number) RETURN

  ckms = clight/1E13

  tspec = spec

  !convolve at fixed sigma_velocity
  IF (smooth_velocity.EQ.1) THEN

     !compute smoothing the fast (and slightly less accurate) way
     IF (smoothspec_fast.EQ.1.OR.PRESENT(ires)) THEN
        
        spec = 0.0

        DO i=1,nspec
           
           IF (lambda(i).LT.minl.OR.lambda(i).GT.maxl) THEN
              spec(i)=tspec(i)
              CYCLE
           ENDIF
           
           IF (PRESENT(ires)) THEN
              sigmal = ires(i)
              IF (sigmal.LE.tiny_number) THEN
                 spec(i)=tspec(i)
                 CYCLE
              ENDIF
           ELSE
              sigmal = sigma
           ENDIF

           xmax = lambda(i)*(m*sigmal/ckms+1)
           ih   = MIN(locate(lambda(1:nspec),xmax),nspec)
           il   = MAX(2*i-ih,1)
           
           IF (il.EQ.ih) THEN
              spec(i) = tspec(i)
           ELSE
              vel(il:ih)  = (lambda(i)/lambda(il:ih)-1)*ckms
              func(il:ih) =  1/SQRT(2*mypi)/sigmal * &
                   EXP(-vel(il:ih)**2/2./sigmal**2)
              !normalize the weights to integrate to unity
              func(il:ih) = func(il:ih) / TSUM(vel(il:ih),func(il:ih))
              spec(i) = TSUM(vel(il:ih),func(il:ih)*tspec(il:ih))
          ENDIF
           
        ENDDO

     !compute smoothing the correct (and somewhat slower) way
     !NB: the accuracy of this approach depends strongly on the 
     !min/max wavelength parameters through the density of the lnlam grid
     ELSE
        
        dlstep = (LOG(maxl)-LOG(minl))/nspec
        DO i=1,nspec
           lnlam(i) = i*dlstep+LOG(minl)
        ENDDO
        
        tspec = linterparr(LOG(lambda(1:nspec)),spec(1:nspec),lnlam)
        
        fwhm   = sigma*2.35482/ckms/dlstep
        psig   = fwhm/2.0/SQRT(-2.0*LOG(0.5)) ! equivalent sigma for kernel
        grange = FLOOR(m*psig)	              ! range for kernel (-range:range)
        
        DO i=1,2*grange+1
           psf(i) = 1.d0/SQRT(2.*mypi)/psig*EXP(-((i-grange-1)/psig)**2/2.)
        ENDDO
        psf(1:2*grange+1) = psf(1:2*grange+1) / SUM(psf(1:2*grange+1))
        
        DO i=grange+1,nspec-grange
           tnspec(i) = SUM( psf(1:2*grange+1)*tspec(i-grange:i+grange) )
        ENDDO
        
        !the ends are not smoothed
        DO i=1,grange
           tnspec(i)=tspec(i)
        ENDDO
        DO i=nspec-grange+1,nspec
           tnspec(i)=tspec(i)
        ENDDO
        
        !interpolate back to the main array
        il = locate(lambda,minl)
        ih = locate(lambda,maxl)
        spec(il:ih) = linterparr(EXP(lnlam),tnspec,lambda(il:ih))
        
     ENDIF

  !convolve at fixed sigma_wavelength
  ELSE
     
     DO i=1,nspec

        IF (lambda(i).LT.minl.OR.lambda(i).GT.maxl) THEN
           spec(i)=tspec(i)
           CYCLE
        ENDIF

        xmax = lambda(i)*(m*sigma+1)
        ih   = MIN(locate(lambda(1:nspec),xmax),nspec)
        il   = MAX(2*i-ih,1)

        IF (il.EQ.ih) THEN
           spec(i) = tspec(i)
        ELSE
           func(il:ih) =  1/SQRT(2*mypi)/sigma * &
                EXP(-(lambda(il:ih)-lambda(i))**2/2./sigma**2)
           !normalize the weights to integrate to unity
           func(il:ih) = func(il:ih) / TSUM(lambda(il:ih),func(il:ih))
           spec(i) = TSUM(lambda(il:ih),func(il:ih)*tspec(il:ih))
        ENDIF

        !gauss = 1/SQRT(2*mypi)/sigma*EXP(-(lambda-lambda(i))**2/2/sigma**2)
        !DO j=1,nspec
        !   spec(i) = spec(i) + gauss(j)*tspec(j)
        !ENDDO
 
    ENDDO
     
  ENDIF

  RETURN

END SUBROUTINE SMOOTHSPEC
