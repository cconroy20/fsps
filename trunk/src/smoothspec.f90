SUBROUTINE SMOOTHSPEC(lambda,spec,sigma,minl,maxl)

  !routine to compute velocity broadening of an input spectrum
  !this is only approximate b/c we are ignoring the variation
  !in dlambda with lambda within each integration
  !integration is truncated at +/-4*sigma

  USE sps_vars; USE sps_utils, ONLY : locate,linterp,tsum,linterparr
  IMPLICIT NONE
  
  REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec
  REAL(SP), INTENT(in), DIMENSION(nspec) :: lambda
  REAL(SP), INTENT(in) :: sigma,minl,maxl
  REAL(SP), DIMENSION(nspec) :: tspec,tnspec,vel,func,gauss,psf,lnlam
  REAL(SP) :: ckms,cg,xmax,xmin,fwhm,psig,dlstep
  INTEGER :: i,j,il,ih,m=4,grange
  !convolve at fixed sigma_velocity
  !if set to 0, then the convolution is at fixed sigma_wavelength
  INTEGER, PARAMETER :: velocity=1

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
 
  IF (sigma.EQ.0) RETURN

  ckms = clight/1E13

  tspec = spec
   
  !compute smoothing the fast (and slightly less accurate) way
  IF (smoothspec_fast.EQ.1) THEN

     spec  = 0.0

     !convolve at fixed sigma_velocity
     IF (velocity.EQ.1) THEN
        
        DO i=1,nspec
           
           IF (lambda(i).LT.minl.OR.lambda(i).GT.maxl) THEN
              spec(i)=tspec(i)
              CYCLE
           ENDIF
           
           xmax = lambda(i)*(m*sigma/ckms+1)
           ih   = MIN(locate(lambda(1:nspec),xmax),nspec)
           il   = MAX(2*i-ih,1)
           
           IF (il.EQ.ih) THEN
              spec(i) = tspec(i)
           ELSE
              vel(il:ih)  = (lambda(i)/lambda(il:ih)-1)*ckms
              func(il:ih) =  1/SQRT(2*mypi)/sigma * &
                   EXP(-vel(il:ih)**2/2./sigma**2)
              !normalize the weights to integrate to unity
              func(il:ih) = func(il:ih) / TSUM(vel(il:ih),func(il:ih))
              spec(i) = TSUM(vel(il:ih),func(il:ih)*tspec(il:ih))
          ENDIF
           
        ENDDO

     !convolve at fixed sigma_wavelength
     ELSE
        
        DO i=1,nspec
           IF (lambda(i).LT.minl.OR.lambda(i).GT.maxl) THEN
              spec(i)=tspec(i)
              CYCLE
           ENDIF
           gauss = 1/SQRT(2*mypi)/sigma*EXP(-(lambda-lambda(i))**2/2/sigma**2)
           DO j=1,nspec
              spec(i) = spec(i) + gauss(j)*tspec(j)
           ENDDO
        ENDDO
        
     ENDIF

  !compute smoothing the correct (and somewhat slower) way
  !NB: the accuracy of this approach depends on the min/max wavelength
  !parameters through the density of the lnlam grid
  ELSE

     dlstep = (LOG(maxl)-LOG(minl))/nspec
     DO i=1,nspec
        lnlam(i) = i*dlstep+LOG(minl)
     ENDDO

     tspec = linterparr(LOG(lambda(1:nspec)),spec(1:nspec),lnlam)
  
     fwhm   = sigma*2.35482/ckms/dlstep
     psig   = fwhm/2.0/SQRT(-2.0*LOG(0.5)) ! equivalent sigma for kernel
     grange = FLOOR(m*psig)	           ! range for kernel (-range:range)

     DO i=1,2*grange+1
        psf(i) = 1.d0/SQRT(2.*mypi)/psig*EXP(-((i-grange-1)/psig)**2/2.)
     ENDDO
     psf(1:2*grange+1) = psf(1:2*grange+1) / SUM(psf(1:2*grange+1))
     
     DO i=grange+1,nspec-grange
        tnspec(i) = SUM( psf(1:2*grange+1)*tspec(i-grange:i+grange) )
     ENDDO
     
     !interpolate back to the main array
     il = locate(lambda,minl)
     ih = locate(lambda,maxl)
     spec(il:ih) = linterparr(EXP(lnlam),tnspec,lambda(il:ih))

  ENDIF


  RETURN

END SUBROUTINE SMOOTHSPEC
