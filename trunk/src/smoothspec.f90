SUBROUTINE SMOOTHSPEC(lambda,spec,sigma,minl,maxl)

  !routine to compute velocity broadening of an input spectrum
  !this is only approximate b/c we are ignoring the variation
  !in dlambda with lambda within each integration
  !integration is truncated at +/-4*sigma

  USE sps_vars; USE sps_utils, ONLY : locate
  IMPLICIT NONE
  
  REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec
  REAL(SP), INTENT(in), DIMENSION(nspec) :: lambda
  REAL(SP), INTENT(in) :: sigma,minl,maxl
  REAL(SP), DIMENSION(nspec) :: tspec,vel,func,gauss
  REAL(SP) :: c,cg,xmax
  INTEGER :: i,j,il,ih,m=4
  !convolve at fixed sigma_velocity
  !if set to 0, then the convolution is at fixed sigma_wavelength
  INTEGER, PARAMETER :: velocity=1

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
 
  IF (sigma.EQ.0) RETURN

  c  = 2.998E5
  cg = 1/SQRT(2*mypi)/sigma

  tspec = spec
  spec  = 0.0
   
  !convolve at fixed sigma_velocity
  IF (velocity.EQ.1) THEN

     DO i=1,nspec
        
        IF (lambda(i).LT.minl.OR.lambda(i).GT.maxl) THEN
           spec(i)=tspec(i)
           CYCLE
        ENDIF

        xmax = lambda(i)*(m*sigma/c+1)
        ih = MIN(locate(lambda(1:nspec),xmax),nspec)
        il = MAX(2*i-ih,1)
        
        vel(il:ih)  = (lambda(i)/lambda(il:ih)-1)*c
        func(il:ih) = tspec(il:ih) * cg*EXP(-vel(il:ih)**2/2./sigma**2)

        IF (il.EQ.ih) THEN
           spec(i) = tspec(i)
        ELSE
           spec(i) = SUM( ABS(vel(il+1:ih)-vel(il:ih-1))*&
                (func(il+1:ih)+func(il:ih-1))/2. )
        ENDIF
       
     ENDDO

  !convolve at fixed sigma_wavelength
  ELSE

     DO i=1,nspec
        IF (lambda(i).LT.minl.OR.lambda(i).GT.maxl) THEN
           spec(i)=tspec(i)
           CYCLE
        ENDIF
        gauss = cg*EXP(-(lambda-lambda(i))**2/2/sigma**2)
        DO j=1,nspec
           spec(i) = spec(i) + gauss(j)*tspec(j)
        ENDDO
     ENDDO
     
  ENDIF

  RETURN

END SUBROUTINE SMOOTHSPEC
