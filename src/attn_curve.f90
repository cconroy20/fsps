FUNCTION ATTN_CURVE(lambda,dtype,pset)

  ! Routine to generate and return the attenuation curve for a chosen 
  ! dust type.  The V-band optical depth (dust2) is passed via the 
  ! pset variable.  Option 3 (Witt & Gordon models) do not use the
  ! dust2 option as those grids fully predict the actual attn curves.

  USE sps_vars
  USE sps_utils, ONLY : locate
  IMPLICIT NONE

  INTEGER, INTENT(in) :: dtype
  TYPE(PARAMS), INTENT(in) :: pset
  INTEGER  :: w63,w1,w2
  REAL(SP) :: eb,zero=0.0,dd63=6300.00,lamv=5500.0,dlam=350.0,lamuvb=2175.0
  REAL(SP), INTENT(in), DIMENSION(nspec) :: lambda
  REAL(SP), DIMENSION(nspec) :: x,a,b,y,fa,fb,hack,cal00,reddy
  REAL(SP), DIMENSION(nspec) :: attn_curve,drude,tmp
 
  !---------------------------------------------------------------------!

  attn_curve = 0.0

  IF (dtype.LT.0.OR.dtype.GT.6) THEN
     WRITE(*,*) 'ATTN_CURVE ERROR: dust_type out of range:',dtype
     STOP
  ENDIF
  
  !-------------------------power-law attenuation-----------------------!

  !power-law attenuation
  IF (dtype.EQ.0) THEN 

     attn_curve = (lambda/lamv)**pset%dust_index * pset%dust2

  !-----------------CCM89 MW curve w/ variable UV bump------------------!

  !MW extinction law w/ UV bump
  ELSE IF (dtype.EQ.1) THEN 

     tmp = 0.0

     !use CCM89 extinction curve parameterization
     x = 1E4/lambda
     y = x-1.82

     !IR
     tmp(mwdindex(2):mwdindex(1)) = &
          (0.574*x(mwdindex(2):mwdindex(1))**1.61) + &
          (-0.527*x(mwdindex(2):mwdindex(1))**1.61)/pset%mwr

     !optical+near-IR
     tmp(mwdindex(3):mwdindex(2)) = &
          (1+0.17699*y(mwdindex(3):mwdindex(2))-&
          0.50447*y(mwdindex(3):mwdindex(2))**2-&
          0.02427*y(mwdindex(3):mwdindex(2))**3+&
          0.72085*y(mwdindex(3):mwdindex(2))**4+&
          0.01979*y(mwdindex(3):mwdindex(2))**5-&
          0.77530*y(mwdindex(3):mwdindex(2))**6+&
          0.32999*y(mwdindex(3):mwdindex(2))**7)+&
          (1.41338*y(mwdindex(3):mwdindex(2))+&
          2.28305*y(mwdindex(3):mwdindex(2))**2+&
          1.07233*y(mwdindex(3):mwdindex(2))**3-&
          5.38434*y(mwdindex(3):mwdindex(2))**4-&
          0.62251*y(mwdindex(3):mwdindex(2))**5+&
          5.3026*y(mwdindex(3):mwdindex(2))**6-&
          2.09002*y(mwdindex(3):mwdindex(2))**7)/pset%mwr

     !near-UV
     a = 1.752-0.316*x-0.104/((x-4.67)**2+0.341)*pset%uvb
     b = -3.09+1.825*x+1.206/((x-4.62)**2+0.263)*pset%uvb
     !this hack parameter is not in the original CCM89 
     !parameterization.  It is a hack designed to result in 
     !a smooth profile in the presence of a variable UVB strength
     hack = (x(mwdindex(3))/x)**6.*(tmp(mwdindex(3))-&
          (a(mwdindex(3))+b(mwdindex(3))/pset%mwr))
     tmp(mwdindex(4):mwdindex(3)) = &
          a(mwdindex(4):mwdindex(3)) + &
          b(mwdindex(4):mwdindex(3))/pset%mwr + &
          hack(mwdindex(4):mwdindex(3))

     !mid-UV
     fa = -0.04473*(x-5.9)**2-0.009779*(x-5.9)**3
     fb =  0.2130*(x-5.9)**2+0.1207*(x-5.9)**3
     a  = 1.752-0.316*x-0.104/((x-4.67)**2+0.341)*pset%uvb+fa
     b  = -3.09+1.825*x+1.206/((x-4.62)**2+0.263)*pset%uvb+fb
     tmp(mwdindex(5):mwdindex(4)) = &
          a(mwdindex(5):mwdindex(4)) + b(mwdindex(5):mwdindex(4))/pset%mwr

     !far-UV
     a = -1.073-0.628*(x-8.)+0.137*(x-8.)**2-0.070*(x-8.)**3
     b = 13.67+4.257*(x-8.)-0.42*(x-8.)**2+0.374*(x-8.)**3
     tmp(mwdindex(6):mwdindex(5)) = &
          a(mwdindex(6):mwdindex(5)) + b(mwdindex(6):mwdindex(5))/pset%mwr

     !set to a constant at lambda < 12 um^-1
     a = -1.073-0.628*(12.-8.)+0.137*(12.-8.)**2-0.070*(12.-8.)**3
     b = 13.67+4.257*(12.-8.)-0.42*(12.-8.)**2+0.374*(12.-8.)**3
     tmp(:mwdindex(6)) = &
          a(:mwdindex(6)) + b(:mwdindex(6))/pset%mwr

     attn_curve = pset%dust2*tmp


  !------------------Calzetti et al. 2000 attenuation-------------------!

  ELSE IF (dtype.EQ.2) THEN 

     w63   = locate(lambda,dd63)
     cal00 = 0.0
     cal00(w63+1:) = 1.17*( -1.857+1.04*(1E4/lambda(w63+1:)) ) + 1.78
     cal00(1:w63)  = 1.17*(-2.156+1.509*(1E4/lambda(1:w63))-&
          0.198*(1E4/lambda(1:w63))**2 + &
          0.011*(1E4/lambda(1:w63))**3) + 1.78
     cal00 = cal00/0.44/4.05  !R=4.05
     w63   = locate(cal00,zero)
     IF (w63.NE.nspec) THEN
        cal00(w63+1:) = 0.0
     ENDIF
 
     attn_curve = cal00 * pset%dust2

  !------------------Witt & Gordon 2000 attenuation--------------------!

  ELSE IF (dtype.EQ.3) THEN
  
     attn_curve = wgdust(:,pset%wgp1,pset%wgp2,pset%wgp3)

  !------------------Kriek & Conroy 2013 attenuation-------------------!

  ELSE IF (dtype.EQ.4) THEN

     !Calzetti curve
     w63   = locate(lambda,dd63)
     cal00 = 0.0
     cal00(w63+1:) = 1.17*( -1.857+1.04*(1E4/lambda(w63+1:)) ) + 1.78
     cal00(1:w63)  = 1.17*(-2.156+1.509*(1E4/lambda(1:w63))-&
          0.198*(1E4/lambda(1:w63))**2 + &
          0.011*(1E4/lambda(1:w63))**3) + 1.78
     !R=4.05 NB: I'm not sure I have this normalization correct...
     cal00 = cal00/0.44/4.05 
     w63   = locate(cal00,zero)
     IF (w63.NE.nspec) THEN
        cal00(w63+1:) = 0.0
     ENDIF

     eb = 0.85 - 1.9 * pset%dust_index  !KC13 Eqn 3

     !Drude profile for 2175A bump
     drude = eb*(lambda*dlam)**2 / &
          ( (lambda**2-lamuvb**2)**2 + (lambda*dlam)**2 )

     attn_curve = pset%dust2*(cal00+drude/4.05)*&
          (lambda/lamv)**pset%dust_index


  !-----------------Gordon et al. (2003) SMC exctincion----------------!

  ELSE IF (dtype.EQ.5) THEN

     attn_curve = pset%dust2*g03smcextn
     
  !------------------Reddy et al. (2015) attenuation-------------------!

  ELSE IF (dtype.EQ.6) THEN

     reddy = 0.0

     ! see Eqn. 8 in Reddy et al. (2015)
     
     w1   = locate(lambda,1500.d0)
     w2   = locate(lambda,6000.d0)
     reddy(w1:w2) = -5.726 + 4.004/(lambda(w1:w2)/1E4) - 0.525/(lambda(w1:w2)/1E4)**2 + &
          0.029/(lambda(w1:w2)/1E4)**3 + 2.505

     reddy(1:w1) = reddy(w1) ! constant extrapolation blueward
     
     w1   = locate(lambda,6000.d0)
     w2   = locate(lambda,28500.d0)
     ! note the last term is not in Reddy et al. but was included to make the
     ! two functions continuous at 0.6um
     reddy(w1:w2) = -2.672 - 0.010/(lambda(w1:w2)/1E4) + 1.532/(lambda(w1:w2)/1E4)**2 + &
          -0.412/(lambda(w1:w2)/1E4)**3 + 2.505 - 0.036221981

     ! convert k_lam to A_lam/A_V assuming Rv=2.505
     reddy = reddy/2.505

     attn_curve = pset%dust2*reddy
     
  ENDIF


END FUNCTION ATTN_CURVE
