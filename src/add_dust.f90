SUBROUTINE ADD_DUST(pset,csp1,csp2,specdust)

  USE nrtype; USE sps_vars; USE nr, ONLY : locate
  IMPLICIT NONE

  REAL(SP), DIMENSION(nspec), INTENT(in) :: csp1,csp2
  TYPE(PARAMS), INTENT(in) :: pset
  REAL(SP), DIMENSION(nspec), INTENT(out) :: specdust
  !number of points to use for dust clump integration
  !(has been tested for convergence)
  INTEGER, PARAMETER :: nclump=50
  INTEGER :: w63, i
  REAL(SP), DIMENSION(nspec)  :: diffuse_dust,tau2,katt,cspi
  REAL(SP), DIMENSION(nspec)  :: x,a,b,y,fa,fb,boost
  REAL(SP), DIMENSION(nclump) :: nden,wclump
  REAL(SP) :: clump_ave

  !-----------------------------------------------------!
  !-----------------Test input params-------------------!
  !-----------------------------------------------------!

  IF (dust_type.NE.0.AND.dust_type.NE.1.AND.dust_type.NE.2&
       .AND.dust_type.NE.3) THEN
     WRITE(*,*) 'ADD_DUST ERROR: unknown dust_type: ', dust_type
     STOP
  ENDIF

  IF (pset%uvb.LT.0) THEN
     WRITE(*,*) 'ADD_DUST ERROR: pset%uvb<0!'
     STOP
  ENDIF

  IF (pset%wgp1.LT.1.OR.pset%wgp2.LT.1) THEN
     WRITE(*,*) 'ADD_DUST ERROR: pset%wgp1 and/or pset%wgp2 '//&
          'out of bounds: ',pset%wgp1,pset%wgp2
     STOP
  ENDIF

  !-----------------------------------------------------!
  !---------------------Add dust------------------------!
  !-----------------------------------------------------!
 
  !set the attenuation curve for the cirrus dust

  !MW extinction law w/ UV bump
  IF (dust_type.EQ.1) THEN 

     tau2 = 0.0
     !use CCM89 extinction curve parameterization
     x = 1E4/spec_lambda
     y = x-1.82

     !IR
     tau2(mwdindex(2):mwdindex(1)) = &
          (0.574*x(mwdindex(2):mwdindex(1))**1.61) + &
          (-0.527*x(mwdindex(2):mwdindex(1))**1.61)/pset%mwr

     !optical+near-IR
     tau2(mwdindex(3):mwdindex(2)) = &
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
     boost = (x(mwdindex(3))/x)**6.*(tau2(mwdindex(3))-&
          (a(mwdindex(3))+b(mwdindex(3))/pset%mwr))
     tau2(mwdindex(4):mwdindex(3)) = &
          a(mwdindex(4):mwdindex(3)) + &
          b(mwdindex(4):mwdindex(3))/pset%mwr + &
          boost(mwdindex(4):mwdindex(3))

     !mid-UV
     fa = -0.04473*(x-5.9)**2-0.009779*(x-5.9)**3
     fb =  0.2130*(x-5.9)**2+0.1207*(x-5.9)**3
        a  = 1.752-0.316*x-0.104/((x-4.67)**2+0.341)*pset%uvb+fa
        b  = -3.09+1.825*x+1.206/((x-4.62)**2+0.263)*pset%uvb+fb
     tau2(mwdindex(5):mwdindex(4)) = &
          a(mwdindex(5):mwdindex(4)) + b(mwdindex(5):mwdindex(4))/pset%mwr

     !far-UV
     a = -1.073-0.628*(x-8.)+0.137*(x-8.)**2
     b = 13.67+4.257*(x-8.)-0.42*(x-8.)**2+0.374*(x-8.)**3
     tau2(mwdindex(6):mwdindex(5)) = &
          a(mwdindex(6):mwdindex(5)) + b(mwdindex(6):mwdindex(5))/pset%mwr

     a = -1.073-0.628*(12.-8.)+0.137*(12.-8.)**2
     b = 13.67+4.257*(12.-8.)-0.42*(12.-8.)**2+0.374*(12.-8.)**3
     tau2(:mwdindex(6)) = &
          a(:mwdindex(6)) + b(:mwdindex(6))/pset%mwr

     tau2 = pset%dust2*tau2

  ELSE IF (dust_type.EQ.0) THEN
     !power-law attenuation
     tau2 = (spec_lambda/5500.)**pset%dust_index * pset%dust2
  ELSE IF (dust_type.EQ.3) THEN
     !Witt & Gordon dust attenuation
     tau2 = wgdust(:,pset%wgp1,pset%wgp2,pset%wgp3)
  ENDIF

  !dust model with separate attenuation for young and old stars
  IF (dust_type.NE.2) THEN
     
     !compute dust extinction for standard uniform screen
     IF (pset%dust_clumps.LE.0.0.OR.dust_type.EQ.3) THEN
        
        diffuse_dust = EXP(-tau2)
        
     !compute dust extinction for a lognormal clump distribution   
     ELSE 
        
        !set up weights for integration
        wclump                  = 1.
        wclump(1:3)             = (/3./8,7./6,23./24/)
        wclump(nclump-2:nclump) = (/23./24,7./6,3./8/)

        !relation between the log mean and sigma of the Gaussian PDF
        !set in order to satisfy physical constraints
        clump_ave = 0.0709*pset%dust_clumps-1.68*pset%dust_clumps**2+&
             1.56*pset%dust_clumps**3-1.96*pset%dust_clumps**4+&
             0.886*pset%dust_clumps**5

        !this is logarithmic density
        DO i=1,nclump
           nden(i) = (i-1)*8.*pset%dust_clumps/nclump
        ENDDO
        nden = nden - SUM(nden)/nclump + clump_ave

        !integrate over clump density, weighted by PDF
        diffuse_dust = 0.0
        DO i=1,nclump
           diffuse_dust = diffuse_dust + &
                (nden(2)-nden(1))*wclump(i)/pset%dust_clumps/SQRT(2*mypi)*&
                EXP(-(nden(i)-clump_ave)**2/2/pset%dust_clumps**2) * &
                EXP(-10**nden(i)*tau2)
        ENDDO

     ENDIF

     !compute total spectrum, including both dust components
     cspi = csp1 * EXP(-pset%dust1*(spec_lambda/5500.)**(pset%dust1_index))*&
             (1-pset%frac_obrun) + csp1*pset%frac_obrun + csp2
     IF (pset%frac_nodust.GE.0.0) THEN
        specdust  = cspi*diffuse_dust*(1-pset%frac_nodust) + &
             cspi*pset%frac_nodust
     ELSE
        specdust  = cspi*diffuse_dust*(1+pset%frac_nodust) - &
             (csp1+csp2)*pset%frac_nodust
     ENDIF

  ELSE

     !Calzetti et al. 2000 attenuation is applied to the entire spectrum 
     w63 = locate(spec_lambda,6300.00)
     katt = 0.0
     katt(w63+1:) = 1.17*( -1.857+1.04*(1E4/spec_lambda(w63+1:)) ) + 1.78
     katt(1:w63)  = 1.17*(-2.156+1.509*(1E4/spec_lambda(1:w63))-&
          0.198*(1E4/spec_lambda(1:w63))**2 + &
          0.011*(1E4/spec_lambda(1:w63))**3) + 1.78
     !R=4.05
     katt = katt/0.44/4.05 * pset%dust2
     w63 = locate(katt,0.0)
     IF (w63.NE.nspec) THEN
        katt(w63+1:) = 0.0
     ENDIF
     specdust = (csp1+csp2)*EXP(-katt)*(1-pset%frac_nodust) + &
          (csp1+csp2)*pset%frac_nodust

  ENDIF

END SUBROUTINE ADD_DUST
