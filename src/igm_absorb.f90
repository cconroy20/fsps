FUNCTION IGM_ABSORB(lam,spec,zz,factor)

  !routine to include IGM absorption via Madau (1995)
  !this routine includes a fudge factor (accessed by pset%igm_factor)
  !that allows the user to scale the IGM optical depth

  USE sps_vars; USE sps_utils, ONLY : locate
  IMPLICIT NONE

  REAL(SP), DIMENSION(nspec), INTENT(in) :: lam,spec
  REAL(SP), DIMENSION(nspec) :: igm_absorb,lobs,xc,tau
  REAL(SP), INTENT(in) :: zz,factor
  REAL(SP) :: z1,lylim,a_metal
  INTEGER, PARAMETER  :: nly=17
  INTEGER :: vv,i
  REAL(SP), DIMENSION(nly) :: lyw,lycoeff

  !--------------------------------------------------------------!

  lyw = (/1215.67, 1025.72, 972.537, 949.743, 937.803,&
          930.748, 926.226, 923.150, 920.963, 919.352,&
          918.129, 917.181, 916.429, 915.824, 915.329,&
          914.919, 914.576/)

  lycoeff = (/0.0036,0.0017,0.0011846,0.0009410,0.0007960,&
       0.0006967,0.0006236,0.0005665,0.0005200,0.0004817,&
       0.0004487,0.0004200,0.0003947,0.000372,0.000352,&
       0.0003334,0.00031644/)

  lylim   = 911.75
  a_metal = 0.0017

  z1   = 1+zz
  lobs = lam*z1
  xc   = lobs/lylim

  tau = 0.0

  !Ly series line blanketing
  DO i=1,nly
     IF (lam(1).GT.lyw(i)) CONTINUE
     vv = MIN(MAX(locate(lam,lyw(i)),1),nspec)
     tau(1:vv) = tau(1:vv) + lycoeff(i) * &
          (lobs(1:vv)/lyw(i))**3.46
     !add metal blanketing (this has ~no effect)
     IF (i.EQ.1) THEN 
        tau(1:vv) = tau(1:vv)+a_metal*(lobs(1:vv)/lyw(i))**1.68
     ENDIF
  ENDDO

  !LyC absorption
  IF (lam(1).LT.lylim) THEN
     vv = MIN(MAX(locate(lam,lylim),1),nspec)
     !approximation to Eqn 16 in Madau (1995); see his footnote 3
     tau(1:vv) = tau(1:vv) + &
          (0.25*xc(1:vv)**3*(z1**0.46-xc(1:vv)**0.46)) + &
          (9.4*xc(1:vv)**1.5*(z1**0.18-xc(1:vv)**0.18)) - &
          (0.7*xc(1:vv)**3*(xc(1:vv)**(-1.32)-z1**(-1.32))) - &
          (0.023*(z1**1.68-xc(1:vv)**1.68))
  ENDIF

  !the LyC fitting function seems to fall apart at really short 
  !wavelengths, so when tau starts to decrease, cap it at the max.
  vv = MAXLOC(tau,1)
  tau(1:vv) = tau(vv)

  !attenuate the input spectrum by the IGM
  !include a fudge factor to dial up/down the strength
  igm_absorb = spec * EXP(-tau*factor)


END FUNCTION IGM_ABSORB
