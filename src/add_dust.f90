SUBROUTINE ADD_DUST(pset,csp1,csp2,specdust,mdust)

  USE sps_vars
  USE sps_utils, ONLY : tsum, locate, attn_curve
  IMPLICIT NONE

  REAL(SP), DIMENSION(nspec), INTENT(in) :: csp1,csp2
  TYPE(PARAMS), INTENT(in) :: pset
  REAL(SP), DIMENSION(nspec), INTENT(out) :: specdust
  REAL(SP), INTENT(out) :: mdust
  INTEGER :: i,qlo,ulo,iself=0
  REAL(SP), DIMENSION(nspec)  :: diff_dust,tau2,cspi
  REAL(SP), DIMENSION(nspec)  :: nu,dumin,dumax
  REAL(SP), DIMENSION(nspec)  :: mduste,duste,oduste,sduste,tduste
  REAL(SP), DIMENSION(7) :: qpaharr
  REAL(SP) :: clump_ave,lboln,lbold,gamma,norm,dq,qpah,umin,du
  REAL(SP), DIMENSION(numin_dl07) :: uminarr

  !---------------------------------------------------------------!
  !----------------------Test input params------------------------!
  !---------------------------------------------------------------!

  IF (dust_type.NE.0.AND.dust_type.NE.1.AND.dust_type.NE.2&
       .AND.dust_type.NE.3.AND.dust_type.NE.4) THEN
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

  IF (pset%dust_clumps.GT.tiny_number) THEN
     WRITE(*,*) 'ADD_DUST ERROR: dust_clumps feature no longer supported'
     STOP
  ENDIF

  IF (pset%frac_obrun.LT.0.0.OR.pset%frac_obrun.GT.1.0.OR.&
       pset%frac_nodust.LT.0.0.OR.pset%frac_nodust.GT.1.0) THEN
     WRITE(*,*) 'ADD_DUST ERROR: frac_obrun and/or frac_nodust out of bounds'
     STOP
  ENDIF

  !---------------------------------------------------------------!
  !----------------------Add dust absorption----------------------!
  !---------------------------------------------------------------!

  !compute attenuation curve for diffuse dust
  diff_dust = EXP(-attn_curve(dust_type,pset))
 
  !combine old and young stars, attenuating the young
  !with a fixed power-law attn curve, and allowing a fraction
  !of the young stars to be dust-free ("OB runaways")
  cspi = csp1 * EXP(-pset%dust1*(spec_lambda/5500.)**(pset%dust1_index))*&
       (1-pset%frac_obrun) + csp1*pset%frac_obrun + csp2
  
  !allow a fraction of the diffuse dust spectrum to be dust-free
  specdust  = cspi*diff_dust*(1-pset%frac_nodust) + &
       cspi*pset%frac_nodust

  !---------------------------------------------------------------!
  !----------------------Add dust emission------------------------!
  !---------------------------------------------------------------!
  
  ! The dust spectrum is computed according to Draine & Li 2007
  ! we are computing the integral of the dust spectrum over P(U)dU
  ! by considering two components, a delta function at Umin and 
  ! a power-law distribution from Umin to Umax=1E6 and alpha=2
  ! the relative weights of the two components are given by gamma
  
  IF (add_dust_emission.EQ.1) THEN 

     !only add dust emission if there is absorption
     IF (pset%dust2.GT.tiny_number.OR.pset%dust1.GT.tiny_number) THEN

        iself=0
        !compute Lbol both before and after dust attenuation
        !this will determine the normalization of the dust emission
        nu    = clight/spec_lambda
        lbold = TSUM(nu,specdust)
        lboln = TSUM(nu,csp1+csp2)
     
        !set up qpah interpolation
        qpaharr = (/0.47,1.12,1.77,2.50,3.19,3.90,4.58/)
        !set limits to qpah: 0.0<qpah<10.0
        qpah = MAX(MIN(pset%duste_qpah,10.0),0.0)
        qlo  = MAX(MIN(locate(qpaharr,qpah),6),1)
        dq   = (qpah-qpaharr(qlo))/(qpaharr(qlo+1)-qpaharr(qlo))

        !set up Umin interpolation
        uminarr = (/0.1,0.15,0.2,0.3,0.4,0.5,0.7,0.8,1.0,1.2,1.5,2.0,&
             2.5,3.0,4.0,5.0,7.0,8.0,12.0,15.0,20.0,25.0/)
        !set limits on Umin: 0.1<Umin<25.0
        umin = MAX(MIN(pset%duste_umin,25.0),0.1)
        ulo  = MAX(MIN(locate(uminarr,umin),numin_dl07),1)
        du   = (umin-uminarr(ulo))/(uminarr(ulo+1)-uminarr(ulo))

        !set limits to gamma (gamma is a fraction)
        gamma = MAX(MIN(pset%duste_gamma,1.0),0.0)

        !bi-linear interpolation over qpah and Umin
        dumin = (1-dq)*(1-du)*dustem2_dl07(:,qlo,2*ulo-1) + &
             dq*(1-du)*dustem2_dl07(:,qlo+1,2*ulo-1) + &
             dq*du*dustem2_dl07(:,qlo+1,2*(ulo+1)-1) + &
             (1-dq)*du*dustem2_dl07(:,qlo,2*(ulo+1)-1)
        dumax = (1-dq)*(1-du)*dustem2_dl07(:,qlo,2*ulo) + &
             dq*(1-du)*dustem2_dl07(:,qlo+1,2*ulo) + &
             dq*du*dustem2_dl07(:,qlo+1,2*(ulo+1)) + &
             (1-dq)*du*dustem2_dl07(:,qlo,2*(ulo+1))

        !combine both parts of P(U)dU
        mduste = (1-gamma)*dumin + gamma*dumax
        mduste = MAX(mduste,tiny_number)

        !normalize the dust emission to the luminosity absorbed by 
        !the dust, i.e., demand that Lbol remains the same
        norm   = TSUM(nu,mduste)
        duste  = mduste/norm * (lboln-lbold)
        duste  = MAX(duste,tiny_number)

        !include dust self-absorption
        !we need to iterate because the dust
        !will re-emit and be re-absorbed, etc.
        tduste = 0.0
        DO WHILE (((lboln-lbold).GT.1E-2).OR.iself.EQ.0)

           oduste = duste
           duste  = duste * diff_dust
           tduste = tduste + duste
           
           lbold = TSUM(nu,duste)  !after  self-abs
           lboln = TSUM(nu,oduste) !before self-abs

           duste = MAX(mduste/norm*(lboln-lbold),tiny_number)
           
           iself=1

        ENDDO

        !this factor assumes Md/Mh=0.01 (appropriate for the 
        !MW3.1 models), and includes conversion factors from 
        !Jy -> Lsun and the hydrogen mass in solar units
        mdust = 3.21E-3 / 4/mypi * (lboln-lbold)/norm

        !add the dust emission to the stellar spectrum
        specdust = specdust + tduste
        
     ELSE
        mdust = tiny_number
     ENDIF

  ENDIF

END SUBROUTINE ADD_DUST
