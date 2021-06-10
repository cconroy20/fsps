SUBROUTINE ADD_DUST(pset,csp1,csp2,specdust,mdust,ncsp1,ncsp2,nebdust)

  USE sps_vars
  USE sps_utils, ONLY : tsum, locate, attn_curve, linterparr
  IMPLICIT NONE

  REAL(SP), DIMENSION(nspec), INTENT(in) :: csp1,csp2
  TYPE(PARAMS), INTENT(in) :: pset
  REAL(SP), DIMENSION(nspec), INTENT(out) :: specdust
  REAL(SP), INTENT(out) :: mdust
  REAL(SP), DIMENSION(nemline), INTENT(in) :: ncsp1,ncsp2
  REAL(SP), DIMENSION(nemline), INTENT(out) :: nebdust
  INTEGER :: i,qlo,ulo,iself=0
  REAL(SP), DIMENSION(nspec)  :: diff_dust,tau2,cspi
  REAL(SP), DIMENSION(nemline)  :: diff_dust_neb,ncspi
  REAL(SP), DIMENSION(nspec)  :: nu,dumin,dumax
  REAL(SP), DIMENSION(nspec)  :: mduste,duste,oduste,sduste,tduste
  REAL(SP) :: clump_ave,lboln,lbold,labs,gamma,norm,dq,du

  !---------------------------------------------------------------!
  !----------------------Test input params------------------------!
  !---------------------------------------------------------------!

  IF (dust_type.LT.0.OR.dust_type.GT.6) THEN
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
  diff_dust = EXP(-attn_curve(spec_lambda,dust_type,pset))
 
  !combine old and young stars, attenuating the young
  !with a fixed power-law attn curve, and allowing a fraction
  !of the young stars to be dust-free ("OB runaways")
  cspi = csp1 * EXP(-pset%dust1*(spec_lambda/5500.)**(pset%dust1_index))*&
       (1-pset%frac_obrun) + csp1*pset%frac_obrun + csp2
  
  !allow a fraction of the diffuse dust spectrum to be dust-free
  specdust  = cspi*diff_dust*(1-pset%frac_nodust) + &
       cspi*pset%frac_nodust

  !as above, for nebular line luminosities
  diff_dust_neb = linterparr(spec_lambda,diff_dust,nebem_line_pos)
  ncspi = ncsp1 * EXP(-pset%dust1*(nebem_line_pos/5500.)**(pset%dust1_index))*&
       (1-pset%frac_obrun) + ncsp1*pset%frac_obrun + ncsp2
  nebdust  = ncspi*diff_dust_neb*(1-pset%frac_nodust) + &
       ncspi*pset%frac_nodust

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
        IF (nebemlineinspec.EQ.0) THEN
           lboln = lboln + SUM(ncsp1) + SUM(ncsp2) - SUM(nebdust)
        ENDIF

        !set up qpah interpolation
        qlo = MAX(MIN(locate(qpaharr,pset%duste_qpah),nqpah_dustem-1),1)
        dq  = (pset%duste_qpah-qpaharr(qlo))/(qpaharr(qlo+1)-qpaharr(qlo))
        dq  = MIN(MAX(dq,0.0),1.0)  !no extrapolation
        
        !set up Umin interpolation
        !set limits on Umin: 0.1<Umin<25.0
        ulo = MAX(MIN(locate(uminarr,pset%duste_umin),numin_dustem),1)
        du  = (pset%duste_umin-uminarr(ulo))/(uminarr(ulo+1)-uminarr(ulo))
        du  = MIN(MAX(du,0.0),1.0)  !no extrapolation

        !set limits to gamma (gamma is a fraction)
        gamma = MAX(MIN(pset%duste_gamma,1.0),0.0)

        !bi-linear interpolation over qpah and Umin
        dumin = (1-dq)*(1-du)*dustem2_dustem(:,qlo,2*ulo-1) + &
             dq*(1-du)*dustem2_dustem(:,qlo+1,2*ulo-1) + &
             dq*du*dustem2_dustem(:,qlo+1,2*(ulo+1)-1) + &
             (1-dq)*du*dustem2_dustem(:,qlo,2*(ulo+1)-1)
        dumax = (1-dq)*(1-du)*dustem2_dustem(:,qlo,2*ulo) + &
             dq*(1-du)*dustem2_dustem(:,qlo+1,2*ulo) + &
             dq*du*dustem2_dustem(:,qlo+1,2*(ulo+1)) + &
             (1-dq)*du*dustem2_dustem(:,qlo,2*(ulo+1))

        !combine both parts of P(U)dU
        mduste = (1-gamma)*dumin + gamma*dumax
        mduste = MAX(mduste,tiny_number)

        !normalize the dust emission to the luminosity absorbed by 
        !the dust, i.e., demand that Lbol remains the same
        labs = (lboln-lbold)
        norm   = TSUM(nu,mduste)
        duste  = mduste/norm * labs
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
        mdust = 3.21E-3 / 4/mypi * labs/norm

        !add the dust emission to the stellar spectrum
        specdust = specdust + tduste
        
     ELSE
        mdust = tiny_number
     ENDIF

  ENDIF

END SUBROUTINE ADD_DUST
