SUBROUTINE INTSPEC(pset,nti,spec_ssp,csp,mass_ssp,lbol_ssp,&
     mass,lbol,specb,massb,lbolb,deltb,sfstart,tau,const,&
     sftrunc,tmax,mdust,tweight)

  !routine to perform integration of SSP over a SFH.
  !Dust absorption and re-radiation is also included here

  !note that integrals over the SFH are all formally wrong
  !basically were assuming int(SFR*spec) = <spec>*int(SFR)
  !we used to do something a priori more sensible but it turns
  !out that analytically calculating int(SFR) is more often than
  !not the better thing to do, esp when SFR is changing rapidly

  USE sps_vars
  USE sps_utils, ONLY : add_dust,intsfr,locate
  IMPLICIT NONE
  
  INTEGER :: i,imax,indsf,wtesc
  REAL(SP) :: t1,t2
  REAL(SP), DIMENSION(ntfull) :: isfr,time
  INTEGER, intent(in), optional :: tweight
  INTEGER,  INTENT(in)    :: nti
  REAL(SP), INTENT(in)    :: massb,lbolb,deltb,sfstart,tau,const,sftrunc,tmax
  REAL(SP), INTENT(inout) :: mass, lbol, mdust
  REAL(SP), INTENT(in), DIMENSION(ntfull) :: mass_ssp,lbol_ssp
  REAL(SP), INTENT(in), DIMENSION(nspec,ntfull) :: spec_ssp
  REAL(SP), INTENT(in), DIMENSION(nspec)    :: specb
  REAL(SP), INTENT(inout), DIMENSION(nspec) :: csp
  REAL(SP), DIMENSION(nspec)  :: csp1,csp2
  TYPE(PARAMS), INTENT(in)    :: pset

  !-----------------------------------------------------!
  !-----------------------------------------------------!

  time = 10**time_full

  csp  = 0.0
  csp1 = 0.0
  csp2 = 0.0
  mass = 0.0
  lbol = 0.0

  !only compute things if SF has "started"
  IF (time(nti).GT.sfstart) THEN

     !find where t~dust_tesc
     wtesc = MIN(MAX(locate(time_full,&
          LOG10(10**pset%dust_tesc+sfstart)),1),ntfull)
     IF (pset%dust1.EQ.0.0.AND.pset%dust2.EQ.0.0) wtesc = ntfull

     indsf = locate(time(nti)-time(1:nti),sfstart)
     IF (sfstart.LE.tiny_number) indsf = nti
     indsf = MIN(indsf,ntfull-1)
     imax  = MIN(MIN(wtesc,nti),indsf)

     !set up the integrated SFR array
     DO i=1,indsf

        t1   = time(nti)-sfstart-time(i)+time(1)
        IF (i.EQ.indsf) t2 = 0.0
        IF (i.NE.indsf) t2 = time(nti)-sfstart-time(i+1)+time(1)

        IF (PRESENT(tweight)) THEN
           isfr(i) = intsfr(pset%sfh,tau,const,sftrunc,sfstart,&
                pset%sf_slope,tmax,t1,t2,1)
        ELSE
           isfr(i) = intsfr(pset%sfh,tau,const,sftrunc,sfstart,&
                pset%sf_slope,tmax,t1,t2)
        ENDIF

     ENDDO

     IF (indsf.EQ.1) THEN
        csp1 = isfr(1)*spec_ssp(:,1)
     ELSE
        !t<tesc
        DO i=1,imax-1
           csp1 = csp1 + isfr(i)*0.5*(spec_ssp(:,i+1)+spec_ssp(:,i))
        ENDDO
        !t>tesc
        IF (nti.GT.wtesc) THEN
           DO i=wtesc,indsf-1
              csp2 = csp2 + isfr(i)*0.5*(spec_ssp(:,i+1)+spec_ssp(:,i))
           ENDDO
           csp2 = csp2 + isfr(indsf)*spec_ssp(:,indsf)
        ELSE
           csp1 = csp1 + isfr(indsf)*spec_ssp(:,indsf)
        ENDIF
     ENDIF

     !compute weighted mass and lbol
     IF (indsf.EQ.1) THEN
        mass = isfr(1)*mass_ssp(1)
        lbol = LOG10( isfr(1)*10**lbol_ssp(1)+tiny_number )
     ELSE
        mass = SUM(isfr(1:indsf-1)/2.*&
             (mass_ssp(1:indsf-1)+mass_ssp(2:indsf)))
        lbol = SUM(isfr(1:indsf-1)/2.*&
             (10**lbol_ssp(1:indsf-1)+10**lbol_ssp(2:indsf)))
        mass = mass + isfr(indsf)*mass_ssp(indsf)
        lbol = lbol + isfr(indsf)*10**lbol_ssp(indsf) 
        lbol = LOG10(lbol+tiny_number)
     ENDIF
  
     !add in an instantaneous burst
     IF (pset%fburst.GT.tiny_number.AND.&
          (pset%sfh.EQ.1.OR.pset%sfh.EQ.4)) THEN
        IF (deltb.LT.10**pset%dust_tesc) THEN
           csp1 = csp1*(1-pset%fburst) + specb*pset%fburst
           csp2 = csp2*(1-pset%fburst)
        ELSE
           csp1 = csp1*(1-pset%fburst)
           csp2 = csp2*(1-pset%fburst) + specb*pset%fburst
        ENDIF
        mass = mass*(1-pset%fburst) + massb*pset%fburst
        lbol = LOG10( 10**lbol*(1-pset%fburst) + &
             10**lbolb*pset%fburst )
     ENDIF

     !add dust, combine young and old csp
     CALL ADD_DUST(pset,csp1,csp2,csp,mdust)

  ENDIF

END SUBROUTINE INTSPEC
