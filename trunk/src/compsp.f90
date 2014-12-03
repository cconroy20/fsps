!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Program to compute magnitudes and spectra for a composite !    
!  stellar population.  Returns a file with the following    !
!  info: age, mass, Lbol, mags in various filters.  SFR      !
!  units are Msun/yr.                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION INTSFR(sfh,tau,const,sftrunc,sfstart,sftheta,t1,t2,tweight)

  !simple routine to integrate the SFH from t1 to t2

  USE sps_vars; USE sps_utils, ONLY : tsum, locate
  IMPLICIT NONE

  INTEGER, INTENT(in)  :: sfh
  REAL(SP), INTENT(in) :: t1,t2,tau,const,sftrunc,sfstart,sftheta
  INTEGER, INTENT(in), OPTIONAL :: tweight
  INTEGER  :: lo,hi
  REAL(SP) :: s1,s2,dt,intsfr,tt1

  !-----------------------------------------------------!
  !-----------------------------------------------------!

  IF (t1.GT.sftrunc.AND.t2.LT.sftrunc) THEN
     tt1 = sftrunc
  ELSE
     tt1 = t1
  ENDIF

  IF (sfh.EQ.1) THEN

     IF (PRESENT(tweight)) THEN
        !in this case the integral is int(t*sfr)
        intsfr = tau*( (t2/tau/1E9+1)*EXP(-t2/tau/1E9) - &
             (tt1/tau/1E9+1)*EXP(-tt1/tau/1E9) ) / &
             (1-EXP(-(sftrunc-sfstart)/1E9/tau)) 
        intsfr = intsfr*(1-const) + const*(tt1**2-t2**2)/2E9/(sftrunc-sfstart)
     ELSE
        intsfr = (EXP(-t2/tau/1E9) - EXP(-tt1/tau/1E9))/&
             (1-EXP(-(sftrunc-sfstart)/1E9/tau)) 
        intsfr = intsfr*(1-const) + const*(tt1-t2)/(sftrunc-sfstart)
     ENDIF

     IF (t2.GT.sftrunc) THEN
        intsfr = 0.0
     ENDIF

  ELSE IF (sfh.EQ.4) THEN

     IF (PRESENT(tweight)) THEN
        !in this case the integral is int(t*sfr)
        intsfr = tau*(EXP(-t2/tau/1E9)*(2+t2/tau/1E9*(t2/tau/1E9+2)) - &
                  EXP(-tt1/tau/1E9)*(2+tt1/tau/1E9*(tt1/tau/1E9+2))) / &
                  (1-EXP(-(sftrunc-sfstart)/1E9/tau)*((sftrunc-sfstart)/1E9/tau+1)) 
        intsfr = intsfr*(1-const) + const*(tt1**2-t2**2)/2E9/(sftrunc-sfstart)
     ELSE
        intsfr = (EXP(-t2/tau/1E9)*(1+t2/tau/1E9) - &
             EXP(-tt1/tau/1E9)*(1+tt1/tau/1E9)) / &
             (1-EXP(-(sftrunc-sfstart)/1E9/tau)*((sftrunc-sfstart)/1E9/tau+1))  
        intsfr = intsfr*(1-const) + const*(tt1-t2)/(sftrunc-sfstart)
     ENDIF

     IF (t2.GT.sftrunc) THEN
        intsfr = 0.0
     ENDIF

  ELSE IF (sfh.EQ.5) THEN
     
     IF ((t1+sfstart).LT.sftrunc) THEN
        intsfr = (EXP(-t2/tau/1E9)*(1+t2/tau/1E9) - &
             EXP(-t1/tau/1E9)*(1+t1/tau/1E9)) * (tau*1E9)**2 
     ELSE
       intsfr =  TAN(sftheta)*&
             (0.5*(t1**2-t2**2)-(t1-t2)*(sftrunc-sfstart))
     ENDIF
     intsfr = MAX(intsfr,0.0) / 1E10

  ELSE IF (sfh.EQ.2.OR.sfh.EQ.3) THEN

     !handle the edges separately
     hi = MIN(MAX(locate(sfh_tab(1,1:ntabsfh),t1),1),ntabsfh-1)
     dt = (t1-sfh_tab(1,hi))/(sfh_tab(1,hi+1)-sfh_tab(1,hi))
     s1 = MAX((1-dt)*sfh_tab(2,hi) + dt*sfh_tab(2,hi+1),0.0)

     lo = MIN(MAX(locate(sfh_tab(1,1:ntabsfh),t2),1),ntabsfh-1)
     dt = (t2-sfh_tab(1,lo))/(sfh_tab(1,lo+1)-sfh_tab(1,lo))
     s2 = MAX((1-dt)*sfh_tab(2,lo) + dt*sfh_tab(2,lo+1),0.0)

     IF ((hi-lo).LT.2) THEN
        intsfr = (t1-t2)*0.5*(s1+s2)
     ELSE
        lo = lo+1
        intsfr = TSUM(sfh_tab(1,lo:hi),sfh_tab(2,lo:hi))
        intsfr = intsfr + (t1-sfh_tab(1,hi))*0.5*(s1+sfh_tab(2,hi))
        intsfr = intsfr + (sfh_tab(1,lo)-t2)*0.5*(s2+sfh_tab(2,lo))
     ENDIF

  ENDIF
  
  
END FUNCTION INTSFR

!------------------------------------------------------------!
!------------------------------------------------------------!

SUBROUTINE INTSPEC(pset,nti,spec_ssp,csp,mass_ssp,lbol_ssp,&
     mass,lbol,specb,massb,lbolb,deltb,sfstart,tau,const,&
     sftrunc,mdust,tweight)

  !routine to perform integration of SSP over a SFH,
  !including attenuation by dust.

  USE sps_vars
  USE sps_utils, ONLY : add_dust,intsfr,locate
  IMPLICIT NONE
  
  INTEGER :: i,imax,indsf,wtesc
  REAL(SP) :: t1,t2
  REAL(SP), DIMENSION(ntfull) :: isfr,time
  INTEGER, intent(in), optional :: tweight
  INTEGER,  INTENT(in)    :: nti
  REAL(SP), INTENT(in)    :: massb,lbolb,deltb,sfstart,tau,const,sftrunc
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
                pset%sf_theta,t1,t2,1)
        ELSE
           isfr(i) = intsfr(pset%sfh,tau,const,sftrunc,sfstart,&
                pset%sf_theta,t1,t2)
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

!-------------------------------------------------------------------!
!-------------------------Main Routine------------------------------!
!-------------------------------------------------------------------!

SUBROUTINE COMPSP(write_compsp,nzin,outfile,mass_ssp,&
     lbol_ssp,tspec_ssp,pset,ocompsp)

  !sfh=1: tau model
  !sfh=2: tabulated SFH (from file)
  !sfh=3: tabulated SFH (stored in sfhtab arr)
  !sfh=4: delayed tau model
  !sfh=5: custom SFH (see Simha et al. 2013)

  !for sfh=1 or 4:
  !If tage >0  -> run only one integration to t=tage
  !If tage<=0  -> produce outputs from tmin<t<tmax

  USE sps_vars
  USE sps_utils, ONLY : getmags,add_dust,linterp,intspec,&
       smoothspec,locate,getindx,write_isochrone,vactoair,igm_absorb
  IMPLICIT NONE
 
  !write_compsp = 1->write mags, 2->write spectra
  !               3->write mags+spec, 4->write indices
  !               5->write CMDs
  INTEGER, INTENT(in) :: write_compsp,nzin
  REAL(SP), INTENT(in), DIMENSION(ntfull,nzin) :: lbol_ssp,mass_ssp
  REAL(SP), INTENT(in), DIMENSION(nspec,ntfull,nzin) :: tspec_ssp
  REAL(SP), DIMENSION(nspec,ntfull,nzin) :: spec_ssp
  CHARACTER(100), INTENT(in) :: outfile

  INTEGER  :: i,j,n,k,stat,klo,jlo,ilo,imin,imax,indsf,indsft
  REAL(SP) :: tau,const,maxtime,psfr,sfstart,zhist,sftrunc
  REAL(SP) :: mass_csp,lbol_csp,dtb,dt,dz,zred=0.,t1,t2,mdust
  REAL(SP) :: mass_burst=0.0,lbol_burst=0.0,delt_burst=0.0,zero=0.0
  REAL(SP), DIMENSION(nbands)  :: mags
  REAL(SP), DIMENSION(nindx)   :: indx
  REAL(SP), DIMENSION(nspec,ntfull) :: ispec
  REAL(SP), DIMENSION(nspec)   :: spec_burst=0.0,csp1,csp2,spec1,spec_csp
  REAL(SP), DIMENSION(ntfull)  :: imass,ilbol,powtime
  REAL(SP), DIMENSION(ntfull)  :: sfr,tsfr
  REAL(SP), DIMENSION(ntabmax) :: tlb
  TYPE(PARAMS), INTENT(in) :: pset
  TYPE(COMPSPOUT), INTENT(inout), DIMENSION(ntfull) :: ocompsp

  !-------------------------------------------------------------!
  !-------------------------------------------------------------!

  !dump the input SSPs into a temporary array so that we can
  !edit the SSPs (this is necessary in order to add nebular em)
  spec_ssp = tspec_ssp

  IF (check_sps_setup.EQ.0) THEN
     WRITE(*,*) 'COMPSP ERROR: '//&
          'SPS_SETUP must be run before calling COMPSP. '
     STOP
  ENDIF

  !-------------------------------------------------------------!
  !-----------------Write the CMDs and exit---------------------!
  !-------------------------------------------------------------!

  IF (write_compsp.EQ.5) THEN
     CALL WRITE_ISOCHRONE(outfile,pset)
     RETURN
  ENDIF

  !-------------------------------------------------------------!
  !---------------------Add Nebular Emission--------------------!
  !-------------------------------------------------------------!

  IF (add_neb_emission.EQ.1) THEN
     IF (nzin.GT.1) THEN
        WRITE(*,*) 'COMPSP ERROR: cannot handle both nebular '//&
             'emission and mult-metallicity SSPs in compsp'
        STOP
     ENDIF
     CALL ADD_NEBULAR(pset,tspec_ssp(:,:,1),spec_ssp(:,:,1))
  ENDIF


  !-------------------------------------------------------------!
  !---------------------Basic CSP Setup-------------------------!
  !-------------------------------------------------------------!

  dtb        = 0.0
  spec_burst = 0.0
  mass_burst = 0.0
  lbol_burst = 0.0
  sfstart    = 0.0
 
  powtime = 10**time_full
  maxtime = powtime(ntfull)
  imin    = 1
  imax    = ntfull

  !SFH-specific setup 
  IF (pset%sfh.EQ.2) THEN

     IF (TRIM(pset%sfh_filename).EQ.'') THEN
        OPEN(3,FILE=TRIM(SPS_HOME)//'/data/sfh.dat',ACTION='READ',STATUS='OLD')
     ELSE
        OPEN(3,FILE=TRIM(SPS_HOME)//'/data/'//TRIM(pset%sfh_filename),&
             ACTION='READ',STATUS='OLD')
     ENDIF
     DO n=1,ntabmax
        IF (nzin.EQ.nz) THEN
           READ(3,*,IOSTAT=stat) sfh_tab(1,n),sfh_tab(2,n),sfh_tab(3,n)
        ELSE
           READ(3,*,IOSTAT=stat) sfh_tab(1,n),sfh_tab(2,n)
           sfh_tab(3,n)=0.0
        ENDIF
        IF (stat.NE.0) GOTO 29
     ENDDO
     WRITE(*,*) 'COMPSP ERROR: didnt finish reading in the sfh file,'
     WRITE(*,*) '     increase ntabmax variable in sps_vars.f90 file'
     STOP
29   CONTINUE
     CLOSE(3)
     ntabsfh = n-1
     sfh_tab(1,1:ntabsfh) = sfh_tab(1,1:ntabsfh)*1E9 !convert to yrs
     !special switch to compute only the last time output
     !in the tabulated file
     IF (pset%tage.EQ.-99.) imin=imax

  ELSE IF (pset%sfh.EQ.3) THEN 

     !sfh_tab array is supposed to already be filled in, check that it is
     IF (ntabsfh.EQ.0) THEN 
        WRITE(*,*) 'COMPSP ERROR: sfh=3 but sfh_tab array not initialized!'
        STOP
     ENDIF

  ELSE

     !set up maxtime variable
     !if tage > 0 then only output one age=tage,
     !otherwise output ages from 0<t<maxtime
     IF (pset%tage.GT.tiny_number) THEN
        maxtime = pset%tage*1E9
        imin    = MIN(MAX(locate(powtime,maxtime),1),ntfull-1)
        imax    = imin+1
     ENDIF

     !find sf_start in the time grid
     IF (pset%sf_start.GT.tiny_number) THEN
        sfstart = pset%sf_start*1E9 !convert to yrs
        indsf = MIN(MAX(locate(powtime,sfstart),1),ntfull-1)
     ELSE
        indsf   = 1
        sfstart = 0.0
     ENDIF

     !find sf_trunc in the time grid
     IF (pset%sf_trunc.GT.tiny_number) THEN
        sftrunc = pset%sf_trunc*1E9 !convert to yrs
        indsft = MIN(MAX(locate(powtime,sftrunc),1),ntfull)
     ELSE
        indsft  = ntfull
        sftrunc = maxtime
     ENDIF

     !set limits on the parameters tau and const
     IF (pset%sfh.EQ.1.OR.pset%sfh.EQ.4.OR.pset%sfh.EQ.5) THEN
        tau   = MIN(MAX(pset%tau,0.1),100.) !tau in Gyr
        const = MIN(MAX(pset%const,0.0),1.0)
     ENDIF

  ENDIF

  !make sure various variables are set correctly
  CALL COMPSP_WARNING(maxtime,pset,nzin,write_compsp)

  !setup output files
  IF (write_compsp.GT.0) &
       CALL COMPSP_SETUP_OUTPUT(write_compsp,pset,outfile,imin,imax)

  !tsfr is a function of the age of Universe
  !these are only used for writing the SFH to file
  IF (pset%sfh.EQ.2.OR.pset%sfh.EQ.3) THEN

     !linearly interpolate the tabulated SFH to the internal time grid
     DO j=1,ntfull
        IF (powtime(j).GT.sfh_tab(1,ntabsfh).OR.&
             powtime(j).LT.sfh_tab(1,1)) THEN
           tsfr(j)   = 0.0
        ELSE
           jlo    = MAX(MIN(locate(LOG10(sfh_tab(1,1:ntabsfh)),&
                time_full(j)),ntabsfh-1),1)
           dt = (powtime(j)-(sfh_tab(1,jlo))) / &
                (sfh_tab(1,jlo+1)-sfh_tab(1,jlo))
           tsfr(j)   = (1-dt)*sfh_tab(2,jlo)+dt*sfh_tab(2,jlo+1)
        ENDIF
     ENDDO

  ELSE IF (pset%sfh.EQ.1.OR.pset%sfh.EQ.4.OR.pset%sfh.EQ.5) THEN

     tsfr  = 0.0
     IF (pset%sfh.EQ.1) &
          tsfr(indsf:indsft)  = EXP(-(powtime(indsf:indsft)-sfstart)/tau/1E9 )/&
          tau/1E9 / (1-EXP(-(sftrunc-sfstart)/1E9/tau))
     IF (pset%sfh.EQ.4) &
          tsfr(indsf:indsft)  = ((powtime(indsf:indsft)-sfstart)/tau/1E9)*&
          EXP(-(powtime(indsf:indsft)-sfstart)/tau/1E9 )/tau/1E9 / &
          (1-EXP(-(sftrunc-sfstart)/1E9/tau)*((sftrunc-sfstart)/1E9/tau+1))
     tsfr(indsf:indsft) = tsfr(indsf:indsft)*(1-const) + const/(sftrunc-sfstart)

     IF (pset%sfh.EQ.5) THEN
        tsfr(indsf:indsft)  = (powtime(indsf:indsft)-sfstart)*&
             EXP(-(powtime(indsf:indsft)-sfstart)/tau/1E9 )
        IF (indsft.LT.ntfull) tsfr(indsft+1:) = tsfr(indsft) + &
             TAN(pset%sf_theta)*(powtime(indsft+1:)-sftrunc)
        tsfr = MAX(tsfr,0.0) ! set SFR=0.0 if SFR<0
     ENDIF

  ELSE IF (pset%sfh.EQ.0) THEN
     tsfr = 0.0
  ENDIF


  !-------------------------------------------------------------!
  !-------------Generate composite spectra and mags-------------!
  !-------------------------------------------------------------!

   !calculate mags at each time step
   DO i=imin,imax

      !Set up tabulated SFH
      IF (pset%sfh.EQ.2.OR.pset%sfh.EQ.3) THEN
         
         IF (nzin.EQ.nz) THEN

            ilbol = 0.0
            imass = 0.0
            ispec = 0.0
            tlb   = 0.0
            ilo = MAX(MIN(locate(LOG10(sfh_tab(1,1:ntabsfh)),&
                 time_full(i)),ntabsfh-1),1)
            tlb(1:ilo) = LOG10(sfh_tab(1,ilo) - sfh_tab(1,1:ilo) + powtime(1))

            DO j=1,i
               !interpolation in time (in logarithmic units)
               jlo = MAX(MIN(locate(tlb(1:ilo),time_full(j)),ilo-1),1)
               dt  = (time_full(j)-tlb(jlo+1)) / (tlb(jlo)-tlb(jlo+1))
               dt  = MAX(MIN(dt,1.0),-1.0) !no extrapolation
               zhist = (1-dt)*sfh_tab(3,jlo+1)+dt*sfh_tab(3,jlo)
               !interpolation over zhist
               klo = MAX(MIN(locate(zlegend,zhist),nz-1),1)
               dz  = (LOG10(zhist)-LOG10(zlegend(klo))) / &
                  (LOG10(zlegend(klo+1))-LOG10(zlegend(klo)))
               dz = MAX(MIN(dz,1.0),-1.0) !don't extrapolate
               ispec(:,j) = (1-dz)*spec_ssp(:,j,klo)+dz*spec_ssp(:,j,klo+1)
               ilbol(j)   = (1-dz)*lbol_ssp(j,klo)  +dz*lbol_ssp(j,klo+1)
               imass(j)   = (1-dz)*mass_ssp(j,klo)  +dz*mass_ssp(j,klo+1)
            ENDDO

         ELSE
            ispec = spec_ssp(:,:,1)
            ilbol = lbol_ssp(:,1)
            imass = mass_ssp(:,1)
         ENDIF

      ENDIF

      !set up an instantaneous burst
      IF ((pset%sfh.EQ.1.OR.pset%sfh.EQ.4).AND.&
           pset%fburst.GT.tiny_number) THEN

         IF ((powtime(i)-pset%tburst*1E9).GT.tiny_number) THEN
            delt_burst = powtime(i)-pset%tburst*1E9
            klo = MAX(MIN(locate(time_full,LOG10(delt_burst)),ntfull-1),1)
            dtb = (LOG10(delt_burst)-time_full(klo))/&
                 (time_full(klo+1)-time_full(klo))
            spec_burst = (1-dtb)*spec_ssp(:,klo,1)+dtb*spec_ssp(:,klo+1,1)
            mass_burst = (1-dtb)*mass_ssp(klo,1)  +dtb*mass_ssp(klo+1,1)
            lbol_burst = (1-dtb)*lbol_ssp(klo,1)  +dtb*lbol_ssp(klo+1,1)
         ENDIF

      ENDIF

      !compute composite spectra, mass, lbol
      IF (pset%sfh.EQ.0) THEN

         csp1 = 0.0
         csp2 = 0.0
         IF (time_full(i).LT.pset%dust_tesc) THEN
            csp1 = spec_ssp(:,i,1)
         ELSE
            csp2 = spec_ssp(:,i,1)
         ENDIF
         !add dust and combine young and old csp
         CALL ADD_DUST(pset,csp1,csp2,spec_csp,mdust)
         mass_csp = mass_ssp(i,1)
         lbol_csp = lbol_ssp(i,1)

      ELSE IF (pset%sfh.EQ.1.OR.pset%sfh.EQ.4.OR.pset%sfh.EQ.5) THEN

         CALL INTSPEC(pset,i,spec_ssp,csp1,mass_ssp,lbol_ssp,&
              mass_csp,lbol_csp,spec_burst,mass_burst,&
              lbol_burst,delt_burst,sfstart,tau,const,sftrunc,mdust)
         IF (compute_light_ages.EQ.1) THEN
            CALL INTSPEC(pset,i,spec_ssp,csp2,mass_ssp,lbol_ssp,&
                 mass_csp,lbol_csp,spec_burst,mass_burst,&
                 lbol_burst,delt_burst,sfstart,tau,const,sftrunc,mdust,1)
            spec_csp = 10**time_full(i)/1E9 - csp2/csp1 - sfstart/1E9
         ELSE
            spec_csp = csp1
         ENDIF

      ELSE IF (pset%sfh.EQ.2.OR.pset%sfh.EQ.3) THEN

         CALL INTSPEC(pset,i,ispec,spec_csp,imass,ilbol,mass_csp,&
              lbol_csp,spec_burst,mass_burst,lbol_burst,&
              delt_burst,sfstart,tau,const,sftrunc,mdust)

      ENDIF

      
      !smooth the spectrum
      IF (pset%sigma_smooth.GT.0.0) THEN
         CALL SMOOTHSPEC(spec_lambda,spec_csp,pset%sigma_smooth,&
              pset%min_wave_smooth,pset%max_wave_smooth)
      ENDIF

      !add IGM absorption
      IF (add_igm_absorption.EQ.1.AND.pset%zred.GT.0.0) THEN
         spec_csp = igm_absorb(spec_lambda,spec_csp,pset%zred,&
              pset%igm_factor)
      ENDIF

      !compute spectral indices
      IF (write_compsp.EQ.4) THEN
         CALL GETINDX(spec_lambda,spec_csp,indx)
      ELSE
         indx=0.0
      ENDIF

      !redshift spectrum; calculate mags
      IF (redshift_colors.EQ.0) THEN
         CALL GETMAGS(pset%zred,spec_csp,mags,pset%mag_compute)
      ELSE
         !here we compute the redshift at the corresponding age
         zred = MIN(MAX(linterp(cosmospl(:,2),cosmospl(:,1),&
              powtime(i)/1E9),0.0),20.0)
         CALL GETMAGS(zred,spec_csp,mags,pset%mag_compute)
      ENDIF

      !only save results if computing all ages
      IF (imax-imin.GT.1.OR.pset%tage.EQ.-99.0) THEN
         CALL SAVE_COMPSP(write_compsp,ocompsp(i),time_full(i),&
              mass_csp,lbol_csp,tsfr(i),mags,spec_csp,mdust,indx)
      ELSE
         !save results temporarily for later interpolation
         ocompsp(i)%mass_csp = mass_csp
         ocompsp(i)%lbol_csp = lbol_csp
         ocompsp(i)%mags     = mags
         ocompsp(i)%indx     = indx
         ocompsp(i)%spec     = spec_csp
      ENDIF

   ENDDO

   !interpolate to maxtime, if tage is set
   IF (imax-imin.EQ.1) THEN
      dt = (LOG10(maxtime)-time_full(imin))/&
           (time_full(imax)-time_full(imin))
      mass_csp = (1-dt)*ocompsp(imin)%mass_csp + &
           dt*ocompsp(imax)%mass_csp
      lbol_csp = (1-dt)*ocompsp(imin)%lbol_csp + &
           dt*ocompsp(imax)%lbol_csp
      DO i=1,nbands
         mags(i) = (1-dt)*ocompsp(imin)%mags(i) + &
              dt*ocompsp(imax)%mags(i)
      ENDDO
      DO i=1,nindx
         indx(i) = (1-dt)*ocompsp(imin)%indx(i) + &
              dt*ocompsp(imax)%indx(i)
      ENDDO
      spec_csp = (1-dt)*ocompsp(imin)%spec + &
           dt*ocompsp(imax)%spec
      CALL SAVE_COMPSP(write_compsp,ocompsp(1),&
           LOG10(maxtime),mass_csp,lbol_csp,zero,mags,&
           spec_csp,mdust,indx)
   ENDIF

   IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) CLOSE(10)
   IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) CLOSE(20)
   

END SUBROUTINE COMPSP

!-------------------------------------------------------------------!
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!
 
SUBROUTINE COMPSP_WARNING(maxtime,pset,nzin,write_compsp)

  !check that variables are properly set

  USE sps_vars
  IMPLICIT NONE
  INTEGER, INTENT(in) :: nzin, write_compsp
  REAL(SP), INTENT(in) :: maxtime
  TYPE(PARAMS), INTENT(in) :: pset

  !-----------------------------------------------------!

  IF (maxtime.LE.1E8.AND.pset%sfh.NE.0) THEN
     WRITE(*,*) 'COMPSP ERROR, maxtime too small:',maxtime
     STOP
  ENDIF

  !the isochrones don't go past 10**10.15 yrs, so warn the user
  !that this will be an extrapolation
  IF (maxtime.GT.10**10.2.AND.isoc_type.NE.'mist') THEN
     WRITE(*,*) 'COMPSP WARNING: log(Tmax)>10.2 yrs -'//&
          ' linear extrapolation beyond this point for log(Tmax)=:',&
          LOG10(maxtime)
  ENDIF
  IF (maxtime.GT.10**10.35.AND.isoc_type.EQ.'mist') THEN
     WRITE(*,*) 'COMPSP WARNING: log(Tmax)>10.35 yrs -'//&
          ' linear extrapolation beyond this point for log(Tmax)=:',&
          LOG10(maxtime)
  ENDIF

  !warn the user about an out-of-bounds burst component
  IF (pset%tburst*1E9.GT.maxtime.AND.pset%fburst.GT.tiny_number.AND.&
       (pset%sfh.EQ.1.OR.pset%sfh.EQ.4)) THEN
     WRITE(*,*) 'COMPSP WARNING: burst time > age of system....'//&
          ' the burst component will NOT be added'
  ENDIF

  IF (pset%tage.LT.pset%sf_start.AND.pset%tage.GT.tiny_number) THEN
     WRITE(*,*) 'COMPSP ERROR: tage<sf_start  stopping...'
     STOP
  ENDIF

  IF (pset%sf_start.LT.0.0) THEN
     WRITE(*,*) 'COMPSP ERROR: sf_start<0.  stopping...'
     STOP
  ENDIF

  IF (pset%sf_start*1E9.GT.maxtime) THEN
     WRITE(*,*) 'COMPSP ERROR: sf_start>maxtime  stopping...'
     STOP
  ENDIF

  !set limits on the parameters tau and const
  IF (pset%sfh.EQ.1.OR.pset%sfh.EQ.4) THEN
     IF (pset%tau.LE.0.1.AND.pset%tau.GE.0.0) THEN
        IF (verbose.EQ.1) THEN
           WRITE(*,*) 'COMPSP WARNING: tau <0.1, setting tau=0.1'
        ENDIF
     ELSE IF (pset%tau.GE.1E2) THEN
        IF (verbose.EQ.1) THEN
           WRITE(*,*) 'COMPSP WARNING: tau >1E2, setting tau=1E2'
        ENDIF
     ENDIF

     IF (pset%const.GT.1.0.OR.pset%const.LT.0.0) THEN
        WRITE(*,*) 'COMPSP ERROR: const out of bounds:',pset%const
        STOP
     ENDIF
  ENDIF

  IF (pset%dust_tesc.LE.5.5) THEN
     WRITE(*,*) 'COMPSP ERROR: pset%dust_tesc<=5.5, you need to set'//&
          ' dust_tesc to a value >5.5; currently it is: ',pset%dust_tesc
     STOP
  ENDIF

  IF (pset%duste_qpah.LT.0.0) THEN
     WRITE(*,*) 'COMPSP WARNING: pset%duste_qpah<0.0, '//&
          'the allowable range is 0-5 and will be set to 0.0'
  ENDIF

  IF (pset%duste_gamma.LT.0.0) THEN
     WRITE(*,*) 'COMPSP WARNING: pset%duste_gamma<0.0, '//&
          'the allowable range is >0, and will be set to 0.0'
  ENDIF
  IF ((pset%sfh.EQ.0.OR.pset%sfh.EQ.1.OR.pset%sfh.EQ.4).AND.nzin.NE.1) THEN
     WRITE(*,*) 'COMPSP_ERROR: sfh=0,1,or,4 but nzin NE 1'
     STOP
  ENDIF

  IF ((pset%sfh.EQ.2.OR.pset%sfh.EQ.3).AND.(nzin.NE.nz.AND.nzin.NE.1)) THEN
     WRITE(*,*) 'COMPSP_ERROR: sfh=2 or 3 but nzin NE (nz OR 1)'
     STOP
  ENDIF

  IF (nzin.NE.1.AND.nzin.NE.nz) THEN
     WRITE(*,*) 'COMPSP_ERROR: nzin NE 1 and nzin NE nz:',nz
     STOP
  ENDIF

  IF (write_compsp.NE.0.AND.write_compsp.NE.1 &
       .AND.write_compsp.NE.2.AND.write_compsp.NE.3 &
       .AND.write_compsp.NE.4.AND.write_compsp.NE.5) THEN
     WRITE(*,*) 'COMPSP ERROR: invalid write_compsp value:', &
          write_compsp
     STOP
  ENDIF

  IF ((pset%sfh.NE.1.AND.pset%sfh.NE.4).AND.&
       compute_light_ages.EQ.1) THEN
     WRITE(*,*) 'COMPSP ERROR: compute_light_ages only works with SFH=1 or 4'
     STOP
  ENDIF
     

END SUBROUTINE COMPSP_WARNING

!------------------------------------------------------------!
!------------------------------------------------------------!

SUBROUTINE COMPSP_SETUP_OUTPUT(write_compsp,pset,outfile,imin,imax)

  USE sps_vars
  USE sps_utils, ONLY : vactoair
  IMPLICIT NONE
  INTEGER, INTENT(in) :: imin,imax,write_compsp
  REAL(SP) :: writeage
  TYPE(PARAMS), INTENT(in) :: pset
  CHARACTER(100), INTENT(in) :: outfile

  !-----------------------------------------------------!

  !open output file for magnitudes
  IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) THEN
     OPEN(10,FILE=TRIM(SPS_HOME)//'/OUTPUTS/'//TRIM(outfile)//'.mags',&
          STATUS='REPLACE')
     CALL COMPSP_HEADER(10,pset)
  ENDIF
  
  !open output file for spectra
  IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) THEN
     OPEN(20,FILE=TRIM(SPS_HOME)//'/OUTPUTS/'//TRIM(outfile)//'.spec',&
          STATUS='REPLACE')
     CALL COMPSP_HEADER(20,pset)
  ENDIF

  !open output file for indices
  IF (write_compsp.EQ.4) THEN
     OPEN(30,FILE=TRIM(SPS_HOME)//'/OUTPUTS/'//TRIM(outfile)//'.indx',&
          STATUS='REPLACE')
     CALL COMPSP_HEADER(30,pset)
  ENDIF

  IF (pset%sfh.EQ.0) THEN
     IF (verbose.NE.0) WRITE(*,*) '  Processing SSP'
     IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) THEN
        WRITE(10,'("#   Processing SSP")')
        WRITE(10,'("#")') 
        WRITE(10,32) 
     ENDIF
     IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) THEN
        WRITE(20,'("#   Processing SSP")')
        WRITE(20,'("#")') 
        WRITE(20,31) 
        IF (imax-imin.EQ.1) WRITE(20,'(I3,1x,I6)') 1,nspec
        IF (imax-imin.GT.1) WRITE(20,'(I3,1x,I6)') ntfull,nspec
        IF (vactoair_flag.EQ.0) THEN
           WRITE(20,'(50000(F15.4))') spec_lambda
        ELSE
           WRITE(20,'(50000(F15.4))') vactoair(spec_lambda)
        ENDIF
     ENDIF
     IF (write_compsp.EQ.4) THEN
        WRITE(30,'("#   Processing SSP")')
        WRITE(30,'("#")') 
        WRITE(30,34) 
     ENDIF
  ELSE
     IF (pset%sfh.EQ.2.OR.pset%sfh.EQ.3) THEN
        IF (verbose.EQ.1) &
             WRITE(*,30) pset%dust1,pset%dust2
        IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) &
             WRITE(10,30) pset%dust1,pset%dust2
        IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) &
             WRITE(20,30) pset%dust1,pset%dust2
     ELSE
        IF (pset%tage.GT.tiny_number) writeage = pset%tage
        IF (pset%tage.LE.tiny_number) writeage = 10**time_full(ntfull)/1E9
        IF (verbose.EQ.1) &
             WRITE(*,33) writeage,LOG10(pset%tau),pset%const,pset%fburst,&
             pset%tburst,pset%sf_start,pset%dust1,pset%dust2
        IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) &
             WRITE(10,33) writeage,LOG10(pset%tau),pset%const,pset%fburst,&
             pset%tburst,pset%sf_start,pset%dust1,pset%dust2
        IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) &
             WRITE(20,33) writeage,LOG10(pset%tau),pset%const,pset%fburst,&
             pset%tburst,pset%sf_start,pset%dust1,pset%dust2
        IF (write_compsp.EQ.4) &
             WRITE(30,33) writeage,LOG10(pset%tau),pset%const,pset%fburst,&
             pset%tburst,pset%sf_start,pset%dust1,pset%dust2
     ENDIF
     IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) THEN 
        WRITE(10,'("#")') 
        WRITE(10,32) 
     ENDIF
     IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) THEN
        WRITE(20,'("#")') 
        WRITE(20,31) 
        IF (imax-imin.EQ.1) WRITE(20,'(I3,1x,I6)') 1,nspec
        IF (imax-imin.GT.1) WRITE(20,'(I3,1x,I6)') ntfull,nspec
        IF (vactoair_flag.EQ.0) THEN
           WRITE(20,'(50000(F15.4))') spec_lambda
        ELSE
           WRITE(20,'(50000(F15.4))') vactoair(spec_lambda)
        ENDIF
       ENDIF
       IF (write_compsp.EQ.4) THEN
          WRITE(20,'("#")') 
          WRITE(20,34) 
       ENDIF
   ENDIF

   !formats
30 FORMAT('#   SFH: tabulated input, dust=(',F6.2,','F6.2,')')
31 FORMAT('#   log(age) log(mass) Log(lbol) log(SFR) spectra')
32 FORMAT('#   log(age) log(mass) Log(lbol) log(SFR) mags (see FILTER_LIST)')
33 FORMAT('#   SFH: Tage=',F6.2,' Gyr, log(tau/Gyr)= ',F6.3,&
        ', const= ',F6.3,', fb= ',F6.3,', tb= ',F6.2,&
        ' Gyr, sf_start= 'F6.3,', dust=(',F6.2,','F6.2,')')
34 FORMAT('#   log(age) indices (see allindices.dat)')


END SUBROUTINE COMPSP_SETUP_OUTPUT

!------------------------------------------------------------!
!------------------------------------------------------------!

SUBROUTINE COMPSP_HEADER(unit,pset)

  !writes headers for the .mag, .spec, .indx files

  USE sps_vars
  IMPLICIT NONE
  INTEGER, INTENT(in) :: unit
  TYPE(PARAMS), INTENT(in) :: pset

  !-----------------------------------------------------!

  IF (pset%sfh.NE.2) THEN
     WRITE(unit,'("#   Log(Z/Zsol): ",F6.3)') &
          LOG10(zlegend(pset%zmet)/zsol)
  ELSE
     WRITE(unit,'("#   Log(Z/Zsol): tabulated")')
  ENDIF
  WRITE(unit,'("#   Fraction of blue HB stars: ",F6.3,'//&
       '"; Ratio of BS to HB stars: ",F6.3)') pset%fbhb, pset%sbss
  WRITE(unit,'("#   Shift to TP-AGB [log(Teff),log(Lbol)]: ",F5.2,1x,F5.2)') &
       pset%delt, pset%dell
  IF (imf_type.EQ.2) THEN
     WRITE(unit,'("#   IMF: ",I1,", slopes= ",3F4.1)') &
          imf_type,pset%imf1,pset%imf2,pset%imf3
  ELSE IF (imf_type.EQ.3) THEN
     WRITE(unit,'("#   IMF: ",I1,", cut-off= ",F4.2)') imf_type,pset%vdmc
  ELSE
     WRITE(unit,'("#   IMF: ",I1)') imf_type
  ENDIF
  IF (compute_vega_mags.EQ.1) THEN
     WRITE(unit,'("#   Mag Zero Point: Vega (not relevant for spec/indx files)")')
  ELSE
     WRITE(unit,'("#   Mag Zero Point: AB (not relevant for spec/indx files)")')
  ENDIF

END SUBROUTINE COMPSP_HEADER

!------------------------------------------------------------!
!------------------------------------------------------------!

SUBROUTINE SAVE_COMPSP(write_compsp,cspo,time,mass,&
     lbol,sfr,mags,spec,mdust,indx)

  !routine to print and save outputs

  USE sps_vars
  IMPLICIT NONE
  INTEGER, INTENT(in) :: write_compsp
  REAL(SP), INTENT(in)    :: time,mass,lbol,sfr,mdust
  REAL(SP), DIMENSION(nspec), INTENT(in)  :: spec
  REAL(SP), DIMENSION(nbands), INTENT(in) :: mags
  REAL(SP), DIMENSION(nindx), INTENT(in)  :: indx
  TYPE(COMPSPOUT), INTENT(inout) :: cspo
  CHARACTER(34) :: fmt

  !-----------------------------------------------------!

  fmt = '(F7.4,1x,3(F8.4,1x),000(F7.3,1x))'
  WRITE(fmt(21:23),'(I3,1x,I4)') nbands

  !dump info into output structure
  cspo%age      = time
  cspo%mass_csp = mass
  cspo%lbol_csp = lbol
  cspo%sfr      = sfr
  cspo%mags     = mags
  cspo%spec     = spec
  cspo%mdust    = mdust
  cspo%indx     = indx

  !write to mags file
  IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) &
       WRITE(10,fmt) time,LOG10(mass+tiny_number),&
       lbol,LOG10(sfr+tiny_number),mags
 
  !write to spectra file
  IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) THEN
     WRITE(20,'(4(F8.4,1x))') time,&
          LOG10(mass+tiny_number),lbol,LOG10(sfr+tiny_number)
     WRITE(20,'(50000(E14.6))') spec
  ENDIF

  !write to indx file
  IF (write_compsp.EQ.4) &
       WRITE(30,'(F8.4,99(F7.3,1x))') time,indx
  
END SUBROUTINE SAVE_COMPSP
