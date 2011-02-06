!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Routine to compute magnitudes and spectra for a composite !    
!  stellar population.  Returns a file with the following    !
!  info: age, mass, Lbol, mags in various filters.  SFR      !
!  units are Msun/yr.                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INTSPEC(pset,nti,wtesc,spec_ssp,sfr,&
     csp,mass_ssp,lbol_ssp,mass,lbol,time,specb,&
     massb,lbolb,deltb)

  !routine to perform integration of SSP over a SFH,
  !including attenuation by dust.

  USE nrtype; USE sps_vars
  USE sps_utils, ONLY : add_dust
  IMPLICIT NONE
  
  INTEGER :: i, ntii,imax
  INTEGER,  INTENT(in)    :: nti,wtesc
  REAL(SP), INTENT(in)    :: massb,lbolb,deltb
  REAL(SP), INTENT(inout) :: mass, lbol
  REAL(SP), INTENT(in), DIMENSION(ntfull) :: mass_ssp,lbol_ssp,sfr,time
  REAL(SP), INTENT(in), DIMENSION(ntfull,nspec) :: spec_ssp
  REAL(SP), INTENT(in), DIMENSION(nspec)    :: specb
  REAL(SP), INTENT(inout), DIMENSION(nspec) :: csp
  REAL(SP), DIMENSION(ntfull) :: weights
  REAL(SP), DIMENSION(nspec) :: csp1,csp2
  TYPE(PARAMS), INTENT(in)   :: pset

  !-----------------------------------------------------!
  !-------------Compute CSP w/o dust--------------------!
  !-----------------------------------------------------!

  !t<tesc
  csp1 = 0.
  imax = MAX(MIN(wtesc,nti),2)
  DO i=1,imax-1
     csp1 = csp1 + (time(i+1)-time(i)) * &
          0.5*(sfr(i+1)*spec_ssp(i+1,:)+sfr(i)*spec_ssp(i,:))
  ENDDO

  !t>tesc
  csp2 = 0.0
  IF (nti.GT.wtesc) THEN
     DO i=wtesc,nti-1
        csp2 = csp2 + (time(i+1)-time(i)) * &
             0.5*(sfr(i+1)*spec_ssp(i+1,:)+sfr(i)*spec_ssp(i,:))
     ENDDO
     ntii = nti
  ELSE
     ntii = imax
  ENDIF

  !compute weighted mass and lbol
  mass = 0.0
  lbol = 0.0
  DO i=1,ntii-1
     mass = mass + (time(i+1)-time(i))*0.5* &
          (sfr(i+1)*mass_ssp(i+1) + &
          sfr(i)*mass_ssp(i))
     lbol = lbol + (time(i+1)-time(i))*0.5* &
          (sfr(i+1)*10**lbol_ssp(i+1) + &
          sfr(i)*10**lbol_ssp(i))
  ENDDO
  lbol = LOG10(lbol)
     
  !add in an instantaneous burst
  IF (pset%fburst.GT.0.0.AND.(pset%sfh.EQ.1.OR.pset%sfh.EQ.4)) THEN
     
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
  CALL ADD_DUST(pset,csp1,csp2,csp)

END SUBROUTINE INTSPEC

!------------------------------------------------------------!
!---------------------Main Routine---------------------------!
!------------------------------------------------------------!

SUBROUTINE COMPSP(write_compsp,nzin,outfile,mass_ssp,&
     lbol_ssp,spec_ssp,pset,ocompsp)

  !We want 1 solar mass of stars to have formed between 
  !Tage and Tuniv except or tabulated SFHs, where the 
  !user specifies the normalization.

  !for sfh=1:
  !If tage >0  -> run only one integration to t=tage
  !If tage<=0  -> produce outputs from tmin<t<tuniv

  USE sps_vars; USE nrtype
  USE sps_utils, ONLY : getmags, add_dust
  USE nrutil, ONLY : assert_eq 
  USE nr, ONLY : locate, splint
  IMPLICIT NONE
 
  !specify type of I/O; number of metallicity inputs
  !write = (1->write mags), (2->write spectra), (3->write both) 
  INTEGER, INTENT(in)                          :: write_compsp,nzin
  !time-dependent mass and Lbol for the SSP
  REAL(SP), INTENT(in), DIMENSION(nzin,ntfull) :: lbol_ssp,mass_ssp
  !time-dependent spectrum for SSP
  REAL(SP), INTENT(in), DIMENSION(nzin,ntfull,nspec) :: spec_ssp
  CHARACTER(100), INTENT(in) :: outfile

  INTEGER  :: i=0,j=1,n=1,stat=0,klo,jlo,wtesc,imin,imax
  REAL(SP) :: tau,const,maxtime,t,writeage=0.,psfr
  REAL(SP) :: mass_csp, lbol_csp,dtb,dt,zred=0.
  REAL(SP) :: mass_burst=0.0,lbol_burst=0.0,delt_burst=0.0
  REAL(SP), DIMENSION(nbands)       :: mags=0.0
  REAL(SP), DIMENSION(ntfull,nspec) :: ispec=0.0
  REAL(SP), DIMENSION(nspec)        :: csp_spec=0.0,spec_burst=0.0,csp1,csp2
  REAL(SP), DIMENSION(ntfull) :: imass=0.0, ilbol=0.0
  REAL(SP), DIMENSION(ntfull) :: powtime=0.
  REAL(SP), DIMENSION(ntfull) :: sfr=0.,tsfr=0.,tzhist=0.,zhist=0.
  CHARACTER(33) :: fmt

  !(TYPE objects defined in sps_vars.f90)
  TYPE(PARAMS), INTENT(in) :: pset
  TYPE(COMPSPOUT), INTENT(inout), DIMENSION(ntfull) :: ocompsp

  !-----------------------------------------------------!
  !--------------------Basic Setup----------------------!
  !-----------------------------------------------------!

  IF (check_sps_setup.EQ.0) THEN
     WRITE(*,*) 'COMPSP ERROR0: '//&
          'SPS_SETUP must be run before calling COMPSP. '
     STOP
  ENDIF
 
  powtime = 10**time_full
  maxtime = powtime(ntfull)
  imin = 1
  imax = ntfull

  !if a tabulated SFH is passed, read in sfh.dat file
  IF (pset%sfh.EQ.2) THEN

     OPEN(3,FILE='sfh.dat',ACTION='READ',STATUS='OLD')
     DO n=1,10000
        READ(3,*,IOSTAT=stat) sfh_tab(1,n),sfh_tab(2,n),sfh_tab(3,n)
        IF (stat.NE.0) GOTO 29
     ENDDO
     WRITE(*,*) 'COMPSP ERROR: didnt finish reading in the sfh file'
     STOP
29   CONTINUE
     CLOSE(3)
     ntabsfh = n-1

     sfh_tab(1,1:ntabsfh) = sfh_tab(1,1:ntabsfh)*1E9

  ELSE IF (pset%sfh.EQ.3) THEN 

     !sfh_tab array is supposed to already be filled in, check that it is
     IF (ntabsfh.EQ.0) THEN 
        WRITE(*,*) 'COMPSP ERROR: sfh=3 but sfh_tab array not initialized!'
        STOP
     ENDIF

  ELSE

     !set up maxtime variable
     IF (pset%tage.GT.0) THEN
        maxtime = pset%tage*1E9
        imin = locate(powtime,maxtime)
        imax = imin+1
     ENDIF

  ENDIF

  !find where t~dust_tesc
  wtesc = MAX(locate(time_full,pset%dust_tesc),1)

  !set limits on the parameters tau and const
  IF (pset%sfh.EQ.1.OR.pset%sfh.EQ.4) THEN
     tau   = MAX(MIN(pset%tau,0.1),100.)
     const = pset%const
  ENDIF

  !make sure various variables are set correctly
  CALL COMPSP_WARNING(maxtime,pset,nzin,write_compsp)

  !tsfr is SFR as a fcn of age of the Universe
  !tzhist is metallicity history as a fcn of age of Universe
  !for analytic SFHs, these are only used for output
  IF (pset%sfh.EQ.2.OR.pset%sfh.EQ.3) THEN
     DO j=1,ntfull
        jlo    = MAX(MIN(locate(LOG10(sfh_tab(1,1:ntabsfh)),&
             time_full(j)),ntabsfh-1),1)
        IF (time_full(j).GT.LOG10(sfh_tab(1,ntabsfh))) THEN
           tsfr(j) = 0.0
        ELSE
           dt = (10**time_full(j)-(sfh_tab(1,jlo))) / &
                (sfh_tab(1,jlo+1)-sfh_tab(1,jlo))
           tsfr(j)   = (1-dt)*sfh_tab(2,jlo)+dt*sfh_tab(2,jlo+1)
           tzhist(j) = (1-dt)*sfh_tab(3,jlo)+dt*sfh_tab(3,jlo+1)
        ENDIF
     ENDDO
  ELSE IF (pset%sfh.EQ.1) THEN
     tsfr  = EXP(-powtime/tau/1E9 )/tau/1E9 / &
          (1-EXP(-maxtime/1E9/tau))*(1-const) &
          + const/maxtime
  ELSE IF (pset%sfh.EQ.4) THEN
     tsfr  = (powtime/tau/1E9)*EXP(-powtime/tau/1E9 )/tau/1E9 / &
          (1-EXP(-maxtime/1E9/tau)*(maxtime/1E9/tau+1))*(1-const) + &
          const/maxtime
  ELSE IF (pset%sfh.EQ.0) THEN
     tsfr = 0.0
  ENDIF

  !-----------------------------------------------------!
  !---------------Setup output files--------------------!
  !-----------------------------------------------------!

  fmt = '(F7.4,1x,3(F8.4,1x),00(F7.3,1x))'
  WRITE(fmt(21:22),'(I2)') nbands

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
        WRITE(20,'(I3)') ntfull
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
        IF (pset%tage.GT.0.) writeage = pset%tage
        IF (pset%tage.LE.0.) writeage = tuniv
        IF (verbose.EQ.1) &
             WRITE(*,33) writeage,LOG10(tau),const,pset%fburst,&
             pset%tburst,pset%dust1,pset%dust2
        IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) &
             WRITE(10,33) writeage,LOG10(tau),const,pset%fburst,&
             pset%tburst,pset%dust1,pset%dust2
        IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) &
             WRITE(20,33) writeage,LOG10(tau),const,pset%fburst,&
             pset%tburst,pset%dust1,pset%dust2
     ENDIF
     IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) THEN 
        WRITE(10,'("#")') 
        WRITE(10,32) 
     ENDIF
     IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) THEN
        WRITE(20,'("#")') 
        WRITE(20,31) 
        WRITE(20,'(I3)') ntfull
       ENDIF
   ENDIF

  !-----------------------------------------------------!
  !---------Generate composite spectra and mags---------!
  !-----------------------------------------------------!

   !calculate mags at each time step
   DO i=imin,imax

      !Set up analytic SFH, burst, tabulated SFH, etc.
      IF (pset%sfh.EQ.2.OR.pset%sfh.EQ.3) THEN
         
         sfr   = 0.0
         ilbol = 0.0
         imass = 0.0
         ispec = 0.0
         DO j=1,i
            jlo = MAX(MIN(locate(powtime(i)-powtime(1:i),&
                 powtime(j)),ntfull-1),1)
            dt = (powtime(j)-(powtime(i)-powtime(jlo))) / &
                 (powtime(jlo+1)-powtime(jlo))
            sfr(j)   = MAX((1-dt)*tsfr(jlo)+dt*tsfr(jlo+1),0.0)
            zhist(j) = (1-dt)*tzhist(jlo)+dt*tzhist(jlo+1)
            
            !interpolation over zhist
            klo = MAX(MIN(locate(zlegend,zhist(j)),nz-1),1)
            t   = (zhist(j)-zlegend(klo)) / &
                 (zlegend(klo+1)-zlegend(klo))
            ispec(j,:) = (1-t)*spec_ssp(klo,j,:)+t*spec_ssp(klo+1,j,:)
            ilbol(j)   = (1-t)*lbol_ssp(klo,j)+t*lbol_ssp(klo+1,j)
            imass(j)   = (1-t)*mass_ssp(klo,j)+t*mass_ssp(klo+1,j)
         ENDDO
         
      ELSE IF (pset%sfh.EQ.1.OR.pset%sfh.EQ.4) THEN
         
         sfr = 0.0
         IF (pset%sfh.EQ.1) &
              sfr(1:i) = EXP(-(powtime(i)-powtime(1:i))/tau/1E9) / &
              tau/1E9/(1-EXP(-maxtime/1E9/tau)) * (1-const) &
              + const/maxtime
         IF (pset%sfh.EQ.4) &
              sfr(1:i)  = ((powtime(i)-powtime(1:i))/tau/1E9)*&
              EXP(-(powtime(i)-powtime(1:i))/tau/1E9 )/tau/1E9 / &
              (1-EXP(-maxtime/1E9/tau)*(maxtime/1E9/tau+1)) * (1-const) &
              + const/maxtime
         
         !set up an instantaneous burst
         IF ((powtime(i)-pset%tburst*1E9).GT.0.0.AND.&
              pset%fburst.GT.0.0) THEN            
            delt_burst = powtime(i)-pset%tburst*1E9
            klo = MAX(MIN(locate(time_full,LOG10(delt_burst)),ntfull-1),1)
            dtb = (LOG10(delt_burst)-time_full(klo))/&
                 (time_full(klo+1)-time_full(klo))
            spec_burst = (1-dtb)*spec_ssp(nzin,klo,:)+dtb*spec_ssp(nzin,klo+1,:)
            mass_burst = (1-dtb)*mass_ssp(nzin,klo)+dtb*mass_ssp(nzin,klo+1)
            lbol_burst = (1-dtb)*lbol_ssp(nzin,klo)+dtb*lbol_ssp(nzin,klo+1)
         ENDIF

      ENDIF

      IF (pset%sfh.EQ.0) THEN

         csp1 = 0.0
         csp2 = 0.0
         IF (time_full(i).LT.pset%dust_tesc) THEN
            csp1 = spec_ssp(nzin,i,:)
         ELSE
            csp2 = spec_ssp(nzin,i,:)
         ENDIF
         !add dust and combine young and old csp
         CALL ADD_DUST(pset,csp1,csp2,csp_spec)
         mass_csp = mass_ssp(nzin,i)
         lbol_csp = lbol_ssp(nzin,i)
         
      ELSE IF (pset%sfh.EQ.1.OR.pset%sfh.EQ.4) THEN
         
         CALL INTSPEC(pset,i,wtesc,spec_ssp,sfr,&
              csp_spec,mass_ssp,lbol_ssp,mass_csp,lbol_csp,powtime,&
              spec_burst,mass_burst,lbol_burst,delt_burst)
         
      ELSE IF (pset%sfh.EQ.2.OR.pset%sfh.EQ.3) THEN

         CALL INTSPEC(pset,i,wtesc,ispec,sfr,&
              csp_spec,imass,ilbol,mass_csp,lbol_csp,powtime,&
              spec_burst,mass_burst,lbol_burst,delt_burst)
         
      ENDIF
      
      !calculate mags at each timestep
      IF (redshift_colors.EQ.0) THEN
         CALL GETMAGS(pset%zred,csp_spec,mags)
      ELSE
         !here we compute the redshift at the corresponding
         !age.  Allows for computation of evolution of 
         !observed magnitudes with redshift
         zred = splint(zagespl(:,2),zagespl(:,1),&
              zagespl(:,3),10**time_full(i)/1E9)
         zred = MIN(MAX(zred,0.0),20.0)
         CALL GETMAGS(zred,csp_spec,mags)
      ENDIF
      
      !dump info into output structure
      ocompsp(i)%age      = time_full(i)
      ocompsp(i)%mass_csp = mass_csp
      ocompsp(i)%lbol_csp = lbol_csp
      ocompsp(i)%sfr      = tsfr(i)
      ocompsp(i)%mags     = mags
      ocompsp(i)%spec     = csp_spec
      
      psfr = LOG10(tsfr(i)+1E-99_dp)
      
      !write to mags file
      IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) &
           WRITE(10,fmt) time_full(i),LOG10(mass_csp),&
           lbol_csp,psfr,mags
      
      !write to spectra file
      IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) THEN
         WRITE(20,'(4(F8.4,1x))') time_full(i),&
              LOG10(mass_csp),lbol_csp,psfr
         WRITE(20,*) csp_spec
      ENDIF
      
   ENDDO

   IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) CLOSE(10)
   IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) CLOSE(20)
   
   !formats
30 FORMAT('#   SFH: tabulated input, dust=(',F6.2,','F6.2,')')
31 FORMAT('#   log(age) log(mass) Log(lbol) log(SFR) spectra')
32 FORMAT('#   log(age) log(mass) Log(lbol) log(SFR) mags (see FILTER_LIST)')
33 FORMAT('#   SFH: Tage=',F6.2,' Gyr, log(tau)= ',F6.3,&
        ' Gyr, const= ',F5.3,', fb= ',F5.3,', tb= ',F6.2,&
        ' Gyr, dust=(',F6.2,','F6.2,')')

END SUBROUTINE COMPSP 

!------------------------------------------------------------!
!------------------------------------------------------------!

SUBROUTINE COMPSP_WARNING(maxtime, pset, nzin, write_compsp)

  USE sps_vars
  IMPLICIT NONE
  INTEGER, INTENT(in) :: nzin, write_compsp
  REAL, INTENT(in) :: maxtime
  TYPE(PARAMS), INTENT(inout) :: pset

  !-----------------------------------------------------!

  IF (maxtime.LE.1E8) THEN
     WRITE(*,*) 'COMPSP ERROR, maxtime too small:',maxtime
     STOP
  ENDIF

  IF (tuniv.LE.0.0) THEN
     WRITE(*,*) 'COMPSP ERROR: tuniv<=0.0'
     STOP
  ENDIF

  !the isochrones don't go past 10**10.15 yrs, so warn the user
  !that this will be an extrapolation
  IF (maxtime.GT.10**10.20) THEN
     WRITE(*,*) 'COMPSP WARNING: Tmax>10^10.2 yrs -'//&
          ' linear extrapolation beyond this point for log(Tmax)=:',&
          LOG10(maxtime)
  ENDIF

  !warn the user about an out-of-bounds burst component
  IF (pset%tburst*1E9.GT.maxtime.AND.pset%fburst.GT.0.0.AND.&
       (pset%sfh.EQ.1.OR.pset%sfh.EQ.4)) THEN
     WRITE(*,*) 'COMPSP WARNING: burst time > age of system....'//&
          ' the burst component will NOT be added'
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

  IF (pset%sfh.EQ.0.AND.nzin.NE.1) THEN
     WRITE(*,*) 'COMPSP_ERROR: sfh=0 but nzin NE 1'
     STOP
  ENDIF

  IF (nzin.NE.1.AND.nzin.NE.nz) THEN
     WRITE(*,*) 'COMPSP_ERROR: nzin NE 1 and nzin NE nz:',nz
     STOP
  ENDIF

  IF (write_compsp.NE.1.AND.write_compsp.NE.2 &
       .AND.write_compsp.NE.3.AND.write_compsp.NE.0) THEN
     WRITE(*,*) 'COMPSP ERROR: invalid write_compsp value:', &
          write_compsp
     STOP
  ENDIF

END SUBROUTINE COMPSP_WARNING

!------------------------------------------------------------!
!------------------------------------------------------------!

SUBROUTINE COMPSP_HEADER(unit,pset)

  USE sps_vars
  IMPLICIT NONE
  INTEGER, INTENT(in) :: unit
  TYPE(PARAMS), INTENT(inout) :: pset

  !-----------------------------------------------------!

  IF (pset%sfh.NE.2) THEN
     WRITE(unit,'("#   Log(Z/Zsol): ",F6.3)') &
          LOG10(zlegend(pset%zmet)/0.0190)
  ELSE
     WRITE(unit,'("#   Log(Z/Zsol): tabulated")')
  ENDIF
  WRITE(unit,'("#   Fraction of blue HB stars: ",F6.3)') pset%fbhb
  WRITE(unit,'("#   Ratio of BS to HB stars  : ",F6.3)') pset%sbss
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

END SUBROUTINE COMPSP_HEADER
