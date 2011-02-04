!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Routine to compute magnitudes and spectra for a composite !    
!  stellar population.  Returns a file with the following    !
!  info: age, mass, Lbol, mags in various filters.  SFR      !
!  units are Msun/yr.                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INTSPEC(pset,nti,wtesc,dlogt,ssp,sfr,&
     csp,mass_tmp,lbol_tmp,mass,lbol,time,spec_burst,&
     mass_burst,lbol_burst,delt_burst)

  !routine to perform integration of SSP over a SFH,
  !including attenuation by dust.

  USE nrtype; USE sps_vars
  USE sps_utils, ONLY : add_dust
  IMPLICIT NONE
  
  INTEGER :: i
  INTEGER,  INTENT(in)    :: nti,wtesc
  REAL(SP), INTENT(in)    :: dlogt, mass_burst,lbol_burst,delt_burst
  REAL(SP), INTENT(inout) :: mass, lbol
  REAL(SP), INTENT(in), DIMENSION(ntall) :: mass_tmp,lbol_tmp,sfr,time
  REAL(SP), INTENT(in), DIMENSION(ntall,nspec) :: ssp
  REAL(SP), INTENT(in), DIMENSION(nspec)    :: spec_burst
  REAL(SP), INTENT(inout), DIMENSION(nspec) :: csp
  REAL(SP), DIMENSION(ntall) :: weights
  REAL(SP), DIMENSION(nspec) :: csp1,csp2
  TYPE(PARAMS), INTENT(in)   :: pset

  !-----------------------------------------------------!
  !-------------Compute CSP w/o dust--------------------!
  !-----------------------------------------------------!

  !t<tesc
  csp1 = 0.
  IF (wtesc.GT.1) THEN
     IF (wtesc.GE.6) THEN
        !do 4th order integration
        weights                = 1.
        weights(1:3)           = (/3./8,7./6,23./24/)
        weights(wtesc-2:wtesc) = (/23./24,7./6,3./8/)
        DO i=1,wtesc
           csp1 = csp1 + 2.3026*dlogt * &
                sfr(i)*ssp(i,:)*weights(i)*time(i)
        ENDDO
     ELSE
        !do trapezoidal integration
        DO i=1,wtesc-1
           csp1 = csp1 + 2.3026*dlogt * &
                0.5*(sfr(i+1)*ssp(i+1,:)*time(i+1) + &
                sfr(i)*ssp(i,:)*time(i))
        ENDDO
     ENDIF
  ENDIF

  !t>tesc, 4th order integration
  csp2 = 0.0
  weights                = 1.
  weights(wtesc:wtesc+2) = (/3./8,7./6,23./24/)
  weights(nti-2:nti)     = (/23./24,7./6,3./8/)
  DO i=wtesc,nti
     csp2 = csp2 + 2.3026*dlogt * &
          sfr(i)*ssp(i,:)*weights(i)*time(i)
  ENDDO
  
  !compute weighted mass and lbol (4th order integration)
  weights            = 1.
  weights(1:3)       = (/3./8,7./6,23./24/)
  weights(nti-2:nti) = (/23./24,7./6,3./8/)
  mass =        2.3026*dlogt*SUM(sfr*mass_tmp*weights*time)
  lbol = LOG10( 2.3026*dlogt*SUM(sfr*10**lbol_tmp*weights*time) )

  !add in an instantaneous burst
  IF (pset%fburst.GT.0.0.AND.pset%sfh.EQ.1) THEN
     
     IF (delt_burst.LT.10**pset%dust_tesc) THEN
        csp1 = csp1*(1-pset%fburst) + spec_burst*pset%fburst
        csp2 = csp2*(1-pset%fburst)
     ELSE
        csp1 = csp1*(1-pset%fburst)
        csp2 = csp2*(1-pset%fburst)   + spec_burst*pset%fburst
     ENDIF

     mass = mass*(1-pset%fburst) + mass_burst*pset%fburst
     lbol = LOG10( 10**lbol*(1-pset%fburst) + 10**lbol_burst*pset%fburst )
     
  ENDIF

  !add dust
  CALL ADD_DUST(pset,csp1,csp2,csp)

END SUBROUTINE INTSPEC

!-------------------------------------------------------!
!-------------------------------------------------------!
!-------------------------------------------------------!

SUBROUTINE COMPSP(write_compsp,nzin,outfile,mass_ssp,&
     lbol_ssp,spec_ssp,pset,ocompsp)

  !results have been compared to BC03 and M05 models for 
  !SSPs and CSPs with tau and dust.

  !We want 1 solar mass of stars to have formed between 
  !Tage and Tuniv except or tabulated SFHs, where the 
  !user specifies the normalization.

  !for sfh=1:
  !If tage >0  -> run only one integration to t=tage
  !If tage<=0  -> produce outputs from tmin<t<tuniv

  !this code currently only includes metallicity histories for
  !tabulated SFHs

  USE sps_vars; USE nrtype; USE sps_utils, ONLY : getmags, add_dust
  USE nrutil, ONLY : assert_eq; USE nr, ONLY : locate, splint
  IMPLICIT NONE
 
  !specify type of I/O; number of metallicity inputs
  !write = (1->write mags), (2->write spectra), (3->write both) 
  INTEGER, INTENT(in)                            :: write_compsp,nzin
  !time-dependent mass and Lbol for the SSP
  REAL(SP), INTENT(in), DIMENSION(nzin,nt)       :: lbol_ssp,mass_ssp
  !time-dependent spectrum for SSP
  REAL(SP), INTENT(in), DIMENSION(nzin,nt,nspec) :: spec_ssp
  CHARACTER(100), INTENT(in)                     :: outfile

  INTEGER  :: i=0,j=1,n=1,stat=0,imin,klo,jlo
  INTEGER  :: imax=1,wtesc,istep
  REAL(SP) :: tau=-1.,const=-1.,maxtime=0.0,t,u,writeage=0.,psfr,zred=0.
  REAL(SP) :: mass_csp, lbol_csp, dlogt, dt2_burst=0.
  REAL(SP) :: mass_burst=0.0,lbol_burst=0.0,delt_burst=0.0
  REAL(SP), DIMENSION(nbands)      :: mags=0.0
  REAL(SP), DIMENSION(ntall,nspec) :: ispec=0.0
  REAL(SP), DIMENSION(nspec)       :: csp_spec=0.0,spec_burst=0.0,csp1,csp2
  REAL(SP), DIMENSION(ntall)   :: mass_ssp2=0.0, lbol_ssp2=0.0
  REAL(SP), DIMENSION(nt)      :: time=0.
  REAL(SP), DIMENSION(ntall)   :: itime=0.,powitime=0.
  REAL(SP), DIMENSION(ntall)   :: sfr=0.,tsfr=0.,tzhist=0.,zhist=0.
  CHARACTER(33) :: fmt

  !(TYPE objects defined in sps_vars.f90)
  TYPE(PARAMS), INTENT(in) :: pset
  TYPE(COMPSPOUT), INTENT(inout), DIMENSION(nt) :: ocompsp

  !-----------------------------------------------------!
  !----------Setup vars and read in sfh.dat-------------!
  !-----------------------------------------------------!

  IF (check_sps_setup.EQ.0) THEN
     WRITE(*,*) 'COMPSP ERROR0: '//&
          'SPS_SETUP must be run before calling COMPSP. '
     STOP
  ENDIF
 
  sfr  = 0.0
  tsfr = 0.0

  IF (write_compsp.NE.1.AND.write_compsp.NE.2 &
       .AND.write_compsp.NE.3.AND.write_compsp.NE.0) THEN
     WRITE(*,*) 'COMPSP ERROR: invalid write_compsp value:', &
          write_compsp
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

  !to balance numerical accuracy and speed, the resolution
  !of the time grid must be a function of tau
  IF (pset%sfh.EQ.1) THEN
     ntime = MIN(MAX(NINT(10**(0.3-LOG10(pset%tau))),3),15*20)
  ELSE
     ntime = 15
  ENDIF

  !if a tabulated SFH is passed, read in sfh.dat file
  IF (pset%sfh.EQ.2) THEN

     OPEN(3,FILE='sfh.dat',ACTION='READ',STATUS='OLD')
     DO n=1,1E4 
        READ(3,*,IOSTAT=stat) sfh_tab(1,n),sfh_tab(2,n),sfh_tab(3,n)
        IF (stat.NE.0) GOTO 29
     ENDDO
     WRITE(*,*) 'COMPSP ERROR: didnt finish reading in the sfh file'
     STOP
29   CONTINUE
     CLOSE(3)
     ntabsfh = n-1

     sfh_tab(1,1:ntabsfh) = sfh_tab(1,1:ntabsfh)*1E9   !convert to years
     maxtime = sfh_tab(1,ntabsfh)

  ELSE IF (pset%sfh.EQ.3) THEN 

     !assume sfh_tab array has already been filled
     IF (ntabsfh.EQ.0) THEN 
        WRITE(*,*) 'COMPSP ERROR: sfh=3 but sfh_tab array not initialized!'
        STOP
     ENDIF
     maxtime = sfh_tab(1,ntabsfh)

  ELSE

     !set up maxtime variable
     IF (pset%tage.LE.0) THEN
        IF (tuniv.LE.0.0) THEN
           WRITE(*,*) 'COMPSP ERROR: tuniv<=0.0'
           STOP
        ENDIF
        maxtime = tuniv*1E9
     ELSE
        maxtime = pset%tage*1E9
        IF (maxtime.LE.1E8) THEN
           WRITE(*,*) 'COMPSP ERROR, maxtime too small:',maxtime
           STOP
        ENDIF
     ENDIF

  ENDIF

  !----------------------------------------------------------!
  !-----------------Set up the time array--------------------!
  !----------------------------------------------------------!
  
  !transfer age array into temp var, units of log(yrs)
  time  = timestep_isoc(pset%zmet,:)

  !set up expanded time array
  DO i=1,nt*ntime
     itime(i) = 5.5+(i-1)*(LOG10(maxtime)-5.5)/(ntime*nt-1)
  ENDDO
  dlogt    = itime(2)-itime(1)
  powitime = 10**itime
  imax     = nt*ntime

  !find where t~dust_tesc
  wtesc = MAX(locate(itime(1:nt*ntime),pset%dust_tesc),1)

  IF (pset%dust_tesc.LT.5.5) THEN
     WRITE(*,*) 'ERROR: pset%dust_tesc<5.5, you need to set dust_tesc '//&
          'to a value >=5.5; currently it is: ',pset%dust_tesc
     STOP
  ENDIF

  IF (pset%sfh.EQ.2.OR.pset%sfh.EQ.3) THEN

     istep = ntime
     imin  = locate(powitime(1:nt*ntime),sfh_tab(1,1))
     imin  = imax - INT((imax-imin)/istep)*istep
     !istep = INT((imax-imin)/nt)
     IF (imin.EQ.1) imin = istep
     !NB: this switch
     IF (pset%tage.EQ.-99.) imin=imax

  ELSE IF (pset%sfh.EQ.1) THEN
     
     istep = ntime
     !set a lower-limit to the integration
     IF (pset%tage.LE.0) THEN
        imin = imax - INT((imax-1)/ntime)*ntime
        IF (imin.EQ.1) imin = istep
     ELSE
        imin = imax
     ENDIF

  ENDIF

  !linearly interpolate the SSPs onto the time grid
  IF (pset%sfh.NE.0.AND.nzin.EQ.1) THEN
     DO j=1,nt*ntime
        !simple constant extrapolation to earlier times
        !this is not great, but its better than doing nothing
        IF (itime(j).LT.6.6) THEN
           ispec(j,:) = spec_ssp(1,1,:)
           lbol_ssp2(j) = lbol_ssp(1,1)
           mass_ssp2(j) = mass_ssp(1,1)
        ELSE
           klo = MAX(MIN(locate(time,itime(j)),nt-1),1)
           ispec(j,:) = spec_ssp(1,klo,:) + (itime(j)-time(klo))*&
                (spec_ssp(1,klo+1,:)-spec_ssp(1,klo,:))/0.05
           lbol_ssp2(j) = lbol_ssp(1,klo) + &
                (itime(j)-time(klo))*(lbol_ssp(1,klo+1)-lbol_ssp(1,klo))/0.05
           mass_ssp2(j) = mass_ssp(1,klo) + &
                (itime(j)-time(klo))*(mass_ssp(1,klo+1)-mass_ssp(1,klo))/0.05
        ENDIF
     ENDDO
  ENDIF

  !the isochrones don't go past 10**10.15 yrs, so warn the user
  !that this will be an extrapolation
  IF (maxtime.GT.1.42E10) THEN
     WRITE(*,*) 'COMPSP WARNING: Tmax>10^10.15 yrs -'//&
          ' linear extrapolation beyond this point for log(Tmax)=:',&
          LOG10(maxtime)
  ENDIF

  !warn the user about an out-of-bounds burst component
  IF (pset%tburst*1E9.GT.maxtime.AND.pset%fburst.GT.0.0.AND.&
       pset%sfh.EQ.1) THEN
     WRITE(*,*) 'COMPSP WARNING: burst time > age of system....'//&
          ' the burst component will NOT be added'
  ENDIF

  !set limits on the parameters tau and const
  IF (pset%sfh.EQ.1) THEN
     IF (pset%tau.LE.0.01.AND.pset%tau.GE.0.0) THEN
        IF (verbose.EQ.1) THEN
           WRITE(*,*) 'COMPSP WARNING: tau <0.01, setting tau=0.01'
        ENDIF
        tau = 0.01
     ELSE IF (pset%tau.GE.1E2) THEN
        IF (verbose.EQ.1) THEN
           WRITE(*,*) 'COMPSP WARNING: tau >1E2, setting tau=1E2'
        ENDIF
        tau = 1E2
     ELSE
        tau = pset%tau
     ENDIF

     const = pset%const     
     IF (const.LE.1E-5.AND.const.NE.0.0) THEN
        IF (verbose.EQ.1) THEN
           WRITE(*,*) 'COMPSP WARNING: const <1E-5, setting const=0.0'
        ENDIF
        const = 0.0
     ENDIF
     IF (const.GT.1.0) THEN
        WRITE(*,*) 'COMPSP ERROR: const >1.0, returning...'
        STOP
     ENDIF
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
     !write a header
     IF (pset%sfh.NE.2) THEN
        WRITE(10,'("#   Log(Z/Zsol): ",F6.3)') LOG10(zlegend(pset%zmet)/0.0190)
     ELSE
        WRITE(10,'("#   Log(Z/Zsol): tabulated")')
     ENDIF
     WRITE(10,'("#   Fraction of blue HB stars: ",F6.3)') pset%fbhb
     WRITE(10,'("#   Ratio of BS to HB stars  : ",F6.3)') pset%sbss
     WRITE(10,'("#   Shift to TP-AGB [log(Teff),log(Lbol)]: ",F5.2,1x,F5.2)') &
          pset%delt, pset%dell
     IF (imf_type.EQ.2) THEN
        WRITE(10,'("#   IMF: ",I1,", slopes= ",3F4.1)') &
             imf_type,pset%imf1,pset%imf2,pset%imf3
     ELSE IF (imf_type.EQ.3) THEN
        WRITE(10,'("#   IMF: ",I1,", cut-off= ",F4.2)') imf_type,pset%vdmc
     ELSE
        WRITE(10,'("#   IMF: ",I1)') imf_type
     ENDIF
  ENDIF
  
  !open output file for spectra
  IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) THEN
     
     OPEN(20,FILE=TRIM(SPS_HOME)//'/OUTPUTS/'//TRIM(outfile)//'.spec',&
          STATUS='REPLACE')
     
     !write a header
     IF (pset%sfh.NE.2) THEN
        WRITE(20,'("#   Log(Z/Zsol): ",F6.3)') LOG10(zlegend(pset%zmet)/0.0190)
     ELSE
        WRITE(20,'("#   Log(Z/Zsol): tabulated")')
     ENDIF
     WRITE(20,'("#   Fraction of blue HB stars: ",F6.3)') pset%fbhb
     WRITE(20,'("#   Ratio of BS to HB stars  : ",F6.3)') pset%sbss
     WRITE(20,'("#   Shift to TP-AGB [log(Teff),log(Lbol)]: ",F5.2,1x,F5.2)') &
          pset%delt, pset%dell
     IF (imf_type.EQ.2) THEN
        WRITE(20,'("#   IMF: ",I1,", slopes= ",3F4.1)') &
             imf_type,pset%imf1,pset%imf2,pset%imf3
     ELSE IF (imf_type.EQ.3) THEN
        WRITE(20,'("#   IMF: ",I1,", cut-off= ",F4.2)') imf_type,pset%vdmc
     ELSE
        WRITE(20,'("#   IMF: ",I1)') imf_type
     ENDIF
  ENDIF

  IF (pset%sfh.EQ.0) THEN
     
     IF (verbose.NE.0) WRITE(*,*) '  Processing SSP'
     IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) THEN
        WRITE(10,'("#   Processing SSP")')
        WRITE(10,'("#")') 
        WRITE(10,'("#   log(age) log(mass) Log(lbol) '//&
             'log(SFR) mags (see FILTER_LIST)")') 
     ENDIF
     IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) THEN
        WRITE(20,'("#   Processing SSP")')
        WRITE(20,'("#")') 
        WRITE(20,'("#   log(age) log(mass) Log(lbol) '//&
             'log(SFR) spectra")') 
        WRITE(20,'(I3)') nt
     ENDIF

  ELSE

    IF (pset%sfh.EQ.2.OR.pset%sfh.EQ.3) THEN
        IF (verbose.EQ.1) &
             WRITE(*,'("   SFH: tabulated input, dust=(",F6.2,","F6.2,")")'),&
             pset%dust1,pset%dust2
        IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) &
             WRITE(10,'("#   SFH: tabulated input, dust=(",F6.2,","F6.2,")")'),&
             pset%dust1,pset%dust2
        IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) &
             WRITE(20,'("#   SFH: tabulated input, dust=(",F6.2,","F6.2,")")'),&
             pset%dust1,pset%dust2
     ELSE
        
        IF (pset%tage.GT.0) THEN
           writeage = pset%tage
        ELSE
           writeage = tuniv
        ENDIF
        IF (verbose.EQ.1) &
             WRITE(*,'("   SFH: Tage=",F6.2," Gyr, log(tau)= ",F6.3,'//&
             '" Gyr, const= ",F5.3,", fb= ",F5.3,", tb= ",F6.2,'//&
             '" Gyr, dust=(",F6.2,","F6.2,")",", log(Z/Zsol)=",F6.3)')&
             writeage,LOG10(tau),const,pset%fburst,pset%tburst,&
             pset%dust1,pset%dust2,LOG10(zlegend(pset%zmet)/0.0190)
        IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) &
             WRITE(10,'("#   SFH: Tage=",F6.2," Gyr, log(tau)= ",F6.3,'//&
             '" Gyr, const= ",F5.3,", fb= ",F5.3,", tb= ",F6.2,'//&
             '" Gyr, dust=(",F6.2,","F6.2,")")')&
             writeage,LOG10(tau),const,pset%fburst,pset%tburst,&
             pset%dust1,pset%dust2
        IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) &
             WRITE(20,'("#   SFH: Tage=",F6.2," Gyr, log(tau)= ",F6.3,'//&
             '" Gyr, const= ",F5.3,", fb= ",F5.3,", tb= ",F6.2,'//&
             '" Gyr, dust=(",F6.2,","F6.2,")")')&
             writeage,LOG10(tau),const,pset%fburst,pset%tburst,&
             pset%dust1,pset%dust2
     ENDIF

     IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) THEN 
        WRITE(10,'("#")') 
        WRITE(10,'("#   log(age) log(mass) Log(lbol) '//&
             'log(SFR) mags (see FILTER_LIST)")') 
     ENDIF

     IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) THEN
        WRITE(20,'("#")') 
        WRITE(20,'("#   log(age) log(mass) Log(lbol) '//&
             'log(SFR) spectra")') 
          WRITE(20,'(I3)') (imax-imin)/istep+1
       ENDIF

     IF (pset%dust_clumps.NE.-99.AND.verbose.EQ.1) &
          WRITE(*,'("Adding clumpy dust with sigma=",F6.2)') pset%dust_clumps

  ENDIF

  !-----------------------------------------------------!
  !---------Generate composite spectra and mags---------!
  !-----------------------------------------------------!

  !just process the SSP
  IF (pset%sfh.EQ.0) THEN

     !calculate mags at each time-step
     DO i=1,nt

        csp1 = 0.0
        csp2 = 0.0
        IF (time(i).LT.pset%dust_tesc) THEN
           csp1 = spec_ssp(nzin,i,:)
        ELSE
           csp2 = spec_ssp(nzin,i,:)
        ENDIF

        CALL ADD_DUST(pset,csp1,csp2,csp_spec)

        IF (redshift_colors.EQ.0) THEN
           CALL GETMAGS(pset%zred,csp_spec,mags)
        ELSE
           zred = splint(zagespl(:,2),zagespl(:,1),&
                zagespl(:,3),10**time(i)/1E9)
           CALL GETMAGS(zred,csp_spec,mags)
        ENDIF

        !dump info into output structure
        ocompsp(i)%age      = time(i)
        ocompsp(i)%mass_csp = mass_ssp(nzin,i)
        ocompsp(i)%lbol_csp = lbol_ssp(nzin,i)
        ocompsp(i)%sfr      = 0.0
        ocompsp(i)%mags     = mags
        ocompsp(i)%spec     = csp_spec

        !write to mags file
        IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) &
             WRITE(10,fmt) time(i),LOG10(mass_ssp(nzin,i)),&
             lbol_ssp(nzin,i),-99.,mags

        !write to spectra file
        IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) THEN
           WRITE(20,'(4(F8.4,1x))') time(i),&
                LOG10(mass_ssp(nzin,i)),lbol_ssp(nzin,i),-99.0
           WRITE(20,*) csp_spec
        ENDIF

     ENDDO

  ELSE 

     !tsfr is SFR as a fcn of age of the Universe
     !tzhist is metallicity history as a fcn of age of Universe
     IF (pset%sfh.EQ.2.OR.pset%sfh.EQ.3) THEN
        DO j=1,nt*ntime
           jlo    = MAX(MIN(locate(LOG10(sfh_tab(1,1:ntabsfh)),&
                itime(j)),ntabsfh-1),1)
           tsfr(j) = sfh_tab(2,jlo) + (10**itime(j)-(sfh_tab(1,jlo)))*&
                (sfh_tab(2,jlo+1)-sfh_tab(2,jlo))/&
                (sfh_tab(1,jlo+1)-sfh_tab(1,jlo))
           tzhist(j) = sfh_tab(3,jlo) + (10**itime(j)-(sfh_tab(1,jlo)))*&
                (sfh_tab(3,jlo+1)-sfh_tab(3,jlo))/&
                (sfh_tab(1,jlo+1)-sfh_tab(1,jlo))
        ENDDO
     ELSE
        tsfr  = EXP(-powitime/tau/1E9 )/tau/1E9 / &
             (1-EXP(-maxtime/1E9/tau))*(1-const) &
             + const/maxtime
     ENDIF

     !actually generate the composite spectra
     !sfr and zhist are functions of the age of the stellar system
     DO i=imin,imax,istep

        IF (pset%sfh.EQ.2.OR.pset%sfh.EQ.3) THEN

           DO j=1,i
              jlo      = MAX(MIN(locate(powitime(i)-powitime(1:i),&
                   powitime(j)),nt*ntime-1),1)
              sfr(j)   = tsfr(jlo) + (powitime(j)-(powitime(i)-powitime(jlo)))&
                   * (tsfr(jlo+1)-tsfr(jlo))/(powitime(jlo+1)-powitime(jlo))
              zhist(j) = tzhist(jlo) + (powitime(j)-(powitime(i)-powitime(jlo)))&
                   * (tzhist(jlo+1)-tzhist(jlo))/(powitime(jlo+1)-powitime(jlo))
              IF (sfr(j).LT.0.0) sfr(j) = 0.0

              jlo = MAX(MIN(locate(time,itime(j)),nt-1),1)
              klo = MAX(MIN(locate(zlegend,zhist(j)),nz-1),1)
              t   = (zhist(j)-zlegend(klo)) / &
                   (zlegend(klo+1)-zlegend(klo))
              u   = (itime(j)-time(jlo)) / 0.05

              !bilinear interpolation over zhist and time
              ispec(j,:) = (1-t)*(1-u)*spec_ssp(klo,jlo,:) + &
                   t*(1-u)*spec_ssp(klo+1,jlo,:) + &
                   t*u*spec_ssp(klo+1,jlo+1,:) + &
                   (1-t)*u*spec_ssp(klo,jlo+1,:)
              lbol_ssp2(j) = (1-t)*(1-u)*lbol_ssp(klo,jlo) + &
                   t*(1-u)*lbol_ssp(klo+1,jlo) + &
                   t*u*lbol_ssp(klo+1,jlo+1) + &
                   (1-t)*u*lbol_ssp(klo,jlo+1)
              mass_ssp2(j) = (1-t)*(1-u)*mass_ssp(klo,jlo) + &
                   t*(1-u)*mass_ssp(klo+1,jlo) + &
                   t*u*mass_ssp(klo+1,jlo+1) + &
                   (1-t)*u*mass_ssp(klo,jlo+1)
           ENDDO

        ELSE

           sfr = 0.0
           sfr(1:i) = EXP(-(powitime(i)-powitime(1:i))/tau/1E9) / &
                tau/1E9/(1-EXP(-maxtime/1E9/tau)) * (1-const) &
                + const/maxtime
  
           !Add an instantaneous burst
           IF ((powitime(i)-pset%tburst*1E9).GT.0.0.AND.&
                pset%fburst.GT.0.0.AND.pset%sfh.EQ.1) THEN
              
              delt_burst = powitime(i)-pset%tburst*1E9
              klo        = MAX(MIN(locate(itime(1:nt*ntime),&
                   LOG10(delt_burst)),nt*ntime-1),1)
              dt2_burst  = (LOG10(delt_burst)-itime(klo))/&
                   (itime(klo+1)-itime(klo))
              spec_burst = ispec(klo,:)+&
                   dt2_burst*(ispec(klo+1,:)-ispec(klo,:))
              mass_burst = mass_ssp2(klo)+&
                   dt2_burst*(mass_ssp2(klo+1)-mass_ssp2(klo))
              lbol_burst = lbol_ssp2(klo)+&
                   dt2_burst*(lbol_ssp2(klo+1)-lbol_ssp2(klo))

           ENDIF

        ENDIF

        !note that we're intergrating over the ages
        !of the SSPs, not cosmic time
        CALL INTSPEC(pset,i,wtesc,dlogt,ispec,sfr,&
             csp_spec,mass_ssp2,lbol_ssp2,&
             mass_csp,lbol_csp,powitime,&
             spec_burst,mass_burst,lbol_burst,delt_burst)

        !calculate mags at each timestep
        IF (redshift_colors.EQ.0) THEN
           CALL GETMAGS(pset%zred,csp_spec,mags)
        ELSE
           zred = splint(zagespl(:,2),zagespl(:,1),&
                zagespl(:,3),10**itime(i)/1E9)
           zred = MAX(zred,0.0)
           CALL GETMAGS(zred,csp_spec,mags)
        ENDIF

        !dump info into output structure
        ocompsp((i-imin)/istep+1)%age      = itime(i)
        ocompsp((i-imin)/istep+1)%mass_csp = mass_csp
        ocompsp((i-imin)/istep+1)%lbol_csp = lbol_csp
        ocompsp((i-imin)/istep+1)%sfr      = tsfr(i)
        ocompsp((i-imin)/istep+1)%mags     = mags
        ocompsp((i-imin)/istep+1)%spec     = csp_spec

        IF (tsfr(i).EQ.0) THEN
           psfr = -99.
        ELSE
           psfr = LOG10(tsfr(i))
        ENDIF

        !write to mags file
        IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) &
             WRITE(10,fmt) itime(i),LOG10(mass_csp),&
             lbol_csp,psfr,mags

        !write to spectra file
        IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) THEN
           WRITE(20,'(4(F8.4,1x))') itime(i),&
                LOG10(mass_csp),lbol_csp,psfr
           WRITE(20,*) csp_spec
        ENDIF
        
     ENDDO

  ENDIF

  IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) CLOSE(10)
  IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) CLOSE(20)

END SUBROUTINE COMPSP 
