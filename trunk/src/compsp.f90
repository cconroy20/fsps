!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Routine to compute magnitudes and spectra for a composite !    
!  stellar population.  Returns a file with the following    !
!  info: age, mass, Lbol, mags in various filters.  SFR      !
!  units are Msun/yr.                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL FUNCTION INTSFR(sfh,tau,const,maxtime,sfstart,t1,t2)

  !routine to integrate an analytic SFH from t1 to t2

  USE nrtype; USE sps_vars
  USE nr, ONLY : locate
  IMPLICIT NONE

  INTEGER, INTENT(in)  :: sfh
  REAL(SP), INTENT(in) :: t1,t2,tau,const,maxtime,sfstart
  INTEGER  :: lo,hi
  REAL(SP) :: s1,s2,dt

  IF (sfh.EQ.1) THEN

     intsfr = (EXP(-t2/tau/1E9) - EXP(-t1/tau/1E9))/&
          (1-EXP(-(maxtime-sfstart)/1E9/tau)) 
     intsfr = intsfr* (1-const) + const*(t1-t2)/(maxtime-sfstart)

  ELSE IF (sfh.EQ.4) THEN

     intsfr = (EXP(-t2/tau/1E9)*(1+t2/tau/1E9) - &
          EXP(-t1/tau/1E9)*(1+t1/tau/1E9)) / &
          (1-EXP(-(maxtime-sfstart)/1E9/tau)*((maxtime-sfstart)/1E9/tau+1))  
     intsfr = intsfr* (1-const) + const*(t1-t2)/(maxtime-sfstart)
   
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
        intsfr = SUM((sfh_tab(1,lo+1:hi)-sfh_tab(1,lo:hi-1))/2.*&
             (sfh_tab(2,lo:hi-1)+sfh_tab(2,lo+1:hi)))
        intsfr = intsfr + (t1-sfh_tab(1,hi))*0.5*(s1+sfh_tab(2,hi))
        intsfr = intsfr + (sfh_tab(1,lo)-t2)*0.5*(s2+sfh_tab(2,lo))
     ENDIF

  ENDIF
  
  
END FUNCTION INTSFR

!------------------------------------------------------------!
!------------------------------------------------------------!

SUBROUTINE INTSPEC(pset,nti,spec_ssp,csp,mass_ssp,lbol_ssp,&
     mass,lbol,specb,massb,lbolb,deltb,sfstart,tau,const,maxtime)

  !routine to perform integration of SSP over a SFH,
  !including attenuation by dust.

  USE nrtype; USE sps_vars
  USE nr, ONLY : locate
  USE sps_utils, ONLY : add_dust
  IMPLICIT NONE
  
  INTEGER :: i,imax,indsf,wtesc
  REAL(SP) :: intsfr,t1,t2
  REAL(SP), DIMENSION(ntfull) :: isfr,time
  INTEGER,  INTENT(in)    :: nti
  REAL(SP), INTENT(in)    :: massb,lbolb,deltb,sfstart,tau,const,maxtime
  REAL(SP), INTENT(inout) :: mass, lbol
  REAL(SP), INTENT(in), DIMENSION(ntfull) :: mass_ssp,lbol_ssp
  REAL(SP), INTENT(in), DIMENSION(ntfull,nspec) :: spec_ssp
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
        
        IF (pset%sfh.EQ.5) THEN
           IF ((t1+sfstart).LT.pset%sf_trunc*1E9) THEN
              isfr(i) = (EXP(-t2/tau/1E9)*(1+t2/tau/1E9) - &
                   EXP(-t1/tau/1E9)*(1+t1/tau/1E9)) * (tau*1E9)**2
           ELSE
              isfr(i) = (pset%sf_trunc*1E9-sfstart)*&
                   EXP(-(pset%sf_trunc-sfstart/1E9)/tau)*(t1-t2)
              isfr(i) = isfr(i) + TAN(pset%sf_theta)*&
                   (0.5*(t1**2-t2**2)-(t1-t2)*(pset%sf_trunc*1E9-sfstart))
           ENDIF
           isfr(i) = MAX(isfr(i),0.0) / 1E10
        ELSE
           isfr(i) = intsfr(pset%sfh,tau,const,maxtime,sfstart,t1,t2)
        ENDIF

     ENDDO

     IF (indsf.EQ.1) THEN
        csp1 = isfr(1)*spec_ssp(1,:)
     ELSE
        !t<tesc
        DO i=1,imax-1
           csp1 = csp1 + isfr(i)*0.5*(spec_ssp(i+1,:)+spec_ssp(i,:))
        ENDDO
        !t>tesc
        IF (nti.GT.wtesc) THEN
           DO i=wtesc,indsf-1
              csp2 = csp2 + isfr(i)*0.5*(spec_ssp(i+1,:)+spec_ssp(i,:))
           ENDDO
           csp2 = csp2 + isfr(indsf)*spec_ssp(indsf,:)
        ELSE
           csp1 = csp1 + isfr(indsf)*spec_ssp(indsf,:)
        ENDIF
     ENDIF

     !compute weighted mass and lbol
     IF (indsf.EQ.1) THEN
        mass = isfr(1)*mass_ssp(1)
        lbol = LOG10( isfr(1)*10**lbol_ssp(1) )
     ELSE
        mass = SUM(isfr(1:indsf-1)/2.*&
             (mass_ssp(1:indsf-1)+mass_ssp(2:indsf)))
        lbol = SUM(isfr(1:indsf-1)/2.*&
             (10**lbol_ssp(1:indsf-1)+10**lbol_ssp(2:indsf)))
        mass = mass + isfr(indsf)*mass_ssp(indsf)
        lbol = lbol + isfr(indsf)*10**lbol_ssp(indsf) 
        lbol = LOG10(lbol)
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
     CALL ADD_DUST(pset,csp1,csp2,csp)

  ENDIF

END SUBROUTINE INTSPEC

!------------------------------------------------------------!
!---------------------Main Routine---------------------------!
!------------------------------------------------------------!

SUBROUTINE COMPSP(write_compsp,nzin,outfile,mass_ssp,&
     lbol_ssp,spec_ssp,pset,ocompsp)

  !for sfh=1 or 4:
  !If tage >0  -> run only one integration to t=tage
  !If tage<=0  -> produce outputs from tmin<t<tmax

  USE sps_vars; USE nrtype
  USE sps_utils, ONLY : getmags, add_dust
  USE nrutil, ONLY : assert_eq 
  USE nr, ONLY : locate, splint, spline
  IMPLICIT NONE
 
  !write_compsp = (1->write mags), (2->write spectra), (3->write both) 
  INTEGER, INTENT(in) :: write_compsp,nzin
  REAL(SP), INTENT(in), DIMENSION(nzin,ntfull) :: lbol_ssp,mass_ssp
  REAL(SP), INTENT(in), DIMENSION(nzin,ntfull,nspec) :: spec_ssp
  CHARACTER(100), INTENT(in) :: outfile

  INTEGER  :: i,j,n,k,stat,klo,jlo,ilo,imin,imax,indsf,indsft
  REAL(SP) :: tau,const,maxtime,writeage,psfr,sfstart,zhist,sftrunc
  REAL(SP) :: mass_csp,lbol_csp,dtb,dt,dz,zred=0.,t1,t2
  REAL(SP) :: mass_burst=0.0,lbol_burst=0.0,delt_burst=0.0
  REAL(SP), DIMENSION(nbands)       :: mags
  REAL(SP), DIMENSION(ntfull,nspec) :: ispec
  REAL(SP), DIMENSION(nspec)  :: spec_csp,spec_burst=0.0,csp1,csp2,spec1
  REAL(SP), DIMENSION(ntfull) :: imass,ilbol,powtime
  REAL(SP), DIMENSION(ntfull) :: sfr,tsfr,tzhist
  REAL(SP), DIMENSION(ntabmax) :: tlb

  !(TYPE objects defined in sps_vars.f90)
  TYPE(PARAMS), INTENT(in) :: pset
  TYPE(COMPSPOUT), INTENT(inout), DIMENSION(ntfull) :: ocompsp

  !-------------------------------------------------------------!
  !------------------------Basic Setup--------------------------!
  !-------------------------------------------------------------!

  dtb        = 0.0
  spec_burst = 0.0
  mass_burst = 0.0
  lbol_burst = 0.0
  sfstart    = 0.0

  IF (check_sps_setup.EQ.0) THEN
     WRITE(*,*) 'COMPSP ERROR: '//&
          'SPS_SETUP must be run before calling COMPSP. '
     STOP
  ENDIF
 
  powtime = 10**time_full
  maxtime = powtime(ntfull)
  imin    = 1
  imax    = ntfull

  !if a tabulated SFH is passed, read in sfh.dat file
  IF (pset%sfh.EQ.2) THEN

     OPEN(3,FILE=TRIM(SPS_HOME)//'/data/sfh.dat',ACTION='READ',STATUS='OLD')
     DO n=1,ntabmax
        READ(3,*,IOSTAT=stat) sfh_tab(1,n),sfh_tab(2,n),sfh_tab(3,n)
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

  ENDIF

  !set limits on the parameters tau and const
  IF (pset%sfh.EQ.1.OR.pset%sfh.EQ.4.OR.pset%sfh.EQ.5) THEN
     tau   = MIN(MAX(pset%tau,0.1),100.)
     const = pset%const
  ENDIF

  !make sure various variables are set correctly
  CALL COMPSP_WARNING(maxtime,pset,nzin,write_compsp)

  !tsfr and tzhist are functions of the age of Universe
  IF (pset%sfh.EQ.2.OR.pset%sfh.EQ.3) THEN

     !linearly interpolate the tabulated SFH to the internal time grid
     DO j=1,ntfull
        IF (time_full(j).GT.LOG10(sfh_tab(1,ntabsfh))) THEN
           tsfr(j)   = 0.0
           tzhist(j) = 0.02
        ELSE
           jlo    = MAX(MIN(locate(LOG10(sfh_tab(1,1:ntabsfh)),&
                time_full(j)),ntabsfh-1),1)
           dt = (powtime(j)-(sfh_tab(1,jlo))) / &
                (sfh_tab(1,jlo+1)-sfh_tab(1,jlo))
           tsfr(j)   = (1-dt)*sfh_tab(2,jlo)+dt*sfh_tab(2,jlo+1)
           tzhist(j) = (1-dt)*sfh_tab(3,jlo)+dt*sfh_tab(3,jlo+1)
        ENDIF
     ENDDO

  ELSE IF (pset%sfh.EQ.1.OR.pset%sfh.EQ.4.OR.pset%sfh.EQ.5) THEN

     tsfr  = 0.0
     IF (pset%sfh.EQ.1) &
          tsfr(indsf:)  = EXP(-(powtime(indsf:)-sfstart)/tau/1E9 )/&
          tau/1E9 / (1-EXP(-(maxtime-sfstart)/1E9/tau))
     IF (pset%sfh.EQ.4) &
          tsfr(indsf:)  = ((powtime(indsf:)-sfstart)/tau/1E9)*&
          EXP(-(powtime(indsf:)-sfstart)/tau/1E9 )/tau/1E9 / &
          (1-EXP(-(maxtime-sfstart)/1E9/tau)*((maxtime-sfstart)/1E9/tau+1))
     tsfr(indsf:) = tsfr(indsf:)*(1-const) + const/(maxtime-sfstart)

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
  !-------------------Setup output files------------------------!
  !-------------------------------------------------------------!

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
        IF (imax-imin.EQ.1) WRITE(20,'(I3)') 1
        IF (imax-imin.GT.1) WRITE(20,'(I3)') ntfull
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
        IF (pset%tage.LE.tiny_number) writeage = time_full(ntfull)
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
        IF (imax-imin.EQ.1) WRITE(20,'(I3)') 1
        IF (imax-imin.GT.1) WRITE(20,'(I3)') ntfull
       ENDIF
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
               dt  = MAX(MIN(dt,1.0),-1.0) !don't extrapolation
               zhist = (1-dt)*sfh_tab(3,jlo+1)+dt*sfh_tab(3,jlo)
               !interpolation over zhist
               klo = MAX(MIN(locate(zlegend,zhist),nz-1),1)
               dz  = (LOG10(zhist)-LOG10(zlegend(klo))) / &
                  (LOG10(zlegend(klo+1))-LOG10(zlegend(klo)))
               dz = MAX(MIN(dz,1.0),-1.0) !don't extrapolate
               ispec(j,:) = (1-dz)*spec_ssp(klo,j,:)+dz*spec_ssp(klo+1,j,:)
               ilbol(j)   = (1-dz)*lbol_ssp(klo,j)  +dz*lbol_ssp(klo+1,j)
               imass(j)   = (1-dz)*mass_ssp(klo,j)  +dz*mass_ssp(klo+1,j)
            ENDDO

         ELSE
            ispec = spec_ssp(1,:,:)
            ilbol = lbol_ssp(1,:)
            imass = mass_ssp(1,:)
         ENDIF

      !set up an instantaneous burst
      ELSE IF ((pset%sfh.EQ.1.OR.pset%sfh.EQ.4).AND.&
           pset%fburst.GT.tiny_number) THEN

         IF ((powtime(i)-pset%tburst*1E9).GT.tiny_number) THEN
            delt_burst = powtime(i)-pset%tburst*1E9
            klo = MAX(MIN(locate(time_full,LOG10(delt_burst)),ntfull-1),1)
            dtb = (LOG10(delt_burst)-time_full(klo))/&
                 (time_full(klo+1)-time_full(klo))
            spec_burst = (1-dtb)*spec_ssp(1,klo,:)+dtb*spec_ssp(1,klo+1,:)
            mass_burst = (1-dtb)*mass_ssp(1,klo)  +dtb*mass_ssp(1,klo+1)
            lbol_burst = (1-dtb)*lbol_ssp(1,klo)  +dtb*lbol_ssp(1,klo+1)
         ENDIF

      ENDIF

      !compute composite spectra, mass, lbol
      IF (pset%sfh.EQ.0) THEN
         csp1 = 0.0
         csp2 = 0.0
         IF (time_full(i).LT.pset%dust_tesc) THEN
            csp1 = spec_ssp(1,i,:)
         ELSE
            csp2 = spec_ssp(1,i,:)
         ENDIF
         !add dust and combine young and old csp
         CALL ADD_DUST(pset,csp1,csp2,spec_csp)
         mass_csp = mass_ssp(1,i)
         lbol_csp = lbol_ssp(1,i)
      ELSE IF (pset%sfh.EQ.1.OR.pset%sfh.EQ.4.OR.pset%sfh.EQ.5) THEN
         CALL INTSPEC(pset,i,spec_ssp,spec_csp,mass_ssp,lbol_ssp,&
              mass_csp,lbol_csp,spec_burst,mass_burst,&
              lbol_burst,delt_burst,sfstart,tau,const,maxtime)
      ELSE IF (pset%sfh.EQ.2.OR.pset%sfh.EQ.3) THEN
         CALL INTSPEC(pset,i,ispec,spec_csp,imass,ilbol,mass_csp,&
              lbol_csp,spec_burst,mass_burst,lbol_burst,&
              delt_burst,sfstart,tau,const,maxtime)
      ENDIF
      
      !calculate mags
      IF (redshift_colors.EQ.0) THEN
         CALL GETMAGS(pset%zred,spec_csp,mags)
      ELSE
         !here we compute the redshift at the corresponding age
         zred = MIN(MAX(splint(zagespl(:,2),zagespl(:,1),&
              zagespl(:,3),powtime(i)/1E9),0.0),20.0)
         CALL GETMAGS(zred,spec_csp,mags)
      ENDIF

      !only save results if computing all ages
      IF (imax-imin.GT.1) THEN
         CALL SAVE_COMPSP(write_compsp,ocompsp(i),time_full(i),&
              mass_csp,lbol_csp,tsfr(i),mags,spec_csp)
      ELSE
         !save results temporarily for later interpolation
         ocompsp(i)%mass_csp = mass_csp
         ocompsp(i)%lbol_csp = lbol_csp
         ocompsp(i)%mags     = mags
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
      spec_csp = (1-dt)*ocompsp(imin)%spec + &
           dt*ocompsp(imax)%spec
      CALL SAVE_COMPSP(write_compsp,ocompsp(1),&
           LOG10(maxtime),mass_csp,lbol_csp,0.0,mags,&
           spec_csp)
   ENDIF

   IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) CLOSE(10)
   IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) CLOSE(20)
   

   !formats
30 FORMAT('#   SFH: tabulated input, dust=(',F6.2,','F6.2,')')
31 FORMAT('#   log(age) log(mass) Log(lbol) log(SFR) spectra')
32 FORMAT('#   log(age) log(mass) Log(lbol) log(SFR) mags (see FILTER_LIST)')
33 FORMAT('#   SFH: Tage=',F6.2,' Gyr, log(tau/Gyr)= ',F6.3,&
        ', const= ',F5.3,', fb= ',F5.3,', tb= ',F6.2,&
        ' Gyr, dust=(',F6.2,','F6.2,')')

END SUBROUTINE COMPSP

!------------------------------------------------------------!
!------------------------------------------------------------!
 
SUBROUTINE COMPSP_WARNING(maxtime, pset, nzin, write_compsp)

  !check that variables are properly set

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

  !the isochrones don't go past 10**10.15 yrs, so warn the user
  !that this will be an extrapolation
  IF (maxtime.GT.10**10.20) THEN
     WRITE(*,*) 'COMPSP WARNING: Tmax>10^10.2 yrs -'//&
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

  !writes headers for the .mag and .spec files

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

!------------------------------------------------------------!
!------------------------------------------------------------!

SUBROUTINE SAVE_COMPSP(write_compsp,cspo,time,mass,&
     lbol,sfr,mags,spec)

  !routine to print and save outputs

  USE sps_vars; USE nrtype
  IMPLICIT NONE
  INTEGER, INTENT(in) :: write_compsp
  REAL, INTENT(in)    :: time,mass,lbol,sfr
  REAL, DIMENSION(nspec), INTENT(in)  :: spec
  REAL, DIMENSION(nbands), INTENT(in) :: mags
  TYPE(COMPSPOUT), INTENT(inout) :: cspo
  CHARACTER(34) :: fmt

  !-----------------------------------------------------!

  fmt = '(F7.4,1x,3(F8.4,1x),000(F7.3,1x))'
  WRITE(fmt(21:23),'(I3)') nbands

  !dump info into output structure
  cspo%age      = time
  cspo%mass_csp = mass
  cspo%lbol_csp = lbol
  cspo%sfr      = sfr
  cspo%mags     = mags
  cspo%spec     = spec
     
  !write to mags file
  IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) &
       WRITE(10,fmt) time,LOG10(mass+tiny_number),&
       lbol,LOG10(sfr+tiny_number),mags
  
  !write to spectra file
  IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) THEN
     WRITE(20,'(4(F8.4,1x))') time,&
          LOG10(mass+tiny_number),lbol,LOG10(sfr+tiny_number)
     WRITE(20,*) spec
  ENDIF

END SUBROUTINE SAVE_COMPSP