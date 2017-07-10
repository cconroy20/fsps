SUBROUTINE COMPSP(write_compsp, nzin, outfile,&
                  mass_ssp, lbol_ssp, tspec_ssp,&
                  pset, ocompsp)
  !
  !
  !N.B. variables not otherwise defined come from sps_vars.f90
  use sps_vars
  use sps_utils, only: write_isochrone, add_nebular, setup_tabular_sfh, &
                       csp_gen, sfhinfo, linterp, agn_dust, &
                       smoothspec, igm_absorb, getindx, getmags
                       
  implicit none

  INTEGER, INTENT(in) :: write_compsp,nzin
  CHARACTER(100), INTENT(in) :: outfile
  REAL(SP), INTENT(in), DIMENSION(ntfull, nzin) :: lbol_ssp,mass_ssp
  REAL(SP), INTENT(in), DIMENSION(nspec, ntfull, nzin) :: tspec_ssp
  TYPE(PARAMS), intent(in) :: pset

  TYPE(COMPSPOUT), INTENT(inout), DIMENSION(ntfull) :: ocompsp

  REAL(SP), DIMENSION(nspec, ntfull, nzin)   :: spec_ssp
  REAL(SP), DIMENSION(nemline, ntfull, nzin) :: emlin_ssp
  REAL(SP), DIMENSION(nemline) :: emlin_csp
  REAL(SP), dimension(ntfull) :: mdust_ssp
  REAL(SP) :: lbol_csp, mass_csp, mdust_csp
  REAL(SP) :: age, mdust, mass_frac, tsfr, zred, frac_linear, maxtime
  REAL(SP), DIMENSION(nspec) :: csp1, csp2, spec_dusty, spec_csp
  REAL(SP), DIMENSION(nbands)  :: mags
  REAL(SP), DIMENSION(nindx)   :: indx
  INTEGER :: i, nage

  ! ------ Various checks and setup ------

  IF (check_sps_setup.EQ.0) THEN
     WRITE(*,*) 'COMPSP ERROR: '//&
          'SPS_SETUP must be run before calling COMPSP. '
     STOP
  ENDIF

  IF ((nzin.GT.1).and.(pset%sfh.ne.2).and.(pset%sfh.ne.3)) THEN
     WRITE(*,*) 'COMPSP ERROR: '//&
          'nzin > 1 no longer supported for non-tabular SFH.'
     STOP
  ENDIF

  call setup_tabular_sfh(pset, nzin)

  ! Make sure various variables are set correctly
  IF (pset%tage.GT.tiny_number) THEN
     maxtime = pset%tage * 1e9
  else
     maxtime = 10**time_full(ntfull)
  endif
  
  CALL COMPSP_WARNING(maxtime, pset, nzin, write_compsp)

  ! Setup output files
  if ((pset%tage.gt.0).or.&
      ((pset%tage.eq.-99).and.((pset%sfh.eq.2).or.(pset%sfh.eq.3)))) then
     nage = 1
  else
     nage = ntfull
  endif
  IF (write_compsp.GT.0) &
       CALL COMPSP_SETUP_OUTPUT(write_compsp, pset, outfile, 1, nage)

  ! Isochrone case just writes the CMDs and exits
  IF (write_compsp.EQ.5) THEN
     CALL WRITE_ISOCHRONE(outfile, pset)
     RETURN
  ENDIF

  ! ------ Prepare SSPs ------
  ! Only doing nebular emission at the moment, dust is added in csp_gen.
  ! We should probably only do this for ages up to tage, if it is set.
  ! Also we will operate on copies of the spectra

  spec_ssp = tspec_ssp

  ! Add nebular emission
  if (add_neb_emission.EQ.1) then
     if (nzin.GT.1) then
        WRITE(*,*) 'COMPSP ERROR: cannot handle both nebular '//&
             'emission and mult-metallicity SSPs in compsp'
        STOP
     endif
     call add_nebular(pset, tspec_ssp(:,:,1), spec_ssp(:,:,1), emlin_ssp(:,:,1))
  else
     emlin_ssp = 0.
  endif


  ! --- Get CSP spectra -------
  
  ! Loop over output ages.
  do i=1,ntfull
     ! ------
     ! First decide what mode we are in.
     if (pset%tage.gt.0) then
        ! A specific age was asked for, so we will only compute one spectrum at
        ! that age.
        age = pset%tage
     else if ((pset%tage.eq.-99).and.((pset%sfh.eq.2).or.(pset%sfh.eq.3))) then
        ! Special switch to just do the last time in the tabular file
        age = maxval(sfh_tab(1, 1:ntabsfh)) / 1E9
     else
        ! Otherwise we will calculate composite spectra for every SSP age.
        age = 10**(time_full(i)-9.)
     endif

     ! -----
     ! Get the spectrum for this age.  Note this is always normalized to one
     ! solar mass formed, so we actually need to renormalize if computing all
     ! ages, which is done using info from `sfhinfo`
     call csp_gen(mass_ssp, lbol_ssp, spec_ssp, &
          pset, age, nzin, mass_csp, lbol_csp, spec_csp,&
          mdust_csp,emlin_ssp,emlin_csp)
     
     call sfhinfo(pset, age, mass_frac, tsfr, frac_linear)
     if (pset%tage.le.0) then   
        mass_csp  = mass_csp * mass_frac
        lbol_csp  = log10(10**lbol_csp * mass_frac)
        spec_csp  = spec_csp * mass_frac
        mdust_csp = mdust_csp * mass_frac
        emlin_csp = emlin_csp * mass_frac
     else
        ! Renormalize the SFR to be appropriate for one solar mass formed.
        tsfr = tsfr / mass_frac
        mass_frac = 1.0
     endif

     ! -------
     ! Now do a bunch of stuff with the spectrum
     ! Smooth the spectrum
     if (pset%sigma_smooth.GT.0.0) then
        call smoothspec(spec_lambda, spec_csp, pset%sigma_smooth,&
                        pset%min_wave_smooth, pset%max_wave_smooth)
     endif
     ! Add IGM absorption
     if (add_igm_absorption.EQ.1.AND.pset%zred.GT.tiny_number) then
        spec_csp = igm_absorb(spec_lambda,spec_csp, pset%zred,&
                              pset%igm_factor)
     endif
     !add AGN dust
     IF (add_agn_dust.EQ.1.AND.pset%fagn.GT.tiny_number) THEN
        spec_csp = agn_dust(spec_lambda, spec_csp, pset, lbol_csp)
     ENDIF
     ! Compute spectral indices
     if (write_compsp.EQ.4) then
        call getindx(spec_lambda, spec_csp, indx)
     else
        indx = 0.0
     endif
     ! Compute mags
     if (redshift_colors.EQ.0) then
        call getmags(pset%zred, spec_csp, mags, pset%mag_compute)
     else
        ! here we compute the redshift at the corresponding age
        zred = min(max(linterp(cosmospl(:,2), cosmospl(:,1), age),&
                       0.0), 20.0)
        write(33,*) zred
        call getmags(zred, spec_csp, mags, pset%mag_compute)
     endif

     ! ---------
     ! Store the spectrum and write....
     call save_compsp(write_compsp, ocompsp(i), log10(age)+9,&
          mass_csp, lbol_csp, tsfr, mags, spec_csp, mdust_csp, mass_frac,&
          indx,emlin_csp)

     ! Terminate the loop if a single specific tage was requested
     if (nage.eq.1) then
        exit
     endif

  enddo

  if (write_compsp.EQ.1.OR.write_compsp.EQ.3) CLOSE(10)
  if (write_compsp.EQ.2.OR.write_compsp.EQ.3) CLOSE(20)

end subroutine compsp

SUBROUTINE COMPSP_WARNING(maxtime,pset,nzin,write_compsp)

  !check that variables are properly set

  USE sps_vars
  IMPLICIT NONE
  INTEGER, INTENT(in) :: nzin, write_compsp
  REAL(SP), INTENT(in) :: maxtime
  TYPE(PARAMS), INTENT(in) :: pset

  !-----------------------------------------------------!

  IF (maxtime.LE.0.AND.pset%sfh.NE.0) THEN
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
          ' the burst component will NOT be added.'
  ENDIF

  IF (pset%sf_start.LT.0.0) THEN
     WRITE(*,*) 'COMPSP ERROR: sf_start<0.  stopping...'
     STOP
  ENDIF

  IF (pset%sf_start*1E9.GT.maxtime) THEN
     WRITE(*,*) 'COMPSP ERROR: sf_start>maxtime  stopping...'
     STOP
  ENDIF

  IF ((pset%sf_trunc.LT.pset%sf_start).AND.(pset%sf_trunc.GT.tiny_number)) THEN
     WRITE(*,*) 'COMPSP WARNING: sf_trunc<sf_start....'//&
          ' sf_trunc will be ignored.'
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

     IF ((pset%const + pset%fburst).GT.1.0) THEN
        WRITE(*,*) 'COMPSP ERROR: const + fburst > 1', pset%const + pset%fburst
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

  if ((pset%dust1.gt.tiny_number).and.(compute_light_ages.eq.1)) then
     WRITE(*,*) 'COMPSP WARNING: compute_light_ages does not take into'//&
          ' account age-dependent dust (dust1 > 0)'
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
     lbol,sfr,mags,spec,mdust,mformed,indx,emlines)

  !routine to print and save outputs

  USE sps_vars
  IMPLICIT NONE
  INTEGER, INTENT(in) :: write_compsp
  REAL(SP), INTENT(in)    :: time,mass,lbol,sfr,mdust,mformed
  REAL(SP), DIMENSION(nspec), INTENT(in)    :: spec
  REAL(SP), DIMENSION(nbands), INTENT(in)   :: mags
  REAL(SP), DIMENSION(nindx), INTENT(in)    :: indx
  REAL(SP), DIMENSION(nemline), INTENT(in)  :: emlines
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
  cspo%spec     = MAX(spec,tiny_number)
  cspo%mdust    = mdust
  cspo%mformed  = mformed
  cspo%indx     = indx
  cspo%emlines  = emlines

  !write to mags file
  IF (write_compsp.EQ.1.OR.write_compsp.EQ.3) &
       WRITE(10,fmt) time,LOG10(mass+tiny_number),&
       lbol,LOG10(sfr+tiny_number),mags
 
  !write to spectra file
  IF (write_compsp.EQ.2.OR.write_compsp.EQ.3) THEN
     WRITE(20,'(4(F8.4,1x))') time,&
          LOG10(mass+tiny_number),lbol,LOG10(sfr+tiny_number)
     WRITE(20,'(50000(E14.6))') MAX(spec,tiny_number)
  ENDIF

  !write to indx file
  IF (write_compsp.EQ.4) &
       WRITE(30,'(F8.4,99(F7.3,1x))') time,indx
  
END SUBROUTINE SAVE_COMPSP
