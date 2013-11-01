SUBROUTINE SPS_SETUP(zin)

  !read in isochrones and spectral libraries for all metallicities.
  !read in band-pass info and the spectrum for Vega.
  !Arrays are stored in a common block defined in sps_vars.f90

  !If zin=-1 then all metallicities are read in, otherwise only the
  !metallicity corresponding to zin in the look-up table zlegend.dat
  !is read.

  USE sps_vars; USE sps_utils
  IMPLICIT NONE
  INTEGER, INTENT(in) :: zin
  INTEGER :: stat=1,n,i,j,m,jj,k
  INTEGER :: n_isoc,z,zmin,zmax,nlam
  CHARACTER(1) :: char,sqpah
  CHARACTER(6) :: zstype
  REAL(SP) :: dumr1,d1,d2,logage,x,a, zero=0.0
  REAL(SP), DIMENSION(nspec) :: tspec
  REAL(SP), DIMENSION(1221) :: tvega_lam,tvega_spec,tsun_lam,tsun_spec
  REAL(SP), DIMENSION(50000) :: readlamb,readband
  REAL(SP), DIMENSION(25) :: wglam
  REAL(SP), DIMENSION(25,18,2,6) :: wgtmp
  REAL(SP), DIMENSION(10000) :: lambda_dagb, fluxin_dagb

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  IF (verbose.EQ.1) THEN 
     WRITE(*,*) 
     WRITE(*,*) '    Setting up SPS...'
  ENDIF

  IF (zin.GT.nz) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: zin GT nz', zin,nz
     STOP
  ENDIF

  !clean out all the common block arrays
  mini_isoc     = 0.
  mact_isoc     = 0.
  logl_isoc     = 0.
  logt_isoc     = 0.
  logg_isoc     = 0.
  ffco_isoc     = 0.
  phase_isoc    = 0.
  nmass_isoc    = 0
  timestep_isoc = 0.
  spec_lambda   = 0.
  speclib       = 0.
  basel_logg    = 0.
  basel_logt    = 0.
  agb_spec_o    = 0.
  agb_logt_o    = 0.
  agb_spec_c    = 0.
  agb_logt_c    = 0.
  n_isoc        = 0
  m             = 1

  !----------------------------------------------------------------!
  !--------------Confirm that variables are properly set-----------!
  !----------------------------------------------------------------!
 
  CALL getenv('SPS_HOME',SPS_HOME)
  IF (LEN_TRIM(SPS_HOME).EQ.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: spsdir environment variable not set!'
     STOP
  ENDIF

  IF (isoc_type.NE.'pdva'.AND.isoc_type.NE.'bsti'.AND.&
       isoc_type.NE.'mist') THEN
     WRITE(*,*) 'SPS_SETUP ERROR: isoc_type var set to invalid type: ',isoc_type
     STOP
  ENDIF

  IF (basel_str.NE.'pdva'.AND.basel_str.NE.'wlbc') THEN
     WRITE(*,*) 'SPS_SETUP ERROR: basel_str var set to invalid type: ',basel_str
     STOP
  ENDIF

  IF (spec_type.NE.'basel'.AND.spec_type.NE.'miles'.AND.spec_type.NE.'picks'&
       .AND.spec_type(1:5).NE.'hrlib'.AND.spec_type.NE.'rrlib') THEN
     WRITE(*,*) 'SPS_SETUP ERROR: spec_type var set to invalid type: ',spec_type
     STOP
  ENDIF

  IF (isoc_type.EQ.'pdva'.AND.spec_type.EQ.'basel'.AND.nz.NE.22) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: isoc_type="pdva", spec_type="basel", but nz NE 22!'
     STOP
  ENDIF

  IF (isoc_type.EQ.'bsti'.AND.spec_type.EQ.'basel'.AND.nz.NE.10) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: isoc_type="bsti", spec_type="basel", but nz NE 10!'
     STOP
  ENDIF

  IF (spec_type.EQ.'miles'.AND.nz.NE.5) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: spec_type="miles", but nz NE 5!'
     STOP
  ENDIF


  !----------------------------------------------------------------!
  !----------------Read in metallicity values----------------------!
  !----------------------------------------------------------------!
  
  !units are simply metal fraction by mass (e.g. Z=0.0190 for Zsun)
  IF (isoc_type.EQ.'pdva') THEN
     OPEN(90,FILE=TRIM(SPS_HOME)//'/ISOCHRONES/Padova/Padova2007/zlegend_'//&
          spec_type(1:5)//'.dat',STATUS='OLD',iostat=stat,ACTION='READ')
  ELSE IF (isoc_type.EQ.'bsti') THEN
     OPEN(90,FILE=TRIM(SPS_HOME)//'/ISOCHRONES/BaSTI/zlegend_'//&
          spec_type(1:5)//'.dat',STATUS='OLD',iostat=stat,ACTION='READ')
 ELSE IF (isoc_type.EQ.'mist') THEN
     OPEN(90,FILE=TRIM(SPS_HOME)//'/ISOCHRONES/MIST/zlegend_'//&
          spec_type(1:5)//'.dat',STATUS='OLD',iostat=stat,ACTION='READ')
  ENDIF
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: zlegend_*.dat '//&
          'cannot be opened'
     STOP 
  END IF
  DO z=1,nz
     READ(90,'(F6.4)') zlegend(z)
  ENDDO
  CLOSE(90)

  IF (zin.LE.0) THEN
     zmin = 1
     zmax = nz
  ELSE
     zmin = zin
     zmax = zin
  ENDIF

  !force the Z to be Zsol when using Pickles library
  IF (spec_type.EQ.'picks') THEN
     zmin = 1
     zmax = 1
  ENDIF

  !----------------------------------------------------------------!
  !-----------------Read in spectral libraries---------------------!
  !----------------------------------------------------------------!

  !read in wavelength array
  IF (spec_type.EQ.'basel') THEN
     OPEN(91,FILE=TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/basel.lambda',&
          STATUS='OLD',iostat=stat,ACTION='READ')
  ELSE IF (spec_type.EQ.'miles') THEN
     OPEN(91,FILE=TRIM(SPS_HOME)//'/SPECTRA/MILES/miles.lambda',&
          STATUS='OLD',iostat=stat,ACTION='READ')
  ELSE IF (spec_type.EQ.'picks') THEN
     OPEN(91,FILE=TRIM(SPS_HOME)//'/SPECTRA/Pickles/pickles.lambda',&
          STATUS='OLD',iostat=stat,ACTION='READ')
  ELSE IF (spec_type(1:5).EQ.'hrlib'.OR.spec_type.EQ.'rrlib') THEN
     OPEN(91,FILE=TRIM(SPS_HOME)//'/SPECTRA/HRLIB/'//spec_type(1:5)//&
          '.lambda',STATUS='OLD',iostat=stat,ACTION='READ')
  ENDIF
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: wavelength grid '//&
          'cannot be opened'
     STOP 
  ENDIF
  DO i=1,nspec
     READ(91,*) spec_lambda(i)
  ENDDO
  CLOSE(91)
  
  !read in basel logg and logt arrays
  !NB: these are the same for all spectral libraries
  OPEN(91,FILE=TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/basel_logt.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,ndim_logt
     READ(91,*) basel_logt(i)
  ENDDO
  CLOSE(91)
  OPEN(91,FILE=TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/basel_logg.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,ndim_logg
     READ(91,*) basel_logg(i)
  ENDDO
  CLOSE(91)

  !read in each metallicity
  DO z=zmin,zmax

     WRITE(zstype,'(F6.4)') zlegend(z)

     !read in the spectral library
     IF (spec_type.EQ.'basel') THEN
        OPEN(92,FILE=TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/basel_'//basel_str//'_z'&
             //zstype//'.spectra.bin',FORM='UNFORMATTED',&
             STATUS='OLD',iostat=stat,ACTION='READ',access='direct',&
             recl=nspec*ndim_logg*ndim_logt*4)
     ELSE IF (spec_type.EQ.'miles') THEN
        OPEN(92,FILE=TRIM(SPS_HOME)//'/SPECTRA/MILES/imiles_z'&
             //zstype//'.spectra.bin',FORM='UNFORMATTED',&
             STATUS='OLD',iostat=stat,ACTION='READ',access='direct',&
             recl=nspec*ndim_logg*ndim_logt*4)
     ELSE IF (spec_type.EQ.'picks') THEN
        OPEN(92,FILE=TRIM(SPS_HOME)//'/SPECTRA/Pickles/pickles.'&
             //'spectra.bin',FORM='UNFORMATTED',&
             STATUS='OLD',iostat=stat,ACTION='READ',access='direct',&
             recl=nspec*ndim_logg*ndim_logt*4)
     ELSE IF (spec_type(1:5).EQ.'hrlib'.OR.spec_type.EQ.'rrlib') THEN
        OPEN(92,FILE=TRIM(SPS_HOME)//'/SPECTRA/HRLIB/'//spec_type//'_z'&
             //zstype//'.spectra.bin',FORM='UNFORMATTED',&
             STATUS='OLD',iostat=stat,ACTION='READ',access='direct',&
             recl=nspec*ndim_logg*ndim_logt*4)
     ENDIF
     IF (stat.NE.0) THEN 
        WRITE(*,*) 'SPS_SETUP ERROR: '//spec_type//&
             ' spectral library cannot be opened ', zstype
        STOP 
     ENDIF

     READ(92,rec=1) rspeclib(:,z,:,:)
     CLOSE(92)

  ENDDO

  speclib = rspeclib

  !these stars are only included with BaSeL or MILES libraries
  IF (spec_type.EQ.'miles'.OR.spec_type.EQ.'basel') THEN

     !read in AGB Teff array for O-rich spectra
     OPEN(93,FILE=TRIM(SPS_HOME)//'/SPECTRA/AGB_spectra/Orich_teff_allZ_'//&
       isoc_type//'_'//spec_type//'.dat',STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: /SPECTRA/AGB_spectra/Orich_teff_allZ_'&
             //isoc_type//'_'//spec_type//'.dat cannot be opened'
        STOP 
     ENDIF
     !burn the header
     READ(93,*) char
     DO i=1,n_agb_o
        READ(93,*) dumr1, agb_logt_o(:,i)
     ENDDO
     CLOSE(93)
     agb_logt_o = LOG10(agb_logt_o)
     
     !read in AGB Teff array for C-rich spectra
     OPEN(94,FILE=TRIM(SPS_HOME)//'/SPECTRA/AGB_spectra/Crich_teff_allZ'//&
          '.dat',STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: /SPECTRA/AGB_spectra/'//&
             'Crich_teff_allZ.dat cannot be opened'
        STOP 
     ENDIF
     !burn the header
     READ(94,*) char
     DO i=1,n_agb_c
        READ(94,*) dumr1, agb_logt_c(i)
     ENDDO
     CLOSE(94)
     agb_logt_c = LOG10(agb_logt_c)
     
     !read in TP-AGB O-rich spectra
     OPEN(95,FILE=TRIM(SPS_HOME)//&
          '/SPECTRA/AGB_spectra/Orich_spec_all_'//spec_type//'.dat',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN

        WRITE(*,*) 'SPS_SETUP ERROR: /AGB_spectra/Orich_spec_all_'//&
             spec_type//'.dat '//'cannot be opened'
        STOP 
     ENDIF
     DO i=1,nspec
        READ(95,*) d1,agb_spec_o(i,:)
     ENDDO
     CLOSE(95)
     
     !read in TP-AGB C-rich spectra
     OPEN(96,FILE=TRIM(SPS_HOME)//&
          '/SPECTRA/AGB_spectra/Crich_spec_all_'//spec_type//'.dat',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: /AGB_spectra/SPECTRA/'//&
             'Crich_spec_all_'//spec_type//'.dat '//'cannot be opened'
        STOP 
     ENDIF
     DO i=1,nspec
        READ(96,*) d1,agb_spec_c(i,:)
     ENDDO
     CLOSE(96)
     
     !read in post-AGB Teff array
     OPEN(94,FILE=TRIM(SPS_HOME)//'/SPECTRA/Hot_spectra/ipagb.teff',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: Hot_spectra/ipagb.teff cannot be opened'
        STOP 
     ENDIF
     DO i=1,ndim_pagb
        READ(94,*) pagb_logt(i)
     ENDDO
     CLOSE(94)
     pagb_logt = LOG10(pagb_logt)

     !read in post-AGB spectra from Rauch 2003
     OPEN(97,FILE=TRIM(SPS_HOME)//'/SPECTRA/Hot_spectra/ipagb_'//&
          spec_type//'.spec_solar',STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: /SPECTRA/Hot_spectra/'//&
             'ipagb.spec_solar cannot be opened'
        STOP 
     ENDIF
     DO i=1,nspec
        READ(97,*) d1,pagb_spec(i,:,2)
     ENDDO
     CLOSE(97)
     OPEN(97,FILE=TRIM(SPS_HOME)//'/SPECTRA/Hot_spectra/ipagb_'//&
          spec_type//'.spec_halo',STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: /SPECTRA/Hot_spectra/'//&
             'ipagb.spec_halo cannot be opened'
        STOP 
     ENDIF
     DO i=1,nspec
        READ(97,*) d1,pagb_spec(i,:,1)
     ENDDO
     CLOSE(97)

     !read in WR Teff array
     OPEN(94,FILE=TRIM(SPS_HOME)//'/SPECTRA/Hot_spectra/CMFGEN_WN.teff',&
          STATUS='OLD',iostat=stat,ACTION= 'READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: Hot_spectra/iwr.teff cannot be opened'
        STOP 
     ENDIF
     DO i=1,ndim_wr
        READ(94,*) wr_logt(i)
     ENDDO
     CLOSE(94)

     !read in WR spectra from Smith et al. 2002
     OPEN(97,FILE=TRIM(SPS_HOME)//'/SPECTRA/Hot_spectra/CMFGEN_WN_Z020_'&
          //spec_type//'.spec',STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: Hot_spectra/iwr.spec '//&
             'cannot be opened'
        STOP 
     ENDIF
     DO i=1,nspec
        READ(97,*) d1,wr_spec(i,:)
     ENDDO
     CLOSE(97)

  ENDIF

  !----------------------------------------------------------------!
  !--------------------Read in isochrones--------------------------!
  !----------------------------------------------------------------!

  !read in all metallicities
  DO z=zmin,zmax

     n_isoc = 0
     WRITE(zstype,'(F6.4)') zlegend(z)

     !open Padova isochrones
     IF (isoc_type.EQ.'pdva') &
          OPEN(97,FILE=TRIM(SPS_HOME)//&
          '/ISOCHRONES/Padova/Padova2007/isoc_z'//&
          zstype//'.dat',STATUS='OLD', IOSTAT=stat,ACTION='READ')
     !open MIST isochrones
     IF (isoc_type.EQ.'mist') &
          OPEN(97,FILE=TRIM(SPS_HOME)//&
          'ISOCHRONES/MIST/isoc_z'//zstype//'.dat',STATUS='OLD',&
          IOSTAT=stat,ACTION='READ')
     !open BaSTI isochrones
     IF (isoc_type.EQ.'bsti') &
          OPEN(97,FILE=TRIM(SPS_HOME)//&
          'ISOCHRONES/BaSTI/NOVER/isoc_z'//zstype//'.dat',STATUS='OLD',&
          IOSTAT=stat,ACTION='READ')

     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: isochrone files cannot be opened'
        STOP 
     END IF
        
     DO i=1,nlines
           
        READ(97,*,IOSTAT=stat) char
        IF (stat.NE.0) GOTO 20
        
        IF (char.EQ.'#') THEN 
           m = 1
        ELSE 
           
           IF (m.EQ.1) n_isoc = n_isoc+1
           BACKSPACE(97)
           IF (m.GT.nm) THEN
              WRITE(*,*) 'SPS_SETUP ERROR: number of mass points GT nm'
              STOP
           ENDIF
           READ(97,*,IOSTAT=stat) logage,mini_isoc(z,n_isoc,m),&
                mact_isoc(z,n_isoc,m),logl_isoc(z,n_isoc,m),&
                logt_isoc(z,n_isoc,m),logg_isoc(z,n_isoc,m),&
                ffco_isoc(z,n_isoc,m),phase_isoc(z,n_isoc,m)
           IF (stat.NE.0) GOTO 20
           IF (m.EQ.1) timestep_isoc(z,n_isoc) = logage
           nmass_isoc(z,n_isoc) = nmass_isoc(z,n_isoc)+1
           m = m+1
           
        ENDIF
        
     ENDDO

     WRITE(*,*) 'SPS_SETUP ERROR: didnt finish reading in the isochrones!'
     STOP
     
20   CONTINUE
     CLOSE(97)

     IF (n_isoc.NE.nt) THEN 
        WRITE(*,*) 'SPS_SETUP ERROR: number of isochrones NE nt',n_isoc,nt
        STOP
     ENDIF

  ENDDO

  !----------------------------------------------------------------!
  !--------Read in dust emission spectra from Draine & Li----------!
  !----------------------------------------------------------------!
  
  !NB: we need to do this before other things b/c the wavelength
  !    array is modified here.

  DO k=0,6

     WRITE(sqpah,'(I1)') k
     OPEN(99,FILE=TRIM(SPS_HOME)//'/dust/DL07/DL07_MW3.1_'//sqpah//'0.dat',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN 
        WRITE(*,*) 'SPS_SETUP ERROR: error opening dust emission file'
        STOP
     ENDIF
     DO i=1,2  !burn the header
        READ(99,*)
     ENDDO
     DO i=1,ndim_dl07
        READ(99,*,IOSTAT=stat) lambda_dl07(i),dustem_dl07(i,:)
        IF (stat.NE.0) THEN 
           WRITE(*,*) 'SPS_SETUP ERROR: error during dust emission read'
           STOP
        ENDIF
     ENDDO
     CLOSE(99)
     lambda_dl07 = lambda_dl07*1E4  !convert to Ang
     
     !expand the wavelength array
     !i=1
     !DO WHILE (spec_lambda(i)/1E4.LT.10.)
     !   i=i+1
     !ENDDO
     !j=1
     !DO WHILE (lambda_dl07(j)/1E4.LE.10.) 
     !   j=j+1
     !ENDDO
     !spec_lambda(i:) = lambda_dl07(j-1:)
     
     !now interpolate the dust spectra onto the master wavelength array
     DO j=1,numin_dl07*2
        DO i=1,nspec
           !the dust models only extend to 1um
           IF (spec_lambda(i)/1E4.GT.1) THEN
              dustem2_dl07(i,k+1,j) = linterp(lambda_dl07,&
                   dustem_dl07(:,j),spec_lambda(i))
           ENDIF
        ENDDO
     ENDDO

  ENDDO

  !logarithmic wavelength grid used in smoothspec.f90
  dlstep = (LOG(spec_lambda(max_lam_smooth))-&
       LOG(spec_lambda(min_lam_smooth)))/nspec
  DO i=1,nspec
     lnlam(i) = i*dlstep+LOG(spec_lambda(min_lam_smooth))
  ENDDO

  !----------------------------------------------------------------!
  !-------------Read in circumstellar AGB dust models--------------!
  !----------------------------------------------------------------!

  !O-rich spectra
  OPEN(99,FILE=TRIM(SPS_HOME)//'/dust/dusty/Orich_dusty.spec',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN 
     WRITE(*,*) 'SPS_SETUP ERROR: error opening dusty models'
     STOP
  ENDIF

  !number of wavelength points in the AGB grid
  READ(99,*) nlam
  !read in the wavelength grid
  READ(99,*,IOSTAT=stat) lambda_dagb(1:nlam)

  lambda_dagb(1:nlam) = lambda_dagb(1:nlam)

  DO i=1,nteff_dagb
     DO j=1,ntau_dagb
        READ(99,*,IOSTAT=stat) teff_dagb(1,i), tau1_dagb(j)
        READ(99,*,IOSTAT=stat) fluxin_dagb(1:nlam)
        IF (stat.NE.0) THEN 
           WRITE(*,*) 'SPS_SETUP ERROR: error reading dusty models'
           STOP
        ENDIF
        !interpolate the dust spectra onto the master wavelength array
        DO k=1,nspec
           IF (spec_lambda(k).GT.lambda_dagb(1)) THEN
              flux_dagb(k,1,i,j) = linterp(lambda_dagb(1:nlam),&
                   fluxin_dagb(1:nlam),spec_lambda(k))
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  CLOSE(99)

  !C-rich spectra
  OPEN(99,FILE=TRIM(SPS_HOME)//'/dust/dusty/Crich_dusty.spec',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN 
     WRITE(*,*) 'SPS_SETUP ERROR: error opening dusty models'
     STOP
  ENDIF

  !number of wavelength points in the AGB grid
  READ(99,*) nlam
  !read in the wavelength grid
  READ(99,*,IOSTAT=stat) lambda_dagb(1:nlam)

  DO i=1,nteff_dagb
     DO j=1,ntau_dagb
        READ(99,*,IOSTAT=stat) teff_dagb(2,i), tau1_dagb(j)
        READ(99,*,IOSTAT=stat) fluxin_dagb(1:nlam)
        IF (stat.NE.0) THEN 
           WRITE(*,*) 'SPS_SETUP ERROR: error reading dusty models'
           STOP
        ENDIF
        !interpolate the dust spectra onto the master wavelength array
        DO k=1,nspec
           IF (spec_lambda(k).GT.lambda_dagb(1)) THEN
              flux_dagb(k,2,i,j) = linterp(lambda_dagb(1:nlam),&
                   fluxin_dagb(1:nlam),spec_lambda(k))
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  CLOSE(99)



  !----------------------------------------------------------------!
  !-------------------Set up magnitude info------------------------!
  !----------------------------------------------------------------!

  !read in Vega-like star (lambda, Flambda)
  !(this is actually a Kurucz (1992) model for Vega)
  OPEN(98,FILE=TRIM(SPS_HOME)//'/SPECTRA/A0V_KURUCZ_92.SED',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: SPECTRA/A0V_KURUCZ_92.SED cannot be opened'
     STOP 
  END IF  
  !burn the header
  READ(98,*)
  DO i=1,1221
     READ(98,*) tvega_lam(i), tvega_spec(i)
  ENDDO
  CLOSE(98)

  !interpolate the Vega spectrum onto the wavelength grid
  DO i=1,nspec
     vega_spec(i) = 10**linterp(LOG10(tvega_lam),&
          LOG10(tvega_spec+tiny_number),LOG10(spec_lambda(i)))
     !convert to fnu
     vega_spec(i) =  vega_spec(i)*spec_lambda(i)**2
     IF (spec_lambda(i).GT.tvega_lam(1221)) vega_spec(i)=tiny_number
  ENDDO

  !read in Solar spectrum; units are fnu, flux is appropriate for
  !deriving absolute magnitudes.  spectrum from STScI, extrapolated 
  !beyond 2.5um with a blackbody.
  OPEN(98,FILE=TRIM(SPS_HOME)//'/SPECTRA/SUN_STScI.SED',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: SPECTRA/SUN_STScI.SED cannot be opened'
     STOP 
  END IF  
  DO i=1,1221
     READ(98,*) tsun_lam(i), tsun_spec(i)
  ENDDO
  CLOSE(98)

  !interpolate the Solar spectrum onto the wavelength grid
  DO i=1,nspec
     sun_spec(i) = 10**linterp(LOG10(tsun_lam),LOG10(tsun_spec+tiny_number),&
          LOG10(spec_lambda(i)))
     IF (spec_lambda(i).GT.tsun_lam(1221)) sun_spec(i)=tiny_number
  ENDDO

  !read in and set up band-pass filters
  OPEN(99,FILE=TRIM(SPS_HOME)//'/data/allfilters.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  READ(99,*)

  DO i=1,nbands

     jj=0
     readlamb = 0.0
     readband = 0.0
     DO j=1,50000
        READ(99,*,iostat=stat) d1,d2
        IF (stat.NE.0) GOTO 909
        IF (jj.EQ.0) THEN
           jj = jj+1
           readlamb(jj) = d1
           !force the transmission to be GE 0
           readband(jj) = MAX(d2,0.0)
        ELSE
           !only read unique lambda points
           IF (readlamb(jj).NE.d1) THEN
              jj = jj+1
              readlamb(jj) = d1
              !force the transmission to be GE 0
              readband(jj) = MAX(d2,0.0)
           ENDIF
        ENDIF
     ENDDO

909  CONTINUE
     IF (j.EQ.50000) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: did not finish reading in filters!'
        STOP
     ENDIF
     IF (jj.EQ.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: error during filter definition read-in',i
        STOP
     ENDIF

     !interpolate the filter onto the master wavelength array
     DO j=1,nspec
        IF (spec_lambda(j).GE.readlamb(1).AND.&
             spec_lambda(j).LE.readlamb(jj)) &
             bands(i,j) = linterp(readlamb(1:jj),readband(1:jj),spec_lambda(j))
     ENDDO

     !normalize
     dumr1 = TSUM(spec_lambda,bands(i,:)/spec_lambda)
     !in this case the band is entirely outside the wavelength array
     IF (dumr1.EQ.0.0) dumr1=1.0 
     bands(i,:) = bands(i,:) / dumr1

     !compute absolute magnitude of the Sun in all filters
     magsun(i) = TSUM(spec_lambda,sun_spec*bands(i,:)/spec_lambda)
     IF (magsun(i).LT.2*tiny_number) THEN
        magsun(i) = 99.0
     ELSE
        magsun(i) = -2.5*LOG10(magsun(i)) - 48.60
     ENDIF

     !compute mags of Vega
     magvega(i) = TSUM(spec_lambda,vega_spec*bands(i,:)/spec_lambda)
     IF (magvega(i).LE.tiny_number) THEN
        magvega(i) = 99.0 
     ELSE
        magvega(i) = -2.5 * LOG10(magvega(i)) - 48.60
     ENDIF
     
     !put Sun magnitudes in the Vega system if keyword is set
     IF (compute_vega_mags.EQ.1.AND.magsun(i).NE.99.0) &
          magsun(i) = (magsun(i)-magsun(1)) - &
          (magvega(i)-magvega(1)) + magsun(1)

  ENDDO
  CLOSE(99)

  !----------------------------------------------------------------!
  !---------------Set up extinction curve indices------------------!
  !----------------------------------------------------------------!

  DO j=1,nspec
     x = 1E4/spec_lambda(j)
     IF (x.GT.12.) mwdindex(6)=j
     IF (x.GE.8.)  mwdindex(5)=j
     IF (x.GE.5.9) mwdindex(4)=j
     IF (x.GE.3.3) mwdindex(3)=j
     IF (x.GE.1.1) mwdindex(2)=j
     IF (x.GE.0.1) mwdindex(1)=j
  ENDDO

  !----------------------------------------------------------------!
  !----------Set up Witt & Gordon 2000 attenuation curves----------!
  !----------------------------------------------------------------!

  !read in the WG00 dust model attenuation curves
  OPEN(99,FILE=TRIM(SPS_HOME)//'/dust/alldirty_h.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  READ(99,*)
  READ(99,*)
  READ(99,*)
  DO i=1,18
     DO j=1,25
        READ(99,*) wglam(j), d1,wgtmp(j,i,1,:)
     ENDDO
  ENDDO
  CLOSE(99)
  OPEN(99,FILE=TRIM(SPS_HOME)//'/dust/alldirty_c.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  READ(99,*)
  READ(99,*)
  READ(99,*)
  DO i=1,18
     DO j=1,25
        READ(99,*) wglam(j), d1,wgtmp(j,i,2,:)
     ENDDO
  ENDDO
  CLOSE(99)

  !interpolate the WG00 models onto the spectral grid
  DO k=1,2
     DO i=1,18
        DO j=1,6
           DO n=1,nspec
              IF (spec_lambda(n).GT.wglam(25)) THEN
                 wgdust(n,i,j,k)=0.0
              ELSE IF (spec_lambda(n).LT.wglam(1)) THEN
                 wgdust(n,i,j,k) = wgtmp(1,i,k,j)
              ELSE
                 wgdust(n,i,j,k) = linterp(wglam,wgtmp(:,i,k,j),&
                      spec_lambda(n))
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  !----------------------------------------------------------------!
  !---------set up the redshift-age relation (age in Gyr)----------!
  !----------------------------------------------------------------!

  DO i=1,500
     a = (i-1)/499.*(1-1/1001.)+1/1001.
     zagespl(i,1) = 1/a-1  !redshift
     zagespl(i,2) = get_tuniv(zagespl(i,1))  !Tuniv(z)
  ENDDO

  !set Tuniv
  tuniv = get_tuniv(zero)

  !----------------------------------------------------------------!
  !-----------------read in index definitions----------------------!
  !----------------------------------------------------------------!

  OPEN(99,FILE=TRIM(SPS_HOME)//'/data/allindices.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,4
     READ(99,*)
  ENDDO
  DO i=1,nindsps
     READ(99,*,IOSTAT=stat) indexdefined(:,i)
     IF (stat.NE.0) THEN 
        WRITE(*,*) 'SPS_SETUP ERROR: error during index defintion read'
        STOP
     ENDIF
  ENDDO
  CLOSE(99)

  !----------------------------------------------------------------!
  !-----------------set up expanded time array---------------------!
  !----------------------------------------------------------------!

  DO i=1,ntfull
     IF (MOD(i-1,time_res_incr).EQ.0) THEN
        time_full(i) = timestep_isoc(zmin,(i-1)/time_res_incr+1)
     ELSE
        IF ((i-1)/time_res_incr+2.LT.nt) THEN
           d1 = (timestep_isoc(zmin,(i-1)/time_res_incr+2)-&
                timestep_isoc(zmin,(i-1)/time_res_incr+1))/time_res_incr 
        ENDIF
        time_full(i) = timestep_isoc(zmin,(i-1)/time_res_incr+1)+d1
        time_full(i) = time_full(i-1)+d1
     ENDIF
  ENDDO

  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!

  !set flag indicating that sps_setup has been run, initializing 
  !important common block vars/arrays
  check_sps_setup = 1

  IF (verbose.EQ.1) THEN
     WRITE(*,*) '      ...done'
     WRITE(*,*)
  ENDIF

END SUBROUTINE SPS_SETUP
