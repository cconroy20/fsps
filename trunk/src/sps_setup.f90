SUBROUTINE SPS_SETUP(zin)

  !read in isochrones and spectral libraries for all metallicities.
  !read in band-pass info and the spectrum for Vega.
  !Arrays are stored in a common block defined in sps_vars.f90

  !If zin=-1 then all metallicities are read in, otherwise only the
  !metallicity corresponding to zin in the look-up table zlegend.dat
  !is read.

  USE sps_vars; USE nrtype; USE sps_utils
  USE nr, ONLY : spline,splint
  IMPLICIT NONE
  INTEGER, INTENT(in) :: zin
  INTEGER :: stat,n,dumi1,i,j,m,jj,k
  INTEGER :: n_isoc,z,zmin,zmax
  CHARACTER(1) :: char
  CHARACTER(6) :: zstype
  REAL(SP) :: dumr1,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,logage
  REAL(SP) :: x,y,a,b,fa,fb,r
  REAL(SP), DIMENSION(1221) :: tvega_lam, tvega_spec
  REAL(SP), DIMENSION(50000) :: readlamb,readband,spl
  REAL(SP), DIMENSION(25) :: wglam
  REAL(SP), DIMENSION(25,18,2,6) :: wgtmp

  !----------------------------------------------------------------!

  IF (verbose.EQ.1) THEN 
     WRITE(*,*) 
     WRITE(*,*) '    Setting up SPS...'
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

  IF (isoc_type.NE.'pdva'.AND.isoc_type.NE.'bsti') THEN
     WRITE(*,*) 'SPS_SETUP ERROR: isoc_type var set to invalid type: ',isoc_type
     STOP
  ENDIF

  IF (basel_str.NE.'pdva'.AND.basel_str.NE.'wlbc') THEN
     WRITE(*,*) 'SPS_SETUP ERROR: basel_str var set to invalid type: ',basel_str
     STOP
  ENDIF

  IF (spec_type.NE.'basel'.AND.spec_type.NE.'miles') THEN
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

  IF (spec_type.EQ.'basel'.AND.nspec.NE.1221) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: spec_type = "basel" but nspec NE 1221!'
     STOP
  ENDIF

  !----------------------------------------------------------------!
  !----------------Read in metallicity values----------------------!
  !----------------------------------------------------------------!
  
  !units are simply metal fraction by mass (e.g. Z=0.0190 for Zsun)
  IF (isoc_type.EQ.'pdva') THEN
     OPEN(90,FILE=TRIM(SPS_HOME)//'/Padova/Padova2007/zlegend_'//&
          spec_type//'.dat',STATUS='OLD',iostat=stat,ACTION='READ')
  ELSE IF (isoc_type.EQ.'bsti') THEN
     OPEN(90,FILE=TRIM(SPS_HOME)//'/BaSTI/zlegend_'//&
          spec_type//'.dat',STATUS='OLD',iostat=stat,ACTION='READ')
  ENDIF
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR0: zlegend_*.dat '//&
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

  !----------------------------------------------------------------!
  !-----------------Read in HR diagram spectra---------------------!
  !----------------------------------------------------------------!

  !read in wavelength array
  IF (spec_type.EQ.'basel') THEN
     OPEN(91,FILE=TRIM(SPS_HOME)//'/BaSeL3.1/basel.lambda',&
          STATUS='OLD',iostat=stat,ACTION='READ')
  ELSE
     OPEN(91,FILE=TRIM(SPS_HOME)//'/MILES/miles.lambda',&
          STATUS='OLD',iostat=stat,ACTION='READ')
  ENDIF

  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR1: wavelength grid '//&
          'cannot be opened'
     STOP 
  END IF

  DO n=1,nspec
     READ(91,*) spec_lambda(n)
  ENDDO
  CLOSE(91)
  
  !read in basel logg and logt arrays
  !NB: these are the same as for the MILES library
  OPEN(91,FILE=TRIM(SPS_HOME)//'/BaSeL3.1/basel_logt.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,ndim_logt
     READ(91,*) basel_logt(i)
  ENDDO
  CLOSE(91)
  OPEN(91,FILE=TRIM(SPS_HOME)//'/BaSeL3.1/basel_logg.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,ndim_logg
     READ(91,*) basel_logg(i)
  ENDDO
  CLOSE(91)

  !read in each metallicity
  DO z=zmin,zmax

     WRITE(zstype,'(F6.4)') zlegend(z)

     !read in spectral library
     IF (spec_type.EQ.'basel') THEN
        OPEN(92,FILE=TRIM(SPS_HOME)//'/BaSeL3.1/basel_'//basel_str//'_z'&
             //zstype//'.spectra.bin',FORM='UNFORMATTED',&
             STATUS='OLD',iostat=stat,ACTION='READ',access='direct',&
             recl=nspec*ndim_logg*ndim_logt*4)
     ELSE
        OPEN(92,FILE=TRIM(SPS_HOME)//'/MILES/imiles_z'&
             //zstype//'.spectra.bin',FORM='UNFORMATTED',&
             STATUS='OLD',iostat=stat,ACTION='READ',access='direct',&
             recl=nspec*ndim_logg*ndim_logt*4)
     ENDIF
     IF (stat.NE.0) THEN 
        WRITE(*,*) 'SPS_SETUP ERROR2: '//spec_type//&
             ' spectral library cannot be opened', zstype
        STOP 
     ENDIF

     READ(92,rec=1) speclib(:,z,:,:)
     CLOSE(92)

  ENDDO
  
  !read in AGB Teff array for O-rich spectra
  OPEN(93,FILE=TRIM(SPS_HOME)//'/AGB_spectra/Orich_teff_allZ_'//&
       isoc_type//'.dat',STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR3: /AGB_spectra/Orich_teff_allZ_'&
          //isoc_type//'.dat cannot be opened'
     STOP 
  ENDIF
  !burn the header
  READ(93,*,IOSTAT=stat) char
  DO i=1,n_agb_o
     READ(93,*) dumr1, agb_logt_o(:,i)
  ENDDO  
  CLOSE(93)
  agb_logt_o = LOG10(agb_logt_o)

  !read in AGB Teff array for C-rich spectra
  OPEN(94,FILE=TRIM(SPS_HOME)//'/AGB_spectra/Crich_teff_allZ_'//&
       isoc_type//'.dat',STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR4: /AGB_spectra/Crich_teff_allZ_'//&
          isoc_type//'.dat cannot be opened'
     STOP 
  ENDIF
  !burn the header
  READ(94,*,IOSTAT=stat) char
  DO i=1,n_agb_c
     READ(94,*) dumr1, agb_logt_c(:,i)
  ENDDO  
  CLOSE(94)
  agb_logt_c = LOG10(agb_logt_c)

  !read in AGB O-rich spectra
  OPEN(95,FILE=TRIM(SPS_HOME)//&
       '/AGB_spectra/Orich_spec_all_'//spec_type//'.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR5: /AGB_spectra/Orich_spec_all_'//&
          spec_type//'.dat '//'cannot be opened'
     STOP 
  ENDIF
  DO i=1,nspec
     READ(95,*) d1,agb_spec_o(i,1),agb_spec_o(i,2),agb_spec_o(i,3),&
          agb_spec_o(i,4),agb_spec_o(i,5),agb_spec_o(i,6),&
          agb_spec_o(i,7),agb_spec_o(i,8),agb_spec_o(i,9)
  ENDDO
  CLOSE(95)
  
  !read in TP-AGB C-rich spectra
  OPEN(96,FILE=TRIM(SPS_HOME)//&
       '/AGB_spectra/Crich_spec_all_'//spec_type//'.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR6: /AGB_spectra/Crich_spec_all_'//&
          spec_type//'.dat '//'cannot be opened'
     STOP 
  ENDIF
  DO i=1,nspec
     READ(96,*) d1,agb_spec_c(i,1),agb_spec_c(i,2),agb_spec_c(i,3),&
          agb_spec_c(i,4),agb_spec_c(i,5)
  ENDDO
  CLOSE(96)

  !----------------------------------------------------------------!
  !--------------------Read in isochrones--------------------------!
  !----------------------------------------------------------------!

  !read in all metallicities
  DO z=zmin,zmax

     n_isoc = 0
     WRITE(zstype,'(F6.4)') zlegend(z)

     ! Add in non-evolving low mass stars (cf. Baraffe et al. 1998)
     ! NB: this is only strictly appropriate for solar metallicity, 
     ! but should be OK since they only contribute mass.
  
     DO i=1,nt
        
        mini_isoc(z,i,1:3) = (/0.10,0.11,0.13/)
        mact_isoc(z,i,1:3) = (/0.10,0.11,0.13/)
        logl_isoc(z,i,1:3) = (/-3.07,-2.93,-2.73/)
        logt_isoc(z,i,1:3) = LOG10((/2813.,2921.,3056./))
        logg_isoc(z,i,1:3) = (/5.251,5.218,5.169/)
        nmass_isoc(z,i)    = nmass_isoc(z,i)+3
        
     ENDDO

     !read in Padova isochrones
     IF (isoc_type.EQ.'pdva') THEN

        OPEN(97,FILE=TRIM(SPS_HOME)//&
             '/Padova/Padova2007/isoc_z'//&
             zstype//'.dat',STATUS='OLD', IOSTAT=stat,ACTION='READ')
        IF (stat.NE.0) THEN
           WRITE(*,*) 'SPS_SETUP ERROR7: Padova isochrone files cannot be opened'
           STOP 
        END IF
        
        DO i=1,nlines
           
           READ(97,*,IOSTAT=stat) char
           IF (stat.NE.0) GOTO 20
           
           IF (char.EQ.'#') THEN 
              m = 4
           ELSE 
              
              IF (m.EQ.4) n_isoc = n_isoc+1
              BACKSPACE(97)
              READ(97,*,IOSTAT=stat) logage,mini_isoc(z,n_isoc,m),&
                   mact_isoc(z,n_isoc,m),logl_isoc(z,n_isoc,m),&
                   logt_isoc(z,n_isoc,m),logg_isoc(z,n_isoc,m),&
                   ffco_isoc(z,n_isoc,m),phase_isoc(z,n_isoc,m)
              IF (stat.NE.0) GOTO 20
              IF (m.EQ.4) timestep_isoc(z,n_isoc) = logage
              nmass_isoc(z,n_isoc) = nmass_isoc(z,n_isoc)+1
              m = m+1
              
           ENDIF
        
        ENDDO

        WRITE(*,*) 'SPS_SETUP ERROR8: didnt finish reading in the isochrones!'
        STOP
     
20      CONTINUE
        CLOSE(97)

     !read in BaSTI isochrones
     ELSE IF (isoc_type.EQ.'bsti') THEN

        OPEN(97,FILE=TRIM(SPS_HOME)//&
             '/BaSTI/isoc_z'//zstype//'.dat',STATUS='OLD',&
             IOSTAT=stat,ACTION='READ')
        IF (stat.NE.0) THEN
           WRITE(*,*) 'SPS_SETUP ERROR7: BaSTI isochrone files cannot be opened'
           STOP 
        END IF

        DO i=1,nlines
           
           READ(97,*,IOSTAT=stat) char
           IF (stat.NE.0) GOTO 21
           
           IF (char.EQ.'#') THEN 
              m = 4
           ELSE
              
              IF (m.EQ.4) n_isoc = n_isoc+1
              BACKSPACE(97)
              READ(97,*,IOSTAT=stat) logage,mini_isoc(z,n_isoc,m),&
                   mact_isoc(z,n_isoc,m),logl_isoc(z,n_isoc,m),&
                   logt_isoc(z,n_isoc,m),phase_isoc(z,n_isoc,m)
              IF (stat.NE.0) GOTO 21
              IF (m.EQ.4) timestep_isoc(z,n_isoc) = logage
              nmass_isoc(z,n_isoc) = nmass_isoc(z,n_isoc)+1
              m = m+1
              
           ENDIF
        
        ENDDO

        WRITE(*,*) 'SPS_SETUP ERROR8: didnt finish reading in the isochrones!'
        STOP
     
21      CONTINUE
        CLOSE(97)

     ENDIF
     
  ENDDO

  !----------------------------------------------------------------!
  !-------------------Set up magnitude info------------------------!
  !----------------------------------------------------------------!

  !read in Vega-like star (lambda, Flambda)
  OPEN(98,FILE=TRIM(SPS_HOME)//'/src/A0V_KURUCZ_92.SED',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR9: A0V_KURUCZ_92.SED cannot be opened'
     STOP 
  END IF  
  !burn the header
  READ(98,*)

  DO i=1,1221
     READ(98,*) tvega_lam(i), tvega_spec(i)
  ENDDO
  CLOSE(98)

  !need to interpolate the Vega spectrum onto the MILES wavelength grid
  !for coding reasons, we need to "interpolate" even onto the BaSeL grid
  CALL spline(tvega_lam,tvega_spec,1.0e30_sp,&
       1.0e30_sp,spl(1:1221))
  DO i=1,nspec
     vega_spec(i) = splint(tvega_lam,tvega_spec,spl(1:1221),spec_lambda(i))
  ENDDO

  !the normalization is not important here, only the shape.
  vega_spec = vega_spec / clight * spec_lambda**2

  !read in and set up band-pass filters
  OPEN(99,FILE=TRIM(SPS_HOME)//'/src/allfilters.dat',&
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
           IF (d2.GE.0.0) readband(jj) = d2
        ELSE
           IF (readlamb(jj).NE.d1) THEN
              jj = jj+1
              readlamb(jj) = d1
              IF (d2.GE.0.0) readband(jj) = d2
           ENDIF
        ENDIF
     ENDDO

909  CONTINUE
     IF (j.EQ.50000) THEN
        WRITE(*,*) 'SPS_SETUP ERROR10: did not finish reading in filters!'
        STOP
     ENDIF

     !interpolate the filter onto the master wavelength array
     CALL spline(readlamb(1:jj),readband(1:jj),1.0e30_sp,&
       1.0e30_sp,spl(1:jj))
     DO j=1,nspec
        IF (spec_lambda(j).GE.readlamb(1).AND.&
             spec_lambda(j).LE.readlamb(jj)) &
             bands(i,j) = &
             splint(readlamb(1:jj),readband(1:jj),spl(1:jj),spec_lambda(j))
     ENDDO

     !normalize
     dumr1 = 0.0
     DO j=1,nspec-1
        dumr1 = dumr1 + (spec_lambda(j+1)-spec_lambda(j))*&
             (bands(i,j+1)/spec_lambda(j+1)+bands(i,j)/spec_lambda(j))/2.
     ENDDO
     bands(i,:) = bands(i,:) / dumr1

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

  !spline the WG00 models onto the spectral grid
  DO k=1,2
     DO i=1,18
        DO j=1,6
           CALL spline(wglam,wgtmp(:,i,k,j),1.0e30_sp,&
                1.0e30_sp,spl(1:25))
           DO n=1,nspec
              IF (spec_lambda(n).GT.wglam(25)) THEN
                 wgdust(n,i,j,k)=0.0
              ELSE IF (spec_lambda(n).LT.wglam(1)) THEN
                 wgdust(n,i,j,k) = wgtmp(1,i,k,j)
              ELSE
                 wgdust(n,i,j,k) = splint(wglam,wgtmp(:,i,k,j),spl(1:25),&
                      spec_lambda(n))
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!

  !set up the redshift-age relation (age in Gyr)
  DO i=1,500
     a = (i-1)/499.*(1-1/1001.)+1/1001.
     zagespl(i,1) = 1/a-1
     zagespl(i,2) = get_tuniv(zagespl(i,1))
  ENDDO
  CALL spline(zagespl(:,2),zagespl(:,1),1.0e30_sp,&
             1.0e30_sp,zagespl(:,3))

  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!

  !read in index definitions
  OPEN(99,FILE=TRIM(SPS_HOME)//'/src/allindices.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,7
     READ(99,*)
  ENDDO
  DO i=1,nindsps
     READ(99,*) indexdefined(:,i)
  ENDDO
  CLOSE(99)

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
