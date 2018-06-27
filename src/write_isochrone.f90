SUBROUTINE WRITE_ISOCHRONE(outfile,pset)

  !routine to write all isochrones and CMDs at a given metallicity
  !note that the output age grid is the native spacing, not boosted
  !by the parameter time_res_incr

  USE sps_vars
  USE sps_utils, ONLY : getmags,getspec,imf_weight,mod_hb,mod_gb,add_bs
  IMPLICIT NONE

  INTEGER :: i,tt,zz
  TYPE(PARAMS), INTENT(in) :: pset
  CHARACTER(100), INTENT(in)  :: outfile
  CHARACTER(60)  :: fmt
  REAL(SP) :: dz=0.0,loggi,hb_wght
  REAL(SP), DIMENSION(nspec)  :: spec
  REAL(SP), DIMENSION(nm)     :: wght
  REAL(SP), DIMENSION(nbands) :: mags
  !temp arrays for the isochrone data
  REAL(SP), DIMENSION(nt,nm)  :: mini,mact,logl,logt,logg,&
       ffco,phase,lmdot
  INTEGER, DIMENSION(nt)      :: nmass

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  hb_wght = 0.0
  wght    = 0.0
  zz      = pset%zmet

  fmt = '(F7.4,1x,F8.4,1x,F14.9,1x,F14.9,1x,7(F8.4,1x),000(F7.3,1x))'
  WRITE(fmt(47:49),'(I3,1x,I4)') nbands

  OPEN(40,FILE=TRIM(SPS_HOME)//'/OUTPUTS/'//TRIM(outfile)//'.cmd',&
       STATUS='REPLACE')
  WRITE(40,*) '# age log(Z) mini mact logl logt logg '//&
       'phase composition log(weight) log(mdot) mags'
       
  !transfer isochrones into temporary arrays
  mini  = mini_isoc(zz,:,:)  !initial mass
  mact  = mact_isoc(zz,:,:)  !actual (present) mass
  logl  = logl_isoc(zz,:,:)  !log(Lbol)
  logt  = logt_isoc(zz,:,:)  !log(Teff)
  logg  = logg_isoc(zz,:,:)  !log(g)
  ffco  = ffco_isoc(zz,:,:)  !is the TP-AGB star C-rich or O-rich?
  lmdot = lmdot_isoc(zz,:,:) !log Mdot
  phase = phase_isoc(zz,:,:) !flag indicating phase of evolution
  nmass = nmass_isoc(zz,:)   !number of elements per isochrone

  DO tt=1,nt

     !compute IMF-based weights
     CALL IMF_WEIGHT(mini(tt,:),wght,nmass(tt))

     !modify the horizontal branch
     !need the hb weight for the blue stragglers too
     IF (pset%fbhb.GT.0.0.OR.pset%sbss.GT.1E-3) &
          CALL MOD_HB(pset%fbhb,tt,mini,mact,logl,logt,logg,phase,&
          wght,hb_wght,nmass,timestep_isoc(zz,tt))

     !add in blue stragglers
     IF (timestep_isoc(zz,tt).GE.bhb_sbs_time.AND.pset%sbss.GT.1E-3) &
          CALL ADD_BS(pset%sbss,tt,mini,mact,logl,logt,logg,phase,&
          wght,hb_wght,nmass)

     !modify the RGB and/or AGB stars
     CALL MOD_GB(zz,tt,timestep_isoc(zz,:),pset%delt,&
          pset%dell,pset%pagb,pset%redgb,pset%agb,nmass(tt),logl,logt,phase,wght)

     DO i=1,nmass(tt)
        
        !get the spectrum
        CALL GETSPEC(pset,mact(tt,i),logt(tt,i),10**logl(tt,i),&
             logg(tt,i),phase(tt,i),ffco(tt,i),lmdot(tt,i),wght(i),spec)
        !calculate magnitudes
        CALL GETMAGS(dz,spec,mags)

        IF (isoc_type.EQ.'bsti') THEN
           loggi = LOG10( gsig4pi*mact_isoc(zz,tt,i)/&
                logl(tt,i) ) + 4*logt(tt,i)
        ELSE
           loggi = logg(tt,i)
        ENDIF

        !write results to file
        WRITE(40,fmt) timestep_isoc(zz,tt),LOG10(zlegend(zz)),&
             mini(tt,i),mact(tt,i),logl(tt,i),logt(tt,i),loggi,phase(tt,i),&
             ffco(tt,i),LOG10(wght(i)),lmdot(tt,i),mags
        
     ENDDO

  ENDDO

  CLOSE(40)

END SUBROUTINE WRITE_ISOCHRONE
