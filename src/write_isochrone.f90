SUBROUTINE WRITE_ISOCHRONE(outfile,pset)

  !routine to write all isochrones and CMDs at a given metallicity

  USE sps_vars; USE sps_utils, ONLY : getmags,getspec,imf_weight
  IMPLICIT NONE

  INTEGER :: i,tt,zz
  TYPE(PARAMS), INTENT(in) :: pset
  CHARACTER(100), INTENT(in)  :: outfile
  CHARACTER(51)  :: fmt
  REAL(SP) :: dz=0.0,loggi,hb_wght
  REAL(SP), DIMENSION(nspec)  :: spec
  REAL(SP), DIMENSION(nm)     :: wght
  REAL(SP), DIMENSION(nbands) :: mags
  !temp arrays for the isochrone data
  REAL(SP), DIMENSION(nt,nm)  :: mini,mact,logl,logt,logg,ffco,phase
  INTEGER, DIMENSION(nt)      :: nmass

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  hb_wght = 0.0
  wght    = 0.0
  zz      = pset%zmet

  fmt = '(F7.4,1x,F8.4,1x,F14.9,1x,6(F8.4,1x),000(F7.3,1x))'
  WRITE(fmt(38:40),'(I3,1x,I4)') nbands

  OPEN(40,FILE=TRIM(SPS_HOME)//'/OUTPUTS/'//TRIM(outfile)//'.cmd',&
       STATUS='REPLACE')
  WRITE(40,*) '# age log(Z) mass logl logt logg '//&
       'phase composition log(weight) mags'
       
  !transfer isochrones into temporary arrays
  mini  = mini_isoc(zz,:,:)  !initial mass
  mact  = mact_isoc(zz,:,:)  !actual (present) mass
  logl  = logl_isoc(zz,:,:)  !log(Lbol)
  logt  = logt_isoc(zz,:,:)  !log(Teff)
  logg  = logg_isoc(zz,:,:)  !log(g)
  ffco  = ffco_isoc(zz,:,:)  !is the TP-AGB star C-rich or O-rich?
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

     !modify the TP-AGB stars and Post-AGB stars
     CALL MOD_GB(zz,tt,timestep_isoc(zz,:),pset%delt,&
          pset%dell,pset%pagb,pset%redgb,nmass(tt),logl,logt,phase,wght)


     DO i=1,nmass_isoc(zz,tt)
        
        !get the spectrum
        CALL GETSPEC(pset,mact(tt,i),logt(tt,i),10**logl(tt,i),&
             logg(tt,i),phase(tt,i),ffco(tt,i),spec)
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
             mini(tt,i),logl(tt,i),logt(tt,i),loggi,phase(tt,i),&
             ffco(tt,i),LOG10(wght(i)),mags
        
     ENDDO
        
  ENDDO

  CLOSE(40)

END SUBROUTINE WRITE_ISOCHRONE
