SUBROUTINE SBF(pset,outfile)

  !routine to read in an isochrone (logt,logl,z) and produce 
  !SBFs for each point.  SBF magnitudes are a light-weighted average
  !of the stellar luminosities over stellar mass

  USE sps_vars
  USE sps_utils, ONLY : imf_weight,mod_hb,add_bs,mod_gb,getmags,getspec
  IMPLICIT NONE

  CHARACTER(100), INTENT(in) :: outfile
  TYPE(PARAMS), INTENT(in)   :: pset
  INTEGER       :: i,j
  CHARACTER(34) :: fmt
  REAL(SP)      :: zero=0.0,hb_wght
  REAL(SP), DIMENSION(nspec)  :: tspec,tspec2,spec1,spec2
  REAL(SP), DIMENSION(nbands) :: mags
  REAL(SP), DIMENSION(nm)     :: wght
  REAL(SP), DIMENSION(nt,nm)  :: mini,mact,logl,logt,logg,ffco,phase,lmdot
  INTEGER, DIMENSION(nt)      :: nmass
  REAL(SP), DIMENSION(nt)     :: time

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  !set up the format 
  fmt = '(F7.4,1x,3(F8.4,1x),000(F7.3,1x))'
  WRITE(fmt(21:23),'(I3,1x,I4)') nbands

  !reset arrays
  hb_wght  = 0.
  wght     = 0.
 
  OPEN(56,FILE=TRIM(SPS_HOME)//'/OUTPUTS/'//TRIM(outfile)//&
       '.mags',STATUS='REPLACE')
  DO i=1,8 !write a dummy header
     WRITE(56,*) '#'
  ENDDO

  !transfer isochrones into temporary arrays
  mini  = mini_isoc(pset%zmet,:,:)  !initial mass
  mact  = mact_isoc(pset%zmet,:,:)  !actual (present) mass
  logl  = logl_isoc(pset%zmet,:,:)  !log(Lbol)
  logt  = logt_isoc(pset%zmet,:,:)  !log(Teff)
  logg  = logg_isoc(pset%zmet,:,:)  !log(g)
  ffco  = ffco_isoc(pset%zmet,:,:)  !is the TP-AGB star C-rich or O-rich?
  lmdot = lmdot_isoc(pset%zmet,:,:) !log Mdot
  phase = phase_isoc(pset%zmet,:,:) !flag indicating phase of evolution
  nmass = nmass_isoc(pset%zmet,:)   !number of elements per isochrone
  time  = timestep_isoc(pset%zmet,:)!age of each isochrone in log(yr)

  DO i=1,nt

     !compute IMF-based weights
     CALL IMF_WEIGHT(mini(i,:),wght,nmass(i))
     
     !modify the horizontal branch
     !need the hb weight for the blue stragglers too
     IF (pset%fbhb.GT.0.0.OR.pset%sbss.GT.1E-3) &
          CALL MOD_HB(pset%fbhb,i,mini,mact,logl,logt,logg,phase,&
          wght,hb_wght,nmass,time(i))

     !add in blue stragglers
     IF (time(i).GE.bhb_sbs_time.AND.pset%sbss.GT.1E-3) &
          CALL ADD_BS(pset%sbss,i,mini,mact,logl,logt,logg,phase,&
          wght,hb_wght,nmass)

     !modify the TP-AGB stars and Post-AGB stars
     CALL MOD_GB(pset%zmet,i,time,pset%delt,pset%dell,pset%pagb,&
          pset%redgb,pset%agb,nmass(i),logl,logt,phase,wght)
 
     spec1 = 0.0
     spec2 = 0.0
     DO j=1,nmass(i)
           
        !get spectrum of ith star
        CALL GETSPEC(pset,mact(i,j),logt(i,j),10**logl(i,j),logg(i,j),&
             phase(i,j),ffco(i,j),lmdot(i,j),wght(j),tspec)

        !compute first and second moments of flux for
        !all stars and also by evolutionary phase
        spec2 = spec2 + wght(j)*tspec**2
        spec1 = spec1 + wght(j)*tspec

     ENDDO

     !compute the SBF
     tspec2 = spec2/spec1

     CALL GETMAGS(zero,tspec2,mags)
     WRITE(56,fmt) time(i),0.0,0.0,0.0,mags

  ENDDO

  CLOSE(56)

END SUBROUTINE SBF
 
