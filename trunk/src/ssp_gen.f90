!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Routine to calculate the evolution of a single stellar          !
!  population from a set of input theoretical isochrones and a     !
!  heterogeneous library of stellar spectra.  The code also allows !
!  for variation in the horizontal branch morphology, TP-AGB       !
!  phase, and the blue straggler population.  The output is a      !
!  time-dependent spectrum from the far-UV to the IR.              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! COMPARISONS/TESTS:
!   1. Mass(t) for both Salpeter & Chabrier IMFs are within <3%
!      of those from B&C03.
!   2. Color evolution in many bands brackets the BC03 and 
!      Maraston 05 models for all metallicities.
!   3. Spectral evolution is bracketed by the BC03 and M05 models
!
! PARAMETER RANGES:
!   1.  if (fbhb,sbs)<1E-3 then (fbhb,sbs)=0.0
!   2.  if abs(delt)>0.5 then abs(delt)=0.5
!-----------------------------------------------------------!
!-----------------------------------------------------------!

SUBROUTINE STITCH_PAGB(t,mini,mact,logl,logt,logg,wght,nmass,phase)

  !routine to add a few dummy stars to stitch the post-AGB tracks
  !onto the full isochrones

  USE sps_vars; USE nrtype
  IMPLICIT NONE

  INTEGER  :: i,nn,j, nexp=20
  REAL(SP), INTENT(inout), DIMENSION(nt,nm) :: mini,mact,logl,&
       logt,phase,logg
  REAL(SP), INTENT(inout), DIMENSION(nm) :: wght
  INTEGER, INTENT(inout),  DIMENSION(nt) :: nmass
  INTEGER, INTENT(in) :: t

  nn = nmass(t)

  pagb: DO i=1,nn
     IF (phase(t,i).EQ.6.0.AND.logg(t,i).EQ.-99.) THEN
        wght(i) = wght(i)/REAL(nexp+1)
        DO j=1,nexp
           logt(t,nn+j)  = logt(t,i-1)+j/REAL(nexp+1)*&
                (logt(t,i)-logt(t,i-1))
           logl(t,nn+j)  = logl(t,i-1)
           mini(t,nn+j)  = mini(t,i)
           mact(t,nn+j)  = mact(t,i)
           wght(nn+j)    = wght(i)
           phase(t,nn+j) = 6.0
        ENDDO
        nmass(t) = nmass(t)+nexp
        EXIT pagb
     ENDIF
  ENDDO pagb

END SUBROUTINE STITCH_PAGB

!-----------------------------------------------------------!
!---------------------Driver Routine------------------------!
!-----------------------------------------------------------!

SUBROUTINE SSP_GEN(pset,mass_ssp,lbol_ssp,spec_ssp)

  USE sps_vars; USE nrtype
  USE sps_utils
  IMPLICIT NONE
  
  INTEGER :: i=1, j=1
  !weight given to the entire horizontal branch
  REAL(SP) :: hb_wght   
  !array of IMF weights
  REAL(SP), DIMENSION(nm) :: wght
  !SSP spectrum
  REAL(SP), INTENT(inout), DIMENSION(nt,nspec) :: spec_ssp
  !Mass and Lbol info
  REAL(SP), INTENT(inout), DIMENSION(nt) :: mass_ssp, lbol_ssp

  !temp arrays for the isochrone data
  REAL(SP), DIMENSION(nt,nm) :: mini,mact,logl,logt,logg,ffco,phase
  !arrays holding the number of mass elements for each isochrone
  !and the age of each isochrone
  INTEGER, DIMENSION(nt)  :: nmass
  REAL(SP), DIMENSION(nt) :: time_ssp
  REAL(SP), DIMENSION(nspec) :: tspec

  !structure containing all necessary parameters
  !(TYPE objects defined in sps_vars.f90)
  TYPE(PARAMS), INTENT(in) :: pset

  CHARACTER(2) :: istr,istr2

  !-----------------------------------------------------------!
  !-------------------Read in user input----------------------!
  !-----------------------------------------------------------!

  IF (check_sps_setup.EQ.0) THEN
     WRITE(*,*) 'SSP_GEN ERROR0: '//&
          'SPS_SETUP must be run once before calling SSP_GEN. '
     STOP
  ENDIF

  !reset arrays
  hb_wght  = 0.
  wght     = 0.
  spec_ssp = 0.
  mass_ssp = 0.
  lbol_ssp = 0.

  !set up metallicity
  IF (pset%zmet.LT.1.OR.pset%zmet.GT.nz) THEN
     WRITE(*,*) 'SSP_GEN ERROR: metallicity outside of range',pset%zmet
     STOP
  ENDIF

  IF (imf_type.NE.0.AND.imf_type.NE.1.AND.imf_type.NE.2.&
       .AND.imf_type.NE.3) THEN
     WRITE(*,*) 'SSP_GEN ERROR: IMF type outside of range',imf_type
     STOP
  ENDIF

  !Dump Padova isochrones into temporary arrays
  mini     = mini_isoc(pset%zmet,:,:)  !initial mass
  mact     = mact_isoc(pset%zmet,:,:)  !actual (present) mass
  logl     = logl_isoc(pset%zmet,:,:)  !log(Lbol)
  logt     = logt_isoc(pset%zmet,:,:)  !log(Teff)
  logg     = logg_isoc(pset%zmet,:,:)  !log(g)
  ffco     = ffco_isoc(pset%zmet,:,:)  !is the TP-AGB star C-rich or O-rich?
  phase    = phase_isoc(pset%zmet,:,:) !flag indicating phase of evolution
  nmass    = nmass_isoc(pset%zmet,:)   !number of elements per isochrone
  time_ssp = timestep_isoc(pset%zmet,:)!age of each isochrone in log(yr)

  !dump IMF parameters into common block
  imf_alpha(1) = pset%imf1
  imf_alpha(2) = pset%imf2
  imf_alpha(3) = pset%imf3
  vdmc         = pset%vdmc

  !write for control
  IF (verbose.NE.0) THEN
     WRITE(*,*)
     WRITE(*,'("   Log(Z/Zsol): ",F6.3)') LOG10(zlegend(pset%zmet)/0.019)
     WRITE(*,'("   Fraction of blue HB stars: ",F6.3)') pset%fbhb
     WRITE(*,'("   Ratio of BS to HB stars  : ",F6.3)') pset%sbss
     WRITE(*,'("   Shift to TP-AGB [log(Teff),log(Lbol)]: ",F5.2,1x,F5.2)') &
          pset%delt, pset%dell
     IF (imf_type.EQ.2) THEN
        WRITE(*,'("   IMF: ",I1,", slopes= ",3F4.1)') &
             imf_type,imf_alpha
     ELSE IF (imf_type.EQ.3) THEN
        WRITE(*,'("   IMF: ",I1,", cut-off= ",F4.2)') imf_type,vdmc
     ELSE
        WRITE(*,'("   IMF: ",I1)') imf_type
     ENDIF
  ENDIF

  !-----------------------------------------------------------!
  !---------------------Generate SSPs-------------------------!
  !-----------------------------------------------------------!

  !loop over each isochrone
  DO i=1,nt

     !compute IMF-based weights
     CALL IMF_WEIGHT(mini(i,:),wght,nmass(i))
     
     !modify the horizontal branch
     !need the hb weight for the blue stragglers too
     IF (pset%fbhb.GT.0.0.OR.pset%sbss.GT.1E-3) &
          CALL MOD_HB(pset%fbhb,i,mini,mact,logl,logt,phase,&
          wght,hb_wght,nmass,time_ssp(i))

     !add in blue stragglers
     IF (time_ssp(i).GE.bhb_sbs_time.AND.pset%sbss.GT.1E-3) &
          CALL ADD_BS(pset%sbss,i,mini,mact,logl,logt,phase,&
          wght,hb_wght,nmass)

     !modify the TP-AGB stars and Post-AGB stars
     CALL MOD_AGB(pset%zmet,i,time_ssp,pset%delt,pset%dell,pset%pagb,&
          pset%redgb,nmass(i),logl,logt,phase,wght)

     !stitch the post-AGB tracks smoothly onto the full isochrone
     !IF (pset%pagb.NE.0.AND.MAXVAL(phase(i,:)).GE.5.0) &
     !     CALL STITCH_PAGB(i,mini,mact,logl,logt,logg,wght,nmass,phase)

     !write out isochrones for testing purposes
     IF (i.GE.10.AND.MOD(i,5).EQ.0.AND.1.EQ.0) THEN
        WRITE(istr,'(I2)') i
        WRITE(istr2,'(I2)') 10+pset%zmet
        OPEN(89,file='isoc_basti_Z'//istr2//'_t'//istr//'.dat')
        DO j=1,nmass(i)
           WRITE(89,*) time_ssp(i),mini(i,j),logl(i,j),logt(i,j),&
                wght(j),phase(i,j)
        ENDDO
        CLOSE(89)
     ENDIF

     !compute IMF-weighted mass of the SSP
     mass_ssp(i) = SUM(wght(1:nmass(i))*mact(i,1:nmass(i)))
     !add in remant masses
     CALL ADD_REMNANTS(mass_ssp(i),MAXVAL(mini(i,:)))

     !compute IMF-weighted bolometric luminosity (actually log(Lbol))
     lbol_ssp(i) = LOG10(SUM(wght(1:nmass(i))*10**logl(i,1:nmass(i))))

     !compute SSP spectrum
     spec_ssp(i,:) = 0.
     DO j=1,nmass(i)
        CALL GETSPEC(pset%zmet,mini(i,j),mact(i,j),logt(i,j),&
             10**logl(i,j),phase(i,j),ffco(i,j),tspec)
        spec_ssp(i,:) = wght(j)*tspec + spec_ssp(i,:)
     ENDDO

  ENDDO
  
END SUBROUTINE SSP_GEN


