SUBROUTINE MOD_HB(f_bhb,t,mini,mact,logl,logt,logg,phase, &
     wght,hb_wght,nmass,hbtime)

  !routine to modify the horizontal branch to include bluer 
  !stars.  The fiducial HB stars are identified and then a fraction
  !f_bhb uniformly redistributed from the red clump to 10^4 K 
  !We also want to be able to call this routine when f_bhb=0 
  !in order to count the weight of stars on the HB.

  !see e.g. Sarajedini et al. 2007.  This treatment of BHB 
  !produces spectra similar to the Maraston 2005 models
  !see also Jimenez et al. 2004

  !Note that the parameter bhb_sbs_time, set in sps_vars.f90,
  !sets the turn-on time for this modification.

  USE sps_vars
  IMPLICIT NONE

  REAL(SP), INTENT(inout), DIMENSION(nt,nm) :: mini,mact,&
       logl,logt,logg,phase
  REAL(SP), INTENT(inout), DIMENSION(nm) :: wght
  REAL(SP), DIMENSION(nm) :: tphase=0.0
  INTEGER, INTENT(inout), DIMENSION(nt) :: nmass
  REAL(SP), INTENT(inout) :: hb_wght
  INTEGER, INTENT(in) :: t
  REAL(SP), INTENT(in) :: f_bhb, hbtime

  !number of blue HB to add per HB star
  !(not important b/c their total weight remains fixed)
  INTEGER, PARAMETER :: nhb=10
  INTEGER :: j, i, flip=0, tnhb
  REAL(SP) :: tgrad=0., hblum=-999.,minteff=1E6
  REAL(SP), DIMENSION(nhb) :: dumarr=0.

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  hblum   = -999.
  flip    = 0
  hb_wght = 0.
  tphase  = phase(t,:)

  !we need to count the total number of HB stars in 
  !these isochrones.  Also the minimum Teff for the HB
  IF (isoc_type.EQ.'bsti'.OR.isoc_type.EQ.'mist') THEN
     tnhb = 0
     minteff=1E6
     DO i=1,nm
        IF (tphase(i).EQ.3) THEN
           tnhb=tnhb+1
           !need this delta(mass) cut to remove the stars 
           !descending from the TRGB in the MIST models
           IF (logt(t,i).LT.minteff.AND.(mini(t,i+1)-mini(t,i)).GT.1E-6)  &
                minteff=logt(t,i)
        ENDIF
     ENDDO
     i=1
  ENDIF

  DO j=2,nm

     IF (isoc_type.EQ.'pdva') THEN

        IF (flip.NE.-1) THEN 
        
           tgrad = (logl(t,j)-logl(t,j-1)) / (mini(t,j)-mini(t,j-1))
           
           !if the lum jump is negative and large, we're on the HB
           IF (tgrad.LE.-5E2.AND.logl(t,j-1).GT.2.5) THEN
              flip = 1
              hblum = logl(t,j)
           ENDIF
        
           IF (flip.EQ.1) THEN

              !modify the HB
              !once lum increases by 0.1 dex we've left the HB
              IF (ABS(logl(t,j)-hblum).LT.0.1) THEN 

                 !keep track of total HB weight
                 hb_wght  = hb_wght+wght(j)
                 
                 !Blue HB stars have to be old
                 IF (f_bhb.GT.1E-3.AND.hbtime.GE.bhb_sbs_time) THEN

                    !add blue HB stars (their mass and Lbol remain the same)
                    mini(t,nmass(t)+1:nmass(t)+nhb)  = dumarr+mini(t,j)
                    mact(t,nmass(t)+1:nmass(t)+nhb)  = dumarr+mact(t,j)
                    logl(t,nmass(t)+1:nmass(t)+nhb)  = dumarr+logl(t,j)
                    phase(t,nmass(t)+1:nmass(t)+nhb) = dumarr+8.
                    DO i=1,nhb
                       !distribute Teff uniformly to high T
                       logt(t,nmass(t)+i) = logt(t,j)+(4.2-logt(t,j))*i/REAL(nhb)
                       !logt(t,nmass(t)+i) = 4.0+(4.3-4.0)*i/REAL(nhb)
                       !compute logg
                       logg(t,nmass(t)+i) = LOG10( gsig4pi*mact(t,nmass(t)+i)/&
                            10**logl(t,nmass(t)+i) ) + 4*logt(t,nmass(t)+i) 

                    ENDDO
                    wght(nmass(t)+1:nmass(t)+nhb)   = dumarr + &
                         f_bhb*wght(j)/nhb
                    !modify the weight of the existing HB stars
                    wght(j)   = wght(j) * (1-f_bhb)                 
                    !update number of stars in the isochrone
                    nmass(t) = nmass(t)+nhb
                 ENDIF

              ELSE 
                 flip = -1
              ENDIF
           
           ENDIF
        
        ENDIF

     ELSE IF ((isoc_type.EQ.'bsti'.OR.isoc_type.EQ.'mist') &
          .AND.tphase(j).EQ.3) THEN

        !keep track of total HB weight
        hb_wght  = hb_wght+wght(j)
        
        !Blue HB stars have to be old
        !here, we're adding one additional BHB stars per CHeB star
        !we're also putting the original BHB star to the RC, so that
        !if f_bhb is small, then the actual BHB contribution is small
        !if f_bhb<1E-4, then the default MIST BHB is used
        IF (f_bhb.GT.1E-4.AND.hbtime.GE.bhb_sbs_time) THEN

           !update number of stars in the isochrone
           nmass(t) = nmass(t)+1
           !add blue HB stars (their mass and Lbol remain the same)
           mini(t,nmass(t))  = mini(t,j)
           mact(t,nmass(t))  = mact(t,j)
           logl(t,nmass(t))  = logl(t,j)
           phase(t,nmass(t)) = 8.
           logt(t,j) = minteff 
           !distribute Teff uniformly to high T
           logt(t,nmass(t)) = logt(t,j)+(4.5-logt(t,j))*i/REAL(tnhb)
           i=i+1
           wght(nmass(t))   = f_bhb*wght(j)
           !modify the weight of the existing HB stars
           wght(j)   = wght(j) * (1-f_bhb)                
           !compute logg
           logg(t,nmass(t)) = LOG10( gsig4pi*mact(t,nmass(t))/&
                10**logl(t,nmass(t)) ) + 4*logt(t,nmass(t)) 
              
        ENDIF

     ENDIF

  ENDDO

  IF (nmass(t).GT.nm) THEN
     WRITE(*,*) 'MOD_HB ERROR: number of mass points GT nm'
     STOP
  ENDIF

  RETURN

END SUBROUTINE MOD_HB
