SUBROUTINE GETSPEC(pset,mact,logt,lbol,logg,phase,ffco,lmdot,wght,spec)

  ! Routine to return the spectrum of a star with an input logg,logt
  ! the phase flag determines if the star is a WR, P-AGB, or TP-AGB star
  ! ffco specifies the chemical composition of WR and TP-AGB stars

  ! This subroutine is a major bottleneck.  The spectra must be
  ! recomputed each time the IMF or isochrone parameters change.

  USE sps_vars; USE sps_utils, ONLY: locate
  IMPLICIT NONE

  REAL(SP), INTENT(in) :: mact,logt,lbol,logg,phase,ffco,wght,lmdot
  TYPE(PARAMS), INTENT(in) :: pset
  REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec  
  REAL(SP), DIMENSION(nspec) :: ispec
  REAL(SP) :: t,u,r2,test1,test2,test3,test4,loggi,teffi,rwr,twr,logt_cut
  INTEGER  :: klo,jlo,flag

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  spec  = tiny_number
  ispec = tiny_number
  flag  = 0
  logt_cut = MAX(wmb_logt(1),logt_wmb_hot)

  !compute radius squared (cm^2) 
  !we need to re-compute logg becuase we modify the logt,logl of AGB stars
  loggi = LOG10( gsig4pi*mact/lbol ) + 4*logt
  r2    = mact*msun*newton/10**loggi

  !post-AGB non-LTE model spectra from Rauch 2003
  !the H-Ni composition spectra are used here.
  !this library has two Zs, Solar and 0.1Solar, simply use one or the other
  IF (phase.EQ.6.0.AND.logt.GE.4.699) THEN
    
     flag = flag+1
     jlo = MIN(MAX(locate(pagb_logt,logt),1),ndim_pagb-1)
     t   = (logt-pagb_logt(jlo)) / &
          (pagb_logt(jlo+1)-pagb_logt(jlo))
     t = MIN(MAX(t,0.0),1.0) !no extrapolation
     klo = 1
     IF (zlegend(pset%zmet)/zsol.GT.0.5) klo=2 
     !the post-agb library is normalized to unity
     spec = lbol*((1-t)*pagb_spec(:,jlo,klo)+t*pagb_spec(:,jlo+1,klo))

  !WR library from Smith et al. 2002 (CMFGEN)
  !NB: there is currently no log(g) dependence in the WR spectra
  ELSE IF (phase.EQ.9.0.AND.use_wr_spectra.EQ.1) THEN
     
     IF (isoc_type.EQ.'mist') THEN
        !below is from Maeder (1990)
        !terminal velocity = 3000 km/s = 3E8 cm/s
        !kappa_es = 0.2*(1+X) cm2/g; assume X=0.73 (Solar)
        rwr = SQRT(r2)+3*0.2*(1+0.73)*ABS(10**lmdot)/8./mypi/3E8/yr2sc*msun
        twr = 10**logt*SQRT(SQRT(r2)/rwr)
        !this is from Smith et al. 2002
        twr = LOG10( 0.6*10**logt + 0.4*twr )
     ELSE
        twr = logt
     ENDIF

     IF (ffco.LT.10) THEN

        !WN spectra
        flag = flag+1
        jlo  = MIN(MAX(locate(wrn_logt,twr),1),ndim_wr-1)
        t    = (twr-wrn_logt(jlo))/(wrn_logt(jlo+1)-wrn_logt(jlo))
        t    = MIN(MAX(t,0.0),1.0) !no extrapolation
        !the WR library is normalized to unity
        spec = lbol*((1-t)*wrn_spec(:,jlo,pset%zmet)+&
             t*wrn_spec(:,jlo+1,pset%zmet))

     ELSE IF (ffco.GE.10) THEN
     
        !WC spectra
        flag = flag+1
        jlo  = MIN(MAX(locate(wrc_logt,twr),1),ndim_wr-1)
        t    = (twr-wrc_logt(jlo))/(wrc_logt(jlo+1)-wrc_logt(jlo))
        t    = MIN(MAX(t,0.0),1.0) !no extrapolation
        !the WR library is normalized to unity
        spec = lbol*((1-t)*wrc_spec(:,jlo,pset%zmet)+&
             t*wrc_spec(:,jlo+1,pset%zmet))

     ENDIF

  !O-rich TP-AGB spectra, from Lancon & Mouhcine 2002
  ELSE IF (phase.EQ.5.0.AND.logt.LT.3.6.AND.ffco.LE.1.0) THEN
     
     flag = flag+1
     jlo  = MAX(MIN(locate(agb_logt_o(pset%zmet,:),logt),n_agb_o-1),1)
     t    = (logt - agb_logt_o(pset%zmet,jlo)) / &
          (agb_logt_o(pset%zmet,jlo+1)-agb_logt_o(pset%zmet,jlo))
     t    = MIN(MAX(t,0.0),1.0) !no extrapolation

     !The spectra are Fdlambda, need to convert to Fdnu and 
     !interpolate in Teff.
     spec = lbol*spec_lambda*spec_lambda/clight * &
          ( (1-t)*agb_spec_o(:,jlo) + t*(agb_spec_o(:,jlo+1)) )

  !C-rich TP-AGB spectra
  ELSE IF (phase.EQ.5.0.AND.logt.LT.3.6.AND.ffco.GT.1.0) THEN
     
     flag = flag+1

     !use Aringer et al. (2009) synthetic spectra
     IF (cstar_aringer.EQ.1) THEN 
        jlo  = MAX(MIN(locate(agb_logt_car,logt),n_agb_car-1),1)
        t    = (logt - agb_logt_car(jlo)) / &
             (agb_logt_car(jlo+1)-agb_logt_car(jlo))
        t    = MIN(MAX(t,0.0),1.0) !no extrapolation
        !The spectra are Fnu/Lbol
        spec = lbol*( (1-t)*agb_spec_car(:,jlo) + t*(agb_spec_car(:,jlo+1)) )

     !use LW02 empirical spectra
     ELSE
        jlo  = MAX(MIN(locate(agb_logt_c,logt),n_agb_c-1),1)
        t    = (logt - agb_logt_c(jlo)) / &
             (agb_logt_c(jlo+1)-agb_logt_c(jlo))
        t    = MIN(MAX(t,0.0),1.0)
        !The spectra are Fdlambda, need to convert to Fdnu and 
        !interpolate in Teff.
        spec = lbol*spec_lambda*spec_lambda/clight * &
             ( (1-t)*agb_spec_c(:,jlo) + t*(agb_spec_c(:,jlo+1)) )
     ENDIF

  !use WMBasic grid from JJ Eldridge for T>25,000K MS stars
  ELSE IF (phase.EQ.0.0.AND.logt.GT.logt_cut) THEN

     flag = flag+1

     jlo = MIN(MAX(locate(wmb_logt,logt),1),ndim_wmb_logt-1)
     klo = MIN(MAX(locate(wmb_logg,loggi),1),ndim_wmb_logg-1)
     t   = (logt-wmb_logt(jlo)) / (wmb_logt(jlo+1)-wmb_logt(jlo))
     t   = MIN(MAX(t,0.0),1.0) !no extrapolation (this means >50K -> 50K)
     u   = (loggi-wmb_logg(klo))  / (wmb_logg(klo+1)-wmb_logg(klo))
     u   = MIN(MAX(u,0.0),1.0) !no extrapolation in logg

     spec = (1-t)*(1-u)*wmb_spec(:,pset%zmet,jlo,klo) + &
          t*(1-u)*wmb_spec(:,pset%zmet,jlo+1,klo) + &
          t*u*wmb_spec(:,pset%zmet,jlo+1,klo+1) + &
          (1-t)*u*wmb_spec(:,pset%zmet,jlo,klo+1)
     !the WMBasic library is normalized to unity
     spec = spec*lbol

  !use the primary library for the rest of the isochrone
  ELSE

     flag = flag+1

     !find the subgrid containing point i 
     jlo = MIN(MAX(locate(speclib_logt,logt),1),ndim_logt-1)
     klo = MIN(MAX(locate(speclib_logg,loggi),1),ndim_logg-1)
     t   = (logt-speclib_logt(jlo))/(speclib_logt(jlo+1)-speclib_logt(jlo))
     t   = MIN(MAX(t,0.0),1.0) !no extrapolation
     u   = (loggi-speclib_logg(klo))/(speclib_logg(klo+1)-speclib_logg(klo))
     u   = MIN(MAX(u,0.0),1.0) !no extrapolation

     test1 = speclib(whlam5000,pset%zmet,jlo,klo)
     test2 = speclib(whlam5000,pset%zmet,jlo+1,klo)
     test3 = speclib(whlam5000,pset%zmet,jlo,klo+1)
     test4 = speclib(whlam5000,pset%zmet,jlo+1,klo+1)
     
     !if all four components are zero, set the flag to zero
     IF ((test1.LE.tiny30.AND.test2.LE.tiny30.AND.&
          test3.LE.tiny30.AND.test4.LE.tiny30)) flag=0

     !catch stars that fall off part of the grid
     !the flux at 5000A should never be zero unless a spec is missing
     IF ((test1.LE.tiny30.OR.test2.LE.tiny30.OR.&
          test3.LE.tiny30.OR.test4.LE.tiny30).AND.flag.EQ.1) THEN

        IF (verbose.EQ.99) & 
             WRITE(*,'(" GETSPEC WARNING: Part of the '//&
             'point is off the grid: Z=",I2,'//&
             '" logT=",F5.2," logg=",F5.2," phase=",I2," lg IMF*L=",F5.2)') &
             pset%zmet,logt,loggi,INT(phase),LOG10(wght*lbol)

        !this is a very crude hack.  just pick one of the spectra
        IF (test1.GT.tiny30) ispec = speclib(:,pset%zmet,jlo,klo)
        IF (test2.GT.tiny30) ispec = speclib(:,pset%zmet,jlo+1,klo)
        IF (test3.GT.tiny30) ispec = speclib(:,pset%zmet,jlo,klo+1)
        IF (test4.GT.tiny30) ispec = speclib(:,pset%zmet,jlo+1,klo+1)

     ELSE

        !bilinear interpolation
        ispec = (1-t)*(1-u)*speclib(:,pset%zmet,jlo,klo) + &
             t*(1-u)*speclib(:,pset%zmet,jlo+1,klo) + &
             t*u*speclib(:,pset%zmet,jlo+1,klo+1) + &
             (1-t)*u*speclib(:,pset%zmet,jlo,klo+1)

     ENDIF

     !at long last the extra factor of 4pi (below) has been found!
     !see p244-245 of Collins' "Fundamentals of Stellar Astrophysics"
     spec = 4*mypi*4*mypi*r2/lsun * ispec

  ENDIF
 
  !make sure the spectrum never has any zeros or negative numbers
  spec = MAX(spec,tiny_number)
  
  IF (verbose.EQ.1) THEN
     !IF (flag.EQ.0.AND.(spec_type.EQ.'basel'.OR.spec_type.EQ.'ckc14').AND.&
     !     phase.NE.6.AND.phase.NE.9) THEN
     IF (flag.EQ.0..AND.wght.GT.0.0) THEN
        WRITE(*,'(" GETSPEC WARNING: point entirely off the grid: Z=",I2,'//&
             '" logT=",F5.2," logg=",F5.2," phase=",I2," lg IMF*L=",F5.2)') &
            pset%zmet,logt,loggi,INT(phase),LOG10(wght*lbol)
     ENDIF
  ENDIF

  IF (flag.GT.1) THEN
     WRITE(*,'(" GETSPEC ERROR: isochrone point assigned *two* spectra!")') 
  ENDIF

  !add circumstellar dust around AGB stars
  IF ((phase.EQ.4.OR.phase.EQ.5) &
       .AND.add_agb_dust_model.EQ.1.AND.pset%agb_dust.GT.tiny_number) THEN
     CALL ADD_AGB_DUST(pset%agb_dust,spec,mact,&
          logt,LOG10(lbol),logg,zlegend(pset%zmet),ffco,lmdot)
  ENDIF


  !pure blackbody; no longer used but kept here for posterity
  !teffi = 10**logt
  !spec = 15/mypi*lbol/clight*&
  !     (hck/teffi)*(hck/teffi)*(hck/teffi)*(hck/teffi) / &
  !     (spec_lambda*spec_lambda*spec_lambda)/ ( &
  !     EXP(hck/spec_lambda/teffi)-1 )

  
END SUBROUTINE GETSPEC
