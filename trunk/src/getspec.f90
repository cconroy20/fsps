SUBROUTINE GETSPEC(zz,mini,mact,logt,lbol,phase,ffco,spec)

  ! Use the theoretical/empirical stellar spectra library
  ! BaSeL to compute IMF-weighted spectra of SSPs.
  ! For very hot stars, assume a Blackbody.  For TP-AGB
  ! stars, use the empirical library of Lancon & Mouhcine 2002.
  ! Luminosities are in solar units; Fluxes are in Fnu units (Lsun/Hz)
  ! BS spectra from M67 look like MS stars (Liu et al. 2008)

  !This subroutine is a major bottleneck.  We must re-compute
  !spectra each time the IMF or isochrone parameters change.

  USE sps_vars; USE nrtype
  USE nr, ONLY: locate, ran
  IMPLICIT NONE

  REAL(SP), INTENT(in) :: mini,mact,logt,lbol,phase,ffco
  INTEGER,  INTENT(in) :: zz
  REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec  
  REAL(SP) :: loggi,t,u,r2,tpagbtdiff,sumtest,tco
  INTEGER  :: klo,jlo,idum=-1,flag

  !-------------------------------------------------------------!

  spec = 0.0
  flag = 0

  tco = ffco
  !allows us to dilute the C star fraction
  !IF (tco.GT.1.0) THEN
  !   IF (ran(idum).LT.1.0) tco = 0.0
  !ENDIF

  !convert lum, temp, and mass into surface gravity
  !we only really need to do this for the added BHB and SBS stars
  loggi = LOG10( gsig4pi*mact/lbol ) + 4*logt
  !compute radius squared (cm^2) 
  r2    = mact*msun*newton/10**loggi

 
  !post-AGB non-LTE model spectra from Rauch 2003
  !the H-Ni composition spectra are used here.
  !this library has two Zs, Solar and 0.1Solar, simply use one or the other
  IF (phase.EQ.6.0.AND.logt.GE.4.699.AND.spec_type.NE.'miles') THEN
    
     flag = 1
     jlo = MIN(MAX(locate(pagb_logt,logt),1),ndim_pagb-1)
     t   = (logt-pagb_logt(jlo)) / &
          (pagb_logt(jlo+1)-pagb_logt(jlo))
     klo = 1
     IF (zlegend(zz)/0.0190.GT.0.5) klo=2 
     !the post-agb library is normalized to unity
     spec = lbol * ( pagb_spec(:,jlo,klo) + &
          t*(pagb_spec(:,jlo+1,klo)-pagb_spec(:,jlo,klo)) )

  !WN WR stars (N-rich), from Smith et al. 2002
  !also use this library for very hot stars that are not post-AGB,
  !such stars come from the Padova isochrones, for which no 
  !phase information is available
  ELSE IF ((phase.EQ.9.0.AND.tco.LT.10.AND.spec_type.NE.'miles').OR.&
       (phase.NE.6.0.AND.logt.GE.4.699)) THEN

     !NB: there is currently no Z or log(g) dependence in the WR spectra
     
     flag = 1
     jlo = MIN(MAX(locate(wr_logt,logt),1),ndim_wr-1)
     t   = (logt-wr_logt(jlo)) / &
          (wr_logt(jlo+1)-wr_logt(jlo))
     t = MIN(MAX(t,-1.0),1.0) !no extrapolation
     !the WR library is normalized to unity
     spec = lbol * ( wr_spec(:,jlo) + &
          t*(wr_spec(:,jlo+1)-wr_spec(:,jlo)) )

  !WC WR stars (C-rich), from Smith et al. 2002
  ELSE IF (phase.EQ.9.0.AND.tco.GE.10.AND.spec_type.NE.'miles') THEN
     
     flag = 1
     jlo = MIN(MAX(locate(wr_logt,logt),1),ndim_wr-1)
     t   = (logt-wr_logt(jlo)) / &
          (wr_logt(jlo+1)-wr_logt(jlo))
     t = MIN(MAX(t,-1.0),1.0) !no extrapolation
     !the WR library is normalized to unity
     spec = lbol * ( wr_spec(:,jlo) + &
          t*(wr_spec(:,jlo+1)-wr_spec(:,jlo)) )

  !O-rich TP-AGB spectra, from Lancon & Mouhcine 2002
  ELSE IF (phase.EQ.5.0.AND.logt.LT.3.6.AND.tco.LE.1.0) THEN
     
     flag = 1
     jlo = MAX(MIN(locate(agb_logt_o(zz,:),logt),n_agb_o-1),1)
     !never let the TP-AGB spectra be an extrapolation
     IF (logt.LT.agb_logt_o(zz,1)) THEN
        tpagbtdiff = 0.0
     ELSE IF (logt.GT.agb_logt_o(zz,n_agb_o)) THEN
        tpagbtdiff = agb_logt_o(zz,jlo+1)-agb_logt_o(zz,jlo)
     ELSE
        tpagbtdiff = logt - agb_logt_o(zz,jlo)
     ENDIF

     !The spectra are Fdlambda, need to convert to Fdnu and 
     !interpolate in Teff.
     spec = lbol*spec_lambda*spec_lambda/clight * &
          ( agb_spec_o(:,jlo) + tpagbtdiff * (agb_spec_o(:,jlo+1) - &
          agb_spec_o(:,jlo))/(agb_logt_o(zz,jlo+1)-agb_logt_o(zz,jlo)) )

  !C-rich TP-AGB spectra, from Lancon & Mouhcine 2002
  ELSE IF (phase.EQ.5.0.AND.logt.LT.3.6.AND.tco.GT.1.0) THEN
     
     flag = 1
     jlo = MAX(MIN(locate(agb_logt_c,logt),n_agb_c-1),1)
     !never let the TP-AGB spectra be an extrapolation
     IF (logt.LT.agb_logt_c(1)) THEN
        tpagbtdiff = 0.0
     ELSE IF (logt.GT.agb_logt_c(n_agb_c)) THEN
        tpagbtdiff = agb_logt_c(jlo+1)-agb_logt_c(jlo)
     ELSE
        tpagbtdiff = logt - agb_logt_c(jlo)
     ENDIF

     !The spectra are Fdlambda, need to convert to Fdnu and 
     !interpolate in Teff.
     spec = lbol*spec_lambda*spec_lambda/clight * &
          ( agb_spec_c(:,jlo) + tpagbtdiff * (agb_spec_c(:,jlo+1) - &
          agb_spec_c(:,jlo))/(agb_logt_c(jlo+1)-agb_logt_c(jlo)) )

  !use the primary library for the rest of the isochrone
  ELSE IF (logt.LT.4.699) THEN
     
     flag = 1
     !find the subgrid containing point i 
     klo = MIN(MAX(locate(basel_logg,loggi),1),ndim_logg-1)
     jlo = MIN(MAX(locate(basel_logt,logt),1),ndim_logt-1)

     t   = (logt-basel_logt(jlo)) / &
          (basel_logt(jlo+1)-basel_logt(jlo))
     u   = (loggi-basel_logg(klo))   / &
          (basel_logg(klo+1)-basel_logg(klo))
     
     !catch stars that fall off the grid
     sumtest = SUM(speclib(:,zz,jlo:jlo+1,klo:klo+1))
     IF (sumtest.EQ.0.0) write(*,*) 'GETSPEC WARNING: you '//&
          'have somehow managed to place part of an isochrone '//&
          'off of the spectral grid.  This is not good.  Z=',&
          zz,'(logt,logg)=',logt,loggi
     
     !bilinear interpolation over every spectral element
     !NB: extra factor of 4pi, that I can't explain,
     !but its needed for everything to work out.
     spec = 4*mypi*4*mypi*r2/lsun* ( &
          (1-t)*(1-u)*speclib(:,zz,jlo,klo) + &
          t*(1-u)*speclib(:,zz,jlo+1,klo) + &
          t*u*speclib(:,zz,jlo+1,klo+1) + &
          (1-t)*u*speclib(:,zz,jlo,klo+1) )
 
  ENDIF
 
  !make sure the spectrum never has any zero's
  spec = MAX(spec,0.0)
  
  IF (flag.EQ.0.AND.spec_type.NE.'miles') THEN
     WRITE(*,*) 'GETSPEC ERROR: isochrone point not assigned a spectrum',&
          logt,loggi,phase
  ELSE IF (flag.GT.1) THEN
     WRITE(*,*) 'GETSPEC ERROR: isochrone point assigned *two* spectra!'
  ENDIF

  !pure blackbody; no longer used but kept here for posterity
  !teffi = 10**logt
  !spec = 15/mypi*lbol/clight*&
  !     (hck/teffi)*(hck/teffi)*(hck/teffi)*(hck/teffi) / &
  !     (spec_lambda*spec_lambda*spec_lambda)/ ( &
  !     EXP(hck/spec_lambda/teffi)-1 )

  
END SUBROUTINE GETSPEC
