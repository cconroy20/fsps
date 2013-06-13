SUBROUTINE GETSPEC(zz,mact,logt,lbol,logg,phase,ffco,spec)

  ! Routine to return the spectrum of a star with an input logg,logt
  ! the phase flag determines if the star is a WR, P-AGB, or TP-AGB star
  ! ffco specifies the chemical composition of WR and TP-AGB stars

  !This subroutine is a major bottleneck.  The spectra must be
  !recomputed each time the IMF or isochrone parameters change.

  USE sps_vars; USE nr, ONLY: locate
  IMPLICIT NONE

  REAL(SP), INTENT(in) :: mact,logt,lbol,logg,phase,ffco
  INTEGER,  INTENT(in) :: zz
  REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec  
  REAL(SP) :: t,u,r2,tpagbtdiff,sum1,sum2,sum3,sum4,loggi
  INTEGER  :: klo,jlo,flag

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  spec = 0.0
  flag = 0

  !compute radius squared (cm^2) 
  !we need to re-compute logg becuase we modify the logt,logl of AGB stars
  loggi = LOG10( gsig4pi*mact/lbol ) + 4*logt
  r2    = mact*msun*newton/10**loggi

  !post-AGB non-LTE model spectra from Rauch 2003
  !the H-Ni composition spectra are used here.
  !this library has two Zs, Solar and 0.1Solar, simply use one or the other
  IF (phase.EQ.6.0.AND.logt.GE.4.699.AND.spec_type.EQ.'basel') THEN
    
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
  !also use this library for very hot stars that are not labeled as 
  !post-AGB. Such stars come from the Padova isochrones, for which no 
  !phase information is available
  !NB: there is currently no Z or log(g) dependence in the WR spectra
  !NB: notice also that currently the WN and WC libraries are the same
  ELSE IF ((phase.EQ.9.0.AND.ffco.LT.10.AND.spec_type.EQ.'basel')&
       .OR.(phase.NE.6.0.AND.logt.GE.4.699)) THEN
     
     flag = 1
     jlo = MIN(MAX(locate(wr_logt,logt),1),ndim_wr-1)
     t   = (logt-wr_logt(jlo)) / &
          (wr_logt(jlo+1)-wr_logt(jlo))
     t = MIN(MAX(t,-1.0),1.0) !no extrapolation
     !the WR library is normalized to unity
     spec = lbol * ( wr_spec(:,jlo) + &
          t*(wr_spec(:,jlo+1)-wr_spec(:,jlo)) )

  !WC WR stars (C-rich), from Smith et al. 2002
  ELSE IF (phase.EQ.9.0.AND.ffco.GE.10.AND.spec_type.EQ.'basel') THEN
     
     flag = 1
     jlo = MIN(MAX(locate(wr_logt,logt),1),ndim_wr-1)
     t   = (logt-wr_logt(jlo)) / &
          (wr_logt(jlo+1)-wr_logt(jlo))
     t = MIN(MAX(t,-1.0),1.0) !no extrapolation
     !the WR library is normalized to unity
     spec = lbol * ( wr_spec(:,jlo) + &
          t*(wr_spec(:,jlo+1)-wr_spec(:,jlo)) )

  !O-rich TP-AGB spectra, from Lancon & Mouhcine 2002
  ELSE IF (phase.EQ.5.0.AND.logt.LT.3.6.AND.ffco.LE.1.0) THEN
     
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
  ELSE IF (phase.EQ.5.0.AND.logt.LT.3.6.AND.ffco.GT.1.0) THEN
     
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
     
     IF (verbose.EQ.1) THEN 
        !catch stars that fall off part of the grid
        sum1 = SUM(speclib(:,zz,jlo,klo))
        sum2 = SUM(speclib(:,zz,jlo+1,klo))
        sum3 = SUM(speclib(:,zz,jlo,klo+1))
        sum4 = SUM(speclib(:,zz,jlo+1,klo+1))
        IF ((sum1.EQ.0.0.OR.sum2.EQ.0.OR.sum3.EQ.0.OR.sum4.EQ.0)&
             .AND.phase.NE.6.0) &
             write(*,'("GETSPEC WARNING: A '//&
             'star is off the grid: Z=",I2,'//&
             '" logT=",F5.2," logg=",F5.2," phase=",F2.0)') &
             zz,logt,loggi,phase
     ENDIF
     
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
  spec = MAX(spec,tiny_number)
  
  IF (flag.EQ.0.AND.spec_type.EQ.'basel') THEN
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
