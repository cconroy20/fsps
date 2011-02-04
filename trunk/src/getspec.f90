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
  INTEGER, DIMENSION(5) :: mizz
  REAL(SP) :: loggi,t,u,r2,tpagbtdiff,teffi, sumtest,tco
  INTEGER  :: klo,jlo,j,tzz,idum=-1

  spec = 0.0

  IF (spec_type.EQ.'miles') THEN
     IF (isoc_type.EQ.'pdva') mizz = (/6,12,17,20,22/)
     IF (isoc_type.EQ.'bsti') mizz = (/2,5,7,8,9/)
     tzz = mizz(zz)
  ELSE
     tzz = zz
  ENDIF

  tco = ffco
!  IF (tco.GT.1.0) THEN
!     IF (ran(idum).LT.0.3) tco = 0.0
!  ENDIF

  !convert lum, temp, and mass into surface gravity
  !we only really need to do this for the added BHB and SBS stars
  loggi = LOG10( gsig4pi*mact/lbol ) + 4*logt
  !compute radius squared [cm^2] 
  r2    = mact*msun*newton/10**loggi

  !post-AGB stars have very high logg.  For T<5E4 we need to use
  !the BaSeL library, which has maximum logg=5.5
  IF (loggi.GT.5.5) loggi = 5.5
  
  !stars hotter than 5E4 are assumed to be pure blackbodies
  IF (logt.GT.4.6990) THEN
     
     teffi = 10**logt
     !this correctly integrates to the total Lbol
     spec = 15/mypi*lbol/clight*&
          (hck/teffi)*(hck/teffi)*(hck/teffi)*(hck/teffi) / &
          (spec_lambda*spec_lambda*spec_lambda)/ ( &
          EXP(hck/spec_lambda/teffi)-1 )

  !O-rich TP-AGB spectra, from Lancon & Mouhcine 2002
  ELSE IF (phase.EQ.5.0.AND.logt.LT.3.6.AND.tco.LE.1.0) THEN
     
     jlo = MAX(MIN(locate(agb_logt_o(tzz,:),logt),n_agb_o-1),1)

     !never let the TP-AGB spectra be an extrapolation
     IF (logt.LT.agb_logt_o(tzz,1)) THEN
        tpagbtdiff = 0.0
     ELSE IF (logt.GT.agb_logt_o(tzz,n_agb_o)) THEN
        tpagbtdiff = agb_logt_o(tzz,jlo+1)-agb_logt_o(tzz,jlo)
     ELSE
        tpagbtdiff = logt - agb_logt_o(tzz,jlo)
     ENDIF

     !The spectra are Fdlambda, need to convert to Fdnu and 
     !interpolate in Teff.
     spec = lbol*spec_lambda*spec_lambda/clight * &
          ( agb_spec_o(:,jlo) + tpagbtdiff * (agb_spec_o(:,jlo+1) - &
          agb_spec_o(:,jlo))/(agb_logt_o(tzz,jlo+1)-agb_logt_o(tzz,jlo)) )
     
  !C-rich TP-AGB spectra, from Lancon & Mouhcine 2002
  ELSE IF (phase.EQ.5.0.AND.logt.LT.3.6.AND.tco.GT.1.0) THEN
     
     jlo = MAX(MIN(locate(agb_logt_c(tzz,:),logt),n_agb_c-1),1)
     
     !never let the TP-AGB spectra be an extrapolation
     IF (logt.LT.agb_logt_c(tzz,1)) THEN
        tpagbtdiff = 0.0
     ELSE IF (logt.GT.agb_logt_c(tzz,n_agb_c)) THEN
        tpagbtdiff = agb_logt_c(tzz,jlo+1)-agb_logt_c(tzz,jlo)
     ELSE
        tpagbtdiff = logt - agb_logt_c(tzz,jlo)
     ENDIF
     
     !The spectra are Fdlambda, need to convert to Fdnu and 
     !interpolate in Teff.
     spec = lbol*spec_lambda*spec_lambda/clight * &
          ( agb_spec_c(:,jlo) + tpagbtdiff * (agb_spec_c(:,jlo+1) - &
          agb_spec_c(:,jlo))/(agb_logt_c(tzz,jlo+1)-agb_logt_c(tzz,jlo)) )
     
  !use the primary library for the rest of the isochrone
  ELSE
     
     !find the subgrid containing point i 
     klo = MIN(MAX(locate(basel_logg,loggi),1),ndim_logg-1)
     jlo = MIN(MAX(locate(basel_logt,logt),1),ndim_logt-1)

     t   = (logt-basel_logt(jlo)) / &
          (basel_logt(jlo+1)-basel_logt(jlo))
     u   = (loggi-basel_logg(klo))   / &
          (basel_logg(klo+1)-basel_logg(klo))
     
     !catch stars that fall off the BaSeL grid
     sumtest = SUM(speclib(:,zz,jlo:jlo+1,klo:klo+1))
     IF (sumtest.EQ.0.0) write(*,*) 'GETSPEC WARNING: you '//&
          'have somehow managed to place (part of) an isochrone '//&
          'off of the spectral grid.  This is not good.  Z=',&
          zz,'(logt,logg)=',logt,loggi
     
     !bilinear interpolation over every spectral element
     !NB: extra factor of 4pi, that I can't really account for, 
     !but its needed for everything to work out.
     spec = 4*mypi*4*mypi*r2/lsun* ( &
          (1-t)*(1-u)*speclib(:,zz,jlo,klo) + &
            t*(1-u)*speclib(:,zz,jlo+1,klo) + &
              t*u*speclib(:,zz,jlo+1,klo+1) + &
            (1-t)*u*speclib(:,zz,jlo,klo+1) )
     
  ENDIF
  
END SUBROUTINE GETSPEC
