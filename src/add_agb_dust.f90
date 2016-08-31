FUNCTION COMPUTE_TAU1(cstar,mact,logt,logl,logg,zz,lmdot)

  !routine to compute tau at 1um from input
  !stellar parameters.  See Villaume et al. (2015)
  
  USE sps_vars
  IMPLICIT NONE
  INTEGER, INTENT(in) :: cstar
  REAL(SP), INTENT(in)  :: mact,logt,logl,logg, zz,lmdot
  REAL(SP) :: compute_tau1, radius, period, vexp, vexp_max
  REAL(SP) :: rin, delta, kappa, delta_agb, mdot

  !-----------------------------------------------------------!
  !write(*,*) 'in function: ', cstar
  !dust-to-gas ratio depends on C/O
  IF (cstar.EQ.1) THEN
     delta_agb = 0.0025
  ELSE
     delta_agb = 0.01
  ENDIF

  !extinction coefficient depends on grain composition (C/O)
  IF (cstar.EQ.1) THEN
     !sum of AmC and SiC (90%+10%)
     kappa = 3200.
  ELSE
     !Sil
     kappa = 3000.
  ENDIF
  
  !stellar radius in solar units
  radius = SQRT(mact*msun*newton/10**logg) / rsun
 
  !fundamental pulsation period in days
  period = 10**( -2.07 + 1.94*LOG10(radius)-0.9*LOG10(mact) )
  
  !expansion velocity in km/s
  vexp     = -13.5 + 0.056*period
  !metallicity dependent maximum (Bressan et al. 1998)
  !this implicitly introduces a (strong) metallicity
  !dependence on the dust-to-gas ratio (see below)
  !vexp_max = 6.5*zz/0.008 + 0.00226*period
  vexp     = MAX(MIN(vexp,15.0),3.0)

  !mass-loss rate in Msun/yr
  IF (use_isoc_mdot.EQ.1) THEN
     IF (isoc_type.NE.'mist') THEN
        WRITE(*,*) 'ADD_AGB_DUST ERROR: use_isoc_mist=1 but isoc_type NE MIST!'
        STOP
     ENDIF
     !use Mdot from the isochrone files
     mdot = MIN(10**lmdot,1E-4)
  ELSE
     !see Vassiliadis & Wood (1993) for details
     IF (period.LT.500.) THEN
        IF (mact.LT.2.5) THEN
           mdot = 10**(-11.4+0.0123*period)
        ELSE
           mdot = 10**(-11.4+0.0125*(period-100*(mact-2.5)))
        ENDIF
     ELSE
        !superwind phase
        mdot = 10**logl/vexp*1.93E3*yr2sc/clight
     ENDIF
  ENDIF

  !inner radius (assuming Td=1100K for C-rich and 
  !Td=700 for O-rich) in cm
  IF (cstar.EQ.1) THEN
    rin = 1.92E12 * (10**logl)**0.5
  ELSE 
    rin = 4.74E12 * (10**logl)**0.5
  ENDIF

  !dust-to-gas ratio
  delta = delta_agb * vexp**2 / 225. * (10**logl/1E4)**(-0.6)

  !finally, compute tau (and convert to cgs)
  compute_tau1 = kappa * delta * (mdot*msun/yr2sc) &
       / rin / (4*mypi) / (vexp*1E5)

  !WRITE(*,'(13ES11.2)') mact,logl,logg,radius,period,mdot,rin,&
  !     vexp,delta,compute_tau1

END FUNCTION COMPUTE_TAU1

!------------------------------------------------------------!
!------------------------------------------------------------!
!------------------------------------------------------------!

SUBROUTINE ADD_AGB_DUST(weight,tspec,mact,logt,logl,logg,zz,&
     tco,lmdot)

  !routine to add a circumstellar dust shell to AGB stars
  !computes the optical depth at 1um from stellar parameters
  !then looks up the corresponding DUSTY model given the C/O
  !ratio and Teff.
  
  USE sps_vars; USE sps_utils, ONLY: locate, smoothspec
  IMPLICIT NONE

  REAL(SP), DIMENSION(nspec), INTENT(inout) :: tspec
  REAL(SP), INTENT(in)  :: weight,mact,logt,logl,logg,zz,tco,lmdot
  INTEGER :: cstar,jlo,klo
  REAL(SP) :: tau1,dj,dk, compute_tau1, loggi
  REAL(SP), DIMENSION(nspec) :: dusty
  
  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  IF (tco.GT.1) THEN
     cstar=1 
  ELSE 
     cstar=0
  ENDIF

  !we need to recompute loggi here because the BaSTI
  !isochrones do not come with logg pre-tabulated
    IF (isoc_type.EQ.'bsti') THEN
        loggi = LOG10( gsig4pi*mact/10**logl ) + 4*logt
    ELSE
        loggi = logg
    ENDIF
   !compute tau1 based on input stellar parameters
   tau1 = compute_tau1(cstar,mact,logt,logl,loggi,zz,lmdot)

  !allow the user to manually adjust the tau1 value
  !by an overall scale factor
  tau1 = tau1*weight

  IF (tau1.LE.tiny_number) RETURN

  !find dusty model given tau1,tco,Teff
  jlo = MIN(MAX(locate(teff_dagb(cstar+1,:),10**logt),1),&
       nteff_dagb-1)
  klo = MIN(MAX(locate(tau1_dagb(cstar+1,:),LOG10(tau1)),1),ntau_dagb-1)

  dj   = (10**logt-teff_dagb(cstar+1,jlo)) / &
       (teff_dagb(cstar+1,jlo+1)-teff_dagb(cstar+1,jlo))
  dk   = (LOG10(tau1)-tau1_dagb(cstar+1,klo)) / &
       (tau1_dagb(cstar+1,klo+1)-tau1_dagb(cstar+1,klo))

  !don't allow extrapolation off the grid
  dj = MIN(MAX(dj,-1.0),1.0)
  dk = MIN(MAX(dk,-1.0),1.0)

  !bilinear interpolation
  dusty = (1-dj)*(1-dk)*flux_dagb(:,cstar+1,jlo,klo) + &
       dj*(1-dk)*flux_dagb(:,cstar+1,jlo+1,klo) + &
       dj*dk*flux_dagb(:,cstar+1,jlo+1,klo+1) + &
       (1-dj)*dk*flux_dagb(:,cstar+1,jlo,klo+1) 

  !implement the dusty spectra (which are in units of 
  !flux_out/flux_in) into the AGB spectra
!  tspec = tspec * dusty
  CALL SMOOTHSPEC(spec_lambda,tspec,10000.d0,30000.d0,100000000.d0)
  tspec = tspec * dusty

END SUBROUTINE ADD_AGB_DUST
