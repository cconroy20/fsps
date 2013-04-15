MODULE SPS_VARS

  ! module to set up most arrays and variables 

  IMPLICIT NONE
  SAVE

!define either BaSeL/Kurucz, MILES, or HR CaT spectral library
#define BASEL 1
#define MILES 0
!define either Padova or BaSTI isochrones
#define PADOVA 1
#define BASTI 0

  INTEGER, PARAMETER :: SP = KIND(1.0D0)
  INTEGER, PARAMETER :: RSP = KIND(1.0)

  !------Common parameters that may be altered by the user-------!
  
  !setup cosmology (WMAP7).  Used only for z(t) relation.
  REAL(SP) :: om0=0.26, ol0=0.74, H0=72.
  
  !controls the level of output (0 = no output to screen)
  INTEGER, PARAMETER :: verbose=0

  !Turn-on time for BHB and SBS phases, time is in log(yrs)
  REAL(SP), PARAMETER :: bhb_sbs_time=9.3
  
  !turn on/off convolution of SSP with P(Z)
  !NB: P(Z) convolution has not been tested in some time, use with caution
  INTEGER, PARAMETER :: pzcon=0
  
  !the factor by which we increase the time array
  INTEGER, PARAMETER :: time_res_incr=2

  !turn on/off computation of light-weighted stellar ages
  !NB: This is only partially implemented, and only for straight tau models
  INTEGER, PARAMETER :: compute_light_ages=0

  !turn on/off the Draine & Li 2007 dust emission model 
  INTEGER, PARAMETER :: add_dust_emission=1

  !turn on/off the addition of stellar remnants to the 
  !computation of stellar masses
  INTEGER, PARAMETER :: add_stellar_remnants=1

  !set attenuation-law for the diffuse ISM
  !0 - power-law attenuation.  See dust_index variable below
  !1 - MW extinction law, parameterized by Cardelli et al. 1989,
  !    with a UV bump strength parameterized by uvb (see params below)
  !2 - Calzetti attenuation law
  !3 - Witt & Gordon 2000 attenuation curve models
  INTEGER :: dust_type=0

  !IMF definition
  !0 = Salpeter (parameters defined above)
  !1 = Chabrier 2003 (parameters defined above)
  !2 = Kroupa 2001 (three slopes must be specified in imf_alpha)
  !3 = van Dokkum 2008 (parameter must be specified in imf_vdmc)
  !4 = Dave 2008 (parameter specified in imf_mdave)
  !5 = user-defined piece-wise power-law, specified in imf.dat
  INTEGER :: imf_type=2

  !flag specifying zero-point of magnitudes
  !0 - AB system
  !1 - Vega system
  INTEGER  :: compute_vega_mags=0

  !flag indicating whether or not the output colors
  !will be redshifted to the age of the Universe corresponding
  !to the age of the SSP or CSP 
  !(only works when using compsp.f90 to compute mags)
  !0 - colors redshifted to a fixed redshift, specified in parameter set
  !1 - colors redshifted according to the age of the SSP or CSP
  INTEGER :: redshift_colors=0

  !------------Pre-compiler defintions------------!
  
  !flag indicating type of isochrones to use
  !'bsti' = BaSTI, 'pdva' = Padova 2007
#if (BASTI)
  CHARACTER(4), PARAMETER :: isoc_type = 'bsti'
#else
  CHARACTER(4), PARAMETER :: isoc_type = 'pdva'
#endif

  !flag indicating type of spectral library to use
#if (MILES)
  CHARACTER(5), PARAMETER :: spec_type = 'miles'
#else
  CHARACTER(5), PARAMETER :: spec_type = 'basel'
#endif

  !number of metallicities in the isochrones
  !number of elements per stellar spectrum
#if (MILES)
  INTEGER, PARAMETER :: nz=5
  INTEGER, PARAMETER :: nspec=5252
#else
  INTEGER, PARAMETER :: nspec=1963
#if (BASTI)
  INTEGER, PARAMETER :: nz=10
#else   
  INTEGER, PARAMETER :: nz=22
#endif
#endif

  !flag indicating the type of normalization used in the BaSeL library
  !pdva = normalized to Padova isochrones
  !wlbc = normalized to Teff-color relations
  !NB: currently only the wlbc option is included in the public release
  CHARACTER(4), PARAMETER :: basel_str = 'wlbc'

  !---------Dimensions of various arrays----------!

  !You must change the number of bands here if
  !filters are added to allfilters.dat
  INTEGER, PARAMETER :: nbands=86
  !number of indices defined in allindices.dat
  INTEGER, PARAMETER :: nindsps=30
  
  !The following parameters should never be changed
   !unless you are changing the libraries

  !max dimension of mass and time arrays for isochrones
  INTEGER, PARAMETER :: nm=1500, nt=94
  !dimension of Kurucz/BaSeL atlas
  INTEGER, PARAMETER :: ndim_basel=1221
  !max number of lines to read in
  INTEGER, PARAMETER ::  nlines=100000
  !max number of lines in tabulated SFH
  INTEGER, PARAMETER :: ntabmax=10000
  !dimensions of BaSeL library
  INTEGER, PARAMETER :: ndim_logt=68, ndim_logg=19
  !number of O-rich AGB spectra
  INTEGER, PARAMETER :: n_agb_o=9
  !number of C-rich AGB spectra
  INTEGER, PARAMETER :: n_agb_c=5
  !number of post-AGB spectra
  INTEGER, PARAMETER :: ndim_pagb=14
  !number of WR spectra (WN only)
  INTEGER, PARAMETER :: ndim_wr=12
  !wavelength dimension of the Draine & Li 2007 dust model
  INTEGER, PARAMETER :: ndim_dl07=1001
  !number of Umin models from Drain & Li 2007 dust model
  INTEGER, PARAMETER :: numin_dl07=22

  !------------IMF-related Constants--------------!
  
  !Salpeter IMF index
  REAL(SP) :: salp_ind= 2.35
  !min/max masses for the IMF
  REAL(SP) :: imf_lower_limit = 0.08, imf_upper_limit=120.
  !Chabrier 2003 IMF parameters
  REAL(SP), PARAMETER :: chab_mc=0.08, chab_sigma2=0.69*0.69,&
       chab_ind=1.3
  !van Dokkum 2008 IMF parameters
  REAL(SP), PARAMETER :: vd_sigma2=0.69*0.69, vd_ah=0.0443,&
       vd_ind=1.3, vd_al=0.14, vd_nc=25.
  REAL(SP) :: mlim_bh=40.0, mlim_ns=8.5

  !-------------Physical Constants---------------!
  !-------in cgs units where applicable----------!

  !constant such that g = C MT^4/L
  REAL(SP), PARAMETER :: gsig4pi = 1/4.13E10
  !pi
  REAL(SP), PARAMETER :: mypi    = 3.14159265
  !hc/k (Ang*K)
  REAL(SP), PARAMETER :: hck     = 1.43878E8
  !speed of light (Ang/s)
  REAL(SP), PARAMETER :: clight  = 2.9979E18
  !hc^2/sigma_SB
  REAL(SP), PARAMETER :: hc2sig  = 0.105021
  !Solar mass in grams
  REAL(SP), PARAMETER :: msun    = 1.989E33
  !Solar luminosity in erg/s
  REAL(SP), PARAMETER :: lsun    = 3.839E33
  !Newton's constant
  REAL(SP), PARAMETER :: newton  = 6.67428E-8
  !cm in a pc
  REAL(SP), PARAMETER :: pc2cm   = 3.08568E18
  
  !other important parameters
  REAL(SP), PARAMETER :: huge_number = 1E33
  REAL(SP), PARAMETER :: tiny_number = 1E-33
  
  !---------------Common Block-------------------!
    
  INTEGER :: check_sps_setup = 0
  
  !IMF parameters for Kroupa 2001 IMF
  !the user does not set these vars explicitly.  They
  !are set in the PARAMS structure below and are 
  !copied internally
  REAL(SP), DIMENSION(3) :: imf_alpha=1.3
  !IMF cut-off for van Dokkum parameterization
  !the user does not set this var explicitly.  It
  !is set in the PARAMS structure below and 
  !copied internally
  REAL(SP) :: imf_vdmc  = 0.08
  !IMF transition mass for Dave parameterization
  !the user does not set this var explicitly.  It
  !is set in the PARAMS structure below and 
  !copied internally
  REAL(SP) :: imf_mdave = 0.5
  !parameters for user-defined IMF
  INTEGER :: n_user_imf = 0
  REAL(SP), DIMENSION(3,100) :: imf_user_alpha=0.0

  !environment variable for SPS home directory
  CHARACTER(250) :: SPS_HOME=''

  !Age of Universe in Gyr (set in sps_setup.f90, based on cosmo params)
  REAL(SP) :: tuniv=0.0
  
  !this specifies the size of the full time grid
  INTEGER, PARAMETER :: ntfull = time_res_incr*nt

  !array of index definitions
  REAL(SP), DIMENSION(7,nindsps) :: indexdefined=0.0

  !array holding MW extinction curve indices
  INTEGER, DIMENSION(6) :: mwdindex=0

  !array holding Witt & Gordon dust models
  !wgdust(lam,tau,model,homo/clump)
  REAL(SP), DIMENSION(nspec,18,6,2) :: wgdust=0.0

  !Index for P(Z) distribution.  1=closed box;
  !P(Z) = z^zpow*exp(-z/pmetals)  (see pz_convol.f90)
  !pmetals set in PARAMS structure
  REAL(SP) :: zpow=1.0

  !array holding redshift-age relation
  REAL(SP), DIMENSION(500,2) :: zagespl=0.0
 
  !array holding tabulated SFH 
  REAL(SP), DIMENSION(3,ntabmax) :: sfh_tab=0.0
  INTEGER :: ntabsfh=0

  !bandpass filters 
  REAL(SP), DIMENSION(nbands,nspec) :: bands
  !magnitude of the Sun in all filters
  REAL(SP), DIMENSION(nbands) :: magsun,magvega
  !Vega-like star spectrum for Vega magnitude zero-point
  !spectrum of Sun, for absolute mags of Sun
  REAL(SP), DIMENSION(nspec)  :: vega_spec=0.,sun_spec=0.0
  !common wavelength array
  REAL(SP), DIMENSION(nspec)  :: spec_lambda=0.0

  !arrays for stellar spectral information in HR diagram
  REAL(SP), DIMENSION(ndim_logt) :: basel_logt=0.0
  REAL(SP), DIMENSION(ndim_logg) :: basel_logg=0.0
  REAL(SP), DIMENSION(nspec,nz,ndim_logt,ndim_logg)  :: speclib=0.0
  REAL(RSP), DIMENSION(nspec,nz,ndim_logt,ndim_logg) :: rspeclib=0.0
  
  !AGB library (Lancon & Mouhcine 2002)
  REAL(SP), DIMENSION(nspec,n_agb_o) :: agb_spec_o=0.
  REAL(SP), DIMENSION(nz,n_agb_o)    :: agb_logt_o=0.
  REAL(SP), DIMENSION(nspec,n_agb_c) :: agb_spec_c=0.
  REAL(SP), DIMENSION(n_agb_c)       :: agb_logt_c=0.
  
  !post-AGB library (Rauch 2003)
  REAL(SP), DIMENSION(nspec,ndim_pagb,2) :: pagb_spec=0.0
  REAL(SP), DIMENSION(ndim_pagb)         :: pagb_logt=0.0

  !WR library (Smith et al. 2002)
  REAL(SP), DIMENSION(nspec,ndim_wr) :: wr_spec=0.0
  REAL(SP), DIMENSION(ndim_wr)       :: wr_logt=0.0

  !dust emission model (Draine & Li 2007)
  REAL(SP), DIMENSION(ndim_dl07)              :: lambda_dl07=0.0
  REAL(SP), DIMENSION(ndim_dl07,numin_dl07*2) :: dustem_dl07=0.0
  REAL(SP), DIMENSION(nspec,7,numin_dl07*2)   :: dustem2_dl07=0.0

  !arrays for the isochrone data
  REAL(SP), DIMENSION(nz,nt,nm) :: mact_isoc=0.,logl_isoc=0.,&
       logt_isoc=0.,logg_isoc=0.,ffco_isoc=0.,phase_isoc=0.
  REAL(SP), DIMENSION(nz,nt,nm) :: mini_isoc=0.

  !arrays holding the number of mass elements for each isochrone,
  !the age of each isochrone, and the metallicity of each isochrone
  INTEGER, DIMENSION(nz,nt)  :: nmass_isoc=0
  REAL(SP), DIMENSION(nz,nt) :: timestep_isoc=0.0
  REAL(SP), DIMENSION(nz)    :: zlegend=-99.
  
  !arrays for the full Z-dep SSP spectra
  REAL(SP), DIMENSION(nz,ntfull,nspec) :: spec_ssp_zz=0.
  REAL(SP), DIMENSION(nz,ntfull)       :: mass_ssp_zz=0.,lbol_ssp_zz=0.
  
  REAL(SP), DIMENSION(ntfull) :: time_full=0.

  !------------Define TYPE structures-------------!
  
  !structure for the set of parameters necessary to generate a model
  TYPE PARAMS
     REAL(SP) :: pagb=1.0,dell=0.0,delt=0.0,fbhb=0.0,sbss=0.0,tau=1.0,&
          const=0.0,tage=0.0,fburst=0.0,tburst=11.0,dust1=0.0,dust2=0.0,&
          logzsol=-0.2,zred=0.0,pmetals=0.02,imf1=1.3,imf2=2.3,imf3=2.3,&
          vdmc=0.08,dust_clumps=-99.,frac_nodust=0.0,dust_index=-0.7,&
          dust_tesc=7.0,frac_obrun=0.0,uvb=1.0,mwr=3.1,redgb=1.0,&
          dust1_index=-1.0,mdave=0.5,sf_start=0.0,sf_trunc=0.0,sf_theta=0.0,&
          duste_gamma=0.01,duste_umin=1.0,duste_qpah=3.5,fcstar=1.0,&
          masscut=150.0
     INTEGER :: zmet=1,sfh=0,wgp1=1,wgp2=1,wgp3=1,evtype=-1
  END TYPE PARAMS
  
  !structure for the output of the compsp routine
  TYPE COMPSPOUT
     REAL(SP) :: age=0.,mass_csp=0.,lbol_csp=0.,sfr=0.,mdust=0.0
     REAL(SP), DIMENSION(nbands) :: mags=0.
     REAL(SP), DIMENSION(nspec) :: spec=0.
  END TYPE COMPSPOUT
  
  !-----the following structures are not used in the public code-----!
  !--they are included here because some users of FSPS utilize them--!

  !structure for observational data
  TYPE OBSDAT
     REAL(SP)                    :: zred=0.0,logsmass=0.0
     REAL(SP), DIMENSION(nbands) :: mags=0.0,magerr=0.0
     REAL(SP), DIMENSION(nspec)  :: spec=0.0, specerr=99.
  END TYPE OBSDAT

  !structure for using P(z) in chi2
  INTEGER, PARAMETER :: npzphot   = 200
  TYPE TPZPHOT
     REAL(SP), DIMENSION(npzphot) :: zz=0.0,pz=0.0
  END TYPE TPZPHOT

  !used for Powell minimization
  TYPE(OBSDAT) :: powell_data

END MODULE SPS_VARS
