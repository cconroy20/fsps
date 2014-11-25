MODULE SPS_VARS

  ! module to set up most arrays and variables 

  IMPLICIT NONE
  SAVE

!-------set the spectral library------!
#define BASEL 1
#define MILES 0
! "CKC14" currently under development.  do not use!
#define CKC14 0

!------set the isochrone library------!
#define PADOVA 1
#define BASTI 0
! "MIST" currently under development.  do not use!
#define MIST 0


  !--------------------------------------------------------------!
  !--------------------------------------------------------------!

  !note that "SP" actually means double precision; this is a hack
  !to turn the nr routines into DP
  INTEGER, PARAMETER :: SP = KIND(1.d0)

  !------Common parameters that may be altered by the user-------!
  
  !setup cosmology (WMAP7).  Used only for z(t) relation.
  REAL(SP) :: om0=0.27, ol0=0.73, H0=72.
  !define solar metallicity
  REAL(SP), PARAMETER :: zsol=0.0190 
  
  !controls the level of output
  !0 = minimal output to screen.
  !1 = lots of output to screen.  useful for debugging.
  INTEGER, PARAMETER :: verbose=0

  !flag specifying TP-AGB normalization scheme
  !0 = default Padova 2007 isochrones
  !1 = Conroy & Gunn 2010 normalization
  !2 = Villaume, Conroy, Johnson 2014 normalization
  INTEGER :: tpagb_norm_type=2

  !turn-on time for BHB and SBS phases, time is in log(yrs)
  REAL(SP), PARAMETER :: bhb_sbs_time=9.3
  
  !turn on/off convolution of SSP with P(Z) (pz_convol.f90)
  !NB: pz_convol.f90 has not been tested in some time, use with caution
  INTEGER :: pzcon=0
  
  !the factor by which we increase the time array
  INTEGER, PARAMETER :: time_res_incr=2

  !turn on/off computation of light-weighted stellar ages
  !NB: currently only works with sfh=1,4 options
  INTEGER :: compute_light_ages=0

  !turn on/off the Draine & Li 2007 dust emission model 
  INTEGER :: add_dust_emission=1

  !turn on/off the AGB circumstellar dust model
  !see Villaume et al. (2014) for details
  INTEGER :: add_agb_dust_model=1

  !turn on/off a Cloudy-based nebular emission model (cont+lines)
  INTEGER :: add_neb_emission=0
  !turn on/off the nebular continuum component (automatically 
  !turned off if the above is set to 0)
  INTEGER  :: add_neb_continuum=1
  !minimum resolution (in velocity) for nebular lines, based 
  !on the resolution of the spectral libraries.
  REAL(SP) :: neb_res_min=1.0

  !turn on/off IGM absorption a la Madau (1995)
  INTEGER :: add_igm_absorption=0

  !turn on/off the addition of stellar remnants to the 
  !computation of stellar masses
  INTEGER :: add_stellar_remnants=1

  !if set, use a simpler, algorithm to smooth
  !the spectra.  Accurate to ~0.1% and somewhat faster than the 
  !correct approach.  NB: one should be careful when choosing
  !to run the slow version, as the accuracy depends on the min/max
  !wavelength parameters.  Contact me if you are intersted in this feature.
  INTEGER :: smoothspec_fast=1

  !if set, smooth the spectrum in velocity space, otherwise
  !smooth in Angstrom space (in all cases the width of the 
  !kernel is a sigma, not FWHM)
  INTEGER :: smooth_velocity=1

  !set attenuation-law for the diffuse ISM
  !0 = power-law attenuation.  See dust_index variable below
  !1 = MW extinction law, parameterized by Cardelli et al. 1989,
  !    with a UV bump strength parameterized by uvb (see params below)
  !2 = Calzetti attenuation law
  !3 = Witt & Gordon 2000 attenuation curve models
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
  !0 = AB system
  !1 = Vega system
  INTEGER  :: compute_vega_mags=0

  !output wavelengths in air (rather than vac) if set to 1
  INTEGER :: vactoair_flag=0

  !flag indicating whether or not the output colors
  !will be redshifted to the age of the Universe corresponding
  !to the age of the SSP or CSP 
  !(only works when using compsp.f90 to compute mags)
  !0 = colors redshifted to a fixed redshift, specified in parameter set
  !1 = colors redshifted according to the age of the SSP or CSP
  INTEGER :: redshift_colors=0

  !------------Pre-compiler defintions------------!
  
  !flag indicating type of isochrones to use
  !'bsti' = BaSTI, 'pdva' = Padova 2007
#if (BASTI)
  CHARACTER(4), PARAMETER :: isoc_type = 'bsti'
  INTEGER, PARAMETER :: nt=94
#elif (MIST)
  CHARACTER(4), PARAMETER :: isoc_type = 'mist'
  INTEGER, PARAMETER :: nt=107
#else
  CHARACTER(4), PARAMETER :: isoc_type = 'pdva'
  INTEGER, PARAMETER :: nt=94
#endif

  !flag indicating type of spectral library to use
#if (MILES)
  CHARACTER(5), PARAMETER :: spec_type = 'miles'
#elif (CKC14)
  !CHARACTER(5), PARAMETER :: spec_type = 'ckc14'
  CHARACTER(6), PARAMETER :: spec_type = 'ckc14z'
#else
  CHARACTER(5), PARAMETER :: spec_type = 'basel'
#endif

  !number of metallicities in the isochrones
  !number of elements per stellar spectrum
#if (MILES)
  INTEGER, PARAMETER :: nz=5
  INTEGER, PARAMETER :: nspec=5994
#elif (CKC14)
  INTEGER, PARAMETER :: nz=1
  INTEGER, PARAMETER :: nspec=47378   ! 47378, 26500
#else
  INTEGER, PARAMETER :: nspec=1963
#if (BASTI)
  INTEGER, PARAMETER :: nz=10
#elif (MIST)
  INTEGER, PARAMETER :: nz=1
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
  !kr06=61, kr02,kr03,kr04=101, kr01,kr11=102, normal=122
  INTEGER, PARAMETER :: nbands=122
  !number of indices defined in allindices.dat
  INTEGER, PARAMETER :: nindx=30
  
  !The following parameters should never be changed
  !unless you are changing the libraries

  !max dimension of array for each isochrone
  INTEGER, PARAMETER :: nm=1500
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
  !number of WR spectra
  INTEGER, PARAMETER :: ndim_wr=12
  !wavelength dimension of the Draine & Li 2007 dust model
  INTEGER, PARAMETER :: ndim_dl07=1001
  !number of Umin models from Drain & Li 2007 dust model
  INTEGER, PARAMETER :: numin_dl07=22
  !parameters for circumstellar dust models
  INTEGER, PARAMETER :: ntau_dagb=50, nteff_dagb=6
  !number of emission lines and continuum emission points
  INTEGER, PARAMETER :: nemline=108, nlam_nebcont=1963
  !number of metallicity, age, and ionization parameter points
  INTEGER, PARAMETER :: nebnz=7, nebnage=12, nebnip=7

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
  !mass limits for BH and neutron star initial-mass mass relations
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
  !Solar radius in cm
  REAL(SP), PARAMETER :: rsun    = 6.955E10
  !Solar luminosity in erg/s
  REAL(SP), PARAMETER :: lsun    = 3.839E33
  !Newton's constant
  REAL(SP), PARAMETER :: newton  = 6.67428E-8
  !cm in a pc
  REAL(SP), PARAMETER :: pc2cm   = 3.08568E18
  !seconds per year
  REAL(SP), PARAMETER :: yr2sc   = 3.15569E7
  !Planck's constant
  REAL(SP), PARAMETER :: hplank  = 6.6261E-27
  !constant to convert mags into propert units (see getmags.f90)
  REAL(SP), PARAMETER :: mag2cgs = LOG10(lsun/4.0/mypi/(pc2cm*pc2cm)/100.0)

  !define large and small numbers.  numbers whose abs values
  !are less than tiny_number are treated as equal to 0.0
  REAL(SP), PARAMETER :: huge_number = 10**(70.d0)
  REAL(SP), PARAMETER :: tiny_number = 10**(-70.d0)
  REAL(SP), PARAMETER :: tiny30      = 10**(-30.0)
  
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
  REAL(SP), DIMENSION(3,100) :: imf_user_alpha=0.

  !environment variable for SPS home directory
  CHARACTER(250) :: SPS_HOME=''
  !name of the filter file, if blank it defaults to allfilters.dat
  CHARACTER(30)  :: alt_filter_file=''

  !Age of Universe in Gyr (set in sps_setup.f90)
  REAL(SP) :: tuniv=0.

  !index in the wavelength array where lambda=5000A, Ly_lim
  INTEGER :: whlam5000,whlylim
  
  !this specifies the size of the full time grid
  INTEGER, PARAMETER :: ntfull = time_res_incr*nt

  !array of index definitions
  REAL(SP), DIMENSION(7,nindx) :: indexdefined=0.

  !array holding MW extinction curve indices
  INTEGER, DIMENSION(6) :: mwdindex=0

  !array holding Witt & Gordon dust models
  !wgdust(lam,tau,model,homo/clump)
  REAL(SP), DIMENSION(nspec,18,6,2) :: wgdust=0.

  !Index for P(Z) distribution.  1=closed box;
  !P(Z) = z^zpow*exp(-z/pmetals)  (see pz_convol.f90)
  !pmetals set in PARAMS structure
  REAL(SP) :: zpow=1.0

  !array holding redshift-age relation
  REAL(SP), DIMENSION(500,2) :: zagespl=0.
 
  !array holding tabulated SFH 
  REAL(SP), DIMENSION(3,ntabmax) :: sfh_tab=0.
  INTEGER :: ntabsfh=0

  !array of bandpass filters
  REAL(SP), DIMENSION(nspec,nbands) :: bands
  !magnitude of the Sun in all filters
  REAL(SP), DIMENSION(nbands) :: magsun,magvega,filter_leff
  !Vega-like star spectrum for Vega magnitude zero-point
  !spectrum of Sun, for absolute mags of Sun
  REAL(SP), DIMENSION(nspec)  :: vega_spec=0.,sun_spec=0.
  !common wavelength and frequench arrays
  REAL(SP), DIMENSION(nspec)  :: spec_lambda=0.,spec_nu=0.0

  !arrays for stellar spectral information in HR diagram
  REAL(SP), DIMENSION(ndim_logt) :: speclib_logt=0.
  REAL(SP), DIMENSION(ndim_logg) :: speclib_logg=0.
  REAL(KIND(1.0)), DIMENSION(nspec,nz,ndim_logt,ndim_logg) :: speclib=0.
  
  !AGB library (Lancon & Mouhcine 2002)
  REAL(SP), DIMENSION(nspec,n_agb_o) :: agb_spec_o=0.
  REAL(SP), DIMENSION(nz,n_agb_o)    :: agb_logt_o=0.
  REAL(SP), DIMENSION(nspec,n_agb_c) :: agb_spec_c=0.
  REAL(SP), DIMENSION(n_agb_c)       :: agb_logt_c=0.
  
  !post-AGB library (Rauch 2003)
  REAL(SP), DIMENSION(nspec,ndim_pagb,2) :: pagb_spec=0.
  REAL(SP), DIMENSION(ndim_pagb)         :: pagb_logt=0.

  !WR library (Smith et al. 2002)
  REAL(SP), DIMENSION(nspec,ndim_wr,nz) :: wrn_spec=0.,wrc_spec=0.
  REAL(SP), DIMENSION(ndim_wr)          :: wrn_logt=0.,wrc_logt=0.

  !dust emission model (Draine & Li 2007)
  REAL(SP), DIMENSION(ndim_dl07)              :: lambda_dl07=0.
  REAL(SP), DIMENSION(ndim_dl07,numin_dl07*2) :: dustem_dl07=0.
  REAL(SP), DIMENSION(nspec,7,numin_dl07*2)   :: dustem2_dl07=0.

  !circumstellar AGB dust model (Villaume et al. in prep)
  REAL(SP), DIMENSION(nspec,2,nteff_dagb,ntau_dagb) :: flux_dagb=0.
  REAL(SP), DIMENSION(2,ntau_dagb)                  :: tau1_dagb=0.
  REAL(SP), DIMENSION(2,nteff_dagb)                 :: teff_dagb=0.

  !nebular emission model
  REAL(SP), DIMENSION(nemline) :: nebem_line_pos=0.
  REAL(SP), DIMENSION(nemline,nebnz,nebnage,nebnip) :: nebem_line=0.
  REAL(SP), DIMENSION(nspec,nebnz,nebnage,nebnip) :: nebem_cont=0.
  REAL(SP), DIMENSION(nebnz)   :: nebem_logz=0.
  REAL(SP), DIMENSION(nebnage) :: nebem_age=0.
  REAL(SP), DIMENSION(nebnip)  :: nebem_logu=0.

  !arrays for the isochrone data
  REAL(SP), DIMENSION(nz,nt,nm) :: mact_isoc=0.,logl_isoc=0.,&
       logt_isoc=0.,logg_isoc=0.,ffco_isoc=0.,phase_isoc=0.
  REAL(SP), DIMENSION(nz,nt,nm) :: mini_isoc=0.

  !arrays holding the number of mass elements for each isochrone,
  !the age of each isochrone, and the metallicity of each isochrone
  INTEGER, DIMENSION(nz,nt)  :: nmass_isoc=0
  REAL(SP), DIMENSION(nz,nt) :: timestep_isoc=0.
  REAL(SP), DIMENSION(nz)    :: zlegend=-99.
  
  !arrays for the full Z-dep SSP spectra
  REAL(SP), DIMENSION(nspec,ntfull,nz) :: spec_ssp_zz=0.
  REAL(SP), DIMENSION(ntfull,nz)       :: mass_ssp_zz=0.,lbol_ssp_zz=0.
  
  REAL(SP), DIMENSION(ntfull) :: time_full=0.

  !------------Define TYPE structures-------------!
  
  !structure for the set of parameters necessary to generate a model
  TYPE PARAMS
     REAL(SP) :: pagb=1.0,dell=0.,delt=0.,fbhb=0.,sbss=0.,tau=1.0,&
          const=0.,tage=0.,fburst=0.,tburst=11.0,dust1=0.,dust2=0.,&
          logzsol=0.,zred=0.,pmetals=0.02,imf1=1.3,imf2=2.3,imf3=2.3,&
          vdmc=0.08,dust_clumps=-99.,frac_nodust=0.,dust_index=-0.7,&
          dust_tesc=7.0,frac_obrun=0.,uvb=1.0,mwr=3.1,redgb=1.0,&
          dust1_index=-1.0,mdave=0.5,sf_start=0.,sf_trunc=0.,sf_theta=0.,&
          duste_gamma=0.01,duste_umin=1.0,duste_qpah=3.5,fcstar=1.0,&
          masscut=150.0,sigma_smooth=0.,agb_dust=1.0,min_wave_smooth=1E3,&
          max_wave_smooth=1E4,gas_logu=-2.0,gas_logz=0.,igm_factor=1.0
     INTEGER :: zmet=1,sfh=0,wgp1=1,wgp2=1,wgp3=1,evtype=-1
     INTEGER, DIMENSION(nbands) :: mag_compute=1
     CHARACTER(50) :: imf_filename='', sfh_filename=''
  END TYPE PARAMS
  
  !structure for the output of the compsp routine
  TYPE COMPSPOUT
     REAL(SP) :: age=0.,mass_csp=0.,lbol_csp=0.,sfr=0.,mdust=0.
     REAL(SP), DIMENSION(nbands) :: mags=0.
     REAL(SP), DIMENSION(nspec)  :: spec=0.
     REAL(SP), DIMENSION(nindx)  :: indx=0.
  END TYPE COMPSPOUT
  
  !-----the following structures are not used in the public code-----!
  !--they are included here because some users of FSPS utilize them--!

  !structure for observational data
  TYPE OBSDAT
     REAL(SP)                    :: zred=0.,logsmass=0.
     REAL(SP), DIMENSION(nbands) :: mags=0.,magerr=0.
     REAL(SP), DIMENSION(nspec)  :: spec=0.,specerr=99.
  END TYPE OBSDAT

  !structure for using P(z) in chi2
  INTEGER, PARAMETER :: npzphot   = 200
  TYPE TPZPHOT
     REAL(SP), DIMENSION(npzphot) :: zz=0.,pz=0.
  END TYPE TPZPHOT

  !used for Powell minimization
  TYPE(OBSDAT) :: powell_data, sedfit_data

END MODULE SPS_VARS
