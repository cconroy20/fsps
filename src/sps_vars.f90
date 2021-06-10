MODULE SPS_VARS

  ! module to set up most arrays and variables

  IMPLICIT NONE
  SAVE

!-------set the spectral library------!
#ifndef MILES
#define MILES 1
#endif

#ifndef BASEL
#define BASEL 0
#endif

#ifndef C3K
#define C3K 0
#endif

!------set the isochrone library------!
#ifndef MIST
#define MIST 1
#endif

#ifndef PADOVA
#define PADOVA 0
#endif

#ifndef PARSEC
#define PARSEC 0
#endif

#ifndef BASTI
#define BASTI 0
#endif

#ifndef GENEVA
#define GENEVA 0
#endif

!note that in the case of BPASS the SSPs are already pre-computed
!so the spectral library, IMF, etc. is fixed in this case.  
#ifndef BPASS
#define BPASS 0
#endif

!------set the dust emission model------!
#ifndef DL07
#define DL07 1
#endif

#ifndef THEMIS
#define THEMIS 0
#endif
  
  !--------------------------------------------------------------!
  !--------------------------------------------------------------!

  !note that "SP" actually means double precision; this is a hack
  !to turn the nr routines into DP
  INTEGER, PARAMETER :: SP = KIND(1.d0)

  !------Common parameters that may be altered by the user-------!

  !setup cosmology (WMAP7).  Used only for z(t) relation.
  REAL(SP) :: om0=0.27, ol0=0.73, H0=72.

  !controls the level of output
  !0 = minimal output to screen.
  !1 = lots of output to screen.  useful for debugging.
  INTEGER, PARAMETER :: verbose=0

  !flag specifying TP-AGB normalization scheme
  !0 = default Padova 2007 isochrones
  !1 = Conroy & Gunn 2010 normalization
  !2 = Villaume, Conroy, Johnson 2015 normalization
  INTEGER :: tpagb_norm_type=2

  !turn-on time for BHB and SBS phases, time is in log(yrs)
  REAL(SP), PARAMETER :: bhb_sbs_time=9.5

  !turn on/off convolution of SSP with P(Z) (pz_convol.f90)
  !NB: pz_convol.f90 has not been tested in some time, use with caution
  INTEGER :: pzcon=0

  !the factor by which we increase the time array
  !this should no longer need to be set to anything other than 1
  INTEGER, PARAMETER :: time_res_incr=1

  !whether to interpolate the SSPs in logt (0) or t (1)
  integer :: interpolation_type = 0

  !The log of the minimum age to use when computing CSPs.  The spectrum for
  !this age is taken from the youngest available SSP.  Should be less than ~3
  real(SP) :: tiny_logt = 0.0

  !Use Aringer et al. (2009) Carbon star library if set
  !otherwise use Lancon & Wood (2002) empirical spectra
  INTEGER, PARAMETER :: cstar_aringer=1

  !turn on/off computation of light-weighted stellar ages
  !NB: currently only works with sfh=1,4 options
  INTEGER :: compute_light_ages=0

  !turn on/off the Draine & Li 2007 dust emission model
  INTEGER :: add_dust_emission=1

  !turn on/off the Nenkova et al. 2008 AGN torus dust model
  INTEGER :: add_agn_dust=1

  !turn on/off the AGB circumstellar dust model
  !see Villaume et al. (2014) for details
  INTEGER :: add_agb_dust_model=1

  !turn on/off the WR spectral library
  !if off (0), will use the main default library instead
  INTEGER :: use_wr_spectra=1

  !Use Eldridge 2017 WMBasic library for stars hotter than 25,000 K
  !or this value, whichever is larger
  real(SP) :: logt_wmb_hot = 0.0

  !turn on/off a Cloudy-based nebular emission model (cont+lines)
  !if set to 2, then the nebular emission lines are added at the SSP
  !level, which may be useful if the nebular parameters are fixed
  INTEGER :: add_neb_emission=0
  !turn on/off the nebular continuum component (automatically
  !turned off if the above is set to 0)
  INTEGER  :: add_neb_continuum=1
  !include dust in the Cloudy tables or not
  INTEGER :: cloudy_dust=0

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

  !if set, smooth the SSPs within ssp_gen by an instrumental
  !LSF that is specified in data/lsf.dat
  INTEGER :: smooth_lsf=0

  !set attenuation-law for the diffuse ISM
  !0 = power-law attenuation.  See dust_index variable below
  !1 = MW extinction law, parameterized by Cardelli et al. 1989,
  !    with a UV bump strength parameterized by uvb (see params below)
  !2 = Calzetti attenuation law
  !3 = Witt & Gordon 2000 attenuation curve models
  !4 = Kriek & Conroy (2013) attenuation model
  !5 = Gordon et al. (2003) SMC bar extinction
  !6 = Reddy et al. (2015) attenuation
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

  !flag indicating whether to use the Mdot tabulated in the isochrone
  !files (if available) for the AGB dust model.  Note: only use this
  !feature with isochrone files that include Mdot (e.g., MIST)
  INTEGER :: use_isoc_mdot=0

  !flag indicating if the Gaussians used for implementing
  !nebular emission lines should be set up on initialization
  INTEGER :: setup_nebular_gaussians=0

  !Width of Gaussian kernels for initial nebular smoothing
  !if setup_nebular_gaussians=1 (units=km/s if smooth_velocity=1)
  REAL(SP) :: nebular_smooth_init=100.

  !flag to include emission lines in the spectrum
  !if not set, the line luminosities are still computed
  INTEGER :: nebemlineinspec=1

  !------------Pre-compiler defintions------------!

  !flag indicating type of isochrones to use
  !and number of metallicities in the set
#if (BASTI)
  REAL(SP), PARAMETER :: zsol = 0.020
  CHARACTER(4), PARAMETER :: isoc_type = 'bsti'
  INTEGER, PARAMETER :: nt=94
  INTEGER, PARAMETER :: nz=10
#elif (GENEVA)
  REAL(SP), PARAMETER :: zsol = 0.020
  CHARACTER(4), PARAMETER :: isoc_type = 'gnva'
  INTEGER, PARAMETER :: nt=51
  INTEGER, PARAMETER :: nz=5
#elif (MIST)
  REAL(SP), PARAMETER :: zsol = 0.0142
  CHARACTER(4), PARAMETER :: isoc_type = 'mist'
  INTEGER, PARAMETER :: nt=107
  INTEGER, PARAMETER :: nz=12
#elif (PARSEC)
  REAL(SP), PARAMETER :: zsol = 0.01524
  CHARACTER(4), PARAMETER :: isoc_type = 'prsc'
  INTEGER, PARAMETER :: nt=93
  INTEGER, PARAMETER :: nz=15
#elif (PADOVA)
  REAL(SP), PARAMETER :: zsol = 0.019
  CHARACTER(4), PARAMETER :: isoc_type = 'pdva'
  INTEGER, PARAMETER :: nt=94
  INTEGER, PARAMETER :: nz=22
#elif (BPASS)
  REAL(SP), PARAMETER :: zsol = 0.020
  CHARACTER(4), PARAMETER :: isoc_type = 'bpss'
  INTEGER, PARAMETER :: nt=43
  INTEGER, PARAMETER :: nz=12
#endif

  !flag indicating type of spectral library to use
  !and number of elements per stellar spectrum
#if (BPASS)
  REAL(SP), PARAMETER :: zsol_spec = 0.020
  CHARACTER(5), PARAMETER :: spec_type = 'bpass'
  INTEGER, PARAMETER :: nzinit=1
  INTEGER, PARAMETER :: nspec=15000
#else
#if (MILES)
  REAL(SP), PARAMETER :: zsol_spec = 0.019
  CHARACTER(5), PARAMETER :: spec_type = 'miles'
  INTEGER, PARAMETER :: nzinit=5
  INTEGER, PARAMETER :: nspec=5994
#elif (C3K)
  REAL(SP), PARAMETER :: zsol_spec = 0.0134
  CHARACTER(11), PARAMETER :: spec_type = 'c3k_afe+0.0'
  INTEGER, PARAMETER :: nzinit=11
  INTEGER, PARAMETER :: nspec=11149
#elif (BASEL)
  REAL(SP), PARAMETER :: zsol_spec = 0.020
  CHARACTER(5), PARAMETER :: spec_type = 'basel'
  INTEGER, PARAMETER :: nzinit=6
  INTEGER, PARAMETER :: nspec=1963
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
  INTEGER, PARAMETER :: nbands=159
  !number of indices defined in allindices.dat
  INTEGER, PARAMETER :: nindx=30

  !The following parameters should never be changed
  !unless you are changing the libraries

  !max dimension of array for each isochrone
  INTEGER, PARAMETER :: nm=2000  !10000
  !max number of lines to read in
  INTEGER, PARAMETER ::  nlines=1000000
  !max number of lines in tabulated SFH, LSF
  INTEGER, PARAMETER :: ntabmax=20000
  !dimensions of BaSeL library
  INTEGER, PARAMETER :: ndim_logt=68, ndim_logg=19
  !number of O-rich, C-rich AGB spectra (and Aringer C-rich spec)
  INTEGER, PARAMETER :: n_agb_o=9, n_agb_c=5, n_agb_car=9
  !number of post-AGB spectra
  INTEGER, PARAMETER :: ndim_pagb=14
  !number of WR spectra
  INTEGER, PARAMETER :: ndim_wr=12
  !dimensions of WMBasic grid
  INTEGER, PARAMETER :: ndim_wmb_logt=11,ndim_wmb_logg=3
  !parameters for circumstellar dust models
  INTEGER, PARAMETER :: ntau_dagb=50, nteff_dagb=6
  !number of emission lines and continuum emission points
  INTEGER, PARAMETER :: nemline=128, nlam_nebcont=1963
  !number of metallicity, age, and ionization parameter points
  INTEGER, PARAMETER :: nebnz=11, nebnage=10, nebnip=7
  !number of optical depths for AGN dust models
  INTEGER, PARAMETER :: nagndust=9
  !number of spectral points in the input library
  INTEGER, PARAMETER :: nagndust_spec=125

  !------------IMF-related Constants--------------!

  !Salpeter IMF index
  REAL(SP) :: salp_ind= 2.35
  !min/max masses for the IMF
  REAL(SP) :: imf_lower_limit = 0.08, imf_upper_limit=120.
  REAL(SP) :: imf_lower_bound
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

  !array holding the Gordon et al. (2003) SMC extinction
  REAL(SP), DIMENSION(nspec) :: g03smcextn=0.
  
  !Index for P(Z) distribution.  1=closed box;
  !P(Z) = z^zpow*exp(-z/pmetals)  (see pz_convol.f90)
  !pmetals set in PARAMS structure
  REAL(SP) :: zpow2=1.0

  !array holding redshift-age-DL relations
  REAL(SP), DIMENSION(500,3) :: cosmospl=0.

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

  !arrays for the WMBasic grid
  REAL(SP), DIMENSION(ndim_wmb_logt) :: wmb_logt=0.
  REAL(SP), DIMENSION(ndim_wmb_logg) :: wmb_logg=0.
  REAL(KIND(1.0)), DIMENSION(nspec,nz,ndim_wmb_logt,ndim_wmb_logg) :: wmb_spec=0.

  !AGB library (Lancon & Mouhcine 2002)
  REAL(SP), DIMENSION(nspec,n_agb_o) :: agb_spec_o=0.
  REAL(SP), DIMENSION(nz,n_agb_o)    :: agb_logt_o=0.
  REAL(SP), DIMENSION(nspec,n_agb_c) :: agb_spec_c=0.
  REAL(SP), DIMENSION(n_agb_c)       :: agb_logt_c=0.
  !C-rich library (Aringer et al. 2009)
  REAL(SP), DIMENSION(n_agb_car)       :: agb_logt_car=0.
  REAL(SP), DIMENSION(nspec,n_agb_car) :: agb_spec_car=0.

  !post-AGB library (Rauch 2003)
  REAL(SP), DIMENSION(nspec,ndim_pagb,2) :: pagb_spec=0.
  REAL(SP), DIMENSION(ndim_pagb)         :: pagb_logt=0.

  !WR library (Smith et al. 2002)
  REAL(SP), DIMENSION(nspec,ndim_wr,nz) :: wrn_spec=0.,wrc_spec=0.
  REAL(SP), DIMENSION(ndim_wr)          :: wrn_logt=0.,wrc_logt=0.

#if (DL07)
  !dust emission model (Draine & Li 2007)
  INTEGER, PARAMETER :: ndim_dustem=1001
  INTEGER, PARAMETER :: numin_dustem=22, nqpah_dustem=7
  CHARACTER(6), PARAMETER :: str_dustem='DL07'
  REAL(SP), DIMENSION(nqpah_dustem), PARAMETER :: &
       qpaharr = (/0.47,1.12,1.77,2.50,3.19,3.90,4.58/)
  REAL(SP), DIMENSION(numin_dustem) :: uminarr = &
       (/0.1,0.15,0.2,0.3,0.4,0.5,0.7,0.8,1.0,1.2,1.5,2.0,&
       2.5,3.0,4.0,5.0,7.0,8.0,12.0,15.0,20.0,25.0/)
#elif (THEMIS)
  !dust emission model (THEMIS; Jones et al. 2013, 2017)
  INTEGER, PARAMETER :: ndim_dustem=576
  INTEGER, PARAMETER :: numin_dustem=37, nqpah_dustem=11
  CHARACTER(6), PARAMETER :: str_dustem='THEMIS'
  REAL(SP), DIMENSION(nqpah_dustem), PARAMETER :: &
       qpaharr = (/0.02,0.06,0.10,0.14,0.17,0.20,0.24,0.28,0.32,0.36,0.40/)/2.2*100
  REAL(SP), DIMENSION(numin_dustem) :: uminarr = &
       (/0.1,0.12,0.15,0.17,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,1.0,&
       1.2,1.5,1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0,&
       12.0, 15.0, 17.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 80.0/)
#endif

  REAL(SP), DIMENSION(ndim_dustem)                :: lambda_dustem=0.
  REAL(SP), DIMENSION(ndim_dustem,numin_dustem*2) :: dustem_dustem=0.
  REAL(SP), DIMENSION(nspec,nqpah_dustem,numin_dustem*2) :: dustem2_dustem=0.

  
  !circumstellar AGB dust model (Villaume et al. 2015)
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
  !minimum resolution for nebular lines, based
  !on the resolution of the spectral libraries.
  REAL(SP), DIMENSION(nspec)   :: neb_res_min=0.0
  REAL(SP), DIMENSION(nspec,nemline) :: gaussnebarr=0.0

  !arrays for AGN dust
  REAL(SP), DIMENSION(nagndust)       :: agndust_tau=0.
  REAL(SP), DIMENSION(nspec,nagndust) :: agndust_spec=0.

  !arrays for the isochrone data
  REAL(SP), DIMENSION(nz,nt,nm) :: mact_isoc=0.,logl_isoc=0.,&
       logt_isoc=0.,logg_isoc=0.,ffco_isoc=0.,phase_isoc=0.,&
       mini_isoc=0.,lmdot_isoc=0.

  !arrays holding the number of mass elements for each isochrone,
  !the age of each isochrone, and the metallicity of each isochrone
  INTEGER, DIMENSION(nz,nt)  :: nmass_isoc=0
  REAL(SP), DIMENSION(nz,nt) :: timestep_isoc=0.
  REAL(SP), DIMENSION(nz)    :: zlegend=-99.
  REAL(SP), DIMENSION(nzinit):: zlegendinit=-99.

  !arrays for the full Z-dep SSP spectra
  REAL(SP), DIMENSION(nspec,ntfull,nz) :: spec_ssp_zz=0.
  REAL(SP), DIMENSION(ntfull,nz)       :: mass_ssp_zz=0.,lbol_ssp_zz=0.

  REAL(SP), DIMENSION(ntfull) :: time_full=0.

  !array for full BPASS SSPs
  REAL(SP), DIMENSION(nspec,nt,nz) :: bpass_spec_ssp=0.
  REAL(SP), DIMENSION(nt,nz)       :: bpass_mass_ssp=0.
  
  !------------Define TYPE structures-------------!

  !structure for the set of parameters necessary to generate a model
  TYPE PARAMS
     REAL(SP) :: pagb=1.0,dell=0.,delt=0.,fbhb=0.,sbss=0.,tau=1.0,&
          const=0.,tage=0.,fburst=0.,tburst=11.0,dust1=0.,dust2=0.,&
          logzsol=0.,zred=0.,pmetals=0.02,imf1=1.3,imf2=2.3,imf3=2.3,&
          vdmc=0.08,dust_clumps=-99.,frac_nodust=0.,dust_index=-0.7,&
          dust_tesc=7.0,frac_obrun=0.,uvb=1.0,mwr=3.1,redgb=1.0,agb=1.0,&
          dust1_index=-1.0,mdave=0.5,sf_start=0.,sf_trunc=0.,sf_slope=0.,&
          duste_gamma=0.01,duste_umin=1.0,duste_qpah=3.5,fcstar=1.0,&
          masscut=150.0,sigma_smooth=0.,agb_dust=1.0,min_wave_smooth=1E3,&
          max_wave_smooth=1E4,gas_logu=-2.0,gas_logz=0.,igm_factor=1.0,&
          fagn=0.0,agn_tau=10.0
     INTEGER :: zmet=1,sfh=0,wgp1=1,wgp2=1,wgp3=1,evtype=-1
     INTEGER, DIMENSION(nbands) :: mag_compute=1
     INTEGER, DIMENSION(nt) :: ssp_gen_age=1
     CHARACTER(50) :: imf_filename='', sfh_filename=''
  END TYPE PARAMS

  !structure for the output of the compsp routine
  TYPE COMPSPOUT
     REAL(SP) :: age=0.,mass_csp=0.,lbol_csp=0.,sfr=0.,mdust=0.,mformed=0.
     REAL(SP), DIMENSION(nbands)  :: mags=0.
     REAL(SP), DIMENSION(nspec)   :: spec=0.
     REAL(SP), DIMENSION(nindx)   :: indx=0.
     REAL(SP), DIMENSION(nemline) :: emlines=0.
  END TYPE COMPSPOUT

  ! A structure to hold SFH params converted to intrinsic units
  TYPE SFHPARAMS
     REAL(SP) :: tau=1.0,tage=0.,tburst=0.,sf_trunc=0.,sf_slope=0.,&
          tq=0.,t0=0.,tb=0.
     INTEGER :: type=0,use_simha_limits=0
  END TYPE SFHPARAMS

  TYPE TLSF
     REAL(SP), DIMENSION(nspec) :: lsf=0.
     REAL(SP) :: minlam=0.,maxlam=0.
  END TYPE TLSF

  TYPE(TLSF) :: lsfinfo

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

  !used for creating a pre-tabulated grid of CSPs as a function
  !of tau and metallicity
  !INTEGER, PARAMETER :: ntaugrid=20
  !REAL(SP), DIMENSION(ntaugrid) :: taugrid=0.0
  !REAL, DIMENSION(nspec,ntfull,ntaugrid,nz) :: csp_grid=0.0
  !INTEGER :: csp_grid_flag=0


END MODULE SPS_VARS
