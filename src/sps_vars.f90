MODULE SPS_VARS

  ! module to set up most arrays and variables

  IMPLICIT NONE
  SAVE

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

  !turn on/off the X-ray binary (ULX) model from Garofali et al. (in prep)
  INTEGER :: add_xrb_emission=0

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
  ! Replaced by runtime variables

  !flag indicating type of isochrones to use
  !and number of metallicities in the set
  REAL(SP) :: zsol = 0.0
  CHARACTER(LEN=:), ALLOCATABLE :: isoc_type
  INTEGER :: nt=0
  INTEGER :: nz=0

  !flag indicating type of spectral library to use
  !and number of elements per stellar spectrum
  REAL(SP) :: zsol_spec = 0.0
  CHARACTER(LEN=:), ALLOCATABLE :: spec_type
  INTEGER :: nzinit=0
  INTEGER :: nspec=0

  !flag indicating the type of normalization used in the BaSeL library
  !pdva = normalized to Padova isochrones
  !wlbc = normalized to Teff-color relations
  !NB: currently only the wlbc option is included in the public release
  CHARACTER(4), PARAMETER :: basel_str = 'wlbc'

  !---------Dimensions of various arrays----------!

  !You must change the number of bands here if
  !filters are added to allfilters.dat
  INTEGER :: nbands=0
  !number of indices defined in allindices.dat
  INTEGER :: nindx=0

  !The following parameters should never be changed
  !unless you are changing the libraries

  !max dimension of array for each isochrone
  INTEGER, PARAMETER :: nm=2000
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
  INTEGER, PARAMETER :: nemline=166, nlam_nebcont=1963
  !number of metallicity, age, and ionization parameter points
  INTEGER, PARAMETER :: nebnz=11, nebnage=10, nebnip=7
  !number of optical depths for AGN dust models
  INTEGER, PARAMETER :: nagndust=9
  !number of spectral points in the input library
  INTEGER, PARAMETER :: nagndust_spec=125

  INTEGER :: nspec_xrb=0
  INTEGER :: nt_xrb=0
  INTEGER :: nz_xrb=0

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
  INTEGER :: ntfull = 0

  !array of index definitions
  REAL(SP), ALLOCATABLE :: indexdefined(:,:)

  !array holding MW extinction curve indices
  INTEGER, DIMENSION(6) :: mwdindex=0

  !array holding Witt & Gordon dust models
  !wgdust(lam,tau,model,homo/clump)
  REAL(SP), ALLOCATABLE :: wgdust(:,:,:,:)

  !array holding the Gordon et al. (2003) SMC extinction
  REAL(SP), ALLOCATABLE :: g03smcextn(:)

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
  REAL(SP), ALLOCATABLE :: bands(:,:)
  !magnitude of the Sun in all filters
  REAL(SP), ALLOCATABLE :: magsun(:),magvega(:),filter_leff(:)
  !Vega-like star spectrum for Vega magnitude zero-point
  !spectrum of Sun, for absolute mags of Sun
  REAL(SP), ALLOCATABLE  :: vega_spec(:),sun_spec(:)
  !common wavelength and frequench arrays
  REAL(SP), ALLOCATABLE  :: spec_lambda(:),spec_nu(:)
  !common wavelength and frequency arrays for dummy resolution files
  REAL(SP), ALLOCATABLE  :: spec_res(:)

  !arrays for stellar spectral information in HR diagram
  REAL(SP), DIMENSION(ndim_logt) :: speclib_logt=0.
  REAL(SP), DIMENSION(ndim_logg) :: speclib_logg=0.
  REAL(KIND(1.0)), ALLOCATABLE :: speclib(:,:,:,:)

  !arrays for the WMBasic grid
  REAL(SP), DIMENSION(ndim_wmb_logt) :: wmb_logt=0.
  REAL(SP), DIMENSION(ndim_wmb_logg) :: wmb_logg=0.
  REAL(KIND(1.0)), ALLOCATABLE :: wmb_spec(:,:,:,:)

  !AGB library (Lancon & Mouhcine 2002)
  REAL(SP), ALLOCATABLE :: agb_spec_o(:,:)
  REAL(SP), ALLOCATABLE :: agb_logt_o(:,:)
  REAL(SP), ALLOCATABLE :: agb_spec_c(:,:)
  REAL(SP), ALLOCATABLE :: agb_logt_c(:)
  !C-rich library (Aringer et al. 2009)
  REAL(SP), DIMENSION(n_agb_car)       :: agb_logt_car=0.
  REAL(SP), ALLOCATABLE :: agb_spec_car(:,:)

  !post-AGB library (Rauch 2003)
  REAL(SP), ALLOCATABLE :: pagb_spec(:,:,:)
  REAL(SP), DIMENSION(ndim_pagb)         :: pagb_logt=0.

  !WR library (Smith et al. 2002)
  REAL(SP), ALLOCATABLE :: wrn_spec(:,:,:),wrc_spec(:,:,:)
  REAL(SP), DIMENSION(ndim_wr)          :: wrn_logt=0.,wrc_logt=0.

  !dust emission model (Draine & Li 2007 or THEMIS)
  INTEGER :: ndim_dustem=0
  INTEGER :: numin_dustem=0, nqpah_dustem=0
  CHARACTER(6) :: str_dustem='DL07'
  
  REAL(SP), ALLOCATABLE :: qpaharr(:)
  REAL(SP), ALLOCATABLE :: uminarr(:)

  REAL(SP), ALLOCATABLE :: lambda_dustem(:)
  REAL(SP), ALLOCATABLE :: dustem_dustem(:,:)
  REAL(SP), ALLOCATABLE :: dustem2_dustem(:,:,:)

  !circumstellar AGB dust model (Villaume et al. 2015)
  REAL(SP), ALLOCATABLE :: flux_dagb(:,:,:,:)
  REAL(SP), DIMENSION(2,ntau_dagb)                  :: tau1_dagb=0.
  REAL(SP), DIMENSION(2,nteff_dagb)                 :: teff_dagb=0.

  !nebular emission model
  REAL(SP), DIMENSION(nemline) :: nebem_line_pos=0.
  REAL(SP), DIMENSION(nemline,nebnz,nebnage,nebnip) :: nebem_line=0.,xnebem_line=0.
  REAL(SP), ALLOCATABLE :: nebem_cont(:,:,:,:),xnebem_cont(:,:,:,:)
  REAL(SP), DIMENSION(nebnz)   :: nebem_logz=0.
  REAL(SP), DIMENSION(nebnage) :: nebem_age=0.
  REAL(SP), DIMENSION(nebnip)  :: nebem_logu=0.
  !minimum resolution for nebular lines, based
  !on the resolution of the spectral libraries.
  REAL(SP), ALLOCATABLE   :: neb_res_min(:)
  REAL(SP), ALLOCATABLE :: gaussnebarr(:,:)

  !arrays for AGN dust
  REAL(SP), DIMENSION(nagndust)       :: agndust_tau=0.
  REAL(SP), ALLOCATABLE :: agndust_spec(:,:)

  !arrays for the isochrone data
  REAL(SP), ALLOCATABLE :: mact_isoc(:,:,:),logl_isoc(:,:,:),&
       logt_isoc(:,:,:),logg_isoc(:,:,:),ffco_isoc(:,:,:),phase_isoc(:,:,:),&
       mini_isoc(:,:,:),lmdot_isoc(:,:,:)

  !arrays holding the number of mass elements for each isochrone,
  !the age of each isochrone, and the metallicity of each isochrone
  INTEGER, ALLOCATABLE  :: nmass_isoc(:,:)
  REAL(SP), ALLOCATABLE :: timestep_isoc(:,:)
  REAL(SP), ALLOCATABLE    :: zlegend(:)
  REAL(SP), ALLOCATABLE :: zlegendinit(:)

  !arrays for the full Z-dep SSP spectra
  REAL(SP), ALLOCATABLE :: spec_ssp_zz(:,:,:)
  REAL(SP), ALLOCATABLE :: mass_ssp_zz(:,:),lbol_ssp_zz(:,:)
  REAL(SP), ALLOCATABLE :: time_full(:)

  !array for ssp weights
  REAL(SP), ALLOCATABLE :: weight_ssp(:,:)

  !array for young and old ages
  REAL(SP), ALLOCATABLE :: spec_young(:),spec_old(:)

  !array for full BPASS SSPs
  REAL(SP), ALLOCATABLE :: bpass_spec_ssp(:,:,:)
  REAL(SP), ALLOCATABLE :: bpass_mass_ssp(:,:)

  !arrays for X-ray binaries
  REAL(SP), ALLOCATABLE :: lam_xrb(:)
  REAL(SP), ALLOCATABLE :: spec_xrb(:,:,:)
  REAL(SP), ALLOCATABLE :: ages_xrb(:)
  REAL(SP), ALLOCATABLE :: zmet_xrb(:)
  
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
          fagn=0.0,agn_tau=10.0,frac_xrb=1.0,dust3=0.
     INTEGER :: zmet=1,sfh=0,wgp1=1,wgp2=1,wgp3=1,evtype=-1
     INTEGER, ALLOCATABLE :: mag_compute(:)
     INTEGER, ALLOCATABLE :: ssp_gen_age(:)
     CHARACTER(50) :: imf_filename='', sfh_filename=''
  END TYPE PARAMS

  !structure for the output of the compsp routine
  TYPE COMPSPOUT
     REAL(SP) :: age=0.,mass_csp=0.,lbol_csp=0.,sfr=0.,mdust=0.,mformed=0.
     REAL(SP), ALLOCATABLE  :: mags(:)
     REAL(SP), ALLOCATABLE   :: spec(:)
     REAL(SP), ALLOCATABLE   :: indx(:)
     REAL(SP), ALLOCATABLE :: emlines(:)
  END TYPE COMPSPOUT

  ! A structure to hold SFH params converted to intrinsic units
  TYPE SFHPARAMS
     REAL(SP) :: tau=1.0,tage=0.,tburst=0.,sf_trunc=0.,sf_slope=0.,&
          tq=0.,t0=0.,tb=0.
     INTEGER :: type=0,use_simha_limits=0
  END TYPE SFHPARAMS

  TYPE TLSF
     REAL(SP), ALLOCATABLE :: lsf(:)
     REAL(SP) :: minlam=0.,maxlam=0.
  END TYPE TLSF

  TYPE(TLSF) :: lsfinfo

  !-----the following structures are not used in the public code-----!
  !--they are included here because some users of FSPS utilize them--!

  !structure for observational data
  TYPE OBSDAT
     REAL(SP)                    :: zred=0.,logsmass=0.
     REAL(SP), ALLOCATABLE :: mags(:),magerr(:)
     REAL(SP), ALLOCATABLE  :: spec(:),specerr(:)
  END TYPE OBSDAT

  !structure for using P(z) in chi2
  INTEGER, PARAMETER :: npzphot   = 200
  TYPE TPZPHOT
     REAL(SP), DIMENSION(npzphot) :: zz=0.,pz=0.
  END TYPE TPZPHOT

  !used for Powell minimization
  TYPE(OBSDAT) :: powell_data, sedfit_data

END MODULE SPS_VARS
