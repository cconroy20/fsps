 PROGRAM SIMPLE

  !set up modules
  USE sps_vars; USE sps_utils  
  IMPLICIT NONE

  !NB: the various structure types are defined in sps_vars.f90
  !    variables not explicitly defined here are defined in sps_vars.f90

  INTEGER :: i
  !define variable for SSP spectrum
  REAL(SP), DIMENSION(ntfull,nspec)  :: spec_ssp
  !define variables for Mass and Lbol info
  REAL(SP), DIMENSION(ntfull)    :: mass_ssp,lbol_ssp
  CHARACTER(100) :: file1='', file2=''
  !structure containing all necessary parameters
  TYPE(PARAMS) :: pset
  !define structure for CSP spectrum
  TYPE(COMPSPOUT), DIMENSION(ntfull) :: ocompsp
  REAL(SP) :: ssfr6,ssfr7,ssfr8,ave_age

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
  
  ! Lets compute an SSP, solar metallicity, with a Chabrier IMF
  ! with no dust, and the 'default' assumptions regarding the 
  ! locations of the isochrones

  imf_type  = 0             !define the IMF (1=Chabrier 2003)
                            !see sps_vars.f90 for details of this var
  pset%zmet = 10            !define the metallicity (see the manual)
                            !20 = solar metallacity

  CALL SPS_SETUP(pset%zmet) !read in the isochrones and spectral libraries

  !define the parameter set.  These are the default values, specified 
  !in sps_vars.f90, but are explicitly included here for transparency
  pset%sfh   = 0     !set SFH to "SSP"
  pset%const = 1.
  pset%zred  = 0.0   !redshift  
  pset%dust1 = 0.0   !dust parameter 1
  pset%dust2 = 0.0   !dust parameter 2

  pset%dell  = 0.0   !shift in log(L) for TP-AGB stars
  pset%delt  = 0.0   !shift in log(Teff) for TP-AGB stars
  pset%fbhb  = 0.0   !fraction of blue HB stars
  pset%sbss  = 0.0   !specific frequency of BS stars

  !compute the SSP
  CALL SSP_GEN(pset,mass_ssp,lbol_ssp,spec_ssp)
  !compute mags and write out mags and spec for SSP
  file1 = 'SSP_BPASS.out'
  CALL COMPSP(3,1,file1,mass_ssp,lbol_ssp,spec_ssp,pset,ocompsp)



  ! Now lets compute a 1 Gyr tau model SFH with a Salpeter IMF
  ! with a simple dust model, at a particular time 
  ! (rather than outputing all the time info)

  imf_type  = 0                !define the IMF (0=Salpeter)
                               !see sps_vars.f90 for details of this var

  !NB: you only need to re-run SPS_SETUP if you have changed the metallicity
  !    or, even better (but slower), you can call SPS_SETUP(-1) and this will
  !    set up all the metallicities at once.
  CALL SPS_SETUP(pset%zmet)    !read in the isochrones and spectral libraries

  !define the parameter set. 
  pset%sfh   = 1     !set SFH to "CSP"; sfh=1 means a normal tau model
  pset%tau   = 2.0   !tau units are Gyr
  pset%dust1 = 1.0   !dust parameter 1
  pset%dust2 = 0.3   !dust parameter 2
  !pset%tage  = 12.5  !age at which we want the mags

  !compute the CSP
  CALL SSP_GEN(pset,mass_ssp,lbol_ssp,spec_ssp)
  !compute mags, and write out mags and spec for CSP
  file1 = 'CSP.out'
  CALL COMPSP(3,1,file1,mass_ssp,lbol_ssp,spec_ssp,pset,ocompsp)

  !compute basic SFH statistics for the last entry in the ocompsp array
  !results are returned in the variables ssfr6,...,ave_age
  CALL SFHSTAT(pset,ocompsp(1),ssfr6,ssfr7,ssfr8,ave_age)


END PROGRAM SIMPLE
