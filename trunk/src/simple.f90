PROGRAM SIMPLE

  !set up modules
  USE sps_vars; USE nrtype; USE sps_utils
  
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

  !-----------------------------------------------------------!
  
  ! Lets compute an SSP, solar metallicity, with a Chabrier IMF
  ! with no dust, and the 'default' assumptions regarding the 
  ! locations of the isochrones

  imf_type  = 1                !define the IMF (1=Chabrier 2003)
                               !see sps_vars.f90 for details of this var
  pset%zmet = 20               !define the metallicity (see the lookup table)
                               !20 = solar metallacity

  CALL SPS_SETUP(pset%zmet)    !read in the isochrones and spectral libraries

  !define the parameter set.  These are the default values, specified 
  !in sps_vars.f90, but are explicitly included here for transparency
  pset%sfh   = 0     !set SFH to "SSP"
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
  file1 = 'SSP.out'
  CALL COMPSP(3,1,file1,mass_ssp,lbol_ssp,spec_ssp,pset,ocompsp)

  ! Now lets compute a 1 Gyr tau model SFH with a van Dokkum 2008 IMF,
  ! with a simple dust model, for a Universe that is 10 Gyr old

  imf_type  = 3                !define the IMF (3=van Dokkum 2003)
                               !see sps_vars.f90 for details of this var
  pset%zmet = 20               !define the metallicity (see the lookup table)
                               !20 = solar metallacity

  CALL SPS_SETUP(pset%zmet)    !read in the isochrones and spectral libraries

  !define the parameter set. 
  pset%sfh   = 1     !set SFH to "CSP"
  pset%tau   = 1.0
  pset%tage  = 0.0   !set this flag to compute ages from 0<t<tuniv
  tuniv      = 10.0  !units in Gyr
  pset%zred  = 0.0   !redshift  
  pset%dust1 = 1.0   !dust parameter 1
  pset%dust2 = 0.3   !dust parameter 2

  !compute the CSP
  CALL SSP_GEN(pset,mass_ssp,lbol_ssp,spec_ssp)
  !compute mags, and write out mags and spec for CSP
  file1 = 'CSP.out'
  CALL COMPSP(3,1,file1,mass_ssp,lbol_ssp,spec_ssp,pset,ocompsp)

END PROGRAM SIMPLE
