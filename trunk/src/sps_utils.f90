MODULE SPS_UTILS

  INTERFACE
     SUBROUTINE VELBROAD(lambda,spec,sigma)
       USE sps_vars; USE nr, ONLY : locate; USE nrtype
       REAL(SP), INTENT(inout), DIMENSION(:) :: spec,lambda
       REAL(SP), INTENT(in) :: sigma
     END SUBROUTINE VELBROAD
  END INTERFACE
  
  INTERFACE
     REAL FUNCTION INTIND(lam,func,lo,hi)
       USE sps_vars; USE nrtype; USE nr, ONLY : locate
       REAL(SP), INTENT(in), DIMENSION(:) :: lam,func
       REAL(SP), INTENT(in) :: lo,hi
     END FUNCTION INTIND
  END INTERFACE

  INTERFACE
     SUBROUTINE GETZMET(smass,pos)
       USE nrtype; USE sps_vars
       IMPLICIT NONE
       REAL(SP), INTENT(in) :: smass
       TYPE(PARAMS), INTENT(inout) :: pos
     END SUBROUTINE GETZMET
  END INTERFACE

  INTERFACE
     SUBROUTINE GETINDX(lambda,spec,indices)
       USE sps_vars; USE nrtype
       REAL(SP), INTENT(in), DIMENSION(:) :: spec,lambda
       REAL(SP), INTENT(inout), DIMENSION(nindsps) :: indices
     END SUBROUTINE GETINDX
  END INTERFACE

  INTERFACE
     SUBROUTINE FITGAL_INIT(switch,pos,powell_pos)
       USE sps_vars; USE nrtype
       IMPLICIT NONE
       INTEGER, INTENT(in) :: switch
       TYPE(PARAMS), INTENT(inout) :: pos
       REAL(SP), OPTIONAL, DIMENSION(:), INTENT(inout) :: powell_pos
     END SUBROUTINE FITGAL_INIT
  END INTERFACE

  INTERFACE
     SUBROUTINE FITFAST_INIT(switch,pos,powell_pos)
       USE sps_vars; USE nrtype
       INTEGER, INTENT(in) :: switch
       TYPE(PARAMS), INTENT(inout) :: pos
       REAL(SP), OPTIONAL, DIMENSION(:), INTENT(inout) :: powell_pos
     END SUBROUTINE FITFAST_INIT
    END INTERFACE

  INTERFACE
     SUBROUTINE ADD_REMNANTS(mass,maxmass)
       USE nrtype
       IMPLICIT NONE
       REAL(SP), INTENT(inout) :: mass
       REAL(SP), INTENT(in) :: maxmass
     END SUBROUTINE ADD_REMNANTS
  END INTERFACE

  INTERFACE
     SUBROUTINE MOD_HB(f_bhb,t,mini,mact,logl,logt,phase, &
          wght,hb_wght,nmass,hbtime)
       USE sps_vars; USE nrtype
       IMPLICIT NONE
       REAL(SP), INTENT(inout), DIMENSION(nt,nm) :: mini,mact,&
            logl,logt,phase
       REAL(SP), INTENT(inout), DIMENSION(nm) :: wght
       REAL(SP), DIMENSION(nm) :: tphase=0.0
       INTEGER, INTENT(inout), DIMENSION(nt) :: nmass
       REAL(SP), INTENT(inout) :: hb_wght
       INTEGER, INTENT(in) :: t
       REAL(SP), INTENT(in) :: f_bhb, hbtime
     END SUBROUTINE MOD_HB
  END INTERFACE

  INTERFACE
     SUBROUTINE MOD_AGB(zz,t,age,delt,dell,pagb,redgb,&
          nn,logl,logt,phase,wght)
       USE sps_vars; USE nrtype
       IMPLICIT NONE
       INTEGER,  INTENT(in) :: t, nn,zz
       REAL(SP), INTENT(inout), DIMENSION(nt,nm) :: logl,logt
       REAL(SP), INTENT(in), DIMENSION(nt,nm)    :: phase
       REAL(SP), INTENT(inout), DIMENSION(nm)    :: wght
       REAL(SP), INTENT(in) :: delt, dell, pagb,redgb
       REAL(SP), INTENT(in), DIMENSION(nt) :: age
     END SUBROUTINE MOD_AGB
  END INTERFACE

  INTERFACE
     SUBROUTINE ADD_BS(s_bs,t,mini,mact,logl,logt,phase, &
          wght,hb_wght,nmass)
       USE nrtype; USE sps_vars
       IMPLICIT NONE
       REAL(SP), INTENT(inout), DIMENSION(nt,nm) :: mini,mact,&
            logl,logt,phase
       REAL(SP), INTENT(inout), DIMENSION(nm) :: wght
       REAL(SP), INTENT(in) :: hb_wght,s_bs
       INTEGER, INTENT(in)  :: t
       INTEGER, INTENT(inout), DIMENSION(nt)  :: nmass
     END SUBROUTINE ADD_BS
  END INTERFACE

  INTERFACE
     SUBROUTINE GETSPEC(zz,mini,mact,logt,lbol,phase,ffco,spec)
       USE sps_vars; USE nrtype
       IMPLICIT NONE
       REAL(SP), INTENT(in) :: mini,mact,logt,lbol,phase,ffco
       INTEGER,  INTENT(in) :: zz
       REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec 
     END SUBROUTINE GETSPEC
  END INTERFACE

  INTERFACE
     SUBROUTINE ADD_DUST(pset,csp1,csp2,specdust)
       USE nrtype; USE sps_vars
       IMPLICIT NONE
       REAL(SP), DIMENSION(nspec), INTENT(in) :: csp1,csp2
       TYPE(PARAMS), INTENT(in) :: pset
       REAL(SP), DIMENSION(nspec), INTENT(out) :: specdust
     END SUBROUTINE ADD_DUST
  END INTERFACE

  INTERFACE
     SUBROUTINE IMF_WEIGHT(mini,wght,nmass)
       USE nrtype; USE sps_vars
       IMPLICIT NONE
       REAL(SP), INTENT(inout), DIMENSION(nm) :: wght
       REAL(SP), INTENT(in), DIMENSION(nm)    :: mini
       INTEGER, INTENT(in) :: nmass
     END SUBROUTINE IMF_WEIGHT
  END INTERFACE

  INTERFACE
     FUNCTION IMF(mass)
       USE sps_vars; USE nrtype
       IMPLICIT NONE
       REAL(SP), DIMENSION(:), INTENT(in) :: mass
       REAL(SP), DIMENSION(size(mass)) :: imf
     END FUNCTION IMF
  END INTERFACE  

  INTERFACE
     SUBROUTINE READ_SPEC(file,lambda,spec,time,mass,lbol,n_isoc)
       USE sps_vars; USE nrtype
       IMPLICIT NONE
       INTEGER, INTENT(out) :: n_isoc
       CHARACTER(60), INTENT(in) :: file
       REAL(SP), INTENT(out), DIMENSION(nspec) :: lambda
       REAL(SP), INTENT(out), DIMENSION(nt,nspec) :: spec
       REAL(SP), INTENT(out), DIMENSION(nt) :: time, mass, lbol
     END SUBROUTINE READ_SPEC
  END INTERFACE
  
  INTERFACE
     SUBROUTINE GETMAGS(zred,spec,mags) 
       USE sps_vars 
       USE nrtype; USE nrutil, ONLY : assert_eq
       USE nr, ONLY : spline, splint
       IMPLICIT NONE
       REAL(SP), INTENT(in), DIMENSION(nspec) :: spec
       REAL(SP), INTENT(in) :: zred
       REAL(SP), DIMENSION(nbands) :: mags
     END SUBROUTINE GETMAGS
  END INTERFACE
  
  INTERFACE
     REAL FUNCTION SPS_PRIORS(pos)
       USE sps_vars
       IMPLICIT NONE
       TYPE(PARAMS), INTENT(in) :: pos
     END FUNCTION SPS_PRIORS
  END INTERFACE
  
  INTERFACE
     SUBROUTINE ZINTERP(zpos,spec,lbol,mass)
       USE sps_vars; USE nrtype
       USE nrutil, ONLY : assert_eq; USE nr, ONLY : locate
       IMPLICIT NONE
       REAL(SP),INTENT(in) :: zpos
       REAL(SP),INTENT(inout),DIMENSION(nt) :: mass, lbol
       REAL(SP),INTENT(inout),DIMENSION(nt,nspec) :: spec
     END SUBROUTINE ZINTERP
  END INTERFACE

  INTERFACE
     FUNCTION SPS_CHI2(data,pos,csp)
       USE sps_vars; USE nrtype
       IMPLICIT NONE
       TYPE(OBSDAT), INTENT(in)       :: data
       TYPE(COMPSPOUT), INTENT(inout) :: csp
       TYPE(PARAMS), INTENT(in)       :: pos
       REAL(SP) :: sps_chi2
     END FUNCTION SPS_CHI2
  END INTERFACE

  INTERFACE
     FUNCTION SPS_CHI2_PHOTZ(data,pos,csp,pzphot)
       USE sps_vars; USE nrtype
       IMPLICIT NONE
       TYPE(OBSDAT), INTENT(in)       :: data
       TYPE(COMPSPOUT), INTENT(inout) :: csp
       TYPE(PARAMS), INTENT(in)       :: pos
       TYPE(TPZPHOT), INTENT(in)       :: pzphot
       REAL(SP) :: sps_chi2_photz
     END FUNCTION SPS_CHI2_PHOTZ
  END INTERFACE

  INTERFACE
     FUNCTION GET_TUNIV(z)
       USE nrtype; USE sps_vars
       IMPLICIT NONE
       REAL, INTENT(in) :: z
       REAL(SP) :: get_tuniv
     END FUNCTION GET_TUNIV
  END INTERFACE

END MODULE SPS_UTILS
