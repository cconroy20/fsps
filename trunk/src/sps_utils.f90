MODULE SPS_UTILS

  INTERFACE
     SUBROUTINE ADD_BS(s_bs,t,mini,mact,logl,logt,logg,phase, &
          wght,hb_wght,nmass)
       USE sps_vars
       REAL(SP), INTENT(inout), DIMENSION(nt,nm) :: mini,mact,&
            logl,logt,logg,phase
       REAL(SP), INTENT(inout), DIMENSION(nm) :: wght
       REAL(SP), INTENT(in) :: hb_wght,s_bs
       INTEGER, INTENT(in)  :: t
       INTEGER, INTENT(inout), DIMENSION(nt)  :: nmass
     END SUBROUTINE ADD_BS
  END INTERFACE

  INTERFACE
     SUBROUTINE ADD_DUST(pset,csp1,csp2,specdust,mdust)
       USE sps_vars
       REAL(SP), INTENT(out) :: mdust
       REAL(SP), DIMENSION(nspec), INTENT(in) :: csp1,csp2
       TYPE(PARAMS), INTENT(in) :: pset
       REAL(SP), DIMENSION(nspec), INTENT(out) :: specdust
     END SUBROUTINE ADD_DUST
  END INTERFACE

  INTERFACE
     SUBROUTINE ADD_REMNANTS(mass,maxmass)
       USE sps_vars
       REAL(SP), INTENT(inout) :: mass
       REAL(SP), INTENT(in) :: maxmass
     END SUBROUTINE ADD_REMNANTS
  END INTERFACE

  INTERFACE
     SUBROUTINE COMPSP(write_compsp,nzin,outfile,mass_ssp,&
          lbol_ssp,spec_ssp,pset,ocompsp)
       USE sps_vars
       INTEGER, INTENT(in) :: write_compsp,nzin
       REAL(SP), INTENT(in), DIMENSION(nzin,ntfull) :: lbol_ssp,mass_ssp
       REAL(SP), INTENT(in), DIMENSION(nzin,ntfull,nspec) :: spec_ssp
       CHARACTER(100), INTENT(in) :: outfile
       TYPE(PARAMS), INTENT(in)   :: pset
       TYPE(COMPSPOUT), INTENT(inout), DIMENSION(ntfull) :: ocompsp
     END SUBROUTINE COMPSP
  END INTERFACE

  INTERFACE
     SUBROUTINE GETZMET(smass,pos)
       USE sps_vars
       REAL(SP), INTENT(in) :: smass
       TYPE(PARAMS), INTENT(inout) :: pos
     END SUBROUTINE GETZMET
  END INTERFACE

  INTERFACE
     SUBROUTINE GETINDX(lambda,spec,indices)
       USE sps_vars
       REAL(SP), INTENT(in), DIMENSION(nspec) :: spec,lambda
       REAL(SP), INTENT(inout), DIMENSION(nindsps) :: indices
     END SUBROUTINE GETINDX
  END INTERFACE

  !this is a private routine not included in the public release
  INTERFACE
     SUBROUTINE FITGAL_INIT(switch,pos,powell_pos)
       USE sps_vars
       INTEGER, INTENT(in) :: switch
       TYPE(PARAMS), INTENT(inout) :: pos
       REAL(SP), OPTIONAL, DIMENSION(:), INTENT(inout) :: powell_pos
     END SUBROUTINE FITGAL_INIT
  END INTERFACE

  !this is a private routine not included in the public release
  INTERFACE
     SUBROUTINE FITFAST_INIT(switch,pos,powell_pos)
       USE sps_vars
       INTEGER, INTENT(in) :: switch
       TYPE(PARAMS), INTENT(inout) :: pos
       REAL(SP), OPTIONAL, DIMENSION(:), INTENT(inout) :: powell_pos
     END SUBROUTINE FITFAST_INIT
    END INTERFACE

  INTERFACE
     FUNCTION GET_TUNIV(z)
       USE sps_vars
       REAL(SP), INTENT(in) :: z
       REAL(SP) :: get_tuniv
     END FUNCTION GET_TUNIV
  END INTERFACE
  
  INTERFACE
     FUNCTION LINTERP(xin,yin,xout)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(in) :: xin,yin
       REAL(SP), INTENT(in)  :: xout
       REAL(SP) :: linterp
     END FUNCTION LINTERP
  END INTERFACE

  INTERFACE
     SUBROUTINE GETMAGS(zred,spec,mags) 
       USE sps_vars
       REAL(SP), INTENT(in), DIMENSION(nspec) :: spec
       REAL(SP), INTENT(in) :: zred
       REAL(SP), DIMENSION(nbands) :: mags
     END SUBROUTINE GETMAGS
  END INTERFACE
  
  INTERFACE
     SUBROUTINE GETSPEC(zz,mact,logt,lbol,logg,phase,ffco,spec)
       USE sps_vars
       REAL(SP), INTENT(in) :: mact,logt,lbol,logg,phase,ffco
       INTEGER,  INTENT(in) :: zz
       REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec 
     END SUBROUTINE GETSPEC
  END INTERFACE

  INTERFACE
     FUNCTION INTIND(lam,func,lo,hi)
       USE sps_vars; USE nr, ONLY : locate
       REAL(SP), INTENT(in), DIMENSION(nspec) :: lam,func
       REAL(SP), INTENT(in) :: lo,hi
       REAL(SP) :: intind
     END FUNCTION INTIND
  END INTERFACE

  INTERFACE
     FUNCTION INTSFR(sfh,tau,const,maxtime,sfstart,t1,t2,tweight)
       USE sps_vars
       INTEGER, INTENT(in)  :: sfh
       REAL(SP), INTENT(in) :: t1,t2,tau,const,maxtime,sfstart
       REAL(SP) :: intsfr
       INTEGER, intent(in), optional :: tweight
     END FUNCTION INTSFR
  END INTERFACE

  INTERFACE
     SUBROUTINE INTSPEC(pset,nti,spec_ssp,csp,mass_ssp,lbol_ssp,&
          mass,lbol,specb,massb,lbolb,deltb,sfstart,tau,const,maxtime,&
          mdust,tweight)
       USE sps_vars
       INTEGER, intent(in), optional :: tweight
       INTEGER,  INTENT(in)    :: nti
       REAL(SP), INTENT(in)    :: massb,lbolb,deltb,sfstart,tau,const,maxtime
       REAL(SP), INTENT(inout) :: mass, lbol, mdust
       REAL(SP), INTENT(in), DIMENSION(ntfull) :: mass_ssp,lbol_ssp
       REAL(SP), INTENT(in), DIMENSION(ntfull,nspec) :: spec_ssp
       REAL(SP), INTENT(in), DIMENSION(nspec)    :: specb
       REAL(SP), INTENT(inout), DIMENSION(nspec) :: csp
       REAL(SP), DIMENSION(nspec)  :: csp1,csp2
       TYPE(PARAMS), INTENT(in)    :: pset
     END SUBROUTINE INTSPEC
  END INTERFACE

  INTERFACE
     FUNCTION IMF(mass)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(in) :: mass
       REAL(SP), DIMENSION(size(mass)) :: imf
     END FUNCTION IMF
  END INTERFACE 

  INTERFACE
     SUBROUTINE IMF_WEIGHT(mini,wght,nmass)
       USE sps_vars
       REAL(SP), INTENT(inout), DIMENSION(nm) :: wght
       REAL(SP), INTENT(in), DIMENSION(nm)    :: mini
       INTEGER, INTENT(in) :: nmass
     END SUBROUTINE IMF_WEIGHT
  END INTERFACE 

  INTERFACE
     SUBROUTINE MOD_GB(zz,t,age,delt,dell,pagb,redgb,&
          nn,logl,logt,phase,wght)
       USE sps_vars
       INTEGER,  INTENT(in) :: t, nn,zz
       REAL(SP), INTENT(inout), DIMENSION(nt,nm) :: logl,logt
       REAL(SP), INTENT(in), DIMENSION(nt,nm)    :: phase
       REAL(SP), INTENT(inout), DIMENSION(nm)    :: wght
       REAL(SP), INTENT(in) :: delt, dell, pagb,redgb
       REAL(SP), INTENT(in), DIMENSION(nt) :: age
     END SUBROUTINE MOD_GB
  END INTERFACE

  INTERFACE
     SUBROUTINE MOD_HB(f_bhb,t,mini,mact,logl,logt,logg,phase, &
          wght,hb_wght,nmass,hbtime)
       USE sps_vars
       REAL(SP), INTENT(inout), DIMENSION(nt,nm) :: mini,mact,&
            logl,logt,logg,phase
       REAL(SP), INTENT(inout), DIMENSION(nm) :: wght
       REAL(SP), DIMENSION(nm) :: tphase=0.0
       INTEGER, INTENT(inout), DIMENSION(nt) :: nmass
       REAL(SP), INTENT(inout) :: hb_wght
       INTEGER, INTENT(in) :: t
       REAL(SP), INTENT(in) :: f_bhb, hbtime
     END SUBROUTINE MOD_HB
  END INTERFACE

  !this is a private routine not included in the public release
  INTERFACE
     SUBROUTINE READ_SPEC(file,lambda,spec,time,mass,lbol,n_isoc)
       USE sps_vars
       INTEGER, INTENT(out) :: n_isoc
       CHARACTER(60), INTENT(in) :: file
       REAL(SP), INTENT(out), DIMENSION(nspec) :: lambda
       REAL(SP), INTENT(out), DIMENSION(nt,nspec) :: spec
       REAL(SP), INTENT(out), DIMENSION(nt) :: time, mass, lbol
     END SUBROUTINE READ_SPEC
  END INTERFACE
  
  !this is a private routine not included in the public release
  INTERFACE
     FUNCTION SPS_PRIORS(pos)
       USE sps_vars
       TYPE(PARAMS), INTENT(in) :: pos
       REAL(SP) :: sps_priors
     END FUNCTION SPS_PRIORS
  END INTERFACE
  
  !this is a private routine not included in the public release
  INTERFACE
     FUNCTION SPS_CHI2(data,pos,csp,c)
       USE sps_vars
       TYPE(OBSDAT), INTENT(in)       :: data
       TYPE(COMPSPOUT), INTENT(inout) :: csp
       TYPE(PARAMS), INTENT(in)       :: pos
       REAL(SP), INTENT(out), OPTIONAL :: c
       REAL(SP) :: sps_chi2
     END FUNCTION SPS_CHI2
  END INTERFACE

  !this is a private routine not included in the public release
  INTERFACE
     FUNCTION SPS_CHI2_PHOTZ(data,pos,csp,pzphot)
       USE sps_vars
       TYPE(OBSDAT), INTENT(in)       :: data
       TYPE(COMPSPOUT), INTENT(inout) :: csp
       TYPE(PARAMS), INTENT(in)       :: pos
       TYPE(TPZPHOT), INTENT(in)       :: pzphot
       REAL(SP) :: sps_chi2_photz
     END FUNCTION SPS_CHI2_PHOTZ
  END INTERFACE

  INTERFACE
     SUBROUTINE VELBROAD(lambda,spec,sigma)
       USE sps_vars
       REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec
       REAL(SP), INTENT(in), DIMENSION(nspec) :: lambda
       REAL(SP), INTENT(in) :: sigma
     END SUBROUTINE VELBROAD
  END INTERFACE

  INTERFACE
     SUBROUTINE WRITE_ISOCHRONE(file,zz)
       USE sps_vars
       INTEGER, INTENT(in) :: zz
       CHARACTER(100), INTENT(in)  :: file
     END SUBROUTINE WRITE_ISOCHRONE
  END INTERFACE

  !this is a private routine not included in the public release
  INTERFACE
     SUBROUTINE ZINTERP(zpos,spec,lbol,mass)
       USE sps_vars
       REAL(SP),INTENT(in) :: zpos
       REAL(SP),INTENT(inout),DIMENSION(ntfull) :: mass, lbol
       REAL(SP),INTENT(inout),DIMENSION(ntfull,nspec) :: spec
     END SUBROUTINE ZINTERP
  END INTERFACE

END MODULE SPS_UTILS
