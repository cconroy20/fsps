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
       REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec
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

  INTERFACE
     FUNCTION MYRAN()
       USE sps_vars
       REAL(SP) :: myran
     END FUNCTION MYRAN
  END INTERFACE

  INTERFACE
     FUNCTION TSUM(xin,yin)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(in) :: xin,yin
       REAL(SP) :: tsum
     END FUNCTION TSUM
  END INTERFACE

  INTERFACE
     SUBROUTINE SMOOTHSPEC(lambda,spec,sigma)
       USE sps_vars
       REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec
       REAL(SP), INTENT(in), DIMENSION(nspec) :: lambda
       REAL(SP), INTENT(in) :: sigma
     END SUBROUTINE SMOOTHSPEC
  END INTERFACE

  INTERFACE
     SUBROUTINE WRITE_ISOCHRONE(file,zz)
       USE sps_vars
       INTEGER, INTENT(in) :: zz
       CHARACTER(100), INTENT(in)  :: file
     END SUBROUTINE WRITE_ISOCHRONE
  END INTERFACE

END MODULE SPS_UTILS
