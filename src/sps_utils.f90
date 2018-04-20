MODULE SPS_UTILS

  INTERFACE
     SUBROUTINE ADD_AGB_DUST(weight,tspec,mact,logt,logl,logg,&
          zz,tco,lmdot)
       USE sps_vars
       REAL(SP), DIMENSION(nspec), INTENT(out) :: tspec
       REAL(SP), INTENT(in)  :: weight,mact,logt,logl,logg,zz,tco,lmdot
     END SUBROUTINE ADD_AGB_DUST
  END INTERFACE

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
     SUBROUTINE ADD_DUST(pset,csp1,csp2,specdust,mdust,ncsp1,ncsp2,nebdust)
       USE sps_vars
       REAL(SP), INTENT(out) :: mdust
       REAL(SP), DIMENSION(nspec), INTENT(in) :: csp1,csp2
       TYPE(PARAMS), INTENT(in) :: pset
       REAL(SP), DIMENSION(nspec), INTENT(out) :: specdust
       REAL(SP), DIMENSION(nemline), INTENT(in) :: ncsp1,ncsp2
       REAL(SP), DIMENSION(nemline), INTENT(out) :: nebdust
     END SUBROUTINE ADD_DUST
  END INTERFACE

  INTERFACE
     SUBROUTINE ADD_NEBULAR(pset,sspi,sspo,nebemline)
       USE sps_vars
       TYPE(PARAMS), INTENT(in) :: pset
       REAL(SP), INTENT(in), DIMENSION(nspec,ntfull)    :: sspi
       REAL(SP), INTENT(inout), DIMENSION(nspec,ntfull) :: sspo
       REAL(SP), INTENT(inout), DIMENSION(nemline,ntfull), OPTIONAL :: nebemline

     END SUBROUTINE ADD_NEBULAR
  END INTERFACE

  INTERFACE
     SUBROUTINE ADD_REMNANTS(mass,maxmass)
       USE sps_vars
       REAL(SP), INTENT(inout) :: mass
       REAL(SP), INTENT(in) :: maxmass
     END SUBROUTINE ADD_REMNANTS
  END INTERFACE

  INTERFACE
     FUNCTION AGN_DUST(lam,spec,pset,lbol_csp)
       USE sps_vars
       REAL(SP), DIMENSION(nspec), INTENT(in) :: lam,spec
       REAL(SP), INTENT(in) :: lbol_csp
       TYPE(PARAMS), INTENT(in) :: pset
       REAL(SP), DIMENSION(nspec) :: agn_dust
     END FUNCTION AGN_DUST
  END INTERFACE

  INTERFACE
     FUNCTION AIRTOVAC(lam)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(in) :: lam
       REAL(SP), DIMENSION(SIZE(lam)) :: airtovac
     END FUNCTION AIRTOVAC
  END INTERFACE

 INTERFACE
     FUNCTION ATTN_CURVE(lambda,dtype,pset)
       USE sps_vars
       INTEGER, INTENT(in) :: dtype
       REAL(SP), INTENT(in), DIMENSION(nspec) :: lambda
       TYPE(PARAMS), INTENT(in) :: pset
       REAL(SP), DIMENSION(nspec) :: attn_curve
     END FUNCTION ATTN_CURVE
  END INTERFACE

  INTERFACE
     SUBROUTINE COMPSP(write_compsp,nzin,outfile,mass_ssp,&
          lbol_ssp,spec_ssp,pset,ocompsp)
       USE sps_vars
       INTEGER, INTENT(in) :: write_compsp,nzin
       REAL(SP), INTENT(in), DIMENSION(ntfull,nzin) :: lbol_ssp,mass_ssp
       REAL(SP), INTENT(in), DIMENSION(nspec,ntfull,nzin) :: spec_ssp
       CHARACTER(100), INTENT(in) :: outfile
       TYPE(PARAMS), INTENT(in)   :: pset
       TYPE(COMPSPOUT), INTENT(inout), DIMENSION(ntfull) :: ocompsp
     END SUBROUTINE COMPSP
  END INTERFACE

  INTERFACE
     SUBROUTINE COMPSP_GRID(pset,nti,specout)
       USE sps_vars
       TYPE(PARAMS), INTENT(in) :: pset
       INTEGER, INTENT(in) :: nti
       REAL(SP), DIMENSION(nspec), INTENT(inout) :: specout
     END SUBROUTINE COMPSP_GRID
  END INTERFACE

  INTERFACE
     SUBROUTINE CSP_GEN(mass_ssp, lbol_ssp, spec_ssp, pset, tage, nzin,&
                        mass_csp, lbol_csp, spec_csp, mdust_csp,emlin_ssp,emlin_csp)
       USE sps_vars
       REAL(SP), DIMENSION(ntfull), INTENT(in) :: mass_ssp, lbol_ssp
       REAL(SP), DIMENSION(nspec, ntfull), INTENT(in) :: spec_ssp
       TYPE(PARAMS), intent(in) :: pset
       REAL(SP), INTENT(in)  :: tage
       INTEGER, INTENT(IN) :: nzin
       REAL(SP), INTENT(out) :: mass_csp, lbol_csp, mdust_csp
       REAL(SP), INTENT(out), DIMENSION(nspec) :: spec_csp
       REAL(SP), DIMENSION(nemline, ntfull, nzin), intent(in) :: emlin_ssp
       REAL(SP), DIMENSION(nemline), intent(out) :: emlin_csp
     END SUBROUTINE CSP_GEN
  END INTERFACE

  INTERFACE
     FUNCTION FUNCINT(func,a,b)
       USE sps_vars
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP) :: funcint
       INTERFACE
          FUNCTION func(x)
            USE sps_vars
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(SIZE(x)) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION FUNCINT
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
       REAL(SP), INTENT(inout), DIMENSION(nindx) :: indices
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
     FUNCTION GET_LUMDIST(z)
       USE sps_vars
       REAL(SP), INTENT(in) :: z
       REAL(SP) :: get_lumdist
     END FUNCTION GET_LUMDIST
  END INTERFACE
  
  INTERFACE
     SUBROUTINE GETMAGS(zred,spec,mags,mag_compute)
       USE sps_vars
       REAL(SP), INTENT(in) :: zred
       REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec
       REAL(SP), DIMENSION(nbands) :: mags
       INTEGER, DIMENSION(nbands), INTENT(in), OPTIONAL  :: mag_compute
     END SUBROUTINE GETMAGS
  END INTERFACE
  
  INTERFACE
     SUBROUTINE GETSPEC(pset,mact,logt,lbol,logg,phase,ffco,lmdot,wght,spec)
       USE sps_vars
       REAL(SP), INTENT(in) :: mact,logt,lbol,logg,phase,ffco,wght,lmdot
       TYPE(PARAMS), INTENT(in) :: pset
       REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec 
     END SUBROUTINE GETSPEC
  END INTERFACE

  INTERFACE
     FUNCTION IGM_ABSORB(lam,spec,zz,factor)
       USE sps_vars
       REAL(SP), DIMENSION(nspec), INTENT(in) :: lam,spec
       REAL(SP), INTENT(in) :: zz,factor
       REAL(SP), DIMENSION(nspec) :: igm_absorb
     END FUNCTION IGM_ABSORB
  END INTERFACE

  INTERFACE
     FUNCTION INTIND(lam,func,lo,hi)
       USE sps_vars
       REAL(SP), INTENT(in), DIMENSION(nspec) :: lam,func
       REAL(SP), INTENT(in) :: lo,hi
       REAL(SP) :: intind
     END FUNCTION INTIND
  END INTERFACE

  INTERFACE
     FUNCTION INTSFWGHT(sspind, logt, sfh)
       USE sps_vars
       TYPE(SFHPARAMS), INTENT(in) :: sfh
       INTEGER, INTENT(in) :: sspind
       REAL(SP), DIMENSION(2), INTENT(in) :: logt
       REAL(SP) :: intsfwght
     END FUNCTION INTSFWGHT
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
     FUNCTION LINTERP(xin,yin,xout)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(in) :: xin,yin
       REAL(SP), INTENT(in)  :: xout
       REAL(SP) :: linterp
     END FUNCTION LINTERP
  END INTERFACE

  INTERFACE
     FUNCTION LINTERPARR(xin,yin,xout)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(in) :: xin,yin
       REAL(SP), INTENT(in), DIMENSION(:) :: xout
       REAL(SP), DIMENSION(SIZE(xout)) :: linterparr
     END FUNCTION LINTERPARR
  END INTERFACE

  INTERFACE
     FUNCTION LOCATE(xx,x)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(IN) :: xx
       REAL(SP), INTENT(IN) :: x
       INTEGER :: locate
     END FUNCTION LOCATE
  END INTERFACE

  INTERFACE
     SUBROUTINE MOD_GB(zz,t,age,delt,dell,pagb,redgb,agb,&
          nn,logl,logt,phase,wght)
       USE sps_vars
       INTEGER,  INTENT(in) :: t, nn,zz
       REAL(SP), INTENT(inout), DIMENSION(nt,nm) :: logl,logt
       REAL(SP), INTENT(in), DIMENSION(nt,nm)    :: phase
       REAL(SP), INTENT(inout), DIMENSION(nm)    :: wght
       REAL(SP), INTENT(in) :: delt, dell, pagb,redgb, agb
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
     SUBROUTINE SBF(pset,outfile)
       USE sps_vars
       CHARACTER(100), INTENT(in) :: outfile
       TYPE(PARAMS), INTENT(in)    :: pset
     END SUBROUTINE SBF
  END INTERFACE

  INTERFACE
     SUBROUTINE SETUP_TABULAR_SFH(pset, nzin)
       USE sps_vars
       TYPE(PARAMS), INTENT(in) :: pset
       INTEGER, INTENT(in) :: nzin
     END SUBROUTINE SETUP_TABULAR_SFH
  END INTERFACE

  INTERFACE
     SUBROUTINE SFHINFO(pset, age, mfrac, sfr, frac_linear)
       USE sps_vars
       TYPE(PARAMS), INTENT(in) :: pset
       REAL(SP), INTENT(in) :: age
       REAL(SP), INTENT(out) :: mfrac, sfr, frac_linear
     END SUBROUTINE SFHINFO
  END INTERFACE

  INTERFACE
     FUNCTION SFHLIMIT(tlim, sfh)
       USE sps_vars
       TYPE(SFHPARAMS), INTENT(in) :: sfh
       REAL(SP), INTENT(in) :: tlim
       REAL(SP) :: sfhlimit
     END FUNCTION SFHLIMIT
  END INTERFACE

  INTERFACE
     FUNCTION SFH_WEIGHT(sfh, imin, imax)
       USE sps_vars
       TYPE(SFHPARAMS), INTENT(in) :: sfh
       INTEGER, INTENT(in) :: imin, imax
       REAL(SP), DIMENSION(ntfull) :: sfh_weight
     END FUNCTION SFH_WEIGHT
  END INTERFACE

  INTERFACE
     FUNCTION TSUM(xin,yin)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(in) :: xin,yin
       REAL(SP) :: tsum
     END FUNCTION TSUM
  END INTERFACE

  INTERFACE
     SUBROUTINE SMOOTHSPEC(lambda,spec,sigma,minl,maxl,ires)
       USE sps_vars
       REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec
       REAL(SP), INTENT(in), DIMENSION(nspec) :: lambda
       REAL(SP), INTENT(in), DIMENSION(nspec), OPTIONAL :: ires
       REAL(SP), INTENT(in) :: sigma,minl,maxl
     END SUBROUTINE SMOOTHSPEC
  END INTERFACE

  INTERFACE
     FUNCTION VACTOAIR(lam)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(in) :: lam
       REAL(SP), DIMENSION(SIZE(lam)) :: vactoair
     END FUNCTION VACTOAIR
  END INTERFACE

  INTERFACE
     SUBROUTINE WRITE_ISOCHRONE(outfile,pset)
       USE sps_vars
       TYPE(PARAMS), INTENT(in) :: pset
       CHARACTER(100), INTENT(in)  :: outfile
     END SUBROUTINE WRITE_ISOCHRONE
  END INTERFACE

  INTERFACE
     SUBROUTINE ZTINTERP(zpos,spec,lbol,mass,tpos,zpow)
       USE sps_vars
       REAL(SP),INTENT(in) :: zpos
       REAL(SP),INTENT(in), OPTIONAL :: tpos,zpow
       REAL(SP),INTENT(inout),DIMENSION(:) :: mass, lbol
       REAL(SP),INTENT(inout),DIMENSION(:,:) :: spec
     END SUBROUTINE ZTINTERP
  END INTERFACE

END MODULE SPS_UTILS
