SUBROUTINE ADD_BS(s_bs,t,mini,mact,logl,logt,phase, &
     wght,hb_wght,nmass)

  !routine to add blue straggler stars into older 
  !stellar populations.  The BS stars are just an extension
  !of the main sequence.
  !They are given a weight relative to the horizontal branch, 
  !S_bs, a metric which is common in the literature but is 
  !otherwise not terribly physical.

  !They are also given masses appropriate for MS stars with 
  !their luminosity and Teff.  This is suggested by the recent
  !results of M67 (Liu et al. 2008), but should nontheless be 
  !considered an uncertainty.

  !Note that the parameter bhb_sbs_time, set in sps_vars.f90,
  !sets the turn-on time for this modification

  USE sps_vars; USE nr, ONLY : spline,splint
  USE nrtype; USE nrutil, ONLY : assert_eq
  IMPLICIT NONE

  REAL(SP), INTENT(inout), DIMENSION(nt,nm) :: mini,mact,&
       logl,logt,phase
  REAL(SP), INTENT(inout), DIMENSION(nm) :: wght
  REAL(SP), INTENT(in) :: hb_wght, s_bs
  INTEGER, INTENT(in)  :: t
  INTEGER, INTENT(inout), DIMENSION(nt) :: nmass

  REAL(SP), DIMENSION(nm) :: t2lspl=0., l2tspl=0., l2mspl=0.
  !number of BS stars to add (not important b/c their 
  !total weight remains fixed)
  INTEGER, PARAMETER :: nbs = 20
  !weight given to total BS population
  REAL(SP) :: bs_wght=0.
  REAL(SP) :: tol=0., msspl=0.
  INTEGER  :: maxt=0, i, k, idum=-1

  tol = 0.0
  bs_wght = s_bs * hb_wght
  
  !find the extent of the t~0 MS
  maxt = 1
  DO WHILE(logl(1,maxt).LT.3.5)
     maxt = maxt+1
  ENDDO

  !set up the Lbol, Teff, and Mass splines along the t~0 MS
  CALL spline(logt(1,1:maxt),logl(1,1:maxt),1.0e30_sp,&
       1.0e30_sp,t2lspl(1:maxt))
  CALL spline(logl(1,1:maxt),logt(1,1:maxt),1.0e30_sp,&
       1.0e30_sp,l2tspl(1:maxt))
  CALL spline(logl(1,1:maxt),mini(1,1:maxt),1.0e30_sp,&
       1.0e30_sp,l2mspl(1:maxt))

  !find the MS turn-off at the current age
  i=0
  DO WHILE (tol.LT.0.2)
     i = i+1
     msspl = splint(logt(1,1:maxt),logl(1,1:maxt),&
             t2lspl(1:maxt),logt(t,i))
     tol   = ABS(msspl-logl(t,i))
  ENDDO

  !now add the BS stars
  DO k=1,nbs

     !distribute them uniformly from L_TO to L_TO+0.75 dex
     !luminosity offset by 0.2 dex to better match NGC 5466
     logl(t,nmass(t)+k) = k*0.75/nbs + logl(t,i-1) + 0.2
     mini(t,nmass(t)+k) = splint(logl(1,1:maxt),mini(1,1:maxt),&
          l2mspl(1:maxt),logl(t,nmass(t)+k))
     mact(t,nmass(t)+k) = mini(t,nmass(t)+k)
     logt(t,nmass(t)+k) = splint(logl(1,1:maxt),logt(1,1:maxt),&
          l2tspl(1:maxt),logl(t,nmass(t)+k))
     phase(t,nmass(t)+k) = 7.
     wght(nmass(t)+k)    = 1./nbs*bs_wght
     
  ENDDO

  !update number of stars in the isochrone
  nmass(t) = nmass(t) + nbs

  RETURN

END SUBROUTINE ADD_BS
