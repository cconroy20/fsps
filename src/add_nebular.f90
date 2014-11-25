SUBROUTINE ADD_NEBULAR(pset,sspi,sspo)

  !routine to add nebular emission (both line and continuum)
  !to input SSPs (sspi).  Returns SSPs as output (sspo).

  USE sps_vars; USE sps_utils, ONLY : locate,tsum
  IMPLICIT NONE

  INTEGER :: t,i,nti,a1,z1,u1
  REAL(SP) :: da,dz,du,sigma,dlam,qq,compute_q
  TYPE(PARAMS), INTENT(in) :: pset
  REAL(SP), INTENT(in), DIMENSION(nspec,ntfull)    :: sspi
  REAL(SP), INTENT(inout), DIMENSION(nspec,ntfull) :: sspo
  REAL(SP), DIMENSION(nemline) :: tmpnebline
  REAL(SP), DIMENSION(nspec)   :: tmpnebcont,nu
  REAL(SP), DIMENSION(nspec,nemline) :: tmparr

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  !locate the maximum nebular age point in the full time array
  !nti = locate(time_full,nebem_age(nebnage))
  nti = locate(time_full,7.d0)

  !set limits on the velocity dispersion for broadening
  IF (smooth_velocity.EQ.1) THEN
     sigma = MAX(pset%sigma_smooth,neb_res_min)
  ELSE
     dlam  = MAX(pset%sigma_smooth,1.0)
  ENDIF

  !set up the interpolation variables for logZ and logU
  z1 = MAX(MIN(locate(nebem_logz,pset%gas_logz),nebnz-1),1)
  dz = (pset%gas_logz-nebem_logz(z1))/(nebem_logz(z1+1)-nebem_logz(z1))
  dz = MAX(MIN(dz,1.0),0.0) !no extrapolations
  u1 = MAX(MIN(locate(nebem_logu,pset%gas_logu),nebnip-1),1)
  du = (pset%gas_logu-nebem_logu(u1))/(nebem_logu(u1+1)-nebem_logu(u1))
  du = MAX(MIN(du,1.0),0.0) !no extrapolations

  !set up a "master" array of normalized Gaussians
  !this makes the code run *much* faster
  DO i=1,nemline
     IF (smooth_velocity.EQ.1) THEN
        dlam = nebem_line_pos(i) * sigma/clight*1E13
     ENDIF
     tmparr(:,i) = 1/SQRT(2*mypi)/dlam*&
          EXP(-(spec_lambda-nebem_line_pos(i))**2/2/dlam**2)  / &
          clight*nebem_line_pos(i)**2
  ENDDO

  sspo = sspi

  DO t=1,nti

     !remove ionizing photons from the stellar source
     sspo(1:whlylim,t) = sspi(1:whlylim,t) * pset%frac_obrun

     !the number of ionizing photons is computed here
     !some fraction of the stars are "runaways" which means
     !that they are not embedded in the HII region
     qq = tsum(spec_nu(:whlylim),sspi(:whlylim,t)/spec_nu(:whlylim))/&
          hplank*lsun 
     qq = qq * (1-pset%frac_obrun)

     !set up age interpolant
     a1 = MAX(MIN(locate(nebem_age,time_full(t)),nebnage-1),1)
     da = (time_full(t)-nebem_age(a1))/(nebem_age(a1+1)-nebem_age(a1))
     da = MAX(MIN(da,1.0),0.0) !no extrapolations
    
     !add nebular continuum emission
     IF (add_neb_continuum.EQ.1) THEN
        tmpnebcont = &   !interpolate in Zgas, logU, age
             (1-dz)*(1-da)*(1-du)* nebem_cont(:,z1,a1,u1)+&
             (1-dz)*(1-da)*(du)*   nebem_cont(:,z1,a1,u1+1)+&
             (1-dz)*(da)*(1-du)*   nebem_cont(:,z1,a1+1,u1)+&
             (1-dz)*(da)*(du)*     nebem_cont(:,z1,a1+1,u1+1)+&
             (dz)*(1-da)*(1-du)*   nebem_cont(:,z1+1,a1,u1)+&
             (dz)*(1-da)*(du)*     nebem_cont(:,z1+1,a1,u1+1)+&
             (dz)*(da)*(1-du)*     nebem_cont(:,z1+1,a1+1,u1)+&
             (dz)*(da)*(du)*       nebem_cont(:,z1+1,a1+1,u1+1)
        sspo(:,t) = sspo(:,t) + 10**tmpnebcont * qq
     ENDIF

     !add line emission
     tmpnebline = &    !interpolate in Zgas, logU, age
          (1-dz)*(1-da)*(1-du)* nebem_line(:,z1,a1,u1)+&
          (1-dz)*(1-da)*(du)*   nebem_line(:,z1,a1,u1+1)+&
          (1-dz)*(da)*(1-du)*   nebem_line(:,z1,a1+1,u1)+&
          (1-dz)*(da)*(du)*     nebem_line(:,z1,a1+1,u1+1)+&
          (dz)*(1-da)*(1-du)*   nebem_line(:,z1+1,a1,u1)+&
          (dz)*(1-da)*(du)*     nebem_line(:,z1+1,a1,u1+1)+&
          (dz)*(da)*(1-du)*     nebem_line(:,z1+1,a1+1,u1)+&
          (dz)*(da)*(du)*       nebem_line(:,z1+1,a1+1,u1+1)

     DO i=1,nemline
        sspo(:,t) = sspo(:,t) + 10**tmpnebline(i)*qq*tmparr(:,i)
     ENDDO

  ENDDO


END SUBROUTINE ADD_NEBULAR
