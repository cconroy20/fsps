SUBROUTINE ADD_NEBULAR(pset,sspi,sspo,nebemline)

  !routine to add nebular emission (both line and continuum)
  !to input SSPs (sspi).  Returns SSPs as output (sspo).

  USE sps_vars; USE sps_utils, ONLY : locate,tsum
  IMPLICIT NONE

  INTEGER :: t,i,nti,a1,z1,u1
  REAL(SP) :: da,dz,du,sigma,dlam,qq,compute_q
  TYPE(PARAMS), INTENT(in) :: pset
  REAL(SP), INTENT(in), DIMENSION(nspec,ntfull)    :: sspi
  REAL(SP), INTENT(inout), DIMENSION(nspec,ntfull) :: sspo
  REAL(SP), INTENT(inout), DIMENSION(nemline,ntfull), OPTIONAL :: nebemline
  REAL(SP), DIMENSION(nemline) :: tmpnebline
  REAL(SP), DIMENSION(nspec)   :: tmpnebcont,nu

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  !locate the maximum nebular age point in the full time array
  nti = locate(time_full,nebem_age(nebnage))
  !right now we only include nebular emission for ages<2x10^7 yr
  !nti = locate(time_full,7.d30)

  !set up the interpolation variables for logZ and logU
  z1 = MAX(MIN(locate(nebem_logz,pset%gas_logz),nebnz-1),1)
  dz = (pset%gas_logz-nebem_logz(z1))/(nebem_logz(z1+1)-nebem_logz(z1))
  dz = MAX(MIN(dz,1.0),0.0) !no extrapolation
  u1 = MAX(MIN(locate(nebem_logu,pset%gas_logu),nebnip-1),1)
  du = (pset%gas_logu-nebem_logu(u1))/(nebem_logu(u1+1)-nebem_logu(u1))
  du = MAX(MIN(du,1.0),0.0) !no extrapolation

  !set up a "master" array of normalized Gaussians
  !in sps_setup.f90 this makes the code much faster
  IF (setup_nebular_gaussians.EQ.0.AND.nebemlineinspec.EQ.1) THEN
     DO i=1,nemline
        IF (smooth_velocity.EQ.1) THEN
           !smoothing variable is km/s
           dlam = nebem_line_pos(i)*pset%sigma_smooth/clight*1E13
        ELSE
           !smoothing variable is A
           dlam = pset%sigma_smooth
        ENDIF
        !broaden the line to at least the resolution element 
        !of the spectrum (x2).
        dlam = MAX(dlam,neb_res_min(i)*2)
        gaussnebarr(:,i) = 1/SQRT(2*mypi)/dlam*&
             EXP(-(spec_lambda-nebem_line_pos(i))**2/2/dlam**2)  / &
             clight*nebem_line_pos(i)**2
     ENDDO
  ENDIF

  sspo = sspi
  nebemline = 0.0
  
  DO t=1,nti

     !remove ionizing photons from the stellar source
     sspo(1:whlylim,t) = sspi(1:whlylim,t)*MAX(MIN(pset%frac_obrun,1.0),0.0)

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

     IF (PRESENT(nebemline)) nebemline(:,t) = 10**tmpnebline*qq

     IF (nebemlineinspec.EQ.1) THEN
        DO i=1,nemline
           sspo(:,t) = sspo(:,t) + 10**tmpnebline(i)*qq*gaussnebarr(:,i)
        ENDDO
     ENDIF
        
  ENDDO


END SUBROUTINE ADD_NEBULAR
