SUBROUTINE ADD_NEBULAR(pset,sspi,sspo)

  USE sps_vars; USE sps_utils, ONLY : locate
  IMPLICIT NONE

  INTEGER :: t,i,nti,a1,z1,u1
  REAL(SP) :: da,dz,du,sigma,dlam
  TYPE(PARAMS), INTENT(in) :: pset
  REAL(SP), INTENT(in), DIMENSION(nspec,ntfull)    :: sspi
  REAL(SP), INTENT(inout), DIMENSION(nspec,ntfull) :: sspo

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  !locate the maximum nebular age point in the full time array
  nti = locate(time_full,nebem_age(nebnage))

  !set the velocity dispersion for broadening
  IF (pset%sigma_smooth.GT.tiny_number) THEN 
     sigma = pset%sigma_smooth
  ELSE
     sigma=10.
  ENDIF

  DO t=1,nti

     !interpolate in Zgas, logU, age
     a1 = MAX(MIN(locate(nebem_age,time_full(t)),nebnage-1),1)
     z1 = MAX(MIN(locate(nebem_zgas,10**pset%gas_logzsol*zsol),&
          nebnz-1),1)
     u1 = MAX(MIN(locate(nebem_logu,pset%gas_logu),nebnip-1),1)

     !remove ionizing photons from the stellar source
     sspo(:,t) = sspi(:,t)
     
     !add nebular continuum emission
     !note that the continuum tables are stored per unit U.
   !  sspo(:,t) = sspo(:,t) + nebem_cont(:,1,1) !* 10**pset%gas_logu

     !add line emission
     DO i=1,nemline
        dlam = nebem_line_pos(i) * sigma/clight*1E13
        !sspo(:,t) = sspo(:,t) + nebem_line(i,1,1,1)/SQRT(2*mypi)/dlam*&
        !     EXP(-(spec_lambda-nebem_line_pos(i))**2/2/dlam**2)
     ENDDO

  ENDDO



END SUBROUTINE ADD_NEBULAR
