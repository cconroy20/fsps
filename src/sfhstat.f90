SUBROUTINE SFHSTAT(pos,model,ssfr6,ssfr7,ssfr8,ave_age)

  !compute basic statistics given a parameterized star formation history.
  !required inputs are the parameter set and a single element output
  !structure from compsp
  
  USE sps_vars
  IMPLICIT NONE

  TYPE(PARAMS), INTENT(in)    :: pos
  TYPE(COMPSPOUT), INTENT(in) :: model
  REAL(SP), INTENT(inout)     :: ssfr6,ssfr7,ssfr8,ave_age
  REAL(SP) :: dt
  
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  dt = 10**model%age/1E9 - pos%sf_start

  IF (dt.LT.0.) THEN
     WRITE(*,*) 'SFHSTAT ERROR: dt<0.0, stopping....'
     STOP
  ENDIF

  !compute mass-weighted stellar age of tau component
  IF (pos%sfh.EQ.1) THEN
     ave_age = pos%tau*(1.-EXP(-dt/pos%tau)*(dt/pos%tau+1.))/&
          (1.-EXP(-dt/pos%tau))
     ssfr6 = (EXP(-(dt-1E-3)/pos%tau)-EXP(-dt/pos%tau))/(1.-EXP(-dt/pos%tau))
     ssfr7 = (EXP(-(dt-1E-2)/pos%tau)-EXP(-dt/pos%tau))/(1.-EXP(-dt/pos%tau))
     ssfr8 = (EXP(-(dt-1E-1)/pos%tau)-EXP(-dt/pos%tau))/(1.-EXP(-dt/pos%tau))
  ELSE IF (pos%sfh.EQ.4) THEN
     ave_age = (2.-EXP(-dt/pos%tau)*(dt/pos%tau*(dt/pos%tau+2.)+2.))*&
          pos%tau/(1.-EXP(-dt/pos%tau)*(dt/pos%tau+1))
     ssfr6 = (EXP(-(dt-1E-3)/pos%tau)*((dt-1E-3)/pos%tau)-&
          EXP(-dt/pos%tau)*(dt/pos%tau))/&
          (1.-EXP(-dt/pos%tau)*(dt/pos%tau+1))
     ssfr7 = (EXP(-(dt-1E-2)/pos%tau)*((dt-1E-2)/pos%tau)-&
          EXP(-dt/pos%tau)*(dt/pos%tau))/&
          (1.-EXP(-dt/pos%tau)*(dt/pos%tau+1))
     ssfr8 = (EXP(-(dt-1E-1)/pos%tau)*((dt-1E-1)/pos%tau)-&
          EXP(-dt/pos%tau)*(dt/pos%tau))/&
          (1.-EXP(-dt/pos%tau)*(dt/pos%tau+1))
  ELSE
     WRITE(*,*) 'SFHSTAT ERROR: you should not be calling sfhstat '//&
          'for sfh types NE 1 or 4, stopping....'
     STOP
  ENDIF

  !add constant component
  ave_age = ave_age*(1.-pos%const) + pos%const*dt/2.
  ssfr6   = ssfr6  *(1.-pos%const) + pos%const*1E-3/dt
  ssfr7   = ssfr7  *(1.-pos%const) + pos%const*1E-2/dt
  ssfr8   = ssfr8  *(1.-pos%const) + pos%const*1E-1/dt

  !convert to lookback time
  ave_age = dt - ave_age

  !only add the burst if the burst has happened
  IF (10**model%age/1E9.GT.pos%tburst) &
       ave_age = (1.-pos%fburst)*ave_age + pos%fburst*pos%tburst
  IF (dt-pos%tburst.LE.1E-3) ssfr6 = ssfr6 + pos%fburst
  IF (dt-pos%tburst.LE.1E-2) ssfr7 = ssfr7 + pos%fburst
  IF (dt-pos%tburst.LE.1E-1) ssfr8 = ssfr8 + pos%fburst


  !convert from integral(SFR) to log(<SSFR>), in 1/Gyr
  ssfr6 = LOG10(MAX(ssfr6/model%mass_csp/1E-3,tiny_number))
  ssfr7 = LOG10(MAX(ssfr7/model%mass_csp/1E-2,tiny_number))
  ssfr8 = LOG10(MAX(ssfr8/model%mass_csp/1E-1,tiny_number)) 


END SUBROUTINE SFHSTAT
