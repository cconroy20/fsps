SUBROUTINE SFHSTAT(pos,model,ssfr6,ssfr7,ssfr8,ave_age)

  !compute basic statistics given a parameterized star formation history.
  !required inputs are the parameter set and a single element output
  !structure from compsp
  
  USE sps_vars; USE nrtype
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
  ELSE IF (pos%sfh.EQ.4) THEN
     ave_age = (2.-EXP(-dt/pos%tau)*(dt/pos%tau*(dt/pos%tau+2.)+2.))*&
          pos%tau/(1.-EXP(-dt/pos%tau)*(dt/pos%tau+1))
  ELSE
     WRITE(*,*) 'SFHSTAT ERROR: you should not be calling sfhstat '//&
          'for sfh types NE 1 or 4, stopping....'
     STOP
  ENDIF

  !add constant component
  ave_age = ave_age*(1.-pos%const) + pos%const*dt/2.

  !convert to lookback time
  ave_age = dt - ave_age

  !only add the burst if the burst has happened
  IF (10**model%age/1E9.GT.pos%tburst) &
       ave_age = (1.-pos%fburst)*ave_age + pos%fburst*pos%tburst


  ssfr6 = pos%const/pos%tage + (1.-pos%const-pos%fburst)*&
       EXP(-(pos%tage-pos%sf_start)/pos%tau)

  ssfr7 = pos%const/pos%tage + (1.-pos%const-pos%fburst)*&
       EXP(-(pos%tage-pos%sf_start)/pos%tau)

  ssfr8 = pos%const/pos%tage + (1.-pos%const-pos%fburst)*&
       EXP(-(pos%tage-pos%sf_start)/pos%tau)

  !convert from integral(SFR) to log(<SSFR>), in 1/Gyr
  ssfr6 = LOG10(ssfr6/model%mass_csp / 1E-3)
  ssfr7 = LOG10(ssfr7/model%mass_csp / 1E-2)
  ssfr8 = LOG10(ssfr8/model%mass_csp / 1E-1) 


  WRITE(*,*) dt,dt/pos%tau,ave_age,ssfr6,ssfr7,ssfr8


END SUBROUTINE SFHSTAT
