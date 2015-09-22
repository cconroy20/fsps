FUNCTION INTSFR(sfh,tau,const,sftrunc,sfstart,sfslope,tmax,t1,t2,tweight)

  !simple routine to integrate the SFH from t1 to t2
  !the edges should now be handled correctly (as of 7/15)

  USE sps_vars; USE sps_utils, ONLY : tsum, locate
  IMPLICIT NONE

  INTEGER, INTENT(in)  :: sfh
  REAL(SP), INTENT(in) :: t1,t2,tau,const,sftrunc,sfstart,sfslope,tmax
  INTEGER, INTENT(in), OPTIONAL :: tweight
  INTEGER  :: lo,hi
  REAL(SP) :: s1,s2,dt,intsfr,tt1,norm,sft

  !-----------------------------------------------------!
  !-----------------------------------------------------!

  intsfr=0.0

  !changed to sftrunc-sfstart (7/15)
  IF (t1.GT.(sftrunc-sfstart)) THEN
     tt1 = sftrunc-sfstart
  ELSE
     tt1 = t1
  ENDIF

  IF (sfh.EQ.1) THEN

     IF (PRESENT(tweight)) THEN
        !in this case the integral is int(t*sfr)
        intsfr = tau*( (t2/tau/1E9+1)*EXP(-t2/tau/1E9) - &
             (tt1/tau/1E9+1)*EXP(-tt1/tau/1E9) ) / &
             (1-EXP(-(sftrunc-sfstart)/1E9/tau)) 
        intsfr = intsfr*(1-const) + const*(tt1**2-t2**2)/2E9/(sftrunc-sfstart)
     ELSE
        intsfr = (EXP(-t2/tau/1E9) - EXP(-tt1/tau/1E9))/&
             (1-EXP(-(sftrunc-sfstart)/1E9/tau)) 
        intsfr = intsfr*(1-const) + const*(tt1-t2)/(sftrunc-sfstart)
     ENDIF

     IF (t2.GT.(sftrunc-sfstart)) intsfr = 0.0

  ELSE IF (sfh.EQ.4) THEN

     IF (PRESENT(tweight)) THEN
        !in this case the integral is int(t*sfr)
        intsfr = tau*(EXP(-t2/tau/1E9)*(2+t2/tau/1E9*(t2/tau/1E9+2)) - &
                  EXP(-tt1/tau/1E9)*(2+tt1/tau/1E9*(tt1/tau/1E9+2))) / &
                  (1-EXP(-(sftrunc-sfstart)/1E9/tau)*((sftrunc-sfstart)/1E9/tau+1)) 
        intsfr = intsfr*(1-const) + const*(tt1**2-t2**2)/2E9/(sftrunc-sfstart)
     ELSE
        intsfr = (EXP(-t2/tau/1E9)*(1+t2/tau/1E9) - &
             EXP(-tt1/tau/1E9)*(1+tt1/tau/1E9)) / &
             (1-EXP(-(sftrunc-sfstart)/1E9/tau)*((sftrunc-sfstart)/1E9/tau+1))
        intsfr = intsfr*(1-const) + const*(tt1-t2)/(sftrunc-sfstart)
     ENDIF

     IF (t2.GT.(sftrunc-sfstart)) intsfr = 0.0

  ELSE IF (sfh.EQ.5) THEN

     !SFR at the transition time
     sft = ((sftrunc-sfstart)/tau/1E9)*&
          EXP(-(sftrunc-sfstart)/tau/1E9 )/tau/1E9
     norm = (1-EXP(-(sftrunc-sfstart)/1E9/tau)*((sftrunc-sfstart)/1E9/tau+1))
     !add the normalization due to the linearly declining comp.
     norm = norm + sft*(tmax-sftrunc)*(1-sfslope*sftrunc/1e9)+&
          sft/1E9*sfslope*0.5*(tmax**2-sftrunc**2)
  
     IF ((t1+sfstart).LE.sftrunc) THEN
        intsfr = (EXP(-t2/tau/1E9)*(1+t2/tau/1E9) - &
             EXP(-t1/tau/1E9)*(1+t1/tau/1E9))
     ELSE IF ((t1+sfstart).GT.sftrunc.AND.(t1+sfstart).LE.tmax) THEN
        IF ((t2+sfstart).GT.sftrunc) THEN
           !both t1 and t2 are >sftrunc
           intsfr = sft*(t1-t2)*(1-sfslope*sftrunc/1e9)+&
                sft/1E9*sfslope*0.5*((t1+sfstart)**2-(t2+sfstart)**2)
        ELSE
           !in this case t1>sftrunc but t2<sftrunc, so break integral in two
           intsfr = (EXP(-t2/tau/1E9)*(1+t2/tau/1E9) - &
                EXP(-(sftrunc-sfstart)/tau/1E9)*(1+(sftrunc-sfstart)/tau/1E9))
           intsfr = intsfr + sft*(t1+sfstart-sftrunc)*(1-sfslope*sftrunc/1e9)+&
                sft/1E9*sfslope*0.5*((t1+sfstart)**2-(sftrunc)**2)
        ENDIF
     ELSE
        intsfr = 0.0
     ENDIF

     write(*,*) t1/1E9,t2/1E9,intsfr,sftrunc/1E9,tmax/1E9

     intsfr = MAX(intsfr/norm,0.0)

  ELSE IF (sfh.EQ.2.OR.sfh.EQ.3) THEN

     !handle the edges separately
     hi = MIN(MAX(locate(sfh_tab(1,1:ntabsfh),t1),1),ntabsfh-1)
     dt = (t1-sfh_tab(1,hi))/(sfh_tab(1,hi+1)-sfh_tab(1,hi))
     s1 = MAX((1-dt)*sfh_tab(2,hi) + dt*sfh_tab(2,hi+1),0.0)

     lo = MIN(MAX(locate(sfh_tab(1,1:ntabsfh),t2),1),ntabsfh-1)
     dt = (t2-sfh_tab(1,lo))/(sfh_tab(1,lo+1)-sfh_tab(1,lo))
     s2 = MAX((1-dt)*sfh_tab(2,lo) + dt*sfh_tab(2,lo+1),0.0)

     IF ((hi-lo).LT.2) THEN
        intsfr = (t1-t2)*0.5*(s1+s2)
     ELSE
        lo = lo+1
        intsfr = TSUM(sfh_tab(1,lo:hi),sfh_tab(2,lo:hi))
        intsfr = intsfr + (t1-sfh_tab(1,hi))*0.5*(s1+sfh_tab(2,hi))
        intsfr = intsfr + (sfh_tab(1,lo)-t2)*0.5*(s2+sfh_tab(2,lo))
     ENDIF

  ENDIF
  
  
END FUNCTION INTSFR
