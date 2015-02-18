FUNCTION INTSFR(sfh,tau,const,sftrunc,sfstart,sftheta,t1,t2,tweight)

  !simple routine to integrate the SFH from t1 to t2

  USE sps_vars; USE sps_utils, ONLY : tsum, locate
  IMPLICIT NONE

  INTEGER, INTENT(in)  :: sfh
  REAL(SP), INTENT(in) :: t1,t2,tau,const,sftrunc,sfstart,sftheta
  INTEGER, INTENT(in), OPTIONAL :: tweight
  INTEGER  :: lo,hi
  REAL(SP) :: s1,s2,dt,intsfr,tt1

  !-----------------------------------------------------!
  !-----------------------------------------------------!

  IF (t1.GT.sftrunc.AND.t2.LT.sftrunc) THEN
     tt1 = sftrunc
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

     IF (t2.GT.sftrunc) THEN
        intsfr = 0.0
     ENDIF

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

     IF (t2.GT.sftrunc) THEN
        intsfr = 0.0
     ENDIF

  ELSE IF (sfh.EQ.5) THEN
     
     IF ((t1+sfstart).LT.sftrunc) THEN
        intsfr = (EXP(-t2/tau/1E9)*(1+t2/tau/1E9) - &
             EXP(-t1/tau/1E9)*(1+t1/tau/1E9)) * (tau*1E9)**2 
     ELSE
       intsfr =  TAN(sftheta)*&
             (0.5*(t1**2-t2**2)-(t1-t2)*(sftrunc-sfstart))
     ENDIF
     intsfr = MAX(intsfr,0.0) / 1E10

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
