! ------------------------------------
! Routines used to evaluate the indefinite integrals that arise when exactly
! integrating interpolated spectra weighted by particular SFHs.
!
! Adding new SFHs is as easy as adding new sfh%type cases in both sfhint_*
! functions below (and making sure the sfhlimits.f90 make sense)
! ------------------------------------

function intsfwght(sspind, logt, sfh)
  ! Wrapper on the sfhint_* routines to choose the correct interpolation type
  ! and calculate the definite integral of the weighted SFR between the given
  ! limits.
  !
  ! Inputs
  ! -------
  !
  ! sspind:
  !    The index of the SSP forming the bracket with the SSP you are getting a
  !    weight for.
  !
  ! logt:
  !    the limits of the integral, 2-element array (lo, hi) of
  !    log(lookback time) (in yrs).
  !
  ! sfh:
  !    Structure containing the SFH parameters in units of years, including an
  !    integer `type` specifiying the form of SFR(t).
  !
  ! Outputs
  !--------
  !  intsfwght:
  !    The exact definite integral between the limits specified in `logt`.
  
  use sps_vars, only: interpolation_type, SFHPARAMS, SP
  implicit none
  
  integer, intent(in) :: sspind
  real(SP), intent(in), dimension(2) :: logt
  type(SFHPARAMS), intent(in) :: sfh

  real(SP) :: intsfwght, sfwght_log, sfwght_lin

  if (interpolation_type.eq.0) then
     intsfwght =  sfwght_log(sspind, logt(2), sfh) - sfwght_log(sspind, logt(1), sfh)
  else if (interpolation_type.eq.1) then
     intsfwght = (sfwght_lin(sspind, 10**logt(2), sfh) - sfwght_lin(sspind, 10**logt(1), sfh))
  endif
  
end function intsfwght
   

function sfwght_log(sspind, logt, sfh)
  ! Evaluates the indefinite integral of the interpolation weight in log time,
  ! weighted by the SFH, and evaluated at `logt`.  In detail, this function
  ! returns:
  !   (\int dt \, \mathrm{SFR}(t) \, (x - \log t) )|_{\texttt{logt}}
  ! where x = time_full(ssp_ind) and `t` is *lookback* time.
  !
  ! Inputs
  ! -------
  !
  ! sspind:
  !    The index of the SSP forming the bracket with the SSP you are getting a
  !    weight for.
  !
  ! logt:
  !    Where the indefinite integral is evaluated, log(t_lookback) (in yrs).
  !
  ! sfh:
  !    Structure containing the SFH parameters in units of years, including an
  !    integer `type` specifiying the form of SFR(t).
  !
  ! Outputs
  !--------
  !  sfwght_log:
  !    The exact indefinite integral, evaluated at `logt`.  Scalar float.
  
  use sps_vars, only: time_full, tiny_logt, SFHPARAMS, SP
  implicit none
  
  integer, intent(in) :: sspind
  real(SP), intent(in) :: logt
  type(SFHPARAMS), intent(in) :: sfh

  real(SP) :: sfwght_log

  real(SP) :: loge, ei
  real(SP) :: logage, tprime ! intermediate time variables
  real(SP) :: a, b, c ! dummy variables used to break up long expressions

  REAL(SP), PARAMETER :: e=2.7182818284590452353602874713526624977572_sp
  loge = log10(e)

  ! Zero index means use ~0 age.
  if (sspind.gt.0) then
     logage = time_full(sspind)
  else
     logage = tiny_logt
  endif
     
  if (sfh%type.eq.0) then
     ! SFR = Constant ~ 1
     sfwght_log = 10**logt * (logage - logt + loge)
     
  else if (sfh%type.eq.1) then
     ! SFR = exponential ~ exp(-T/tau)
     tprime = 10**logt / sfh%tau
     sfwght_log = (logage - logt) * exp(tprime) + loge * ei(tprime)

  else if (sfh%type.eq.4) then
     ! SFR = delayed exponential ~ T/tau exp(-T/tau)
     tprime = 10**logt / sfh%tau ! t/tau
     a = (10**logt - sfh%tage - sfh%tau) * (logt - logage)
     b = sfh%tau * loge
     c = (sfh%tage + sfh%tau) * loge
     sfwght_log = (a - b) * exp(tprime) + c * ei(tprime)
     
  else if (sfh%type.eq.5) then
     !SFR = linear ~ (1 - sf_slope * (T - T_trunc)), T > T_trunc
     tprime = max(0.0, sfh%tage - sfh%sf_trunc) !t_q
     a = 1 - sfh%sf_slope * tprime
     b = a * 10**logt * (logage - logt + loge)
     c = sfh%sf_slope * (10**logt)**2 / 2 * (logage - logt + loge / 2)
     sfwght_log = b + c
     
  endif
  
end function sfwght_log


function sfwght_lin(sspind, t, sfh)
  ! Evaluates the indefinite integral of the interpolation weight in linear
  ! time, weighted by the SFH, and evaluated at `logt`.  In detail, this
  ! function returns:
  !   (\int dt \, \mathrm{SFR}(t) \, (x - t) )|_{\texttt{t}}
  ! where x = 10**time_full(ssp_ind), and `t` is *lookback* time.
  !
  ! Inputs
  ! -------
  !
  ! sspind:
  !    The index of the SSP forming the bracket with the SSP you are getting a
  !    weight for.
  !
  ! t:
  !    Where the indefinite integral is evaluated, linear years (lookback time).
  !
  ! sfh:
  !    Structure containing the SFH parameters in units of years, including an
  !    integer `type` specifiying the form of SFR(t).
  !
  ! Outputs
  !--------
  !  sfwght_lin:
  !    The indefinite integral, evaluated at `t`

  use sps_vars, only: time_full, tiny_logt, SFHPARAMS, SP
  implicit none
  
  integer, intent(in) :: sspind 
  real(SP), intent(in) :: t
  type(SFHPARAMS), intent(in) :: sfh

  real(SP) :: sfwght_lin
  
  real(SP) :: age, tprime
  real(SP) :: a

  ! Convert from log(age_ssp) to age_ssp,
  ! accounting for the case sspind=0 (where age~0)
  if (sspind.gt.0) then
     age = 10**time_full(sspind)
  else
     age = 10**tiny_logt
  endif

  if (sfh%type.eq.0) then
     ! SFR = Constant ~ 1
     sfwght_lin = age * t - t**2 / 2

  else if (sfh%type.eq.1) then
     ! SFR = exponential ~ 1/tau * exp(-T/tau)
     tprime = t / sfh%tau
     sfwght_lin = (age - t + sfh%tau) * exp(tprime)

  else if (sfh%type.eq.4) then
     ! SFR = delayed exponential ~ T/tau**2 * exp(-T/tau)
     tprime = t / sfh%tau
     a = sfh%tage * age - (sfh%tage + age) * (t - sfh%tau) + &
          t**2 - 2*t*sfh%tau + 2*sfh%tau**2
     sfwght_lin = a * exp(tprime)
     
  else if (sfh%type.eq.5) then
     ! SFR = linear ~ (1 - sf_slope * (T - T_trunc)), T > T_trunc
     tprime = max(0.0, sfh%tage - sfh%sf_trunc) !t_q
     a = 1 - sfh%sf_slope * tprime
     sfwght_lin = a * age * t + (sfh%sf_slope*age - a) * t**2 / 2 - sfh%sf_slope * t**3 / 3

  endif

end function sfwght_lin


FUNCTION EI(X)
!
! ============================================
! Purpose: Compute exponential integral Ei(x)
! Input :  x  --- Argument of Ei(x)
! Output:  EI --- Ei(x)
!
! Shanjie Zhang and Jianming Jin
!
!  Copyrighted but permission granted to use code in programs.
!  Buy their book "Computation of Special Functions", 1996, John Wiley & Sons, Inc.
!
!       ============================================
!
  implicit none
  INTEGER, PARAMETER :: SP = KIND(1.d0)
  REAL(SP) :: X, EI

  INTEGER :: K, MAXIT=1000
  REAL(SP), PARAMETER :: GA=0.5772156649015328D0, eps=1.0D-20
  REAL(SP) :: R
  
  IF (X.EQ.0.0) THEN
     EI=-1.0D+300
  ELSE IF (X .LT. 0) THEN
     WRITE(*,*) "EI: X < 0."
     STOP
  ELSE IF (DABS(X).LE.40.0) THEN
     ! Power series around x=0
     EI=1.0D0
     R=1.0D0
     DO K = 1, MAXIT
        R = R*K*X / (K+1.0D0)**2
        EI = EI+R
        IF (DABS(R/EI).LE.eps) THEN
           CONTINUE
        ELSE IF (K.eq.MAXIT) then
           write(*,*) "EI: Series failed to converge."
           STOP
        ENDIF
           
     ENDDO
     EI=GA+DLOG(X)+X*EI

  ELSE
     ! Asymptotic expansion (the series is not convergent)
     EI=1.0D0
     R=1.0D0
     DO K=1,20
        R=R*K/X
        EI=EI+R
     ENDDO
     EI=DEXP(X)/X*EI
  ENDIF
  RETURN
END FUNCTION EI
