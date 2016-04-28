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
     intsfwght = (sfwght_log(sspind, logt(2), sfh) - sfwght_log(sspind, logt(1), sfh))
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

  real(SP) :: loge, expi
  real(SP) :: logage, tprime ! intermediate time variables
  real(SP) :: a, b, c ! dummy variables used to break up long expressions

  loge = log10(exp(1.0))

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
     sfwght_log = (logage - logt) * exp(tprime) + loge * expi(tprime)
     
  else if (sfh%type.eq.4) then
     ! SFR = delayed exponential ~ T/tau exp(-T/tau)
     tprime = 10**logt / sfh%tau ! t/tau
     a = (10**logt - sfh%tage - sfh%tau) * (logt - logage)
     b = sfh%tau * loge
     c = (sfh%tage + sfh%tau) * loge
     sfwght_log = (a - b) * exp(tprime) + c * expi(tprime)
     
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
  real(SP) :: loge, a

  loge = log10(exp(1.0))

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


FUNCTION expi(arg)
  ! Computes the exponential integral Ei(x) for x > 0.
  ! Parameters:
  !   `eps` is the relative error, or absolute error near the zero of Ei at x=0.3725;
  !   `eul` is Euler’s constant
  !   `maxit` is the maximum number of iterations allowed;
  !   `fmin` is a number near the smallest representable floating-point number.
  ! Adapted from the NR *public domain* code `ei` (http://numerical.recipes/pubdom/nr.f90.txt)
  use sps_vars, only: SP
  implicit none
  
  INTEGER, PARAMETER :: maxit=1000
  REAL(SP) :: expi, arg
  REAL(SP), PARAMETER :: eps=6.e-8, eul=.57721566, fmin=1.d-70
  INTEGER :: k
  REAL(SP) :: fact,prev,sum,term
  
  if (arg.le.0.) then
     write(*,*) "EXPI: arg < 0"
     STOP
  endif
  
  ! Special case: avoid failure of convergence test because of underflow.
  if (arg.lt.fmin) then   
     expi = log(arg) + eul
  !else if (arg.le.-log(eps)) then ! Use power series.
  else if (.true.) then ! always use power series....
     sum = 0.
     fact = 1.
     do k=1, maxit
        fact = fact * arg / k
        term = fact / k
        sum = sum + term
        if (term.lt.(eps*sum)) then
           continue
        else if (k.eq.maxit) then
           write(*,*) 'EXPI: Series failed to converge.'
           STOP
        endif
     enddo
     expi = sum + log(arg) + eul
  else ! Use asymptotic series.
     sum = 0. ! Start with second term.
     term = 1.
     do k=1, maxit
        prev = term
        term = term * k / arg
        ! Since final sum is greater than one, term itself approximates the relative error.
        if (term.lt.eps) continue
        if (term.lt.prev) then
           sum = sum + term ! Still converging: add new term.
        else
           sum = sum - prev ! Diverging: subtract previous term and exit.
           continue
        endif
     enddo
     expi = exp(arg) * (1. + sum) / arg
  endif
  return
END FUNCTION expi
