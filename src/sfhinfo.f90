subroutine sfhinfo(pset, age, mfrac, sfr, frac_linear)
  ! Get the SFR integrated from T=0 to T=age, normalized by the SFR integrated
  ! from T=0 to T=Tmax, where Tmax is the maximum isochrone/SSP age.
  !
  ! Inputs
  ! ---------
  !
  ! pset:
  !   A `PARAMS` structure containing the SFH parameters
  !
  ! age:
  !   The age (in Gyr) at which to calculate the SFR and mass fractions
  !
  ! Outputs:
  ! ----------
  !
  ! mfrac:
  !    Fraction of mass formed by `age` \equiv mass_formed(age) /
  !    mass_formed(maxtime) where maxtime is the age of the oldest isochrone.
  !
  ! sfr:
  !    The SFR at `age`, normalized s.t. the integral of the SFR from 0 to
  !    maxtime is 1.
  !
  ! frac_linear:
  !   If pset%sfh=5, this gives the fraction of m_formed(age) that was formed
  !   in the linear portion.
  !
  use sps_vars, only: time_full, ntfull, PARAMS, SP
  implicit none
  type(PARAMS), intent(in) :: pset
  real(SP), intent(in) :: age

  real(SP), intent(out) :: mfrac, sfr, frac_linear
  
  real(SP) :: Tmax, Tprime, Tz, Ttrunc, Thi
  real(SP) :: m, gammainc
  real(SP) :: mass_tau, mass_linear, mfrac_burst
  real(SP) :: total_mass_tau, total_mass_linear
  real(SP) :: sfr_tau, sfr_linear, sfr_trunc
  integer :: power

  ! Defaults, used for SFH=2,3
  mfrac = 0.
  sfr = 0.
  frac_linear = 0.
  
  if (pset%sfh.eq.0) then
     mfrac = 1.0
     sfr = 0. ! actually the sfr is infinity

  else if ((pset%sfh.eq.1).or.(pset%sfh.eq.4).or.(pset%sfh.eq.5)) then
     ! Compute tau model component, for SFH=1,4,5
     !
     ! Integration limits are from 0 to Tmax and 0 to Tprime, where
     !   Tmax is the maximum isochrone age, and
     !   Tprime is the given `age`,
     !   both adjusted for sf_start and clipped to Ttrunc
     Tmax = 10**(time_full(ntfull) - 9) - pset%sf_start
     Tprime = age - pset%sf_start
     Ttrunc = pset%sf_trunc - pset%sf_start
     ! Deal with truncation that happens after sf_start but before Tmax and/or Tprime.
     if (Ttrunc.gt.0) then
        Tmax = min(Tmax, Ttrunc)
        Tprime = min(Tprime, Ttrunc)
     endif

     ! Now integrate to get mass formed by Tprime and by Tmax
     if (pset%sfh.eq.1) power = 1
     if ((pset%sfh.eq.4).or.(pset%sfh.eq.5)) power = 2
     total_mass_tau = pset%tau * gammainc(power, Tmax/pset%tau)
     mass_tau = pset%tau * gammainc(power, Tprime/pset%tau)
     ! The SFR at Tprime (unnormalized)
     sfr_tau = (Tprime/pset%tau)**(power-1.) * exp(Tprime/pset%tau)
  endif
  
  ! Add the constant and burst portions, for SFH=1,4.
  if ((pset%sfh.eq.1).or.(pset%sfh.eq.4)) then
     
     !Fraction of the burst mass formed by Tprime.
     if (Tprime.gt.(pset%tburst-pset%sf_start)) then
        mfrac_burst = 1.0
     else
        mfrac_burst = 0.0
     endif
     ! Add formed mass fractions for each component, weighted by component fractions.
     ! Fraction of the constant mass formed by Tprime is just Tprime/Tmax
     mfrac = (1. - pset%const - pset%fburst) * mass_tau / total_mass_tau + &
          pset%const * Tprime / Tmax + &
          pset%fburst * mfrac_burst
     ! N.B. for Tprime = tburst, sfr is infinite, but we ignore that case.
     ! For constant SFR=1, total_mass_const = Tmax
     if (Tprime.eq.Ttrunc) then
        ! We've truncated.
        sfr = 0
     else
        sfr = (1. - pset%const - pset%fburst) * sfr_tau / total_mass_tau + &
              pset%const * 1.0 / Tmax
     endif

  ! Add the linear portion, for Simha, SFH=5.
  ! This is the integral of sfr_trunc*(1 - m * (T - Ttrunc)) from Ttrunc to Tz
  else if (pset%sfh.eq.5) then
     ! Need to recalculate Tprime and Tmax, since we don't want them clipped by
     ! Ttrunc anymore
     Tprime = age - pset%sf_start
     Tmax = 10**(time_full(ntfull) - 9) - pset%sf_start
     m = -pset%sf_slope
     if (m.gt.0) then
        ! find time at which SFR=0, if m>0
        Tz = Ttrunc + 1.0 / m
     else
        ! m <= 0 will never reach SFR=0
        Tz = Tmax
     endif
     
     ! Logic for Linear portion
     if ((Ttrunc.le.0).or.(Ttrunc.gt.Tmax)) then
        ! Truncation does not occur during the SFH.
        total_mass_linear = 0.
        mass_linear = 0.
        sfr = sfr_tau / total_mass_tau
        frac_linear = 0.
     else
        ! Truncation does occur, integrate linear to zero crossing or Tmax.
        Thi = min(Tz, Tmax)
        sfr_trunc = (Ttrunc/pset%tau) * exp(Ttrunc/pset%tau) / total_mass_tau
        total_mass_linear = sfr_trunc * ((Thi - Ttrunc) - m/2.*(Thi**2 + Ttrunc**2) + &
                                         m * Thi * Ttrunc)

        if (Ttrunc.ge.Tprime) then
           ! But Tprime has not reached truncation yet.
           mass_linear = 0.0
           sfr = sfr_tau / total_mass_tau
        else
           ! Tprime is past truncation how much mass formed?
           Thi = min(Tz, Tprime)
           mass_linear = sfr_trunc * ((Thi - Ttrunc) - m/2.*(Thi**2 + Ttrunc**2) + &
                                      m * Thi * Ttrunc)
           sfr = max(sfr_trunc * (1 - m * (Tprime - Ttrunc)), 0.)
        endif
        
     endif

     ! Fractions including tau and linear portions
     mfrac = (mass_tau + mass_linear) / (total_mass_tau + total_mass_linear)
     ! Mass fraction of the simha SFH at Tprime that formed in the linear portion.
     frac_linear = mass_linear / (mass_linear + mass_tau)
     
  endif

end subroutine sfhinfo

function gammainc(power, arg)
  !
  ! Calculate incomplete gamma for a = 1 or 2

  use sps_vars, only: SP
  implicit none
  integer, intent(in) :: power
  real(SP), intent(in) :: arg

  real(SP) :: gammainc

  if (power.eq.2) then
     gammainc = 1.0 - exp(-arg) - arg * exp(-arg)
  else if (power.eq.1) then
     gammainc = 1.0 - exp(-arg)
  else
     write(*,*) "gammainc: power must be 1 or 2"
     STOP
  endif

end function gammainc
