subroutine sfhinfo(pset, age, mfrac, sfr, frac_linear)
  ! Get the SFR integrated from T=0 to T=age, normalized by the SFR integrated
  ! from T=0 to T=Tmax, where Tmax is the maximum isochrone/SSP age.
  !
  ! Inputs
  ! ---------
  !
  ! pset:
  !   A `PARAMS` structure containing the SFH parameters.
  !
  ! age:
  !   The age (in forward Gyr) at which to calculate the SFR and mass
  !   fractions.
  !
  ! Outputs:
  ! ----------
  !
  ! mfrac:
  !    Fraction of mass formed by `age` \equiv mass_formed(age) /
  !    mass_formed(maxtime) where maxtime is the age of the oldest isochrone.
  !
  ! sfr:
  !    The SFR at `age`, normalized s.t. the integral of the SFR from 0 to the
  !    oldest isochrone age is 1.  Infinite SFR from bursts is not included.
  !
  ! frac_linear:
  !   If pset%sfh=5, this gives the fraction of m_formed(age) that was formed
  !   in the linear portion.
  !
  use sps_vars, only: time_full, ntfull, tiny_number, &
                      sfh_tab, ntabsfh, &
                      PARAMS, SP
  use sps_utils, only: locate
  implicit none

  type(PARAMS), intent(in) :: pset
  real(SP), intent(in) :: age

  real(SP), intent(out) :: mfrac, sfr, frac_linear
  
  real(SP) :: Tmax, Tprime, Tz, Ttrunc, Thi
  real(SP) :: m, gammainc
  real(SP) :: mass_tau, mass_linear, mfrac_burst
  real(SP) :: total_mass_tau, total_mass_linear
  real(SP) :: sfr_tau, sfr_linear, sfr_trunc, sfr_const
  integer :: power, itab

  ! Defaults
  mfrac = 1.
  sfr = 0.
  frac_linear = 0.
  
  if (pset%sfh.eq.0) then
     ! SSPs
     mfrac = 1.0
     sfr = 0. ! actually the sfr is infinity
     
  else if ((pset%sfh.eq.1).or.(pset%sfh.eq.4).or.(pset%sfh.eq.5)) then
     ! Compute tau model component, for SFH=1,4,5
     !
     ! Integration limits are from 0 to Tmax and 0 to Tprime, where
     !   - Tmax is the maximum isochrone age, and
     !   - Tprime is the given `age`,
     ! both adjusted for sf_start, and we'll clip them to Ttrunc.
     Tmax = 10**(time_full(ntfull) - 9) - pset%sf_start
     Tprime = age - pset%sf_start
     ! Ttrunc only matters if nonzero and not before sf_start.
     if ((pset%sf_trunc.lt.tiny_number).or.(pset%sf_trunc.lt.pset%sf_start)) then
        Ttrunc = Tmax
     else
        Ttrunc = pset%sf_trunc - pset%sf_start
     endif

     ! Now integrate to get mass formed by Tprime and by Tmax, dealing with
     ! truncation that happens after sf_start but before Tmax and/or Tprime.
     if (pset%sfh.eq.1) power = 1
     if ((pset%sfh.eq.4).or.(pset%sfh.eq.5)) power = 2
     total_mass_tau = pset%tau * gammainc(power, min(Tmax, Ttrunc) / pset%tau)
     if (Tprime.lt.0) then
        ! Deal with sf that hasn't started yet
        mass_tau = 0.
        sfr_tau = 0.
     else
        mass_tau = pset%tau * gammainc(power, min(Tprime, Ttrunc) / pset%tau)
        ! The SFR at Tprime (unnormalized)
        sfr_tau = (min(Tprime, Ttrunc) / pset%tau)**(power-1.) * exp(-min(Tprime, Ttrunc) / pset%tau)
     endif

  endif

  ! Add the constant and burst portions, for SFH=1,4.
  if ((pset%sfh.eq.1).or.(pset%sfh.eq.4)) then

     ! Fraction of the burst mass formed by Tprime
     if (Tprime.gt.(pset%tburst-pset%sf_start)) then
        mfrac_burst = 1.0
     else
        mfrac_burst = 0.0
     endif
     ! SFR from constant portion at Tprime
     if (Tprime.gt.0) then
        sfr_const = 1.0 / min(Tmax, Ttrunc)
     else
        sfr_const = 0.0
     endif

     ! Add formed mass fractions for each component, weighted by component fractions.
     ! Fraction of the constant mass formed by Tprime is just Tprime/Tmax
     mfrac = (1. - pset%const - pset%fburst) * mass_tau / total_mass_tau + &
          pset%const * max(min(Tprime, Ttrunc), 0.0) * sfr_const + &
          pset%fburst * mfrac_burst

     ! N.B. for Tprime = tburst, sfr is infinite, but we ignore that case.
     ! For constant SFR=1, total_mass_const = Tmax
     if (Tprime.gt.Ttrunc) then
        ! We've truncated.
        sfr = 0.
     else
        sfr = (1. - pset%const - pset%fburst) * sfr_tau / total_mass_tau + &
              pset%const * sfr_const
     endif

  ! Add the linear portion, for Simha, SFH=5.
  ! This is the integral of sfr_trunc*(1 - m * (T - Ttrunc)) from Ttrunc to Tz
  else if (pset%sfh.eq.5) then
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
        sfr_trunc = (Ttrunc/pset%tau) * exp(-Ttrunc/pset%tau)
        total_mass_linear = sfr_trunc * ((Thi - Ttrunc) - m/2.*(Thi**2 + Ttrunc**2) + &
                                         m * Thi * Ttrunc)

        if (Ttrunc.ge.Tprime) then
           ! But Tprime has not reached truncation yet.
           mass_linear = 0.0
           sfr = sfr_tau / (total_mass_tau + total_mass_linear)
        else
           ! Tprime is past truncation, how much mass formed?
           Thi = min(Tz, Tprime)
           mass_linear = sfr_trunc * ((Thi - Ttrunc) - m/2.*(Thi**2 + Ttrunc**2) + &
                                      m * Thi * Ttrunc)
           sfr = max(sfr_trunc * (1 - m * (Tprime - Ttrunc)), 0.) / (total_mass_tau + total_mass_linear)
        endif

     endif

     ! Fractions including tau and linear portions, wrapped in if/else to avoid NaNs
     if (Tprime.ge.0) then
        mfrac = (mass_tau + mass_linear) / (total_mass_tau + total_mass_linear)
        ! Mass fraction of the simha SFH at Tprime that formed in the linear portion.
        frac_linear = mass_linear / (mass_linear + mass_tau)
     else
        mfrac = 0.0
        frac_linear = 0.0
        sfr = 0.0
     endif

  endif

  ! Tabular.  Simple linear interpolation to get the sfr.
  ! The table is in units of yrs of forward time and M_sun/yr.
  if ((pset%sfh.eq.2).or.(pset%sfh.eq.3)) then
     itab = max(min(locate(sfh_tab(1, 1:ntabsfh), age*1e9), ntabsfh-1), 1)
     m = (sfh_tab(2, itab+1) - sfh_tab(2, itab)) / (sfh_tab(1, itab+1) - sfh_tab(1, itab))
     sfr = sfh_tab(2, itab) + m * (age*1e9 - sfh_tab(1, itab))
     sfr = max(sfr, 0.0) * 1e9 ! convert to per Gyr
  endif
  
  ! Convert SFR from per Gyr to per year
  sfr = sfr / 1e9
  
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
