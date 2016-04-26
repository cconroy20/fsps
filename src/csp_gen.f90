subroutine csp_gen(mass_ssp, lbol_ssp, spec_ssp, mdust_ssp, &
                   pset, tage, &
                   mass_csp, lbol_csp, spec_csp, mdust_csp)
  !
  ! Return the spectrum (and mass and lbol) of a composite stellar population.
  !
  ! Inputs
  ! ---------
  !
  ! mass_ssp, lbol_ssp, spec_ssp:
  !   The (surviving) stellar masses, bolometric luminosities, and spectra of
  !   the SSPs.
  !
  ! pset:
  !   A `PARAMS` structure containing the SFH parameters
  !
  ! tage:
  !   The age (in Gyr) at which the spectrum is desired.  Note that this can be
  !   different than pset%tage if the latter is 0.
  !
  ! Outputs
  ! ---------
  !
  ! mass_csp, lbol_csp, spec_csp:
  !   The (surviving) stellar masses, bolometric luminosity, and spectrum of
  !   the composite stellar population at tage, normalized to 1 M_sun *formed*.
  
  use sps_vars, only: ntfull, nspec, time_full, tiny_number, &
                      sfhparams, params, SP
  use sps_utils, only: locate, sfh_weight, sfhinfo
  implicit none

  real(SP), intent(in), dimension(ntfull) :: mass_ssp, lbol_ssp, mdust_ssp
  real(SP), intent(in), dimension(nspec, ntfull) :: spec_ssp
  type(PARAMS), intent(in) :: pset
  real(SP), intent(in) :: tage
  
  real(SP), intent(out) :: mass_csp, lbol_csp, mdust_csp
  real(SP), intent(out), dimension(nspec) :: spec_csp


  real(SP), dimension(ntfull) :: total_weights=0., w1=0., w2=0.
  integer :: i, j, imin=0, imax=ntfull, ntabsfh
  type(SFHPARAMS) :: sfhpars
  real(SP) :: mass, m1, m2, frac_linear, mfrac, sfr
  
  ! Build a structure containing useful units, numbers, and switches for the
  ! weight calculations.
  call convert_sfhparams(pset, tage, sfhpars)
  ! Only calculate SFH weights for SSPs up to tage (plus the next one).
  imin = 0
  imax = min(max(locate(time_full, log10(sfhpars%tage)) + 2, 1), ntfull)
  
  ! ----- Get SFH weights -----


  ! SSP.
  if (pset%sfh.eq.0) then
     ! Make sure to use SSP weighting scheme
     sfhpars%type = -1
     ! Use tage as the burst lookback time, instead of tage-tburst.
     sfhpars%tb = sfhpars%tage
     ! Only need weights at two SSP points.
     imin = min(max(locate(time_full, log10(sfhpars%tage)), 0), ntfull)
     imax = min(imin+1, ntfull)
     ! These come pre-normalized to 1 Msun
     total_weights = sfh_weight(sfhpars, imin, imax)
  endif


  ! Tau and delayed-tau.
  if ((pset%sfh.eq.1).or.(pset%sfh.eq.4)) then
     imin = 0
     sfhpars%type = pset%sfh
     total_weights = sfh_weight(sfhpars, imin, imax)
     ! Could save some loops by having proper normalization analytically from sfh_weight
     m1 = sum(total_weights(1:imax))
     if (m1.lt.tiny_number) m1 = 1.0
     total_weights = total_weights / m1
  endif
  
  ! Add constant and burst weights for SFH=1,4
  if (((pset%sfh.eq.1).or.(pset%sfh.eq.4)).and.&
       ((pset%const.gt.0).or.(pset%fburst.gt.tiny_number))) then
     imin = 0
     ! Constant
     sfhpars%type = 0
     w1 = sfh_weight(sfhpars, imin, imax)
     m1 = sum(w1(1:imax))
     ! Burst.  These weights come pre-normalized to 1 Msun.
     sfhpars%type = -1
     w2 = sfh_weight(sfhpars, imin, imax)
     ! Sum with proper relative normalization.  Beware divide by zero.
     if (m1.lt.tiny_number) m1 = 1.0
     total_weights = (1 - pset%const - pset%fburst) * total_weights + &
                      pset%const * (w1 / m1) + &
                      pset%fburst * w2
  endif

  
  ! Simha
  if (pset%sfh.eq.5) then
     imin = 0
     ! Delayed-tau model portion.
     sfhpars%type = 4
     w1 = sfh_weight(sfhpars, imin, imax)
     m1 = sum(w1(1:imax))
     ! Linear portion.  Need to set `use_simha_limits` flag to get correct limits.
     sfhpars%type = 5
     sfhpars%use_simha_limits = 1
     w2 = sfh_weight(sfhpars, imax, imax)
     sfhpars%use_simha_limits = 0
     m2 = sum(w2(1:imax))
     ! Normalize and sum.  need to be careful of divide by zero here, if all
     ! linear or all delay-tau s.t. weights sum to zero for a component.
     if (m1.lt.tiny_number) m1 = 1.0
     if (m2.lt.tiny_number) m2 = 1.0
     ! need to get the fraction of the mass formed in the linear portion
     call sfhinfo(pset, tage, mfrac, sfr, frac_linear)
     total_weights = (w1 / m1) * (1 - frac_linear) + frac_linear * (w2 / m2)
  endif


  ! Tabular.  Time units in sfhtab are assumed to be linear years of lookback time.
  if (pset%sfh.eq.2.or.pset%sfh.eq.3) then
     total_weights = 0.
!     call setup_tabular()
!     ! Linearly interpolate in the bins.
!     sfhpars%type = 5
!     ! Loop over each bin.
!     do i=1,ntabsfh-1
!        ! mass formed in this bin assuming linear
!        mass = (sfhtab(i,2) + sfhtab(i+1, 2)) * (sfhtab(i+1,1) - sfhtab(i, 1)) / 2
!        ! min and max ssps to consider
!        imin = min(max(locate(time_full, log10(sfhtab(i, 1))) - 1, 0), ntfull)
!        imax = min(max(locate(time_full, log10(sfhtab(i+1, 1))) + 2, 0), ntfull)
!        ! set integration limits
!        sfhpars%tq = sfhtab(i, 1)
!        sfhpars%tage = sfhtab(i+1, 1)
!        sfhpars%sf_slope = (sfhtab(i, 2) - sfhtab(i+1, 2)) / (sfhtab(i+1, 1) - sfhtab(i, 1))
!        
!        ! get the weights for this bin in the tabulated sfh and add to the
!        ! total weight, after normalizing
!        w1 = sfh_weight(sfhpars, imin, imax)
!        m1 = sum(w1)
!        if (m1.lt.tiny_number) m1 = 1.0
!        total_weights = total_weights + w1 * (mass / m1)
!     enddo
     imin = 0
     imax = ntfull
  endif

  ! Now weight each SSP by `total_weight` and sum.
  ! This matrix multiply could probably be optimized!!!!
  spec_csp = 0.
  do j=max(imin, 1), imax
     if (total_weights(j).gt.tiny_number) then
        spec_csp = spec_csp + total_weights(j) * spec_ssp(:, j)
     endif
  enddo
  mass_csp = sum(mass_ssp * total_weights)
  lbol_csp = sum(lbol_ssp * total_weights)
  mdust_csp = sum(mdust_ssp * total_weights)

end subroutine csp_gen


subroutine convert_sfhparams(pset, tage, sfh)
  ! Convert the pset values to yrs, subtract sf_start from the relevant times,
  ! pre-calculate some useful lookback times, and store in a structure.  Note
  ! that this subroutine uses but does not alter the pset values.
  !
  ! Inputs
  ! ----------
  !
  ! pset:
  !    The parameter set.
  !
  ! tage:
  !    The age (in forward time, not lookback time) at which you are
  !    calculating the composite spectrum, in Gyr.
  !
  ! Outputs
  ! ----------
  ! sfh:
  !    An SFHPARAMS structure, containing SFH parameters in yrs, with sf_start
  !    subtracted, and some useful lookback times.
  !       - `tq` is the truncation time, in lookback time
  !       - `t0` is the zero crossing time for a linear SFH, in lookback time.
  !       - `tb` is the burst time, in lookback time.
  !
  use sps_vars, only: tiny_number, SFHPARAMS, PARAMS, SP
  implicit none
  
  type(PARAMS), intent(in) :: pset
  real(SP), intent(in) :: tage
  
  type(SFHPARAMS), intent(inout) :: sfh

  real(SP) :: start
  ! Convert units from Gyr to yr, subtract sf_start
  start = pset%sf_start * 1e9
  sfh%tage = tage * 1e9 - start
  sfh%tburst = pset%tburst * 1e9 - start
  sfh%sf_trunc = pset%sf_trunc * 1e9 - start
  sfh%tau = pset%tau * 1e9
  ! Note the sign flip here!
  sfh%sf_slope = -pset%sf_slope / 1e9

  ! convert tburst to lookback time
  sfh%tb = sfh%tage - sfh%tburst
  
  ! convert sf_trunc to to lookback time
  if ((sfh%sf_trunc.le.0).or.(sfh%sf_trunc.gt.sfh%tage)) then
     sfh%tq = 0.
  else
     sfh%tq = sfh%tage - sfh%sf_trunc
  endif
  ! For simha get zero crossing time (in lookback time), avoiding divison by zero
  ! Note that only positive slopes have a chance to hit zero SFR.
  if (sfh%sf_slope.gt.tiny_number) then
     sfh%t0 = sfh%tq - 1. / sfh%sf_slope
  else
     sfh%t0 = 0.
  endif
  ! If the zero crossing time is outside the linear regime, set it to zero.
  if ((sfh%t0.gt.sfh%tq).or.(sfh%t0.le.0)) then
     sfh%t0 = 0.
  endif

end subroutine convert_sfhparams
