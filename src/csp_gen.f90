subroutine csp_gen(mass_ssp, lbol_ssp, spec_ssp, &
     pset, tage, nzin, mass_csp, lbol_csp, spec_csp, &
     mdust_csp,emlin_ssp,emlin_csp)
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
  !   The age (in Gyr of forward time) at which the spectrum is desired.  Note
  !   that this can be different than pset%tage if the latter is 0.
  !
  ! Outputs
  ! ---------
  !
  ! mass_csp, lbol_csp, spec_csp:
  !   The (surviving) stellar masses, bolometric luminosity, and spectrum of
  !   the composite stellar population at tage, normalized to 1 M_sun *formed*.
  
  use sps_vars, only: ntfull, nspec, time_full, tiny_number, tiny_logt, &
                      zlegend, nz, sfh_tab, ntabsfh, compute_light_ages, &
                      SFHPARAMS, PARAMS, SP, nemline, dust_type
  use sps_utils, only: locate, sfh_weight, sfhinfo, add_dust
  implicit none

  real(SP), intent(in), dimension(ntfull, nzin) :: mass_ssp, lbol_ssp
  real(SP), intent(in), dimension(nspec, ntfull, nzin) :: spec_ssp
  type(PARAMS), intent(in) :: pset
  real(SP), intent(in) :: tage
  integer, intent(in) :: nzin
  
  real(SP), intent(out) :: mass_csp, lbol_csp, mdust_csp
  real(SP), intent(out), dimension(nspec) :: spec_csp

  real(SP), DIMENSION(nemline, ntfull, nzin), intent(in) :: emlin_ssp
  real(SP), DIMENSION(nemline), intent(out) :: emlin_csp

  real(SP), dimension(nspec) :: csp1, csp2, lw_age
  real(SP), dimension(nemline) :: ncsp1, ncsp2, nlw_age
  real(SP), dimension(ntfull, nzin) :: total_weights
  real(SP), dimension(ntfull) :: w1=0., w2=0.
  integer :: i, j, k, imin, imax, i_tesc
  type(SFHPARAMS) :: sfhpars
  real(SP) :: m1, m2, frac_linear, mfrac, sfr, fburst
  real(SP) :: t1, t2, dt, zbin, dz  ! for tabular calculations
  real(SP) :: lbol_age, mass_age ! for mass and lbol weighted ages

  ! ------- Setup ----------

  ! Build a structure containing useful units, numbers, and switches for the
  ! weight calculations. The units of the parameters in the `sfhparams`
  ! structure are years of lookback time.
  call convert_sfhparams(pset, tage, sfhpars)
  ! Only calculate SFH weights for SSPs up to tage
  ! (plus the next couple, to bracket and be safe).
  imin = 0
  imax = min(max(locate(time_full, log10(sfhpars%tage)) + 2, 1), ntfull)

  ! ----- Get SFH weights -----

  total_weights = 0.

  ! SSP.
  if (pset%sfh.eq.0) then
     ! Make sure to use SSP weighting scheme
     sfhpars%type = -1
     ! Use tage as the burst lookback time, instead of tage-tburst.
     sfhpars%tb = sfhpars%tage
     ! Only need weights at two SSP points. Though in practice this doesn't
     ! matter, as the appropriate ages are located within sfh_weight.
     ! But it speeds up the matrix multiply later
     imin = max(imax - 2, 1)
     ! These come pre-normalized to 1 Msun
     total_weights(:, 1) = sfh_weight(sfhpars, imin, imax)
  endif


  ! Tau and delayed-tau.
  if ((pset%sfh.eq.1).or.(pset%sfh.eq.4)) then
     imin = 0
     sfhpars%type = pset%sfh
     total_weights(:, 1) = sfh_weight(sfhpars, imin, imax)
     ! Could save some loops by having proper normalization analytically from sfh_weight
     m1 = sum(total_weights(1:imax, 1))
     if (m1.lt.tiny_number) m1 = 1.0
     total_weights(:, 1) = total_weights(:, 1) / m1
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
     ! If burst happens after age of system then we kill it entirely.
     sfhpars%type = -1
     if (sfhpars%tb.lt.0) then
        w2 = 0.
        fburst = 0.
     else
        w2 = sfh_weight(sfhpars, imin, imax)
        fburst = pset%fburst
        ! We'll need to add any early bursts when summing later
        imax = max(imax, min(max(locate(time_full, log10(sfhpars%tb)) + 2, 1), ntfull))
     endif
     ! Sum with proper relative normalization.  Beware divide by zero.
     if (m1.lt.tiny_number) m1 = 1.0
     total_weights(:,1) = (1 - pset%const - fburst) * total_weights(:,1) + &
                          pset%const * (w1 / m1) + &
                          fburst * w2
  endif

  
  ! Simha
  if (pset%sfh.eq.5) then
     imin = 0
     ! Delayed-tau model portion, could set imin here to be just before
     ! sfhpars%tq, for small speed increase.
     ! imin = min(max(locate(time_full, log10(sfhpars%tq)), 0), ntfull-1)
     sfhpars%type = 4
     w1 = sfh_weight(sfhpars, imin, imax)
     m1 = sum(w1(1:imax))
     ! Linear portion.  Need to set `use_simha_limits` flag to get correct
     ! integration limits. Could set imax here to be just after sfhpars%tq,
     ! but with imin=0, for small speed increase.
     ! imax = imin + 1
     ! imin = 0
     sfhpars%type = 5
     sfhpars%use_simha_limits = 1
     w2 = sfh_weight(sfhpars, imin, imax)
     sfhpars%use_simha_limits = 0
     m2 = sum(w2(1:imax))
     ! Normalize and sum.  need to be careful of divide by zero here, if all
     ! linear or all delay-tau s.t. weights sum to zero for a component.
     if (m1.lt.tiny_number) m1 = 1.0
     if (m2.lt.tiny_number) m2 = 1.0
     ! need to get the fraction of the mass formed in the linear portion
     call sfhinfo(pset, tage, mfrac, sfr, frac_linear)
     total_weights(:, 1) = (w1 / m1) * (1 - frac_linear) + frac_linear * (w2 / m2)
     ! imax = min(max(locate(time_full, log10(sfhpars%tage)) + 2, 1), ntfull)
  endif


  ! Tabular.  Time units in sfh_tab are assumed to be linear years of time
  ! since big bang (forward time).  We are going to treat this as a sum of
  ! linear SFHs, one for each bin in the table
  if (pset%sfh.eq.2.or.pset%sfh.eq.3) then
     total_weights = 0.

     ! Assume linear SFH within the bins
     sfhpars%type = 5
     ! Loop over each bin in the table.
     do j=1, ntabsfh-1
        ! Edges of the bin in lookback time. Note that the order of sfhtab gets
        ! flipped, since it is given in forward time and then we convert to
        ! lookback time.  So j=0 is the `oldest` in terms of lookback time
        t1 = tage*1e9 - sfh_tab(1, j+1)
        t2 = tage*1e9 - sfh_tab(1, j)
        if (t2.lt.0) then
           ! Entire bin is in the future, skip
           cycle
        endif
        ! Metallicity of the bin. Just a straight average.
        zbin = (sfh_tab(3, j) + sfh_tab(3, j+1)) / 2

        ! Linear slope.  Positive should be sfr *decreasing* in time since big bang.
        sfhpars%sf_slope = -(sfh_tab(2, j+1) - sfh_tab(2, j)) / (t2 - t1) / sfh_tab(2, j+1)
        ! Set integration limits using bin edges clipped to valid times.
        ! That is, don't include any portion of a bin that goes to negative
        ! time, or beyond the oldest isochrone.
        sfhpars%tq = min(max(t1, 10**tiny_logt), 10**time_full(ntfull))  ! lower limit (in lookback time)
        sfhpars%tage = min(max(t2, 10**tiny_logt), 10**time_full(ntfull)) ! upper limit
        sfhpars%sf_trunc = sfhpars%tage - sfhpars%tq
        ! Mass that formed within these valid times.
        dt = (sfhpars%tage - sfhpars%tq)
        m2 = sfh_tab(2, j+1) * (1 + sfhpars%sf_slope/2. * (sfhpars%tage + sfhpars%tq - 2*t1)) * dt
        ! min and max ssps to consider, being conservative.
        imin = min(max(locate(time_full, log10(t1)) - 1, 0), ntfull)
        imax = min(max(locate(time_full, log10(t2)) + 2, 0), ntfull)

        ! Get the weights for this bin in the tabulated sfh and add to the
        ! total weight, after normalizing.
        w1 = sfh_weight(sfhpars, imin, imax)
        m1 = sum(w1)
        if (m1.lt.tiny_number) m1 = 1.0
        ! This is where we'd assign to specific metallicities, if taking that
        ! into account.  This scheme assumes entire bin is at average of the
        ! two enclosing Z values.
        if (nzin.gt.1) then
           k = max(min(locate(zlegend, zbin), nz-1), 1)
           dz = (log10(zbin) - log10(zlegend(k))) / &
                (log10(zlegend(k+1)) - log10(zlegend(k)))
           dz = max(min(dz, 1.0), -1.0) !don't extrapolate
           total_weights(:, k) = total_weights(:, k) + (1-dz) * w1 * (m2 / m1)
           total_weights(:, k+1) = total_weights(:, k+1) + dz * w1 * (m2 / m1)
        else
           total_weights(:, 1) = total_weights(:, 1) + w1 * (m2 / m1)
        endif
     enddo
     ! Reset imin and imax for the spectral sum.
     imin = 0
     imax = ntfull
  endif

  ! ----- Combine SSPs with dust -------
  ! Now weight each SSP by `total_weight`, assign to young or old, and feed to
  ! add_dust, which does the final sum as well as attenuating.
  ! This matrix multiply could probably be optimized!!!!

  csp1  = 0. ! young
  csp2  = 0. ! old
  ncsp1 = 0. ! young nebular
  ncsp2 = 0. ! old nebular
  if (compute_light_ages.eq.1) then
     nlw_age  = 0.
     lw_age   = 0.
     lbol_age = 0.
     mass_age = 0.
  endif

  ! Dust treatment is not strictly correct, since the age bin older than
  ! dust_tesc will include some contribution from young star dust due to the
  ! interpolation, and changing dust_tesc by values smaller than half the ssp
  ! age grid resolution will have no effect on the output.
  i_tesc = locate(time_full, pset%dust_tesc)
  do i=max(imin, 1), imax
     do k=1,nzin
        if (total_weights(i, k).gt.tiny_number) then
           if (i.le.i_tesc) then
              csp1  = csp1  + total_weights(i, k) * spec_ssp(:, i, k)
              ncsp1 = ncsp1 + total_weights(i, k) * emlin_ssp(:, i, k)
           else
              csp2  = csp2  + total_weights(i, k) * spec_ssp(:, i, k)
              ncsp2 = ncsp2 + total_weights(i, k) * emlin_ssp(:, i, k)
           endif
           ! Now do numerator in case of light and mass weighted ages.
           if (compute_light_ages.eq.1) then
              nlw_age = nlw_age + total_weights(i, k) * emlin_ssp(:, i, k) &
                       * 10**(time_full(i)-9)
              lw_age = lw_age + total_weights(i, k) * spec_ssp(:, i, k) &
                       * 10**(time_full(i)-9)
              lbol_age = lbol_age + 10**lbol_ssp(i, k) * total_weights(i, k) &
                         * 10**(time_full(i)-9)
              mass_age = mass_age + mass_ssp(i, k) * total_weights(i, k) &
                         * 10**(time_full(i)-9)
           endif
        endif
     enddo
  enddo
  mass_csp = sum(mass_ssp * total_weights)
  lbol_csp = log10(sum(10**lbol_ssp * total_weights))

  ! Here we add young and old spectra with dust.
  if (((pset%dust1.gt.tiny_number).or.(pset%dust2.gt.tiny_number).or.(dust_type.eq.3))&
       .and.(compute_light_ages.eq.0)) then
     call add_dust(pset, csp1, csp2, spec_csp, mdust_csp, ncsp1, ncsp2, emlin_csp)

  else
     spec_csp  = csp1 + csp2
     emlin_csp = ncsp1 + ncsp2
     mdust_csp = 0.0
  endif

  if (compute_light_ages.eq.1) then
     spec_csp  = lw_age / spec_csp
     emlin_csp = nlw_age / emlin_csp
     lbol_csp  = lbol_age / 10**lbol_csp
     mass_csp  = mass_age / mass_csp
  endif

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
  !       - `sf_slope` is the fractional change in the SFR in inverse years.  It is
  !          positive for SFR that increases with *lookback* time.
  !
  use sps_vars, only: tiny_number, SFHPARAMS, PARAMS, SP
  implicit none
  
  type(PARAMS), intent(in) :: pset
  real(SP), intent(in) :: tage
  
  type(SFHPARAMS), intent(inout) :: sfh

  real(SP) :: start

  ! Define a starting time iff SFH=1,4,5
  if ((pset%sfh.eq.1).or.(pset%sfh.eq.4).or.(pset%sfh.eq.5)) then
     start = pset%sf_start * 1e9
  else
     start = 0.
  endif
  ! Convert units from Gyr to yr, subtract sf_start
  sfh%tage = tage * 1e9 - start
  sfh%tburst = pset%tburst * 1e9 - start
  sfh%sf_trunc = pset%sf_trunc * 1e9 - start
  sfh%tau = pset%tau * 1e9
  ! Note the sign flip here!  pset%sf_slope is positive for sfr that increases
  ! with forward time, sfh%sf_slope is positive for sfr that increases with
  ! *lookback time*.  Yeah, awesome!
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
