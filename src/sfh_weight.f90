function sfh_weight(sfh, imin, imax)
  ! Function to calculate the weights
  !
  ! Inputs
  ! --------
  !
  ! sfh:
  !   An SFHPARAMS structure, containing SFH parameters in yrs, and some
  !   useful lookback times.
  !
  ! imin:
  !   The index of the youngest SSP to consider.  If zero, weights for the t=0
  !   -> t=time_full(1) bin will be added to sfh_weights(1)
  !
  ! imax:
  !   The oldest SSP to consider.  This should usually be the next oldest SSP
  !   compared to tage (or tage - sf_start)
  !
  ! Outputs
  ! --------
  !
  ! sfh_weight:
  !   array of shape `ntfull` that gives the proper weights for the SSPs to
  !   produce the SFH.

  use sps_vars, only: ntfull, time_full, tiny_logt, &
                      SFHPARAMS, SP
  use sps_utils, only: intsfwght, sfhlimit, locate
  implicit none
  
  type(SFHPARAMS), intent(in) :: sfh
  integer, intent(in) :: imin, imax

  real(SP), dimension(ntfull) ::sfh_weight
  
  integer :: i, istart
  real(SP), dimension(2) :: tlim
  real(SP) :: dt, delta_time, log_tb
  real(SP), dimension(ntfull) :: tmp_wght=0. !left=0., right=0.
  

  ! Check if this is an SSP.  If so, do simple weights and return.
  if (sfh%type.eq.-1) then
     sfh_weight = 0.
     ! Check for burst out of valid limits.
     ! N.B. we allow bursts before sf_start (sfh%tb > sfh%tage)
     if (sfh%tb.lt.0) then
        return
     endif
     log_tb = log10(sfh%tb)
     if (log_tb.le.time_full(1)) then
        ! Just return all weight in the youngest ssp.
        sfh_weight(1) = 1.0
        return
     endif

     istart = min(max(locate(time_full, log_tb), 1), ntfull-1)
     dt = delta_time(time_full(istart), time_full(istart+1))
     sfh_weight(istart) = delta_time(log_tb, time_full(istart+1)) / dt
     sfh_weight(istart+1) = delta_time(time_full(istart), log_tb) / dt
     return
  endif

  ! Loop over each SSP and calculate its weight in the given sfh.
  ! Note we skip i=0 for now, which is a flag for adding in the bin from t=0 to
  ! t=time_full(1).
  tmp_wght = 0.
  istart = max(imin, 1)
  do i=istart, imax
     if (i.gt.1) then
        ! There is a younger (`left`) bin, and we calculate its contribution to
        ! the weight.
        ! First calculate actual limits for the younger bin.  
        tlim(1) = sfhlimit(time_full(i-1), sfh)
        tlim(2) = sfhlimit(time_full(i), sfh)
        ! The elements of `tlim` will be equal if there is no valid SFR in the
        ! younger bin; only proceed if there is a non-zero sfr in the younger bin.
        if (tlim(1).ne.tlim(2)) then
           dt = delta_time(time_full(i-1), time_full(i))
           ! Note sign flip here
           !left(i) = 0. - intsfwght(i-1, tlim, sfh) / dt
           tmp_wght(i) = tmp_wght(i) - intsfwght(i-1, tlim, sfh) / dt
        endif
     endif
     if (i.lt.ntfull) then
        ! There is an older (`right`) bin, we calculate its contribution to the
        ! weight.
        tlim(1) = sfhlimit(time_full(i), sfh)
        tlim(2) = sfhlimit(time_full(i+1), sfh)
        ! The elements of `tlim` will be equal if there is no valid SFR in the
        ! older bin; only proceed if there is a non-zero sfr in the older bin.
        if (tlim(1).ne.tlim(2)) then
           dt = delta_time(time_full(i), time_full(i+1))
           !right(i) = intsfwght(i+1, tlim, sfh) / dt
           tmp_wght(i) = tmp_wght(i) + intsfwght(i+1, tlim, sfh) / dt
        endif
     endif
  enddo

  ! Do we need to add weights from the zero bin to the first SSP?
  !   We assume the t~0 spectrum is the same as the t=10**time_full(1) spectrum
  !   (i.e., nearest neighbor extrapolation), so the t~0 weight gets added to
  !   sfh_weight(1)
  if (imin.eq.0) then
     tlim(1) = sfhlimit(tiny_logt, sfh)
     tlim(2) = sfhlimit(time_full(1), sfh)
     if (tlim(1).ne.tlim(2)) then
        dt = delta_time(tiny_logt, time_full(1))
        ! Contribution of i=1 to younger bin.
        tmp_wght(1) = tmp_wght(1) - intsfwght(0, tlim, sfh) / dt
        ! Contribution of i=0 to older bin
        tmp_wght(1) = tmp_wght(1) + intsfwght(1, tlim, sfh) / dt
     endif
  endif

  sfh_weight = tmp_wght

end function sfh_weight

  
function delta_time(logt1, logt2)
  ! Dumb function to properly calculate dt based on interpolation type.
  !
  ! Returns (logt2 - logt1), or (10**logt2 - 10**logt1)

  use sps_vars, only: interpolation_type, SP
  implicit none

  real(SP), intent(in) :: logt1, logt2
  real(SP) :: delta_time

  if (interpolation_type.eq.0) then
     delta_time = logt2 - logt1
  else if (interpolation_type.eq.1) then
     delta_time = 10**logt2 - 10**logt1
  endif

end function delta_time
