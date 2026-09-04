!----------------------------------------------------------------------
! check_solar_calibration (formerly chkcal)
!----------------------------------------------------------------------
! Checks whether the current model is a converged solar calibration
! (log L and log R, and optionally log Z/X, within tolerance of the
! solar values) and, if not, applies a fixed 2x2 Newton correction to
! the envelope hydrogen fraction (X) and mixing-length alpha using
! pre-determined empirical partial derivatives of log L and log R with
! respect to X and alpha (star%dlum_dx/star%drad_dx/star%dlum_dalpha/star%drad_dalpha).
! Called after every completed calibration TRIPLE (mod(nk,3) == 0;
! the protocol comment in run_yrec's run loop is the reference):
! card 3k+1 rescales, 3k+2 settles to 1e8 yr, 3k+3 evolves to the
! solar age. On a miss, the corrected (X[, Z], alpha) are written
! into the NEXT triple's cards (run_index+1 .. +3), reading the
! previous guess from run_index-2 (this triple's rescale card).
subroutine check_solar_calibration(log_l_lsun, log_r_rsun, run_index, current_zx)

      use star_info_lib, only: star
      use luout_lib
      use math_lib
      implicit none

      double precision, intent(in) :: log_l_lsun, log_r_rsun
      integer, intent(in) :: run_index
      double precision, intent(in) :: current_zx

! --- locals ---
      double precision :: log_zx_mismatch, delta_z

! CHECK TO SEE IF THE MODEL IS A CALIBRATED SOLAR MODEL.
! IF NOT, ESTIMATE CORRECTIONS TO X AND ALPHA.
      if(abs(log_l_lsun).lt.star%ctrl%luminosity_tolerance .and. &
           abs(log_r_rsun).lt.star%ctrl%radius_tolerance)then
! ADD CHECK FOR Z/X
         if(star%ctrl%calibrate_solar_zx .and.star%ctrl%target_solar_zx.gt.0.0d0)then
            log_zx_mismatch = log10(current_zx)-log10(star%ctrl%target_solar_zx)
            if(abs(log_zx_mismatch).lt.star%ctrl%zx_tolerance)then
               star%solar_calibration_active = .true.
               return
            endif
         else
! CALIBRATED SOLAR MODEL.  SET UP OUTPUT FLAGS AND EXIT
            star%solar_calibration_active = .true.
            return
         endif
      endif

! Added code to use pre-determined partial derivatives    LLP  6/28/09
      star%dlum_dx = -3.78d0            ! empirical result:  -3.783    RMS error .070
      star%drad_dx = -0.89d0             ! empirical result:  -0.890    RMS error .048
      star%dlum_dalpha = 0.0139d0            ! empirical result:  +0.139    RMS error .0022
      star%drad_dalpha = -0.050d0            ! empirical result:  -0.0504   RMS error .0059
! (2026 goto campaign: an unconditional `goto 1234` here bypassed the
! original run_index-based partial-derivative code, which was therefore
! dead; it was removed. See git history for the bypassed block.)
!     USE DERIVATIVES OF L AND R WITH RESPECT TO X AND ALPHA TO
!     GET IMPROVED GUESSES FOR ALPHA AND X.
         star%delta_alpha = ((log_l_lsun*star%drad_dx/star%dlum_dx-log_r_rsun)/ &
              (star%drad_dalpha-star%dlum_dalpha*star%drad_dx/star%dlum_dx))
         star%delta_x = -(log_l_lsun + star%dlum_dalpha*star%delta_alpha)/star%dlum_dx
         star%job%mixing_length_array(run_index+1) = &
              star%job%mixing_length_array(run_index-2)+star%delta_alpha
         star%job%mixing_length_array(run_index+2) = star%job%mixing_length_array(run_index+1)
         star%job%mixing_length_array(run_index+3) = star%job%mixing_length_array(run_index+1)
         if(star%ctrl%calibrate_solar_zx)then
            star%job%rescale_params(3,run_index+1) = &
                 star%job%rescale_params(3,run_index-2)*star%ctrl%target_solar_zx/current_zx
            delta_z = star%job%rescale_params(3,run_index+1) - &
                 star%job%rescale_params(3,run_index-2)
            star%job%initial_z_array(run_index+1) = star%job%rescale_params(3,run_index+1)
            star%job%initial_z_array(run_index+2) = star%job%rescale_params(3,run_index+1)
            star%job%initial_z_array(run_index+3) = star%job%rescale_params(3,run_index+1)
         endif
         star%job%rescale_params(2,run_index+1) = star%job%rescale_params(2,run_index-2)+star%delta_x
         star%job%initial_x_array(run_index+1) = star%job%rescale_params(2,run_index+1)
         star%job%initial_x_array(run_index+2) = star%job%rescale_params(2,run_index+1)
         star%job%initial_x_array(run_index+3) = star%job%rescale_params(2,run_index+1)
         write(terminal_unit,*) "New BL, Old BL, Delta BL: ", &
             log_l_lsun, star%log_l_prev, log_l_lsun-star%log_l_prev
         write(terminal_unit,*) "New RL, Old RL, Delta RL: ", &
             log_r_rsun, star%log_r_prev, log_r_rsun-star%log_r_prev
         write(terminal_unit,*) "New X, Old X, DX: ", &
             star%job%rescale_params(2,run_index+1), star%job%rescale_params(2,run_index-2), &
             star%delta_x
         write(terminal_unit,*) "New A, Old A, DA: ", &
             star%job%mixing_length_array(run_index+1), &
             star%job%mixing_length_array(run_index-2), star%delta_alpha
         if(star%ctrl%calibrate_solar_zx)then
            write(terminal_unit,*) "New Z, Old Z, DZ: ", &
             star%job%rescale_params(3,run_index+1), star%job%rescale_params(3,run_index-2), &
             delta_z
         endif
         star%log_l_prev = log_l_lsun
         star%log_r_prev = log_r_rsun
end subroutine check_solar_calibration
