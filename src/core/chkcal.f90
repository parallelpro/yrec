!----------------------------------------------------------------------
! chkcal
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original chkcal.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
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
! (An older header here claimed the pre-5/96 even/odd pair protocol;
! fixed 2026.)
subroutine chkcal(log_l_lsun, log_r_rsun, run_index, current_zx)

      use star_info_lib, only: star
      use luout_lib
      implicit none

      double precision, intent(in) :: log_l_lsun, log_r_rsun
      integer, intent(in) :: run_index
      double precision, intent(in) :: current_zx


!      COMMON/SETT/ENDAGE(50),SETDT(50),LENDAG(50),LSETDT(50)



!      COMMON/CALS2/TOLL,TOLR,LCALS
!     DATA TOLL,TOLR/1.0D-5,1.0D-4/

! locals -- log_zx_mismatch/delta_z/log_zx_mismatch_prev (originally
! ZL/DZ/ZLP); ZL is set only in the early-return Z/X check below but
! (via the blanket SAVE, matching the original) is still read again
! from its previous-call value at the "New Z" report near label 1234.
      double precision :: log_zx_mismatch, delta_z, log_zx_mismatch_prev

! CHECK TO SEE IF THE MODEL IS A CALIBRATED SOLAR MODEL.
! IF NOT, ESTIMATE CORRECTIONS TO X AND ALPHA.
      if(abs(log_l_lsun).lt.star%ctrl%luminosity_tolerance .and. &
           abs(log_r_rsun).lt.star%ctrl%radius_tolerance)then
! ADD CHECK FOR Z/X
         if(star%ctrl%calibrate_solar_zx .and.star%ctrl%target_solar_zx.gt.0.0d0)then
            log_zx_mismatch = log10(current_zx)-log10(star%ctrl%target_solar_zx)
            if(abs(log_zx_mismatch).lt.star%ctrl%zx_tolerance)then
               star%solar_calibration_active = .true.
               continue
               return
            endif
         else
! CALIBRATED SOLAR MODEL.  SET UP OUTPUT FLAGS AND EXIT
            star%solar_calibration_active = .true.
            continue
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
         write(iowr,*) "New BL, Old BL, Delta BL: ", &
             log_l_lsun, star%log_l_prev, log_l_lsun-star%log_l_prev
         write(iowr,*) "New RL, Old RL, Delta RL: ", &
             log_r_rsun, star%log_r_prev, log_r_rsun-star%log_r_prev
         write(iowr,*) "New X, Old X, DX: ", &
             star%job%rescale_params(2,run_index+1), star%job%rescale_params(2,run_index-2), &
             star%delta_x
         write(iowr,*) "New A, Old A, DA: ", &
             star%job%mixing_length_array(run_index+1), &
             star%job%mixing_length_array(run_index-2), star%delta_alpha
         if(star%ctrl%calibrate_solar_zx)then
            write(iowr,*) "New Z, Old Z, DZ: ", &
             star%job%rescale_params(3,run_index+1), star%job%rescale_params(3,run_index-2), &
             delta_z
         endif
         star%log_l_prev = log_l_lsun
         star%log_r_prev = log_r_rsun
         log_zx_mismatch_prev = log_zx_mismatch
         continue
         return
!      ENDIF
      return
end subroutine chkcal
