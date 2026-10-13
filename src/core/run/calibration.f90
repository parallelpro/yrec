!----------------------------------------------------------------------
! calibration
!----------------------------------------------------------------------
! The two halves of the calibration loops: the run-list setup before a
! calibration (setup_solar_calibration, setup_star_calibration) and the
! check after each run that declares convergence or sets up the next
! trial (check_solar_calibration, check_star_calibration).
module calibration
      use star_info_lib, only: star_info
      implicit none
      private
      public :: check_solar_calibration, check_star_calibration, setup_solar_calibration, setup_star_calibration

contains

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! check_solar_calibration (formerly chkcal)
!----------------------------------------------------------------------
! Tests whether the model at the end of a calibration triple is a
! converged solar model: |log(L/Lsun)| < luminosity_tolerance and
! |log(R/Rsun)| < radius_tolerance, and, when calibrate_solar_zx is
! set with a positive target_solar_zx, also
! |log(Z/X) - log(Z/X)_target| < zx_tolerance. If so, it sets
! solar_calibration_active and returns. If not, it takes one Newton
! step in (X, alpha) with fixed empirical partial derivatives of log L
! and log R (star%dlum_dx, drad_dx, dlum_dalpha, drad_dalpha), rescales
! Z by the Z/X mismatch when calibrate_solar_zx is set, and writes the
! new envelope (X, Z, alpha) into the next triple's cards.
!
! Called after every completed triple (mod(nk,3) == 0; the protocol
! comment in run_yrec's run loop is the reference): card 3k+1
! rescales, 3k+2 settles to 1e8 yr, 3k+3 evolves to the solar age.
! The current guess is read from run_index-2 (this triple's rescale
! card), and the corrected one is written into run_index+1 .. +3.
!
! Inputs: log_l_lsun = log(L/Lsun) and log_r_rsun = log(R/Rsun) of the
! model (both zero for the Sun), run_index = nk, and current_zx, the
! surface Z/X.
subroutine check_solar_calibration(s, log_l_lsun, log_r_rsun, run_index, current_zx)

      use io_units
      use math_lib
      implicit none
      type(star_info), intent(inout) :: s

      double precision, intent(in) :: log_l_lsun, log_r_rsun
      integer, intent(in) :: run_index
      double precision, intent(in) :: current_zx

! --- locals ---
      double precision :: log_zx_mismatch, delta_z

! Converged when L and R match the Sun (and Z/X its target, when that
! is calibrated too).
      if(abs(log_l_lsun).lt.s%ctrl%luminosity_tolerance .and. &
           abs(log_r_rsun).lt.s%ctrl%radius_tolerance)then
! The Z/X test.
         if(s%ctrl%calibrate_solar_zx .and.s%ctrl%target_solar_zx.gt.0.0d0)then
            log_zx_mismatch = log10(current_zx)-log10(s%ctrl%target_solar_zx)
            if(abs(log_zx_mismatch).lt.s%ctrl%zx_tolerance)then
               s%solar_calibration_converged = .true.
               return
            endif
         else
! A calibrated solar model: flag it and return.
            s%solar_calibration_converged = .true.
            return
         endif
      endif

! Fixed empirical partial derivatives d log L/dX, d log R/dX,
! d log L/d alpha and d log R/d alpha (LLP 6/28/09). The trailing
! comments give the fitted values and their RMS errors; the fit
! recorded for dlum_dalpha, +0.139, is ten times the value used.
      s%dlum_dx = -3.78d0            ! fit: -3.783, RMS error .070
      s%drad_dx = -0.89d0             ! fit: -0.890, RMS error .048
      s%dlum_dalpha = 0.0139d0            ! fit as recorded: +0.139, RMS error .0022
      s%drad_dalpha = -0.050d0            ! fit: -0.0504, RMS error .0059
! (The original code also had a run_index-based estimate of these
! derivatives, which an unconditional jump always bypassed; it has been
! removed.)
!
! The Newton step: the corrections that zero log L and log R solve
!   dlum_dx dX + dlum_dalpha dalpha = -log L
!   drad_dx dX + drad_dalpha dalpha = -log R,
! which gives
!   dalpha = (log L drad_dx/dlum_dx - log R)
!            / (drad_dalpha - dlum_dalpha drad_dx/dlum_dx),
!   dX     = -(log L + dlum_dalpha dalpha)/dlum_dx.
         s%delta_alpha = ((log_l_lsun*s%drad_dx/s%dlum_dx-log_r_rsun)/ &
              (s%drad_dalpha-s%dlum_dalpha*s%drad_dx/s%dlum_dx))
         s%delta_x = -(log_l_lsun + s%dlum_dalpha*s%delta_alpha)/s%dlum_dx
         s%job%mixing_length_array(run_index+1) = &
              s%job%mixing_length_array(run_index-2)+s%delta_alpha
         s%job%mixing_length_array(run_index+2) = s%job%mixing_length_array(run_index+1)
         s%job%mixing_length_array(run_index+3) = s%job%mixing_length_array(run_index+1)
! Z/X calibration: Z_new = Z_old (Z/X)_target / (Z/X)_current.
         if(s%ctrl%calibrate_solar_zx)then
            s%job%rescale_params(3,run_index+1) = &
                 s%job%rescale_params(3,run_index-2)*s%ctrl%target_solar_zx/current_zx
            delta_z = s%job%rescale_params(3,run_index+1) - &
                 s%job%rescale_params(3,run_index-2)
            s%job%initial_z_array(run_index+1) = s%job%rescale_params(3,run_index+1)
            s%job%initial_z_array(run_index+2) = s%job%rescale_params(3,run_index+1)
            s%job%initial_z_array(run_index+3) = s%job%rescale_params(3,run_index+1)
         endif
         s%job%rescale_params(2,run_index+1) = s%job%rescale_params(2,run_index-2)+s%delta_x
         s%job%initial_x_array(run_index+1) = s%job%rescale_params(2,run_index+1)
         s%job%initial_x_array(run_index+2) = s%job%rescale_params(2,run_index+1)
         s%job%initial_x_array(run_index+3) = s%job%rescale_params(2,run_index+1)
         write(terminal_unit,*) "New BL, Old BL, Delta BL: ", &
             log_l_lsun, s%log_l_prev, log_l_lsun-s%log_l_prev
         write(terminal_unit,*) "New RL, Old RL, Delta RL: ", &
             log_r_rsun, s%log_r_prev, log_r_rsun-s%log_r_prev
         write(terminal_unit,*) "New X, Old X, DX: ", &
             s%job%rescale_params(2,run_index+1), s%job%rescale_params(2,run_index-2), &
             s%delta_x
         write(terminal_unit,*) "New A, Old A, DA: ", &
             s%job%mixing_length_array(run_index+1), &
             s%job%mixing_length_array(run_index-2), s%delta_alpha
         if(s%ctrl%calibrate_solar_zx)then
            write(terminal_unit,*) "New Z, Old Z, DZ: ", &
             s%job%rescale_params(3,run_index+1), s%job%rescale_params(3,run_index-2), &
             delta_z
         endif
         s%log_l_prev = log_l_lsun
         s%log_r_prev = log_r_rsun
end subroutine check_solar_calibration

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! check_star_calibration (formerly chkscal)
!----------------------------------------------------------------------
! The star-calibration check, called for every model after the first
! of an even card (the odd cards rescale, the even cards evolve; see
! check_stop_conditions). It tests whether the track has just crossed
! the target radius R* (in either direction). If not, it stores the
! model's L, age and R and returns. If so, it interpolates log L and the
! age to R*, linearly in R, and then either
!  - finds |L(R*) - L*| <= target_star_luminosity_tolerance (in Lsun,
!    a linear tolerance), sets star_found_flag, and arms one more run,
!    stopped at the interpolated age, in the next pair of cards; or
!  - sets X for the next pair: on the first pass (nk = 2), X - 0.01
!    (Y + 0.01 at fixed Z); after that, a secant step in X against
!    log L(R*).
subroutine check_star_calibration(s, log_l_lsun, log_teff, current_age, run_index)

      use io_units
      use physical_constants
      use math_lib
      implicit none
      type(star_info), intent(inout) :: s

      double precision, intent(in) :: log_l_lsun, log_teff, current_age
      integer, intent(in) :: run_index

! locals
      double precision :: teff_current, r_rsun_current
      double precision :: dlogl_dlogr, dage_dlogr
      double precision :: new_x, prev_x, dx_dlogl
! The step in X that seeds the second calibration run (the "Y=Y+0.01"
! of the messages below).
      double precision, parameter :: x_calibration_step = 0.01d0

! The quantities involved:
!   star%star_found_flag                  true once a track passes
!                                         through both R* and L*
!   star%just_passed_target_radius_flag   true if this model has just
!                                         crossed R*
!   ctrl%target_luminosity_lsun           L*/Lsun, reached by adjusting X
!   ctrl%target_star_luminosity_tolerance the tolerance on L/Lsun at R*
!   job%target_radius_rsun                R*/Rsun (linear, despite the
!                                         log_r_* names of the model
!                                         values)
!   r_rsun_current                        R/Rsun of this model (linear)
!   star%r_rsun_prev_model                 R/Rsun of the previous model
!                                         (linear)
!   current_age, star%age_prev_model      age of this and of the
!                                         previous model (Gyr)
!   star%age_at_target_radius             age at R* (Gyr)
!   log_l_lsun, star%log_l_prev_model     log(L/Lsun) of this and of the
!                                         previous model
!   star%log_l_at_target_radius           log(L/Lsun) at R*
!   star%log_l_at_target_radius_prev_run  the same, from the previous run
!   prev_x                                X of the previous run,
!                                         rescale_params(2, nk-1)
!
! Has the track crossed R*? R/Rsun = sqrt(L/(4 pi sigma))/(Teff^2 Rsun).
! If it has not, store this model's L, age and R and return.
      s%just_passed_target_radius_flag=.false.
      teff_current = exp10(log_teff)
      r_rsun_current = sqrt((exp10(log_l_lsun))*s%solar_luminosity_cgs/ &
           (c4pi*csig))/(teff_current*teff_current*s%solar_radius_cgs)
      if(r_rsun_current.gt.s%r_rsun_prev_model) then
         if(.not.(r_rsun_current.gt.s%job%target_radius_rsun.and. &
              s%r_rsun_prev_model.le.s%job%target_radius_rsun)) then
              s%log_l_prev_model = log_l_lsun
              s%age_prev_model = current_age
              s%r_rsun_prev_model = r_rsun_current
              return
          end if
      else
          if (.not.(r_rsun_current.lt.s%job%target_radius_rsun.and. &
               s%r_rsun_prev_model.ge.s%job%target_radius_rsun)) then
              s%log_l_prev_model = log_l_lsun
              s%age_prev_model = current_age
              s%r_rsun_prev_model = r_rsun_current
              return
           end if
      endif
!
      write(*,*) ' Just passed R*'
      write(run_log_unit,*) '#Just passed R*'
      s%just_passed_target_radius_flag = .true.
!
! The previous model lies on one side of R* and this one on the other:
! interpolate log L and the age to R*, linearly in R/Rsun,
!   log L(R*) = log L + (d log L/dR) (R* - R),
!   age(R*)   = age + (d age/dR) (R* - R).
! (Despite their names, dlogl_dlogr and dage_dlogr are derivatives with
! respect to the linear R/Rsun.)
      dlogl_dlogr = (log_l_lsun-s%log_l_prev_model)/ &
           (r_rsun_current-s%r_rsun_prev_model)
      s%log_l_at_target_radius = log_l_lsun + &
           dlogl_dlogr*(s%job%target_radius_rsun-r_rsun_current)
      dage_dlogr = (current_age-s%age_prev_model)/ &
           (r_rsun_current-s%r_rsun_prev_model)
      s%age_at_target_radius = current_age + &
           dage_dlogr*(s%job%target_radius_rsun-r_rsun_current)
      write(*,*) ' X, LogL/Lsun at R* =', s%job%rescale_params(2,run_index-1), &
           s%log_l_at_target_radius
      write(run_log_unit,*) '#X, LogL/Lsun at R* =', &
           s%job%rescale_params(2,run_index-1), s%log_l_at_target_radius
      if (abs(exp10(s%log_l_at_target_radius)-s%ctrl%target_luminosity_lsun) &
           .le. s%ctrl%target_star_luminosity_tolerance) then
! The track passes through the target L and R. Do one more run with the
! same X, stopped at the interpolated age at R* (converted from Gyr to
! yr).
         s%star_found_flag=.true.
       s%job%end_age_stop_active(run_index+1) = .true.
       s%job%target_end_age(run_index+1) = s%age_at_target_radius*1.0d9
       s%job%end_age_stop_active(run_index+2) = .true.
       s%job%target_end_age(run_index+2) = s%age_at_target_radius*1.0d9
       s%job%rescale_params(2,run_index+1) = s%job%rescale_params(2,run_index-1)
       s%job%initial_x_array(run_index+1) = s%job%rescale_params(2,run_index+1)
       s%job%initial_x_array(run_index+2) = s%job%initial_x_array(run_index+1)
!
      write(*,*) ' Have hit R* & L*, prepare final run to age:', &
        s%age_at_target_radius
      write(run_log_unit,*)'#Have hit R* & L*, prepare final run to age:', &
        s%age_at_target_radius
         return
      else
         if (run_index .eq. 2) then
! First pass: save log L at R*, lower X by x_calibration_step (Y up by
! the same at fixed Z), and start the next run.
            s%log_l_at_target_radius_prev_run = s%log_l_at_target_radius
          prev_x = s%job%rescale_params(2,run_index-1)
          new_x = prev_x - x_calibration_step
          s%job%rescale_params(2,run_index+1) = new_x
          s%job%initial_x_array(run_index+1) = new_x
          s%job%initial_x_array(run_index+2) = new_x
!
      write(*,*) ' NK=2, Y=Y+0.01, Setup next run, X=', new_x
      write(run_log_unit,*) '#NK=2, Y=Y+0.01, Setup next run, X=', new_x
            return
       else
! Later passes (nk = 4, 6, 8, ...): a secant step from this run's and
! the previous run's X and log L at R*,
!   dX/dlog L = (X_n - X_n-1)/(log L(R*)_n - log L(R*)_n-1),
!   X_new = X_n + (dX/dlog L) (log L* - log L(R*)_n),
! with X_n = rescale_params(2, nk-1) and X_n-1 = rescale_params(2, nk-3).
! Save log L(R*) and start the next run.
            new_x = s%job%rescale_params(2, run_index-1)
          prev_x = s%job%rescale_params(2, run_index-3)
            dx_dlogl = (new_x-prev_x)/ &
                 (s%log_l_at_target_radius-s%log_l_at_target_radius_prev_run)
          new_x = dx_dlogl*(log10(s%ctrl%target_luminosity_lsun)- &
               s%log_l_at_target_radius)+new_x
!
      write(*,*) ' Setup next run, NK, X =', run_index, new_x
      write(run_log_unit, *) ' Setup next run, NK, X =', run_index, new_x
          s%log_l_at_target_radius_prev_run = s%log_l_at_target_radius
          s%job%rescale_params(2,run_index+1) = new_x
          s%job%initial_x_array(run_index+1) = new_x
          s%job%initial_x_array(run_index+2) = new_x
       end if
      endif
end subroutine check_star_calibration

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! setup_solar_calibration (legacy setcal.f)
!----------------------------------------------------------------------
! Sets up the run list for an automatic solar calibration from the
! three seed cards (card 1 rescales, card 2 settles, card 3 evolves to
! the solar age):
!  - num_runs becomes max_runs-2 = 48 cards, i.e. 16 triples;
!  - every card takes card 1's X, Z, alpha and envelope-mass request;
!  - cards 4, 7, ..., 46 repeat card 1's rescaling (a fresh model read
!    with card 1's rescale parameters);
!  - cards 5, 8, ..., 47 continue the model in memory (rescale_kind 1,
!    evolve) with card 2's model count, stop age and timestep;
!  - cards 6, 9, ..., 48 do the same with card 3's.
! Card 2 is hardwired to stop at 1e8 yr and card 3 at target_solar_age,
! both times age_scale_factor (the Monte-Carlo age factor, 1 otherwise),
! and their timestep overrides are scaled by the same factor.
! check_solar_calibration writes each new (X, Z, alpha) guess into the
! next triple.
subroutine setup_solar_calibration(s, age_scale_factor)

      use controls_buffer, only: max_runs
      implicit none
      type(star_info), intent(inout) :: s

      double precision, intent(in) :: age_scale_factor
! --- locals ---
      integer :: i, j
! 48 cards: 16 three-card calibration cycles (the run arrays in
! controls_buffer have max_runs = 50 slots).
      integer, parameter :: num_calibration_runs = max_runs-2

! Card 1's composition, mixing length and envelope-mass request go to
! every card.
      s%job%num_runs = num_calibration_runs
      do i = 2,num_calibration_runs
         s%job%initial_x_array(i) = s%job%initial_x_array(1)
         s%job%initial_z_array(i) = s%job%initial_z_array(1)
         s%job%mixing_length_array(i) = s%job%mixing_length_array(1)
         s%job%has_senv0_array(i) = s%job%has_senv0_array(1)
         s%job%senv0_array(i) = s%job%senv0_array(1)
      end do
! The rescale cards repeat card 1.
      do i = 4,46,3
         s%job%rescale_kind(i) = s%job%rescale_kind(1)
         s%job%first_call_flag(i) = .true.
         s%job%num_models(i) = s%job%num_models(1)
         s%job%rsclzc(i) = s%job%rsclzc(1)
         s%job%rsclzm1(i) = s%job%rsclzm1(1)
         s%job%rsclzm2(i) = s%job%rsclzm2(1)
         do j = 1,4
            s%job%rescale_params(j,i) = s%job%rescale_params(j,1)
         end do
      end do
! Card 2 stops at 1e8 yr and card 3 at target_solar_age (legacy
! CALSOLAGE; MHP 06/13), both scaled by age_scale_factor, as are their
! timesteps.
      s%job%target_end_age(2) = 1.0d8*age_scale_factor
      s%job%timestep_override(2) = s%job%timestep_override(2)*age_scale_factor
      s%job%target_end_age(3) = s%ctrl%target_solar_age*age_scale_factor
      s%job%timestep_override(3) = s%job%timestep_override(3)*age_scale_factor
! The settling cards repeat card 2.
      do i = 5,47,3
         s%job%rescale_kind(i) = 1
         s%job%first_call_flag(i) = .false.
         s%job%num_models(i) = s%job%num_models(2)
         s%job%target_end_age(i) = s%job%target_end_age(2)
         s%job%end_age_stop_active(i) = s%job%end_age_stop_active(2)
         s%job%timestep_override(i) = s%job%timestep_override(2)
         s%job%timestep_override_active(i) = s%job%timestep_override_active(2)
      end do
! The cards that evolve to the solar age repeat card 3.
      do i = 6,num_calibration_runs,3
         s%job%rescale_kind(i) = 1
         s%job%first_call_flag(i) = .false.
         s%job%num_models(i) = s%job%num_models(3)
         s%job%target_end_age(i) = s%job%target_end_age(3)
         s%job%end_age_stop_active(i) = s%job%end_age_stop_active(3)
         s%job%timestep_override(i) = s%job%timestep_override(3)
         s%job%timestep_override_active(i) = s%job%timestep_override_active(3)
      end do
      s%solar_calibration_converged = .false.
      return
end subroutine setup_solar_calibration

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! setup_star_calibration (legacy setscal.f)
!----------------------------------------------------------------------
! Sets up the run list to evolve a model to a target point in the HR
! diagram: the target luminosity L* (target_luminosity_lsun) and either
! the target Teff or the target radius R* (specify_teff_flag); the
! other follows from L = 4 pi R^2 sigma Teff^4. num_runs becomes
! max_runs = 50 cards (25 pairs). Every card takes card 1's X, Z, alpha
! and envelope-mass request; the odd cards 3, 5, ..., 49 repeat card
! 1's rescaling, and the even cards 4, 6, ..., 50 continue the model
! with card 2's model count, stop age and timestep, as in
! setup_solar_calibration.
subroutine setup_star_calibration(s)

      use controls_buffer, only: max_runs
      use io_units
      use physical_constants
      use math_lib
      implicit none
      type(star_info), intent(inout) :: s
! --- locals ---
      integer :: i, j

      s%star_found_flag = .false.
      s%just_passed_target_radius_flag = .false.
! R*/Rsun = sqrt(L* Lsun/(4 pi sigma))/(Teff*^2 Rsun), or
! Teff* = [L* Lsun/(4 pi sigma R*^2 Rsun^2)]^(1/4).
      if (s%ctrl%specify_teff_flag) then
         s%job%target_radius_rsun = sqrt(s%ctrl%target_luminosity_lsun* &
              s%solar_luminosity_cgs/(c4pi*csig))/(s%job%target_teff*s%job%target_teff* &
              s%solar_radius_cgs)
      else
         s%job%target_teff = pow((s%ctrl%target_luminosity_lsun*s%solar_luminosity_cgs)/ &
              (c4pi*csig*s%job%target_radius_rsun*s%job%target_radius_rsun* &
              s%solar_radius_cgs*s%solar_radius_cgs), 0.25d0)
      end if
! No previous model yet.
      s%r_rsun_prev_model = 0
! Use all max_runs = 50 cards, copying the seed cards' parameters to
! the 48 calibration cards that follow.
      s%job%num_runs = max_runs
      do i = 2,max_runs
         s%job%initial_x_array(i) = s%job%initial_x_array(1)
         s%job%initial_z_array(i) = s%job%initial_z_array(1)
         s%job%mixing_length_array(i) = s%job%mixing_length_array(1)
         s%job%has_senv0_array(i) = s%job%has_senv0_array(1)
         s%job%senv0_array(i) = s%job%senv0_array(1)
      end do
! The rescale cards repeat card 1.
      do i = 3,max_runs-1,2
         s%job%rescale_kind(i) = s%job%rescale_kind(1)
         s%job%first_call_flag(i) = .true.
         s%job%num_models(i) = s%job%num_models(1)
         s%job%rsclzc(i) = s%job%rsclzc(1)
         s%job%rsclzm1(i) = s%job%rsclzm1(1)
         s%job%rsclzm2(i) = s%job%rsclzm2(1)
         do j = 1,4
            s%job%rescale_params(j,i) = s%job%rescale_params(j,1)
         end do
      end do
! The evolving cards repeat card 2.
      do i = 4,max_runs,2
         s%job%rescale_kind(i) = 1
         s%job%first_call_flag(i) = .false.
         s%job%num_models(i) = s%job%num_models(2)
         s%job%target_end_age(i) = s%job%target_end_age(2)
         s%job%end_age_stop_active(i) = s%job%end_age_stop_active(2)
         s%job%timestep_override(i) = s%job%timestep_override(2)
         s%job%timestep_override_active(i) = s%job%timestep_override_active(2)
      end do
      write(*,*) ' Evolve to R*, L* = ', s%job%target_radius_rsun, &
           s%ctrl%target_luminosity_lsun
      write(run_log_unit,*) '#Evolve to R*, L* = ', s%job%target_radius_rsun, &
           s%ctrl%target_luminosity_lsun
      return
end subroutine setup_star_calibration

end module calibration
