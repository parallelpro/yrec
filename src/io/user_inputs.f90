!----------------------------------------------------------------------
! user_inputs
!----------------------------------------------------------------------
! map_user_inputs (legacy remap.f, MHP 6/2014), a step of the controls
! read (after adopt_canonical_names): the descriptively named controls
! onto the arrays the code uses.
!  - With use_new_nuclear_rates (legacy LNEWNUC), the per-reaction S0
!    values and their first and second derivatives (s0_*, s0p_*,
!    s0pp_*) become the cross-section scales s%cross_section_scale
!    (legacy SSTANDARD), the extra scales s%extra_cross_section_scale,
!    the derivative scales s%s_prime_scale and s%s_double_prime_scale (legacy QS0E,
!    QQS0EE; Be7 + p's are s%be7p_s_prime_scale and s%be7p_s_double_prime_scale)
!    and the two N15 + p branch scales.
!  - With use_new_rate_variables (use_new_rate_variables, which despite its name has
!    nothing to do with reaction rates), the tol_*, flag_* and time_*
!    controls overwrite the rezoning tolerances chi_grid_scale (legacy
!    HPTTOL) and the timestep limiters timestep_factors, and set
!    use_structure_dt_limits, max_domega_global and
!    diffusion_timestep_factor.
! The derivative and N15 + p branch scales are always first reset to
! unity.
!
! The network follows C12, C13, N14 and O16 out of equilibrium. N15
! burning is assumed to follow N14 + p at once, with the branching
! between N15(p,alpha)C12 and N15(p,gamma)O16 set by the relative S0 of
! the two; to change it, change s0_n15_p_c12_branch and
! s0_n15_p_o16_branch. O17 proton capture is assumed to follow O16
! proton capture at once, producing N14. The helium-burning reactions
! (triple alpha, C12 + He4, C13 + He4, N14 + He4) use the Caughlan &
! Fowler (1988) rates; resonances matter for these, so a simple scaling
! by S0 is not enough, and they keep unit scale (updates should use other
! functional forms, e.g. NACRE II).
module user_inputs
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use star_info_lib, only: star_info
      use controls_buffer, only: use_new_nuclear_rates, use_new_rate_variables, &
           s0_pp, s0_he3he3, s0_he3he4, s0_p_c12, s0_p_c13, s0_p_n14, s0_p_o16, &
           s0_pep, s0_be7_electron, s0_be7_p, s0_hep, s0_n15_p_c12_branch, &
           s0_n15_p_o16_branch, s0p_pp, s0p_he3he3, s0p_he3he4, s0p_p_c12, &
           s0p_p_c13, s0p_p_n14, s0p_p_o16, s0p_be7_p, s0pp_p_c12, s0pp_p_c13, &
           s0pp_p_o16, s0pp_be7_p, &
           chi_grid_scale, ichi_dm_min, ichi_dm_max, ichi_dl_max, ichi_dp_core_max, &
           ichi_dp_env_max, ichi_dp_czbase_max, ichi_czbase_fine_width, ichi_dx_max, &
           ichi_dz_max, ichi_flag_dx, ichi_flag_dz, ichi_flag_dw, &
           tol_dm_min, tol_dm_max, tol_dl_max, tol_dp_core_max, tol_dp_env_max, &
           tol_dp_czbase_max, tol_czbase_fine_width, tol_dx_max, tol_dz_max, &
           flag_dx, flag_dz, flag_dw, &
           timestep_factors, itime_core_min, itime_dx_core_tot, itime_dx_core_frac, &
           itime_dx_total, itime_dx_shell, itime_dy_core_tot, itime_dy_core_frac, &
           itime_dy_shell, itime_dy_total, itime_dt, itime_dp, itime_dr, itime_dl, &
           itime_max_dt_frac, &
           time_core_min, time_dx_core_tot, time_dx_core_frac, time_dx_total, &
           time_dx_shell, time_dy_core_tot, time_dy_core_frac, time_dy_shell, &
           time_dy_total, time_dt, time_dp, time_dr, time_dl, time_max_dt_frac, &
           time_dw_global, time_dw_mix, use_structure_time_limits, &
           use_structure_dt_limits, max_domega_global, diffusion_timestep_factor
      use net_lib, only: r_pp, r_he3he3, r_he3he4, r_pc12, r_pc13, r_pn14, r_po16, &
           xs_pep, xs_be7e, xs_be7p, xs_hep, reference_s0, pep_coefficient, &
           be7_ecapture_coefficient, s0_be7p_kevb, s0_hep_kevb, s0_n15pa_kevb, &
           s0_n15pg_kevb
      implicit none
      private
      public :: map_user_inputs

! The reference S0 values (keV b) against which the input S0 values are
! scaled are net's (net_lib): each table reaction's reference_s0, the
! S0 of Bahcall & Ulrich (1988) its fit is normalized to; for pep, the
! proportionality constant of the pep rate relative to pp
! (pep_coefficient), for Be7 + e- the constant of its rate
! (be7_ecapture_coefficient) -- both rates go as T^(-1/2), and the code
! uses T9 rather than T6, so a constant quoted for T6 = T/1e6 K is
! divided by 1000^(1/2); the S0 of Be7 + p and hep's fits (s0_be7p_kevb,
! s0_hep_kevb) and of the N15 + p branches (s0_n15pa_kevb,
! s0_n15pg_kevb). See set_cross_section_scales.
! The Solar Fusion I S0 values (Adelberger et al. 1998), in keV b. The
! rate fits' first- and second-derivative terms are hard-coded as S'/S
! and S''/S from that compilation, so each reaction's input S'/S and
! S''/S are scaled against its own Solar Fusion I ratio. (2026-10: the
! F77 remap.f divided reactions 2-7 by the pp S0 instead, which made
! those scales ~1e-22 and silently removed the q2..q5 terms of He3+He3,
! He3+He4, C12+p, C13+p, N14+p and O16+p.) No derivative scales are
! taken for pep, Be7 + e-, hep or the N15 + p branches.
      real(dp), parameter :: s0_pp_a98 = 4.00d-22, s0_he3he3_a98 = 5.4d3, &
           s0_he3he4_a98 = 5.3d-1, s0_p_c12_a98 = 1.34d0, s0_p_c13_a98 = 7.6d0, &
           s0_p_n14_a98 = 3.5d0, s0_p_o16_a98 = 9.4d0, s0_be7_p_a98 = 0.019d0
! The Solar Fusion I first derivatives S'(0) (Adelberger et al. 1998),
! in b.
      real(dp), parameter :: s_prime_pp_a98 = 4.48d-24, s_prime_he3he3_a98 = -4.1d0, &
           s_prime_he3he4_a98 = -3.0d-4, s_prime_p_c12_a98 = 2.6d-3, &
           s_prime_p_c13_a98 = -7.8d-3, s_prime_p_n14_a98 = -1.28d-2, &
           s_prime_p_o16_a98 = -2.4d-2, s_prime_be7_p_a98 = -1.33d-5
! The Solar Fusion I second derivatives S''(0) (Adelberger et al. 1998),
! in b/keV.
      real(dp), parameter :: s_double_prime_p_c12_a98 = 8.3d-5, s_double_prime_p_c13_a98 = 7.3d-4, &
           s_double_prime_p_o16_a98 = 5.7d-5, s_double_prime_be7_p_a98 = 7.22d-8

contains

subroutine map_user_inputs(s)
      type(star_info), intent(inout) :: s

      call reset_rate_scale_defaults(s)
      if (use_new_nuclear_rates) then
         call set_cross_section_scales(s)
         call set_first_derivative_scales(s)
         call set_second_derivative_scales(s)
      end if
      if (use_new_rate_variables) then
         call map_zoning_tolerances
         call map_timestep_tolerances
      end if
end subroutine map_user_inputs

! ---------------------------------------------------------------
! The rate-derivative scales and the N15 + p branch scales, at unity.
subroutine reset_rate_scale_defaults(s)
      type(star_info), intent(inout) :: s
      s%s_prime_scale = 1.0d0
      s%s_double_prime_scale = 1.0d0
      s%be7p_s_prime_scale = 1.0d0
      s%be7p_s_double_prime_scale = 1.0d0
      s%n15_pa_scale = 1.0d0
      s%n15_pg_scale = 1.0d0
end subroutine reset_rate_scale_defaults

! ---------------------------------------------------------------
! Each reaction's cross-section scale, S0 over net's reference S0, and
! the scales of the two N15 + p branches. The S-factor law
! (net/reaction_rates.f90's compute_s_factor_rate) multiplies each rate
! by cross_section_scale (legacy SSTANDARD), the ratio of the desired S0
! to that of Bahcall & Ulrich (1988); here the scales come from explicit
! S0 values for the individual reactions. The helium-burning reactions
! keep unit scale.
subroutine set_cross_section_scales(s)
      type(star_info), intent(inout) :: s
      s%cross_section_scale = 1.0d0
      s%extra_cross_section_scale = 1.0d0
      s%cross_section_scale(r_pp) = s0_pp/reference_s0(r_pp)
      s%cross_section_scale(r_he3he3) = s0_he3he3/reference_s0(r_he3he3)
      s%cross_section_scale(r_he3he4) = s0_he3he4/reference_s0(r_he3he4)
      s%cross_section_scale(r_pc12) = s0_p_c12/reference_s0(r_pc12)
      s%cross_section_scale(r_pc13) = s0_p_c13/reference_s0(r_pc13)
      s%cross_section_scale(r_pn14) = s0_p_n14/reference_s0(r_pn14)
      s%cross_section_scale(r_po16) = s0_p_o16/reference_s0(r_po16)
      s%extra_cross_section_scale(xs_pep) = s0_pep/pep_coefficient
      s%extra_cross_section_scale(xs_be7e) = s0_be7_electron/be7_ecapture_coefficient
      s%extra_cross_section_scale(xs_be7p) = s0_be7_p/s0_be7p_kevb
      s%extra_cross_section_scale(xs_hep) = s0_hep/s0_hep_kevb
! N15(p,alpha)C12 and N15(p,gamma)O16
      s%n15_pa_scale = s0_n15_p_c12_branch/s0_n15pa_kevb
      s%n15_pg_scale = s0_n15_p_o16_branch/s0_n15pg_kevb
end subroutine set_cross_section_scales

! ---------------------------------------------------------------
! The first-derivative scales. The rate fits (net/reaction_rates.f90;
! legacy ENGEB and NRATES) carry the S-factor derivatives in their q2-q5
! terms: q2 and q3 are proportional to S'/S and q4 and q5 to S''/S, both
! from Solar Fusion I; s_prime_scale multiplies the q2 and q3 terms. Be7 + p,
! outside the reaction table, has its own.
subroutine set_first_derivative_scales(s)
      type(star_info), intent(inout) :: s
      s%s_prime_scale(r_pp) = derivative_scale(s0p_pp, s0_pp, s_prime_pp_a98, s0_pp_a98)
      s%s_prime_scale(r_he3he3) = derivative_scale(s0p_he3he3, s0_he3he3, &
           s_prime_he3he3_a98, s0_he3he3_a98)
      s%s_prime_scale(r_he3he4) = derivative_scale(s0p_he3he4, s0_he3he4, &
           s_prime_he3he4_a98, s0_he3he4_a98)
      s%s_prime_scale(r_pc12) = derivative_scale(s0p_p_c12, s0_p_c12, s_prime_p_c12_a98, s0_p_c12_a98)
      s%s_prime_scale(r_pc13) = derivative_scale(s0p_p_c13, s0_p_c13, s_prime_p_c13_a98, s0_p_c13_a98)
      s%s_prime_scale(r_pn14) = derivative_scale(s0p_p_n14, s0_p_n14, s_prime_p_n14_a98, s0_p_n14_a98)
      s%s_prime_scale(r_po16) = derivative_scale(s0p_p_o16, s0_p_o16, s_prime_p_o16_a98, s0_p_o16_a98)
      s%be7p_s_prime_scale = derivative_scale(s0p_be7_p, s0_be7_p, s_prime_be7_p_a98, s0_be7_p_a98)
end subroutine set_first_derivative_scales

! ---------------------------------------------------------------
! The second-derivative scales: s_double_prime_scale multiplies the q4 and q5
! terms of the rate fits. None are set for pp, He3 + He3, He3 + He4 or
! N14 + p (dropped in the 2003 version); theirs stay at unity.
subroutine set_second_derivative_scales(s)
      type(star_info), intent(inout) :: s
      s%s_double_prime_scale(r_pc12) = derivative_scale(s0pp_p_c12, s0_p_c12, &
           s_double_prime_p_c12_a98, s0_p_c12_a98)
      s%s_double_prime_scale(r_pc13) = derivative_scale(s0pp_p_c13, s0_p_c13, &
           s_double_prime_p_c13_a98, s0_p_c13_a98)
      s%s_double_prime_scale(r_po16) = derivative_scale(s0pp_p_o16, s0_p_o16, &
           s_double_prime_p_o16_a98, s0_p_o16_a98)
      s%be7p_s_double_prime_scale = derivative_scale(s0pp_be7_p, s0_be7_p, &
           s_double_prime_be7_p_a98, s0_be7_p_a98)
end subroutine set_second_derivative_scales

! ---------------------------------------------------------------
! A derivative's scale, (S^(n)/S0)_input / (S^(n)/S0)_SolarFusionI: zero
! where the input S0 is not positive (MHP 4/25), instead of a division
! by zero.
pure function derivative_scale(derivative, s0, derivative_a98, s0_a98) result(scale)
      real(dp), intent(in) :: derivative, s0, derivative_a98, s0_a98
      real(dp) :: scale
      if (s0 > 0.0d0) then
         scale = (derivative/s0)/(derivative_a98/s0_a98)
      else
         scale = 0.0d0
      end if
end function derivative_scale

! ---------------------------------------------------------------
! The named zoning tolerances into chi_grid_scale: the log m, L/L_surface
! and log P spacings (radiative interior, convective envelope, near the
! base of the convection zone), the width in log P of the finely zoned
! region there, the X and Z spacings (usually disabled), and the X, Z and
! log omega (rotating models) differences that flag a shell (to avoid
! numerical diffusion).
subroutine map_zoning_tolerances
      chi_grid_scale(ichi_dm_min) = tol_dm_min
      chi_grid_scale(ichi_dm_max) = tol_dm_max
      chi_grid_scale(ichi_dl_max) = tol_dl_max
      chi_grid_scale(ichi_dp_core_max) = tol_dp_core_max
      chi_grid_scale(ichi_dp_env_max) = tol_dp_env_max
      chi_grid_scale(ichi_dp_czbase_max) = tol_dp_czbase_max
      chi_grid_scale(ichi_czbase_fine_width) = tol_czbase_fine_width
      chi_grid_scale(ichi_dx_max) = tol_dx_max
      chi_grid_scale(ichi_dz_max) = tol_dz_max
      chi_grid_scale(ichi_flag_dx) = flag_dx
      chi_grid_scale(ichi_flag_dz) = flag_dz
      chi_grid_scale(ichi_flag_dw) = flag_dw
end subroutine map_zoning_tolerances

! ---------------------------------------------------------------
! The named timestep tolerances into timestep_factors (the itime_* slots: the core-
! burning criterion's minimum central abundance; the X and Y burned in
! the core, absolutely and as a fraction, in the shell, and in total
! (Msun); the largest changes of T, P, R and L at any shell; the largest
! timestep change from model to model) and the other timestep controls.
! use_structure_dt_limits switches the structure limits (the T, P, R and
! L slots) on or off; they are meant for the pre-main sequence, and the
! code switches them off by itself later on. max_domega_global is the
! largest change of omega per global timestep, diffusion_timestep_factor
! the largest in one of rotational mixing's sub-steps.
subroutine map_timestep_tolerances
      timestep_factors(itime_core_min) = time_core_min
      timestep_factors(itime_dx_core_tot) = time_dx_core_tot
      timestep_factors(itime_dx_core_frac) = time_dx_core_frac
      timestep_factors(itime_dx_total) = time_dx_total
      timestep_factors(itime_dx_shell) = time_dx_shell
      timestep_factors(itime_dy_core_tot) = time_dy_core_tot
      timestep_factors(itime_dy_core_frac) = time_dy_core_frac
      timestep_factors(itime_dy_shell) = time_dy_shell
      timestep_factors(itime_dy_total) = time_dy_total
      timestep_factors(itime_dt) = time_dt
      timestep_factors(itime_dp) = time_dp
      timestep_factors(itime_dr) = time_dr
      timestep_factors(itime_dl) = time_dl
      use_structure_dt_limits = use_structure_time_limits
      max_domega_global = time_dw_global
      diffusion_timestep_factor = time_dw_mix
      timestep_factors(itime_max_dt_frac) = time_max_dt_frac
end subroutine map_timestep_tolerances

end module user_inputs
