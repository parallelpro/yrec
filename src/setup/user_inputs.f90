!----------------------------------------------------------------------
! user_inputs
!----------------------------------------------------------------------
! The mapping of the newer, descriptively named controls onto the
! parameters the code uses (map_user_inputs), during the controls
! read.
module user_inputs
      implicit none
      private
      public :: map_user_inputs

contains

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! map_user_inputs (legacy remap.f, MHP 6/2014)
!----------------------------------------------------------------------
! Maps the newer, descriptively named namelist controls onto the
! arrays the code uses:
!  - with use_new_nuclear_rates (legacy LNEWNUC), the per-reaction S0
!    values and their first and second derivatives (s0_*, s0p_*,
!    s0pp_*) become the cross-section scales star%cross_section_scale
!    (legacy SSTANDARD) and the derivative scales star%qs0e_scale and
!    star%qqs0ee_scale (legacy QS0E, QQS0EE);
!  - with lnewvars (use_new_rate_variables), the tol_*, flag_* and
!    time_* controls overwrite the rezoning tolerances chi_grid_scale
!    (legacy HPTTOL) and the timestep limiters atime, and set
!    use_structure_dt_limits, max_domega_global and
!    diffusion_timestep_factor.
! The derivative and N15 + p branching scales are always first reset
! to unity.
subroutine map_user_inputs
      use star_info_lib, only: star
      use controls_buffer
      use net_lib, only: r_pp, r_he3he3, r_he3he4, r_pc12, r_pc13, r_pn14, r_po16, &
           xs_pep, xs_be7e, xs_be7p, xs_hep, iq_be7p
      implicit none

      integer :: i

! The reference S0 values (keV b) against which the input S0 values are
! scaled: those of Bahcall & Ulrich (1988), the normalization of the
! rates in compute_zone_energy (see set_cross_section_scales).
      double precision :: s0_pp_ref, s0_he3he3_ref, s0_he3he4_ref, &
           s0_p_c12_ref, s0_p_c13_ref
      data s0_pp_ref,s0_he3he3_ref,s0_he3he4_ref,s0_p_c12_ref,s0_p_c13_ref/ &
           4.07D-22,5.15D3,5.4D-1,1.45D0,5.5D0/
! For pep, the reference is the proportionality constant of the pep
! rate relative to pp; for Be7 + e-, the constant of the linear term.
! Both rates go as T^(-1/2), and the code uses T9 rather than T6, so a
! constant quoted for T6 = T/1e6 K must be divided by 1000^(1/2).
      double precision :: s0_p_n14_ref, s0_p_o16_ref, s0_pep_ref, &
           s0_be7_electron_ref
      data s0_p_n14_ref,s0_p_o16_ref,s0_pep_ref,s0_be7_electron_ref/ &
           3.32D0,9.4E0,3.4848D-6,1.752D-10/
      double precision :: s0_be7_p_ref, s0_hep_ref, s0_n15_p_c12_branch_ref, &
           s0_n15_p_o16_branch_ref
      data s0_be7_p_ref,s0_hep_ref,s0_n15_p_c12_branch_ref, &
           s0_n15_p_o16_branch_ref/0.0243D0,8.0D-20,6.75D4,6.4D1/
! The Solar Fusion I S0 values (Adelberger et al. 1998), in keV b. The
! rate fits' first- and second-derivative terms are hard-coded as S'/S
! and S''/S from that compilation, so each reaction's input S'/S and
! S''/S are scaled below against its own Solar Fusion I ratio.
! (2026-10: the F77 remap.f divided reactions 2-7 by the pp S0
! instead, which made those scales ~1e-22 and silently removed the
! q2..q5 terms of He3+He3, He3+He4, C12+p, C13+p, N14+p and O16+p.)
      double precision :: s0_pp_a98, s0_he3he3_a98, s0_he3he4_a98, &
           s0_p_c12_a98, s0_p_c13_a98, s0_p_n14_a98, s0_p_o16_a98
      data s0_pp_a98,s0_he3he3_a98,s0_he3he4_a98,s0_p_c12_a98,s0_p_c13_a98/ &
           4.00D-22,5.4D3,5.3D-1,1.34D0,7.6D0/
      data s0_p_n14_a98,s0_p_o16_a98/3.5D0,9.4D0/
! No derivative scales are taken for pep, Be7 + e-, hep or the N15 + p
! branches, so their Solar Fusion I values are not needed.
      double precision :: s0_be7_p_a98
      data s0_be7_p_a98/0.019D0/
! The Solar Fusion I first derivatives S'(0) (Adelberger et al. 1998),
! in b.
      double precision :: qs0e_pp_a98, qs0e_he3he3_a98, qs0e_he3he4_a98, &
           qs0e_p_c12_a98
      data qs0e_pp_a98,qs0e_he3he3_a98,qs0e_he3he4_a98,qs0e_p_c12_a98/ &
           4.48D-24,-4.1D0,-3.0D-4,2.6D-3/
      double precision :: qs0e_p_c13_a98, qs0e_p_n14_a98, qs0e_p_o16_a98, &
           qs0e_be7_p_a98
      data qs0e_p_c13_a98,qs0e_p_n14_a98,qs0e_p_o16_a98,qs0e_be7_p_a98/ &
           -7.8D-3,-1.28D-2,-2.4D-2,-1.33D-5/
! The Solar Fusion I second derivatives S''(0) (Adelberger et al.
! 1998), in b/keV.
      double precision :: qqs0ee_p_c12_a98, qqs0ee_p_c13_a98, &
           qqs0ee_p_o16_a98, qqs0ee_be7_p_a98
      data qqs0ee_p_c12_a98,qqs0ee_p_c13_a98,qqs0ee_p_o16_a98, &
           qqs0ee_be7_p_a98/8.3D-5,7.3D-4,5.7D-5,7.22D-8/
! The network follows C12, C13, N14 and O16 out of equilibrium. N15
! burning is assumed to follow N14 + p at once, with the branching
! between N15(p,alpha)C12 and N15(p,gamma)O16 set by the relative S0
! of the two; to change it, change s0_n15_p_c12_branch and
! s0_n15_p_o16_branch. O17 proton capture is assumed to follow O16
! proton capture at once, producing N14. The helium-burning reactions
! (triple alpha, C12 + He4, C13 + He4, N14 + He4) use the Caughlan &
! Fowler (1988) rates. Resonances matter for these, so a simple
! scaling by S0 is not enough; updates should use other functional
! forms instead (e.g. NACRE II).
!
! With use_new_nuclear_rates (legacy LNEWNUC), the cross-section
! scales are reset to unity and rebuilt from the individual S0 values
! and their derivatives.
      call reset_rate_scale_defaults
      if (use_new_nuclear_rates) then
         call set_cross_section_scales
         call set_first_derivative_scales
         call set_second_derivative_scales
      end if
! With lnewvars (use_new_rate_variables, which despite its name has
! nothing to do with reaction rates), the named zoning and timestep
! tolerances overwrite chi_grid_scale (legacy HPTTOL) and atime.
      if (lnewvars) then
         call map_zoning_tolerances
         call map_timestep_tolerances
      end if
      return

contains

! ---------------------------------------------------------------
! Reset the rate-derivative scales and the N15 + p branching scales to
! unity.
subroutine reset_rate_scale_defaults
      do i = 1,8
         star%qs0e_scale(i) = 1.0d0
         star%qqs0ee_scale(i) = 1.0d0
      end do
      star%c12_alpha_scale = 1.0d0
      star%o16_gamma_scale = 1.0d0
end subroutine reset_rate_scale_defaults

! ---------------------------------------------------------------
! Set each reaction's cross-section scale, S0/S0_ref, and the scales
! of the two N15 + p branches.
subroutine set_cross_section_scales
         do i = 1,17
            star%cross_section_scale(i) = 1.0d0
         end do
! compute_zone_energy multiplies each rate by cross_section_scale
! (legacy SSTANDARD), the ratio of the desired S0 to that of Bahcall &
! Ulrich (1988). Here the namelist vector is replaced by explicit S0
! values for the individual reactions, and the implied scales are
! derived. The first and second derivatives may be changed separately
! (unlike in the original compute_zone_energy): their inputs are S'/S0
! and S''/S0, taken relative to the Solar Fusion I (Adelberger et al.
! 1998) ratios.
         star%cross_section_scale(r_pp) = s0_pp/s0_pp_ref
         star%cross_section_scale(r_he3he3) = s0_he3he3/s0_he3he3_ref
         star%cross_section_scale(r_he3he4) = s0_he3he4/s0_he3he4_ref
         star%cross_section_scale(r_pc12) = s0_p_c12/s0_p_c12_ref
         star%cross_section_scale(r_pc13) = s0_p_c13/s0_p_c13_ref
         star%cross_section_scale(r_pn14) = s0_p_n14/s0_p_n14_ref
         star%cross_section_scale(r_po16) = s0_p_o16/s0_p_o16_ref
! Slots 8-13 keep unit scale: C13 + He4, C12 + He4, N14 + He4 and
! triple alpha are resonant and not fit by a simple S0, and O16 + He4
! and C12 + C12 are not used.
         star%cross_section_scale(xs_pep) = s0_pep/s0_pep_ref
         star%cross_section_scale(xs_be7e) = s0_be7_electron/s0_be7_electron_ref
         star%cross_section_scale(xs_be7p) = s0_be7_p/s0_be7_p_ref
         star%cross_section_scale(xs_hep) = s0_hep/s0_hep_ref
! The scales of the N15 + p branches, N15(p,alpha)C12 and
! N15(p,gamma)O16.
         star%c12_alpha_scale = s0_n15_p_c12_branch/s0_n15_p_c12_branch_ref
         star%o16_gamma_scale = s0_n15_p_o16_branch/s0_n15_p_o16_branch_ref
end subroutine set_cross_section_scales

! ---------------------------------------------------------------
! Set the first-derivative scales, for each reaction
!   qs0e_scale = (S'/S0)_input / (S'/S0)_SolarFusionI,
! zero where the input S0 is not positive.
subroutine set_first_derivative_scales
! The rate fits (net/reaction_rates.f90; legacy ENGEB and NRATES) carry
! the S-factor derivatives in their q2-q5 terms: q2 and q3 are
! proportional to S'/S and q4 and q5 to S''/S, both from Solar Fusion I
! (Adelberger et al. 1998). qs0e_scale multiplies the q2 and q3 terms.
! A non-positive S0 gives a zero scale (MHP 4/25) instead of a division
! by zero.
         if (s0_pp.gt.0.0d0) then
            star%qs0e_scale(r_pp) = (s0p_pp/s0_pp)/(qs0e_pp_a98/s0_pp_a98)
         else
            star%qs0e_scale(r_pp) = 0.0d0
         end if
         if (s0_he3he3.gt.0.0d0) then
            star%qs0e_scale(r_he3he3) = (s0p_he3he3/s0_he3he3)/(qs0e_he3he3_a98/s0_he3he3_a98)
         else
            star%qs0e_scale(r_he3he3) = 0.0d0
         end if
         if (s0_he3he4.gt.0.0d0) then
            star%qs0e_scale(r_he3he4) = (s0p_he3he4/s0_he3he4)/(qs0e_he3he4_a98/s0_he3he4_a98)
         else
            star%qs0e_scale(r_he3he4) = 0.0d0
         end if
         if (s0_p_c12.gt.0.0d0) then
            star%qs0e_scale(r_pc12) = (s0p_p_c12/s0_p_c12)/(qs0e_p_c12_a98/s0_p_c12_a98)
         else
            star%qs0e_scale(r_pc12) = 0.0d0
         end if
         if (s0_p_c13.gt.0.0d0) then
            star%qs0e_scale(r_pc13) = (s0p_p_c13/s0_p_c13)/(qs0e_p_c13_a98/s0_p_c13_a98)
         else
            star%qs0e_scale(r_pc13) = 0.0d0
         end if
         if (s0_p_n14.gt.0.0d0) then
            star%qs0e_scale(r_pn14) = (s0p_p_n14/s0_p_n14)/(qs0e_p_n14_a98/s0_p_n14_a98)
         else
            star%qs0e_scale(r_pn14) = 0.0d0
         end if
         if (s0_p_o16.gt.0.0d0) then
            star%qs0e_scale(r_po16) = (s0p_p_o16/s0_p_o16)/(qs0e_p_o16_a98/s0_p_o16_a98)
         else
            star%qs0e_scale(r_po16) = 0.0d0
         end if
! Slot iq_be7p (8) of the derivative scales holds Be7 + p.
         if (s0_be7_p.gt.0.0d0) then
            star%qs0e_scale(iq_be7p) = (s0p_be7_p/s0_be7_p)/ &
                 (qs0e_be7_p_a98/s0_be7_p_a98)
         else
            star%qs0e_scale(iq_be7p) = 0.0d0
         end if
end subroutine set_first_derivative_scales

! ---------------------------------------------------------------
! Set the second-derivative scales, for each reaction
!   qqs0ee_scale = (S''/S0)_input / (S''/S0)_SolarFusionI,
! zero where the input S0 is not positive.
subroutine set_second_derivative_scales
! qqs0ee_scale multiplies the q4 and q5 terms of the rate fits. No
! second-derivative scales are set for pp, He3 + He3, He3 + He4 or
! N14 + p (dropped in the 2003 version); theirs stay at unity.
         if (s0_p_c12.gt.0.0d0) then
            star%qqs0ee_scale(r_pc12) = (s0pp_p_c12/s0_p_c12)/ &
                 (qqs0ee_p_c12_a98/s0_p_c12_a98)
         else
            star%qqs0ee_scale(r_pc12) = 0.0d0
         end if
         if (s0_p_c13.gt.0.0d0) then
            star%qqs0ee_scale(r_pc13) = (s0pp_p_c13/s0_p_c13)/ &
                 (qqs0ee_p_c13_a98/s0_p_c13_a98)
         else
            star%qqs0ee_scale(r_pc13) = 0.0d0
         end if
         if (s0_p_o16.gt.0.0d0) then
            star%qqs0ee_scale(r_po16) = (s0pp_p_o16/s0_p_o16)/ &
                 (qqs0ee_p_o16_a98/s0_p_o16_a98)
         else
            star%qqs0ee_scale(r_po16) = 0.0d0
         end if
         if (s0_be7_p.gt.0.0d0) then
            star%qqs0ee_scale(iq_be7p) = (s0pp_be7_p/s0_be7_p)/ &
                 (qqs0ee_be7_p_a98/s0_be7_p_a98)
         else
            star%qqs0ee_scale(iq_be7p) = 0.0d0
         end if
end subroutine set_second_derivative_scales

! ---------------------------------------------------------------
! Copy the named zoning tolerances into chi_grid_scale.
subroutine map_zoning_tolerances
! Minimum log m spacing between shells.
         chi_grid_scale(ichi_dm_min) = tol_dm_min
! Maximum log m spacing between shells.
         chi_grid_scale(ichi_dm_max) = tol_dm_max
! Maximum L/L_surface spacing between shells.
         chi_grid_scale(ichi_dl_max) = tol_dl_max
! Maximum log P spacing between shells in the radiative interior.
         chi_grid_scale(ichi_dp_core_max) = tol_dp_core_max
! Maximum log P spacing between shells in the convective envelope.
         chi_grid_scale(ichi_dp_env_max) = tol_dp_env_max
! Maximum log P spacing between shells near the base of the convection
! zone.
         chi_grid_scale(ichi_dp_czbase_max) = tol_dp_czbase_max
! Width (in log P) of the finely zoned region near the base of the
! convection zone.
         chi_grid_scale(ichi_czbase_fine_width) = tol_czbase_fine_width
! Maximum X spacing between shells (usually disabled).
         chi_grid_scale(ichi_dx_max) = tol_dx_max
! Maximum Z spacing between shells (usually disabled).
         chi_grid_scale(ichi_dz_max) = tol_dz_max
! The X difference between shells that flags a shell (to avoid
! numerical diffusion).
         chi_grid_scale(ichi_flag_dx) = flag_dx
! The Z difference between shells that flags a shell (to avoid
! numerical diffusion).
         chi_grid_scale(ichi_flag_dz) = flag_dz
! The log omega difference between shells that flags a shell (rotating
! models only).
         chi_grid_scale(ichi_flag_dw) = flag_dw
end subroutine map_zoning_tolerances

! ---------------------------------------------------------------
! Copy the named timestep tolerances into atime and the other
! timestep controls.
subroutine map_timestep_tolerances
! The minimum central abundance for the core-burning timestep
! criterion.
         atime(itime_core_min) = time_core_min
! The largest absolute change in X burned in the core per step.
         atime(itime_dx_core_tot) = time_dx_core_tot
! The largest fraction of the core X burned per step.
         atime(itime_dx_core_frac) = time_dx_core_frac
! The largest mass of X (Msun) burned per step.
         atime(itime_dx_total) = time_dx_total
! The largest fraction of X burned in the shell per step.
         atime(itime_dx_shell) = time_dx_shell
! The largest absolute change in Y burned in the core per step.
         atime(itime_dy_core_tot) = time_dy_core_tot
! The largest fraction of the core Y burned per step.
         atime(itime_dy_core_frac) = time_dy_core_frac
! The largest fraction of Y burned in the shell per step.
         atime(itime_dy_shell) = time_dy_shell
! The largest mass of Y (Msun) burned in the shell per step.
         atime(itime_dy_total) = time_dy_total
! The largest changes in T, P, R and L at any shell per step.
         atime(itime_dt) = time_dt
         atime(itime_dp) = time_dp
         atime(itime_dr) = time_dr
         atime(itime_dl) = time_dl
! use_structure_dt_limits switches the structure limits above
! (atime(itime_dt) .. atime(itime_dl), slots 8-11) on or off. They are
! meant for the pre-main sequence; the code switches them off by itself
! later on.
         use_structure_dt_limits = lstruct_time
! The largest change in omega per global timestep.
         max_domega_global = time_dw_global
! Rotational mixing can run in a series of smaller steps within the
! timestep; this is the largest change in omega in one such step.
         diffusion_timestep_factor = time_dw_mix
! The largest factor by which the timestep may change from one model
! to the next.
         atime(itime_max_dt_frac) = time_max_dt_frac
end subroutine map_timestep_tolerances

end subroutine map_user_inputs

end module user_inputs
