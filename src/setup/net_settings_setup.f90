!----------------------------------------------------------------------
! net_settings_setup
!----------------------------------------------------------------------
! The network's settings record, built from the run's controls and rate
! scales by the star layer; run_yrec installs it (net_lib's
! set_net_settings, which checks it and the reaction table) at the start
! of every Monte-Carlo run, after the run's rate factors are applied.
module net_settings_setup
      use net_lib, only: net_settings_t
      implicit none
      private
      public :: net_settings_from_star

contains

!----------------------------------------------------------------------
! net_settings_from_star: the network's settings from the controls and
! the star's cross-section scales (io/user_inputs.f90 builds those
! from the s0_* inputs; the Monte-Carlo runs rescale four of them).
!   star reads: be7p_s_double_prime_scale, be7p_s_prime_scale,
!     cross_section_scale, ctrl%be9_palpha_rate_scale,
!     ctrl%be9_pd_rate_scale, ctrl%be9_pg_rate_scale,
!     ctrl%calc_neutrinos, ctrl%li6_rate_scale, ctrl%li7_rate_scale,
!     ctrl%nuclear_logT_cutoffs, ctrl%use_itoh_neutrino_loss,
!     ctrl%weak_screening_threshold, extra_cross_section_scale,
!     n15_pa_scale, n15_pg_scale, s_double_prime_scale, s_prime_scale
!   star writes: none
pure function net_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(net_settings_t) :: settings

      settings%logT_cutoffs = s%ctrl%nuclear_logT_cutoffs
      settings%weak_screening_threshold = s%ctrl%weak_screening_threshold
      settings%cross_section_scale = s%cross_section_scale
      settings%s_prime_scale = s%s_prime_scale
      settings%s_double_prime_scale = s%s_double_prime_scale
      settings%extra_cross_section_scale = s%extra_cross_section_scale
      settings%be7p_s_prime_scale = s%be7p_s_prime_scale
      settings%be7p_s_double_prime_scale = s%be7p_s_double_prime_scale
      settings%n15_pg_scale = s%n15_pg_scale
      settings%n15_pa_scale = s%n15_pa_scale
      settings%li6_rate_scale = s%ctrl%li6_rate_scale
      settings%li7_rate_scale = s%ctrl%li7_rate_scale
      settings%be9_pg_rate_scale = s%ctrl%be9_pg_rate_scale
      settings%be9_pd_rate_scale = s%ctrl%be9_pd_rate_scale
      settings%be9_palpha_rate_scale = s%ctrl%be9_palpha_rate_scale
      settings%calc_neutrinos = s%ctrl%calc_neutrinos
      settings%use_itoh_neutrino_loss = s%ctrl%use_itoh_neutrino_loss
end function net_settings_from_star

end module net_settings_setup
