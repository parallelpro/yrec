!----------------------------------------------------------------------
! net_settings_builder
!----------------------------------------------------------------------
! The network's settings record, built from the run's controls and rate
! scales by the star layer (setup/README.md); settings_install's
! install_monte_carlo_settings installs it (net_lib's set_net_settings,
! which checks it and the reaction table) at the start of every
! Monte-Carlo run, after begin_monte_carlo_run has applied the run's
! rate factors:
!     net_settings_from_star          net_settings_t
module net_settings_builder
      use net_lib, only: net_settings_t
      implicit none
      private
      public :: net_settings_from_star

contains

!----------------------------------------------------------------------
! net_settings_from_star: the network's settings from the controls and
! the star's cross-section scales (io/user_inputs.f90 builds those
! from the s0_* inputs; the Monte-Carlo runs rescale four of them).
!   star reads: cross_section_scale_in_use,
!     extra_cross_section_scale_in_use,
!     physics%be7p_s_double_prime_scale, physics%be7p_s_prime_scale,
!     physics%be9_palpha_rate_scale, physics%be9_pd_rate_scale,
!     physics%be9_pg_rate_scale, physics%calc_neutrinos,
!     physics%li6_rate_scale, physics%li7_rate_scale,
!     physics%n15_pa_scale, physics%n15_pg_scale,
!     physics%nuclear_logT_cutoffs, physics%s_double_prime_scale,
!     physics%s_prime_scale, physics%use_itoh_neutrino_loss,
!     physics%weak_screening_threshold
!   star writes: none
pure function net_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(net_settings_t) :: settings

      settings%logT_cutoffs = s%physics%nuclear_logT_cutoffs
      settings%weak_screening_threshold = s%physics%weak_screening_threshold
      settings%cross_section_scale = s%cross_section_scale_in_use
      settings%s_prime_scale = s%physics%s_prime_scale
      settings%s_double_prime_scale = s%physics%s_double_prime_scale
      settings%extra_cross_section_scale = s%extra_cross_section_scale_in_use
      settings%be7p_s_prime_scale = s%physics%be7p_s_prime_scale
      settings%be7p_s_double_prime_scale = s%physics%be7p_s_double_prime_scale
      settings%n15_pg_scale = s%physics%n15_pg_scale
      settings%n15_pa_scale = s%physics%n15_pa_scale
      settings%li6_rate_scale = s%physics%li6_rate_scale
      settings%li7_rate_scale = s%physics%li7_rate_scale
      settings%be9_pg_rate_scale = s%physics%be9_pg_rate_scale
      settings%be9_pd_rate_scale = s%physics%be9_pd_rate_scale
      settings%be9_palpha_rate_scale = s%physics%be9_palpha_rate_scale
      settings%calc_neutrinos = s%physics%calc_neutrinos
      settings%use_itoh_neutrino_loss = s%physics%use_itoh_neutrino_loss
end function net_settings_from_star

end module net_settings_builder
