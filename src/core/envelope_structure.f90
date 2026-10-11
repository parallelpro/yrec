!----------------------------------------------------------------------
! envelope_structure
!----------------------------------------------------------------------
! The envelope of the current model, point by point (1..num_env_points):
! log P, log T, log m, log rho, log r, X, Z, the convective flag, the
! gradients, the convective velocity, beta, gamma1, cp, the ionization
! fractions, the opacity and its derivatives, the luminosity and
! d ln rho/d ln T. The envelope integration (envelope_integration_kernel)
! fills it; read_starting_model, envelope_refit, rebuild_envelope and
! the stitched model (and its pulse files) read it.
module envelope_structure
      use star_info_lib, only: json
      implicit none
      private

      type, public :: envelope_structure_state
            double precision :: env_log10_pressure(json), &
                 env_log10_temperature(json)
            double precision :: env_log10_mass(json), env_log10_density(json), &
                 env_log10_radius(json)
            double precision :: env_hydrogen_fraction(json), &
                 env_metal_fraction(json)
            logical :: env_convective_flag(json)
            double precision :: env_gradients(3,json), &
                 env_convective_velocity(json), env_beta(json)
            double precision :: env_gamma1(json), env_specific_heat_cp(json), &
                 env_ion_fraction(3,json)
            double precision :: env_opacity(json), env_luminosity(json), &
                 env_dlnrho_dlnt(json)
! the opacity derivatives the envelope
! integration computes at every point (envelope_derivs: pt_scr%qqod,
! pt_scr%qqot), kept so that the pulse writers can fill the GYRE/GSM
! kap_kap_T and kap_kap_rho columns of the envelope points (they were
! written as zero, which switches off the kappa mechanism of the
! hydrogen ionisation zone in any nonadiabatic calculation on the file).
            double precision :: env_dlnkap_dlnrho(json), env_dlnkap_dlnt(json)
            integer :: num_env_points
      end type envelope_structure_state

      type(envelope_structure_state), public :: env_struct

end module envelope_structure
