!----------------------------------------------------------------------
! transport_mechanisms
!----------------------------------------------------------------------
! The transport mechanisms' velocities, one slot at a time
! (transport_mechanism_layout). compute_mechanism_velocity calls the
! kernel of the mechanism in a slot (instability_velocities) with what
! it reads from the sub-step's interface record (stab), its physics
! (sub_physics) and omega. circulation_velocities calls it for each
! slot in order, then scales, averages and sums the slots without
! knowing which mechanisms they hold. A mechanism writes its own slot
! and, when the layout declares that it corrects another
! (transport_mechanism_corrects), that one's too, which must come
! earlier: the shear corrects the GSF velocity. A star-blind module
! (tools/check_boundaries.py STAR_BLIND_FILES), checked by
! test_rotation_kernels.
module transport_mechanisms
      use composition_layout, only: json
      use transport_mechanism_layout, only: i_mech_eddington_sweet, i_mech_gsf, &
           i_mech_shear, i_mech_solberg_hoiland, n_transport_mechanisms
      use rotation_configuration, only: rotation_settings
      use rotation_scratch_lib, only: stability_state, substep_physics_state
      use instability_velocities, only: compute_eddington_sweet_velocity, &
           inhibit_by_mu_gradient, compute_gsf_velocity, compute_shear_velocities, &
           compute_solberg_hoiland_velocity
      use rotation_messages, only: rotation_note, message_length, unwritten_line
      implicit none
      private
      public :: compute_mechanism_velocity

contains

! ---------------------------------------------------------------
! The velocity (cm/s) of the mechanism in slot `slot` at the interfaces
! zone_min..zone_max, into mechanism_velocity(slot,:), and into the slot
! it corrects, if any. It reads stab (the interfaces' structure, the
! zones' radii, d omega/d ln r, the dynamical-shear limit and the
! coefficients), sub_physics%mu, omega and the zones' log P and
! convective flags, and the settings gsf_inhibition_mode and
! critical_reynolds.
!   i_mech_eddington_sweet  U_ES corrected for the quadrupole moment
!                           (setting stab%circulation_correction_ratio),
!                           then inhibited by the mu gradient;
!   i_mech_gsf              U_GSF (Kippenhahn), with its inhibition;
!   i_mech_shear            the secular and dynamical shear, which also
!                           corrects U_GSF where both act;
!   i_mech_solberg_hoiland  the Solberg-Hoiland instability, which can
!                           act only where the shear is dynamical.
! A slot outside 1..n_transport_mechanisms does nothing.
subroutine compute_mechanism_velocity(slot, zone_min, zone_max, stab, sub_physics, &
     omega, log_pressure, am_transport_convective_flag, mechanism_velocity)
      integer, intent(in) :: slot, zone_min, zone_max
      type(stability_state), intent(inout) :: stab
      type(substep_physics_state), intent(in) :: sub_physics
      double precision, intent(in) :: omega(json), log_pressure(json)
      logical, intent(in) :: am_transport_convective_flag(json)
      double precision, intent(inout) :: mechanism_velocity(n_transport_mechanisms, json)
! where the shear is dynamical, and the inner omega at the limit there
      logical :: dynamical(json)
      double precision :: dynamical_inner_omega(json)
      integer :: i
      character(len=message_length) :: lines(1)

      select case (slot)
      case (i_mech_eddington_sweet)
         call compute_eddington_sweet_velocity(zone_min, zone_max, am_transport_convective_flag, &
              stab%radius, omega, stab%omega_interface, stab%radius_interface, &
              stab%local_gravity, stab%quadrupole_moment, stab%es_velocity_coeff1, &
              stab%es_velocity_coeff2, stab%es_shear_coeff, &
              stab%circulation_correction_ratio, mechanism_velocity(i_mech_eddington_sweet,:))
         call inhibit_by_mu_gradient(zone_min, zone_max, sub_physics%mu, log_pressure, &
              stab%del_grad_diff_interface, mechanism_velocity(i_mech_eddington_sweet,:))
      case (i_mech_gsf)
         call compute_gsf_velocity(zone_min, zone_max, rotation_settings%gsf_inhibition_mode, &
              am_transport_convective_flag, stab%radius, omega, stab%domega_dlnr, &
              stab%dynamical_shear_omega_limit, &
              sub_physics%mu, stab%mean_molecular_weight_interface, stab%omega_interface, &
              stab%kinematic_viscosity_interface, stab%thermal_diffusivity_interface, &
              stab%gravity_factor_interface, stab%gsf_kippenhahn_coeff, &
              mechanism_velocity(i_mech_gsf,:))
      case (i_mech_shear)
         call compute_shear_velocities(zone_min, zone_max, rotation_settings%critical_reynolds, &
              stab%radius, omega, stab%domega_dlnr, stab%dynamical_shear_omega_limit, &
              sub_physics%mu, stab%mean_molecular_weight_interface, stab%radius_interface, &
              stab%kinematic_viscosity_interface, stab%thermal_diffusivity_interface, &
              stab%mu_gradient_richardson_coeff, mechanism_velocity(i_mech_shear,:), &
              mechanism_velocity(i_mech_gsf,:), dynamical, dynamical_inner_omega)
! The report of each dynamically unstable interface (the kernel is
! pure): the zones' omega and the inner omega at the limit.
         do i = zone_min, zone_max
            if (.not. dynamical(i)) cycle
            lines = unwritten_line
            write(lines,9911) i, omega(i), omega(i-1), dynamical_inner_omega(i)
            call rotation_note(lines)
         end do
 9911    format(1x,'DYNAMICAL SHEAR-SHELL',i5,1p,' WTOP',e11.3, &
              ' WBOT',e11.3,' LIMIT',e11.3)
      case (i_mech_solberg_hoiland)
         call compute_solberg_hoiland_velocity(zone_min, zone_max, &
              am_transport_convective_flag, stab%domega_dlnr, &
              stab%dynamical_shear_omega_limit, sub_physics%mu, &
              stab%mean_molecular_weight_interface, stab%omega_interface, &
              stab%radius_interface, stab%thermal_diffusivity_interface, &
              stab%mu_gradient_richardson_coeff, &
              mechanism_velocity(i_mech_solberg_hoiland,:))
      end select
end subroutine compute_mechanism_velocity

end module transport_mechanisms
