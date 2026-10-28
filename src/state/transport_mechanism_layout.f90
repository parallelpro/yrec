!----------------------------------------------------------------------
! transport_mechanism_layout
!----------------------------------------------------------------------
! The slots of rotation's transport mechanisms, and what each declares:
! its name; the slot whose velocity it corrects, if any (0 for none),
! which must be an earlier slot (the shear corrects the GSF velocity
! where both act); and its two controls, a velocity scale and a mixing
! efficiency, which rotation_settings_setup copies into
! rotation_settings%velocity_scale(i_mech_*) and %mixing_scale(i_mech_*)
! (test_rotation_kernels.py checks that both are registry controls and
! are copied so). Each mechanism has a velocity at every interface,
! star%mechanism_velocity(i_mech_*, k) (cm/s), which
! transport_mechanisms' compute_mechanism_velocity computes slot by
! slot. The diffusion coefficients sum the mechanisms in slot order, so
! a new mechanism takes the next slot: the sum of the present ones is
! then unchanged bit for bit. A constants-only module, so that the star
! and rotation's kernels share the slots without either using the other
! (star_info_lib re-exports them).
module transport_mechanism_layout
      implicit none
      private
      integer, parameter, public :: i_mech_eddington_sweet = 1, i_mech_gsf = 2, &
           i_mech_shear = 3, n_transport_mechanisms = 3
      character(len=16), parameter, public :: &
           transport_mechanism_name(n_transport_mechanisms) = &
           [character(len=16) :: 'Eddington-Sweet', 'GSF', 'shear']
      integer, parameter, public :: &
           transport_mechanism_corrects(n_transport_mechanisms) = [0, 0, i_mech_gsf]
      character(len=32), parameter, public :: &
           transport_mechanism_velocity_scale_control(n_transport_mechanisms) = &
           [character(len=32) :: 'es_velocity_scale', 'gsf_velocity_scale', &
           'secular_shear_velocity_scale'], &
           transport_mechanism_mixing_scale_control(n_transport_mechanisms) = &
           [character(len=32) :: 'es_mixing_scale', 'gsf_mixing_scale', &
           'secular_shear_mixing_scale']
end module transport_mechanism_layout
