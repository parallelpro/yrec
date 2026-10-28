!----------------------------------------------------------------------
! transport_mechanism_layout
!----------------------------------------------------------------------
! The slots of rotation's transport mechanisms. Each mechanism has a
! velocity at every interface, star%mechanism_velocity(i_mech_*, k)
! (cm/s), and a velocity scale and a mixing efficiency,
! rotation_settings%velocity_scale(i_mech_*) and %mixing_scale(i_mech_*).
! The diffusion coefficients sum the mechanisms in slot order, so a new
! mechanism takes the next slot: the sum of the present ones is then
! unchanged bit for bit. A constants-only module, so that the star and
! rotation's kernels share the slots without either using the other
! (star_info_lib re-exports them).
module transport_mechanism_layout
      implicit none
      private
      integer, parameter, public :: i_mech_eddington_sweet = 1, i_mech_gsf = 2, &
           i_mech_shear = 3, n_transport_mechanisms = 3
      character(len=16), parameter, public :: &
           transport_mechanism_name(n_transport_mechanisms) = &
           [character(len=16) :: 'Eddington-Sweet', 'GSF', 'shear']
end module transport_mechanism_layout
