!----------------------------------------------------------------------
! species_table_lib
!----------------------------------------------------------------------
! 2026 readability sweep (W3, rotation): the species descriptors that
! the microscopic-diffusion routines used to carry as private
! (weight, charge) lists.  Values are exactly the former literals.
!
! The Thoul et al. (1994) solver works on four species in fixed
! columns: 1 = hydrogen, 2 = helium, 3 = the diffused heavy species
! (iron for the metal/hydrogen runs, one light element at a time for
! lithium diffusion), 4 = electrons.  thoul_col_* name the columns;
! thoul_h1 / thoul_he4 / thoul_fe / thoul_electron hold the descriptor
! values shared by gravitational_settling_setup, microdiff_coefficients
! and microdiff.
!
! Not shared here (they differ value-for-value from the Thoul set and
! stay local): viscos.f90's 11-row kinetic-theory table (H1 1.007825,
! He4 4.0026, Z 1.0/0, ...) and check_composition.f90's four-entry
! mean-weight table (H1 1.007825, He4 4.002603, C12 12.0, He3 3.01603).
!
! Rotation-internal; not a public rotation entry.
module species_table_lib
      use star_info_lib, only: i_li6, i_li7, i_be9
      use phys_const_lib, only: m_electron_amu
      implicit none
      private

! (name, atomic weight in amu, charge in units of e)
      type, public :: species_props
         character(len=4) :: name
         double precision :: weight
         double precision :: charge
      end type species_props

! Thoul column numbers of the four-species vectors
! (atomic_weight(4), atomic_charge(4), species_mass_fraction, concentration,
! the species_fraction(3,json) rows of microdiff).
      integer, parameter, public :: thoul_col_h = 1, thoul_col_he = 2, &
           thoul_col_metal = 3, thoul_col_electron = 4

      type(species_props), parameter, public :: &
           thoul_h1 = species_props('H1', 1.008d0, 1.0d0), &
           thoul_he4 = species_props('He4', 4.004d0, 2.0d0), &
           thoul_fe = species_props('Fe', 55.86d0, 26.0d0), &
           thoul_electron = species_props('e-', m_electron_amu, -1.0d0)

! Light elements diffused one at a time by microdiff when
! ctrl%diffuse_lithium is set, with their composition rows.
      integer, parameter, public :: num_light_diffused = 3
      type(species_props), parameter, public :: &
           light_diffused(num_light_diffused) = [ &
           species_props('Li6', 6.015d0, 3.0d0), &
           species_props('Li7', 7.016d0, 3.0d0), &
           species_props('Be9', 9.012d0, 4.0d0)]
      integer, parameter, public :: &
           light_diffused_row(num_light_diffused) = [i_li6, i_li7, i_be9]

end module species_table_lib
