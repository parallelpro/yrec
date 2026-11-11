!----------------------------------------------------------------------
! species_table
!----------------------------------------------------------------------
! The species of the settling routines: their atomic weights (amu) and
! charges (e).
!
! The Thoul et al. (1994) solve takes four species in fixed columns: 1
! hydrogen, 2 helium, 3 the heavy species (iron for the metals and when
! hydrogen settles; one light element at a time with diffuse_lithium),
! 4 the electrons. thoul_col_* name the columns; thoul_h1, thoul_he4,
! thoul_fe and thoul_electron are the species both pipelines use. The
! light elements are a table: a fourth is a row here (its weight,
! charge and composition slot) and a count.
!
! Not shared here, because their values differ from these: rotation's
! viscosity table (viscos: H1 1.007825, He4 4.0026, ...) and
! check_composition's mean weights (H1 1.007825, He4 4.002603, C12 12,
! He3 3.01603).
module species_table
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: i_li6, i_li7, i_be9
      use physical_constants, only: m_electron_amu
      implicit none
      private

! (name, atomic weight in amu, charge in units of e)
      type, public :: species_props
         character(len=4) :: name
         real(dp) :: weight
         real(dp) :: charge
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

! The light elements microdiff settles one at a time (diffuse_lithium),
! with their composition rows.
      integer, parameter, public :: num_light_diffused = 3
      type(species_props), parameter, public :: &
           light_diffused(num_light_diffused) = [ &
           species_props('Li6', 6.015d0, 3.0d0), &
           species_props('Li7', 7.016d0, 3.0d0), &
           species_props('Be9', 9.012d0, 4.0d0)]
      integer, parameter, public :: &
           light_diffused_row(num_light_diffused) = [i_li6, i_li7, i_be9]

end module species_table
