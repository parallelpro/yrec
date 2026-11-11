!----------------------------------------------------------------------
! eos_table_source
!----------------------------------------------------------------------
! What eos knows of a table source (SCV, OPAL and MHD now; a new source
! is a module that returns one of these, an entry in eos_source_list's
! list, and its settings): its descriptor, eos_table_source_t, and the
! standard procedures a source names in it.
!   name         the source's name (messages, the README)
!   stage        where its result enters (stage_*, below)
!   selected     whether the run's settings select it
!   read_tables  eos_init's read (or open) of its tables, from the paths
!   in_window    a cheap test of whether the point may be in its tables
!                (always_in_window: no such test)
!   lookup       its result at (log T, log P, X, Z) and whether its table
!                covers the point, from a copy of the result so far
!   edge_weight  inside its table, whether the point is near its edge and
!                the weight of its result there (whole_weight_inside: no
!                edge weight, its result is taken whole where it covers)
!   adopt        how its result replaces the result so far inside the
!                table (adopt_table_result, adopt_whole_result)
!   blend        how it is weighted in near the edge (blend_table_result)
! The stages, in the order blend_eos consults them:
!   stage_in_place_of_saha  below the Saha cutoff, in place of the Saha
!                           solution and before the ramp into the fully
!                           ionized gas (SCV)
!   stage_over_the_gas      over the gas's result, after the ramp (OPAL)
!   stage_whole_eos         instead of all of the above: eos_get takes its
!                           result alone (MHD)
! Within a stage the list's order is the order of application, so a
! later source wins where two cover a point.
module eos_table_source
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use eos_results, only: eos_result
      use eos_configuration, only: eos_table_paths_t
      use math_lib
      implicit none
      private
      public :: always_in_window, whole_weight_inside, adopt_table_result, adopt_whole_result, &
           blend_table_result

      integer, parameter, public :: stage_in_place_of_saha = 1, stage_over_the_gas = 2, &
           stage_whole_eos = 3

      abstract interface
         logical function source_selected()
         end function source_selected
         subroutine source_reader(paths, ierr)
            import :: eos_table_paths_t
            type(eos_table_paths_t), intent(in) :: paths
            integer, intent(out) :: ierr
         end subroutine source_reader
         logical function source_window(log10_temperature, r)
            import :: dp, eos_result
            real(dp), intent(in) :: log10_temperature
            type(eos_result), intent(in) :: r
         end function source_window
         subroutine source_lookup(log10_temperature, log10_pressure, hydrogen_fraction, &
              metal_fraction, table_r, covers, ierr)
            import :: dp, eos_result
            real(dp), intent(in) :: log10_temperature, log10_pressure, hydrogen_fraction, &
                 metal_fraction
            type(eos_result), intent(inout) :: table_r
            logical, intent(out) :: covers
            integer, intent(out) :: ierr
         end subroutine source_lookup
         subroutine source_edge_weight(temperature, log10_density, weight, in_table, &
              near_edge, ierr)
            import :: dp
            real(dp), intent(in) :: temperature, log10_density
            real(dp), intent(out) :: weight
            logical, intent(out) :: in_table, near_edge
            integer, intent(out) :: ierr
         end subroutine source_edge_weight
         subroutine source_adopt(table_r, r)
            import :: eos_result
            type(eos_result), intent(in) :: table_r
            type(eos_result), intent(inout) :: r
         end subroutine source_adopt
         subroutine source_blend(table_r, weight, r)
            import :: dp, eos_result
            type(eos_result), intent(in) :: table_r
            real(dp), intent(in) :: weight
            type(eos_result), intent(inout) :: r
         end subroutine source_blend
      end interface

      type, public :: eos_table_source_t
         character(len=16) :: name = ''
         integer :: stage = 0
         procedure(source_selected), pointer, nopass :: selected => null()
         procedure(source_reader), pointer, nopass :: read_tables => null()
         procedure(source_window), pointer, nopass :: in_window => null()
         procedure(source_lookup), pointer, nopass :: lookup => null()
         procedure(source_edge_weight), pointer, nopass :: edge_weight => null()
         procedure(source_adopt), pointer, nopass :: adopt => null()
         procedure(source_blend), pointer, nopass :: blend => null()
      end type eos_table_source_t

contains

! in_window for a source with no cheap window test (its lookup decides)
logical function always_in_window(log10_temperature, r)
      real(dp), intent(in) :: log10_temperature
      type(eos_result), intent(in) :: r
      always_in_window = .true.
end function always_in_window

! edge_weight for a source taken whole wherever its lookup covers the
! point (no edge ramp: SCV, M88)
subroutine whole_weight_inside(temperature, log10_density, weight, in_table, near_edge, ierr)
      real(dp), intent(in) :: temperature, log10_density
      real(dp), intent(out) :: weight
      logical, intent(out) :: in_table, near_edge
      integer, intent(out) :: ierr
      ierr = 0
      weight = 1.0d0
      in_table = .true.
      near_edge = .false.
end subroutine whole_weight_inside

! adopt for a source whose lookup refines a copy of the result so far
! (SCV): its whole record replaces r's
subroutine adopt_whole_result(table_r, r)
      type(eos_result), intent(in) :: table_r
      type(eos_result), intent(inout) :: r
      r = table_r
end subroutine adopt_whole_result

! adopt for a table that gives the density, beta, the mean weights and the
! first derivatives (OPAL): those replace r's, the density from its log;
! the ionization fractions, eta, the second derivatives and the Saha
! solver's fields stay r's.
subroutine adopt_table_result(table_r, r)
      type(eos_result), intent(in) :: table_r
      type(eos_result), intent(inout) :: r
      r%log10_density = table_r%log10_density
      r%density = exp10(r%log10_density)
      r%beta = table_r%beta
      r%beta_inverse = 1.0d0/r%beta
      r%beta14 = 4.0d0*r%beta_inverse*(1.0d0 - r%beta)
      r%specific_gas_constant = table_r%specific_gas_constant
      r%ion_mean_weight_inverse = table_r%ion_mean_weight_inverse
      r%electron_mean_weight_inverse = table_r%electron_mean_weight_inverse
      r%dlnrho_dlnt = table_r%dlnrho_dlnt
      r%dlnrho_dlnp = table_r%dlnrho_dlnp
      r%specific_heat_cp = table_r%specific_heat_cp
      r%adiabatic_gradient = table_r%adiabatic_gradient
end subroutine adopt_table_result

! blend for such a table near its edge: the same fields, r + w (table -
! r), with w the source's edge weight
subroutine blend_table_result(table_r, weight, r)
      type(eos_result), intent(in) :: table_r
      real(dp), intent(in) :: weight
      type(eos_result), intent(inout) :: r
      r%log10_density = r%log10_density + &
           weight*(table_r%log10_density - r%log10_density)
      r%density = exp10(r%log10_density)
      r%beta = r%beta + weight*(table_r%beta - r%beta)
      r%beta_inverse = 1.0d0/r%beta
      r%beta14 = 4.0d0*r%beta_inverse*(1.0d0 - r%beta)
      r%specific_gas_constant = r%specific_gas_constant + &
           weight*(table_r%specific_gas_constant - r%specific_gas_constant)
      r%ion_mean_weight_inverse = r%ion_mean_weight_inverse + &
           weight*(table_r%ion_mean_weight_inverse - r%ion_mean_weight_inverse)
      r%electron_mean_weight_inverse = r%electron_mean_weight_inverse + &
           weight*(table_r%electron_mean_weight_inverse - r%electron_mean_weight_inverse)
      r%dlnrho_dlnt = r%dlnrho_dlnt + &
           weight*(table_r%dlnrho_dlnt - r%dlnrho_dlnt)
      r%dlnrho_dlnp = r%dlnrho_dlnp + &
           weight*(table_r%dlnrho_dlnp - r%dlnrho_dlnp)
      r%specific_heat_cp = r%specific_heat_cp + &
           weight*(table_r%specific_heat_cp - r%specific_heat_cp)
      r%adiabatic_gradient = r%adiabatic_gradient + &
           weight*(table_r%adiabatic_gradient - r%adiabatic_gradient)
end subroutine blend_table_result

end module eos_table_source
