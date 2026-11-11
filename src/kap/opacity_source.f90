!----------------------------------------------------------------------
! opacity_source
!----------------------------------------------------------------------
! What an opacity table family is to kap_eval (docs/TARGET_ARCHITECTURE.md
! 4.5(c)): a kap_source descriptor, the point it is evaluated at
! (kap_point) and what it returns (kap_result). Each family module
! returns its descriptor (opal95_source, ...), with the family's
! procedures behind four interfaces; opacity_families keeps the
! families of each role in precedence order.
!   kap_point     log rho, log T, X and Z (log is log10, CGS), and the
!                 ionized fractions of H, He+ and He++ that the
!                 conductive opacity weights its species by (the eos's)
!   kap_result    kappa (cm^2/g), log10 kappa, dln kappa/dln rho at
!                 constant T and dln kappa/dln T at constant rho
!   kap_source    a family: its name, whether the settings select it,
!                 its window in log T (set by its role), whether it
!                 covers any Z itself or blends two tables in Z (and
!                 their Zs), and its procedures:
!     load        read its tables and build what its lookups need
!     refresh     rebuild its surface-composition tables at a new
!                 envelope X (null when it has none)
!     lookup      the opacity at a point; in_table is .false. where the
!                 family has no value there (Kurucz 1990 outside its
!                 table), and the next role's family is used alone
!     core_lookup an interior family's opacity in a helium-burning core
!                 from its single-Z table, and that table's Z, which
!                 kap_eval blends toward the pure-Z table (null for the
!                 other roles)
module opacity_source
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: kap_settings_t, kap_table_paths_t, opacity_table_state
      implicit none
      private
      public :: kap_point, kap_result, kap_source, source_load, source_refresh, &
           source_lookup, source_core_lookup

      type :: kap_point
         real(dp) :: log10_density, log10_temperature, hydrogen_fraction, &
              metal_fraction
         real(dp) :: ion_fraction(3)
      end type kap_point

      type :: kap_result
         real(dp) :: opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt
      end type kap_result

      abstract interface
         subroutine source_load(settings, paths, envelope_hydrogen_fraction, tables, ierr)
            import :: dp, kap_settings_t, kap_table_paths_t, opacity_table_state
            type(kap_settings_t), intent(in) :: settings
            type(kap_table_paths_t), intent(in) :: paths
            real(dp), intent(in) :: envelope_hydrogen_fraction
            type(opacity_table_state), intent(inout) :: tables
            integer, intent(out) :: ierr
         end subroutine source_load

         subroutine source_refresh(settings, hydrogen_fraction, tables, ierr)
            import :: dp, kap_settings_t, opacity_table_state
            type(kap_settings_t), intent(in) :: settings
            real(dp), intent(in) :: hydrogen_fraction
            type(opacity_table_state), intent(inout) :: tables
            integer, intent(out) :: ierr
         end subroutine source_refresh

         subroutine source_lookup(settings, tables, point, res, in_table, ierr)
            import :: kap_settings_t, opacity_table_state, kap_point, kap_result
            type(kap_settings_t), intent(in) :: settings
            type(opacity_table_state), intent(inout) :: tables
            type(kap_point), intent(in) :: point
            type(kap_result), intent(inout) :: res
            logical, intent(out) :: in_table
            integer, intent(out) :: ierr
         end subroutine source_lookup

         subroutine source_core_lookup(settings, tables, point, res, table_z, ierr)
            import :: dp, kap_settings_t, opacity_table_state, kap_point, kap_result
            type(kap_settings_t), intent(in) :: settings
            type(opacity_table_state), intent(inout) :: tables
            type(kap_point), intent(in) :: point
            type(kap_result), intent(inout) :: res
            real(dp), intent(out) :: table_z
            integer, intent(out) :: ierr
         end subroutine source_core_lookup
      end interface

      type :: kap_source
         character(len=32) :: name = 'none'
         logical :: selected = .false.
! the family's window in log T (its role's, from the ramp's settings)
         real(dp) :: logt_min = -huge(1.0d0), logt_max = huge(1.0d0)
! .true. when it interpolates in Z itself (OPAL 1995), so it needs no
! table at the zone's Z
         logical :: covers_any_z = .false.
! .true. when, with use_two_z_tables, it reads a second table and blends
! the two linearly in Z: their Zs and the controls that set them
         logical :: blends_two_z = .false.
         real(dp) :: table_z(2) = 0.0d0
         character(len=16) :: table_z_controls(2) = ''
         procedure(source_load), pointer, nopass :: load => null()
         procedure(source_refresh), pointer, nopass :: refresh => null()
         procedure(source_lookup), pointer, nopass :: lookup => null()
         procedure(source_core_lookup), pointer, nopass :: core_lookup => null()
      end type kap_source

end module opacity_source
