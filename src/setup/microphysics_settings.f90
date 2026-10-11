!----------------------------------------------------------------------
! microphysics_settings
!----------------------------------------------------------------------
! The settings records of the microphysics domains, built from the
! run's controls by the star layer (2026-10, Phase 3). Each domain
! installs its record at init and reads no other configuration:
!     opacity_settings   kap (kap_init)
module microphysics_settings
      use kap_lib, only: kap_settings_t
      implicit none
      private
      public :: opacity_settings

contains

!----------------------------------------------------------------------
! opacity_settings: kap's settings from the controls.
!   star reads: ctrl%alex_table_z1, ctrl%kurucz_table_z1,
!     ctrl%kurucz_table_z2, ctrl%laol_table_z1, ctrl%laol_table_z2,
!     ctrl%metal_fraction_match_tolerance,
!     ctrl%molecular_opacity_logt_max, ctrl%molecular_opacity_logt_min,
!     ctrl%opal95_single_table_z, ctrl%opal_table_z1,
!     ctrl%opal_table_z2, ctrl%tollaol, ctrl%use_alex06_tables,
!     ctrl%use_alex95_tables, ctrl%use_conductive_opacity,
!     ctrl%use_kurucz90_tables, ctrl%use_laol89_tables,
!     ctrl%use_opal92_tables, ctrl%use_opal95_tables,
!     ctrl%use_pure_z_table, use_two_z_tables
!   star writes: none
function opacity_settings() result(s)
      use star_info_lib, only: star
      type(kap_settings_t) :: s

      s%use_opal95_tables = star%ctrl%use_opal95_tables
      s%use_opal92_tables = star%ctrl%use_opal92_tables
      s%use_laol89_tables = star%ctrl%use_laol89_tables
      s%use_pure_z_table = star%ctrl%use_pure_z_table
      s%use_alex06_tables = star%ctrl%use_alex06_tables
      s%use_alex95_tables = star%ctrl%use_alex95_tables
      s%use_kurucz90_tables = star%ctrl%use_kurucz90_tables
      s%use_conductive_opacity = star%ctrl%use_conductive_opacity
      s%use_two_z_tables = star%use_two_z_tables
      s%opal95_single_table_z = star%ctrl%opal95_single_table_z
      s%opal_table_z1 = star%ctrl%opal_table_z1
      s%opal_table_z2 = star%ctrl%opal_table_z2
      s%laol_table_z1 = star%ctrl%laol_table_z1
      s%laol_table_z2 = star%ctrl%laol_table_z2
      s%kurucz_table_z1 = star%ctrl%kurucz_table_z1
      s%kurucz_table_z2 = star%ctrl%kurucz_table_z2
      s%alex_table_z1 = star%ctrl%alex_table_z1
      s%molecular_opacity_logt_min = star%ctrl%molecular_opacity_logt_min
      s%molecular_opacity_logt_max = star%ctrl%molecular_opacity_logt_max
      s%metal_fraction_match_tolerance = star%ctrl%metal_fraction_match_tolerance
      s%laol_extrapolation_tolerance = star%ctrl%tollaol
end function opacity_settings

end module microphysics_settings
