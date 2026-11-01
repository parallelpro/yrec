!----------------------------------------------------------------------
! kap_settings_builder
!----------------------------------------------------------------------
! kap's settings record and its tables' paths, built from the run's
! controls (setup/README.md). settings_install's install_run_settings
! installs them once per run, with kap_init, which also loads the tables:
!     kap_settings_from_star          kap_settings_t
!     kap_table_paths_from_star       kap_table_paths_t
module kap_settings_builder
      use kap_lib, only: kap_settings_t, kap_table_paths_t
      implicit none
      private
      public :: kap_settings_from_star, kap_table_paths_from_star

contains

!----------------------------------------------------------------------
! kap_settings_from_star: kap's settings from the controls. One field is
! derived: use_two_z_tables, a second table set at another Z, which the
! Z ramp (control%use_z_ramp) or Z diffusion (physics%use_diffusion_z)
! needs.
!   star reads: control%alex_table_z1, control%kurucz_table_z1,
!     control%kurucz_table_z2, control%laol_table_z1,
!     control%laol_table_z2, control%opal95_single_table_z,
!     control%opal_table_z1, control%opal_table_z2,
!     control%use_alex06_tables, control%use_alex95_tables,
!     control%use_kurucz90_tables, control%use_laol89_tables,
!     control%use_opal92_tables, control%use_opal95_tables,
!     control%use_pure_z_table, control%use_z_ramp,
!     physics%laol_extrapolation_tolerance,
!     physics%metal_fraction_match_tolerance,
!     physics%molecular_opacity_logt_max,
!     physics%molecular_opacity_logt_min,
!     physics%use_conductive_opacity, physics%use_diffusion_z
!   star writes: none
pure function kap_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(kap_settings_t) :: settings

      settings%use_opal95_tables = s%control%use_opal95_tables
      settings%use_opal92_tables = s%control%use_opal92_tables
      settings%use_laol89_tables = s%control%use_laol89_tables
      settings%use_pure_z_table = s%control%use_pure_z_table
      settings%use_alex06_tables = s%control%use_alex06_tables
      settings%use_alex95_tables = s%control%use_alex95_tables
      settings%use_kurucz90_tables = s%control%use_kurucz90_tables
      settings%use_conductive_opacity = s%physics%use_conductive_opacity
      settings%use_two_z_tables = s%control%use_z_ramp .or. s%physics%use_diffusion_z
      settings%opal95_single_table_z = s%control%opal95_single_table_z
      settings%opal_table_z1 = s%control%opal_table_z1
      settings%opal_table_z2 = s%control%opal_table_z2
      settings%laol_table_z1 = s%control%laol_table_z1
      settings%laol_table_z2 = s%control%laol_table_z2
      settings%kurucz_table_z1 = s%control%kurucz_table_z1
      settings%kurucz_table_z2 = s%control%kurucz_table_z2
      settings%alex_table_z1 = s%control%alex_table_z1
      settings%molecular_opacity_logt_min = s%physics%molecular_opacity_logt_min
      settings%molecular_opacity_logt_max = s%physics%molecular_opacity_logt_max
      settings%metal_fraction_match_tolerance = s%physics%metal_fraction_match_tolerance
      settings%laol_extrapolation_tolerance = s%physics%laol_extrapolation_tolerance
end function kap_settings_from_star

!----------------------------------------------------------------------
! kap_table_paths_from_star: the tables kap_init reads, from the
! controls.
!   star reads: control%alex06_table_path, control%alex95_table_paths,
!     control%conductive_table_path, control%kurucz_table2_path,
!     control%kurucz_table_path, control%laol_table2_path,
!     control%laol_table_path, control%opal92_table2_path,
!     control%opal92_table_path, control%opal95_table_path,
!     control%pure_z_table_path
!   star writes: none
pure function kap_table_paths_from_star(s) result(paths)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(kap_table_paths_t) :: paths

      paths%opal95 = s%control%opal95_table_path
      paths%opal92 = s%control%opal92_table_path
      paths%opal92_z2 = s%control%opal92_table2_path
      paths%laol = s%control%laol_table_path
      paths%laol_z2 = s%control%laol_table2_path
      paths%pure_z = s%control%pure_z_table_path
      paths%alex06 = s%control%alex06_table_path
      paths%alex94 = s%control%alex95_table_paths
      paths%kurucz = s%control%kurucz_table_path
      paths%kurucz_z2 = s%control%kurucz_table2_path
      paths%conductive = s%control%conductive_table_path
end function kap_table_paths_from_star

end module kap_settings_builder
