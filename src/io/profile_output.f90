!----------------------------------------------------------------------
! profile_output
!----------------------------------------------------------------------
! The profile files, in MESA's profile layout: <outdir>/profile<N>.data,
! zone 1 the surface. The columns are selected by the
! profile_columns_file control (the default selection is
! defaults/profile_columns.list, see output_columns).
!
! The column table and the values are core's (stitched_model: the
! table profile_columns, at the named indices ip_<name>, and the values
! s%stitched%prof over the stitched full-star model: interior, envelope
! and atmosphere). The writer only reads them. The column selection is a
! profile_selection record, part of yrec_output's output_state, which
! run_yrec owns.
module profile_output
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use stitched_model, only: n_prof_cols, profile_columns, ip_zone
      use output_columns, only: max_cols, parse_columns, output_column, &
           integer_column, prof_default_names, column_integer_format, &
           column_real_format, column_numbers_format, column_names_format
      implicit none
      private
      public :: profile_selection, profile_output_init, write_profile

! One job's profile columns (indices into profile_columns).
      type :: profile_selection
         integer :: nsel = 0
         integer :: sel(max_cols) = 0
      end type profile_selection

contains

! ---------------------------------------------------------------
! The job's profile columns, from profile_columns_file.
!   star reads: control%profile_columns_file
!   star writes: none
subroutine profile_output_init(s, profiles, ierr)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(profile_selection), intent(out) :: profiles
      integer, intent(out) :: ierr

      call parse_columns(s%control%profile_columns_file, profile_columns(), &
           prof_default_names, profiles%sel, profiles%nsel, 'profile', ierr)
end subroutine profile_output_init

! ---------------------------------------------------------------
! Profile iprof, <out_dir><profile_data_prefix><iprof>.data, of the
! stitched model (s%stitched).
!   star reads: age_gyr, card_initial_z, log_Teff, model_number,
!     physics%profile_data_prefix, pulsation_mass_msun, star_mass,
!     stitched, timestep_yr
!   star writes: none
subroutine write_profile(s, profiles, out_dir, iprof, ierr)
      use star_info_lib, only: star_info
      use output_files, only: open_output_file, close_output_file
      type(star_info), intent(in) :: s
      type(profile_selection), intent(in) :: profiles
      character(len=*), intent(in) :: out_dir
      integer, intent(in) :: iprof
      integer, intent(inout) :: ierr
      type(output_column) :: cols(n_prof_cols)
      character(len=:), allocatable :: path
      character(len=256) :: msg
      character(len=12) :: numstr
      integer :: u, i, k, kz, ios
      real(dp) :: v

      cols = profile_columns()
      write(numstr, '(i0)') iprof
      path = trim(out_dir) // trim(s%physics%profile_data_prefix) // trim(numstr) // '.data'
      call open_output_file(path, u, ierr)
      if (ierr /= 0) return
! global block (MESA profile shape: numbers / names / values)
      write(u, column_numbers_format, iostat=ios, iomsg=msg) 1, 2, 3, 4, 5, 6, 7, 8
      if (ios == 0) write(u, column_names_format, iostat=ios, iomsg=msg) adjustr('model_number'), &
           adjustr('num_zones'), adjustr('initial_mass'), &
           adjustr('initial_z'), adjustr('star_age'), &
           adjustr('time_step'), adjustr('log_Teff'), &
           adjustr('star_mass')
! initial_mass is the kind card's starting mass (pulsation_mass_msun,
! stamped by begin_kind_card); initial_z is the run's initial Z, the
! first card's starting Z (card_initial_z(1)), as MESA writes it and as
! FGONG's glob(4) is; star_age and time_step are in years, matching
! MESA's profile header.
      if (ios == 0) write(u, '(2(1x,i40),6(1x,es40.16e3))', iostat=ios, iomsg=msg) &
           s%model_number, s%stitched%n_ext, s%pulsation_mass_msun, &
           s%card_initial_z(1), s%age_gyr*1.0d9, &
           s%timestep_yr, s%log_Teff, s%star_mass
      if (ios == 0) write(u, '(a)', iostat=ios, iomsg=msg) ''
      if (ios == 0) write(u, column_numbers_format, iostat=ios, iomsg=msg) (k, k = 1, profiles%nsel)
      if (ios == 0) write(u, column_names_format, iostat=ios, iomsg=msg) &
           (adjustr(trim(cols(profiles%sel(k))%name)), k = 1, profiles%nsel)
! rows: zone 1 = the outermost point (MESA convention). The extended
! model (interior + envelope + atmosphere, from stitched_model)
! is stored inward-to-outward, so walk it in reverse.
      do k = 1, s%stitched%n_ext
         kz = s%stitched%n_ext - k + 1
         do i = 1, profiles%nsel
            if (ios /= 0) exit
            if (profiles%sel(i) == ip_zone) then
               v = dble(k)
            else
               v = s%stitched%prof(profiles%sel(i), kz)
            end if
            if (cols(profiles%sel(i))%kind == integer_column) then
               write(u, column_integer_format, advance='no', iostat=ios, iomsg=msg) nint(v)
            else
               write(u, column_real_format, advance='no', iostat=ios, iomsg=msg) v
            end if
         end do
         if (ios /= 0) exit
         write(u, '(a)', iostat=ios, iomsg=msg) ''
      end do
      call close_output_file(u, path, ios, msg, ierr)
end subroutine write_profile

end module profile_output
