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
           integer_column, prof_default_names
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
subroutine profile_output_init(s, profiles, ierr)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(profile_selection), intent(out) :: profiles
      integer, intent(out) :: ierr

      call parse_columns(s%ctrl%profile_columns_file, profile_columns(), &
           prof_default_names, profiles%sel, profiles%nsel, 'profile', ierr)
end subroutine profile_output_init

! ---------------------------------------------------------------
! Profile iprof, <out_dir><profile_data_prefix><iprof>.data, of the
! stitched model.
subroutine write_profile(s, profiles, out_dir, iprof, stitched, ierr)
      use star_info_lib, only: star_info
      use stitched_model_state, only: stitched_model_record
      use io_units, only: check_io
      type(star_info), intent(in) :: s
      type(profile_selection), intent(in) :: profiles
      character(len=*), intent(in) :: out_dir
      integer, intent(in) :: iprof
! stitched: the stitched model whose columns are written.
      type(stitched_model_record), intent(in) :: stitched
      integer, intent(inout) :: ierr
      type(output_column) :: cols(n_prof_cols)
      character(len=256) :: path, msg
      character(len=12) :: numstr
      integer :: u, i, k, kz, ios
      real(dp) :: v

      cols = profile_columns()
      write(numstr, '(i0)') iprof
      path = trim(out_dir) // trim(s%ctrl%profile_data_prefix) // trim(numstr) // '.data'
      open(newunit=u, file=path, status='REPLACE', action='WRITE', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'write', path, ierr)
      if (ierr /= 0) return
! global block (MESA profile shape: numbers / names / values)
      write(u, '(8(1x,i40))') 1, 2, 3, 4, 5, 6, 7, 8
      write(u, '(8(1x,a40))') adjustr('model_number'), &
           adjustr('num_zones'), adjustr('initial_mass'), &
           adjustr('initial_z'), adjustr('star_age'), &
           adjustr('time_step'), adjustr('log_Teff'), &
           adjustr('star_mass')
! initial_mass/initial_z are the kind card's starting values
! (pulsation_mass_msun is stamped by begin_kind_card); star_age and
! time_step are in years, matching MESA's profile header.
      write(u, '(2(1x,i40),6(1x,es40.16e3))') s%model_number, &
           stitched%n_ext, s%pulsation_mass_msun, &
           s%job%initial_envelope_z, s%age_gyr*1.0d9, &
           s%timestep_yr, s%log_Teff, s%star_mass
      write(u, '(a)') ''
      write(u, '(999(1x,i40))') (k, k = 1, profiles%nsel)
      write(u, '(999(1x,a40))') &
           (adjustr(trim(cols(profiles%sel(k))%name)), k = 1, profiles%nsel)
! rows: zone 1 = the outermost point (MESA convention). The extended
! model (interior + envelope + atmosphere, from stitched_model)
! is stored inward-to-outward, so walk it in reverse.
      do k = 1, stitched%n_ext
         kz = stitched%n_ext - k + 1
         do i = 1, profiles%nsel
            if (profiles%sel(i) == ip_zone) then
               v = dble(k)
            else
               v = stitched%prof(profiles%sel(i), kz)
            end if
            if (cols(profiles%sel(i))%kind == integer_column) then
               write(u, '(1x,i40)', advance='no') nint(v)
            else
               write(u, '(1x,es40.16e3)', advance='no') v
            end if
         end do
         write(u, '(a)') ''
      end do
      close(u)
end subroutine write_profile

end module profile_output
