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
! and atmosphere). The writer only reads them. The module state below is
! per job (the output directory and the selection), reset by
! profile_output_init.
module profile_output
      use stitched_model, only: n_prof_cols, profile_columns, ip_zone
      use output_columns, only: max_cols, parse_columns, output_column, &
           integer_column, prof_default_names
      implicit none
      private
      public :: profile_output_init, write_profile

      character(len=256) :: prof_out_dir = ' '
      integer :: prof_nsel = 0
      integer :: prof_sel(max_cols)

contains

! ---------------------------------------------------------------
subroutine profile_output_init(s, out_dir, ierr)
      use star_info_lib, only: star_info
      type(star_info), intent(inout) :: s
      character(len=*), intent(in) :: out_dir
      integer, intent(out) :: ierr

      prof_out_dir = out_dir
      ierr = 0
      call parse_columns(s%ctrl%profile_columns_file, profile_columns(), &
           prof_default_names, prof_sel, prof_nsel, 'profile', ierr)
end subroutine profile_output_init

! ---------------------------------------------------------------
subroutine write_profile(s, iprof, stitched, ierr)
      use star_info_lib, only: star_info
      use stitched_model_state, only: stitched_model_record
      use io_units, only: check_io
      type(star_info), intent(inout) :: s
      integer, intent(in) :: iprof
! stitched: the stitched model whose columns are written.
      type(stitched_model_record), intent(in) :: stitched
      integer, intent(inout) :: ierr
      type(output_column) :: cols(n_prof_cols)
      character(len=256) :: path, msg
      character(len=12) :: numstr
      integer :: u, i, k, kz, ios
      double precision :: v

      cols = profile_columns()
      write(numstr, '(i0)') iprof
      path = trim(prof_out_dir) // trim(s%ctrl%profile_data_prefix) // trim(numstr) // '.data'
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
      write(u, '(999(1x,i40))') (k, k = 1, prof_nsel)
      write(u, '(999(1x,a40))') &
           (adjustr(trim(cols(prof_sel(k))%name)), k = 1, prof_nsel)
! rows: zone 1 = the outermost point (MESA convention). The extended
! model (interior + envelope + atmosphere, from stitched_model)
! is stored inward-to-outward, so walk it in reverse.
      do k = 1, stitched%n_ext
         kz = stitched%n_ext - k + 1
         do i = 1, prof_nsel
            if (prof_sel(i) == ip_zone) then
               v = dble(k)
            else
               v = stitched%prof(prof_sel(i), kz)
            end if
            if (cols(prof_sel(i))%kind == integer_column) then
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
