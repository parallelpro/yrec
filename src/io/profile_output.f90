!----------------------------------------------------------------------
! profile_output
!----------------------------------------------------------------------
! The MESA-layout profile files (split out of yrec_output, 2026):
! the column name table and the writer for <outdir>/profile<N>.data
! (zone 1 = the surface); columns selectable via the
! profile_columns_file control (default selection from
! defaults/profile_columns.list -- see output_columns_lib).
!
! Writers are PURE READERS of the stitched full-star model
! (interior + envelope + atmosphere), assembled every step by
! core/stitched_model.f90 into the stx_* arrays. Module state below
! is per-job (output directory, parsed selection), reset by
! profile_output_init on every job.
module profile_output
      use stitched_model_lib, only: n_ext, stx_prof, n_prof_cols
      use output_columns_lib, only: max_cols, parse_columns, &
           prof_default_names, n_prof_default
      implicit none
      private
      public :: profile_output_init, write_profile

      character(len=256) :: prof_out_dir = ' '
      integer :: prof_nsel = 0
      integer :: prof_sel(max_cols)

contains

! ---------------------------------------------------------------
subroutine profile_output_init(out_dir, ierr)
      use star_info_lib, only: star
      character(len=*), intent(in) :: out_dir
      integer, intent(out) :: ierr
      character(len=24) :: prof_names(n_prof_cols)

      prof_out_dir = out_dir
      ierr = 0
      call profile_column_names(prof_names)
      call parse_columns(star%ctrl%profile_columns_file, prof_names, &
           n_prof_cols, prof_default_names, n_prof_default, &
           prof_sel, prof_nsel, 'profile', ierr)
end subroutine profile_output_init

! ---------------------------------------------------------------
subroutine profile_column_names(names)
      character(len=24), intent(out) :: names(n_prof_cols)

      names(1) = 'zone'
      names(2) = 'mass'
      names(3) = 'logR'
      names(4) = 'logT'
      names(5) = 'logRho'
      names(6) = 'logP'
      names(7) = 'luminosity'
      names(8) = 'dm'
      names(9) = 'mixing_type'
      names(10) = 'gamma1'
      names(11) = 'opacity'
      names(12) = 'gradr'
      names(13) = 'gradT'
      names(14) = 'grada'
      names(15) = 'conv_vel'
      names(16) = 'beta'
      names(17) = 'eta'
      names(18) = 'grav'
      names(19) = 'eps_ppI'
      names(20) = 'eps_ppII'
      names(21) = 'eps_ppIII'
      names(22) = 'eps_cno'
      names(23) = 'eps_he3'
      names(24) = 'eps_nuc'
      names(25) = 'eps_neu'
      names(26) = 'eps_grav'
      names(27) = 'h1'
      names(28) = 'h2'
      names(29) = 'he3'
      names(30) = 'he4'
      names(31) = 'li6'
      names(32) = 'li7'
      names(33) = 'be9'
      names(34) = 'c12'
      names(35) = 'c13'
      names(36) = 'n14'
      names(37) = 'n15'
      names(38) = 'o16'
      names(39) = 'o17'
      names(40) = 'o18'
      names(41) = 'z'
      names(42) = 'omega'
      names(43) = 'j_rot'
      names(44) = 'i_rot'
      names(45) = 'fp_rot'
      names(46) = 'ft_rot'
      names(47) = 'v_es'
      names(48) = 'v_gsf'
      names(49) = 'v_ss'
      names(50) = 'cp'
      names(51) = 'delta'
      names(52) = 'mu'
      names(53) = 'mu_e_inv'
      names(54) = 'brunt_N2'
      names(55) = 'lamb_S2'
      names(56) = 'gradL'
      names(57) = 'gradr_div_grada'
      names(58) = 'csound'
      names(59) = 'D_omega'
      names(60) = 'D_mix_rot'
end subroutine profile_column_names

! ---------------------------------------------------------------
subroutine write_profile(iprof)
      use star_info_lib, only: star
      integer, intent(in) :: iprof
      character(len=24) :: names(n_prof_cols)
      character(len=256) :: path
      character(len=12) :: numstr
      integer :: u, i, k, kz
      double precision :: v

      call profile_column_names(names)
      write(numstr, '(i0)') iprof
      path = trim(prof_out_dir) // trim(star%ctrl%profile_data_prefix) // trim(numstr) // '.data'
      open(newunit=u, file=path, status='REPLACE', action='WRITE')
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
      write(u, '(2(1x,i40),6(1x,es40.16e3))') star%model_number, &
           n_ext, star%pulsation_mass_msun, &
           star%job%initial_envelope_z, star%dage*1.0d9, &
           star%timestep_yr, star%log_Teff, star%star_mass
      write(u, '(a)') ''
      write(u, '(999(1x,i40))') (k, k = 1, prof_nsel)
      write(u, '(999(1x,a40))') &
           (adjustr(trim(names(prof_sel(k)))), k = 1, prof_nsel)
! rows: zone 1 = the outermost point (MESA convention). The extended
! model (interior + envelope + atmosphere, from stitched_model_lib)
! is stored inward-to-outward, so walk it in reverse.
      do k = 1, n_ext
         kz = n_ext - k + 1
         do i = 1, prof_nsel
            if (prof_sel(i) == 1) then
               v = dble(k)
            else
               v = stx_prof(prof_sel(i), kz)
            end if
            if (prof_sel(i) == 1 .or. prof_sel(i) == 9) then
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
