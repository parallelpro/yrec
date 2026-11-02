!----------------------------------------------------------------------
! pulsation_output
!----------------------------------------------------------------------
! The pulsation-code model files, written by yrec_output. Their column
! definitions are ported from MESA (mesa-26.04.1/star/private/
! pulse_*.f90). write_fgong_pulse and write_gsm_pulse take the extended
! point set from build_pulse_points (interior, envelope and atmosphere,
! centre to the top of the atmosphere; the column schema of pts is
! documented there); M_star, R_star and L_star refer to the photosphere,
! so atmosphere points sit above R_star, as in MESA's add_atmosphere.
! write_gyre_pulse writes the GYRE text file, both the one that goes
! with each profile (the extended point set) and the interior-only one
! every pulse_gyre_interval models (the point set's first nz points;
! yrec_output passes them). GSM output needs HDF5 (make USE_HDF5=1);
! without it, write_gsm_pulse is a stub and gsm_supported is .false.
module pulsation_output
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use run_messages, only: check_io
      use output_files, only: open_replacement_file, close_replacement_file
      implicit none
      private
      public :: write_fgong_pulse, write_gsm_pulse, write_gyre_pulse, gsm_supported

contains

! ---------------------------------------------------------------
! FGONG (ivers 1300: 15 globals, 40 columns, five values per line in
! the 1P,5(X,E26.18E3) layout), points from the surface to the centre.
! Quantities YREC does not track per zone (the central d2P/d2rho terms,
! the dlnGamma1 partials, X_Be7, X_Ne20, r_X) are zero, as MESA writes
! several of them. var(18) = R_star - r is negative for atmosphere
! points.
!   star reads: age_gyr, card_initial_x, card_initial_z, log_Teff,
!     mixing_length_alpha, model_number
!   star writes: none
subroutine write_fgong_pulse(s, n, pts, mstar_g, rstar_cm, lstar_cgs, &
     pulse_path, ierr)
      use star_info_lib, only: star_info
      use physical_constants, only: cgl, ln10
      use math_lib
      use stitched_model, only: n_pulse_cols, ipul_r, ipul_m, &
           ipul_L, ipul_P, ipul_T, ipul_rho, ipul_grad, ipul_N2, &
           ipul_gamma1, ipul_grada, ipul_delta, ipul_kap, ipul_eps, &
           ipul_cp, ipul_mu_e_inv, ipul_h1, ipul_z, &
           ipul_species_base, ipul_eps_grav
      type(star_info), intent(in) :: s
      integer, intent(in) :: n
      real(dp), intent(in) :: pts(:,:)
      real(dp), intent(in) :: mstar_g, rstar_cm, lstar_cgs
      character(len=*), intent(in) :: pulse_path
      integer, intent(inout) :: ierr

      integer, parameter :: iconst = 15, ivar = 40, ivers = 1300
! the globals' and each point's values, five per line
      character(len=*), parameter :: fgong_record_format = '(1P,5(X,E26.18E3))'
      real(dp) :: glob(iconst), var(ivar)
      integer :: u, j, k, i, ios
      character(len=256) :: msg
      real(dp) :: radius_cm, mass_g, grav
      real(dp) :: grav_const_cgs

      call open_replacement_file(pulse_path, u, ierr)
      if (ierr /= 0) return

      grav_const_cgs = exp(ln10*cgl)
      glob = 0.0d0
      glob(1) = mstar_g
      glob(2) = rstar_cm
      glob(3) = lstar_cgs
! the run's initial Z and X, the first card's starting composition
! (card_initial_z(1) and _x(1), the profile header's initial_z too; the
! history header's initial_z and initial_x are the same unless card 1
! rescales the model to another composition), and the model's mixing
! length (the current card's), as MESA writes s%initial_z,
! 1 - initial_y - initial_z and s%mixing_length_alpha
! (mesa-26.04.1 star/private/pulse_fgong.f90
! :171-176; the FGONG notes call glob(5) X_0, the initial hydrogen
! abundance). (A later card's own X and Z can be unset, 0, on an evolve
! card.) glob(11) and glob(12), the central second derivatives, are not
! computed: zero, as tools tolerate.
      glob(4) = s%card_initial_z(1)
      glob(5) = s%card_initial_x(1)
      glob(6) = s%mixing_length_alpha
      glob(13) = s%age_gyr*1.0d9
      glob(14) = exp10(s%log_Teff)
      glob(15) = grav_const_cgs

! MESA-style 4-line comment header (FGONG spec: lines 1-4 free text),
! followed by the standard `nn iconst ivar ivers` record. ivers 1300
! selects the wide E26.18E3 layout, exactly as MESA writes it.
      write(u,'(a)',iostat=ios,iomsg=msg) 'FGONG file'
      if (ios == 0) write(u,'(a)',iostat=ios,iomsg=msg) 'Created by YREC'
      if (ios == 0) write(u,'(a,i8,a,1pe16.9,a,0pf10.6)',iostat=ios,iomsg=msg) &
           'model_number ', s%model_number, '  star_age_yr ', s%age_gyr*1.0d9, &
           '  log_Teff ', s%log_Teff
      if (ios == 0) write(u,'(a)',iostat=ios,iomsg=msg) ''
      if (ios == 0) write(u,'(4I10)',iostat=ios,iomsg=msg) n, iconst, ivar, ivers
      if (ios == 0) write(u,fgong_record_format,iostat=ios,iomsg=msg) (glob(i), i = 1, iconst)

      do j = 1, n
         if (ios /= 0) exit
         k = n - j + 1
         radius_cm = pts(ipul_r,k)
         mass_g = pts(ipul_m,k)

         var = 0.0d0
         var(1) = radius_cm
         if (mass_g > 0.0d0) then
            var(2) = log(mass_g/mstar_g)
         else
            var(2) = -1.0d99
         end if
         var(3) = pts(ipul_T,k)
         var(4) = pts(ipul_P,k)
         var(5) = pts(ipul_rho,k)
         var(6) = pts(ipul_h1,k)
         var(7) = pts(ipul_L,k)
         var(8) = pts(ipul_kap,k)
         var(9) = pts(ipul_eps,k)
         var(10) = pts(ipul_gamma1,k)
         var(11) = pts(ipul_grada,k)
         var(12) = pts(ipul_delta,k)
         var(13) = pts(ipul_cp,k)
         var(14) = pts(ipul_mu_e_inv,k)
         if (radius_cm > 0.0d0) then
            grav = grav_const_cgs*mass_g/(radius_cm*radius_cm)
            var(15) = pts(ipul_N2,k)*radius_cm/grav
         end if
! var(16) r_X: not tracked -> 0
         var(17) = pts(ipul_z,k)
         var(18) = rstar_cm - radius_cm
         var(19) = pts(ipul_eps_grav,k)
! var(20) unused
         var(21) = pts(ipul_species_base+1,k)
         var(22) = pts(ipul_species_base+2,k)
         var(23) = pts(ipul_species_base+3,k)
         var(24) = pts(ipul_species_base+4,k)
         var(25) = pts(ipul_species_base+5,k)
! var(26)-var(28) dlnGamma1 partials: not tracked -> 0
         var(29) = pts(ipul_species_base+6,k)
         var(30) = pts(ipul_species_base+7,k)
         var(31) = pts(ipul_species_base+8,k)
! var(32) X_Be7: YREC tracks Be9, not Be7 -> 0
         var(33) = pts(ipul_species_base+9,k)
         var(34) = pts(ipul_species_base+10,k)
         var(35) = pts(ipul_species_base+11,k)
! var(36) X_Ne20 and var(37)-var(40): not tracked -> 0

         write(u,fgong_record_format,iostat=ios,iomsg=msg) (var(i), i = 1, ivar)
      end do
      call close_replacement_file(u, pulse_path, ios, msg, ierr)
end subroutine write_fgong_pulse

! ---------------------------------------------------------------
! GSM, the GYRE Stellar Model HDF5 format (v1.10): root attributes n,
! M_star, R_star, L_star and version, and one double dataset per column,
! the same 18 columns as the GYRE text file from the same point set.
#ifdef YREC_USE_HDF5
subroutine write_gsm_pulse(n, pts, mstar_g, rstar_cm, lstar_cgs, &
     pulse_path, ierr)
      use hdf5
      use h5lt
      use stitched_model, only: ipul_gyre_last

      integer, intent(in) :: n
      real(dp), intent(in) :: pts(:,:)
      real(dp), intent(in) :: mstar_g, rstar_cm, lstar_cgs
      character(len=*), intent(in) :: pulse_path
      integer, intent(inout) :: ierr

      ! GSM schema version, stamped in the root 'version' attribute.
      ! 110 = GSM v1.10, the schema whose 18-column layout (kap_kap_T,
      ! eps_eps_rho, ... naming) is written below -- matching MESA's
      ! pulse_gsm.f90 case(110).  (101 is the *text* GYRE format's schema
      ! number and is not a valid GSM version; GYRE aborts on it.)
      integer, parameter :: gyre_schema = 110
      integer :: i, herr, close_err
      integer(hid_t) :: file_id
      integer(hsize_t) :: dims1(1)
      character(len=12) :: dset_names(18)

      dset_names = [character(len=12) :: 'r', 'M_r', 'L_r', 'P', 'T', &
           'rho', 'nabla', 'N2', 'Gamma_1', 'nabla_ad', 'delta', 'kap', &
           'kap_kap_T', 'kap_kap_rho', 'eps', 'eps_eps_T', &
           'eps_eps_rho', 'Omega_rot']

! Each call's status (herr, 0 on success) is checked: after the first
! that fails the later writes are skipped, the file and the library are
! closed, and the failure is reported.
      call h5open_f(herr)
      call check_io(herr, 'HDF5 could not start', 'write', pulse_path, ierr)
      if (ierr /= 0) return
      call h5fcreate_f(trim(pulse_path), H5F_ACC_TRUNC_F, file_id, herr)
      if (herr /= 0) then
         call h5close_f(close_err)
         call check_io(herr, 'HDF5 could not create the file', 'write', pulse_path, ierr)
         return
      end if

      call h5ltset_attribute_int_f(file_id, '/', 'n', &
           [n], 1_size_t, herr)
      if (herr == 0) call h5ltset_attribute_double_f(file_id, '/', 'M_star', &
           [mstar_g], 1_size_t, herr)
      if (herr == 0) call h5ltset_attribute_double_f(file_id, '/', 'R_star', &
           [rstar_cm], 1_size_t, herr)
      if (herr == 0) call h5ltset_attribute_double_f(file_id, '/', 'L_star', &
           [lstar_cgs], 1_size_t, herr)
      if (herr == 0) call h5ltset_attribute_int_f(file_id, '/', 'version', &
           [gyre_schema], 1_size_t, herr)

      dims1(1) = int(n, hsize_t)
      do i = 1, ipul_gyre_last
         if (herr /= 0) exit
         call h5ltmake_dataset_double_f(file_id, trim(dset_names(i)), &
              1, dims1, pts(i,1:n), herr)
      end do

      call h5fclose_f(file_id, close_err)
      if (herr == 0) herr = close_err
      call h5close_f(close_err)
      if (herr == 0) herr = close_err
      call check_io(herr, 'HDF5 could not write the file', 'write', pulse_path, ierr)
end subroutine write_gsm_pulse
#else
subroutine write_gsm_pulse(n, pts, mstar_g, rstar_cm, lstar_cgs, &
     pulse_path, ierr)
      use run_messages, only: report_error
      integer, intent(in) :: n
      real(dp), intent(in) :: pts(:,:)
      real(dp), intent(in) :: mstar_g, rstar_cm, lstar_cgs
      character(len=*), intent(in) :: pulse_path
      integer, intent(inout) :: ierr

! (output_open rejects GSM on this build, so this is a last guard)
      call report_error('pulse_format = GSM requires an HDF5-enabled build: ' // &
           'make clean && make USE_HDF5=1 (HDF5_DIR defaults to the MESA SDK)')
      ierr = 1
end subroutine write_gsm_pulse
#endif

! ---------------------------------------------------------------
! Whether this build writes GSM files, so output_open can reject
! pulse_format = 'GSM' at startup (the stub above is a last guard).
logical function gsm_supported()
#ifdef YREC_USE_HDF5
      gsm_supported = .true.
#else
      gsm_supported = .false.
#endif
end function gsm_supported

! ---------------------------------------------------------------
! GYRE text format (schema 101, 18 columns), points 1..n of a pulse
! point set, centre outward: the file written with each profile when
! pulse_format is 'GYRE' (the extended point set, M_star, R_star and
! L_star the photosphere's), and the interior-only file every
! pulse_gyre_interval models (the interior points, M_star, R_star and
! L_star the outermost interior shell's).
subroutine write_gyre_pulse(n, pts, mstar_g, rstar_cm, lstar_cgs, path, ierr)
      use stitched_model, only: ipul_gyre_last
      integer, intent(in) :: n
      real(dp), intent(in) :: pts(:,:)
      real(dp), intent(in) :: mstar_g, rstar_cm, lstar_cgs
      character(len=*), intent(in) :: path
      integer, intent(inout) :: ierr
      integer, parameter :: gyre_schema = 101
      integer :: u, j, i, ios
      character(len=256) :: msg

      call open_replacement_file(path, u, ierr)
      if (ierr /= 0) return
      write(u,'(I6,3(1X,1PE26.16),1X,I6)',iostat=ios,iomsg=msg) n, mstar_g, rstar_cm, &
           lstar_cgs, gyre_schema
      do j = 1, n
         if (ios /= 0) exit
         write(u,'(I6,99(1X,1PE26.16))',iostat=ios,iomsg=msg) j, (pts(i,j), i = 1, ipul_gyre_last)
      end do
      call close_replacement_file(u, path, ios, msg, ierr)
end subroutine write_gyre_pulse

end module pulsation_output
