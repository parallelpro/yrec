!----------------------------------------------------------------------
! write_gsm_pulse
!----------------------------------------------------------------------
! New (2026, MESA-style output): writes the converged model as a GSM
! (GYRE Stellar Model) HDF5 file -- root attributes n / M_star /
! R_star / L_star / version plus one double dataset per column, the
! schema ported from mesa-26.04.1/star/private/pulse_gsm.f90 (data
! schema 101, the same 18 columns yrec_output's GYRE text writer
! emits; both consume the identical build_pulse_points array, so the
! text and HDF5 products are numerically identical).
!
! Takes the extended point set (interior + envelope + atmosphere,
! center to top of atmosphere); M_star / R_star / L_star refer to
! the photosphere, so atmosphere points sit above R_star, matching
! MESA's add_atmosphere convention. Column schema of pts(35,n) is
! documented at build_pulse_points.
!
! HDF5 is an OPTIONAL dependency: build with `make USE_HDF5=1`
! (HDF5_DIR defaults to the MESA SDK, which ships static HDF5 +
! Fortran modules). Without it this file compiles to a stub that
! reports the missing capability and stops -- selecting
! pulse_format = 'GSM' is then a configuration error.
#ifdef YREC_USE_HDF5
subroutine write_gsm_pulse(n, pts, mstar_g, rstar_cm, lstar_cgs, &
     pulse_path)
      use hdf5
      use h5lt
      implicit none

      integer, intent(in) :: n
      double precision, intent(in) :: pts(35, n)
      double precision, intent(in) :: mstar_g, rstar_cm, lstar_cgs
      character(len=*), intent(in) :: pulse_path

      integer, parameter :: gyre_schema = 101
      integer :: i, herr
      integer(hid_t) :: file_id
      integer(hsize_t) :: dims1(1)
      character(len=12) :: dset_names(18)

      dset_names = [character(len=12) :: 'r', 'M_r', 'L_r', 'P', 'T', &
           'rho', 'nabla', 'N2', 'Gamma_1', 'nabla_ad', 'delta', 'kap', &
           'kap_kap_T', 'kap_kap_rho', 'eps', 'eps_eps_T', &
           'eps_eps_rho', 'Omega_rot']

      call h5open_f(herr)
      call h5fcreate_f(trim(pulse_path), H5F_ACC_TRUNC_F, file_id, herr)
      if (herr /= 0) then
         write(*,*) 'write_gsm_pulse: cannot create ', trim(pulse_path)
         return
      end if

      call h5ltset_attribute_int_f(file_id, '/', 'n', &
           [n], 1_size_t, herr)
      call h5ltset_attribute_double_f(file_id, '/', 'M_star', &
           [mstar_g], 1_size_t, herr)
      call h5ltset_attribute_double_f(file_id, '/', 'R_star', &
           [rstar_cm], 1_size_t, herr)
      call h5ltset_attribute_double_f(file_id, '/', 'L_star', &
           [lstar_cgs], 1_size_t, herr)
      call h5ltset_attribute_int_f(file_id, '/', 'version', &
           [gyre_schema], 1_size_t, herr)

      dims1(1) = int(n, hsize_t)
      do i = 1, 18
         call h5ltmake_dataset_double_f(file_id, trim(dset_names(i)), &
              1, dims1, pts(i,1:n), herr)
      end do

      call h5fclose_f(file_id, herr)
      call h5close_f(herr)
      return
end subroutine write_gsm_pulse
#else
subroutine write_gsm_pulse(n, pts, mstar_g, rstar_cm, lstar_cgs, &
     pulse_path)
      use luout_lib
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: pts(35, n)
      double precision, intent(in) :: mstar_g, rstar_cm, lstar_cgs
      character(len=*), intent(in) :: pulse_path

      write(*,*) 'pulse_format = GSM requires an HDF5-enabled build:'
      write(*,*) '  make clean && make USE_HDF5=1'
      write(*,*) '(HDF5_DIR defaults to the MESA SDK; see the Makefile.)'
      write(short_file_unit,*) 'pulse_format = GSM requires USE_HDF5=1 build'
      stop 1
end subroutine write_gsm_pulse
#endif

! ---------------------------------------------------------------
! 2026 io-writer stops -> ierr: capability probe so output_init_mesa
! can reject pulse_format='GSM' at configuration time (through the
! read_controls ierr chain) instead of the stub stopping at the
! first pulse write. The stub's stop above remains as a last-resort
! guard only.
logical function gsm_supported()
#ifdef YREC_USE_HDF5
      gsm_supported = .true.
#else
      gsm_supported = .false.
#endif
end function gsm_supported
