!----------------------------------------------------------------------
! write_gsm_pulse
!----------------------------------------------------------------------
! New (2026, MESA-style output): writes the converged model as a GSM
! (GYRE Stellar Model) HDF5 file -- root attributes n / M_star /
! R_star / L_star / version plus one double dataset per column, the
! schema ported from mesa-26.04.1/star/private/pulse_gsm.f90 (data
! schema 101, the same 18 columns io/write_gyre_pulse.f90 emits as
! text; per-point formulas mirror that writer exactly, so the text
! and HDF5 products are numerically identical).
!
! HDF5 is an OPTIONAL dependency: build with `make USE_HDF5=1`
! (HDF5_DIR defaults to the MESA SDK, which ships static HDF5 +
! Fortran modules). Without it this file compiles to a stub that
! reports the missing capability and stops -- selecting
! pulse_format = 'GSM' is then a configuration error.
#ifdef YREC_USE_HDF5
subroutine write_gsm_pulse(num_shells, model_number, pulse_path)
      use star_info_lib, only: star
      use const_lib
      use hdf5
      use h5lt
      implicit none

      integer, intent(in) :: num_shells, model_number
      character(len=*), intent(in) :: pulse_path

      integer, parameter :: gyre_schema = 101
      double precision, allocatable :: col(:,:)
      integer :: i, k, herr
      integer(hid_t) :: file_id
      integer(hsize_t) :: dims1(1)
      character(len=12) :: dset_names(18)
      double precision :: radius_cm, mass_g, luminosity_erg_s, &
           pressure_cgs, temperature_k, density_cgs, delta, grav, &
           brunt_n2

      dset_names = [character(len=12) :: 'r', 'M_r', 'L_r', 'P', 'T', &
           'rho', 'nabla', 'N2', 'Gamma_1', 'nabla_ad', 'delta', 'kap', &
           'kap_kap_T', 'kap_kap_rho', 'eps', 'eps_eps_T', &
           'eps_eps_rho', 'Omega_rot']

      allocate(col(num_shells, 18))
! Per-point values: same sources and formulas as write_gyre_pulse
! (schema-101 text), center (1) to surface (num_shells).
      do k = 1, num_shells
         radius_cm = exp(ln10*star%logR(k))
         mass_g = star%m(k)
         luminosity_erg_s = star%luminosity_lsun(k)*solar_luminosity_cgs
         pressure_cgs = exp(ln10*star%logP(k))
         temperature_k = exp(ln10*star%logT(k))
         density_cgs = exp(ln10*star%logRho(k))
         delta = -star%pulse%pulse_dlnrho_dlnt(k)
         if (radius_cm > 0.0d0) then
            grav = exp(ln10*cgl)*mass_g/(radius_cm*radius_cm)
            brunt_n2 = grav*grav*(density_cgs/pressure_cgs)*delta* &
                 (star%diag%del_grad(3,k)-star%diag%del_grad(2,k))
         else
            brunt_n2 = 0.0d0
         end if
         col(k,1) = radius_cm
         col(k,2) = mass_g
         col(k,3) = luminosity_erg_s
         col(k,4) = pressure_cgs
         col(k,5) = temperature_k
         col(k,6) = density_cgs
         col(k,7) = star%diag%del_grad(2,k)
         col(k,8) = brunt_n2
         col(k,9) = star%run%adiabatic_index_gamma1(k)
         col(k,10) = star%diag%del_grad(3,k)
         col(k,11) = delta
         col(k,12) = star%diag%so(k)
         col(k,13) = star%pulse%pulse_dlnkap_dlnt(k)
         col(k,14) = star%pulse%pulse_dlnkap_dlnrho(k)
         col(k,15) = star%diag%sesum(k)
         col(k,16) = star%pulse%pulse_dlneps_dlnt(k)
         col(k,17) = star%pulse%pulse_dlneps_dlnrho(k)
         col(k,18) = star%omega(k)
      end do

      call h5open_f(herr)
      call h5fcreate_f(trim(pulse_path), H5F_ACC_TRUNC_F, file_id, herr)
      if (herr /= 0) then
         write(*,*) 'write_gsm_pulse: cannot create ', trim(pulse_path)
         return
      end if

      call h5ltset_attribute_int_f(file_id, '/', 'n', &
           [num_shells], 1_size_t, herr)
      call h5ltset_attribute_double_f(file_id, '/', 'M_star', &
           [star%m(num_shells)], 1_size_t, herr)
      call h5ltset_attribute_double_f(file_id, '/', 'R_star', &
           [exp(ln10*star%logR(num_shells))], 1_size_t, herr)
      call h5ltset_attribute_double_f(file_id, '/', 'L_star', &
           [star%luminosity_lsun(num_shells)*solar_luminosity_cgs], &
           1_size_t, herr)
      call h5ltset_attribute_int_f(file_id, '/', 'version', &
           [gyre_schema], 1_size_t, herr)

      dims1(1) = int(num_shells, hsize_t)
      do i = 1, 18
         call h5ltmake_dataset_double_f(file_id, trim(dset_names(i)), &
              1, dims1, col(:,i), herr)
      end do

      call h5fclose_f(file_id, herr)
      call h5close_f(herr)
      deallocate(col)
      return
end subroutine write_gsm_pulse
#else
subroutine write_gsm_pulse(num_shells, model_number, pulse_path)
      use luout_lib
      implicit none
      integer, intent(in) :: num_shells, model_number
      character(len=*), intent(in) :: pulse_path

      write(*,*) 'pulse_format = GSM requires an HDF5-enabled build:'
      write(*,*) '  make clean && make USE_HDF5=1'
      write(*,*) '(HDF5_DIR defaults to the MESA SDK; see the Makefile.)'
      write(short_file_unit,*) 'pulse_format = GSM requires USE_HDF5=1 build'
      stop 1
end subroutine write_gsm_pulse
#endif
