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
! write_gyre_ext writes the GYRE file that goes with each profile, and
! write_gyre_pulse the interior-only GYRE file every pulse_gyre_interval
! models. GSM output needs HDF5 (make USE_HDF5=1); without it,
! write_gsm_pulse is a stub and gsm_supported is .false.
module pulsation_output
      use io_units, only: check_io
      implicit none
      private
      public :: write_fgong_pulse, write_gsm_pulse, write_gyre_pulse, write_gyre_ext, &
           gsm_supported

contains

! ---------------------------------------------------------------
! FGONG (ivers 1300: 15 globals, 40 columns, five values per line in
! the 1P,5(X,E26.18E3) layout), points from the surface to the centre.
! Quantities YREC does not track per zone (the central d2P/d2rho terms,
! the dlnGamma1 partials, X_Be7, X_Ne20, r_X) are zero, as MESA writes
! several of them. var(18) = R_star - r is negative for atmosphere
! points.
subroutine write_fgong_pulse(n, pts, mstar_g, rstar_cm, lstar_cgs, &
     pulse_path, ierr)
      use star_info_lib, only: star
      use physical_constants
      use math_lib
      use stitched_model, only: n_pulse_cols, ipul_r, ipul_m, &
           ipul_L, ipul_P, ipul_T, ipul_rho, ipul_grad, ipul_N2, &
           ipul_gamma1, ipul_grada, ipul_delta, ipul_kap, ipul_eps, &
           ipul_cp, ipul_mu_e_inv, ipul_h1, ipul_z, &
           ipul_species_base, ipul_eps_grav
      implicit none

      integer, intent(in) :: n
      double precision, intent(in) :: pts(n_pulse_cols, n)
      double precision, intent(in) :: mstar_g, rstar_cm, lstar_cgs
      character(len=*), intent(in) :: pulse_path
      integer, intent(inout) :: ierr

      integer, parameter :: iconst = 15, ivar = 40, ivers = 1300
      double precision :: glob(iconst), var(ivar)
      integer :: u, j, k, i, ios
      character(len=256) :: msg
      double precision :: radius_cm, mass_g, grav
      double precision :: grav_const_cgs

      open(newunit=u, file=pulse_path, status='REPLACE', form='FORMATTED', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'write', pulse_path, ierr)
      if (ierr /= 0) return

      grav_const_cgs = exp(ln10*cgl)
      glob = 0.0d0
      glob(1) = mstar_g
      glob(2) = rstar_cm
      glob(3) = lstar_cgs
      glob(4) = pts(ipul_z, n)
      glob(5) = pts(ipul_h1, n)
! glob(6) (mixing length alpha), glob(11)/glob(12) (central second
! derivatives) not tracked here -- zero, as tools tolerate.
      glob(13) = star%age_gyr*1.0d9
      glob(14) = exp10(star%log_Teff)
      glob(15) = grav_const_cgs

! MESA-style 4-line comment header (FGONG spec: lines 1-4 free text),
! followed by the standard `nn iconst ivar ivers` record. ivers 1300
! selects the wide E26.18E3 layout, exactly as MESA writes it.
      write(u,'(a)') 'FGONG file'
      write(u,'(a)') 'Created by YREC'
      write(u,'(a,i8,a,1pe16.9,a,0pf10.6)') 'model_number ', &
           star%model_number, '  star_age_yr ', star%age_gyr*1.0d9, &
           '  log_Teff ', star%log_Teff
      write(u,'(a)') ''
      write(u,'(4I10)') n, iconst, ivar, ivers
      write(u,'(1P,5(X,E26.18E3))') (glob(i), i = 1, iconst)

      do j = 1, n
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

         write(u,'(1P,5(X,E26.18E3))') (var(i), i = 1, ivar)
      end do

      close(u)
      return
end subroutine write_fgong_pulse

! ---------------------------------------------------------------
! GYRE text format (schema 101, 18 columns) from the interior shells,
! centre (1) to surface, from the star%pulse_* arrays the Henyey
! coefficients fill for every shell.
subroutine write_gyre_pulse(num_shells, mass_coordinate, &
     log_density, log_luminosity, log_pressure, log_radius, &
     log_temperature, omega, pulse_path, ierr)
      use star_info_lib, only: star, json
      use physical_constants
      use math_lib
      implicit none

      integer, intent(in) :: num_shells
      double precision, intent(in) :: mass_coordinate(json), &
           log_density(json), log_luminosity(json), log_pressure(json), &
           log_radius(json), log_temperature(json), omega(json)


      character(len=*), intent(in) :: pulse_path
      integer, intent(inout) :: ierr
      integer :: gyre_unit, i, ios
      character(len=256) :: msg
      integer, parameter :: gyre_schema = 101
      double precision :: radius_cm, mass_g, luminosity_erg_s, &
           pressure_cgs, temperature_k, density_cgs, delta, grav, &
           global_data(3)
      double precision :: brunt_n2(json), dr
      double precision :: grav_const_cgs

      open(newunit=gyre_unit,file=pulse_path,status='UNKNOWN',form='FORMATTED', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'write', pulse_path, ierr)
      if (ierr /= 0) return

      global_data(1) = mass_coordinate(num_shells)
      global_data(2) = exp(ln10*log_radius(num_shells))
      global_data(3) = log_luminosity(num_shells)*star%solar_luminosity_cgs
      write(gyre_unit,100) num_shells,global_data,gyre_schema
 100  format(I6,3(1X,1PE26.16),1X,I6)

! Brunt-Vaisala N^2, exact gradient form (centered differences in r):
!     N^2 = g * [ (1/Gamma_1) dlnP/dr - dlnRho/dr ]
! The thermal-only form g^2 rho delta (grada - gradT)/P previously used
! here is exact only for homogeneous composition -- it omits the Ledoux
! mu-gradient term that dominates N^2 in composition-gradient layers
! (see the matching note in core/diagnostics/stitched_model.f90's
! build_pulse_points, which fixes the same formula for the
! profile-coupled GSM/FGONG/GYRE writers). CONVECTIVE shells (where
! the thermal form is negative) keep the thermal value: the mixture is
! homogeneous there, so the thermal form is exact and smooth, while
! the centered difference of a near-adiabatic stratification is
! cancellation noise of random sign.
      grav_const_cgs = exp(ln10*cgl)
      do i = 2, num_shells - 1
         radius_cm = exp(ln10*log_radius(i))
         dr = exp(ln10*log_radius(i+1)) - exp(ln10*log_radius(i-1))
         if (radius_cm > 0.0d0 .and. dr > 0.0d0 .and. &
              star%adiabatic_index_gamma1(i) > 0.0d0) then
            grav = grav_const_cgs*mass_coordinate(i)/(radius_cm*radius_cm)
            brunt_n2(i) = grav*grav* &
                 (exp(ln10*(log_density(i) - log_pressure(i))))* &
                 (-star%pulse_dlnrho_dlnt(i))* &
                 (star%grada(i) - star%gradT(i))
            if (brunt_n2(i) >= 0.0d0) then
               brunt_n2(i) = grav*ln10*( &
                    (log_pressure(i+1) - log_pressure(i-1)) / &
                    star%adiabatic_index_gamma1(i) - &
                    (log_density(i+1) - log_density(i-1)) ) / dr
            end if
         else
            brunt_n2(i) = 0.0d0
         end if
      end do
      if (num_shells >= 2) then
         brunt_n2(1) = brunt_n2(2)
         brunt_n2(num_shells) = brunt_n2(num_shells-1)
      end if

      do i = 1,num_shells
         radius_cm = exp(ln10*log_radius(i))
         mass_g = mass_coordinate(i)
         luminosity_erg_s = log_luminosity(i)*star%solar_luminosity_cgs
         pressure_cgs = exp(ln10*log_pressure(i))
         temperature_k = exp(ln10*log_temperature(i))
         density_cgs = exp(ln10*log_density(i))
! delta = chiT/chiRho = -star%pulse_dlnrho_dlnt, since chiRho=1/star%pulse_dlnrho_dlnp
! and chiT=-chiRho*star%pulse_dlnrho_dlnt (see core/structure/structure_equations.f90).
         delta = -star%pulse_dlnrho_dlnt(i)
! kap_kap_T/eps_eps_T columns follow the GSM/GYRE convention (MESA
! pulse_gyre.f90): the ABSOLUTE derivatives kap*dlnkap/dlnT and
! eps*dlneps/dlnT = deps/dlnT, not the bare log-derivatives.
! pulse_deps_* already hold deps/dln (see core/diagnostics/stitched_model.f90's
! build_pulse_points), so they are not multiplied by eps.
         write(gyre_unit,110) i,radius_cm,mass_g,luminosity_erg_s, &
              pressure_cgs,temperature_k,density_cgs,star%gradT(i), &
              brunt_n2(i),star%adiabatic_index_gamma1(i),star%grada(i),delta, &
              star%opacity_zone(i), &
              star%opacity_zone(i)*star%pulse_dlnkap_dlnt(i), &
              star%opacity_zone(i)*star%pulse_dlnkap_dlnrho(i), &
              star%eps_total(i), &
              star%pulse_deps_dlnt(i), &
              star%pulse_deps_dlnrho(i),omega(i)
 110     format(I6,99(1X,1PE26.16))
      end do

      close(gyre_unit)

      return
end subroutine write_gyre_pulse

! ---------------------------------------------------------------
! GSM, the GYRE Stellar Model HDF5 format (v1.10): root attributes n,
! M_star, R_star, L_star and version, and one double dataset per column,
! the same 18 columns as the GYRE text file from the same point set.
#ifdef YREC_USE_HDF5
subroutine write_gsm_pulse(n, pts, mstar_g, rstar_cm, lstar_cgs, &
     pulse_path, ierr)
      use hdf5
      use h5lt
      use stitched_model, only: n_pulse_cols, ipul_gyre_last
      implicit none

      integer, intent(in) :: n
      double precision, intent(in) :: pts(n_pulse_cols, n)
      double precision, intent(in) :: mstar_g, rstar_cm, lstar_cgs
      character(len=*), intent(in) :: pulse_path
      integer, intent(inout) :: ierr

      ! GSM schema version, stamped in the root 'version' attribute.
      ! 110 = GSM v1.10, the schema whose 18-column layout (kap_kap_T,
      ! eps_eps_rho, ... naming) is written below -- matching MESA's
      ! pulse_gsm.f90 case(110).  (101 is the *text* GYRE format's schema
      ! number and is not a valid GSM version; GYRE aborts on it.)
      integer, parameter :: gyre_schema = 110
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
      call check_io(herr, 'HDF5 could not create the file', 'write', pulse_path, ierr)
      if (ierr /= 0) return

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
      do i = 1, ipul_gyre_last
         call h5ltmake_dataset_double_f(file_id, trim(dset_names(i)), &
              1, dims1, pts(i,1:n), herr)
      end do

      call h5fclose_f(file_id, herr)
      call h5close_f(herr)
      return
end subroutine write_gsm_pulse
#else
subroutine write_gsm_pulse(n, pts, mstar_g, rstar_cm, lstar_cgs, &
     pulse_path, ierr)
      use io_units
      use stitched_model, only: n_pulse_cols, ipul_gyre_last
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: pts(n_pulse_cols, n)
      double precision, intent(in) :: mstar_g, rstar_cm, lstar_cgs
      character(len=*), intent(in) :: pulse_path
      integer, intent(inout) :: ierr

      write(*,*) 'pulse_format = GSM requires an HDF5-enabled build:'
      write(*,*) '  make clean && make USE_HDF5=1'
      write(*,*) '(HDF5_DIR defaults to the MESA SDK; see the Makefile.)'
      write(run_log_unit,*) 'pulse_format = GSM requires USE_HDF5=1 build'
! 2026 ierr-not-stop: an output-format misconfiguration no longer
! kills the run -- the pulse file is skipped (warned above, once per
! attempt); the evolution itself is unaffected.
      return
end subroutine write_gsm_pulse
#endif

! ---------------------------------------------------------------
! Whether this build writes GSM files, so output_init_mesa can reject
! pulse_format = 'GSM' at startup (the stub above is a last guard).
logical function gsm_supported()
#ifdef YREC_USE_HDF5
      gsm_supported = .true.
#else
      gsm_supported = .false.
#endif
end function gsm_supported

! ---------------------------------------------------------------
! GYRE text format (schema 101) over the extended point set: the file
! written with each profile when pulse_format is 'GYRE'.
subroutine write_gyre_ext(n, pts, mstar_g, rstar_cm, lstar_cgs, path, ierr)
      use stitched_model, only: n_pulse_cols, ipul_gyre_last
      integer, intent(in) :: n
      double precision, intent(in) :: pts(n_pulse_cols, n)
      double precision, intent(in) :: mstar_g, rstar_cm, lstar_cgs
      character(len=*), intent(in) :: path
      integer, intent(inout) :: ierr
      integer, parameter :: gyre_schema = 101
      integer :: u, j, i, ios
      character(len=256) :: msg

      open(newunit=u, file=path, status='REPLACE', form='FORMATTED', iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'write', path, ierr)
      if (ierr /= 0) return
      write(u,'(I6,3(1X,1PE26.16),1X,I6)') n, mstar_g, rstar_cm, &
           lstar_cgs, gyre_schema
      do j = 1, n
         write(u,'(I6,99(1X,1PE26.16))') j, (pts(i,j), i = 1, ipul_gyre_last)
      end do
      close(u)
end subroutine write_gyre_ext

end module pulsation_output
