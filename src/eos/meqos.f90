!----------------------------------------------------------------------
! meqos
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original meqos.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Mihalas, Hummer, Dappen equation of state. Implemented by YC,
! March 7, 1990. Given log10(T), log10(P), X, and Z, calls the MHD
! table interpolator (mhdpx, converted separately; out of scope for
! this pass) and converts its output into the same thermodynamic
! derivative set used by the rest of the EOS machinery.
subroutine meqos(log10_temperature, temperature, log10_pressure, &
     pressure, log10_density, density, hydrogen_fraction, metal_fraction, &
     beta, beta_inverse, beta14, ion_fraction, specific_gas_constant, &
     ion_mean_weight_inverse, electron_mean_weight_inverse, &
     electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
     specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
     dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
     specific_heat_cp_dt, specific_heat_cp_dp)

! LATMO,KSAHA NEEDED FOR EQSAHA
      use luout_lib
      use const_lib
      implicit none

      double precision, intent(in) :: log10_temperature, log10_pressure
      double precision, intent(out) :: temperature, pressure, &
           log10_density, density
      double precision, intent(in) :: hydrogen_fraction, metal_fraction
      double precision, intent(out) :: beta, beta_inverse, beta14
      double precision, intent(out) :: ion_fraction(3)
      double precision, intent(out) :: specific_gas_constant, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, &
           electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
           dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp


      integer, parameter :: ivarx = 25
      double precision, parameter :: cnvs = 0.434294481d0
      double precision, parameter :: zero = 0.0d0

! common/comp/: none of these members are used in this file; declared
! only to preserve the storage layout shared with every other file
! that references common/comp/. Names are chosen to match their usage
! in eqstat2.f90, where they are read.
      double precision :: envelope_hydrogen_fraction, &
           envelope_metal_fraction, zenvm, envelope_amu, &
           envelope_species_fractions(12), xnew, znew, stotal, senv
      common/comp/ envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, envelope_amu, envelope_species_fractions, xnew, znew, &
           stotal, senv

! common/ctlim/: not used in this file; declared only to preserve
! layout. saha_log10t_cutoff is the name used where this member is
! actually read (eqstat2.f90).
      double precision :: atime(14), tcut(5), saha_log10t_cutoff, &
           tenv0, tenv1, tenv, tgcut
      common/ctlim/ atime, tcut, saha_log10t_cutoff, tenv0, tenv1, tenv, tgcut



! common/ccout2/: no member is used anywhere in this batch of files;
! all are unrenamed placeholders preserving the storage layout.
      logical :: ldebug, lcorr, lmilne, ltrack, lstpch
      common/ccout2/ ldebug, lcorr, lmilne, ltrack, lstpch

! common/mhdout/: raw output vector from the MHD table interpolator
! (mhdpx), indexed as documented inline below where each element is
! read.
      double precision :: mhd_output(ivarx)
      common/mhdout/ mhd_output

      save

      integer :: ier_flag
      double precision :: mhdpx_r10
      double precision :: chi_rho, chi_t, log10_specific_heat_cp
      double precision :: specific_gas_constant_check
      integer :: ion_idx

      ier_flag = 0
      temperature = 10.0d0**log10_temperature
      pressure = 10.0d0**log10_pressure
      call mhdpx(log10_pressure, log10_temperature, hydrogen_fraction, &
           mhdpx_r10)
      log10_density = mhd_output(1)
      density = 10.0d0**log10_density
      dlnrho_dlnp = 1.0d0/mhd_output(4)
      chi_rho = mhd_output(4)
      chi_t = mhd_output(5)
      dlnrho_dlnt = -chi_t/chi_rho
      log10_specific_heat_cp = mhd_output(9)
      specific_heat_cp = 10.0d0**log10_specific_heat_cp
      specific_heat_cp_dp = dlnrho_dlnp*mhd_output(12)
      specific_heat_cp_dt = mhd_output(13) + mhd_output(12)*dlnrho_dlnt
      adiabatic_gradient = mhd_output(8)
      adiabatic_gradient_dp = dlnrho_dlnp*mhd_output(10)*cnvs/adiabatic_gradient
      adiabatic_gradient_dt = (mhd_output(11) + mhd_output(10)*dlnrho_dlnt)* &
           cnvs/adiabatic_gradient
      dlnrho_dlnp_dt = dlnrho_dlnt*(adiabatic_gradient_dp - 1.0d0 + &
           dlnrho_dlnp + specific_heat_cp_dp)
      dlnrho_dlnt_dt = dlnrho_dlnt*(adiabatic_gradient_dt + 1.0d0 + &
           dlnrho_dlnt + specific_heat_cp_dt)
      beta = 10.0d0**(mhd_output(20) - mhd_output(2))
      beta_inverse = 1.0d0/beta
      beta14 = 1.0d0 - beta
      do ion_idx = 1, 3
      ion_fraction(ion_idx) = mhd_output(ion_idx+13)
      end do
      electron_degeneracy_parameter = mhd_output(18)
      specific_gas_constant = 10.0d0**mhd_output(20)/(density*temperature)
      call mu(temperature, pressure, density, hydrogen_fraction, &
           metal_fraction, specific_gas_constant_check, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, beta)
      if (abs((specific_gas_constant-specific_gas_constant_check)/ &
           specific_gas_constant).gt.5.0d-7) then
          write(iowr,*)' ERROR(MHD) IN MEAN WEIGHTS ... '
          write(iowr,*) specific_gas_constant, &
               specific_gas_constant_check
        write(iowr,*) 'ERROR (MHD): CHECK MU'
          write(short_file_unit,*)' ERROR(MHD): IN MEAN WEIGHTS ... '
          write(short_file_unit,*) specific_gas_constant, &
               specific_gas_constant_check
        write(short_file_unit,*) 'ERROR (MHD): CHECK MU'
          stop
      end if
      return
!   999 CONTINUE
      write(short_file_unit,*) 'ERROR(MHD):... MHD TABLE FAIL'
      stop
end subroutine meqos
