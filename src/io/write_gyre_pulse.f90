!----------------------------------------------------------------------
! write_gyre_pulse
!----------------------------------------------------------------------
! New (2026) as part of the YREC readability refactor's pulsation-
! output feature: writes the converged model's structure to a GYRE-
! format stellar model file (MESA/GYRE schema 101, 18 data columns),
! independent of YREC's own older path-length-triggered OPAL-format
! pulsation writer (misc/pdist.f90 + io/wrtmod.f90, gated by LPULSE).
! This routine is instead triggered by wrtout.f90 every
! pulse_gyre_interval converged models (common/pulsegyre/,
! new NAMELIST /control/ member, default 0 = off; see core/parmin.f90).
!
! The physics quantities are the same ones YREC already computes for
! every shell in misc/coefft.f90 (common/scrtch/, common/pulse1/,
! common/sound/) -- see the "MHP 8/25 unconditional" note there for
! why those arrays are now always populated instead of only when the
! older LPULSE mechanism is active. Column formulas and the two format
! strings are taken from mesa-26.04.1/star/private/pulse_gyre.f90.
!
! YREC's own shell numbering already runs center (1) to surface
! (num_shells), matching GYRE's expected ordering directly -- no
! reversal needed (unlike MESA, which stores its own arrays surface-
! to-center internally and must reverse before writing).
subroutine write_gyre_pulse(num_shells, model_number, mass_coordinate, &
     log_density, log_luminosity, log_pressure, log_radius, &
     log_temperature, omega)
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: num_shells, model_number
      double precision, intent(in) :: mass_coordinate(json), &
           log_density(json), log_luminosity(json), log_pressure(json), &
           log_radius(json), log_temperature(json), omega(json)

! common/scrtch/: so (opacity) and del_grad (2=actual gradient,
! 3=adiabatic gradient) are used here. Naming matches wrtout.f90/
! coefft.f90.
      double precision :: sesum(json), seg(7,json), sbeta(json), seta(json)
      logical :: locons(json)
      double precision :: so(json), del_grad(3,json), sfxion(3,json), &
           svel(json), scp(json)
      common/scrtch/ sesum, seg, sbeta, seta, locons, so, del_grad, &
           sfxion, svel, scp
! common/sound/: adiabatic_index_gamma1 (Gamma_1), used here. Naming
! matches wrtout.f90.
      double precision :: adiabatic_index_gamma1(json)
      logical :: sound_speed_output_active
      common/sound/ adiabatic_index_gamma1, sound_speed_output_active
! common/pulse1/: the kappa/epsilon logarithmic derivatives, and
! pulse_dlnrho_dlnt/pulse_dlnrho_dlnp (used below to form
! delta=chiT/chiRho), are used here. Naming matches coefft.f90, which
! now computes these for every shell unconditionally.
      double precision :: pulse_dlnrho_dlnp(json), pulse_dlneps_dlnrho(json), &
           pulse_dlneps_dlnt(json), pulse_dlnkap_dlnrho(json), &
           pulse_dlnkap_dlnt(json), pulse_specific_heat(json), &
           pulse_mean_molecular_weight(json), pulse_dlnrho_dlnt(json), &
           pulse_electron_mean_molecular_weight(json)
      logical :: lpumod
      common/pulse1/ pulse_dlnrho_dlnp, pulse_dlneps_dlnrho, pulse_dlneps_dlnt, &
           pulse_dlnkap_dlnrho, pulse_dlnkap_dlnt, pulse_specific_heat, &
           pulse_mean_molecular_weight, pulse_dlnrho_dlnt, &
           pulse_electron_mean_molecular_weight, lpumod
! common/const/: solar_luminosity_cgs, used to convert the per-shell
! linear L/Lsun array (log_luminosity here, despite its name -- see
! core/main.f90's cross-callee naming note on HL/luminosity_lsun) to
! erg/s. Naming matches wrtout.f90.
      double precision :: solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude
      common/const/ solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude

      integer :: gyre_unit, i
      integer, parameter :: gyre_schema = 101
      character(len=5) :: model_suffix
      character(len=64) :: gyre_path
      double precision :: radius_cm, mass_g, luminosity_erg_s, &
           pressure_cgs, temperature_k, density_cgs, delta, grav, &
           brunt_n2, global_data(3)

      write(model_suffix,'(I5.5)') model_number
      gyre_path = 'gyre_profile_'//model_suffix//'.data.GYRE'
      open(newunit=gyre_unit,file=gyre_path,status='UNKNOWN',form='FORMATTED')

      global_data(1) = mass_coordinate(num_shells)
      global_data(2) = exp(ln10*log_radius(num_shells))
      global_data(3) = log_luminosity(num_shells)*solar_luminosity_cgs
      write(gyre_unit,100) num_shells,global_data,gyre_schema
 100  format(I6,3(1X,1PE26.16),1X,I6)

      do i = 1,num_shells
         radius_cm = exp(ln10*log_radius(i))
         mass_g = mass_coordinate(i)
         luminosity_erg_s = log_luminosity(i)*solar_luminosity_cgs
         pressure_cgs = exp(ln10*log_pressure(i))
         temperature_k = exp(ln10*log_temperature(i))
         density_cgs = exp(ln10*log_density(i))
! delta = chiT/chiRho = -pulse_dlnrho_dlnt, since chiRho=1/pulse_dlnrho_dlnp
! and chiT=-chiRho*pulse_dlnrho_dlnt (see misc/coefft.f90 around line 639).
         delta = -pulse_dlnrho_dlnt(i)
         if (radius_cm.gt.0.0d0) then
            grav = exp(ln10*cgl)*mass_g/(radius_cm*radius_cm)
            brunt_n2 = grav*grav*(density_cgs/pressure_cgs)*delta* &
                 (del_grad(3,i)-del_grad(2,i))
         else
            brunt_n2 = 0.0d0
         end if
         write(gyre_unit,110) i,radius_cm,mass_g,luminosity_erg_s, &
              pressure_cgs,temperature_k,density_cgs,del_grad(2,i), &
              brunt_n2,adiabatic_index_gamma1(i),del_grad(3,i),delta, &
              so(i),pulse_dlnkap_dlnt(i),pulse_dlnkap_dlnrho(i),sesum(i), &
              pulse_dlneps_dlnt(i),pulse_dlneps_dlnrho(i),omega(i)
 110     format(I6,99(1X,1PE26.16))
      end do

      close(gyre_unit)

      return
end subroutine write_gyre_pulse
