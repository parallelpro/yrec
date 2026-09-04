!----------------------------------------------------------------------
! write_gyre_pulse
!----------------------------------------------------------------------
! New (2026) as part of the YREC readability refactor's pulsation-
! output feature: writes the converged model's structure to a GYRE-
! format stellar model file (MESA/GYRE schema 101, 18 data columns),
! independent of the profile-coupled GSM/FGONG/GYRE writers in
! yrec_output (write_pulse). This routine is triggered from
! yrec_output's output_write_model every pulse_gyre_interval converged
! models (NAMELIST /control/ member, default 0 = off; see
! io/read_controls.f90).
!
! The physics quantities are the same ones YREC already computes for
! every shell in core/henyey_coefficients.f90 (the star%pulse_* and
! star%grada/gradT/adiabatic_index_gamma1 arrays, always populated).
! Column formulas and the two format strings are taken from
! mesa-26.04.1/star/private/pulse_gyre.f90.
!
! YREC's own shell numbering already runs center (1) to surface
! (num_shells), matching GYRE's expected ordering directly -- no
! reversal needed (unlike MESA, which stores its own arrays surface-
! to-center internally and must reverse before writing).
subroutine write_gyre_pulse(num_shells, mass_coordinate, &
     log_density, log_luminosity, log_pressure, log_radius, &
     log_temperature, omega, pulse_path)
      use star_info_lib, only: star, json
      use phys_const_lib
      use math_lib
      implicit none

      integer, intent(in) :: num_shells
      double precision, intent(in) :: mass_coordinate(json), &
           log_density(json), log_luminosity(json), log_pressure(json), &
           log_radius(json), log_temperature(json), omega(json)


      character(len=*), intent(in) :: pulse_path
      integer :: gyre_unit, i
      integer, parameter :: gyre_schema = 101
      double precision :: radius_cm, mass_g, luminosity_erg_s, &
           pressure_cgs, temperature_k, density_cgs, delta, grav, &
           global_data(3)
      double precision :: brunt_n2(json), dr
      double precision :: grav_const_cgs

      open(newunit=gyre_unit,file=pulse_path,status='UNKNOWN',form='FORMATTED')

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
! (see the matching note in core/stitched_model.f90's
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
! and chiT=-chiRho*star%pulse_dlnrho_dlnt (see core/henyey_coefficients.f90).
         delta = -star%pulse_dlnrho_dlnt(i)
! kap_kap_T/eps_eps_T columns follow the GSM/GYRE convention (MESA
! pulse_gyre.f90): the ABSOLUTE derivatives kap*dlnkap/dlnT and
! eps*dlneps/dlnT, not the bare log-derivatives.
         write(gyre_unit,110) i,radius_cm,mass_g,luminosity_erg_s, &
              pressure_cgs,temperature_k,density_cgs,star%gradT(i), &
              brunt_n2(i),star%adiabatic_index_gamma1(i),star%grada(i),delta, &
              star%opacity_zone(i), &
              star%opacity_zone(i)*star%pulse_dlnkap_dlnt(i), &
              star%opacity_zone(i)*star%pulse_dlnkap_dlnrho(i), &
              star%eps_total(i), &
              star%eps_total(i)*star%pulse_dlneps_dlnt(i), &
              star%eps_total(i)*star%pulse_dlneps_dlnrho(i),omega(i)
 110     format(I6,99(1X,1PE26.16))
      end do

      close(gyre_unit)

      return
end subroutine write_gyre_pulse
