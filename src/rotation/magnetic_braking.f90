!----------------------------------------------------------------------
! magnetic_braking
!----------------------------------------------------------------------
! The wind torque on the surface convection zone: the angular momentum a
! magnetized wind carries away, the boundary condition of the angular
! momentum equation. Two laws:
!   Matt / PMM (use_pmm_wind_law): dOmega/dt = C S Omega
!       min(Omega_s, Omega_sat)**(n-1) f_cen / I_cz, where S (the
!       structure factor) scales with M, R, L, P_phot and tau_cz
!       (Matt et al. 2012, 2015), Omega_s is Omega with the optional
!       Rossby / B-field scaling, and f_cen the centrifugal reduction of
!       Matt et al. (2012, ApJ 754, L26), relative to the Sun;
!   Kawaler (1988): dOmega/dt = C (Mdot/1e-14)**exmd Omega (R/Rsun)**exr
!       (M/Msun)**exm min(Omega, Omega_sat)**(n-1) / I_cz, with Mdot fixed
!       at 2e-14 Msun/yr.
! The constant C, the exponents and Omega_sat are controls.
!     brake_convection_zone       remove the angular momentum lost over
!                                 a step from the convection zone's j/m
!                                 (formerly mwind and the Kawaler wind)
!     compute_solid_body_braking  dOmega at the start and end of a step
!                                 for a solid-body star (formerly mcowind
!                                 and the Kawaler spindown)
! A step whose first-guess loss exceeds Omega is split into sub-steps;
! in each, Omega solves Omega = Omega_0 - dOmega(Omega_avg) by averaging
! the explicit and implicit solutions (MHP 12/91).
module magnetic_braking
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use star_info_lib, only: star, json
      use physical_constants, only: ln10, cgl, log10_radius_from_l_teff
      use math_lib
      implicit none
      private
      public :: brake_convection_zone, compute_solid_body_braking

! The law in force for one call, and what its torque needs.
      type :: braking_law
         logical :: matt = .true.
         real(dp) :: omega_saturation = 0.0d0
         real(dp) :: cz_moment_of_inertia = 0.0d0, total_mass_msun = 0.0d0
! Matt: tau_cz interpolated in the step (zero without Rossby scaling),
! the structure factor, the solar centrifugal term and log R
         real(dp) :: turnover_timescale = 0.0d0, structure_factor = 0.0d0
         real(dp) :: fsun = 0.0d0, log10_radius = 0.0d0
! Kawaler: R and the fixed mass-loss rate
         real(dp) :: total_radius_cm = 0.0d0, mass_loss_rate_msun_yr = 0.0d0
      end type braking_law

contains

!----------------------------------------------------------------------
! brake_convection_zone (formerly mwind and the Kawaler wind)
!----------------------------------------------------------------------
! Spins the surface convection zone (zones start_zone..end_zone, masses
! cz_mass_bottom..cz_mass_top, moment of inertia I_cz) down over
! full_timestep from omega_surface, and subtracts the j/m lost from each
! of its zones. The Kawaler law acts only with wind_loss_active; the
! Matt law always does. ierr = 1: the Kawaler saturation needs tau_cz
! and it is undefined.
subroutine brake_convection_zone(log_luminosity_lsun, full_timestep, cz_mass_bottom, &
     cz_mass_top, start_zone, end_zone, wind_loss_active, omega_surface, &
     total_mass_msun, log_teff, cz_moment_of_inertia, &
     specific_angular_momentum, ierr)
      real(dp), intent(in) :: log_luminosity_lsun, full_timestep, &
           cz_mass_bottom, cz_mass_top
      integer, intent(in) :: start_zone, end_zone
      logical, intent(in) :: wind_loss_active
      real(dp), intent(in) :: omega_surface, total_mass_msun, log_teff, &
           cz_moment_of_inertia
      real(dp), intent(inout) :: specific_angular_momentum(json)
      integer, intent(out) :: ierr
      type(braking_law) :: law
      real(dp) :: domega_test, sub_timestep
      real(dp) :: omega_substep_start, omega_iter, omega_iter_new, &
           domega_relative_change, omega_iter_prev
      integer :: num_substeps, substep_idx, iter_count, zone_idx
      real(dp) :: cz_mass, delta_j_per_mass

      ierr = 0
      law%cz_moment_of_inertia = cz_moment_of_inertia
      law%total_mass_msun = total_mass_msun
      if (star%ctrl%use_pmm_wind_law) then
         call set_matt_law(log_luminosity_lsun, total_mass_msun, log_teff, law)
      else
         call set_kawaler_law(.true., law, ierr)
         if (ierr /= 0) return
         if (.not. wind_loss_active) return
         call set_kawaler_radius(log_luminosity_lsun, log_teff, law)
      end if

! The loss rate grows with Omega; a step too long for the first-guess
! loss to leave Omega positive is split into sub-steps.
      domega_test = domega(law, full_timestep, omega_surface, .true.)
      if (domega_test > omega_surface) then
         num_substeps = int(domega_test/omega_surface)+1
         sub_timestep = full_timestep/real(num_substeps, dp)
      else
         num_substeps = 1
         sub_timestep = full_timestep
      endif
      omega_substep_start = omega_surface
      do substep_idx = 1,num_substeps
! first guess: the loss at the starting Omega; then at the average of
! the starting and the previous ending Omega
         iter_count = 0
         omega_iter = omega_substep_start
         omega_iter_prev = omega_substep_start
         omega_fixed_point: do
            iter_count = iter_count + 1
            omega_iter_new = omega_substep_start - domega(law, sub_timestep, omega_iter, .false.)
            domega_relative_change = 2.0d0*abs((omega_iter_prev-omega_iter_new)/ &
                 (omega_iter_prev+omega_iter_new))
            if (domega_relative_change > 1.0d-6) then
               omega_iter = 0.5d0*(omega_substep_start+omega_iter_new)
               omega_iter_prev = omega_iter_new
               if (iter_count <= 20) cycle omega_fixed_point
            endif
            exit omega_fixed_point
         end do omega_fixed_point
         omega_substep_start = omega_iter_new
      end do
! the j/m lost, from every zone of the convection zone
      cz_mass = cz_mass_top - cz_mass_bottom
      delta_j_per_mass = (omega_surface-omega_iter_new)*cz_moment_of_inertia/ &
           cz_mass
      do zone_idx = start_zone,end_zone
         specific_angular_momentum(zone_idx) = &
              specific_angular_momentum(zone_idx) - delta_j_per_mass
      end do
end subroutine brake_convection_zone

!----------------------------------------------------------------------
! compute_solid_body_braking (formerly mcowind and the Kawaler spindown)
!----------------------------------------------------------------------
! For a star rotating as a solid body: the drop in Omega over
! full_timestep from the torque at the start of the step (at omega_old)
! and at its end (at omega_surface). domega_end is averaged with its
! previous value after the first iteration of the caller's solve.
subroutine compute_solid_body_braking(log_luminosity_lsun, full_timestep, cz_moment_of_inertia, &
     iteration_number, omega_surface, total_mass_msun, log_teff, &
     omega_old, domega_start, domega_end, ierr)
      real(dp), intent(in) :: log_luminosity_lsun, full_timestep, &
           cz_moment_of_inertia
      integer, intent(in) :: iteration_number
      real(dp), intent(in) :: omega_surface, total_mass_msun, log_teff
      real(dp), intent(in) :: omega_old
      real(dp), intent(out) :: domega_start
      real(dp), intent(inout) :: domega_end
      integer, intent(out) :: ierr
      type(braking_law) :: law
      real(dp) :: omega_first, omega_now, wind_coefficient, fcen1, fcen2, &
           domega_end_this_iter

      ierr = 0
      law%cz_moment_of_inertia = cz_moment_of_inertia
      omega_first = omega_old
      omega_now = omega_surface
      if (star%ctrl%use_pmm_wind_law) then
         if (star%ctrl%scale_by_rossby_number) then
            law%turnover_timescale = interpolated_turnover_timescale()
            law%omega_saturation = star%job%wind_saturation_omega*star%ctrl%pmm_solar_turnover_timescale/ &
                 law%turnover_timescale
! with B-field scaling, Omega and Omega_sat both carry tau_cz/tau_cz(Sun)
            if (star%ctrl%scale_by_b_field) then
               omega_first = omega_first*law%turnover_timescale/ &
                    star%ctrl%pmm_solar_turnover_timescale
               omega_now = omega_now*law%turnover_timescale/ &
                    star%ctrl%pmm_solar_turnover_timescale
               law%omega_saturation = law%omega_saturation*law%turnover_timescale/ &
                    star%ctrl%pmm_solar_turnover_timescale
            endif
         else
            law%omega_saturation = star%job%wind_saturation_omega
         endif
         law%structure_factor = matt_structure_factor(total_mass_msun, log_luminosity_lsun, log_teff)
         wind_coefficient = (full_timestep/cz_moment_of_inertia)*star%ctrl%constfactor* &
              law%structure_factor
! the solar centrifugal term, without the 1/G of brake_convection_zone's:
! the original divided by exp(ln10*gl) with gl an unset local (zero)
         law%fsun = 0.5*star%ctrl%pmm_solar_omega**2*star%solar_radius_cgs**3/star%solar_mass_cgs
         law%log10_radius = log10_radius_from_l_teff(log_luminosity_lsun, log_teff, &
              star%log10_solar_luminosity)
         fcen1 = matt_centrifugal_factor(omega_old, law%fsun, law%log10_radius, total_mass_msun)
         fcen2 = matt_centrifugal_factor(omega_surface, law%fsun, law%log10_radius, total_mass_msun)
      else
         call set_kawaler_law(.false., law, ierr)
         if (ierr /= 0) return
         call set_kawaler_radius(log_luminosity_lsun, log_teff, law)
         wind_coefficient = full_timestep/cz_moment_of_inertia*star%ctrl%constfactor* &
              pow((law%mass_loss_rate_msun_yr/1.0d-14), star%ctrl%exmd) &
              *pow((law%total_radius_cm/star%solar_radius_cgs), star%ctrl%exr)*pow(total_mass_msun, star%ctrl%exm)
         fcen1 = 1.0d0
         fcen2 = 1.0d0
      endif
      domega_start = wind_coefficient*pow(min(omega_first,law%omega_saturation), &
           star%ctrl%wind_law_omega_exponent-1.0d0)*omega_old*fcen1
      domega_end_this_iter = wind_coefficient*pow(min(omega_now,law%omega_saturation), &
           star%ctrl%wind_law_omega_exponent-1.0d0)*omega_surface*fcen2
      if (iteration_number == 1) then
         domega_end = domega_end_this_iter
      else
         domega_end = 0.5d0*(domega_end+domega_end_this_iter)
      endif
end subroutine compute_solid_body_braking

!----------------------------------------------------------------------
! set_matt_law
!----------------------------------------------------------------------
! The Matt law's saturation (scaled by the Rossby number when asked:
! Omega_sat(Sun) tau_cz(Sun)/tau_cz, unless the B-field scaling moves
! the factor onto Omega instead), structure factor, solar centrifugal
! term and log R.
subroutine set_matt_law(log_luminosity_lsun, total_mass_msun, log_teff, law)
      real(dp), intent(in) :: log_luminosity_lsun, total_mass_msun, log_teff
      type(braking_law), intent(inout) :: law

      law%matt = .true.
      if (star%ctrl%scale_by_rossby_number) then
         law%turnover_timescale = interpolated_turnover_timescale()
         if (star%ctrl%scale_by_b_field) then
            law%omega_saturation = star%job%wind_saturation_omega
         else
            law%omega_saturation = star%job%wind_saturation_omega*star%ctrl%pmm_solar_turnover_timescale/ &
                 law%turnover_timescale
         endif
      else
         law%omega_saturation = star%job%wind_saturation_omega
      endif
      law%structure_factor = matt_structure_factor(total_mass_msun, log_luminosity_lsun, log_teff)
      law%fsun = 0.5*star%ctrl%pmm_solar_omega**2*star%solar_radius_cgs**3/exp(ln10*cgl)/star%solar_mass_cgs
      law%log10_radius = log10_radius_from_l_teff(log_luminosity_lsun, log_teff, &
           star%log10_solar_luminosity)
end subroutine set_matt_law

!----------------------------------------------------------------------
! set_kawaler_law
!----------------------------------------------------------------------
! The Kawaler law's saturation. MHP 3/09: Omega_sat > 1 means
! Omega_sat(Sun) tau_cz(Sun), so Omega_sat = Omega_sat/tau_cz(star);
! for the convection zone (check_transport) only with angular momentum
! transport on. ierr = 1 when tau_cz is needed and undefined (the
! original stopped; the drivers stop on the error).
subroutine set_kawaler_law(check_transport, law, ierr)
      logical, intent(in) :: check_transport
      type(braking_law), intent(inout) :: law
      integer, intent(out) :: ierr

      ierr = 0
      law%matt = .false.
      if (check_transport .and. .not. star%job%instability_transport_active) then
         law%omega_saturation = star%job%wind_saturation_omega
      else if (star%job%wind_saturation_omega > 1.0d0) then
         if (star%convective_turnover_timescale > 1.0d0) then
            law%omega_saturation = star%job%wind_saturation_omega/star%convective_turnover_timescale
         else
            write(*,911)star%job%wind_saturation_omega,star%convective_turnover_timescale
 911        format('ERROR IN WIND - TAUCZ NOT DEFINED ',1P2E12.3,'STOPPED')
            ierr = 1
            return
         endif
      else
         law%omega_saturation = star%job%wind_saturation_omega
      endif
end subroutine set_kawaler_law

! R from L and Teff, and the Kawaler law's fixed mass-loss rate.
subroutine set_kawaler_radius(log_luminosity_lsun, log_teff, law)
      real(dp), intent(in) :: log_luminosity_lsun, log_teff
      type(braking_law), intent(inout) :: law
      law%total_radius_cm = exp(ln10*log10_radius_from_l_teff(log_luminosity_lsun, log_teff, &
           star%log10_solar_luminosity))
      law%mass_loss_rate_msun_yr = 2.0d-14
end subroutine set_kawaler_radius

!----------------------------------------------------------------------
! domega
!----------------------------------------------------------------------
! The drop in Omega over dt from the law's torque at omega. For the
! Matt law, Omega is scaled by tau_cz/tau_cz(Sun) for the saturation
! test with B-field scaling: at the first guess (initial) only with
! Rossby scaling too, in the iteration always (where tau_cz is zero
! without Rossby scaling, as in the original).
real(dp) function domega(law, dt, omega, initial)
      type(braking_law), intent(in) :: law
      real(dp), intent(in) :: dt, omega
      logical, intent(in) :: initial
      real(dp) :: omega_scaled, fcen

      if (law%matt) then
         omega_scaled = omega
         if (star%ctrl%scale_by_b_field .and. (star%ctrl%scale_by_rossby_number .or. .not. initial)) &
              omega_scaled = omega*law%turnover_timescale/star%ctrl%pmm_solar_turnover_timescale
         fcen = matt_centrifugal_factor(omega, law%fsun, law%log10_radius, law%total_mass_msun)
         domega = (dt/law%cz_moment_of_inertia)*star%ctrl%constfactor* &
              law%structure_factor*omega &
              *pow(min(omega_scaled,law%omega_saturation), (star%ctrl%wind_law_omega_exponent-1.0d0))*fcen
      else
         domega = (dt/law%cz_moment_of_inertia)*star%ctrl%constfactor* &
              pow((law%mass_loss_rate_msun_yr/1.0d-14), star%ctrl%exmd) &
              *omega*pow((law%total_radius_cm/star%solar_radius_cgs), star%ctrl%exr)* &
              pow(law%total_mass_msun, star%ctrl%exm) &
              *pow(min(omega,law%omega_saturation), (star%ctrl%wind_law_omega_exponent-1.0d0))
      end if
end function domega

!----------------------------------------------------------------------
! matt_structure_factor
!----------------------------------------------------------------------
! G Somers 3/17: the Matt law's dependence on the star,
! M**exm (R/Rsun)**exr (L/Lsun)**exl (P_phot/P_phot,Sun)**expr
! (tau_cz/tau_cz,Sun)**extau (the last only with Rossby scaling), with
! P_phot and tau_cz interpolated in the step.
real(dp) function matt_structure_factor(total_mass_msun, log_luminosity_lsun, log_teff)
      real(dp), intent(in) :: total_mass_msun, log_luminosity_lsun, log_teff
      real(dp) :: radius_rsun, luminosity_lsun, photospheric_pressure_ratio, &
           turnover_ratio

      radius_rsun = exp(ln10*log10_radius_from_l_teff(log_luminosity_lsun, log_teff, &
           star%log10_solar_luminosity))/star%solar_radius_cgs
      luminosity_lsun = exp10(log_luminosity_lsun)
      photospheric_pressure_ratio = exp10((star%pphot0+star%fracstep*(star%pphot-star%pphot0)))/ &
           (exp10(star%ctrl%pmm_solar_pressure))
      if (star%ctrl%scale_by_rossby_number) then
         turnover_ratio = interpolated_turnover_timescale()/star%ctrl%pmm_solar_turnover_timescale
      else
         turnover_ratio = 1.
      endif
      matt_structure_factor = pow(total_mass_msun, star%ctrl%exm) * pow(radius_rsun, star%ctrl%exr) * &
           pow(luminosity_lsun, star%ctrl%exl) * pow(photospheric_pressure_ratio, star%ctrl%expr) &
           * pow(turnover_ratio, star%ctrl%extau)
end function matt_structure_factor

!----------------------------------------------------------------------
! matt_centrifugal_factor
!----------------------------------------------------------------------
! The centrifugal reduction of Matt et al. (2012), relative to the Sun:
! ((c_2**2 + f_sun)/(c_2**2 + f))**excen with f = omega**2 R**3/(2 G M).
real(dp) function matt_centrifugal_factor(omega, fsun, log10_radius, total_mass_msun)
      real(dp), intent(in) :: omega, fsun, log10_radius, total_mass_msun
      real(dp) :: f
      f = 0.5*omega**2*exp(ln10*(3.0*log10_radius-cgl))/ &
           total_mass_msun/star%solar_mass_cgs
      matt_centrifugal_factor = pow(((star%ctrl%c_2**2+fsun)/(star%ctrl%c_2**2+f)), star%ctrl%excen)
end function matt_centrifugal_factor

! tau_cz interpolated in the step (MHP 8/17)
real(dp) function interpolated_turnover_timescale()
      interpolated_turnover_timescale = star%convective_turnover_timescale_old+ &
           star%fracstep*(star%convective_turnover_timescale-star%convective_turnover_timescale_old)
end function interpolated_turnover_timescale

end module magnetic_braking
