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
! The constant C, the exponents and Omega_sat follow from the controls
! (derive_wind_law).
!     derive_wind_law             the law's constant, exponents and
!                                 saturation from its controls
!     brake_convection_zone       remove the angular momentum lost over
!                                 a step from the convection zone's j/m
!                                 (formerly mwind and the Kawaler wind)
!     compute_solid_body_braking  dOmega at the start and end of a step
!                                 for a solid-body star (formerly mcowind
!                                 and the Kawaler spindown)
!     outer_point_wind_terms      the wind's terms for the outermost point
!                                 of a transport region's grid, from
!                                 compute_solid_body_braking
! A step whose first-guess loss exceeds Omega is split into sub-steps;
! in each, Omega solves Omega = Omega_0 - dOmega(Omega_avg) by averaging
! the explicit and implicit solutions (MHP 12/91).
! The laws read the controls and the solar units from rotation_settings
! and the step's state from a wind_conditions record; never the star.
module magnetic_braking
      use rotation_configuration, only: rotation_settings
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: json
      use physical_constants, only: ln10, log10_radius_from_l_teff
      use math_lib
      use rotation_messages, only: rotation_report, message_length, unwritten_line
      implicit none
      private
      public :: brake_convection_zone, compute_solid_body_braking, outer_point_wind_terms
      public :: derive_wind_law

! The wind law's controls (setup's rotation_settings_from_star fills
! them): the Matt law or the Kawaler law; the Matt law's name (K97, V13
! or CUS); the Kawaler index alfa (wind_law_index) and the constant
! factor fk of both laws (wind_law_constant); the Matt exponents a, b,
! c, d, m as read (CUS keeps them); the solar torque (erg) and the solar
! omega (rad/s) that normalize the Matt law; the saturation omega
! (rad/s) and the same in units of the solar omega (solar_omega_sat,
! which wins below solar_units_saturation_limit); the field scaling as
! read.
      type, public :: wind_law_inputs
         logical :: use_pmm_wind_law = .false.
         character(len=3) :: law_name = '   '
         real(dp) :: law_index = 0.0d0, law_constant = 0.0d0
         real(dp) :: exponent_a = 0.0d0, exponent_b = 0.0d0, exponent_c = 0.0d0, &
              exponent_d = 0.0d0, exponent_m = 0.0d0
         real(dp) :: solar_torque = 0.0d0, solar_omega = 0.0d0
         real(dp) :: saturation_omega = 0.0d0, saturation_omega_solar = 0.0d0
         logical :: scale_by_b_field = .false.
      end type wind_law_inputs

! The law in force: the constant C and the exponents of omega, the
! mass-loss rate, R, M, L, the photospheric pressure, the turnover time
! and the centrifugal term, the centrifugal constant, the saturation
! omega (rad/s) and the field scaling (rotation_settings' wind fields).
      type, public :: wind_law_coefficients
         real(dp) :: constant = 0.0d0, omega_exponent = 0.0d0, &
              mass_loss_exponent = 0.0d0, radius_exponent = 0.0d0, &
              mass_exponent = 0.0d0, luminosity_exponent = 0.0d0, &
              pressure_exponent = 0.0d0, turnover_exponent = 0.0d0, &
              centrifugal_exponent = 0.0d0, centrifugal_constant = 0.0d0, &
              saturation_omega = 0.0d0
         logical :: scale_by_b_field = .false.
      end type wind_law_coefficients

! The Kawaler (1988) law's normalization: Jdot = fk kawaler_torque_scale
! kawaler_omega_scale**alfa ... (cgs). The Matt law's centrifugal
! constant c_2 (Matt et al. 2012, ApJ 754, L26); V13's torque exponent m
! (van Saders & Pinsonneault 2013). solar_omega_sat below
! solar_units_saturation_limit gives the saturation omega in units of
! the solar omega.
      real(dp), parameter :: kawaler_torque_scale = 2.036d33, kawaler_omega_scale = 1.452d9
      real(dp), parameter :: matt_centrifugal_constant = 0.0506d0, v13_exponent_m = 0.22d0
      real(dp), parameter :: solar_units_saturation_limit = 1.0d3

! The step's state the laws read (the caller fills it from the star):
! the convective turnover time tau_cz (s) at the start and the end of the
! step, log10 of the photospheric pressure (dyn/cm^2) at the start and
! the end, and the fraction of the step elapsed at the current rotation
! sub-step, at which the Matt law interpolates both.
      type, public :: wind_conditions
         real(dp) :: turnover_timescale_start = 0.0d0, turnover_timescale_end = 0.0d0
         real(dp) :: log_pphot_start = 0.0d0, log_pphot_end = 0.0d0
         real(dp) :: step_fraction = 0.0d0
      end type wind_conditions

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
! the step's state (the caller's)
         type(wind_conditions) :: conditions
      end type braking_law

! brake_convection_zone's implicit midpoint iteration stops when two
! iterates of omega differ by at most brake_iteration_tolerance,
! relatively, or after brake_iteration_limit corrections. The Kawaler
! law's mass-loss rate is fixed at kawaler_mass_loss_rate_msun_yr and
! enters as (Mdot/kawaler_mass_loss_unit_msun_yr)**exmd.
      double precision, parameter :: brake_iteration_tolerance = 1.0d-6
      integer, parameter :: brake_iteration_limit = 20
      double precision, parameter :: kawaler_mass_loss_rate_msun_yr = 2.0d-14, &
           kawaler_mass_loss_unit_msun_yr = 1.0d-14
contains

!----------------------------------------------------------------------
! derive_wind_law
!----------------------------------------------------------------------
! The wind law's coefficients from its controls (formerly computed by
! the controls read, derive_rotation_controls).
! The Kawaler (1988) law, for an index alfa:
!   Jdot = fk 2.036e33 1.452e9**alfa (Mdot/1e-14)**(1-2alfa/3)
!          Omega**(1+4alfa/3) (R/Rsun)**(2-alfa) (M/Msun)**(-alfa/3),
! so the exponents of Mdot, Omega, R and M; fk = 1 reproduces the
! observed solar angular momentum loss.
! The Matt (PMM) law: K97 is Krishnamurthi et al. (1997) (m = alfa/3,
! a = 2, b = 1, c = 2, d = 0, no field scaling), V13 van Saders &
! Pinsonneault (2013, after Matt et al. 2008, 2012) (m = 0.22, a = 2,
! b = 1, c = 0, d = 1, field scaling), and CUS keeps the exponents
! given. With every term scaled to solar,
!   Jdot = w Bmag**4m Mdot**(1-2m) R**(5m+2) M**(-m) Fcen**md
!   Bmag = Pphot**0.5 (w R**c)**b Tcz**d
!   Mdot = R M**-1 Lbol w**a  (from GM/R Mdot = Lx, Lx/Lbol ~ w**a)
! so Jdot ~ w**(1+a-2ma+4mb) Tcz**(4md) Lbol**(1-2m) Pphot**(2m) Fcen**md,
! and, ignoring the GM/R term in Mdot, R**(2+5m+4mbc) M**(-m). The
! centrifugal term {K2**2/(K2**2 + 0.5 w**2 R**3/GM)} has the exponent
! m d; d switches the Tcz and centrifugal terms off. The Tcz term is
! computed in the braking itself, so its exponent here is 0. The
! constant makes the torque at the solar omega fk times the solar
! torque.
pure function derive_wind_law(inputs) result(law)
      type(wind_law_inputs), intent(in) :: inputs
      type(wind_law_coefficients) :: law
      real(dp), parameter :: one_third = 1.0d0/3.0d0, two_thirds = 2.0d0/3.0d0
      real(dp) :: a, b, c, d, m, alfa

      alfa = inputs%law_index
      law%saturation_omega = inputs%saturation_omega
      law%scale_by_b_field = inputs%scale_by_b_field
      if (.not. inputs%use_pmm_wind_law) then
         law%constant = inputs%law_constant*kawaler_torque_scale*pow(kawaler_omega_scale, alfa)
         law%mass_loss_exponent = 1.0d0 - two_thirds*alfa
         law%omega_exponent = 1.0d0 + 2.0d0*two_thirds*alfa
         law%radius_exponent = 2.0d0 - alfa
         law%mass_exponent = -one_third*alfa
         return
      end if
      a = inputs%exponent_a
      b = inputs%exponent_b
      c = inputs%exponent_c
      d = inputs%exponent_d
      m = inputs%exponent_m
      select case (inputs%law_name)
      case ('K97')
         m = alfa/3.0d0
         a = 2.0d0
         b = 1.0d0
         c = 2.0d0
         d = 0.0d0
         law%scale_by_b_field = .false.
      case ('V13')
         m = v13_exponent_m
         a = 2.0d0
         b = 1.0d0
         c = 0.0d0
         d = 1.0d0
         law%scale_by_b_field = .true.
      end select
      law%omega_exponent = 1.0d0 + a - 2.0d0*a*m + 4.0d0*m*b
      law%turnover_exponent = 0.0d0
      law%radius_exponent = 2.0d0 + 5.0d0*m - 4.0d0*m*c
      law%mass_exponent = -m
      law%luminosity_exponent = 1.0d0 - 2.0d0*m
      law%pressure_exponent = 2.0d0*m*d
      law%centrifugal_exponent = m*d
      law%centrifugal_constant = matt_centrifugal_constant
      law%constant = inputs%law_constant*inputs%solar_torque/pow(inputs%solar_omega, law%omega_exponent)
      if (inputs%saturation_omega_solar < solar_units_saturation_limit) &
           law%saturation_omega = inputs%saturation_omega_solar*inputs%solar_omega
end function derive_wind_law

!----------------------------------------------------------------------
! brake_convection_zone (formerly mwind and the Kawaler wind)
!----------------------------------------------------------------------
! Spins the surface convection zone (zones start_zone..end_zone, masses
! cz_mass_bottom..cz_mass_top, moment of inertia I_cz) down over
! full_timestep from omega_surface, and subtracts the j/m lost from each
! of its zones. The Kawaler law acts only with wind_loss_active; the
! Matt law always does. conditions: the step's tau_cz, P_phot and elapsed
! fraction. ierr = 1: the Kawaler saturation needs tau_cz and it is
! undefined.
subroutine brake_convection_zone(log_luminosity_lsun, full_timestep, cz_mass_bottom, &
     cz_mass_top, start_zone, end_zone, wind_loss_active, omega_surface, &
     total_mass_msun, log_teff, cz_moment_of_inertia, &
     specific_angular_momentum, conditions, ierr)
      real(dp), intent(in) :: log_luminosity_lsun, full_timestep, &
           cz_mass_bottom, cz_mass_top
      integer, intent(in) :: start_zone, end_zone
      logical, intent(in) :: wind_loss_active
      real(dp), intent(in) :: omega_surface, total_mass_msun, log_teff, &
           cz_moment_of_inertia
      real(dp), intent(inout) :: specific_angular_momentum(json)
      type(wind_conditions), intent(in) :: conditions
      integer, intent(out) :: ierr
      type(braking_law) :: law
      real(dp) :: domega_test, sub_timestep
      real(dp) :: omega_substep_start, omega_iter, omega_iter_new, &
           domega_relative_change, omega_iter_prev
      integer :: num_substeps, substep_idx, iter_count, zone_idx
      real(dp) :: cz_mass, delta_j_per_mass

      ierr = 0
      law%conditions = conditions
      law%cz_moment_of_inertia = cz_moment_of_inertia
      law%total_mass_msun = total_mass_msun
! Neither law removes angular momentum with the wind off (the Matt law
! once did; docs/bug_sweep.md O17).
      if (rotation_settings%use_pmm_wind_law) then
         if (.not. wind_loss_active) return
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
            if (domega_relative_change > brake_iteration_tolerance) then
               omega_iter = 0.5d0*(omega_substep_start+omega_iter_new)
               omega_iter_prev = omega_iter_new
               if (iter_count <= brake_iteration_limit) cycle omega_fixed_point
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
! conditions: the step's tau_cz, P_phot and elapsed fraction.
subroutine compute_solid_body_braking(log_luminosity_lsun, full_timestep, cz_moment_of_inertia, &
     iteration_number, omega_surface, total_mass_msun, log_teff, &
     omega_old, domega_start, domega_end, conditions, ierr)
      real(dp), intent(in) :: log_luminosity_lsun, full_timestep, &
           cz_moment_of_inertia
      integer, intent(in) :: iteration_number
      real(dp), intent(in) :: omega_surface, total_mass_msun, log_teff
      real(dp), intent(in) :: omega_old
      real(dp), intent(out) :: domega_start
      real(dp), intent(inout) :: domega_end
      type(wind_conditions), intent(in) :: conditions
      integer, intent(out) :: ierr
      type(braking_law) :: law
      real(dp) :: omega_first, omega_now, wind_coefficient, fcen1, fcen2, &
           domega_end_this_iter

      ierr = 0
      law%conditions = conditions
      law%cz_moment_of_inertia = cz_moment_of_inertia
      omega_first = omega_old
      omega_now = omega_surface
      if (rotation_settings%use_pmm_wind_law) then
         if (rotation_settings%scale_by_rossby_number .or. rotation_settings%scale_by_b_field) &
              law%turnover_timescale = interpolated_turnover_timescale(law%conditions)
         if (rotation_settings%scale_by_rossby_number) then
            law%omega_saturation = rotation_settings%wind_saturation_omega*rotation_settings%pmm_solar_turnover_timescale/ &
                 law%turnover_timescale
! with B-field scaling, Omega and Omega_sat both carry tau_cz/tau_cz(Sun)
            if (rotation_settings%scale_by_b_field) then
               omega_first = omega_first*law%turnover_timescale/ &
                    rotation_settings%pmm_solar_turnover_timescale
               omega_now = omega_now*law%turnover_timescale/ &
                    rotation_settings%pmm_solar_turnover_timescale
               law%omega_saturation = law%omega_saturation*law%turnover_timescale/ &
                    rotation_settings%pmm_solar_turnover_timescale
            endif
         else
! without Rossby scaling the saturation is fixed; B-field scaling
! still scales Omega (as brake_convection_zone does; docs/bug_sweep.md
! O17)
            law%omega_saturation = rotation_settings%wind_saturation_omega
            if (rotation_settings%scale_by_b_field) then
               omega_first = omega_first*law%turnover_timescale/ &
                    rotation_settings%pmm_solar_turnover_timescale
               omega_now = omega_now*law%turnover_timescale/ &
                    rotation_settings%pmm_solar_turnover_timescale
            endif
         endif
         law%structure_factor = matt_structure_factor(total_mass_msun, log_luminosity_lsun, log_teff, &
              law%conditions)
         wind_coefficient = (full_timestep/cz_moment_of_inertia)*rotation_settings%wind_constant* &
              law%structure_factor
! the solar centrifugal term, (1/2) Omega_sun^2 R_sun^3/(G M_sun), as
! brake_convection_zone's (the original divided by exp(ln10*gl) with gl
! an unset local, zero, so it lacked the 1/G; docs/bug_sweep.md O17)
         law%fsun = 0.5*rotation_settings%pmm_solar_omega**2*rotation_settings%solar_radius_cgs**3 &
              /exp(ln10*rotation_settings%log10_gravitational_constant)/rotation_settings%solar_mass_cgs
         law%log10_radius = log10_radius_from_l_teff(log_luminosity_lsun, log_teff, &
              rotation_settings%log10_solar_luminosity)
         fcen1 = matt_centrifugal_factor(omega_old, law%fsun, law%log10_radius, total_mass_msun)
         fcen2 = matt_centrifugal_factor(omega_surface, law%fsun, law%log10_radius, total_mass_msun)
      else
         call set_kawaler_law(.false., law, ierr)
         if (ierr /= 0) return
         call set_kawaler_radius(log_luminosity_lsun, log_teff, law)
         wind_coefficient = full_timestep/cz_moment_of_inertia*rotation_settings%wind_constant* &
              pow((law%mass_loss_rate_msun_yr/kawaler_mass_loss_unit_msun_yr), rotation_settings%wind_mass_loss_exponent) &
              *pow((law%total_radius_cm/rotation_settings%solar_radius_cgs), rotation_settings%wind_radius_exponent)*pow(total_mass_msun, rotation_settings%wind_mass_exponent)
         fcen1 = 1.0d0
         fcen2 = 1.0d0
      endif
      domega_start = wind_coefficient*pow(min(omega_first,law%omega_saturation), &
           rotation_settings%wind_law_omega_exponent-1.0d0)*omega_old*fcen1
      domega_end_this_iter = wind_coefficient*pow(min(omega_now,law%omega_saturation), &
           rotation_settings%wind_law_omega_exponent-1.0d0)*omega_surface*fcen2
      if (iteration_number == 1) then
         domega_end = domega_end_this_iter
      else
         domega_end = 0.5d0*(domega_end+domega_end_this_iter)
      endif
end subroutine compute_solid_body_braking

! ---------------------------------------------------------------
! outer_point_wind_terms
! ---------------------------------------------------------------
! The wind's terms for the outermost point of a transport region's grid,
! the surface boundary condition of am_diffusion_coeffs. They act only
! when the region reaches the surface convection zone
! (reaches_surface_cz). Under disk locking (disk_locked) the point keeps
! its omega (fix_omega_at_surface) and loses nothing. Otherwise, with
! use_wind_torque, wind_loss_explicit and wind_loss_implicit are the
! drops of the surface omega over the sub-step from the torque at its
! start (at omega_start) and at its current end (at omega_now), for the
! point's moment of inertia outer_point_inertia
! (compute_solid_body_braking); the implicit one is averaged with the
! surface region's previous one, surface_implicit, when an earlier
! iteration set it (surface_implicit_set), and becomes the new
! surface_implicit. cz_moment_of_inertia and omega_surface are then set
! to the inertia and the omega the torque used. Elsewhere both terms are
! zero and surface_implicit is left alone. conditions: the step's
! tau_cz, P_phot and elapsed fraction. ierr from
! compute_solid_body_braking.
subroutine outer_point_wind_terms(reaches_surface_cz, disk_locked, iteration, &
     log_luminosity_lsun, sub_timestep, outer_point_inertia, omega_now, &
     total_mass_msun, log_teff, omega_start, conditions, &
     surface_implicit, surface_implicit_set, &
     wind_loss_explicit, wind_loss_implicit, fix_omega_at_surface, &
     cz_moment_of_inertia, omega_surface, ierr)
      logical, intent(in) :: reaches_surface_cz, disk_locked
      integer, intent(in) :: iteration
      real(dp), intent(in) :: log_luminosity_lsun, sub_timestep, outer_point_inertia, &
           omega_now, total_mass_msun, log_teff, omega_start
      type(wind_conditions), intent(in) :: conditions
      real(dp), intent(inout) :: surface_implicit
      logical, intent(inout) :: surface_implicit_set
      real(dp), intent(out) :: wind_loss_explicit, wind_loss_implicit
      logical, intent(out) :: fix_omega_at_surface
      real(dp), intent(inout) :: cz_moment_of_inertia, omega_surface
      integer, intent(out) :: ierr

      ierr = 0
      fix_omega_at_surface = .false.
      if (reaches_surface_cz .and. disk_locked) then
         wind_loss_explicit = 0.0d0
         wind_loss_implicit = 0.0d0
         fix_omega_at_surface = .true.
      else if (reaches_surface_cz .and. rotation_settings%use_wind_torque) then
         cz_moment_of_inertia = outer_point_inertia
         omega_surface = omega_now
! compute_solid_body_braking averages with its domega_end on entry
! from the second iteration on: start afresh when no earlier iteration
! gave the surface region a term.
         wind_loss_implicit = surface_implicit
         call compute_solid_body_braking(log_luminosity_lsun, sub_timestep, &
              cz_moment_of_inertia, merge(iteration, 1, surface_implicit_set), &
              omega_surface, total_mass_msun, log_teff, &
              omega_start, wind_loss_explicit, wind_loss_implicit, conditions, ierr)
         surface_implicit = wind_loss_implicit
         surface_implicit_set = .true.
      else
         wind_loss_explicit = 0.0d0
         wind_loss_implicit = 0.0d0
      end if
end subroutine outer_point_wind_terms

! The Matt law's saturation (scaled by the Rossby number when asked:
! Omega_sat(Sun) tau_cz(Sun)/tau_cz, unless the B-field scaling moves
! the factor onto Omega instead), structure factor, solar centrifugal
! term and log R.
subroutine set_matt_law(log_luminosity_lsun, total_mass_msun, log_teff, law)
      real(dp), intent(in) :: log_luminosity_lsun, total_mass_msun, log_teff
      type(braking_law), intent(inout) :: law

      law%matt = .true.
! B-field scaling needs tau_cz too, with or without Rossby scaling (it
! was left zero without, which zeroed the torque; docs/bug_sweep.md O17)
      if (rotation_settings%scale_by_rossby_number .or. rotation_settings%scale_by_b_field) &
           law%turnover_timescale = interpolated_turnover_timescale(law%conditions)
      if (rotation_settings%scale_by_rossby_number) then
         if (rotation_settings%scale_by_b_field) then
            law%omega_saturation = rotation_settings%wind_saturation_omega
         else
            law%omega_saturation = rotation_settings%wind_saturation_omega*rotation_settings%pmm_solar_turnover_timescale/ &
                 law%turnover_timescale
         endif
      else
         law%omega_saturation = rotation_settings%wind_saturation_omega
      endif
      law%structure_factor = matt_structure_factor(total_mass_msun, log_luminosity_lsun, log_teff, &
              law%conditions)
      law%fsun = 0.5*rotation_settings%pmm_solar_omega**2*rotation_settings%solar_radius_cgs**3 &
           /exp(ln10*rotation_settings%log10_gravitational_constant)/rotation_settings%solar_mass_cgs
      law%log10_radius = log10_radius_from_l_teff(log_luminosity_lsun, log_teff, &
           rotation_settings%log10_solar_luminosity)
end subroutine set_matt_law

! The Kawaler law's saturation. MHP 3/09: Omega_sat > 1 means
! Omega_sat(Sun) tau_cz(Sun), so Omega_sat = Omega_sat/tau_cz(star);
! for the convection zone (check_transport) only with angular momentum
! transport on. ierr = 1 when tau_cz is needed and undefined (the
! original stopped; the drivers stop on the error).
subroutine set_kawaler_law(check_transport, law, ierr)
      logical, intent(in) :: check_transport
      type(braking_law), intent(inout) :: law
      integer, intent(out) :: ierr

      character(len=message_length) :: lines(1)
      ierr = 0
      law%matt = .false.
      if (check_transport .and. .not. rotation_settings%instability_transport_active) then
         law%omega_saturation = rotation_settings%wind_saturation_omega
      else if (rotation_settings%wind_saturation_omega > 1.0d0) then
         if (law%conditions%turnover_timescale_end > 1.0d0) then
            law%omega_saturation = rotation_settings%wind_saturation_omega/law%conditions%turnover_timescale_end
         else
            lines = unwritten_line
            write(lines,911)rotation_settings%wind_saturation_omega,law%conditions%turnover_timescale_end
            call rotation_report(lines)
 911        format('ERROR IN WIND - TAUCZ NOT DEFINED ',1P2E12.3,'STOPPED')
            ierr = 1
            return
         endif
      else
         law%omega_saturation = rotation_settings%wind_saturation_omega
      endif
end subroutine set_kawaler_law

! R from L and Teff, and the Kawaler law's fixed mass-loss rate.
subroutine set_kawaler_radius(log_luminosity_lsun, log_teff, law)
      real(dp), intent(in) :: log_luminosity_lsun, log_teff
      type(braking_law), intent(inout) :: law
      law%total_radius_cm = exp(ln10*log10_radius_from_l_teff(log_luminosity_lsun, log_teff, &
           rotation_settings%log10_solar_luminosity))
      law%mass_loss_rate_msun_yr = kawaler_mass_loss_rate_msun_yr
end subroutine set_kawaler_radius

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
         if (rotation_settings%scale_by_b_field .and. (rotation_settings%scale_by_rossby_number .or. .not. initial)) &
              omega_scaled = omega*law%turnover_timescale/rotation_settings%pmm_solar_turnover_timescale
         fcen = matt_centrifugal_factor(omega, law%fsun, law%log10_radius, law%total_mass_msun)
         domega = (dt/law%cz_moment_of_inertia)*rotation_settings%wind_constant* &
              law%structure_factor*omega &
              *pow(min(omega_scaled,law%omega_saturation), (rotation_settings%wind_law_omega_exponent-1.0d0))*fcen
      else
         domega = (dt/law%cz_moment_of_inertia)*rotation_settings%wind_constant* &
              pow((law%mass_loss_rate_msun_yr/kawaler_mass_loss_unit_msun_yr), rotation_settings%wind_mass_loss_exponent) &
              *omega*pow((law%total_radius_cm/rotation_settings%solar_radius_cgs), rotation_settings%wind_radius_exponent)* &
              pow(law%total_mass_msun, rotation_settings%wind_mass_exponent) &
              *pow(min(omega,law%omega_saturation), (rotation_settings%wind_law_omega_exponent-1.0d0))
      end if
end function domega

! G Somers 3/17: the Matt law's dependence on the star,
! M**exm (R/Rsun)**exr (L/Lsun)**exl (P_phot/P_phot,Sun)**expr
! (tau_cz/tau_cz,Sun)**extau (the last only with Rossby scaling), with
! P_phot and tau_cz interpolated in the step.
real(dp) function matt_structure_factor(total_mass_msun, log_luminosity_lsun, log_teff, &
     conditions)
      real(dp), intent(in) :: total_mass_msun, log_luminosity_lsun, log_teff
      type(wind_conditions), intent(in) :: conditions
      real(dp) :: radius_rsun, luminosity_lsun, photospheric_pressure_ratio, &
           turnover_ratio

      radius_rsun = exp(ln10*log10_radius_from_l_teff(log_luminosity_lsun, log_teff, &
           rotation_settings%log10_solar_luminosity))/rotation_settings%solar_radius_cgs
      luminosity_lsun = exp10(log_luminosity_lsun)
      photospheric_pressure_ratio = exp10((conditions%log_pphot_start+conditions%step_fraction* &
           (conditions%log_pphot_end-conditions%log_pphot_start)))/ &
           (exp10(rotation_settings%pmm_solar_pressure))
      if (rotation_settings%scale_by_rossby_number) then
         turnover_ratio = interpolated_turnover_timescale(conditions)/rotation_settings%pmm_solar_turnover_timescale
      else
         turnover_ratio = 1.
      endif
      matt_structure_factor = pow(total_mass_msun, rotation_settings%wind_mass_exponent) * pow(radius_rsun, rotation_settings%wind_radius_exponent) * &
           pow(luminosity_lsun, rotation_settings%wind_luminosity_exponent) * pow(photospheric_pressure_ratio, rotation_settings%wind_pressure_exponent) &
           * pow(turnover_ratio, rotation_settings%wind_turnover_exponent)
end function matt_structure_factor

! The centrifugal reduction of Matt et al. (2012), relative to the Sun:
! ((c_2**2 + f_sun)/(c_2**2 + f))**excen with f = omega**2 R**3/(2 G M).
real(dp) function matt_centrifugal_factor(omega, fsun, log10_radius, total_mass_msun)
      real(dp), intent(in) :: omega, fsun, log10_radius, total_mass_msun
      real(dp) :: f
      f = 0.5*omega**2*exp(ln10*(3.0*log10_radius-rotation_settings%log10_gravitational_constant))/ &
           total_mass_msun/rotation_settings%solar_mass_cgs
      matt_centrifugal_factor = pow(((rotation_settings%wind_centrifugal_constant**2+fsun)/(rotation_settings%wind_centrifugal_constant**2+f)), rotation_settings%wind_centrifugal_exponent)
end function matt_centrifugal_factor

! tau_cz interpolated in the step (MHP 8/17)
real(dp) function interpolated_turnover_timescale(conditions)
      type(wind_conditions), intent(in) :: conditions
      interpolated_turnover_timescale = conditions%turnover_timescale_start+ &
           conditions%step_fraction*(conditions%turnover_timescale_end-conditions%turnover_timescale_start)
end function interpolated_turnover_timescale

end module magnetic_braking
