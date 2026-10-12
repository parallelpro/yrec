!----------------------------------------------------------------------
! secular_transport
!----------------------------------------------------------------------
! SECULR CALCULATES THE TRANSFER OF ANGULAR MOMENTUM(AND THE ASSOCIATED
! TRANSPORT OF COMPOSITION) DUE TO SECULAR ROTATIONAL INSTABILITIES.
!
! A few local/dummy quantities are named differently by the various
! callees this file dispatches to. Where that
! happens the name used here is the one shared by the most callees (or
! the more descriptive of the two), and the alternate name(s) are
! noted alongside the declaration, following henyey_iteration's convention:
!  - dynamical_shear_omega_limit (QWRMAX): circulation_velocities calls the same
!    quantity dlnomega_dlnr_max.
!  - radius_unlogged (HRU): circulation_velocities/zahn_coupling_factor/compute_quadrupole/am_transport_grid
!    call it radius; diffusion_velocity_scales calls it radius_mid_prev.
!  - mixing_diffusion_coeff (COD2, dummy arg): diffuse_composition_driver calls it
!    diffusion_coeff.
!  - diffusion_velocity (HV, dummy arg): circulation_velocities calls it
!    total_circulation_velocity.
!  - eq_mixing_diffusion_coeff (ECOD2): diffuse_composition_driver calls it
!    equally_spaced_diffusion_coeff.
!  - eq_angular_momentum (EJ): equal_grid_to_model calls it angular_momentum.
!  - eq_am_diffusion_coeff (ECOD): am_diffusion_coeffs calls it diffusion_coeff.
!  - eq_delta_angular_momentum (DJ): equal_grid_to_model calls it
!    delta_angular_momentum; solve_angular_momentum_transport calls it dj.
!  - eq_mass (EM): diffuse_composition_driver calls it equally_spaced_mass.
!  - diffusion_converged (LOKAD): check_angular_momentum calls it
!    already_converged_flag.
!  - wind_loss_explicit/wind_loss_implicit (WIND1/WIND2):
!    magnetic_braking's compute_solid_body_braking calls the same pair
!    domega_start/domega_end.
! The midpoint structure arrays live in rot_scr (rotation_scratch_lib, filled by
! mid_timestep_model each sub-step) and the star-model arrays/scalars
! (log_mass/m/dm/j_rot/xa/nz/log_L/log_total_mass/star_mass/log_Teff)
! are read from star directly. What remains is genuine per-call data:
! the sub-step, the saved-at-substep-start angular momentum, the
! unstable-domain bounds, the surface-CZ terms, and the outputs.
! Wrapped as a module procedure so the compiler enforces the call.
module secular_transport
      use rotation_configuration, only: rotation_settings
      use transport_checks, only: check_angular_momentum, check_composition
      use am_grid, only: am_transport_grid, equal_grid_to_model
      use rotational_mixing, only: diffuse_composition_driver
      use circulation, only: circulation_velocities, diffusion_velocity_scales, zahn_coupling_factor, compute_quadrupole
      use am_transport_coefficients, only: am_diffusion_coeffs, am_advection_diffusion_coeffs
      use omega_profile, only: solid_body_omega
      use stability_coefficients, only: rotation_stability_setup
      implicit none
contains

subroutine apply_secular_transport(sub_timestep, specific_angular_momentum_saved, &
     zone_min, zone_max, redo_flag, cut_count, &
     cz_moment_of_inertia, cz_mass_bottom, cz_mass_top, omega_surface, &
     surface_cz_active, mixing_diffusion_coeff, diffusion_velocity, &
     diffusion_solve_ok, ierr)
      use rotation_scratch_lib
      use magnetic_braking, only: brake_convection_zone, compute_solid_body_braking
      use convection_lib, only: homogenize_convection_zones
      use star_info_lib, only: star, json
      use physical_constants
      use numerics_lib
      use math_lib
      implicit none
      double precision, intent(inout) :: sub_timestep
      double precision, intent(in) :: specific_angular_momentum_saved(json)
      integer, intent(in) :: zone_min, zone_max
      logical, intent(out) :: redo_flag
      integer, intent(inout) :: cut_count
      double precision, intent(out) :: cz_moment_of_inertia
      double precision, intent(in) :: cz_mass_bottom
      double precision, intent(in) :: cz_mass_top
      double precision, intent(inout) :: omega_surface
      logical, intent(in) :: surface_cz_active
      double precision, intent(inout) :: mixing_diffusion_coeff(json)
      double precision, intent(inout) :: diffusion_velocity(json)
      logical, intent(inout) :: diffusion_solve_ok
      integer, intent(out) :: ierr

! --- locals (JSON-dimensioned) ---
      double precision :: dynamical_shear_omega_limit(json), eq_mass(json), &
           am_diffusion_coeff(json), eq_mixing_diffusion_coeff(json), &
           eq_moment_of_inertia(json), eq_angular_momentum(json), &
           eq_am_diffusion_coeff(json), radius_unlogged(json), radius_mid(json), &
           dlnomega_dlnr(json), eq_delta_angular_momentum(json), &
           specific_angular_momentum_prev(json), omega_start(json), &
           eq_omega(json)

! Tridiagonal-solve work arrays (Thomas algorithm) threaded from
! am_diffusion_coeffs (fills them) to solve_angular_momentum_transport
! (solves). Its solution(:) (the new omega distribution) is never read back here --
! this call site only uses eq_delta_angular_momentum/
! sum_delta_angular_momentum -- so it is captured in an unused local.
! surface_wind_loss_term is am_diffusion_coeffs' surface wind
! angular-momentum-loss term, consumed by solve_angular_momentum_transport.
      double precision :: sub_diag(json), diag(json), super_diag(json), &
           rhs(json), unused_tridia_solution(json), surface_wind_loss_term

! --- other locals ---
      logical :: disk_lock_active
      logical :: lcz_first_zone, lcz_last_zone
      double precision :: total_luminosity
      integer :: i, j
      double precision :: log_radius_center
! dlnr_weight (originally DR1) is the sum of the two one-sided inverse
! log-radius spacings used to build the centered dlnomega/dlnr
! estimate; not otherwise named elsewhere.
      double precision :: dlnr_weight
      integer :: iteration
      double precision :: omega_surface_start
      logical :: wind_loss_active
      integer :: solid_body_zone_start, solid_body_zone_end
      logical :: in_unstable_region
      integer :: scan_start_zone
      logical :: unstable_zone_found
      integer :: zone_begin, zone_end
      double precision :: grid_spacing
      logical :: single_interface_flag
      logical :: fix_omega_at_surface
      logical :: diffusion_converged
      double precision :: wind_loss_explicit, wind_loss_implicit
      double precision :: sum_delta_angular_momentum
      double precision :: solid_cz_mass_bottom, solid_cz_mass_top
      integer :: solid_start_zone
      integer :: species_begin, species_end

      ierr = 0

! The advective transport option (LDIFAD,
! use_diffusion_advection_transport) has never been functional. The
! equal-grid coefficients am_advection_diffusion_coeffs reads
! (eq_velocity_coeff*, geometric_factor, am_2nd/3rd_deriv_coeff,
! dchi_dr_edge, mixing_geometric_factor, equatorial_radius) are written
! nowhere, and neither were their F77 originals (EV*, FGEOM, ECOD5/6,
! /DIFAD4/, /EGRIDCHI/ -- only dadcoeft.f referenced them), so the path
! divides by a zero geometric factor. The flag is inert unless rotation
! and instability transport are on, which is why the many non-rotating
! decks that set it run fine; here, where it would take effect, the run
! stops instead.
      if (rotation_settings%use_diffusion_advection_transport) then
         write(*,'(a)') 'apply_secular_transport: use_diffusion_advection_transport' &
              //' (LDIFAD) is not functional -- its advective transport' &
              //' coefficients are never computed. Set it .false.'
         ierr = 1
         return
      end if

! MHP 9/94
! DISK LOCKING CHECKED
      disk_lock_active = .false.
      if(disk_locking_engaged()) disk_lock_active = .true.
!  FIRST AND LAST SHELL BOUNDARY CONDITIONS REQUIRE THAT THESE SHELLS BE
!  TREATED AS 'CONVECTIVE' FOR DIFFUSION PURPOSES. FIX THIS.
      lcz_first_zone = rot_scr%am_transport_convective_flag_mid(1)
      lcz_last_zone = rot_scr%am_transport_convective_flag_mid(star%nz)
!  SET UP FACTORS USED TO COMPUTE DIFFUSION VELOCITIES.
      call rotation_stability_setup(rot_scr%log_density_mid,rot_scr%mean_gravity_mid,rot_scr%luminosity_lsun_mid,rot_scr%log_pressure_mid, &
           rot_scr%log_radius_mid,star%m,rot_scr%log_temperature_mid,zone_min,zone_max, &
           star%nz,radius_unlogged,dynamical_shear_omega_limit)
      total_luminosity = star%solar_luminosity_cgs*rot_scr%luminosity_lsun_mid(star%nz)
      call compute_omega_gradients()
      call compute_quadrupole(rot_scr%log_density_mid,rot_scr%mean_gravity_mid,radius_unlogged,rot_scr%omega_mid,star%nz,ierr)
      if (ierr /= 0) return
      do i = 1,star%nz
         star%vfc(i) = 0.0D0
      end do
!  CHECK STABILITY OF THE MODEL AND COMPUTE CIRCULATION VELOCITIES.
!  MECHANISMS CONSIDERED:SECULAR AND DYNAMICAL SHEAR,EDDINGTON CIRCULATION,
!  AND THE GSF INSTABILITY.
      iteration = 1
      call compute_circulation()
!  STORE INITIAL ANGULAR MOMENTUM DISTRIBUTION.
      do i =1,star%nz
         specific_angular_momentum_prev(i) = specific_angular_momentum_saved(i)
         omega_start(i) = rot_scr%omega_mid(i)
      end do
!  STORE INITIAL SURFACE ANGULAR VELOCITY FOR USE IN ANGULAR MOMENTUM
!  LOSS CALCULATIONS.
      omega_surface_start = rot_scr%omega_substep_start(star%nz)
      diffusion_solve_ok = .false.
!  ON THE FIRST LEVEL OF ITERATION, THE UNPERTURBED MODEL IS USED TO
!  CALCULATE THE DIFFUSION VELOCITIES. ON THE SECOND AND SUBSEQUENT
!  ITERATIONS, THE ANGULAR MOMENTUM DISTRIBUTION FROM THE PREVIOUS
!  ITERATION IS USED TO GET THE VELOCITIES.  THIS 'NEW' VELOCITY IS THEN
!  AVERAGED WITH THE VELOCITY FOUND IN THE PREVIOUS ITERATION TO GET A
!  CORRECTED V AND THUS A MORE ACCURATE RUN OF DIFFUSION COEFFICIENTS.
      do iteration = 1,rotation_settings%max_diffusion_iters
         omega_surface = rot_scr%omega_mid(star%nz)
         if(iteration.gt.1) call restart_iteration()
         call spin_down_uncoupled_surface()
         if (ierr /= 0) return
!  IF LDO=F,NO INSTABILITIES OCCUR (STABLE AGAINST ALL MECHANISMS).
         if(.not.unstable_zone_found) return
!  TREAT CENTRAL AND SURFACE ZONES AS ALWAYS CONVECTIVE
!  (SHOULD BE FIXED TO GIVE BETTER CENTRAL/SURFACE B.C.)
         rot_scr%am_transport_convective_flag_mid(1) = .true.
         rot_scr%am_transport_convective_flag_mid(star%nz) = .true.
!  UNSTABLE REGION EXISTS.
!  FIND DIFFUSION COEFFICIENTS(COD) FOR ALL UNSTABLE REGIONS.
         call diffusion_velocity_scales(radius_unlogged,star%nz,radius_mid, &
              am_diffusion_coeff,mixing_diffusion_coeff)
! The model-grid transport coefficients are profile columns (D_omega,
! D_mix); the step's last sub-step wins.
         do i = 1, star%nz
            star%am_diffusion_coeff(i) = am_diffusion_coeff(i)
            star%mixing_diffusion_coeff(i) = mixing_diffusion_coeff(i)
         end do
!  EACH UNSTABLE REGION IS SOLVED SEPARATELY STARTING HERE.
         in_unstable_region = .false.
         scan_start_zone = zone_min
         region_loop: do
         call find_unstable_region()
!  IF NO NON-ZERO V'S ENCOUNTERED, EXIT.
         if(.not.unstable_zone_found) exit region_loop
!  TRANSFORM TO EQUAL GRID SPACING IN CHI FOR THE REGION.
! CHI IS A NORMALIZED SUM OF THE VARIABLES USED TO PLACE POINTS
! IN THE HENYEY SCHEME, CHOSEN SUCH THAT THE GRID USED IN THE
! ANGULAR MOMENTUM EVOLUTION IS CLOSE TO THE GRID STORED FOR
! THE STRUCTURAL EVOLUTION.
         call am_transport_grid(am_diffusion_coeff,mixing_diffusion_coeff,rot_scr%log_density_mid, &
              rot_scr%moment_of_inertia_mid,specific_angular_momentum_saved,rot_scr%luminosity_lsun_mid, &
              rot_scr%log_pressure_mid,rot_scr%log_radius_mid,radius_unlogged,star%log_mass,star%m, &
              star%dm,star%log_total_mass,zone_begin,zone_end, &
              rot_scr%am_transport_convective_flag_mid,star%nz,omega_start,grid_spacing, &
              eq_am_diffusion_coeff,eq_mixing_diffusion_coeff, &
              eq_moment_of_inertia,eq_angular_momentum,eq_mass,eq_omega, &
              single_interface_flag)
!  LDUM2=T IF TWO C.Z.'S ARE SEPARATED BY ONE RADIATIVE ZONE;
!  SKIP IF THIS OCCURS.
         if(single_interface_flag) then
            if(scan_start_zone.le.zone_max) then
               cycle region_loop
            else
               exit region_loop
            endif
         endif
! MHP 3/09 SKIP ANGULAR MOMENTUM EVOLUTION FOR SOLID BODY MODELS
! JNT 09/25 FOR 05/15 IMPJMOD=1 IS THE SAME AS LSOLID
         if(.not.rotation_settings%force_solid_body_rotation .and. &
              (rotation_settings%solid_body_mode_flag.ne.1))then
            call transport_angular_momentum_in_region()
            if (ierr /= 0) return
         endif
!  PERFORM COMPOSITION DIFFUSION.
!  UNTIL THE FINAL ITERATION, ONLY COMPOSITION DIFFUSION OF SPECIES WHICH
!  AFFECT GRADIENTS IN MEAN MOLECULAR WEIGHT IS COMPUTED (H,HE3,HE4).
!  ON THE FINAL ITERATION, DIFFUSION OF ALL SPECIES IS PERFORMED.
         species_begin = 1
         species_end = 4
         call mix_composition()
         if (ierr /= 0) return
!  RETURN FOR NEXT REGION IF APPLICABLE
         if(scan_start_zone.gt.zone_max) exit region_loop
         end do region_loop
! CHECK SOLUTION,UPDATE OMEGA,AND SEE IF ANOTHER ITERATION IS NEEDED.
         rot_scr%am_transport_convective_flag_mid(1) = lcz_first_zone
         rot_scr%am_transport_convective_flag_mid(star%nz) = lcz_last_zone
         call check_angular_momentum(rot_scr%log_density_mid,specific_angular_momentum_prev, &
              specific_angular_momentum_saved,rot_scr%log_radius_mid,star%log_mass,star%dm, &
              iteration, &
              rot_scr%am_transport_convective_flag_mid,star%nz,sub_timestep,rot_scr%eta_squared_mid, &
              rot_scr%moment_of_inertia_mid,star%j_rot,cut_count, &
              diffusion_solve_ok,redo_flag,rot_scr%omega_mid,rot_scr%di_domega_mid,rot_scr%mean_radius_mid, &
              diffusion_converged, ierr)
         if (ierr /= 0) return
! CHECK COMPOSITION DIFFUSION AND RECOMPUTE MEAN MOLECULAR WEIGHT.
         if(.not.redo_flag)call check_composition(star%xa,iteration, &
              star%nz,sub_timestep,cut_count, &
              diffusion_solve_ok,redo_flag, ierr)
         if (ierr /= 0) return
! MHP 9/93
         if(rotation_settings%no_am_transport_in_core)diffusion_solve_ok = .true.
! IF LOK=T,CONVERGED.
         if(diffusion_solve_ok)exit
! IF LREDO=T, A PROBLEM REQUIRES TIMESTEP CUTTING.
         if(redo_flag)return
      end do
! IF THE ITERATION LIMIT WAS REACHED WITHOUT CONVERGENCE, THE LAST
! ITERATE IS ACCEPTED (NO REDO, NO ERROR) AND DIFFUSION OF THE
! REMAINING SPECIES PROCEEDS; check_angular_momentum ONLY PRINTS.
      diffusion_solve_ok = .true.
! PERFORM COMPOSITION DIFFUSION OF REMAINING SPECIES.
      rot_scr%am_transport_convective_flag_mid(1) = .true.
      rot_scr%am_transport_convective_flag_mid(star%nz) = .true.
      species_begin = 5
      if(rotation_settings%use_extended_composition)then
         species_end = 15
      else
         species_end = 11
      endif
      call mix_composition()
      if (ierr /= 0) return
      rot_scr%am_transport_convective_flag_mid(1) = lcz_first_zone
      rot_scr%am_transport_convective_flag_mid(star%nz) = lcz_last_zone

contains

! ---------------------------------------------------------------
! d omega/d ln r at each interface of the transport region, centered
! (weights from the one-sided log-radius spacings).
subroutine compute_omega_gradients()
      do i = zone_min,zone_max
         log_radius_center = log10(rot_scr%radius_interface(i))
         dlnr_weight = 1.0D0/ln10/(rot_scr%log_radius_mid(i)-log_radius_center)+ &
              1.0D0/ln10/(log_radius_center-rot_scr%log_radius_mid(i-1))
         dlnomega_dlnr(i) = 0.25D0*(rot_scr%omega_mid(i)-rot_scr%omega_mid(i-1))*dlnr_weight
      end do
end subroutine compute_omega_gradients

! ---------------------------------------------------------------
! The circulation and instability velocities of this iteration (and
! whether any transport acts), and the Zahn coupling factor under lvfc
! (as in the F77 seculr.f; a 2026 change had gated it on ldifad, which
! left star%vfc zero and killed the rotational mixing in
! diffusion_velocity_scales, fixed in the 2026 bug sweep).
subroutine compute_circulation()
      call circulation_velocities(rot_scr%log_radius_mid,radius_unlogged,zone_min,zone_max,iteration, &
           rot_scr%am_transport_convective_flag_mid,star%nz,rot_scr%omega_mid,unstable_zone_found, &
           dlnomega_dlnr,dynamical_shear_omega_limit,diffusion_velocity, &
           total_luminosity,sub_timestep,rot_scr%log_pressure_mid)
      if(rotation_settings%use_zahn_coupling_factor) &
           call zahn_coupling_factor(rot_scr%log_density_mid,radius_unlogged,diffusion_velocity, &
           zone_min,zone_max,rot_scr%omega_mid)
end subroutine compute_circulation

! ---------------------------------------------------------------
! The second and later iterations: omega averaged with its value at the
! start, the velocities recomputed from it, then the angular momentum
! and composition reset to their start and the convection zones remixed
! to their proper depth (otherwise a deepening convection zone takes the
! local abundance at its old base; MHP 10/91).
subroutine restart_iteration()
      do i = 1,star%nz
         rot_scr%omega_mid(i) = 0.5D0*(rot_scr%omega_mid(i)+omega_start(i))
      end do
      call compute_omega_gradients()
      call compute_circulation()
      do i = 1,star%nz
         star%j_rot(i) = specific_angular_momentum_saved(i)
         do j = 1,4
            star%xa(j,i) = rot_scr%composition_snapshot(j,i)
         end do
      end do
      call homogenize_convection_zones(star%xa,star%dm,rot_scr%am_transport_convective_flag_mid, &
           star%nz)
end subroutine restart_iteration

! ---------------------------------------------------------------
! The wind on a surface convection zone the diffusion does not couple to
! the interior (no velocity at its base), or on the whole star when it
! rotates as a solid body (solid_body_mode_flag = 1 acts like
! force_solid_body_rotation).
subroutine spin_down_uncoupled_surface()
      if(surface_cz_active.and.diffusion_velocity(zone_max).eq.0.0D0 &
           .and. .not.disk_lock_active)then
! the moment of inertia of the surface convection zone
         cz_moment_of_inertia = 0.0D0
         do i = zone_max,star%nz
            cz_moment_of_inertia = cz_moment_of_inertia + rot_scr%moment_of_inertia_mid(i)
         end do
         wind_loss_active = rotation_settings%use_wind_torque
         call brake_convection_zone(star%log_L,sub_timestep,cz_mass_bottom, &
              cz_mass_top,zone_max,star%nz,wind_loss_active,omega_surface, &
              star%star_mass,star%log_Teff,cz_moment_of_inertia, &
              star%j_rot, ierr)
         if (ierr /= 0) return
      else if(surface_cz_active .and. (rotation_settings%force_solid_body_rotation .or. &
           (rotation_settings%solid_body_mode_flag.eq.1)))then
         wind_loss_active = rotation_settings%use_wind_torque
         solid_cz_mass_bottom = 0.0D0
         solid_cz_mass_top = exp(ln10*star%log_total_mass)
         solid_start_zone = 1
         call brake_convection_zone(star%log_L,sub_timestep,solid_cz_mass_bottom, &
              solid_cz_mass_top,solid_start_zone,star%nz,wind_loss_active, &
              omega_surface,star%star_mass,star%log_Teff,cz_moment_of_inertia, &
              star%j_rot, ierr)
         if (ierr /= 0) return
         solid_body_zone_start = 1
         solid_body_zone_end = star%nz
         call solid_body_omega(rot_scr%log_density_mid,star%j_rot,rot_scr%log_radius_mid, &
              star%log_mass,star%dm,solid_body_zone_start,solid_body_zone_end, &
              rot_scr%eta_squared_mid,rot_scr%moment_of_inertia_mid,rot_scr%omega_mid,rot_scr%di_domega_mid,rot_scr%mean_radius_mid,star%nz)
      endif
end subroutine spin_down_uncoupled_surface

! ---------------------------------------------------------------
! The next unstable region, from scan_start_zone on: zone_begin is the
! zone below its first nonzero velocity, zone_end the zone above its
! last; unstable_zone_found is false when there is none.
subroutine find_unstable_region()
      unstable_zone_found = .false.
      do j = scan_start_zone,zone_max
         if(diffusion_velocity(j).gt.0.0D0) then
            unstable_zone_found = .true.
            if(.not.in_unstable_region) then
               in_unstable_region = .true.
               zone_begin = j - 1
            endif
         else if(in_unstable_region) then
            zone_end = j - 1
            in_unstable_region = .false.
            scan_start_zone = j + 1
            exit
         endif
      end do
      if (j .gt. zone_max) then
! the last interface is unstable: the region ends at zone_max
      if(in_unstable_region) zone_end = zone_max
      scan_start_zone = zone_max + 1
      end if
end subroutine find_unstable_region

! ---------------------------------------------------------------
! Angular momentum transport across the region: the wind terms when the
! surface convection zone is part of it (none under disk locking, which
! fixes the surface omega), the diffusion equation on the equally spaced
! grid and its solve (the advection-diffusion branch is the never
! functional ldifad path), and the change of angular momentum mapped
! back to the model.
subroutine transport_angular_momentum_in_region()
!  CHECK IF SURFACE C.Z. IS PART OF THE UNSTABLE REGION.
!  IF SO,CALCULATE TERMS FOR DIFFUSION.
      if(surface_cz_active.and.zone_end.eq.zone_max)then
         if(disk_lock_active)then
            wind_loss_explicit = 0.0D0
            wind_loss_implicit = 0.0D0
            fix_omega_at_surface = .true.
         else
            fix_omega_at_surface = .false.
            if(rotation_settings%use_wind_torque)then
               cz_moment_of_inertia = eq_moment_of_inertia(rot_scr%num_eq_points)
               omega_surface = rot_scr%omega_mid(star%nz)
               call compute_solid_body_braking(star%log_L,sub_timestep, &
                    cz_moment_of_inertia,iteration,omega_surface, &
                    star%star_mass,star%log_Teff,omega_surface_start, &
                    wind_loss_explicit,wind_loss_implicit, ierr)
               if (ierr /= 0) return
            else
               wind_loss_explicit = 0.0D0
               wind_loss_implicit = 0.0D0
            endif
         endif
      else
         fix_omega_at_surface = .false.
         wind_loss_explicit = 0.0D0
         wind_loss_implicit = 0.0D0
      endif
      diffusion_converged = .false.
      if(.not.rotation_settings%use_diffusion_advection_transport)then
!  SET UP DIFFUSION EQUATION ARRAYS TO SOLVE FOR OMEGA AT END OF TSTEP
         call am_diffusion_coeffs(eq_am_diffusion_coeff,grid_spacing,sub_timestep, &
              eq_moment_of_inertia,eq_angular_momentum,eq_omega,rot_scr%num_eq_points, &
              wind_loss_explicit,wind_loss_implicit,fix_omega_at_surface, &
              sub_diag,diag,super_diag,rhs,surface_wind_loss_term)
!  SOLVE MATRIX FOR THE RUN OF OMEGA AT THE END OF THE TIMESTEP AT THE
!  EQUALLY SPACED GRID POINTS.
         call solve_angular_momentum_transport(rot_scr%num_eq_points,eq_moment_of_inertia, &
              eq_delta_angular_momentum,sum_delta_angular_momentum, &
              sub_diag,diag,super_diag,rhs,unused_tridia_solution, &
              surface_wind_loss_term, ierr)
         if (ierr /= 0) return
      else
! SOLVE FOR OMEGA AND ITS DERIVATIVES IN A BAND MATRIX
         call am_advection_diffusion_coeffs(grid_spacing,sub_timestep,eq_moment_of_inertia, &
              eq_omega,rot_scr%num_eq_points,wind_loss_explicit,wind_loss_implicit, &
              eq_delta_angular_momentum,eq_mixing_diffusion_coeff, &
              sum_delta_angular_momentum,fix_omega_at_surface, &
              diffusion_converged, ierr)
         if (ierr /= 0) return
      endif
!  TRANSFORM THE NEW ANGULAR MOMENTUM DISTRIBUTION BACK TO THE ORIGINAL GRID
!  POINTS IN THE UNSTABLE REGION.
      call equal_grid_to_model(eq_delta_angular_momentum,eq_angular_momentum, &
           star%dm,zone_begin,zone_end,rot_scr%am_transport_convective_flag_mid, &
           star%nz,sum_delta_angular_momentum,star%j_rot)
end subroutine transport_angular_momentum_in_region

! ---------------------------------------------------------------
! The rotational mixing of species species_begin..species_end across the
! transport region.
subroutine mix_composition()
      call diffuse_composition_driver(sub_timestep,mixing_diffusion_coeff, &
           eq_mixing_diffusion_coeff,eq_mass,rot_scr%log_density_mid,rot_scr%luminosity_lsun_mid, &
           rot_scr%log_pressure_mid,rot_scr%log_radius_mid,star%log_mass,star%m,star%dm, &
           star%log_total_mass,diffusion_velocity,zone_begin,zone_end,zone_max, &
           zone_min,rot_scr%am_transport_convective_flag_mid,diffusion_solve_ok, &
           star%nz,star%xa,species_begin,species_end, ierr)
end subroutine mix_composition

end subroutine apply_secular_transport

!----------------------------------------------------------------------
! solve_angular_momentum_transport (formerly tridia, numerics_lib)
!----------------------------------------------------------------------
! The Thomas-algorithm solve of the angular-momentum diffusion system
! (sub_diag/diag/super_diag, right-hand side rhs = omega at the start)
! for omega at the end (solution), with dj(:) = (omega_new - omega_old)
! * ei(:) per point (ei the moment of inertia) and sumdj their sum.
! dj(n) starts from dj_n_seed, the surface wind loss (formerly smuggled
! in through the shared common block's gamma_elim(n)). Before the back
! substitution, omega(n) is reset to rhs(n)*(1 + (omega(n)-rhs(n))/rhs(n)).
! KC 2025-05-31 removed the unused ej argument.
subroutine solve_angular_momentum_transport(n, ei, dj, sumdj, sub_diag, diag, super_diag, rhs, &
     solution, dj_n_seed, ierr)
      use star_info_lib, only: json

      integer, intent(in) :: n
      double precision, intent(in) :: ei(json)
      double precision, intent(out) :: dj(json)
      double precision, intent(out) :: sumdj
! sub_diag/diag/super_diag are the tridiagonal matrix's three
! diagonals and rhs is the right-hand side, all filled in by the
! caller; solution is solved for here. dj_n_seed is dj(n)'s initial
! value (see header note above).
      double precision, intent(in) :: sub_diag(json), diag(json), &
           super_diag(json), rhs(json)
      double precision, intent(out) :: solution(json)
      double precision, intent(in) :: dj_n_seed
      integer, intent(out) :: ierr

! gamma_elim is solver-internal scratch (no SAVE: nothing here depends
! on it persisting across calls, now that dj_n_seed carries the one
! value that used to be read from it across calls).
      double precision :: gamma_elim(json)

      double precision :: rhs_orig(json)
      integer :: i, j
      double precision :: bet, fj

      ierr = 0
      dj(n) = dj_n_seed
      do i = 1, n
         rhs_orig(i) = rhs(i)
      end do
      bet = diag(1)
      solution(1) = rhs(1)/bet
      do j = 2, n
         gamma_elim(j) = super_diag(j-1)/bet
         bet = diag(j) - sub_diag(j)*gamma_elim(j)
         if (bet.eq.0) then
            write(*,*) 'solve_angular_momentum_transport: singular matrix'
            ierr = 1
            return
         end if
         solution(j) = (rhs(j) - sub_diag(j)*solution(j-1))/bet
      end do
      dj(n) = dj(n)+(solution(n)-rhs_orig(n))*ei(n)
      sumdj = dj(n)
      fj = 1.0d0+(solution(n)-rhs_orig(n))/rhs_orig(n)
      solution(n) = rhs_orig(n)*fj
      do j = n-1, 1, -1
         solution(j) = solution(j) - gamma_elim(j+1)*solution(j+1)
         dj(j) = (solution(j)-rhs_orig(j))*ei(j)
         sumdj = sumdj + dj(j)
      end do

      return
end subroutine solve_angular_momentum_transport

end module secular_transport
