!----------------------------------------------------------------------
! omega_profile
!----------------------------------------------------------------------
! The rotation curve omega(r) from the specific angular momentum j.
! omega_from_j walks the star from the center out: each radiative zone
! is solved alone (solid_body_omega), and each convection zone as a
! whole under the imposed rotation law (enforce_rotation_profile).
module omega_profile
      use rotation_configuration, only: rotation_settings, disk_locking_engaged
      implicit none
      private
      public :: omega_from_j, solid_body_omega, enforce_rotation_profile, &
           default_inertia_tolerance

! The relative tolerance on J of solid_body_omega's Newton iteration
! when the caller gives none: zero, so the iteration always takes its
! 20 corrections. It is the F77 TOLERI, which nothing ever set;
! check_angular_momentum loosens it (inertia_tolerance) for the
! transport's iterations before the last.
      double precision, parameter :: default_inertia_tolerance = 0.0d0

contains

! ---------------------------------------------------------------
! omega_from_j (legacy getrot): the rotation curve for a given
! distribution of specific angular momentum j (cm^2/s).
! Each radiative zone (am_transport_convective_flag false) gets the
! omega at which I(omega) omega = j dm (solid_body_omega). Each run of
! convective zones (overshoot included) goes to enforce_rotation_profile,
! which imposes the convection-zone rotation law on the run as a whole.
! Updates omega (rad/s), j, and each zone's shape and inertia:
! eta_squared (Radau's eta_2), mean_radius (r0), moment_of_inertia (I)
! and qiw (di_domega of zone_moments_of_inertia). disk_gate_age_gyr,
! the age the disk-locking test uses (the last converged model's), and
! inertia_tolerance and shape_tolerance (optional) go to every
! solid_body_omega call, and the latter on to shape.
subroutine omega_from_j(log_density, specific_angular_momentum, log_radius, &
     log_mass, shell_mass, am_transport_convective_flag, num_zones, &
     eta_squared, moment_of_inertia, omega, qiw, mean_radius, &
     disk_gate_age_gyr, inertia_tolerance, &
     shape_tolerance)
      use composition_layout, only: json

      implicit none

      double precision, intent(inout) :: log_density(json), &
           specific_angular_momentum(json), log_radius(json), log_mass(json), &
           shell_mass(json)
      logical, intent(in) :: am_transport_convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: eta_squared(json), &
           moment_of_inertia(json), omega(json), qiw(json), mean_radius(json)
      double precision, intent(in) :: disk_gate_age_gyr
      double precision, intent(in), optional :: inertia_tolerance
      double precision, intent(in), optional :: shape_tolerance
      logical :: in_convection_zone
      integer :: zone_index, zone_start, zone_end

      in_convection_zone = .false.
      do zone_index = 1,num_zones
         if (.not.in_convection_zone) then
            if (.not.am_transport_convective_flag(zone_index)) then
! A radiative zone: solve it alone.
               zone_start = zone_index
               zone_end = zone_index
               call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
                    log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                    moment_of_inertia,omega,qiw,mean_radius,num_zones, &
                    disk_gate_age_gyr, &
                    inertia_tolerance=inertia_tolerance, &
                 shape_tolerance=shape_tolerance)
            else
! The first zone of a convection zone: note where it starts.
               in_convection_zone = .true.
               zone_start = zone_index
            endif
         else if (.not.am_transport_convective_flag(zone_index)) then
! The convection zone ended at the zone below: impose the rotation law
! on all of it.
            zone_end = zone_index - 1
            call enforce_rotation_profile(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones, &
                 disk_gate_age_gyr, &
                 inertia_tolerance=inertia_tolerance, &
                 shape_tolerance=shape_tolerance)
! Then solve the current, radiative, zone alone.
            zone_start = zone_index
            zone_end = zone_index
            call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones, &
                 disk_gate_age_gyr, &
                 inertia_tolerance=inertia_tolerance, &
                 shape_tolerance=shape_tolerance)
            in_convection_zone = .false.
         else if (zone_index.eq.num_zones) then
! The convection zone reaches the surface: impose the rotation law on
! it as above. (A convection zone made of the surface zone alone starts
! at zone_index = num_zones and never reaches this branch, so that
! zone's omega is left as it was.)
            zone_end = num_zones
            call enforce_rotation_profile(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones, &
                 disk_gate_age_gyr, &
                 inertia_tolerance=inertia_tolerance, &
                 shape_tolerance=shape_tolerance)
         endif
      end do
      return
end subroutine omega_from_j

! ---------------------------------------------------------------
! solid_body_omega: the single rotation rate omega of zones
! zone_start..zone_end that carries their total angular momentum
! J = sum_i dm_i j_i. Newton iteration on J(omega) = omega I(omega):
! at each guess, shape gives the zones' r0 (mean_radius) and eta_2
! (eta_squared), and zone_moments_of_inertia their I and di_domega.
! The iteration stops when |J - omega sum_i I_i| <= J times
! inertia_tolerance (default_inertia_tolerance if absent), or after 20
! corrections (then
! the last guess is kept). On return omega, I, di_domega, r0 and eta_2
! belong to the last guess, and j_i = I_i omega / dm_i.
! Disk locking (MHP 9/94): while disk_locking_engaged(disk_gate_age_gyr), a surface
! convection zone of more than one zone is held at
! rotation_settings%disk_omega_rad_s with no iteration, and its j is
! reset to match; the disk takes up the difference in J.
subroutine solid_body_omega(log_density, specific_angular_momentum, log_radius, &
     log_mass, shell_mass, zone_start, zone_end, eta_squared, &
     moment_of_inertia, omega, di_domega, mean_radius, num_zones, &
     disk_gate_age_gyr, inertia_tolerance, &
     shape_tolerance)
      use moments_of_inertia, only: zone_moments_of_inertia
      use rotating_shape, only: shape

      use composition_layout, only: json
      implicit none

      double precision, intent(in) :: log_density(json)
      double precision, intent(inout) :: specific_angular_momentum(json)
      double precision, intent(in) :: log_radius(json), log_mass(json), &
           shell_mass(json)
      integer, intent(in) :: zone_start, zone_end
      double precision, intent(inout) :: eta_squared(json)
      double precision, intent(out) :: moment_of_inertia(json)
      double precision, intent(inout) :: omega(json)
      double precision, intent(out) :: di_domega(json)
      double precision, intent(out) :: mean_radius(json)
      integer, intent(in) :: num_zones
      double precision, intent(in) :: disk_gate_age_gyr
      double precision, intent(in), optional :: inertia_tolerance
      double precision, intent(in), optional :: shape_tolerance
! --- locals ---
      double precision :: tolerance
      double precision :: total_angular_momentum, omega_sum
      integer :: iteration_count
      logical :: disk_locked
      double precision :: omega_guess
      double precision :: total_moment_of_inertia, total_di_domega
      double precision :: new_angular_momentum, delta_angular_momentum, &
           delta_omega
      integer :: zone_idx

      tolerance = default_inertia_tolerance
      if (present(inertia_tolerance)) tolerance = inertia_tolerance
      total_angular_momentum = 0.0d0
      omega_sum = 0.0d0
      iteration_count = 0
! The total angular momentum J of the zones, and the mean of their
! current omega as the first guess.
      do zone_idx = zone_start,zone_end
         total_angular_momentum = total_angular_momentum + &
              shell_mass(zone_idx)*specific_angular_momentum(zone_idx)
         omega_sum = omega_sum + omega(zone_idx)
      end do
! Disk locking (MHP 9/94) applies to the surface convection zone only:
! a run of more than one zone that ends at the surface.
      disk_locked = .false.
      if(zone_end.eq.num_zones .and. zone_start.ne.zone_end .and. &
           disk_locking_engaged(disk_gate_age_gyr))disk_locked = .true.
      if(.not.disk_locked)then
         omega_guess = omega_sum/dfloat(zone_end - zone_start + 1)
      else
         omega_guess = rotation_settings%disk_omega_rad_s
      end if
      omega_iter: do
      do zone_idx = zone_start,zone_end
         omega(zone_idx) = omega_guess
      end do
! The shape and the moments of inertia of the zones at omega_guess.
      call shape(log_density,log_radius,log_mass,zone_start,zone_end,omega, &
           eta_squared,mean_radius,shape_tolerance=shape_tolerance)
      call zone_moments_of_inertia(eta_squared,log_radius,log_mass,shell_mass,zone_start, &
           zone_end,omega,mean_radius,moment_of_inertia,di_domega)
      total_moment_of_inertia = 0.0d0
      total_di_domega = 0.0d0
! Their totals, sum I and sum di_domega.
      do zone_idx = zone_start,zone_end
         total_moment_of_inertia = total_moment_of_inertia + &
              moment_of_inertia(zone_idx)
         total_di_domega = total_di_domega + di_domega(zone_idx)
      end do
      new_angular_momentum = omega_guess*total_moment_of_inertia
! If omega_guess misses J by more than the tolerance, take a Newton step
!   delta_omega = (J - omega I) / (I + omega di_domega)
! with I and di_domega the totals. (di_domega is (omega/2) dI/domega,
! not dI/domega; see zone_moments_of_inertia. For stellar omega the
! second term is negligible and the step is in effect (J - omega I)/I.)
      delta_angular_momentum = total_angular_momentum - new_angular_momentum
      if(.not.disk_locked)then
      if(dabs(delta_angular_momentum/total_angular_momentum).gt. &
           tolerance) then
         if(iteration_count.lt.20) then
            iteration_count = iteration_count + 1
            delta_omega = delta_angular_momentum/ &
                 (total_moment_of_inertia + omega_guess*total_di_domega)
            omega_guess = omega_guess + delta_omega
            cycle omega_iter
         end if
      end if
      end if
      exit omega_iter
      end do omega_iter
! Make j consistent with the final omega and I.
      do zone_idx = zone_start,zone_end
         specific_angular_momentum(zone_idx) = &
              moment_of_inertia(zone_idx)*omega(zone_idx)/shell_mass(zone_idx)
      end do

      return
end subroutine solid_body_omega

! ---------------------------------------------------------------
! enforce_rotation_profile (legacy wczimp; JNT 2025/09 from the 2015/05
! version): imposes the rotation law on the convection zone
! istart..iend, conserving its angular momentum J_cz (in mode 3, the
! star's). With alpha = rotation_settings%cz_rotation_law_exponent:
!   - a solid body (solid_body_omega over the zone) if alpha >= 0,
!     force_solid_body_rotation is set, solid_body_mode_flag = 1, or the
!     zone does not reach the surface (iend < num_zones);
!   - otherwise (the surface convection zone, alpha < 0):
!     solid_body_mode_flag = 2: omega = C r^alpha in the zone, and the
!       interior (zones 1..istart-1) a separate solid body with its own
!       angular momentum;
!     solid_body_mode_flag = 3: omega = C r^alpha in the zone, and the
!       interior a solid body at the rate of the zone's base, with C set
!       by the angular momentum of the whole star;
!     alpha <= -2 and solid_body_mode_flag = 0: constant specific
!       angular momentum, j = J_cz / M_cz;
!     else: omega = C r^alpha, C = J_cz / sum_i r_i^alpha I_i.
! The imposed laws use the I of the previous call for C, and then
! re-solve each zone's omega alone at fixed j (solid_body_omega), as I
! depends on omega; omega follows the law up to that change in I.
! (solid_body_mode_flag is legacy IMPJMOD; legacy wcz lacked modes 2
! and 3.) inertia_tolerance (optional) goes to every solid_body_omega
! call.
subroutine enforce_rotation_profile(log_density, specific_angular_momentum, log_radius, &
     log_mass, shell_mass, istart, iend, eta_squared, moment_of_inertia, &
     omega, qiw, mean_radius, num_zones, disk_gate_age_gyr, inertia_tolerance, &
     shape_tolerance)
      use composition_layout, only: json
      use physical_constants
      use math_lib
      implicit none

      double precision, intent(in) :: log_density(json)
      double precision, intent(inout) :: specific_angular_momentum(json)
      double precision, intent(in) :: log_radius(json), log_mass(json)
      double precision, intent(in) :: shell_mass(json)
      integer, intent(inout) :: istart, iend
      double precision, intent(inout) :: eta_squared(json), &
           moment_of_inertia(json)
      double precision, intent(inout) :: omega(json), qiw(json), &
           mean_radius(json)
      integer, intent(in) :: num_zones
      double precision, intent(in) :: disk_gate_age_gyr
      double precision, intent(in), optional :: inertia_tolerance
      double precision, intent(in), optional :: shape_tolerance
! --- locals ---
! In the power-law branches cz_total_mass accumulates sum r^alpha I,
! the normalization of omega = C r^alpha; only in the constant-j branch
! is it the mass of the zone. power_law_norm is C (or j there).
      double precision :: cz_total_am, cz_total_mass, power_law_norm
      integer :: zone_idx, zone_start, zone_end

! Solid body. alpha >= 0 rather than > 0 (JVS 2026/02/13), so alpha = 0
! takes this branch. A convection zone below the surface always rotates
! as a solid body (JNT 2015/05/07): the imposed law and modes 2 and 3
! act on the surface convection zone only.
      if((rotation_settings%cz_rotation_law_exponent.ge.0.0d0) .or. (rotation_settings%force_solid_body_rotation) .or. &
           (rotation_settings%solid_body_mode_flag.eq.1) .or. (iend.lt.num_zones)) then
         call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
              log_mass,shell_mass,istart,iend,eta_squared, &
              moment_of_inertia,omega,qiw,mean_radius,num_zones, &
              disk_gate_age_gyr, &
              inertia_tolerance=inertia_tolerance, &
                 shape_tolerance=shape_tolerance)
! Mode 2: omega = C r^alpha in the convection zone,
! C = J_cz / sum_i r_i^alpha I_i, and everything below it a solid body
! detached from the convection zone.
      else if(rotation_settings%solid_body_mode_flag.eq.2)then
         cz_total_am = specific_angular_momentum(istart)*shell_mass(istart)
         cz_total_mass = exp(ln10*rotation_settings%cz_rotation_law_exponent*log_radius(istart))* &
              moment_of_inertia(istart)
         do zone_idx = istart+1,iend
            cz_total_am = cz_total_am + &
                 specific_angular_momentum(zone_idx)*shell_mass(zone_idx)
            cz_total_mass = cz_total_mass + &
                 exp(ln10*rotation_settings%cz_rotation_law_exponent*log_radius(zone_idx))* &
                 moment_of_inertia(zone_idx)
         end do

! Impose omega = C r^alpha and j = omega I / dm, then re-solve each
! zone's omega at fixed j.
         power_law_norm = cz_total_am/cz_total_mass
         do zone_idx = istart,iend
            omega(zone_idx) = power_law_norm*exp(ln10*rotation_settings%cz_rotation_law_exponent* &
                 log_radius(zone_idx))
            specific_angular_momentum(zone_idx) = omega(zone_idx)* &
                 moment_of_inertia(zone_idx)/shell_mass(zone_idx)
            zone_start = zone_idx
            zone_end = zone_idx
            call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones, &
                 disk_gate_age_gyr, &
                 inertia_tolerance=inertia_tolerance, &
                 shape_tolerance=shape_tolerance)
         end do
! The interior rotates as one solid body with its own angular momentum.
         if (istart.gt.1) then
            call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,1,istart-1,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones, &
                 disk_gate_age_gyr, &
                 inertia_tolerance=inertia_tolerance, &
                 shape_tolerance=shape_tolerance)
         endif
! Mode 3: the whole star (zones 1..iend) shares its angular momentum
! between the convection zone, omega = C r^alpha, and the interior
! (zones 1..istart), held at the rate of the zone's base,
! C r(istart)^alpha:
!   C = J / (r(istart)^alpha sum_{i<=istart} I_i
!            + sum_{istart<i<=iend} r_i^alpha I_i).
      else if(rotation_settings%solid_body_mode_flag.eq.3) then
         cz_total_am = specific_angular_momentum(1)*shell_mass(1)
         cz_total_mass = exp(ln10*rotation_settings%cz_rotation_law_exponent*log_radius(istart))* &
              moment_of_inertia(1)
         if (istart.gt.1) then
            do zone_idx = 2,istart
               cz_total_am = cz_total_am + &
                    specific_angular_momentum(zone_idx)*shell_mass(zone_idx)
               cz_total_mass = cz_total_mass + &
                    exp(ln10*rotation_settings%cz_rotation_law_exponent*log_radius(istart))* &
                    moment_of_inertia(zone_idx)
            end do
         endif
         do zone_idx = (istart+1),iend
            cz_total_am = cz_total_am + &
                 specific_angular_momentum(zone_idx)*shell_mass(zone_idx)
            cz_total_mass = cz_total_mass + &
                 exp(ln10*rotation_settings%cz_rotation_law_exponent*log_radius(zone_idx))* &
                 moment_of_inertia(zone_idx)
         end do
! Impose omega = C r^alpha and j = omega I / dm in the convection zone,
! then re-solve each zone's omega at fixed j.

         power_law_norm = cz_total_am/cz_total_mass
         do zone_idx = istart,iend
            omega(zone_idx) = power_law_norm*exp(ln10*rotation_settings%cz_rotation_law_exponent* &
                 log_radius(zone_idx))
            specific_angular_momentum(zone_idx) = omega(zone_idx)* &
                 moment_of_inertia(zone_idx)/shell_mass(zone_idx)
            zone_start = zone_idx
            zone_end = zone_idx
            call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones, &
                 disk_gate_age_gyr, &
                 inertia_tolerance=inertia_tolerance, &
                 shape_tolerance=shape_tolerance)
         end do
! Set the interior to the rate of the zone's base and re-solve it as one
! solid body.
         if (istart.gt.1) then
            do zone_idx = 1,(istart-1)
               omega(zone_idx) = omega(istart)
               specific_angular_momentum(zone_idx) = omega(zone_idx)* &
                    moment_of_inertia(zone_idx)/shell_mass(zone_idx)
            end do
            call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,1,istart-1,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones, &
                 disk_gate_age_gyr, &
                 inertia_tolerance=inertia_tolerance, &
                 shape_tolerance=shape_tolerance)
         endif

      else if((rotation_settings%cz_rotation_law_exponent.le.-2.0d0) .and. (rotation_settings%solid_body_mode_flag.eq.0))then
! Constant specific angular momentum in the convection zone: sum the
! zone's angular momentum J_cz and mass M_cz.
         cz_total_am = specific_angular_momentum(istart)*shell_mass(istart)
         cz_total_mass = shell_mass(istart)
         do zone_idx = istart+1,iend
            cz_total_am = cz_total_am + &
                 specific_angular_momentum(zone_idx)*shell_mass(zone_idx)
            cz_total_mass = cz_total_mass + shell_mass(zone_idx)
         end do
! Set j = J_cz / M_cz and omega = j dm / I, then re-solve each zone's
! omega at fixed j.
         power_law_norm = cz_total_am/cz_total_mass
         do zone_idx = istart,iend
            specific_angular_momentum(zone_idx) = power_law_norm
            omega(zone_idx) = power_law_norm*shell_mass(zone_idx)/ &
                 moment_of_inertia(zone_idx)
            zone_start = zone_idx
            zone_end = zone_idx
            call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones, &
                 disk_gate_age_gyr, &
                 inertia_tolerance=inertia_tolerance, &
                 shape_tolerance=shape_tolerance)
         end do
      else
! The general law omega = C r^alpha, C constant over the convection
! zone. Conserving the zone's angular momentum,
!   C = J_cz / sum_i r_i^alpha I_i,
! and then j_i = C r_i^alpha I_i / dm_i.
         cz_total_am = specific_angular_momentum(istart)*shell_mass(istart)
         cz_total_mass = exp(ln10*rotation_settings%cz_rotation_law_exponent*log_radius(istart))* &
              moment_of_inertia(istart)
         do zone_idx = istart+1,iend
            cz_total_am = cz_total_am + &
                 specific_angular_momentum(zone_idx)*shell_mass(zone_idx)
            cz_total_mass = cz_total_mass + &
                 exp(ln10*rotation_settings%cz_rotation_law_exponent*log_radius(zone_idx))* &
                 moment_of_inertia(zone_idx)
         end do
! Impose omega = C r^alpha and j = omega I / dm, then re-solve each
! zone's omega at fixed j.
         power_law_norm = cz_total_am/cz_total_mass
         do zone_idx = istart,iend
            omega(zone_idx) = power_law_norm*exp(ln10*rotation_settings%cz_rotation_law_exponent* &
                 log_radius(zone_idx))
            specific_angular_momentum(zone_idx) = omega(zone_idx)* &
                 moment_of_inertia(zone_idx)/shell_mass(zone_idx)
            zone_start = zone_idx
            zone_end = zone_idx
            call solid_body_omega(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,zone_start,zone_end,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones, &
                 disk_gate_age_gyr, &
                 inertia_tolerance=inertia_tolerance, &
                 shape_tolerance=shape_tolerance)
         end do
      endif
      return
end subroutine enforce_rotation_profile

end module omega_profile
