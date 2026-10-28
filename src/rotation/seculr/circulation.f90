!----------------------------------------------------------------------
! circulation
!----------------------------------------------------------------------
! The velocities that drive rotational mixing: the circulation and
! instability velocities (circulation_velocities), the diffusion
! coefficients formed from them (diffusion_velocity_scales), the Zahn
! (1992) coupling factor (zahn_coupling_factor) and the quadrupole
! moment of the distorted star (compute_quadrupole). Each mechanism's
! velocity is a star-blind kernel in instability_velocities, which
! circulation_velocities calls with the structure at the interfaces.
module circulation
      use rotation_configuration, only: rotation_settings
      use instability_velocities, only: compute_eddington_sweet_velocity, &
           inhibit_by_mu_gradient, compute_gsf_velocity, compute_shear_velocities
      use transport_mechanism_layout, only: i_mech_eddington_sweet, i_mech_gsf, &
           i_mech_shear, n_transport_mechanisms
      implicit none
      private
      public :: circulation_velocities, diffusion_velocity_scales, zahn_coupling_factor, compute_quadrupole

contains

! ---------------------------------------------------------------
! The circulation and instability velocities at each interface
! zone_min..zone_max of the transport region (interface i lies between
! zones i-1 and i). Each stability criterion is a structural coefficient
! (rotation_stability_setup, stability_coefficients) times a gradient
! across the interface, compared with a critical value; where the
! interface is unstable the mechanism gets a velocity:
!  - es_circulation_velocity, the Eddington-Sweet circulation (Zahn
!    1991, 1992), corrected for the quadrupole moment and reduced by the
!    mu gradient (Maeder & Zahn 1998);
!  - gsf_circulation_velocity, the GSF instability (Kippenhahn
!    1980), optionally inhibited by viscosity and mu gradients;
!  - secular_shear_velocity, the secular shear instability, or the
!    dynamical one where Ri < 1/4 (Endal & Sofia 1978).
! The Eddington-Sweet and GSF velocities are skipped at interfaces with
! convective zones on both sides. From the second iteration on, each
! velocity is averaged with the previous iteration's. Their sum is
! total_velocity (legacy HV); any_transport_active is set if
! it is positive anywhere. With no_am_transport_in_core the velocities
! are zeroed inside the energy-generating core (zero_transport_in_core).
! domega_dlnr is d omega/d ln r and dynamical_shear_omega_limit its
! dynamical-shear limit.
! Each velocity is later multiplied by r (diffusion_velocity_scales) to
! give a diffusion coefficient.
! The formulas are the kernels of instability_velocities. Its internal steps
! keep the books: they interpolate the physics to the interfaces
! (stab), and store, clear, scale, average and sum the velocities
! (the arguments and circ_scr).
subroutine circulation_velocities(stab, sub_physics, circ_scr, &
     log_radius, radius, zone_min, zone_max, iteration, &
     am_transport_convective_flag, num_zones, omega, any_transport_active, &
     domega_dlnr, dynamical_shear_omega_limit, total_velocity, &
     total_luminosity, log_pressure, mechanism_velocity)
      use rotation_scratch_lib, only: stability_state, &
           substep_physics_state, circulation_velocity_state

      use star_info_lib, only: json
      use math_lib
      implicit none
      type(stability_state), intent(inout) :: stab
      type(substep_physics_state), intent(in) :: sub_physics
      type(circulation_velocity_state), intent(inout) :: circ_scr

      double precision, intent(in) :: log_radius(json), radius(json)
      integer, intent(in) :: zone_min, zone_max, iteration
      logical, intent(in) :: am_transport_convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(in) :: omega(json)
      logical, intent(out) :: any_transport_active
      double precision, intent(in) :: domega_dlnr(json), &
           dynamical_shear_omega_limit(json)
      double precision, intent(out) :: total_velocity(json)
      double precision, intent(in) :: total_luminosity, log_pressure(json)
! Each mechanism's velocity (cm/s; transport_mechanism_layout's slots;
! the star's, for the output): the previous iteration's on entry, this
! iteration's on return.
      double precision, intent(inout) :: mechanism_velocity(n_transport_mechanisms, json)
! --- locals ---
      integer :: i, i0, i1, m
      double precision :: dell

      any_transport_active = .false.
      if (iteration.gt.1) call store_previous_velocities()
      call clear_velocities()
      call interpolate_to_interfaces()
      call compute_eddington_sweet_velocity(zone_min, zone_max, am_transport_convective_flag, &
           radius, omega, stab%omega_interface, stab%radius_interface, &
           stab%local_gravity, stab%quadrupole_moment, stab%es_velocity_coeff1, &
           stab%es_velocity_coeff2, stab%es_shear_coeff, &
           stab%circulation_correction_ratio, mechanism_velocity(i_mech_eddington_sweet,:))
      call inhibit_by_mu_gradient(zone_min, zone_max, sub_physics%mu, log_pressure, &
           stab%del_grad_diff_interface, mechanism_velocity(i_mech_eddington_sweet,:))
      call compute_gsf_velocity(zone_min, zone_max, rotation_settings%gsf_inhibition_mode, &
           am_transport_convective_flag, radius, omega, domega_dlnr, dynamical_shear_omega_limit, &
           sub_physics%mu, stab%mean_molecular_weight_interface, stab%omega_interface, &
           stab%kinematic_viscosity_interface, stab%thermal_diffusivity_interface, &
           stab%gravity_factor_interface, stab%gsf_kippenhahn_coeff, &
           mechanism_velocity(i_mech_gsf,:))
      call compute_shear_velocities(zone_min, zone_max, rotation_settings%critical_reynolds, &
           log_radius, radius, omega, domega_dlnr, dynamical_shear_omega_limit, &
           sub_physics%mu, stab%mean_molecular_weight_interface, stab%radius_interface, &
           stab%kinematic_viscosity_interface, stab%thermal_diffusivity_interface, &
           stab%mu_gradient_richardson_coeff, mechanism_velocity(i_mech_shear,:), &
           mechanism_velocity(i_mech_gsf,:))
      call apply_velocity_scales()
      if (iteration.gt.1) call average_with_previous()
      call sum_velocities()
      if (rotation_settings%no_am_transport_in_core) call zero_transport_in_core()

contains

! ---------------------------------------------------------------
! Keep the previous iteration's velocities for average_with_previous.
subroutine store_previous_velocities()
      do i = zone_min,zone_max
         do m = 1, n_transport_mechanisms
            circ_scr%velocity_prev(m,i) = mechanism_velocity(m,i)
         end do
      end do
! (commented out: the advection-diffusion transport's velocities)
!      if (rotation_settings%use_diffusion_advection_transport) then
!         do i = zone_min,zone_max
!            rot_scr%es_advective_velocity_prev(i) = rot_scr%es_advective_velocity(i)
!            rot_scr%es_diffusive_velocity_prev(i) = rot_scr%es_diffusive_velocity(i)
!         end do
!      end if
end subroutine store_previous_velocities

! ---------------------------------------------------------------
! Every velocity starts at zero, so an interface a mechanism skips (a
! pair of convective zones, a stable interface) keeps none.
subroutine clear_velocities()
      do i = 1,num_zones
         total_velocity(i) = 0.0d0
         do m = 1, n_transport_mechanisms
            mechanism_velocity(m,i) = 0.0d0
         end do
      end do
end subroutine clear_velocities

! ---------------------------------------------------------------
! The mean molecular weight (and what depends on it) and omega change
! during a diffusion sub-step, so their values at the interfaces, where
! stability is evaluated, are found again each iteration. The kinematic
! viscosity and thermal diffusivity use 4-point Lagrange interpolation
! in ln (log_lagrange4) from zones i-2..i+1.
subroutine interpolate_to_interfaces()
! The first interface (2) uses zones 1-4.
      if (zone_min.lt.3) then
         stab%kinematic_viscosity_interface(2)=log_lagrange4(sub_physics%visc, 1, 2)
         stab%thermal_diffusivity_interface(2)=log_lagrange4(sub_physics%thdif, 1, 2)
         i0=3
      else
         i0=zone_min
      end if
! The last interface (num_zones) uses the last four zones.
      if (zone_max.eq.num_zones) then
         stab%kinematic_viscosity_interface(num_zones)= &
              log_lagrange4(sub_physics%visc, num_zones-3, num_zones)
         stab%thermal_diffusivity_interface(num_zones)= &
              log_lagrange4(sub_physics%thdif, num_zones-3, num_zones)
         i1=num_zones-1
      else
         i1=zone_max
      end if
! The other interfaces.
      do i = i0,i1
         stab%kinematic_viscosity_interface(i)=log_lagrange4(sub_physics%visc, i-2, i)
         stab%thermal_diffusivity_interface(i)=log_lagrange4(sub_physics%thdif, i-2, i)
      end do
! mu and omega: the mean of the two neighbouring zones.
      do i = 2,num_zones
         stab%mean_molecular_weight_interface(i) = 0.5d0*(sub_physics%mu(i)+sub_physics%mu(i-1))
         stab%omega_interface(i) = 0.5d0*(omega(i)+omega(i-1))
      end do
end subroutine interpolate_to_interfaces

! ---------------------------------------------------------------
! Multiply each mechanism's velocity by its efficiency factor
! (rotation_settings%velocity_scale: the controls es_velocity_scale,
! gsf_velocity_scale and secular_shear_velocity_scale, legacy FES, FGSF,
! FSS). The mu-gradient inhibition is already applied. A velocity is a
! magnitude, so the product is taken in absolute value (it always was for
! U_ES; the GSF and shear velocities are positive, and every shipped deck
! sets the scales to 1).
subroutine apply_velocity_scales()
      do i = zone_min,zone_max
         do m = 1, n_transport_mechanisms
            mechanism_velocity(m,i) = abs(rotation_settings%velocity_scale(m)* &
                 mechanism_velocity(m,i))
         end do
      end do
! (commented out: the scales of the advection-diffusion transport's
! parts, and their GSF and shear diffusion coefficients r |U|)
! The same for the advection-diffusion parts (MHP 11/94).
!      if (rotation_settings%use_diffusion_advection_transport) then
!         do i = zone_min,zone_max
!               rot_scr%es_diffusive_velocity(i)=abs(rotation_settings%es_velocity_scale* &
!                    rot_scr%es_diffusive_velocity(i))
!               rot_scr%es_advective_velocity(i)=rotation_settings%es_velocity_scale* &
!                    rot_scr%es_advective_velocity(i)
!         end do
!      end if
!      if (rotation_settings%use_diffusion_advection_transport) then
! In the advection-diffusion transport (MHP 5/02) the GSF and shear
! instabilities act as diffusion. In the Endal & Sofia form a velocity
! U (e.g. U_GSF ~ K dlnwdr omega^2/r) enters d/dr(con U r d omega/dr),
! so its diffusion coefficient is r U. In the Zahn form the equation is
! d/dr(con U omega), and a velocity proportional to d omega/dr is in
! effect diffusion with coefficient U omega/(d omega/dr), which is how
! es_diffusive_velocity is formed. For the GSF and shear the code uses
! the Endal & Sofia coefficient: r |U| of each is added to
! es_diffusive_velocity and kept in shear_diffusion_coeff and
! gsf_diffusion_coeff.
!         do i = zone_min,zone_max
!            rot_scr%es_diffusive_velocity(i) = rot_scr%es_diffusive_velocity(i)+ &
!                 stab%radius_interface(i)*(abs(gsf_circulation_velocity(i))+ &
!                 abs(secular_shear_velocity(i)))
!            rot_scr%shear_diffusion_coeff(i) = stab%radius_interface(i)* &
!                 abs(secular_shear_velocity(i))
!            rot_scr%gsf_diffusion_coeff(i) = stab%radius_interface(i)* &
!                 abs(gsf_circulation_velocity(i))
!         end do
!      end if
end subroutine apply_velocity_scales

! ---------------------------------------------------------------
! From the second iteration on, the mean of the new velocities and the
! previous iteration's (MHP 11/94).
subroutine average_with_previous()
      do i = zone_min,zone_max
         do m = 1, n_transport_mechanisms
            mechanism_velocity(m,i) = 0.5d0*(mechanism_velocity(m,i) &
                 + circ_scr%velocity_prev(m,i))
         end do
      end do
! (commented out: the same for the advection-diffusion parts)
!      if (rotation_settings%use_diffusion_advection_transport) then
!         do i = zone_min,zone_max
!            rot_scr%es_advective_velocity(i) = 0.5d0*(rot_scr%es_advective_velocity(i) &
!                 + rot_scr%es_advective_velocity_prev(i))
!            rot_scr%es_diffusive_velocity(i) = 0.5d0*(rot_scr%es_diffusive_velocity(i) &
!                 + rot_scr%es_diffusive_velocity_prev(i))
!         end do
!      end if
end subroutine average_with_previous

! ---------------------------------------------------------------
! total_velocity, the sum of the mechanisms' velocities at
! each interface (mechanism_sum), zeroed below 1e-20 cm/s;
! any_transport_active if any is positive.
subroutine sum_velocities()
      do i =zone_min,zone_max
         total_velocity(i) = mechanism_sum(mechanism_velocity(:,i))
         if (total_velocity(i).lt.1.0d-20) &
              total_velocity(i)=0.0d0
         if (total_velocity(i).gt.0.0d0) any_transport_active=.true.
      end do
end subroutine sum_velocities

! ---------------------------------------------------------------
! With no_am_transport_in_core (mixing without angular momentum
! transport; 9/93): zero every velocity from zone_min outward up to the
! first interface where L reaches 99% of the total luminosity, to avoid
! numerical problems in the hydrogen-burning shell.
subroutine zero_transport_in_core()
      do i = zone_min,zone_max
         dell = stab%luminosity_interface(i)/total_luminosity
         if (dell.lt.9.9d-1) then
            total_velocity(i) = 0.0d0
            do m = 1, n_transport_mechanisms
               mechanism_velocity(m,i) = 0.0d0
            end do
!            rot_scr%es_advective_velocity(i) = 0.0d0
!            rot_scr%es_diffusive_velocity(i) = 0.0d0
         else
            exit
         end if
      end do
end subroutine zero_transport_in_core


! 4-point Lagrange interpolation of a(k..k+3) to interface i in ln a,
! with the weights stab%lagrange_interp_weights(1:4,i) set by
! rotation_stability_setup (stability_coefficients).
      function log_lagrange4(a, k, i)
      double precision, intent(in) :: a(*)
      integer, intent(in) :: k, i
      double precision :: log_lagrange4
      log_lagrange4 = exp(log(a(k))* &
           stab%lagrange_interp_weights(1,i)+log(a(k+1))* &
           stab%lagrange_interp_weights(2,i)+ &
           log(a(k+2))*stab%lagrange_interp_weights(3,i)+log(a(k+3))* &
           stab%lagrange_interp_weights(4,i))
      end function log_lagrange4

end subroutine circulation_velocities

! ---------------------------------------------------------------
! The sum of the mechanisms' values v(1:n), in slot order:
!   ((v(1) + v(2)) + v(3)) + ...
! The order is part of the results (floating-point addition is not
! associative): a new mechanism goes in the next slot
! (transport_mechanism_layout).
pure double precision function mechanism_sum(v)
      double precision, intent(in) :: v(:)
      integer :: m
      mechanism_sum = v(1)
      do m = 2, size(v)
         mechanism_sum = mechanism_sum + v(m)
      end do
end function mechanism_sum

! ---------------------------------------------------------------
! The diffusion coefficients (velocity times the length scale r) of the
! angular momentum and composition equations, at the interfaces
! 2..num_zones; element i is the interface between zones i-1 and i.
! Inputs: radius, the zone radii (cm), and each mechanism's velocity
! (the star's mechanism_velocity) and the Zahn coupling factor vfc.
! Outputs:
!   radius_mid(i) = (radius(i) + radius(i-1))/2, the interface radius r
!     (radius_mid(1) = 0);
!   am_diffusion_coeff (legacy COD)
!     = 4 pi difad_velocity_scale (sum of the mechanisms' U) r;
!   mixing_diffusion_coeff (legacy COD2)
!     = 4 pi difad_velocity_scale mixing_velocity_scale
!       (sum of each mechanism's mixing_scale times its U) r,
!     each mechanism with its own mixing efficiency (2/91), or, with
!     use_zahn_coupling_factor, am_diffusion_coeff times vfc.
! difad_velocity_scale and mixing_velocity_scale are legacy FW and FC.
! With use_constant_background_diffusion, 4 pi difad_velocity_scale
! constant_background_diffusion_coeff is added to am_diffusion_coeff
! only (MHP 9/14, 8/13). With no_am_transport_in_core,
! am_diffusion_coeff is zero everywhere and only the mixing remains.
subroutine diffusion_velocity_scales(radius, num_zones, radius_mid, &
     am_diffusion_coeff, mixing_diffusion_coeff, mechanism_velocity, vfc)
      use star_info_lib, only: json
      use physical_constants
      implicit none
      double precision, intent(in) :: radius(json)
      integer, intent(in) :: num_zones
      double precision, intent(out) :: radius_mid(json)
      double precision, intent(out) :: am_diffusion_coeff(json), &
           mixing_diffusion_coeff(json)
! the velocities of circulation_velocities (cm/s) and the coupling
! factor of zahn_coupling_factor
      double precision, intent(in) :: mechanism_velocity(n_transport_mechanisms, json), &
           vfc(json)
      integer :: i
      double precision :: con1, con2

      radius_mid(1) = 0.0d0
      do i = 2,num_zones
         radius_mid(i) = 0.5d0*(radius(i)+radius(i-1))
      end do
      con1 = c4pi*rotation_settings%difad_velocity_scale
      if (.not.rotation_settings%use_constant_background_diffusion) then
         if (.not.rotation_settings%use_zahn_coupling_factor) then
            con2 = con1*rotation_settings%mixing_velocity_scale
            do i = 2,num_zones
               am_diffusion_coeff(i) = con1*mechanism_sum(mechanism_velocity(:,i))* &
                    radius_mid(i)
               mixing_diffusion_coeff(i) = con2*mechanism_sum(rotation_settings%mixing_scale* &
                    mechanism_velocity(:,i))*radius_mid(i)
            end do
         else
            do i = 2,num_zones
               am_diffusion_coeff(i) = con1*mechanism_sum(mechanism_velocity(:,i))* &
                    radius_mid(i)
               mixing_diffusion_coeff(i) = am_diffusion_coeff(i)*vfc(i)
            end do
         end if
      else
         if (.not.rotation_settings%use_zahn_coupling_factor) then
            con2 = con1*rotation_settings%mixing_velocity_scale
            do i = 2,num_zones
               am_diffusion_coeff(i) = con1*(rotation_settings%constant_background_diffusion_coeff+ &
                    mechanism_sum(mechanism_velocity(:,i))*radius_mid(i))
               mixing_diffusion_coeff(i) = con2*mechanism_sum(rotation_settings%mixing_scale* &
                    mechanism_velocity(:,i))*radius_mid(i)
            end do
         else
            do i = 2,num_zones
               am_diffusion_coeff(i) = con1*mechanism_sum(mechanism_velocity(:,i))* &
                    radius_mid(i)
               mixing_diffusion_coeff(i) = am_diffusion_coeff(i)*vfc(i)
! The background coefficient is added after the mixing coefficient is
! formed, so it does not mix (MHP 8/13).
               am_diffusion_coeff(i) = am_diffusion_coeff(i)+ &
                    con1*rotation_settings%constant_background_diffusion_coeff
            end do
         end if
      end if
! Suppress angular momentum transport but keep the mixing.
      if (rotation_settings%no_am_transport_in_core) then
         do i = 2,num_zones
            am_diffusion_coeff(i) = 0.0d0
         end do
      end if
      return
end subroutine diffusion_velocity_scales

! ---------------------------------------------------------------
! The coupling factor vfc between the angular momentum and the
! composition diffusion coefficients, with the prescription of Zahn
! (1992), at the interfaces zone_min..zone_max:
!   vfc = min(1, f_c U/|2V - alpha U|),
! where U is the total circulation velocity (total_velocity),
!   V = (1/(6 rho r)) d(rho r^2 U)/dr,
!   alpha = (1/2) d ln(r^2 omega)/d ln r,
! and f_c = mixing_velocity_scale (legacy FC) plays the role of C/30.
! vfc = 0 where U <= 0, and everywhere if f_c <= 0 (no mixing).
subroutine zahn_coupling_factor(stab, log_density, radius, total_velocity, zone_min, &
     zone_max, angular_velocity, vfc)
      use rotation_scratch_lib, only: stability_state

      use star_info_lib, only: json
      use physical_constants
      use math_lib
      implicit none
      type(stability_state), intent(in) :: stab

! log_density: log10 rho (g/cm^3) at the zone centres.
! radius: the unlogged radius (cm) at the zone centres.
! total_velocity: the total velocity U at the interfaces (zone edges).
! zone_min, zone_max: the first and last interface.
! angular_velocity: omega at the zone centres.
      double precision, intent(in) :: log_density(json), radius(json), &
           total_velocity(json)
      integer, intent(in) :: zone_min, zone_max
      double precision, intent(in) :: angular_velocity(json)
! vfc: the coupling factor, set at the interfaces zone_min..zone_max
      double precision, intent(inout) :: vfc(json)

      double precision :: half_dlnj_dlnr(json), circ_velocity(json)
      integer :: zone_index
      double precision :: omega_mid
      double precision :: interp_weight
      double precision :: v_minus, v_plus
      double precision :: density_r2_v_minus, density_r2_v_plus
      double precision :: density_radius, d_density_r2, du_dr
      double precision :: denom_test

! f_c <= 0: no mixing.

      if (rotation_settings%mixing_velocity_scale.le.0.0D0) then
         do zone_index = zone_min, zone_max
            vfc(zone_index) = 0.0D0
         end do
         return
      endif
! alpha = (1/2) d ln(r^2 omega)/d ln r
!       = (1/2) d(omega r^2)/dr/(omega r)
! at each interface, with omega the mean of the two zones. (The F77 code
! lacked the 1/2; corrected in 2026.)
      do zone_index = zone_min, zone_max
         omega_mid = 0.5D0*(angular_velocity(zone_index)+ &
              angular_velocity(zone_index-1))
         half_dlnj_dlnr(zone_index) = 0.5D0*(angular_velocity(zone_index)* &
              radius(zone_index)**2-angular_velocity(zone_index-1)* &
              radius(zone_index-1)**2)/ &
              (radius(zone_index)-radius(zone_index-1))/ &
              (omega_mid*stab%radius_interface(zone_index))
      end do
! 2V = (1/(3 rho r)) d(rho r^2 U)/dr at interface i, from rho r^2 U at
! the zone centres i-1 and i, with U at a zone centre interpolated
! linearly between the interfaces on either side; rho and r at the
! interface are stab's.
      interp_weight = (radius(zone_min)-stab%radius_interface(zone_min))/ &
           (stab%radius_interface(zone_min+1)-stab%radius_interface(zone_min))
      v_minus = total_velocity(zone_min)+interp_weight* &
           (total_velocity(zone_min+1)-total_velocity(zone_min))
      density_r2_v_minus = exp(ln10*log_density(zone_min))* &
           radius(zone_min)**2*v_minus
      do zone_index = zone_min + 1, zone_max - 1
         interp_weight = (radius(zone_index)-stab%radius_interface(zone_index))/ &
              (stab%radius_interface(zone_index+1)-stab%radius_interface(zone_index))
         v_plus = total_velocity(zone_index)+interp_weight* &
              (total_velocity(zone_index+1)-total_velocity(zone_index))
         density_r2_v_plus = exp(ln10*log_density(zone_index))* &
              radius(zone_index)**2*v_plus
         density_radius = stab%density_interface(zone_index)*stab%radius_interface(zone_index)
         circ_velocity(zone_index) = cc13*(density_r2_v_plus-density_r2_v_minus)/ &
              density_radius/(radius(zone_index)-radius(zone_index-1))
         density_r2_v_minus = density_r2_v_plus
      end do
! At the first and last interfaces,
!   2V = (1/3) [U d(rho r^2)/dr/(rho r) + r dU/dr],
! with d(rho r^2)/dr across the interface and dU/dr the one-sided
! difference to the neighbouring interface.
      density_radius = stab%density_interface(zone_min)*stab%radius_interface(zone_min)
      d_density_r2 = (exp(ln10*log_density(zone_min))*radius(zone_min)**2 - &
           exp(ln10*log_density(zone_min-1))*radius(zone_min-1)**2)/ &
           (radius(zone_min) - radius(zone_min-1))
      du_dr = (total_velocity(zone_min+1)-total_velocity(zone_min))/ &
           (stab%radius_interface(zone_min+1)-stab%radius_interface(zone_min))
      circ_velocity(zone_min) = cc13*(total_velocity(zone_min)*d_density_r2/ &
           density_radius + stab%radius_interface(zone_min)*du_dr)
      density_radius = stab%density_interface(zone_max)*stab%radius_interface(zone_max)
      d_density_r2 = (exp(ln10*log_density(zone_max))*radius(zone_max)**2 - &
           exp(ln10*log_density(zone_max-1))*radius(zone_max-1)**2)/ &
           (radius(zone_max) - radius(zone_max-1))
      du_dr = (total_velocity(zone_max)-total_velocity(zone_max-1))/ &
           (stab%radius_interface(zone_max)-stab%radius_interface(zone_max-1))
      circ_velocity(zone_max) = cc13*(total_velocity(zone_max)*d_density_r2/ &
           density_radius + stab%radius_interface(zone_max)*du_dr)
! vfc, the same for all mechanisms and at most 1:
!   vfc = f_c/max(|2V - alpha U|/U, f_c).
      do zone_index = zone_min, zone_max
         if (total_velocity(zone_index).le.0.0D0) then
            vfc(zone_index) = 0.0D0
         else
            denom_test = max(abs(circ_velocity(zone_index)- &
                 half_dlnj_dlnr(zone_index)*total_velocity(zone_index))/ &
                 total_velocity(zone_index),rotation_settings%mixing_velocity_scale)
            vfc(zone_index) = rotation_settings%mixing_velocity_scale/denom_test
         endif
      end do

      return
end subroutine zahn_coupling_factor

! ---------------------------------------------------------------
! The quadrupole moment Phi of the rotating star at the zone centres
! (stab%quadrupole_moment, legacy QUAD), with the Sweet (1950)
! formulation as given in Zahn (1992). Phi solves
!   (1/r) d^2(r Phi)/dr^2 - 6 Phi/r^2 + 4 pi G rho (d rho/dP) Phi
!     = -(4/3) pi G (r^2/g) d(rho omega^2)/dr,
! with g the local mean gravity (gravity, ~G m/r^2). With dP/dr =
! -rho g this is
!   Phi'' + (2/r) Phi' - 6 Phi/r^2 - 4 pi G (d rho/dr/g) Phi
!     = -(4/3) pi G (r^2/g) d(rho omega^2)/dr,
! which the code differences on the zone radii into the tridiagonal
! system
!   sub_diag Phi(i-1) + diag Phi(i) + super_diag Phi(i+1) = rhs.
! Boundary conditions: Phi = 0 at the centre; at the surface no mass
! lies outside, so Phi varies as 1/r^3 (the exterior l = 2 solution).
! gravity is copied to stab%local_gravity for circulation_velocities.
! ierr = 1 if a pivot of the elimination is zero.
subroutine compute_quadrupole(stab, &
     log_density, gravity, radius, angular_velocity, num_zones, ierr)
      use rotation_scratch_lib, only: stability_state
      use star_info_lib, only: json
      use physical_constants
      use math_lib
      implicit none
      type(stability_state), intent(inout) :: stab

      double precision, intent(in) :: log_density(json), gravity(json), &
           radius(json), angular_velocity(json)
      integer, intent(in) :: num_zones
      integer, intent(out) :: ierr

! The tridiagonal system and the work arrays of its solution (Thomas
! algorithm).
      double precision :: sub_diag(json), diag(json), super_diag(json), &
           rhs(json), solution(json), gamma_elim(json)

      double precision :: density_omega2(json)
      integer :: zone_index, matrix_row
      double precision :: four_pi_g
      double precision :: dr, weight_plus, weight_minus, inv_dr2, dr_inv_r
      double precision :: drho_dr
      double precision :: dr_plus, dr_minus, d_density_omega2
      double precision :: radius_plus, outer_boundary_ratio
      double precision :: pivot

      ierr = 0
      do zone_index = 1,num_zones
         density_omega2(zone_index) = exp(ln10*log_density(zone_index))* &
              angular_velocity(zone_index)**2
      end do
      four_pi_g = exp(ln10*(c4pil+cgl))
! The derivatives are centred, with dr the smaller of the left and right
! intervals dr_minus, dr_plus and weights w_plus = dr/dr_plus,
! w_minus = dr/dr_minus:
!   Phi''  = [(Phi(i+1) - Phi(i))/dr_plus - (Phi(i) - Phi(i-1))/dr_minus]/dr,
!   Phi'   = [(Phi(i+1) - Phi(i))/dr_plus + (Phi(i) - Phi(i-1))/dr_minus]/2,
! and d rho/dr and d(rho omega^2)/dr likewise.
! Zone 1: the inner neighbour is the centre, r = 0, where Phi = 0, so it
! drops out (sub_diag(1) = 0). There d rho/dr = (rho(2) - rho(1))/dr and
! d(rho omega^2)/dr is taken as omega^2 d rho/dr.
      dr = min(radius(1),radius(2)-radius(1))
      weight_plus = dr/(radius(2)-radius(1))
      weight_minus = dr/radius(1)
      inv_dr2 = 1.0D0/dr**2
      dr_inv_r = 1.0D0/dr/radius(1)
      sub_diag(1) = 0.0D0
      drho_dr = (exp(ln10*log_density(2))-exp(ln10*log_density(1)))/dr
      diag(1) = -inv_dr2*(weight_plus+weight_minus)+dr_inv_r* &
           (weight_minus-weight_plus)-6.0D0/radius(1)**2 &
           -four_pi_g*drho_dr/gravity(1)
      super_diag(1) = weight_plus*(inv_dr2+dr_inv_r)
      rhs(1) = -four_pi_g*drho_dr*(radius(1)*angular_velocity(1))**2/3.0D0/ &
           gravity(1)
! The interior zones. The Phi'' difference divides by min(dr_minus,
! dr_plus) rather than their mean, so it is exact only on a uniform grid.
      do zone_index = 2,num_zones-1
         dr_plus = radius(zone_index+1)-radius(zone_index)
         dr_minus = radius(zone_index)-radius(zone_index-1)
         dr = min(dr_minus,dr_plus)
         weight_plus = dr/dr_plus
         weight_minus = dr/dr_minus
         inv_dr2 = 1.0D0/dr**2
         dr_inv_r = 1.0D0/dr/radius(zone_index)
         drho_dr = 0.5D0*(weight_plus*exp(ln10*log_density(zone_index+1))+ &
              (weight_minus-weight_plus)* &
              exp(ln10*log_density(zone_index))-weight_minus* &
              exp(ln10*log_density(zone_index-1)))/dr
         sub_diag(zone_index) = weight_minus*(inv_dr2-dr_inv_r)
         diag(zone_index) = -inv_dr2*(weight_plus+weight_minus)+dr_inv_r* &
              (weight_minus-weight_plus)-6.0D0/radius(zone_index)**2 &
              -four_pi_g*drho_dr/gravity(zone_index)
         super_diag(zone_index) = weight_plus*(inv_dr2+dr_inv_r)
         d_density_omega2 = 0.5D0*(weight_plus*density_omega2(zone_index+1) + &
              (weight_minus-weight_plus)*density_omega2(zone_index) - &
              weight_minus*density_omega2(zone_index-1))/dr
         rhs(zone_index) = -d_density_omega2*four_pi_g* &
              radius(zone_index)**2/3.0D0/gravity(zone_index)
      end do
! The surface: Phi' + 3 Phi/r = 0, and outside the star Phi(r) =
! Phi(R) (R/r)^3. Phi'' uses a ghost point at R + dr with that value;
! (2/r) Phi' - 6 Phi/r^2 becomes -12 Phi/r^2; the density and forcing
! terms are dropped (rhs = 0).
      dr = radius(num_zones) - radius(num_zones-1)
      radius_plus = radius(num_zones)+dr
      outer_boundary_ratio = (radius(num_zones)/radius_plus)**3
      inv_dr2 = 1.0D0/dr**2
      sub_diag(num_zones) = inv_dr2
      diag(num_zones) = (outer_boundary_ratio-2.0D0)*inv_dr2-1.2D1/ &
           radius(num_zones)**2
      super_diag(num_zones) = 0.0D0
      rhs(num_zones) = 0.0D0
! Solve by forward elimination and back substitution (the Thomas
! algorithm, tridag of Numerical Recipes).
      if (diag(1).eq.0.0D0) then
         write(*,*) 'compute_quadrupole: singular first pivot'
         ierr = 1
         return
      end if
      pivot = diag(1)
      solution(1) = rhs(1)/pivot
      do matrix_row = 2,num_zones
         gamma_elim(matrix_row) = super_diag(matrix_row-1)/pivot
         pivot = diag(matrix_row) - sub_diag(matrix_row)*gamma_elim(matrix_row)
         if (pivot.eq.0.0D0) then
            write(*,*) 'compute_quadrupole: singular pivot in elimination'
            ierr = 1
            return
         end if
         solution(matrix_row) = (rhs(matrix_row) - &
              sub_diag(matrix_row)*solution(matrix_row-1))/pivot
      end do
! Back substitution.
      do matrix_row = num_zones-1,1,-1
         solution(matrix_row) = solution(matrix_row) - &
              gamma_elim(matrix_row+1)*solution(matrix_row+1)
      end do
      do zone_index = 1,num_zones
         stab%quadrupole_moment(zone_index) = solution(zone_index)
         stab%local_gravity(zone_index) = gravity(zone_index)
      end do

      return
end subroutine compute_quadrupole

end module circulation
