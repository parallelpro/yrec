!----------------------------------------------------------------------
! stability_coefficients
!----------------------------------------------------------------------
! The structure at the zone interfaces, and the parts of the
! circulation and instability velocities that depend only on it.
! secular_transport calls rotation_stability_setup once per sub-step,
! before its iterations; circulation then multiplies these coefficients
! by the current omega and its gradient.
!
! Interface i lies between zones i-1 and i, at the mean of their radii.
! Each quantity there is a 4-point Lagrange interpolation in radius from
! zones i-2..i+1 (the first and last interfaces use the first or last
! four zones). P, T and rho are interpolated in log10 and unlogged; the
! gradients, d ln rho/d ln T, the enclosed mass, eps, L and the local
! gravity are interpolated linearly. The actual gradient at the
! interface is min(grad_rad, grad_ad).
!
! With H_P = P r^2/(rho G m) the pressure scale height and
! Delta = max(grad_ad - grad, floor), the coefficients are:
!   - The Eddington-Sweet circulation (Zahn 1992).
!     circulation forms
!       U_ES = omega^2 (es_velocity_coeff1 + omega^2 es_velocity_coeff2)
!            = (8/3) (omega^2 r^3/(G m)) P/(g Delta c_P rho T)
!              [ (L/m) (1 - omega^2/(2 pi G rho)) - eps ],
!     then corrects it for the quadrupole and subtracts
!     es_shear_coeff * omega d omega/dr, with
!       es_shear_coeff = (L/m) P/(g Delta c_P rho T) r/(3 pi G rho).
!   - The GSF instability (Kippenhahn 1980):
!     gsf_kippenhahn_coeff = (32/3) H_P r^3/(G m Delta). circulation
!     multiplies it by the thermal diffusivity, omega^2/r and the
!     gradient of omega or of the specific angular momentum.
!   - The dynamical shear instability: an interface is stable while
!       Ri = (rho/P) (grad_ad - grad) g^2/(d omega/d ln r)^2 > 1/4,
!     so stab%dynamical_shear_omega_limit = 2 g sqrt((rho/P)
!     max(grad_ad - grad, 1e-3)) is the largest stable |d omega/d ln r|.
!   - The mu-gradient Richardson term: the same criterion with the mu
!     stratification in place of the thermal one,
!       Ri_mu = (rho/P) g^2 (-d ln T/d ln mu) (d ln mu/d ln P)/(d omega/d ln r)^2,
!     -d ln T/d ln mu = (b - 1)/(1 + 3b) with b = a T^4/(3P).
!     mu_gradient_richardson_coeff holds everything but d ln mu.
!   - del_grad_diff_interface = max(grad_ad - grad + d ln mu/d ln P,
!     floor), the margin with the mu gradient (Maeder & Zahn 1998), by
!     which circulation inhibits the circulation (it adds d ln mu/d ln P
!     again, so the mu gradient counts twice, as in the F77 code; see
!     inhibit_by_mu_gradient).
! The floor is 1e-3, except at the transport region's first and last
! interfaces, where it is the larger of 1e-3 and half the neighbouring
! zone's grad_ad - grad (plus d ln mu/d ln P).
!
! The commented-out advection-diffusion transport (2026-10-12, the
! author's decision; see am_transport_coefficients) also needed the
! theta relaxation factor, the velocity_coeff* terms, the second- and
! third-derivative geometric factors and the order d omega/dr shear
! terms (difad_shear_coeff1/2, facd2/facd3), set here; they are
! commented out with it.
module stability_coefficients
      use rotation_configuration, only: rotation_settings, gradient_margin_floor
      implicit none
      private
      public :: rotation_stability_setup

contains

subroutine rotation_stability_setup(stab, sub_physics, &
     log_density, local_gravity, luminosity, log_pressure, &
     log_radius, mass_unlogged, log_temperature, transport_zone_begin, &
     transport_zone_end, num_zones)
      use rotation_scratch_lib, only: stability_state, &
           substep_physics_state
      use star_info_lib, only: json
      use physical_constants
      use math_lib
      implicit none
      type(stability_state), intent(inout) :: stab
      type(substep_physics_state), intent(in) :: sub_physics

      double precision, intent(in) :: log_density(json), local_gravity(json), &
           luminosity(json), log_pressure(json), log_radius(json), &
           mass_unlogged(json), log_temperature(json)
      integer, intent(in) :: transport_zone_begin, transport_zone_end, &
           num_zones
! The structure at the interfaces that only the coefficients below use
! (the rest is stab's): pressure, temperature and mass (cgs), the
! actual and adiabatic gradients, d ln rho/d ln T, eps (erg/g/s), c_P.
      double precision :: pressure_interface(json), temperature_interface(json), &
           mass_interface(json), gradT_interface(json), grada_interface(json), &
           qdt_interface(json), eps_total_interface(json)
      double precision :: specific_heat_interface(json)
      integer :: zone_idx, interior_begin, interior_end
      double precision :: grav_const
! (commented out with the advection-diffusion transport: log10 of
! r^4/(G m) at the interface, the theta term)
! (now in rotation/retired/advection_diffusion.f90.txt)
! d ln mu/d ln P; the floor of Delta and Delta itself
      double precision :: dlnmu_dlnp, delta_grad_floor, delta_grad
! H_P; (4/3) r^3/(G m); P/(g Delta c_P rho T);
! (8/3) r^3/(G m) P/(g Delta c_P rho T); L/m; 1/(pi G)
      double precision :: pressure_scale_height, &
           rotation_factor, thermal_factor, es_prefactor, specific_luminosity, &
           one_over_pi_g
! (commented out with the advection-diffusion transport: H_T = H_P/grad;
! F/rho = L/(4 pi rho r^2); H_T at the zones; the mean d eps/d ln T and
! neutrino_loss_fraction; (2/3) r^4/(G m); d ln chi/d ln T = 3 - d ln
! kappa/d ln T; and P/(g Delta_mu c_P rho T), with Delta_mu =
! del_grad_diff_interface)
! (now in rotation/retired/advection_diffusion.f90.txt)

! the unlogged radii of the zones (stab%radius), and the interfaces
! midway between
      do zone_idx = 1,num_zones
         stab%radius(zone_idx) = exp(ln10*log_radius(zone_idx))
      end do
      do zone_idx = 2,num_zones
         stab%radius_interface(zone_idx) = 0.5d0*(stab%radius(zone_idx) + &
              stab%radius(zone_idx-1))
      end do
! The interpolation weights and the interpolated structure. The
! interface between zones 1 and 2 uses zones 1-4, the one between the
! last two zones the last four, and every other interface i zones
! i-2..i+1.
      if(transport_zone_begin.lt.3)then
         call lagrange_weights(1, 2)
         interior_begin = 3
         call interpolate_to_interface(1, 2)
      else
         interior_begin = transport_zone_begin
      endif
      if(transport_zone_end.eq.num_zones)then
         call lagrange_weights(num_zones-3, num_zones)
         interior_end = num_zones-1
         call interpolate_to_interface(num_zones-3, num_zones)
      else
         interior_end = transport_zone_end
      endif
      do zone_idx = interior_begin,interior_end
         call lagrange_weights(zone_idx-2, zone_idx)
      end do
      grav_const = exp(ln10*rotation_settings%log10_gravitational_constant)
      do zone_idx = interior_begin,interior_end
         call interpolate_to_interface(zone_idx-2, zone_idx)
      end do
! the actual gradient: grad_rad was interpolated
      do zone_idx = transport_zone_begin,transport_zone_end
         gradT_interface(zone_idx) = min(gradT_interface(zone_idx),grada_interface(zone_idx))
      end do
      do zone_idx = 1,num_zones
         stab%es_velocity_coeff1(zone_idx) = 0.0d0
         stab%es_velocity_coeff2(zone_idx) = 0.0d0
         stab%es_shear_coeff(zone_idx) = 0.0d0
         stab%gsf_kippenhahn_coeff(zone_idx) = 0.0d0
         stab%mu_gradient_richardson_coeff(zone_idx) = 0.0d0
! (commented out with the advection-diffusion transport)
! (now in rotation/retired/advection_diffusion.f90.txt)
      end do
! (commented out with the advection-diffusion transport: r^4/(G m) at
! each interface, from the geometric means of r and m on either side,
! the d theta/dt term of Maeder & Zahn 1998)
! (now in rotation/retired/advection_diffusion.f90.txt)
! the coefficients at every interface of the transport region
      one_over_pi_g = 4.0d0/c4pi/grav_const
      associate (p_i => pressure_interface, t_i => temperature_interface, &
           rho_i => stab%density_interface, m_i => mass_interface, &
           r_i => stab%radius_interface, l_i => stab%luminosity_interface, &
           g_i => stab%gravity_factor_interface, eps_i => eps_total_interface, &
           grada_i => grada_interface, grad_i => gradT_interface)
      do zone_idx = transport_zone_begin,transport_zone_end
! (the commented-out advection-diffusion transport set it to zero)
         dlnmu_dlnp = (log10(sub_physics%mu(zone_idx))-log10(sub_physics%mu(zone_idx-1)))/ &
              (log_pressure(zone_idx)-log_pressure(zone_idx-1))
         if(zone_idx.eq.transport_zone_begin)then
            delta_grad_floor=max(gradient_margin_floor,0.5d0*(sub_physics%grada(transport_zone_begin+1)- &
                 sub_physics%gradT(transport_zone_begin+1)) &
                 +dlnmu_dlnp)
         else if(zone_idx.eq.transport_zone_end)then
            delta_grad_floor=max(gradient_margin_floor,0.5d0*(sub_physics%grada(transport_zone_end-1)- &
                 sub_physics%gradT(transport_zone_end-1)) &
                 +dlnmu_dlnp)
         else
            delta_grad_floor = gradient_margin_floor
         endif
         delta_grad = max(grada_i(zone_idx)-grad_i(zone_idx),delta_grad_floor)
         stab%del_grad_diff_interface(zone_idx) = &
              max(grada_i(zone_idx)-grad_i(zone_idx)+dlnmu_dlnp,delta_grad_floor)
! (commented out with the advection-diffusion transport: the theta
! relaxation factor P/(g Delta rho))
! (now in rotation/retired/advection_diffusion.f90.txt)
         pressure_scale_height = p_i(zone_idx)*r_i(zone_idx)**2/rho_i(zone_idx)/m_i(zone_idx)/grav_const
!         temperature_scale_height = pressure_scale_height/grad_i(zone_idx)
         rotation_factor = 2.0d0*cc23*r_i(zone_idx)**3/grav_const/m_i(zone_idx)
         stab%gsf_kippenhahn_coeff(zone_idx)= 8.0d0*pressure_scale_height*rotation_factor/delta_grad
         thermal_factor = p_i(zone_idx)/(g_i(zone_idx)*delta_grad*specific_heat_interface(zone_idx)* &
              rho_i(zone_idx)*t_i(zone_idx))
         es_prefactor = 2.0d0*rotation_factor*thermal_factor
         specific_luminosity = l_i(zone_idx)/m_i(zone_idx)
         stab%es_velocity_coeff1(zone_idx) = es_prefactor*(specific_luminosity-eps_i(zone_idx))
         stab%es_velocity_coeff2(zone_idx)= -0.5d0*es_prefactor*specific_luminosity*one_over_pi_g/rho_i(zone_idx)
         stab%es_shear_coeff(zone_idx) = specific_luminosity*thermal_factor*cc13*one_over_pi_g*r_i(zone_idx)/ &
              rho_i(zone_idx)
! (commented out with the advection-diffusion transport)
! (now in rotation/retired/advection_diffusion.f90.txt)
! the dynamical shear: the |d omega/d ln r| at which Ri = 1/4
         stab%dynamical_shear_omega_limit(zone_idx)=sqrt(rho_i(zone_idx)/p_i(zone_idx)* &
              max(grada_i(zone_idx) - grad_i(zone_idx),gradient_margin_floor)) &
              *2.0d0*g_i(zone_idx)
! the mu-gradient Richardson term, (rho/P) g^2 (-d ln T/d ln mu)/d ln P
         stab%mu_gradient_richardson_coeff(zone_idx) = rho_i(zone_idx)*g_i(zone_idx)**2/ &
              p_i(zone_idx)/ln10/(log_pressure(zone_idx)-log_pressure(zone_idx-1))* &
              (radiation_constant_over_3*t_i(zone_idx)**4/p_i(zone_idx)-1.0d0)/ &
              (1.0d0 + 3.0d0*radiation_constant_over_3*t_i(zone_idx)**4/p_i(zone_idx))
      end do
! (commented out: the advection-diffusion transport's terms of order
! d omega/dr, Maeder & Zahn 1998)
! (now in rotation/retired/advection_diffusion.f90.txt)
      end associate
      return

contains

!----------------------------------------------------------------------
! lagrange_weights
!----------------------------------------------------------------------
! The 4-point Lagrange interpolation weights at interface interface_idx
! (radius stab%radius_interface(interface_idx)) from zones
! first_point..first_point+3, into
! stab%lagrange_interp_weights(1:4,interface_idx) (numerics' interp).
subroutine lagrange_weights(first_point, interface_idx)
      use numerics_lib, only: interp
      integer, intent(in) :: first_point, interface_idx

      call interp(stab%radius(first_point:first_point+3), &
           stab%lagrange_interp_weights(:,interface_idx), &
           x_eval=stab%radius_interface(interface_idx))
end subroutine lagrange_weights

!----------------------------------------------------------------------
! interp4
!----------------------------------------------------------------------
! field at interface interface_idx, with lagrange_weights' weights from
! zones first_point..first_point+3 (numerics' lagrange4).
double precision function interp4(field, first_point, interface_idx)
      use numerics_lib, only: lagrange4
      double precision, intent(in) :: field(json)
      integer, intent(in) :: first_point, interface_idx
      interp4 = lagrange4(stab%lagrange_interp_weights(:,interface_idx), &
           field(first_point:first_point+3))
end function interp4

!----------------------------------------------------------------------
! interpolate_to_interface
!----------------------------------------------------------------------
! The structure at interface interface_idx from zones
! first_point..first_point+3: P, T and rho interpolated in log10 and
! unlogged; grad_rad (into gradT_interface, which the caller then
! limits by grad_ad), grad_ad, d ln rho/d ln T, the enclosed mass, eps,
! c_P and the local gravity linearly; the luminosity converted to cgs.
subroutine interpolate_to_interface(first_point, interface_idx)
      integer, intent(in) :: first_point, interface_idx

      pressure_interface(interface_idx)=exp(ln10*interp4(log_pressure, first_point, interface_idx))
      temperature_interface(interface_idx)=exp(ln10*interp4(log_temperature, first_point, interface_idx))
      stab%density_interface(interface_idx)=exp(ln10*interp4(log_density, first_point, interface_idx))
      gradT_interface(interface_idx)=interp4(sub_physics%gradr, first_point, interface_idx)
      grada_interface(interface_idx)=interp4(sub_physics%grada, first_point, interface_idx)
      qdt_interface(interface_idx)=interp4(sub_physics%qdt, first_point, interface_idx)
      mass_interface(interface_idx)=interp4(mass_unlogged, first_point, interface_idx)
      eps_total_interface(interface_idx)=interp4(sub_physics%eps_total, first_point, interface_idx)
      stab%luminosity_interface(interface_idx)=rotation_settings%solar_luminosity_cgs*(interp4(luminosity, first_point, &
           interface_idx))
      stab%gravity_factor_interface(interface_idx)=interp4(local_gravity, first_point, interface_idx)
      specific_heat_interface(interface_idx)=interp4(sub_physics%cp, first_point, interface_idx)
end subroutine interpolate_to_interface

end subroutine rotation_stability_setup

end module stability_coefficients
