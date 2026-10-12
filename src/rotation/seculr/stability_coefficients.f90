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
!     es_shear_coeff * omega d omega/dr.
!   - The GSF instability (Kippenhahn 1980):
!     gsf_kippenhahn_coeff = (32/3) H_P r^3/(G m Delta). circulation
!     multiplies it by the thermal diffusivity, omega^2/r and the
!     gradient of omega or of the specific angular momentum.
!   - The dynamical shear instability: an interface is stable while
!       Ri = (rho/P) (grad_ad - grad) g^2/(d omega/d ln r)^2 > 1/4,
!     so dynamical_shear_omega_limit = 2 g sqrt((rho/P) max(grad_ad -
!     grad, 1e-3)) is the largest stable |d omega/d ln r|.
!   - The mu-gradient Richardson term: the same criterion with the mu
!     stratification in place of the thermal one,
!       Ri_mu = (rho/P) g^2 (-d ln T/d ln mu) (d ln mu/d ln P)/(d omega/d ln r)^2,
!     -d ln T/d ln mu = (b - 1)/(1 + 3b) with b = a T^4/(3P).
!     mu_gradient_richardson_coeff holds everything but d ln mu.
!   - del_grad_diff_interface = max(grad_ad - grad + d ln mu/d ln P,
!     floor), the margin with the mu gradient (Maeder & Zahn 1998), by
!     which circulation inhibits the circulation.
! The floor is 1e-3, except at the transport region's first and last
! interfaces, where it is the larger of 1e-3 and half the neighbouring
! zone's grad_ad - grad (plus d ln mu/d ln P).
!
! The advection-diffusion transport (use_diffusion_advection_transport;
! Maeder & Zahn 1998) also needs the theta relaxation factor, the
! velocity_coeff* terms, the second- and third-derivative geometric
! factors and the order d omega/dr shear terms (difad_shear_coeff1/2,
! facd2/facd3), set here.
module stability_coefficients
      use rotation_configuration, only: rotation_settings
      implicit none
      private
      public :: rotation_stability_setup

contains

subroutine rotation_stability_setup(log_density, local_gravity, luminosity, log_pressure, &
     log_radius, mass_unlogged, log_temperature, transport_zone_begin, &
     transport_zone_end, num_zones, radius_unlogged, &
     dynamical_shear_omega_limit)
      use rotation_scratch_lib
      use star_info_lib, only: star, json
      use physical_constants
      use math_lib
      implicit none

      double precision, intent(in) :: log_density(json), local_gravity(json), &
           luminosity(json), log_pressure(json), log_radius(json), &
           mass_unlogged(json), log_temperature(json)
      integer, intent(in) :: transport_zone_begin, transport_zone_end, &
           num_zones
      double precision, intent(out) :: radius_unlogged(json)
      double precision, intent(out) :: dynamical_shear_omega_limit(json)
      double precision :: specific_heat_interface(json)
      integer :: zone_idx, interior_begin, interior_end
      double precision :: grav_const
! log10 of r^4/(G m) at the interface (the theta term)
      double precision :: log_theta_mean
! d ln mu/d ln P; the floor of Delta and Delta itself
      double precision :: dlnmu_dlnp, delta_grad_floor, delta_grad
! H_P and H_T = H_P/grad; (4/3) r^3/(G m); P/(g Delta c_P rho T);
! (8/3) r^3/(G m) P/(g Delta c_P rho T); L/m; 1/(pi G)
      double precision :: pressure_scale_height, temperature_scale_height, &
           rotation_factor, thermal_factor, es_prefactor, specific_luminosity, &
           one_over_pi_g
! the advection-diffusion terms: F/rho = L/(4 pi rho r^2); H_T at the
! zones; the mean d eps/d ln T and neutrino_loss_fraction;
! (2/3) r^4/(G m); d ln chi/d ln T = 3 - d ln kappa/d ln T; and
! P/(g Delta_mu c_P rho T), with Delta_mu = del_grad_diff_interface
      double precision :: flux_over_density
      double precision :: temperature_scale_height_below, temperature_scale_height_zone
      double precision :: mean_deps_dlnt, mean_neutrino_fraction
      double precision :: second_deriv_factor, dlnchi_dlnt, thermal_factor_mu

! the unlogged radii of the zones, and the interfaces midway between
      do zone_idx = 1,num_zones
         radius_unlogged(zone_idx) = exp(ln10*log_radius(zone_idx))
      end do
      do zone_idx = 2,num_zones
         rot_scr%radius_interface(zone_idx) = 0.5d0*(radius_unlogged(zone_idx) + &
              radius_unlogged(zone_idx-1))
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
      grav_const = exp(ln10*cgl)
      do zone_idx = interior_begin,interior_end
         call interpolate_to_interface(zone_idx-2, zone_idx)
      end do
! the actual gradient: grad_rad was interpolated
      do zone_idx = transport_zone_begin,transport_zone_end
         rot_scr%gradT_interface(zone_idx) = min(rot_scr%gradT_interface(zone_idx),rot_scr%grada_interface(zone_idx))
      end do
      do zone_idx = 1,num_zones
         rot_scr%es_velocity_coeff1(zone_idx) = 0.0d0
         rot_scr%es_velocity_coeff2(zone_idx) = 0.0d0
         rot_scr%es_shear_coeff(zone_idx) = 0.0d0
         rot_scr%gsf_kippenhahn_coeff(zone_idx) = 0.0d0
         rot_scr%mu_gradient_richardson_coeff(zone_idx) = 0.0d0
         rot_scr%es_relaxation_factor(zone_idx) = 0.0d0
         rot_scr%theta_mean(zone_idx) = 0.0d0
         rot_scr%difad_shear_coeff1(zone_idx) = 0.0d0
         rot_scr%difad_shear_coeff2(zone_idx) = 0.0d0
      end do
! r^4/(G m) at each interface, from the geometric means of r and m on
! either side: the d theta/dt term of Maeder & Zahn (1998)
      do zone_idx = 2,num_zones
         log_theta_mean = 2.0d0*(log_radius(zone_idx)+log_radius(zone_idx-1))-0.5d0* &
              (log10(mass_unlogged(zone_idx))+log10(mass_unlogged(zone_idx-1)))-cgl
         rot_scr%theta_mean(zone_idx) = exp(ln10*log_theta_mean)
      end do
! the coefficients at every interface of the transport region
      one_over_pi_g = 4.0d0/c4pi/grav_const
      associate (p_i => rot_scr%pressure_interface, t_i => rot_scr%temperature_interface, &
           rho_i => rot_scr%density_interface, m_i => rot_scr%mass_interface, &
           r_i => rot_scr%radius_interface, l_i => rot_scr%luminosity_interface, &
           g_i => rot_scr%gravity_factor_interface, eps_i => rot_scr%eps_total_interface, &
           grada_i => rot_scr%grada_interface, grad_i => rot_scr%gradT_interface)
      do zone_idx = transport_zone_begin,transport_zone_end
         if(.not.rotation_settings%use_diffusion_advection_transport)then
            dlnmu_dlnp = (log10(mid_scr%mu(zone_idx))-log10(mid_scr%mu(zone_idx-1)))/ &
                 (log_pressure(zone_idx)-log_pressure(zone_idx-1))
         else
            dlnmu_dlnp = 0.0d0
         endif
         if(zone_idx.eq.transport_zone_begin)then
            delta_grad_floor=max(1.0d-3,0.5d0*(mid_scr%grada(transport_zone_begin+1)- &
                 mid_scr%gradT(transport_zone_begin+1)) &
                 +dlnmu_dlnp)
         else if(zone_idx.eq.transport_zone_end)then
            delta_grad_floor=max(1.0d-3,0.5d0*(mid_scr%grada(transport_zone_end-1)- &
                 mid_scr%gradT(transport_zone_end-1)) &
                 +dlnmu_dlnp)
         else
            delta_grad_floor = 1.0d-3
         endif
         delta_grad = max(grada_i(zone_idx)-grad_i(zone_idx),delta_grad_floor)
         rot_scr%del_grad_diff_interface(zone_idx) = &
              max(grada_i(zone_idx)-grad_i(zone_idx)+dlnmu_dlnp,delta_grad_floor)
! the theta relaxation factor P/(g Delta rho)
         rot_scr%es_relaxation_factor(zone_idx) = p_i(zone_idx)/(g_i(zone_idx)*delta_grad*rho_i(zone_idx))
         pressure_scale_height = p_i(zone_idx)*r_i(zone_idx)**2/rho_i(zone_idx)/m_i(zone_idx)/grav_const
         temperature_scale_height = pressure_scale_height/grad_i(zone_idx)
         rotation_factor = 2.0d0*cc23*r_i(zone_idx)**3/grav_const/m_i(zone_idx)
         rot_scr%gsf_kippenhahn_coeff(zone_idx)= 8.0d0*pressure_scale_height*rotation_factor/delta_grad
         thermal_factor = p_i(zone_idx)/(g_i(zone_idx)*delta_grad*specific_heat_interface(zone_idx)* &
              rho_i(zone_idx)*t_i(zone_idx))
         es_prefactor = 2.0d0*rotation_factor*thermal_factor
         specific_luminosity = l_i(zone_idx)/m_i(zone_idx)
         rot_scr%es_velocity_coeff1(zone_idx) = es_prefactor*(specific_luminosity-eps_i(zone_idx))
         rot_scr%es_velocity_coeff2(zone_idx)= -0.5d0*es_prefactor*specific_luminosity*one_over_pi_g/rho_i(zone_idx)
         rot_scr%es_shear_coeff(zone_idx) = specific_luminosity*thermal_factor*cc13*one_over_pi_g*r_i(zone_idx)/ &
              rho_i(zone_idx)
         rot_scr%velocity_coeff0(zone_idx) = thermal_factor
         rot_scr%velocity_coeff1a(zone_idx) = rot_scr%es_velocity_coeff1(zone_idx)/thermal_factor
         rot_scr%velocity_coeff1b(zone_idx) = rot_scr%es_velocity_coeff2(zone_idx)/thermal_factor
         rot_scr%velocity_coeff2a(zone_idx) = specific_luminosity*one_over_pi_g*3.0d0*temperature_scale_height/ &
              rho_i(zone_idx)
         rot_scr%velocity_coeff2b(zone_idx) = -specific_luminosity*one_over_pi_g*r_i(zone_idx)/rho_i(zone_idx)
! the dynamical shear: the |d omega/d ln r| at which Ri = 1/4
         dynamical_shear_omega_limit(zone_idx)=sqrt(rho_i(zone_idx)/p_i(zone_idx)* &
              max(grada_i(zone_idx) - grad_i(zone_idx),1.0d-3)) &
              *2.0d0*g_i(zone_idx)
! the mu-gradient Richardson term, (rho/P) g^2 (-d ln T/d ln mu)/d ln P
         rot_scr%mu_gradient_richardson_coeff(zone_idx) = rho_i(zone_idx)*g_i(zone_idx)**2/ &
              p_i(zone_idx)/ln10/(log_pressure(zone_idx)-log_pressure(zone_idx-1))* &
              (radiation_constant_over_3*t_i(zone_idx)**4/p_i(zone_idx)-1.0d0)/ &
              (1.0d0 + 3.0d0*radiation_constant_over_3*t_i(zone_idx)**4/p_i(zone_idx))
      end do
! The advection-diffusion transport's terms of order d omega/dr
! (Maeder & Zahn 1998).
      if(rotation_settings%use_diffusion_advection_transport)then
         temperature_scale_height_below = exp(ln10*(log_pressure(transport_zone_begin-1)+ &
              2.0d0*log_radius(transport_zone_begin-1)-log_density(transport_zone_begin-1))) &
              /mass_unlogged(transport_zone_begin-1)/grav_const/mid_scr%gradr(transport_zone_begin-1)
         rot_scr%third_deriv_geom_factor(transport_zone_begin-1) = temperature_scale_height_below
         do zone_idx = transport_zone_begin,transport_zone_end
         pressure_scale_height = p_i(zone_idx)*r_i(zone_idx)**2/rho_i(zone_idx)/m_i(zone_idx)/grav_const
         temperature_scale_height = pressure_scale_height/grad_i(zone_idx)
         thermal_factor_mu = p_i(zone_idx)/(g_i(zone_idx)*rot_scr%del_grad_diff_interface(zone_idx)* &
              specific_heat_interface(zone_idx)*rho_i(zone_idx)*t_i(zone_idx))
         second_deriv_factor = cc23*r_i(zone_idx)**4/grav_const/m_i(zone_idx)
         rot_scr%second_deriv_geom_factor(zone_idx) = second_deriv_factor
         dlnchi_dlnt = 3.0d0 - 0.5d0*(rot_scr%dlnkappa_dlnt(zone_idx)+rot_scr%dlnkappa_dlnt(zone_idx-1))
         temperature_scale_height_zone = exp(ln10*(log_pressure(zone_idx)+2.0d0*log_radius(zone_idx)- &
              log_density(zone_idx)))/mass_unlogged(zone_idx)/grav_const/mid_scr%gradT(zone_idx)
         rot_scr%third_deriv_geom_factor(zone_idx) = temperature_scale_height_zone
         mean_deps_dlnt = 0.5d0*(rot_scr%deps_dlnt(zone_idx)+rot_scr%deps_dlnt(zone_idx-1))
         mean_neutrino_fraction = 0.5d0*(rot_scr%neutrino_loss_fraction(zone_idx)+ &
              rot_scr%neutrino_loss_fraction(zone_idx-1))
         rot_scr%difad_shear_coeff1(zone_idx)= &
         - thermal_factor_mu*second_deriv_factor*(mean_neutrino_fraction*mean_deps_dlnt + &
              eps_i(zone_idx)*(1.0d0 - mean_neutrino_fraction - dlnchi_dlnt))
         rot_scr%velocity_coeff2b(zone_idx) = rot_scr%velocity_coeff2b(zone_idx)+ &
              rot_scr%difad_shear_coeff1(zone_idx)/thermal_factor_mu
         rot_scr%difad_shear_coeff2(zone_idx) = 0.0d0
         flux_over_density = l_i(zone_idx)/c4pi/rho_i(zone_idx)/r_i(zone_idx)**2
         rot_scr%facd2(zone_idx) = thermal_factor_mu*(flux_over_density*dlnchi_dlnt - &
              eps_i(zone_idx)*temperature_scale_height)
         rot_scr%facd3(zone_idx) = -thermal_factor_mu*flux_over_density
         end do
      endif
      end associate
      return

contains

!----------------------------------------------------------------------
! lagrange_weights
!----------------------------------------------------------------------
! The 4-point Lagrange interpolation weights at interface interface_idx
! (radius rot_scr%radius_interface(interface_idx)) from zones
! first_point..first_point+3, into
! rot_scr%lagrange_interp_weights(1:4,interface_idx).
subroutine lagrange_weights(first_point, interface_idx)
      integer, intent(in) :: first_point, interface_idx
      double precision :: dr43, dr42, dr41, dr32, dr31, dr21
      double precision :: lag_denom1, lag_denom2, lag_denom3, lag_denom4
      double precision :: lag_x1, lag_x2, lag_x3, lag_x4

      dr43 = radius_unlogged(first_point+3) - radius_unlogged(first_point+2)
      dr42 = radius_unlogged(first_point+3) - radius_unlogged(first_point+1)
      dr41 = radius_unlogged(first_point+3) - radius_unlogged(first_point)
      dr32 = radius_unlogged(first_point+2) - radius_unlogged(first_point+1)
      dr31 = radius_unlogged(first_point+2) - radius_unlogged(first_point)
      dr21 = radius_unlogged(first_point+1) - radius_unlogged(first_point)
      lag_denom1 = -dr21*dr31*dr41
      lag_denom2 = dr21*dr32*dr42
      lag_denom3 = -dr31*dr32*dr43
      lag_denom4 = dr41*dr42*dr43
      lag_x1 = rot_scr%radius_interface(interface_idx) - radius_unlogged(first_point)
      lag_x2 = rot_scr%radius_interface(interface_idx) - radius_unlogged(first_point+1)
      lag_x3 = rot_scr%radius_interface(interface_idx) - radius_unlogged(first_point+2)
      lag_x4 = rot_scr%radius_interface(interface_idx) - radius_unlogged(first_point+3)
      rot_scr%lagrange_interp_weights(1,interface_idx) = (lag_x2*lag_x3*lag_x4)/lag_denom1
      rot_scr%lagrange_interp_weights(2,interface_idx) = (lag_x1*lag_x3*lag_x4)/lag_denom2
      rot_scr%lagrange_interp_weights(3,interface_idx) = (lag_x1*lag_x2*lag_x4)/lag_denom3
      rot_scr%lagrange_interp_weights(4,interface_idx) = (lag_x1*lag_x2*lag_x3)/lag_denom4
end subroutine lagrange_weights

!----------------------------------------------------------------------
! interp4
!----------------------------------------------------------------------
! field at interface interface_idx, with lagrange_weights' weights from
! zones first_point..first_point+3.
double precision function interp4(field, first_point, interface_idx)
      double precision, intent(in) :: field(json)
      integer, intent(in) :: first_point, interface_idx
      interp4 = field(first_point)*rot_scr%lagrange_interp_weights(1,interface_idx)+ &
           field(first_point+1)*rot_scr%lagrange_interp_weights(2,interface_idx)+ &
           field(first_point+2)*rot_scr%lagrange_interp_weights(3,interface_idx)+ &
           field(first_point+3)*rot_scr%lagrange_interp_weights(4,interface_idx)
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

      rot_scr%pressure_interface(interface_idx)=exp(ln10*interp4(log_pressure, first_point, interface_idx))
      rot_scr%temperature_interface(interface_idx)=exp(ln10*interp4(log_temperature, first_point, interface_idx))
      rot_scr%density_interface(interface_idx)=exp(ln10*interp4(log_density, first_point, interface_idx))
      rot_scr%gradT_interface(interface_idx)=interp4(mid_scr%gradr, first_point, interface_idx)
      rot_scr%grada_interface(interface_idx)=interp4(mid_scr%grada, first_point, interface_idx)
      rot_scr%qdt_interface(interface_idx)=interp4(mid_scr%qdt, first_point, interface_idx)
      rot_scr%mass_interface(interface_idx)=interp4(mass_unlogged, first_point, interface_idx)
      rot_scr%eps_total_interface(interface_idx)=interp4(mid_scr%eps_total, first_point, interface_idx)
      rot_scr%luminosity_interface(interface_idx)=star%solar_luminosity_cgs*(interp4(luminosity, first_point, &
           interface_idx))
      rot_scr%gravity_factor_interface(interface_idx)=interp4(local_gravity, first_point, interface_idx)
      specific_heat_interface(interface_idx)=interp4(mid_scr%cp, first_point, interface_idx)
end subroutine interpolate_to_interface

end subroutine rotation_stability_setup

end module stability_coefficients
