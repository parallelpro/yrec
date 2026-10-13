!----------------------------------------------------------------------
! structure_equations
!----------------------------------------------------------------------
! The structure equations of the Henyey solve, across the model.
! compute_henyey_coefficients gathers each shell from the star, has the
! kernel evaluate_shell_equations (shell_equations.f90) form the four
! equations' terms and partials there, hands each pair of shells to
! henyey_eliminate (henyey_elimination.f90) as it goes, and stores the
! shells' values for output.
module structure_equations
      implicit none
      private
      public :: compute_henyey_coefficients

contains

!----------------------------------------------------------------------
! compute_henyey_coefficients (formerly coefft)
!----------------------------------------------------------------------
! Linearizes the structure equations at every shell, from the centre
! out, and eliminates as it goes. At each shell im it fills a
! shell_state from the star and calls evaluate_shell_equations, which
! returns the terms q and their partials (cur; prev holds the shell
! below) and the shell's values for output; it then calls
! henyey_eliminate for the pair (im-1, im), or sets the central boundary
! conditions at im = 1, and stores the output values. The results are
! star%elim_coeff and star%elim_rhs. It also stores per-zone physics for
! output: star%logRho, eps_total and eps_channels, gradr/gradT/grada,
! adiabatic_index_gamma1, the pulse_* arrays, the luminosity channels and
! the neutrino fluxes. The equations themselves, with their units and
! partials, are in shell_equations.f90's header.
! dlnrho_dlnt and dlnrho_dlnp: on entry the starting point of the first
! shell's eos; on return the last shell's values.
subroutine compute_henyey_coefficients(delta_time, dlnrho_dlnt, dlnrho_dlnp, &
     saha_state, envelope_zone_index, ierr)
      use henyey_elimination
      use shell_equations, only: henyey_shell_terms, shell_state, &
           structure_settings, shell_diagnostics, evaluate_shell_equations
      use star_info_lib, only: star, i_eps_grav, i_eps_neu, i_h1, i_lum_grav, &
           i_lum_he_c, n_lum_channels, n_nu_fluxes
      use eos_lib, only: eos_result
      use kap_lib, only: i_kap, i_dlnkap_dlnrho, i_dlnkap_dlnt
      use physical_constants
      use math_lib
      implicit none

      double precision, intent(in) :: delta_time
      double precision, intent(inout) :: dlnrho_dlnt, dlnrho_dlnp
      integer, intent(inout) :: saha_state
      integer, intent(in) :: envelope_zone_index
      integer, intent(out) :: ierr

! --- locals ---
! The eos result, carried from shell to shell: each shell's eos starts
! from the previous shell's values.
      type(eos_result) :: eos_res
      type(structure_settings) :: cfg
      type(shell_state) :: sh
      type(shell_diagnostics) :: diag
! The equation terms at the current shell (cur) and, carried from the
! previous pass of the loop, at the shell below (prev).
      type(henyey_shell_terms) :: prev, cur
      double precision :: energy_sum_inverse
      integer :: im, im1, j

      ierr = 0

      call start_energy_bookkeeping()
      call set_structure_settings()
      eos_res%dlnrho_dlnt = dlnrho_dlnt
      eos_res%dlnrho_dlnp = dlnrho_dlnp
      do im = 1,star%nz
       call gather_shell_state()
       call evaluate_shell_equations(sh, cfg, eos_res, saha_state, cur, diag, ierr)
       if (ierr /= 0) return
       dlnrho_dlnt = eos_res%dlnrho_dlnt
       dlnrho_dlnp = eos_res%dlnrho_dlnp
       call store_shell_structure()
       call add_shell_energy_budget()
       call eliminate_shell_pair()
       call store_zone_energy_channels()
       call store_zone_thermodynamics()
       call store_pulse_derivatives()
      end do

      return

contains

! ---------------------------------------------------------------
! Zero the solar-neutrino fluxes (star%neutrino_flux_total(j), the total
! flux of neutrino species j made in the star) and the luminosity
! channels.
subroutine start_energy_bookkeeping()
      if (star%ctrl%calc_neutrinos) then
         do j = 1,n_nu_fluxes
            star%neutrino_flux_total(j) = 0.0d0
         end do
      end if
      do j = 1,n_lum_channels
       star%luminosity_breakdown(j) = 0.0d0
      end do
end subroutine start_energy_bookkeeping

! ---------------------------------------------------------------
! The settings every shell shares: the time step, the burning cutoff,
! log Teff (spots), the turnover time (the deuterium cap), ln10/Lsun and
! the rotation switch.
subroutine set_structure_settings()
      cfg%dt = delta_time
      cfg%burning_log_temperature_cutoff = star%ctrl%nuclear_logT_cutoffs(1)
      cfg%log_teff = star%log_Teff
      cfg%turnover_timescale = star%convective_turnover_timescale
      cfg%ln10_over_lsun = star%ln10_over_lsun
      cfg%rotation_active = star%job%rotation_active
end subroutine set_structure_settings

! ---------------------------------------------------------------
! Shell im's state. With accretion, the zones at or above
! envelope_zone_index add the accretion shifts of log T and log P
! (mass_change.f90) to their change over the step.
subroutine gather_shell_state()
       sh%log_mass = star%log_mass(im)
       sh%mass = star%m(im)
       sh%shell_mass = star%dm(im)
       sh%log_temperature = star%logT(im)
       sh%log_pressure = star%logP(im)
       sh%log_radius = star%logR(im)
       sh%luminosity_lsun = star%luminosity_lsun(im)
       sh%log_density = star%logRho(im)
       sh%xa = star%xa(:,im)
       sh%fp = star%fp_rot(im)
       sh%ft = star%ft_rot(im)
       sh%log_temperature_delta = star%log_temperature_delta(im)
       sh%log_pressure_delta = star%log_pressure_delta(im)
       if (star%job%use_mass_accretion.and.star%ctrl%mass_accretion_rate.gt.0.0d0) then
          if (im.ge.envelope_zone_index) then
             sh%log_temperature_delta = star%log_temperature_delta(im)+ &
                  star%delta_log_temperature
             sh%log_pressure_delta = star%log_pressure_delta(im)+ &
                  star%delta_log_pressure
          end if
       end if
       sh%rotational_ke = star%kinetic_energy_rot(im)
       sh%rotational_ke_start = star%kinetic_energy_rot_old(im)
       sh%in_overshoot = im >= star%iov1 .and. im <= star%iov2
       sh%in_surface_cz = im.ge.star%jcz
end subroutine gather_shell_state

! ---------------------------------------------------------------
! The shell's density and convective flag.
subroutine store_shell_structure()
       star%logRho(im) = diag%log_density
       star%convective_flag(im) = diag%convective
end subroutine store_shell_structure

! ---------------------------------------------------------------
! The shell's share of the luminosity channels (dm eps_j/Lsun: pp1, pp2,
! pp3, CNO, triple alpha, the neutrino losses, alpha captures in
! i_lum_he_c, eps_grav in i_lum_grav), its He3+He3 and He3-burning
! luminosities (Lsun), and its solar-neutrino fluxes: burn%neutrino_flux
! is the flux at the Earth per gram of the shell (10^10 cm^-2 s^-1 g^-1,
! nuclear_energy.f90), so times dm it is the shell's flux.
subroutine add_shell_energy_budget()
       if (diag%burned) then
            if (star%ctrl%calc_neutrinos) then
               do j = 1,n_nu_fluxes
                  star%neutrino_flux_total(j) = star%neutrino_flux_total(j) + &
                       diag%burn%neutrino_flux(j)*star%dm(im)
               end do
            end if
            do j = 1,6
               star%luminosity_breakdown(j) = star%luminosity_breakdown(j) + &
                    (star%dm(im)/star%solar_luminosity_cgs)* &
                    diag%energy_gen_component(j)
            end do
            star%he3_he3_luminosity_zone(im) = (star%dm(im)/ &
                 star%solar_luminosity_cgs)*diag%burn%eps_he3_he3
            star%he3_burning_luminosity_zone(im) = (star%dm(im)/ &
                 star%solar_luminosity_cgs)*diag%burn%eps_he3_burning
            star%luminosity_breakdown(i_lum_he_c)=star%luminosity_breakdown(i_lum_he_c)+(star%dm(im)/ &
                 star%solar_luminosity_cgs)*diag%burn%eps_alpha_capture
       end if
       if (diag%has_eps_grav) then
            star%gravitational_luminosity(im) = diag%eps_grav
            star%luminosity_breakdown(i_lum_grav) = star%luminosity_breakdown(i_lum_grav) + &
                 (star%dm(im)/star%solar_luminosity_cgs)*diag%eps_grav
       end if
end subroutine add_shell_energy_budget

! ---------------------------------------------------------------
! Eliminate the pair (im-1, im), or set the central boundary conditions
! at im = 1; carry cur into prev.
subroutine eliminate_shell_pair()
         if (im.gt.1) then
            im1 = im
            call henyey_eliminate(im1, prev, cur)
         else
! The central boundary conditions, rows 3 and 4 at shell 1 in the form
! of henyey_elimination.f90,
!   d y = elim_rhs - elim_coeff(:,1) d log P - elim_coeff(:,2) d log T.
! The innermost point is taken as a uniform sphere, m = (4 pi/3) r^3 rho:
!   log R = (log m - log(4 pi/3) - log rho)/3,
! so elim_rhs(3,1) is that minus log R and elim_coeff(3,:,1) is
! (dln rho/dln P, dln rho/dln T)/3. Its luminosity is L = eps m/Lsun =
! ql/ln10, so elim_rhs(4,1) = ql/ln10 - L and elim_coeff(4,:,1) =
! -(ql_dp, ql_dt).
            star%elim_coeff(3,1,1) = cc13*dlnrho_dlnp
            star%elim_coeff(3,2,1) = cc13*dlnrho_dlnt
            star%elim_rhs(3,1) = -cc13*(c4pi3l + diag%log_density - &
                 sh%log_mass) - sh%log_radius
            star%elim_coeff(4,1,1) = -cur%ql_dp
            star%elim_coeff(4,2,1) = -cur%ql_dt
            star%elim_rhs(4,1) = clni*cur%ql - sh%luminosity_lsun
         end if
! Carry this shell's terms into the next pair's elimination.
         prev = cur
end subroutine eliminate_shell_pair

! ---------------------------------------------------------------
! Store the shell's energy generation for output, for every model:
!   eps_total = eps_pp1 + eps_pp2 + eps_pp3 + eps_cno + eps_3alpha (erg/g/s),
!   eps_channels(1:5) = the fraction of eps_total in each of those,
!   eps_channels(i_eps_neu) = eps_nu and eps_channels(i_eps_grav) =
!   eps_grav (erg/g/s).
! Below the burning cutoff only eps_grav is nonzero.
subroutine store_zone_energy_channels()
         if (star%logT(im).lt.star%ctrl%nuclear_logT_cutoffs(1)) then
            star%eps_total(im) = 0.0d0
            star%eps_channels(i_eps_grav,im) = star%gravitational_luminosity(im)
            do j = 1,6
               star%eps_channels(j,im) = 0.0d0
           end do
         else
            star%eps_total(im) = diag%energy_gen_component(1)+diag%energy_gen_component(2)+ &
                 diag%energy_gen_component(3)+diag%energy_gen_component(4)+ &
                 diag%energy_gen_component(5)
            star%eps_channels(i_eps_neu,im) = diag%energy_gen_component(6)
            star%eps_channels(i_eps_grav,im) = star%gravitational_luminosity(im)
            if (star%eps_total(im).gt.1.0d-22) then
               energy_sum_inverse = 1.0d0/star%eps_total(im)
            else
               energy_sum_inverse = 0.0d0
            end if
            do j = 1,5
               star%eps_channels(j,im) = diag%energy_gen_component(j)*energy_sum_inverse
              end do
         end if
end subroutine store_zone_energy_channels

! ---------------------------------------------------------------
! Store the shell's eos, opacity and gradient values and Gamma1 for output.
subroutine store_zone_thermodynamics()
         star%beta(im) = eos_res%beta
         star%eta(im) = eos_res%electron_degeneracy_parameter
         star%opacity_zone(im) = diag%kap_res(i_kap)
         star%gradr(im) = diag%radiative_gradient
         star%gradT(im) = diag%actual_gradient
         star%grada(im) = eos_res%adiabatic_gradient
         do j = 1,3
            star%fxion_zone(j,im) = eos_res%ion_fraction(j)
         end do
         star%conv_vel(im) = diag%convective_velocity
         star%scp(im) = eos_res%specific_heat_cp
         star%adiabatic_index_gamma1(im) = diag%gamma1

! Disabled with rotation's advection-diffusion transport (2026-10-12),
! their only reader: the copies of dln kappa/dln T, d eps/dln T and
! (eps - eps_nu - eps_grav)/eps for rotation_stability_setup. They come
! back with that transport.
!         if (star%job%rotation_active) then
!            rot_scr%dlnkappa_dlnt(im) = diag%kap_res(i_dlnkap_dlnt)
!            if (star%eps_total(im).gt.0.0d0) then
!               total_energy_sum = star%eps_total(im)
!               neutrino_and_grav_sum = star%eps_channels(i_eps_neu,im)+star%eps_channels(i_eps_grav,im)
!               rot_scr%neutrino_loss_fraction(im) = (total_energy_sum - &
!                    neutrino_and_grav_sum)/total_energy_sum
!            else
!               rot_scr%neutrino_loss_fraction(im) = 0.0d0
!            end if
!            rot_scr%deps_dlnt(im) = diag%deps_dlnt
!         end if
end subroutine store_zone_thermodynamics

! ---------------------------------------------------------------
! Store the shell's derivatives and mean molecular weights for the
! pulsation output, for every model: the GYRE-format writer
! (pulsation_output's write_gyre_pulse) runs every pulse_gyre_interval
! models, independently of the older OPAL pulsation writer's flag, and
! needs them (MHP 8/25). No other output file reads these arrays. Below
! the burning cutoff pulse_deps_dlnrho/dlnt are zero.
subroutine store_pulse_derivatives()
         star%pulse_dlnrho_dlnp(im) = dlnrho_dlnp
         star%pulse_deps_dlnrho(im) = diag%deps_dlnrho
         star%pulse_deps_dlnt(im) = diag%deps_dlnt
         star%pulse_dlnkap_dlnrho(im) = diag%kap_res(i_dlnkap_dlnrho)
         star%pulse_dlnkap_dlnt(im) = diag%kap_res(i_dlnkap_dlnt)
         star%pulse_specific_heat(im) = eos_res%specific_heat_cp
! eos_res%specific_gas_constant is R/mu, so mu = R/specific_gas_constant;
! the electron slot is 1/mu_e straight from the eos.
         if (eos_res%specific_gas_constant > 0.0d0) then
            star%pulse_mean_molecular_weight(im) = gas_constant/eos_res%specific_gas_constant
         else
            star%pulse_mean_molecular_weight(im) = 0.0d0
         end if
         star%pulse_electron_mean_weight_inverse(im) = eos_res%electron_mean_weight_inverse
         star%pulse_dlnrho_dlnt(im) = dlnrho_dlnt
! These MLT pulse-output arrays are always zero: the alfmlt/phmlt/cmxmlt
! scalars they once copied were only ever zeroed, and are gone. Whether
! to keep the arrays awaits an author decision (audit SUMMARY.md 1.1 #10).
         star%valfmlt(im) = 0.0d0
         star%vphmlt(im) = 0.0d0
         star%vcmxmlt(im) = 0.0d0
end subroutine store_pulse_derivatives

end subroutine compute_henyey_coefficients

end module structure_equations
