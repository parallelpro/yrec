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
      use star_info_lib, only: star_info
      implicit none
      private
      public :: compute_henyey_coefficients, luminosity_channel_sum

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
! Each shell's eos starts from the previous shell's result, the first
! from a fresh eos_result. (Until 2026-10 the last pass's dln rho/dln P
! and dln rho/dln T were carried into the next pass through the star, as
! a starting value every eos family overwrites before reading it;
! docs/bug_sweep.md O23.)
subroutine compute_henyey_coefficients(s, delta_time, saha_iterations, &
     envelope_zone_index, ierr)
      use henyey_elimination
      use shell_equations, only: henyey_shell_terms, shell_state, &
           structure_settings, shell_diagnostics, evaluate_shell_equations
      use star_info_lib, only: i_eps_grav, i_eps_neu, i_h1, i_lum_grav, &
           i_lum_neu, i_lum_extra, n_lum_channels, n_nu_fluxes
      use net_lib, only: n_chains
      use eos_lib, only: eos_result
      use kap_lib, only: i_kap, i_dlnkap_dlnrho, i_dlnkap_dlnt
      use physical_constants
      use math_lib
      implicit none
      type(star_info), intent(inout) :: s

      double precision, intent(in) :: delta_time
      integer, intent(inout) :: saha_iterations
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
! the current shell's dln rho/dln P and dln rho/dln T (the central
! boundary condition and the pulse arrays)
      double precision :: dlnrho_dlnt, dlnrho_dlnp
      integer :: im, im1, j, chain

      ierr = 0

      call start_energy_bookkeeping()
      call set_structure_settings()
      do im = 1,s%nz
       call gather_shell_state()
       call evaluate_shell_equations(sh, cfg, eos_res, saha_iterations, cur, diag, ierr)
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

! --------------------------------------------------------------- Zero
! the solar-neutrino fluxes (star%neutrino_flux_total_e10_per_cm2_s(j),
! the total flux of neutrino species j made in the star) and the
! luminosity channels.
subroutine start_energy_bookkeeping()
      if (s%ctrl%calc_neutrinos) then
         do j = 1,n_nu_fluxes
            s%neutrino_flux_total_e10_per_cm2_s(j) = 0.0d0
         end do
      end if
      do j = 1,n_lum_channels
       s%luminosity_breakdown_lsun(j) = 0.0d0
      end do
end subroutine start_energy_bookkeeping

! ---------------------------------------------------------------
! The settings every shell shares: the time step, log Teff (spots), the
! turnover time (the deuterium cap), ln10/Lsun, the rotation switch and
! log10 G.
subroutine set_structure_settings()
      cfg%dt = delta_time
      cfg%log_teff = s%log_Teff
      cfg%turnover_timescale = s%convective_turnover_timescale
      cfg%ln10_over_lsun = s%ln10_over_lsun
      cfg%rotation_active = s%ctrl%rotation_active
      cfg%log10_gravitational_constant = cgl
end subroutine set_structure_settings

! ---------------------------------------------------------------
! Shell im's state. With accretion, the zones at or above
! envelope_zone_index add the accretion shifts of log T and log P
! (mass_change.f90) to their change over the step.
subroutine gather_shell_state()
       sh%log_mass = s%log_mass(im)
       sh%mass = s%m(im)
       sh%shell_mass = s%dm(im)
       sh%log_temperature = s%logT(im)
       sh%log_pressure = s%logP(im)
       sh%log_radius = s%logR(im)
       sh%luminosity_lsun = s%luminosity_lsun(im)
       sh%log_density = s%logRho(im)
       sh%xa = s%xa(:,im)
       sh%fp = s%fp_rot(im)
       sh%ft = s%ft_rot(im)
       sh%log_temperature_delta = s%log_temperature_delta(im)
       sh%log_pressure_delta = s%log_pressure_delta(im)
       if (s%job%use_mass_accretion.and.s%ctrl%mass_accretion_rate.gt.0.0d0) then
          if (im.ge.envelope_zone_index) then
             sh%log_temperature_delta = s%log_temperature_delta(im)+ &
                  s%delta_log_temperature
             sh%log_pressure_delta = s%log_pressure_delta(im)+ &
                  s%delta_log_pressure
          end if
       end if
       sh%rotational_ke = s%kinetic_energy_rot(im)
       sh%rotational_ke_start = s%kinetic_energy_rot_start(im)
       sh%in_overshoot = im >= s%envelope_overshoot_bottom_index .and. im <= s%envelope_overshoot_top_index
       sh%in_surface_cz = im.ge.s%envelope_cz_bottom_index_start
end subroutine gather_shell_state

! ---------------------------------------------------------------
! The shell's density and convective flag.
subroutine store_shell_structure()
       s%logRho(im) = diag%log_density
       s%convective_flag(im) = diag%convective
end subroutine store_shell_structure

! ---------------------------------------------------------------
! The shell's share of the luminosity channels (dm eps_j/Lsun: each
! nuclear chain in its channel, luminosity_channel_of_chain; the neutrino
! losses in i_lum_neu; eps_grav in i_lum_grav), the extra energy
! sources' eps_extra (erg/g/s, in
! s%eps_extra) and its share dm eps_extra/Lsun in i_lum_extra, and its
! solar-neutrino fluxes: burn%neutrino_flux
! is the flux at the Earth per gram of the shell (10^10 cm^-2 s^-1 g^-1,
! nuclear_energy.f90), so times dm it is the shell's flux.
subroutine add_shell_energy_budget()
       if (diag%burn%burning) then
            if (s%ctrl%calc_neutrinos) then
               do j = 1,n_nu_fluxes
                  s%neutrino_flux_total_e10_per_cm2_s(j) = s%neutrino_flux_total_e10_per_cm2_s(j) + &
                       diag%burn%neutrino_flux(j)*s%dm(im)
               end do
            end if
            do chain = 1, n_chains
               j = luminosity_channel_of_chain(chain)
               if (j == 0) cycle
               s%luminosity_breakdown_lsun(j) = s%luminosity_breakdown_lsun(j) + &
                    (s%dm(im)/s%solar_luminosity_cgs)* &
                    diag%burn%eps_by_chain(chain)
            end do
            s%luminosity_breakdown_lsun(i_lum_neu) = s%luminosity_breakdown_lsun(i_lum_neu) + &
                 (s%dm(im)/s%solar_luminosity_cgs)*diag%burn%eps_neutrino_loss
       end if
       if (diag%has_eps_grav) then
            s%eps_grav(im) = diag%eps_grav
            s%luminosity_breakdown_lsun(i_lum_grav) = s%luminosity_breakdown_lsun(i_lum_grav) + &
                 (s%dm(im)/s%solar_luminosity_cgs)*diag%eps_grav
       end if
       s%eps_extra(im) = diag%eps_extra
       s%luminosity_breakdown_lsun(i_lum_extra) = s%luminosity_breakdown_lsun(i_lum_extra) + &
            (s%dm(im)/s%solar_luminosity_cgs)*diag%eps_extra
end subroutine add_shell_energy_budget

! ---------------------------------------------------------------
! Eliminate the pair (im-1, im), or set the central boundary conditions
! at im = 1; carry cur into prev.
subroutine eliminate_shell_pair()
         if (im.gt.1) then
            im1 = im
            call henyey_eliminate(s, im1, prev, cur)
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
            s%elim_coeff(3,1,1) = cc13*dlnrho_dlnp
            s%elim_coeff(3,2,1) = cc13*dlnrho_dlnt
            s%elim_rhs(3,1) = -cc13*(c4pi3l + diag%log_density - &
                 sh%log_mass) - sh%log_radius
            s%elim_coeff(4,1,1) = -cur%ql_dp
            s%elim_coeff(4,2,1) = -cur%ql_dt
            s%elim_rhs(4,1) = clni*cur%ql - sh%luminosity_lsun
         end if
! Carry this shell's terms into the next pair's elimination.
         prev = cur
end subroutine eliminate_shell_pair

! ---------------------------------------------------------------
! Store the shell's energy generation for output, for every model:
!   eps_total = eps_nuc, the sum of the chains' eps, alpha captures
!     included (erg/g/s; d(p,g) is in pp I, the neutrino losses are not
!     in it),
!   eps_channels(i_eps_*) = the fraction of eps_total of each chain
!     that has a slot there (pp1, pp2, pp3, CNO, triple alpha:
!     eps_channel_of_chain); the alpha captures have none, so where they
!     burn the five fractions add up to less than 1,
!   eps_channels(i_eps_neu) = eps_nu and eps_channels(i_eps_grav) =
!   eps_grav (erg/g/s).
! Below net's burning cutoff (diag%burn%burning false) only eps_grav is
! nonzero.
subroutine store_zone_energy_channels()
         if (.not. diag%burn%burning) then
            s%eps_total(im) = 0.0d0
            s%eps_channels(i_eps_grav,im) = s%eps_grav(im)
            do j = 1,6
               s%eps_channels(j,im) = 0.0d0
           end do
         else
            s%eps_total(im) = 0.0d0
            do chain = 1, n_chains
               s%eps_total(im) = s%eps_total(im) + diag%burn%eps_by_chain(chain)
            end do
            s%eps_channels(i_eps_neu,im) = diag%burn%eps_neutrino_loss
            s%eps_channels(i_eps_grav,im) = s%eps_grav(im)
            if (s%eps_total(im).gt.1.0d-22) then
               energy_sum_inverse = 1.0d0/s%eps_total(im)
            else
               energy_sum_inverse = 0.0d0
            end if
            do chain = 1, n_chains
               j = eps_channel_of_chain(chain)
               if (j /= 0) s%eps_channels(j,im) = &
                    diag%burn%eps_by_chain(chain)*energy_sum_inverse
            end do
         end if
end subroutine store_zone_energy_channels

! ---------------------------------------------------------------
! Store the shell's eos, opacity and gradient values and Gamma1 for output.
subroutine store_zone_thermodynamics()
         s%beta(im) = eos_res%beta
         s%eta(im) = eos_res%electron_degeneracy_parameter
         s%opacity_zone(im) = diag%kap_res(i_kap)
         s%gradr(im) = diag%radiative_gradient
         s%gradT(im) = diag%actual_gradient
         s%grada(im) = eos_res%adiabatic_gradient
         s%conv_vel(im) = diag%convective_velocity
         s%adiabatic_index_gamma1(im) = diag%gamma1

! Disabled with rotation's advection-diffusion transport (2026-10-12),
! their only reader: the copies of dln kappa/dln T, d eps/dln T and
! (eps - eps_nu - eps_grav)/eps for rotation_stability_setup. They come
! back with that transport.
!         if (star%ctrl%rotation_active) then
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
! pulsation output, for every model: the stitched model's pulse point
! set (build_pulse_points) reads them, for the pulse files written with
! the profiles and for the interior GYRE file written every
! pulse_gyre_interval models whatever write_pulse_flag says (MHP 8/25).
! No other output file reads these arrays.
! pulse_deps_dlnrho/dlnt are the derivatives of eps_total (eps_nuc,
! without the thermal neutrino loss), the eps the pulsation files write
! beside them (MESA's convention); zero below the burning cutoff.
subroutine store_pulse_derivatives()
         s%pulse_deps_dlnrho(im) = diag%burn%deps_nuc_dlnrho
         s%pulse_deps_dlnt(im) = diag%burn%deps_nuc_dlnt
         s%pulse_dlnkap_dlnrho(im) = diag%kap_res(i_dlnkap_dlnrho)
         s%pulse_dlnkap_dlnt(im) = diag%kap_res(i_dlnkap_dlnt)
         s%pulse_specific_heat(im) = eos_res%specific_heat_cp
! eos_res%specific_gas_constant is R/mu, so mu = R/specific_gas_constant;
! the electron slot is 1/mu_e straight from the eos.
         if (eos_res%specific_gas_constant > 0.0d0) then
            s%pulse_mean_molecular_weight(im) = gas_constant/eos_res%specific_gas_constant
         else
            s%pulse_mean_molecular_weight(im) = 0.0d0
         end if
         s%pulse_electron_mean_weight_inverse(im) = eos_res%electron_mean_weight_inverse
         s%pulse_dlnrho_dlnt(im) = dlnrho_dlnt
end subroutine store_pulse_derivatives

end subroutine compute_henyey_coefficients

!----------------------------------------------------------------------
! luminosity_channel_sum: the sum (Lsun) of the luminosity channels
! breakdown_lsun (star%luminosity_breakdown_lsun) that are renormalized
! to the surface luminosity: every i_lum_* channel in slot order, the
! extra energy sources' (i_lum_extra) only when a program registered any
! (with_extra_sources). henyey_iterate's per-iteration renormalization
! and observables' post-convergence one both take it, so a new channel
! (n_lum_channels + 1) is in both sums.
pure double precision function luminosity_channel_sum(breakdown_lsun, &
     with_extra_sources)
      use star_info_lib, only: i_lum_extra, n_lum_channels
      double precision, intent(in) :: breakdown_lsun(n_lum_channels)
      logical, intent(in) :: with_extra_sources
      integer :: j

      luminosity_channel_sum = 0.0d0
      do j = 1, n_lum_channels
         if (j == i_lum_extra .and. .not. with_extra_sources) cycle
         luminosity_channel_sum = luminosity_channel_sum + breakdown_lsun(j)
      end do
end function luminosity_channel_sum

!----------------------------------------------------------------------
! Where each nuclear chain (net_lib's chain_*, the slots of
! nuclear_energy_result%eps_by_chain) is stored: its luminosity channel
! (star%luminosity_breakdown_lsun's i_lum_* slot), and its eps channel
! (star%eps_channels' i_eps_* slot; 0 for the alpha captures, which have
! none). A chain without a case here maps to 0 and is left out of the
! breakdown (eps_total still has it): a new chain needs its case.
pure integer function luminosity_channel_of_chain(chain)
      use net_lib, only: chain_pp1, chain_pp2, chain_pp3, chain_cno, &
           chain_alpha_capture, chain_triple_alpha
      use star_info_lib, only: i_lum_pp1, i_lum_pp2, i_lum_pp3, i_lum_cno, &
           i_lum_he_c, i_lum_3alpha
      integer, intent(in) :: chain

      select case (chain)
      case (chain_pp1)
         luminosity_channel_of_chain = i_lum_pp1
      case (chain_pp2)
         luminosity_channel_of_chain = i_lum_pp2
      case (chain_pp3)
         luminosity_channel_of_chain = i_lum_pp3
      case (chain_cno)
         luminosity_channel_of_chain = i_lum_cno
      case (chain_alpha_capture)
         luminosity_channel_of_chain = i_lum_he_c
      case (chain_triple_alpha)
         luminosity_channel_of_chain = i_lum_3alpha
      case default
         luminosity_channel_of_chain = 0
      end select
end function luminosity_channel_of_chain

pure integer function eps_channel_of_chain(chain)
      use net_lib, only: chain_pp1, chain_pp2, chain_pp3, chain_cno, &
           chain_triple_alpha
      use star_info_lib, only: i_eps_pp1, i_eps_pp2, i_eps_pp3, i_eps_cno, &
           i_eps_3alpha
      integer, intent(in) :: chain

      select case (chain)
      case (chain_pp1)
         eps_channel_of_chain = i_eps_pp1
      case (chain_pp2)
         eps_channel_of_chain = i_eps_pp2
      case (chain_pp3)
         eps_channel_of_chain = i_eps_pp3
      case (chain_cno)
         eps_channel_of_chain = i_eps_cno
      case (chain_triple_alpha)
         eps_channel_of_chain = i_eps_3alpha
      case default
         eps_channel_of_chain = 0
      end select
end function eps_channel_of_chain

end module structure_equations
