!----------------------------------------------------------------------
! solve_composition (formerly kemcom)
!----------------------------------------------------------------------
! New abundances at the end of the timestep for one radiative shell
! (zone_begin = zone_end) or one convection zone (zone_begin < zone_end,
! homogenized), by an implicit burning step of the network in
! net/network_lib.f90. This routine does only the model bookkeeping
! around it; the reactions, the rate law and the Newton iteration live
! in network_lib (2026-10: they replaced the hand-written 7-species
! equations and Jacobian that stood here).
!
!  log_temperature     run of log10 T; shells below
!                      ctrl%nuclear_logT_cutoffs(1) are only mixed
!  rate_by_zone(rr_*,k) net_lib's rate vectors, per Gyr per amu
!  shell_mass          the mass weights for a convection zone
!  composition         out: the new abundances (star%xa layout); the
!                      solved shell's previous neighbour is the first
!                      guess, so shells are solved in increasing order
!  timestep_years      the step
! The start-of-step abundances are star%xa_start. Z changes only once
! hydrogen is exhausted (X < 5e-7 at the top of the region), as
! 1 - X - Y - He3.
subroutine solve_composition(log_temperature, zone_begin, zone_end, &
     rate_by_zone, shell_mass, composition, timestep_years, ierr)

      use star_info_lib, only: star, json, i_h1, i_he3, i_he4, i_metals, &
           i_c12, i_c13, i_n14, i_o16, i_o18, n_species_basic
      use net_lib, only: num_rate_rows
      use network_lib, only: implicit_burn_step, n_solved, n_network, &
           s_h1, s_he3, s_he4
      use luout_lib
      implicit none

      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: zone_begin, zone_end
      double precision, intent(in) :: rate_by_zone(num_rate_rows, json)
      double precision, intent(in) :: shell_mass(json)
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: timestep_years
      integer, intent(out) :: ierr

! star%xa slot of each network species (network_lib's s_* order)
      integer, parameter :: xa_slot(n_network) = &
           [i_h1, i_he3, i_he4, i_c12, i_c13, i_n14, i_o16, i_o18]
      double precision :: region_start(n_species_basic), rate_vec(num_rate_rows)
      double precision :: x_start(n_network), x(n_network), region_mass, metals
      integer :: k, guess_zone

      ierr = 0

! start-of-step abundances and rates of the region (mass-weighted
! averages over a convection zone)
      if (zone_begin /= zone_end) then
         region_mass = 0.0d0
         region_start = 0.0d0
         rate_vec = 0.0d0
         do k = zone_begin, zone_end
            region_mass = region_mass + shell_mass(k)
            region_start = region_start + star%xa_start(1:n_species_basic,k)*shell_mass(k)
            rate_vec = rate_vec + shell_mass(k)*rate_by_zone(:,k)
         end do
         region_start = region_start/region_mass
         rate_vec = rate_vec/region_mass
      else
         region_start = star%xa_start(1:n_species_basic,zone_begin)
         rate_vec = rate_by_zone(:,zone_begin)
      end if

! too cool to burn: only homogenize
      if (log_temperature(zone_begin) < star%ctrl%nuclear_logT_cutoffs(1)) then
         do k = zone_begin, zone_end
            composition(1:n_species_basic,k) = region_start
         end do
         return
      end if

! first guess: the shell just solved below (or this one, at the centre)
      guess_zone = max(zone_begin - 1, 1)
      x_start = region_start(xa_slot)
      x(1:n_solved) = composition(xa_slot(1:n_solved), guess_zone)
      x(n_solved+1:) = x_start(n_solved+1:)

      call implicit_burn_step(x_start, rate_vec, timestep_years*1.0d-9, &
           star%ctrl%min_abundance, star%ctrl%absolute_tolerance, &
           star%ctrl%relative_tolerance, star%ctrl%max_burn_iterations, x, ierr)
      if (ierr /= 0) then
         write (run_log_unit,1000) zone_begin
 1000    format(1X,39('>'),40('<')/1X,'ERROR IN SUBROUTINE KEMCOM'/ &
              1X,'UNABLE TO SOLVE FOR NEW ABUNDANCES IN SHELL',I4/1X, &
              'RUN STOPPED AFTER 50 ATTEMPTS')
         return
      end if

      if (star%xa_start(i_h1,zone_end) < 5.0d-7) then
         metals = 1.0d0 - x(s_h1) - x(s_he3) - x(s_he4)
      else
         metals = region_start(i_metals)
      end if
      do k = zone_begin, zone_end
         composition(xa_slot,k) = x
         composition(i_metals,k) = metals
      end do
end subroutine solve_composition
