!----------------------------------------------------------------------
! ytime
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ytime.f; only variable names, source form, and comment
! style were updated -- except the He-shell branch (2026, see below).
!
! Stars with a helium luminosity have a timestep limit based on the
! time required to burn atime(4) of Y at the maximum in temperature or
! the fraction atime(5) of Y at this point.
! Test for helium flash burning; high central density is taken as a
! sign that the model is a giant undergoing a He flash rather than a
! horizontal branch star.
subroutine timestep_limit_heburn(composition, log_density, luminosity, &
     enclosed_mass, log_temperature, convective_core_edge_zone, &
     num_points, helium_dt, h_shell_zone_begin)
      use net_lib
      use star_info_lib, only: star, json
      use burn_lib
      implicit none

      double precision, intent(in) :: composition(15,json)
      double precision, intent(in) :: log_density(json), luminosity(json), &
           enclosed_mass(json), log_temperature(json)
      integer, intent(in) :: convective_core_edge_zone, num_points
! helium_dt: assigned on every path (2026 -- the He-shell branch no
! longer reads the incoming value; see the note there).
      double precision, intent(out) :: helium_dt
      integer, intent(in) :: h_shell_zone_begin
      double precision :: max_temp, local_log_density, local_log_temperature
      integer :: max_temp_zone, zone_idx, engeb_zone, shell_zone
      double precision :: hydrogen_fraction, helium_fraction, &
           he3_fraction, c12_fraction, c13_fraction, &
           n14_fraction, o16_fraction, o18_fraction, h2_fraction
! only the helium-burning term (energy_gen_5) of the engeb output is
! used here; the other terms are received and discarded.
      double precision :: energy_gen_1, energy_gen_2, energy_gen_3, &
           energy_gen_4, energy_gen_5, qed_correction, qet_correction, &
           total_energy_gen, helium_energy_gen
      double precision :: core_helium_fraction
! energy released per gram of helium burned (erg/g)
      double precision, parameter :: he_burn_energy_per_gram = 5.85d17
! dt_unlimited: "no limit from this criterion" sentinel (seconds).
      double precision, parameter :: dt_unlimited = 1.0d20

      if (log_density(1).ge.5.0d0) then
!     search for temperature maximum
       max_temp = 0.0d0
       max_temp_zone = 1
       do zone_idx = 1,num_points
          if (log_temperature(zone_idx).gt.max_temp) then
             max_temp_zone = zone_idx
             max_temp = log_temperature(zone_idx)
          endif
       end do
!     calculate helium burning rate at tmax
       local_log_density = log_density(max_temp_zone)
       local_log_temperature = log_temperature(max_temp_zone)
       hydrogen_fraction = composition(1,max_temp_zone)
       helium_fraction = composition(2,max_temp_zone)
       he3_fraction = composition(4,max_temp_zone)
       c12_fraction = composition(5,max_temp_zone)
       c13_fraction = composition(6,max_temp_zone)
       n14_fraction = composition(7,max_temp_zone)
       o16_fraction = composition(9,max_temp_zone)
       o18_fraction = composition(11,max_temp_zone)
! h2_fraction is left at its -finit-local-zero value (0) when the
! extended composition is off, as in the original.
       if(star%job%use_extended_composition) then
          h2_fraction = composition(12,max_temp_zone)
       endif
       engeb_zone = max_temp_zone
         call engeb(energy_gen_1,energy_gen_2,energy_gen_3,energy_gen_4, &
              energy_gen_5,qed_correction,qet_correction,total_energy_gen, &
              local_log_density,local_log_temperature,hydrogen_fraction, &
              helium_fraction,he3_fraction,c12_fraction,c13_fraction, &
              n14_fraction,o16_fraction,o18_fraction,h2_fraction, &
              engeb_zone)
         helium_energy_gen = energy_gen_5
       if(helium_energy_gen.lt.1.d-22) helium_energy_gen = 1.d-22
       helium_dt=1.0d15/helium_energy_gen
      else
!  set limits on he burning for non-helium flash stars
!  core_helium_fraction: computed as 1 - Z(core); in a He-burning core
!  X is normally negligible there, so this approximates the core Y.
       core_helium_fraction = 1d0 - composition(3,convective_core_edge_zone)
       if(core_helium_fraction.ge.star%ctrl%atime(1)) then
            helium_dt = min(star%ctrl%atime(4),star%ctrl%atime(5)*core_helium_fraction)
            helium_dt = (he_burn_energy_per_gram/star%solar_luminosity_cgs)*helium_dt* &
               (enclosed_mass(convective_core_edge_zone)/luminosity(convective_core_edge_zone))
         else
! 2026 (bugsweep sec-10/11): core helium exhausted -- the He-shell
! branch. The original (ytime.f) read DELTSY on the right-hand side
! here, i.e. the PREVIOUS call's timestep (kept alive only by the
! caller's blanket SAVE), and multiplied it by M/L again: a
! dimensionless-nonsense recurrence that gave ~0 s on the first
! model after exhaustion and, once the SAVE was gone, exactly 0 ->
! dt = 0 -> NaN (the run_from_zahb_to_tahb failure at Y_c < atime(1)).
! Compute a fresh limit from the He-shell controls instead, mirroring
! the H-shell branch of timestep_limit_hburn: atime(14) (time_dy_total,
! Msun of He burned per step) and atime(12) (time_dy_shell, fraction
! of Y burned at the shell) -- both were mapped by map_user_inputs but
! never consumed until now. The reference point is the shell just
! below the H-burning shell (the He-rich core), as before.
            shell_zone = max(h_shell_zone_begin-1, 1)
            helium_dt = dt_unlimited
            if (luminosity(shell_zone) .gt. 0.0d0) then
               if (star%ctrl%atime(14) .gt. 0.0d0) helium_dt = min(helium_dt, &
                    (he_burn_energy_per_gram/star%solar_luminosity_cgs)*star%ctrl%atime(14)* &
                    (star%solar_mass_cgs/luminosity(shell_zone)))
               if (star%ctrl%atime(12) .gt. 0.0d0) helium_dt = min(helium_dt, &
                    (he_burn_energy_per_gram/star%solar_luminosity_cgs)*star%ctrl%atime(12)* &
                    composition(2,shell_zone)* &
                    (enclosed_mass(shell_zone)/luminosity(shell_zone)))
            end if
         endif
      endif
      return
end subroutine timestep_limit_heburn
