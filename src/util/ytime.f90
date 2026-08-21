!----------------------------------------------------------------------
! ytime
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ytime.f; only variable names, source form, and comment
! style were updated.
!
! Stars with a helium luminosity have a timestep limit based on the
! time required to burn atime(4) of Y at the maximum in temperature or
! the fraction atime(5) of Y at this point.
! Test for helium flash burning; high central density is taken as a
! sign that the model is a giant undergoing a He flash rather than a
! horizontal branch star.
subroutine ytime(energy_gen_terms, composition, log_density, luminosity, &
     enclosed_mass, log_temperature, convective_core_edge_zone, &
     num_points, helium_dt, rate_pp, rate_he3_he3, rate_he3_he4, &
     rate_c12_p, rate_c13_p, rate_n14_p, rate_o16_p, rate_c13_alpha, &
     rate_zero9, rate_c12_alpha, rate_n14_alpha, rate_triple_alpha, &
     rate_zero13, frac_c12_alpha, frac_be7_electron, h_shell_zone_begin)

      use nuclear_lib
      use engeb_diag_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(out) :: energy_gen_terms(6)
      double precision, intent(in) :: composition(15,json)
      double precision, intent(in) :: log_density(json), luminosity(json), &
           enclosed_mass(json), log_temperature(json)
      integer, intent(in) :: convective_core_edge_zone, num_points
! helium_dt: intent(inout), not intent(out) -- in the "not a helium
! flash, core Y below atime(1)" branch below, the original reads
! DELTSY on the right-hand side before this call has assigned it
! (relying on whatever value the caller's actual argument already
! held). This is preserved exactly, not fixed.
      double precision, intent(inout) :: helium_dt
      double precision, intent(out) :: rate_pp(json), rate_he3_he3(json), &
           rate_he3_he4(json), rate_c12_p(json), rate_c13_p(json), &
           rate_n14_p(json), rate_o16_p(json), rate_c13_alpha(json), &
           rate_zero9(json), rate_c12_alpha(json), rate_n14_alpha(json), &
           rate_triple_alpha(json), rate_zero13(json)
      double precision, intent(out) :: frac_c12_alpha(json), &
           frac_be7_electron(json)
      integer, intent(in) :: h_shell_zone_begin







      save

      double precision :: max_temp, local_log_density, local_log_temperature
      integer :: max_temp_zone, zone_idx, engeb_zone
      double precision :: hydrogen_fraction, helium_fraction, &
           metal_fraction, he3_fraction, c12_fraction, c13_fraction, &
           n14_fraction, n15_fraction, o16_fraction, o17_fraction, &
           o18_fraction, h2_fraction, li6_fraction, li7_fraction, &
           be9_fraction
      double precision :: energy_gen_1, energy_gen_2, energy_gen_3, &
           energy_gen_4, energy_gen_5, qed_correction, qet_correction, &
           total_energy_gen, total_nuclear_energy_gen, &
           dead_alpha_capture_energy
      double precision :: core_helium_fraction

! DATA Q3A,QCA/5.85D17,1.7276D18/
      if (log_density(1).ge.5.0d0) then
!     search for temperature maximum
       max_temp = 0.0d0
       max_temp_zone = 1
       do 210 zone_idx = 1,num_points
          if (log_temperature(zone_idx).gt.max_temp) then
             max_temp_zone = zone_idx
             max_temp = log_temperature(zone_idx)
          endif
  210    continue
!     calculate helium burning rate at tmax
       local_log_density = log_density(max_temp_zone)
       local_log_temperature = log_temperature(max_temp_zone)
       hydrogen_fraction = composition(1,max_temp_zone)
       helium_fraction = composition(2,max_temp_zone)
       metal_fraction = composition(3,max_temp_zone)
       he3_fraction = composition(4,max_temp_zone)
       c12_fraction = composition(5,max_temp_zone)
       c13_fraction = composition(6,max_temp_zone)
       n14_fraction = composition(7,max_temp_zone)
       n15_fraction = composition(8,max_temp_zone)
       o16_fraction = composition(9,max_temp_zone)
       o17_fraction = composition(10,max_temp_zone)
       o18_fraction = composition(11,max_temp_zone)
       if(use_extended_composition) then
          h2_fraction = composition(12,max_temp_zone)
          li6_fraction = composition(13,max_temp_zone)
          li7_fraction = composition(14,max_temp_zone)
          be9_fraction = composition(15,max_temp_zone)
       endif
       engeb_zone = max_temp_zone
         call engeb(energy_gen_1,energy_gen_2,energy_gen_3,energy_gen_4, &
              energy_gen_5,qed_correction,qet_correction,total_energy_gen, &
              local_log_density,local_log_temperature,hydrogen_fraction, &
              helium_fraction,he3_fraction,c12_fraction,c13_fraction, &
              n14_fraction,o16_fraction,o18_fraction,h2_fraction, &
              engeb_zone,rate_pp,rate_he3_he3,rate_he3_he4,rate_c12_p, &
              rate_c13_p,rate_n14_p,rate_o16_p,rate_c13_alpha,rate_zero9, &
              rate_c12_alpha,rate_n14_alpha,rate_triple_alpha,rate_zero13, &
              frac_c12_alpha,frac_be7_electron)
         total_nuclear_energy_gen = total_energy_gen
         energy_gen_terms(1) = energy_gen_1
         energy_gen_terms(2) = energy_gen_2
         energy_gen_terms(3) = energy_gen_3
         energy_gen_terms(4) = energy_gen_4
         energy_gen_terms(5) = energy_gen_5
         energy_gen_terms(6) = engeb_diag%neutrino_loss_rate
         dead_alpha_capture_energy = engeb_diag%alpha_capture_energy
       if(energy_gen_terms(5).lt.1.d-22) energy_gen_terms(5) = 1.d-22
       helium_dt=1.0d15/energy_gen_terms(5)
      else
!  set limits on he burning for non-helium flash stars
!  core_helium_fraction: computed as 1 - Z(core); in a He-burning core
!  X is normally negligible there, so this approximates the core Y.
       core_helium_fraction = 1d0 - composition(3,convective_core_edge_zone)
       if(core_helium_fraction.ge.atime(1)) then

            helium_dt = min(atime(4),atime(5)*core_helium_fraction)

          helium_dt = (5.85d17/solar_luminosity_cgs)*helium_dt* &
               (enclosed_mass(convective_core_edge_zone)/luminosity(convective_core_edge_zone))

            else

            helium_dt = (5.85d17/solar_luminosity_cgs)*helium_dt* &
                 (enclosed_mass(h_shell_zone_begin-1)/luminosity(h_shell_zone_begin-1))

         endif

      endif

      return
end subroutine ytime
