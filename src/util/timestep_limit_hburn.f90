!----------------------------------------------------------------------
! xtime
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original xtime.f; only variable names, source form, and comment
! style were updated.
!
! This subroutine finds the timestep based on hydrogen burning.
! For stars with central X > atime(1) the timestep is the time needed
! to burn the minimum of atime(2) of X at the center or the fraction
! atime(3) of the central X.
! Stars with central X < atime(1) are considered by the program to
! have a hydrogen shell burning source. The timestep is the minimum of
! the time required to burn atime(7) of X at the shell midpoint or to
! burn the mass fraction atime(6) of X in the entire star.
subroutine timestep_limit_hburn(log_density, composition, luminosity, enclosed_mass, &
     shell_mass, log_temperature, hydrogen_luminosity, &
     convective_core_edge_zone, h_shell_midpoint_zone, num_points, &
     hydrogen_dt, rate_pp, rate_he3_he3, rate_he3_he4, rate_c12_p, &
     rate_c13_p, rate_n14_p, rate_o16_p, rate_c13_alpha, rate_zero9, &
     rate_c12_alpha, rate_n14_alpha, rate_triple_alpha, rate_zero13, &
     frac_c12_alpha, frac_be7_electron)
      use star_info_lib, only: star, json

      use phys_const_lib
      use net_lib
      use burn_lib
      implicit none

      double precision, intent(in) :: log_density(json)
      double precision, intent(in) :: composition(15,json)
      double precision, intent(in) :: luminosity(json)
      double precision, intent(in) :: enclosed_mass(json)
      double precision, intent(in) :: shell_mass(json)
      double precision, intent(in) :: log_temperature(json)
      double precision, intent(in) :: hydrogen_luminosity
      integer, intent(in) :: convective_core_edge_zone, &
           h_shell_midpoint_zone, num_points
      double precision, intent(out) :: hydrogen_dt
      double precision, intent(out) :: rate_pp(json), rate_he3_he3(json), &
           rate_he3_he4(json), rate_c12_p(json), rate_c13_p(json), &
           rate_n14_p(json), rate_o16_p(json), rate_c13_alpha(json), &
           rate_zero9(json), rate_c12_alpha(json), rate_n14_alpha(json), &
           rate_triple_alpha(json), rate_zero13(json)
      double precision, intent(out) :: frac_c12_alpha(json), &
           frac_be7_electron(json)
! delta_x: hydrogen mass-fraction change budget for the current
! branch (core-burning limit in the first branch, whole-star shell-
! burning limit in the second -- a single reused scratch variable in
! the original, exactly as here).
      double precision :: delta_x
! local_log_density/local_log_temperature: log10(rho)/log10(T) at the
! shell midpoint, fed to rates().
      double precision :: local_log_density, local_log_temperature
      double precision :: hydrogen_fraction, helium_fraction, &
           metal_fraction, he3_fraction, c12_fraction, c13_fraction, &
           n14_fraction, n15_fraction, o16_fraction, o17_fraction, &
           o18_fraction, h2_fraction, li6_fraction, li7_fraction, &
           be9_fraction
      integer :: zone_begin, zone_end
      double precision :: dc_dt, do_dt, dx_dt, dy_dt, shell_dt_x_depletion

!  h-core burning time criterion
! **note that for a convective core, the timestep is based on the time
!   needed to burn the given fraction of hydrogen in the core and not
!   just in the central shell.
      if(composition(1,1).ge.star%ctrl%atime(1)) then
       delta_x = min(star%ctrl%atime(2),star%ctrl%atime(3)*composition(1,convective_core_edge_zone))
       hydrogen_dt =(6.00d18/star%solar_luminosity_cgs)* &
            (enclosed_mass(convective_core_edge_zone)/luminosity(convective_core_edge_zone))*delta_x
       return
      endif
!  h-shell burning criterion
!  limit total mass of hydrogen burned.
      delta_x = star%ctrl%atime(6)*composition(1,num_points)*(star%solar_mass_cgs/star%solar_luminosity_cgs)
      hydrogen_dt = 6.00d18*delta_x/hydrogen_luminosity
!  limit x-depletion at shell mid-point.
!  call nuclear reaction sr's to find dxdt at the shell midpoint.
      zone_begin=h_shell_midpoint_zone
      zone_end=h_shell_midpoint_zone
         local_log_density = log_density(zone_end)
         local_log_temperature = log_temperature(zone_end)
         hydrogen_fraction = composition(1,zone_end)
         helium_fraction = composition(2,zone_end)
         metal_fraction = composition(3,zone_end)
         he3_fraction = composition(4,zone_end)
         c12_fraction = composition(5,zone_end)
         c13_fraction = composition(6,zone_end)
         n14_fraction = composition(7,zone_end)
         n15_fraction = composition(8,zone_end)
         o16_fraction = composition(9,zone_end)
         o17_fraction = composition(10,zone_end)
         o18_fraction = composition(11,zone_end)
         h2_fraction = composition(12,zone_end)
         li6_fraction = composition(13,zone_end)
         li7_fraction = composition(14,zone_end)
        be9_fraction = composition(15,zone_end)
!  setup nuclear energy terms
      call rates(local_log_density,local_log_temperature,hydrogen_fraction, &
           helium_fraction,he3_fraction,c12_fraction,c13_fraction,n14_fraction, &
           o16_fraction,o18_fraction,zone_end,rate_pp,rate_he3_he3,rate_he3_he4, &
           rate_c12_p,rate_c13_p,rate_n14_p,rate_o16_p,rate_c13_alpha,rate_zero9, &
           rate_c12_alpha,rate_n14_alpha,rate_triple_alpha,rate_zero13, &
           frac_c12_alpha,frac_be7_electron)
      call eqburn(rate_pp,rate_he3_he3,rate_he3_he4,rate_c12_p,rate_c13_p, &
           rate_n14_p,rate_o16_p,rate_c12_alpha,rate_triple_alpha,shell_mass, &
           log_temperature,zone_begin,zone_end,dc_dt,do_dt,dx_dt,dy_dt, &
           c12_fraction,o16_fraction,hydrogen_fraction,metal_fraction)
      if(dx_dt.lt.0.0d0 .and. star%ctrl%atime(7).gt.0.0d0) then
         shell_dt_x_depletion = abs(seconds_per_year*1.0d9*star%ctrl%atime(7)/dx_dt)
         hydrogen_dt = min(hydrogen_dt,shell_dt_x_depletion)
      endif
      return
end subroutine timestep_limit_hburn
