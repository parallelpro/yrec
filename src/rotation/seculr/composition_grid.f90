!----------------------------------------------------------------------
! mixgrid
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mixgrid.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Builds the equally-spaced-grid quantities needed by diffuse_composition.f90's
! diffusive composition solve: given the first/last unstable zones of
! a region (zone_begin/zone_end), calls equal_spaced_grid to lay down the
! equally spaced coordinate rot_scr%chi, then computes the equally spaced grid
! masses (equally_spaced_mass) and the geometrically weighted
! diffusion coefficients (equally_spaced_diffusion_coeff) at the zone
! edges, including the Jacobian factor for the transformation from
! radius to the rot_scr%chi coordinate.
subroutine composition_grid(diffusion_coeff, log_density, luminosity_lsun, &
     log_pressure, log_radius, log_mass, enclosed_mass, shell_mass, &
     log_total_mass, zone_begin, zone_end, convective_flag, num_zones, &
     equally_spaced_diffusion_coeff, equally_spaced_mass, &
     single_interface_flag)
      use rotation_scratch_lib
      use equal_grid_lib
      use star_info_lib, only: json
      use phys_const_lib, only: ln10
      use math_lib
      implicit none

      double precision, intent(in) :: diffusion_coeff(json), &
           log_density(json), luminosity_lsun(json), log_pressure(json), &
           log_radius(json), log_mass(json), enclosed_mass(json), &
           shell_mass(json)
      double precision, intent(in) :: log_total_mass
      integer, intent(in) :: zone_begin, zone_end
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(out) :: equally_spaced_diffusion_coeff(json), &
           equally_spaced_mass(json)
      logical, intent(out) :: single_interface_flag
      integer :: idx, search_idx, ntab, ntabb
      double precision :: em_top, em_bot

! FLAG THE SPECIAL CASE OF A SINGLE UNSTABLE INTERFACE AND EXIT
      if (zone_end - zone_begin.le.1) then
         single_interface_flag = .true.
         return
      else
         single_interface_flag = .false.
      end if
! DEFINE A GRID OF EQUALLY SPACED POINTS.
      call equal_spaced_grid(luminosity_lsun, log_pressure, log_mass, zone_begin, &
           zone_end, num_zones)
! EQUAL_SPACED_GRID HAS DEFINED A SET OF CO-ORDINATES (CHI) AND EQUALLY
! SPACED MASS POINTS.  NOW FIND THE OTHER QUANTITIES OF INTEREST AT ZONE
! CENTERS:
! TOTAL ZONE MASSES
! INTERMEDIATE POINTS
      do idx = 2, rot_scr%ntot-1
         equally_spaced_mass(idx) = 0.5d0*(rot_scr%eq_enclosed_mass(idx+1) - rot_scr%eq_enclosed_mass(idx-1))
      end do
! SPECIAL TREATMENT OF THE BOUNDARIES; CAN BE CONVECTIVE.
! IF CONVECTIVE SUM OVER ALL SHELLS.  CARE IS NEEDED TO DO BOOK-KEEPING
! PROPERLY AT THE EDGES - TOP IS HALFWAY TO EQUALLY SPACED POINT, NOT
! HALFWAY TO EDGE OF UNEQUALLY SPACED ORIGINAL SET OF POINTS.
!
! CENTER
      em_top = 0.5d0*(rot_scr%eq_enclosed_mass(2) + rot_scr%eq_enclosed_mass(1))
      if (zone_begin.gt.1) then
         em_bot = 0.5d0*(enclosed_mass(zone_begin) + &
              enclosed_mass(zone_begin-1))
      else
         em_bot = 0.0d0
      end if
      equally_spaced_mass(1) = em_top - em_bot
! (INCLUDE EVERY CONVECTIVE SHELL BELOW ZONE_BEGIN.)
      if (zone_begin.gt.1) then
         do search_idx = zone_begin-1, 1, -1
            if (.not.convective_flag(search_idx)) exit
            equally_spaced_mass(1) = equally_spaced_mass(1) + &
                 shell_mass(search_idx)
         end do
      end if
! SURFACE
      em_bot = 0.5d0*(rot_scr%eq_enclosed_mass(rot_scr%ntot) + rot_scr%eq_enclosed_mass(rot_scr%ntot-1))
      if (zone_end.lt.num_zones) then
         em_top = 0.5d0*(enclosed_mass(zone_end) + &
              enclosed_mass(zone_end+1))
      else
         em_top = exp(ln10*log_total_mass)
      end if
      equally_spaced_mass(rot_scr%ntot) = em_top - em_bot
! (INCLUDE EVERY CONVECTIVE SHELL ABOVE ZONE_END.)
      if (zone_end.lt.num_zones) then
         do search_idx = zone_end+1, num_zones
            if (.not.convective_flag(search_idx)) exit
            equally_spaced_mass(rot_scr%ntot) = equally_spaced_mass(rot_scr%ntot) + &
                 shell_mass(search_idx)
         end do
      end if
! NOW SOLVE FOR QUANTITIES NEEDED AT THE ZONE EDGES.  THESE ARE
! RELATED TO THE DIFFUSION COEFFICIENTS.  UNLIKE THE EQUALLY SPACED
! GRID IN R, WE NEED TO INCLUDE A JACOBIAN TERM FOR THE TRANSFORMATION
! OF VARIABLES.
      ntab = zone_end - zone_begin + 1
      call edge_grid_abscissae(rot_scr%chi, rot_scr%echi, rot_scr%dchi, ntab, &
           rot_scr%ntot, rot_scr%xtab, rot_scr%xval, ntabb)
! DIFFUSION COEFFICIENT FOR MIXING - ASSUME CONSTANT BELOW
! BOTTOM INTERFACE OR ABOVE TOP INTERFACE
      call interp_edge_coeff(diffusion_coeff, zone_begin, zone_end, ntab, ntabb, &
           rot_scr%ntot, rot_scr%xval, rot_scr%xtab, rot_scr%ytab, &
           equally_spaced_diffusion_coeff)
! PRODUCT OF RHO R^2 BY D CHI/DR (LOG10, ON THE EQUAL GRID IN YVAL)
      call dchi_dr_jacobian(log_density, log_radius, log_mass, log_pressure, &
           enclosed_mass, mix_scr%epsm, luminosity_lsun(num_zones), zone_begin, &
           ntab, rot_scr%ntot, rot_scr%chi, rot_scr%xval, rot_scr%xtab, &
           rot_scr%ytab, rot_scr%yval)
! NOW ADD MULTIPLICATIVE FACTORS TO DIFFUSION COEFFICIENTS
      call multiply_by_exp10(equally_spaced_diffusion_coeff, rot_scr%yval, rot_scr%ntot)
      return
end subroutine composition_grid
