!----------------------------------------------------------------------
! envelope_refit_lib
!----------------------------------------------------------------------
! The per-shell fill that the two envelope refits share (2026 W3).
! After atm_get has integrated a fresh envelope from the surface down
! past the current fitting point and the caller has shifted env_struct
! so that its point 1 lies just below the last interior shell,
! append_envelope_points appends the envelope points that lie below
! the new fitting mass star%senv to the model grid as new shells,
! interpolates the last one onto the fitting mass exactly, and
! advances num_zones.
!
! Callers: core/rebuild_envelope.f90 (fitting point moved during the
! run) and core/read_starting_model.f90 (requested envelope mass
! shallower than the starting model's). Everything the two callers do
! before this loop (surface L/Teff/g, the env_struct shift) and after
! it (shell masses, rotation, the prints) differs and stays in the
! caller. Inside the loop they differ in exactly one test -- whether
! an envelope point sitting exactly at the fitting mass is taken as a
! shell (read_starting_model: .le.) or interpolated (rebuild_envelope:
! .lt.) -- carried by accept_point_at_fit_mass. Each caller prints
! its own message when ierr comes back nonzero (degenerate mass
! interval).
module envelope_refit_lib
      implicit none
      private
      public :: append_envelope_points

contains

subroutine append_envelope_points(accept_point_at_fit_mass, species_end_index, &
     composition, log_density, log_luminosity, log_pressure, log_radius, &
     log_mass, log_temperature, convective_flag, num_zones, ierr)
      use star_info_lib, only: star, json, n_species_extended, &
           i_h1, i_he4, i_metals, i_he3
      use envstruct_lib
      logical, intent(in) :: accept_point_at_fit_mass
      integer, intent(in) :: species_end_index
      double precision, intent(inout) :: composition(n_species_extended,json), &
           log_density(json), log_luminosity(json), log_pressure(json), &
           log_radius(json), log_mass(json), log_temperature(json)
      logical, intent(inout) :: convective_flag(json)
      integer, intent(inout) :: num_zones
      integer, intent(out) :: ierr

      integer :: zone_index, env_point_index, k
      logical :: point_below_fit_mass
      double precision :: mass_interp_x0, mass_interp_x1, mass_interp_x2
      double precision :: interp_fraction

      ierr = 0
! ASSIGN NEW POINTS
      do zone_index = num_zones+1,num_zones+env_struct%num_env_points
         env_point_index = zone_index-num_zones
! LUMINOSITY ASSUMED CONSTANT
         log_luminosity(zone_index) = log_luminosity(num_zones)
! INCLUDE NEW POINTS UP TO THE DIFFERENT DESIRED FITTING POINT
         if (accept_point_at_fit_mass) then
            point_below_fit_mass = env_struct%env_log10_mass(env_point_index).le.star%senv
         else
            point_below_fit_mass = env_struct%env_log10_mass(env_point_index).lt.star%senv
         end if
         if(point_below_fit_mass)then
            log_density(zone_index) = env_struct%env_log10_density(env_point_index)
            log_pressure(zone_index) = env_struct%env_log10_pressure(env_point_index)
            log_radius(zone_index) = env_struct%env_log10_radius(env_point_index)
            log_mass(zone_index) = env_struct%env_log10_mass(env_point_index) + star%stotal
            log_temperature(zone_index) = env_struct%env_log10_temperature(env_point_index)
            composition(i_h1,zone_index) = env_struct%env_hydrogen_fraction(env_point_index)
            composition(i_metals,zone_index) = env_struct%env_metal_fraction(env_point_index)
            do k = i_he3,species_end_index
               composition(k,zone_index) = composition(k,num_zones)
            end do
            composition(i_he4,zone_index)=1.0D0-composition(i_h1,zone_index)- &
                 composition(i_metals,zone_index)-composition(i_he3,zone_index)
            convective_flag(zone_index) = env_struct%env_convective_flag(env_point_index)
         else
! POINTS BEYOND THIS ARE ABOVE THE NEW DESIRED FITTING POINT;
! INTERPOLATE LINEARLY, SET NEW NUMBER OF TOTAL POINTS, AND EXIT
            if(env_point_index.eq.1)then
! INTERPOLATE BETWEEN THE LAST INTERIOR POINT AND THE FIRST ENVELOPE POINT
               mass_interp_x0 = log_mass(num_zones)
               mass_interp_x1 = star%stotal + star%senv
               mass_interp_x2 = env_struct%env_log10_mass(env_point_index) + star%stotal
               if(mass_interp_x2-mass_interp_x0.lt.1.0D-14) then
                  ierr = 1
                  return
               end if
               interp_fraction = (mass_interp_x1-mass_interp_x0)/ &
                    (mass_interp_x2-mass_interp_x0)
               log_density(zone_index) = log_density(num_zones)+interp_fraction* &
                    (env_struct%env_log10_density(env_point_index)-log_density(num_zones))
               log_pressure(zone_index) = log_pressure(num_zones)+interp_fraction* &
                    (env_struct%env_log10_pressure(env_point_index)-log_pressure(num_zones))
               log_radius(zone_index) = log_radius(num_zones)+interp_fraction* &
                    (env_struct%env_log10_radius(env_point_index)-log_radius(num_zones))
               log_mass(zone_index) = mass_interp_x1
               log_temperature(zone_index) = log_temperature(num_zones)+interp_fraction* &
                    (env_struct%env_log10_temperature(env_point_index)-log_temperature(num_zones))
               composition(i_h1,zone_index) = composition(i_h1,num_zones)+interp_fraction* &
                    (composition(i_h1,num_zones)-env_struct%env_hydrogen_fraction(env_point_index))
               composition(i_metals,zone_index) = composition(i_metals,num_zones)+interp_fraction* &
                    (composition(i_metals,num_zones)-env_struct%env_metal_fraction(env_point_index))
               do k = i_he3,species_end_index
                  composition(k,zone_index) = composition(k,num_zones)
               end do
               composition(i_he4,zone_index)=1.0D0-composition(i_h1,zone_index)- &
                    composition(i_metals,zone_index)-composition(i_he3,zone_index)
               if(env_struct%env_convective_flag(env_point_index).or.convective_flag(num_zones))then
                  convective_flag(zone_index) = .true.
               else
                  convective_flag(zone_index) = .false.
               endif
            else
! INTERPOLATE BETWEEN THE LAST 2 ENVELOPE POINTS
               mass_interp_x0 = env_struct%env_log10_mass(env_point_index-1) + star%stotal
               mass_interp_x1 = star%stotal + star%senv
               mass_interp_x2 = env_struct%env_log10_mass(env_point_index) + star%stotal
               if(mass_interp_x2-mass_interp_x0.lt.1.0D-14) then
                  ierr = 1
                  return
               end if
               interp_fraction = (mass_interp_x1-mass_interp_x0)/ &
                    (mass_interp_x2-mass_interp_x0)
               log_density(zone_index) = env_struct%env_log10_density(env_point_index-1)+interp_fraction* &
                    (env_struct%env_log10_density(env_point_index)-env_struct%env_log10_density(env_point_index-1))
               log_pressure(zone_index) = env_struct%env_log10_pressure(env_point_index-1)+interp_fraction* &
                    (env_struct%env_log10_pressure(env_point_index)-env_struct%env_log10_pressure(env_point_index-1))
               log_radius(zone_index) = env_struct%env_log10_radius(env_point_index-1)+interp_fraction* &
                    (env_struct%env_log10_radius(env_point_index)-env_struct%env_log10_radius(env_point_index-1))
               log_mass(zone_index) = mass_interp_x1
               log_temperature(zone_index) = env_struct%env_log10_temperature(env_point_index-1)+interp_fraction* &
                    (env_struct%env_log10_temperature(env_point_index)-env_struct%env_log10_temperature(env_point_index-1))
               composition(i_h1,zone_index) = env_struct%env_hydrogen_fraction(env_point_index-1)+interp_fraction* &
                    (env_struct%env_hydrogen_fraction(env_point_index)-env_struct%env_hydrogen_fraction(env_point_index-1))
               composition(i_metals,zone_index) = env_struct%env_metal_fraction(env_point_index-1)+interp_fraction* &
                    (env_struct%env_metal_fraction(env_point_index)-env_struct%env_metal_fraction(env_point_index-1))
               do k = i_he3,species_end_index
                  composition(k,zone_index) = composition(k,num_zones)
               end do
               composition(i_he4,zone_index)=1.0D0-composition(i_h1,zone_index)- &
                    composition(i_metals,zone_index)-composition(i_he3,zone_index)
               if(env_struct%env_convective_flag(env_point_index).or.env_struct%env_convective_flag(env_point_index-1))then
                  convective_flag(zone_index) = .true.
               else
                  convective_flag(zone_index) = .false.
               endif
            endif
            num_zones = zone_index
            exit
         endif
      end do
! ASSIGN THE BOUNDARY AT THE PHOTOSPHERE FOR ENVELOPE MASS BELOW 1.0D-12.
! (On the exit path above num_zones was just set to zone_index, so this
! guard is false there; on fall-through num_zones is unchanged.)
      if (zone_index .gt. num_zones + env_struct%num_env_points) then
      num_zones = num_zones + env_struct%num_env_points
      end if
end subroutine append_envelope_points

end module envelope_refit_lib
