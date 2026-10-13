!----------------------------------------------------------------------
! envelope_refit
!----------------------------------------------------------------------
! Appending envelope points to the model as shells, the step the two
! envelope refits share. Before the call, atm_get has integrated a
! fresh envelope from the surface down to the current fitting point,
! and the caller has shifted env_points so that its point 1 is the
! first point outside the last interior shell.
!
! Callers: core/envelope/envelope_rebuild.f90 (the fitting point moved
! outward during the run) and core/model/starting_model.f90 (a
! requested envelope thinner than the starting model's). What they do
! before the call (the surface L, Teff and g, the env_points shift)
! and after it (shell masses, rotation, prints) stays with them. Inside
! the loop they differ in one test, carried by
! accept_point_at_fit_mass: whether an envelope point exactly at the
! fitting mass is taken as a shell (read_starting_model, .le.) or
! interpolated (rebuild_envelope, .lt.). Each prints its own message
! when ierr comes back nonzero.
module envelope_refit
      use star_info_lib, only: star_info
      implicit none
      private
      public :: append_envelope_points

contains

! ---------------------------------------------------------------
! Appends the env_points points inside the new fitting mass star%senv
! to the model as zones num_zones+1, num_zones+2, ..., and the first
! point outside it as a last zone interpolated linearly in log m onto
! the fitting mass; advances num_zones. Each new zone keeps the
! luminosity of the last interior shell (log_luminosity holds L in
! Lsun, linear, despite its name); it takes X and Z from the envelope,
! He3 and the heavier tracers from the last interior shell, and
! Y = 1 - X - Z - X(He3). ierr = 1 when the mass interval of the
! interpolation is below 1e-14 dex.
subroutine append_envelope_points(s, accept_point_at_fit_mass, species_end_index, &
     composition, log_density, log_luminosity, log_pressure, log_radius, &
     log_mass, log_temperature, convective_flag, num_zones, ierr)
      use star_info_lib, only: json, n_species_extended, i_h1, i_he4, &
           i_metals, i_he3
      type(star_info), intent(inout) :: s
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
      do zone_index = num_zones+1,num_zones+s%env_points%num_env_points
         env_point_index = zone_index-num_zones
! The luminosity is constant through the envelope.
         log_luminosity(zone_index) = log_luminosity(num_zones)
! Points inside the new fitting mass become zones as they are.
         if (accept_point_at_fit_mass) then
            point_below_fit_mass = s%env_points%env_log10_mass(env_point_index).le.s%senv
         else
            point_below_fit_mass = s%env_points%env_log10_mass(env_point_index).lt.s%senv
         end if
         if(point_below_fit_mass)then
            log_density(zone_index) = s%env_points%env_log10_density(env_point_index)
            log_pressure(zone_index) = s%env_points%env_log10_pressure(env_point_index)
            log_radius(zone_index) = s%env_points%env_log10_radius(env_point_index)
            log_mass(zone_index) = s%env_points%env_log10_mass(env_point_index) + s%stotal
            log_temperature(zone_index) = s%env_points%env_log10_temperature(env_point_index)
            composition(i_h1,zone_index) = s%env_points%env_hydrogen_fraction(env_point_index)
            composition(i_metals,zone_index) = s%env_points%env_metal_fraction(env_point_index)
            do k = i_he3,species_end_index
               composition(k,zone_index) = composition(k,num_zones)
            end do
            composition(i_he4,zone_index)=1.0D0-composition(i_h1,zone_index)- &
                 composition(i_metals,zone_index)-composition(i_he3,zone_index)
            convective_flag(zone_index) = s%env_points%env_convective_flag(env_point_index)
         else
! The first point outside the fitting mass: the last zone is
! interpolated linearly in log m onto the fitting mass,
!   v = v0 + f (v2 - v0),  f = (log m_fit - log m0)/(log m2 - log m0),
! between the zone or point below it (0) and this point (2); then
! num_zones is set and the loop ends.
            if(env_point_index.eq.1)then
! Between the last interior shell and the first envelope point; the
! zone is convective if either is. X and Z use the difference
! X(nz) - X_env, the reverse of the other variables; it vanishes when
! the envelope has the shell's X and Z, as both callers integrate it.
               mass_interp_x0 = log_mass(num_zones)
               mass_interp_x1 = s%stotal + s%senv
               mass_interp_x2 = s%env_points%env_log10_mass(env_point_index) + s%stotal
               if(mass_interp_x2-mass_interp_x0.lt.1.0D-14) then
                  ierr = 1
                  return
               end if
               interp_fraction = (mass_interp_x1-mass_interp_x0)/ &
                    (mass_interp_x2-mass_interp_x0)
               log_density(zone_index) = log_density(num_zones)+interp_fraction* &
                    (s%env_points%env_log10_density(env_point_index)-log_density(num_zones))
               log_pressure(zone_index) = log_pressure(num_zones)+interp_fraction* &
                    (s%env_points%env_log10_pressure(env_point_index)-log_pressure(num_zones))
               log_radius(zone_index) = log_radius(num_zones)+interp_fraction* &
                    (s%env_points%env_log10_radius(env_point_index)-log_radius(num_zones))
               log_mass(zone_index) = mass_interp_x1
               log_temperature(zone_index) = log_temperature(num_zones)+interp_fraction* &
                    (s%env_points%env_log10_temperature(env_point_index)-log_temperature(num_zones))
               composition(i_h1,zone_index) = composition(i_h1,num_zones)+interp_fraction* &
                    (composition(i_h1,num_zones)-s%env_points%env_hydrogen_fraction(env_point_index))
               composition(i_metals,zone_index) = composition(i_metals,num_zones)+interp_fraction* &
                    (composition(i_metals,num_zones)-s%env_points%env_metal_fraction(env_point_index))
               do k = i_he3,species_end_index
                  composition(k,zone_index) = composition(k,num_zones)
               end do
               composition(i_he4,zone_index)=1.0D0-composition(i_h1,zone_index)- &
                    composition(i_metals,zone_index)-composition(i_he3,zone_index)
               if(s%env_points%env_convective_flag(env_point_index).or.convective_flag(num_zones))then
                  convective_flag(zone_index) = .true.
               else
                  convective_flag(zone_index) = .false.
               endif
            else
! Between the last two envelope points; the zone is convective if
! either is.
               mass_interp_x0 = s%env_points%env_log10_mass(env_point_index-1) + s%stotal
               mass_interp_x1 = s%stotal + s%senv
               mass_interp_x2 = s%env_points%env_log10_mass(env_point_index) + s%stotal
               if(mass_interp_x2-mass_interp_x0.lt.1.0D-14) then
                  ierr = 1
                  return
               end if
               interp_fraction = (mass_interp_x1-mass_interp_x0)/ &
                    (mass_interp_x2-mass_interp_x0)
               log_density(zone_index) = s%env_points%env_log10_density(env_point_index-1)+interp_fraction* &
                    (s%env_points%env_log10_density(env_point_index)-s%env_points%env_log10_density(env_point_index-1))
               log_pressure(zone_index) = s%env_points%env_log10_pressure(env_point_index-1)+interp_fraction* &
                    (s%env_points%env_log10_pressure(env_point_index)-s%env_points%env_log10_pressure(env_point_index-1))
               log_radius(zone_index) = s%env_points%env_log10_radius(env_point_index-1)+interp_fraction* &
                    (s%env_points%env_log10_radius(env_point_index)-s%env_points%env_log10_radius(env_point_index-1))
               log_mass(zone_index) = mass_interp_x1
               log_temperature(zone_index) = s%env_points%env_log10_temperature(env_point_index-1)+interp_fraction* &
                    (s%env_points%env_log10_temperature(env_point_index)-s%env_points%env_log10_temperature(env_point_index-1))
               composition(i_h1,zone_index) = s%env_points%env_hydrogen_fraction(env_point_index-1)+interp_fraction* &
                    (s%env_points%env_hydrogen_fraction(env_point_index)-s%env_points%env_hydrogen_fraction(env_point_index-1))
               composition(i_metals,zone_index) = s%env_points%env_metal_fraction(env_point_index-1)+interp_fraction* &
                    (s%env_points%env_metal_fraction(env_point_index)-s%env_points%env_metal_fraction(env_point_index-1))
               do k = i_he3,species_end_index
                  composition(k,zone_index) = composition(k,num_zones)
               end do
               composition(i_he4,zone_index)=1.0D0-composition(i_h1,zone_index)- &
                    composition(i_metals,zone_index)-composition(i_he3,zone_index)
               if(s%env_points%env_convective_flag(env_point_index).or.s%env_points%env_convective_flag(env_point_index-1))then
                  convective_flag(zone_index) = .true.
               else
                  convective_flag(zone_index) = .false.
               endif
            endif
            num_zones = zone_index
            exit
         endif
      end do
! Every envelope point, the photosphere included, lay inside the
! fitting mass (an envelope below about 1e-12 of the star, senv at 0):
! all became zones, and the model now ends at the photosphere. (After
! the exit above, num_zones is already zone_index, and this test is
! false.)
      if (zone_index .gt. num_zones + s%env_points%num_env_points) then
      num_zones = num_zones + s%env_points%num_env_points
      end if
end subroutine append_envelope_points

end module envelope_refit
