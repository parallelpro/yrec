!----------------------------------------------------------------------
! microdiff_driver
!----------------------------------------------------------------------
! The newer settling pipeline (use_new_diffusion_routines), the one
! every shipped deck runs. microdiff settles hydrogen and helium, the
! metals (use_diffusion_z) and the light elements Li6, Li7 and Be9
! (diffuse_lithium) over one step, in the formalism of Bahcall & Loeb
! (1990) with the coefficients of Thoul, Bahcall & Loeb (1994). The
! abundance of each species obeys
!     dX/dt = (1/(rho r^2)) d/dr [ D1 + D2 dX/dr ],
! a drift term D1 (pressure, temperature and hydrogen-concentration
! gradients) and a diffusion term D2. Its stages, over the settling
! region diffuse_elements found:
!   1. convert_to_bahcall_loeb: the model's quantities in Bahcall & Loeb
!      units (bahcall_loeb_units), on copies;
!   2. map_microdiff_to_equal_grid (microdiff_mapping): two equally
!      spaced grids in radius, at the centres and the midpoints, and the
!      structure and the abundances on them;
!   3. run_microdiff (microdiff_run), once for hydrogen (helium follows
!      from it), once for the metals and once per light element: the
!      species' change over the step, the drift term by the two-step
!      Lax-Wendroff scheme and the diffusion term by one implicit solve,
!      with the coefficients of microdiff_coefficients; its messages
!      (report_species_solve);
!   4. microdiff_etm (microdiff_mapping): the changes back on the
!      model's zones, and Y from X + Y + Z + He3 = 1.
! The species is solved for as column 3 of the Thoul et al. vectors, the
! heavy species, except hydrogen (column 1), which takes iron as the
! heavy species: the metals are fully ionized iron (species_table).
module microdiff_driver
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: json
      use element_diffusion_configuration, only: diffusion_settings
      implicit none
      private
      public :: microdiff

contains

! ---------------------------------------------------------------
! One step of timestep seconds over the settling region
! zone_begin..zone_end of a model of num_zones zones (centre first):
! the mass fractions composition(species, zone) change in place. The
! structure: log10 radius (cm), density (g/cm^3) and temperature (K),
! the enclosed mass (g), d ln P/dr (1/cm), the temperature gradient
! del_grad = d ln T/d ln P, and the total mass (g); none of them
! changes. ierr: a numerics routine's status (singular matrix), after
! its message.
subroutine microdiff(timestep, composition, dlnp_dr, log_radius, &
     log_density, enclosed_mass, log_temperature, del_grad, &
     num_zones, total_mass, zone_begin, zone_end, ierr)
      use microdiff_mapping, only: microdiff_grid, map_microdiff_to_equal_grid, &
           microdiff_etm
      use microdiff_run, only: species_solve_report, run_microdiff
      use species_table, only: thoul_fe, thoul_col_h, thoul_col_he, &
           thoul_col_metal, num_light_diffused, light_diffused, light_diffused_row
      real(dp), intent(in) :: timestep
      real(dp), intent(inout) :: composition(:,:)
      real(dp), intent(in) :: dlnp_dr(:), log_radius(:), log_density(:), &
           enclosed_mass(:), log_temperature(:), del_grad(:)
      integer, intent(in) :: num_zones
      real(dp), intent(in) :: total_mass
      integer, intent(in) :: zone_begin, zone_end
      integer, intent(out) :: ierr
! The timestep, d ln P/dr, the enclosed mass and the total mass in
! Bahcall & Loeb units: copies, so that the model's own stay exactly as
! they were (docs/bug_sweep.md O16).
      real(dp) :: timestep_bl, dlnp_dr_bl(json), enclosed_mass_bl(json), &
           total_mass_bl
! the radius in Bahcall & Loeb units; the density (g/cm^3) and the
! temperature (K)
      real(dp) :: radius_bl(json), density(json), temperature(json)
! the equally spaced grids, at the centres and the midpoints
      type(microdiff_grid) :: centre_grid, mid_grid
      integer :: num_eq_points
      real(dp) :: grid_spacing
! the mass fractions of the Thoul columns (hydrogen, helium, the species
! solved for as the heavy one) on the two grids
      real(dp) :: species_fraction(3,json), species_fraction_mid(3,json)
! d ln C_H/dr on the two grids: the hydrogen solve computes it, the
! others read it
      real(dp) :: hydrogen_dlnc_dr(json), hydrogen_dlnc_dr_mid(json)
! each species' change on the centre grid
      real(dp) :: hydrogen_change(json), metal_change(json), &
           light_change(num_light_diffused,json)
! what a species' solve did, for its messages
      type(species_solve_report) :: report
      integer :: k

      ierr = 0
! (no change of the species a setting leaves out)
      metal_change = 0.0d0
      light_change = 0.0d0
      timestep_bl = timestep
      dlnp_dr_bl(1:num_zones) = dlnp_dr(1:num_zones)
      enclosed_mass_bl(1:num_zones) = enclosed_mass(1:num_zones)
      total_mass_bl = total_mass
      call convert_to_bahcall_loeb(timestep_bl, dlnp_dr_bl, log_radius, log_density, &
           enclosed_mass_bl, log_temperature, num_zones, total_mass_bl, radius_bl, &
           density, temperature)
      call map_microdiff_to_equal_grid(num_light_diffused, light_diffused_row, composition, &
           dlnp_dr_bl, radius_bl, enclosed_mass_bl, del_grad, zone_begin, zone_end, &
           num_zones, grid_spacing, num_eq_points, density, temperature, &
           centre_grid, mid_grid)
      call load_column(thoul_col_h, centre_grid%hydrogen, mid_grid%hydrogen, &
           num_eq_points, species_fraction, species_fraction_mid)
      call load_column(thoul_col_he, centre_grid%helium, mid_grid%helium, &
           num_eq_points, species_fraction, species_fraction_mid)

! Hydrogen, iron as the heavy species; then the hydrogen column back to
! the grids' X (the solve leaves the change in it).
      call load_column(thoul_col_metal, centre_grid%metal, mid_grid%metal, &
           num_eq_points, species_fraction, species_fraction_mid)
      call run_microdiff(grid_spacing, timestep_bl, total_mass_bl, &
           num_eq_points, centre_grid, species_fraction, hydrogen_dlnc_dr, mid_grid, &
           species_fraction_mid, hydrogen_dlnc_dr_mid, thoul_fe%weight, &
           thoul_fe%charge, thoul_col_h, report, ierr)
      call report_species_solve(report, ierr)
      if (ierr /= 0) return
      hydrogen_change(1:num_eq_points) = species_fraction(thoul_col_h,1:num_eq_points)
      call load_column(thoul_col_h, centre_grid%hydrogen, mid_grid%hydrogen, &
           num_eq_points, species_fraction, species_fraction_mid)

! The metals, as iron.
      if (diffusion_settings%use_diffusion_z) then
         call load_column(thoul_col_metal, centre_grid%metal, mid_grid%metal, &
              num_eq_points, species_fraction, species_fraction_mid)
         call run_microdiff(grid_spacing, timestep_bl, total_mass_bl, &
              num_eq_points, centre_grid, species_fraction, hydrogen_dlnc_dr, mid_grid, &
              species_fraction_mid, hydrogen_dlnc_dr_mid, thoul_fe%weight, &
              thoul_fe%charge, thoul_col_metal, report, ierr)
         call report_species_solve(report, ierr)
         if (ierr /= 0) return
         metal_change(1:num_eq_points) = species_fraction(thoul_col_metal,1:num_eq_points)
      end if

! The light elements, one at a time.
      if (diffusion_settings%diffuse_lithium) then
         do k = 1, num_light_diffused
            call load_column(thoul_col_metal, centre_grid%light(k,:), mid_grid%light(k,:), &
                 num_eq_points, species_fraction, species_fraction_mid)
            call run_microdiff(grid_spacing, timestep_bl, total_mass_bl, &
                 num_eq_points, centre_grid, species_fraction, hydrogen_dlnc_dr, mid_grid, &
                 species_fraction_mid, hydrogen_dlnc_dr_mid, light_diffused(k)%weight, &
                 light_diffused(k)%charge, thoul_col_metal, report, ierr)
            call report_species_solve(report, ierr)
            if (ierr /= 0) return
            light_change(k,1:num_eq_points) = species_fraction(thoul_col_metal,1:num_eq_points)
         end do
      end if

      call microdiff_etm(centre_grid%radius, hydrogen_change, metal_change, light_change, &
           zone_begin, zone_end, num_eq_points, composition, radius_bl, num_zones, &
           num_light_diffused, light_diffused_row)
end subroutine microdiff

! ---------------------------------------------------------------
! Column col of the Thoul vectors on the two grids: centre_values at the
! num_eq_points centres, mid_values at the num_eq_points - 1 midpoints.
pure subroutine load_column(col, centre_values, mid_values, num_eq_points, &
     species_fraction, species_fraction_mid)
      integer, intent(in) :: col, num_eq_points
      real(dp), intent(in) :: centre_values(:), mid_values(:)
      real(dp), intent(inout) :: species_fraction(:,:), species_fraction_mid(:,:)

      species_fraction(col,1:num_eq_points) = centre_values(1:num_eq_points)
      species_fraction_mid(col,1:num_eq_points-1) = mid_values(1:num_eq_points-1)
end subroutine load_column

! ---------------------------------------------------------------
! The messages of one species' solve (element_diffusion_messages): the
! numerics routine that failed (ierr nonzero); else each pass of the
! implicit solve under solver diagnostics (pass 1's change, then zero),
! and a solve that missed settling_tolerance in its
! settling_num_iterations passes, with the last pass's change.
subroutine report_species_solve(report, ierr)
      use microdiff_run, only: species_solve_report
      use element_diffusion_messages, only: note_solve_pass, report_unconverged, &
           report_numerics_failure
      type(species_solve_report), intent(in) :: report
      integer, intent(in) :: ierr
      integer :: pass

      if (ierr /= 0) then
         call report_numerics_failure(trim(report%failed_routine), ierr)
         return
      end if
      do pass = 1, report%num_passes
         if (pass == 1) then
            call note_solve_pass(pass, report%change, report%change_shell)
         else
            call note_solve_pass(pass, 0.0d0, 0)
         end if
      end do
      if (report%converged) return
      if (report%num_passes <= 1) then
         call report_unconverged('MICRODIFF', 'D', diffusion_settings%settling_tolerance, &
              diffusion_settings%settling_num_iterations, report%change, report%change_shell)
      else
         call report_unconverged('MICRODIFF', 'D', diffusion_settings%settling_tolerance, &
              diffusion_settings%settling_num_iterations, 0.0d0, 0)
      end if
end subroutine report_species_solve

! ---------------------------------------------------------------
! microdiff's first stage: the model's quantities in Bahcall & Loeb
! units (bahcall_loeb_units), on microdiff's copies. timestep (s),
! dlnp_dr (1/cm), enclosed_mass and total_mass (g) are converted in
! place; radius_bl is the radius in Bahcall & Loeb units, and density
! (g/cm^3) and temperature (K) are unlogged, for the equal grids.
pure subroutine convert_to_bahcall_loeb(timestep, dlnp_dr, log_radius, log_density, &
     enclosed_mass, log_temperature, num_zones, total_mass, radius_bl, density, &
     temperature)
      use bahcall_loeb_units, only: bl_radius_scale, bl_mass_scale, bl_time_scale
      use physical_constants, only: ln10
      use math_lib
      real(dp), intent(inout) :: timestep, total_mass
      real(dp), intent(inout) :: dlnp_dr(:), enclosed_mass(:)
      real(dp), intent(in) :: log_radius(:), log_density(:), log_temperature(:)
      integer, intent(in) :: num_zones
      real(dp), intent(inout) :: radius_bl(:), density(:), temperature(:)
      integer :: k

      do k = 1, num_zones
         radius_bl(k) = exp(ln10*log_radius(k))*bl_radius_scale
         enclosed_mass(k) = enclosed_mass(k)*bl_mass_scale
         dlnp_dr(k) = dlnp_dr(k)/bl_radius_scale
      end do
      timestep = timestep/bl_time_scale
      total_mass = total_mass*bl_mass_scale
      do k = 1, num_zones
         density(k) = exp(ln10*log_density(k))
         temperature(k) = exp(ln10*log_temperature(k))
      end do
end subroutine convert_to_bahcall_loeb

end module microdiff_driver

! Sources
! - G. Somers (05/2015) wrote microdiff for the light elements, adapting
!   the original routines (original_settling; MHP and BC, 1989-1994);
!   JvS (01/2026) applied the helium and metal factors
!   (settling_y_factor, settling_z_factor) to it.
! - The two-step Lax-Wendroff scheme: Press et al. (1986), Numerical
!   Recipes, Cambridge University Press, p. 633.
! - 2026: the model's quantities converted on copies (docs/bug_sweep.md
!   O16); the implicit step solved once (M253); the region, the grids,
!   the Coulomb logarithms and the messages shared with the original
!   routines (settling_region, diffusion_grid, coulomb_logarithm,
!   element_diffusion_messages).
