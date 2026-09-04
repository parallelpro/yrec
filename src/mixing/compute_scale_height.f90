!----------------------------------------------------------------------
! hsubp
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original hsubp.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Calculates the local pressure scale height at the edge of a
! convective region (PD Dec. 1984; 2/91 MHP revised for YREC/Mark6).
! Calls the equation of state at the edge zone and returns the
! density-based pressure scale height (pscahe). The inputs are log10
! quantities; the result is in cm.
subroutine compute_scale_height(composition, log_density, log_pressure, &
     log_radius, log_mass, log_temperature, edge_zone, pscahe, ierr)
      use star_info_lib, only: json

      use phys_const_lib
      use eos_lib
      use math_lib
      implicit none

      double precision, intent(in) :: composition(15,json), log_density(json), &
           log_pressure(json), log_radius(json), log_mass(json), &
           log_temperature(json)
      integer, intent(in) :: edge_zone
      double precision, intent(out) :: pscahe
      integer, intent(out) :: ierr

      logical :: lderiv, latmo
      integer :: ksaha
      double precision :: hydrogen_fraction, metal_fraction, log10_pressure, &
           log10_temperature, log10_density
! 2026 named-index results: the former 20-variable eos output soup is
! one result array (see eos_lib's index constants); only the inout
! i_log10_density (the density consistent with P,T,X,Z) is consumed here.
      double precision :: eos_res(num_eos_results)
      double precision :: log10_mass, log10_radius

!  CALL EQUATION OF STATE TO GET THE DENSITY AT THE EDGE ZONE.
      ierr = 0
      lderiv = .false.
      latmo = .true.
      ksaha = 0
      hydrogen_fraction = composition(1,edge_zone)
      metal_fraction = composition(3,edge_zone)
      log10_pressure = log_pressure(edge_zone)
      log10_temperature = log_temperature(edge_zone)
      log10_density = log_density(edge_zone)
      eos_res = 0.0d0
      eos_res(i_log10_density) = log10_density
      call eos_get(log10_temperature, log10_pressure, &
           hydrogen_fraction, metal_fraction, eos_res, lderiv, latmo, &
           ksaha, composition_at_zone=composition(:,edge_zone), ierr=ierr)
      if (ierr /= 0) return
!  COMPUTE PRESSURE SCALE HEIGHT  H_P = P R^2/(G M rho).
      log10_mass = log_mass(edge_zone)
      log10_radius = log_radius(edge_zone)
      pscahe = exp(ln10*(log10_pressure + 2.0d0*log10_radius - &
           eos_res(i_log10_density) - cgl - log10_mass))
      return
end subroutine compute_scale_height
