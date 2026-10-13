!----------------------------------------------------------------------
! surface_boundary
!----------------------------------------------------------------------
! The outer boundary condition of the Henyey solve. Envelope
! integrations at the three vertices of a triangle in the
! (log Teff, log L) plane around the current model each give log P,
! log T and log R at the fitting point; compute_surface_boundary fits
! log R, log L and log Teff there as linear functions of the
! fitting-point log P and log T, which henyey_iterate turns into the
! surface rows of the Henyey system. The triangle's geometry and the fit
! are the star-blind kernels of triangle_fit.f90.
module surface_boundary
      implicit none
      private
      public :: compute_surface_boundary

contains

!----------------------------------------------------------------------
! compute_surface_boundary (legacy surfbc.f)
!----------------------------------------------------------------------
! Maintains the triangle of envelope integrations around the current
! (log Teff, log L) = (log10_teff, log_luminosity) and its fit,
! envelope_coeffs = c:
!   log R    = c(1) log P + c(2) log T + c(3)
!   log L    = c(4) log P + c(5) log T + c(6)
!   log Teff = c(7) log P + c(8) log T + c(9)
! at the fitting point (a plane through the three vertices, for each).
! In steps: switch to or from the gray atmosphere when Teff crosses the
! edge of a tabulated (Kurucz or Allard) T-tau relation; start a new
! triangle, or reflect vertices until the triangle contains the point;
! integrate a new envelope (envelope_integration's atm_get) at every
! vertex that needs one; refit the coefficients. The geometry
! (start_triangle, reflect_until_inside) and the fit
! (fit_triangle_coefficients) are triangle_fit's kernels; this routine
! gathers their inputs, integrates the envelopes, and stores the results.
! tri_teffl and tri_logl hold the vertices, and vtx_logp, vtx_logt and
! vtx_logr each vertex's fitting-point log P, log T and log R.
! start_new_triangle builds three new envelopes about the point;
! reset_triangle redoes all three and retriangulates if needed. Both
! come back false. envelope_recomputed_flag is true if any envelope was
! integrated. (log10_teff, the model's, used to come back as the last
! integrated vertex's log Teff, which the first Newton pass then read;
! each vertex now has its own, vertex_log10_teff.)
! hydrogen_fraction, metal_fraction and the two rotation factors are
! intent(inout) only because atm_get's dummies are (they pass through to
! the Bulirsch-Stoer integrator's callbacks); the caller passes locals.
!   star reads: atm_choice_initial, ctrl%spot_filling_factor,
!     ctrl%spot_temp_contrast, ctrl%tri_delta_logl, ctrl%tri_delta_teffl
!   star writes: job%atm_choice, use_ttau_relation
subroutine compute_surface_boundary(tri_teffl, tri_logl, envelope_coeffs, &
     vtx_logp, vtx_logt, &
     vtx_logr, tri_orientation, stored_vertex_index, &
     stored_envelope_state, &
     start_new_triangle, reset_triangle, saha_state, &
     log10_star_mass, log_luminosity, &
     log10_teff, hydrogen_fraction, metal_fraction, &
     pressure_rotation_factor, temperature_rotation_factor, &
     envelope_recomputed_flag, log10_pressure_limit, convective_flag, &
     zone_index, ierr)

      use envelope_integration, only: atm_get, surface_log10_radius_cm, ambient_log10_teff
      use triangle_fit, only: start_triangle, reflect_until_inside, &
           fit_triangle_coefficients
      use atm_lib, only: atm_allard_teff_max
      use star_info_lib
      use io_units
      use run_log, only: solver_diagnostics
      use physical_constants
      use math_lib
      implicit none
      double precision, intent(inout) :: tri_teffl(3), tri_logl(3), &
           envelope_coeffs(9)
      double precision, intent(inout) :: vtx_logp(3), &
           vtx_logt(3), vtx_logr(3)
      double precision, intent(inout) :: tri_orientation
      integer, intent(inout) :: stored_vertex_index
      double precision, intent(inout) :: stored_envelope_state(4)
      logical, intent(inout) :: start_new_triangle, reset_triangle
      integer, intent(inout) :: saha_state
      double precision, intent(in) :: log10_star_mass, log_luminosity
      double precision, intent(in) :: log10_teff
      double precision, intent(inout) :: hydrogen_fraction, metal_fraction, &
           pressure_rotation_factor, temperature_rotation_factor
      logical, intent(out) :: envelope_recomputed_flag
      double precision, intent(in) :: log10_pressure_limit
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: zone_index

      logical :: tri_vertex_valid(3)

      integer, intent(out) :: ierr
      integer :: jerr_atm
! --- locals ---
      integer :: i, vertex_being_computed, j
      double precision :: tri_err
      logical :: envelope_needs_recompute, print_envelope_flag, &
           save_boundary_flag
      double precision :: b, gl, rl, adjusted_teffl, vertex_log10_teff

      ierr = 0
      call switch_atmosphere_choice()
      call mark_triangle_vertices()
      call locate_point_in_triangle()
      call integrate_vertex_envelopes()
      if (ierr /= 0) return
      call refit_coefficients()
      return

contains

! ---------------------------------------------------------------
! Fall back to the gray atmosphere (atm_choice 0) and start a new
! triangle when log Teff reaches 3.95 with the Kurucz table (atm_choice
! 3) or the Allard table's maximum (atm_choice 4) (MHP 9/01); switch
! back to the initial choice once Teff is below again, with a new
! triangle as well. (The switch back used to keep the triangle, whose
! vertices held gray-atmosphere envelopes, until each was reflected.)
subroutine switch_atmosphere_choice()
      if (star%job%atm_choice.eq.3) then
         if (log10_teff.ge.3.95d0) then
            write(*,*)
            write(*,5) log10_teff
 5          format('LOG TEFF OF ',F7.3,' ABOVE 3.95 - SWITCH' &
                   ,'TO GRAY ATMOSPHERE BOUNDARY CONDITION')
            write(*,*)
            star%job%atm_choice = 0
            start_new_triangle = .true.
! Remember the switch.
            star%use_ttau_relation = .true.
         endif
      endif
      if (star%job%atm_choice.eq.4) then
         if (log10_teff.ge.atm_allard_teff_max()) then
            write(*,*)
            write(*,7) log10_teff, atm_allard_teff_max()
 7          format('LOG TEFF OF ',F7.3,' ABOVE Allard Table max ',F7.3 &
                   ,'  - SWITCH TO GRAY ATMOSPHERE BOUNDARY CONDITION')
            write(*,*)
            star%job%atm_choice = 0
            start_new_triangle = .true.
! Remember the switch.
            star%use_ttau_relation = .true.
         endif
      endif
      if (star%use_ttau_relation) then
         if (star%atm_choice_initial.eq.3.and.log10_teff.lt.3.95d0) then
            star%job%atm_choice = star%atm_choice_initial
            start_new_triangle = .true.
            star%use_ttau_relation = .false.
            write(*,9) log10_teff
 9          format('LOG TEFF OF ',F7.3,' BELOW 3.95 - SWITCH' &
           ,' BACK TO KURUCZ ATMOSPHERE BOUNDARY CONDITION')
         else if (star%atm_choice_initial.eq.4.and.log10_teff.lt.atm_allard_teff_max()) then
            star%job%atm_choice = star%atm_choice_initial
            start_new_triangle = .true.
            star%use_ttau_relation = .false.
            write(*,11) log10_teff, atm_allard_teff_max()
 11         format('LOG TEFF OF ',F7.3,' below Allard Table max ',F7.3 &
           ,'  - SWITCH BACK TO ALLARD ATMOSPHERE BOUNDARY CONDITION')
         endif
      endif
end subroutine switch_atmosphere_choice

! ---------------------------------------------------------------
! No envelope recomputed yet; all three vertices need new envelopes for
! a new or reset triangle, otherwise all three are still valid.
! tri_err, the tolerance of the inside test, is always zero.
subroutine mark_triangle_vertices()
      envelope_recomputed_flag = .false.
      if (start_new_triangle.or.reset_triangle) then
       tri_vertex_valid(1) = .false.
       tri_vertex_valid(2) = .false.
       tri_vertex_valid(3) = .false.
       stored_vertex_index = 0
       tri_err = 0.0d0
       reset_triangle = .false.
      else
       tri_vertex_valid(1) = .true.
       tri_vertex_valid(2) = .true.
       tri_vertex_valid(3) = .true.
       tri_err = 0.0d0
      endif
end subroutine mark_triangle_vertices

! ---------------------------------------------------------------
! Start a triangle about (Teff, L), dT = tri_delta_teffl wide and
! dL = tri_delta_logl high, or reflect vertices of the old one through
! the midpoints of their opposite edges until it contains the point; a
! reflected vertex needs a new envelope.
subroutine locate_point_in_triangle()
      logical :: reflected(3)
      integer :: n_reflections
      if (start_new_triangle) then
       call start_triangle(log10_teff, log_luminosity, &
            star%ctrl%tri_delta_teffl, star%ctrl%tri_delta_logl, &
            tri_teffl, tri_logl, tri_orientation)
       start_new_triangle = .false.
      else
       call reflect_until_inside(tri_teffl, tri_logl, tri_orientation, &
            log10_teff, log_luminosity, tri_err, reflected, n_reflections)
       do i = 1,3
          if (reflected(i)) tri_vertex_valid(i) = .false.
       end do
      endif
end subroutine locate_point_in_triangle

! ---------------------------------------------------------------
! Integrate a new envelope at every vertex that needs one; a failure
! leaves ierr nonzero.
subroutine integrate_vertex_envelopes()
      envelope_needs_recompute = .false.
      do i = 1,3
       if (.not.tri_vertex_valid(i)) then
! The vertex's L (Lsun), its radius from L = 4 pi R^2 sigma Teff^4, and
! its surface gravity, log g = log G + log M - 2 log R.
          envelope_recomputed_flag = .true.
          vertex_log10_teff = tri_teffl(i)
          b = exp(ln10*tri_logl(i))
          rl = surface_log10_radius_cm(tri_logl(i), vertex_log10_teff)
          gl = cgl + log10_star_mass - rl - rl
          print_envelope_flag = .false.
          vertex_being_computed=i
          save_boundary_flag = .true.
! Spotted runs (G. Somers): when the outermost shell is convective and
! spots are on, the atmosphere starts from the Teff of the unspotted
! (ambient) surface.
          if (convective_flag(zone_index).and.star%ctrl%spot_filling_factor.ne.0.0.and.star%ctrl%spot_temp_contrast.ne.1.0) then
               adjusted_teffl = ambient_log10_teff(vertex_log10_teff)
          else
             adjusted_teffl = vertex_log10_teff
          endif
          call atm_get(b,pressure_rotation_factor,temperature_rotation_factor,gl, &
                 log10_star_mass,print_envelope_flag, &
                 save_boundary_flag,log10_pressure_limit,rl, &
                 adjusted_teffl,hydrogen_fraction,metal_fraction, &
                 saha_state,ierr=jerr_atm, &
                 vertex_index=vertex_being_computed, &
                 stored_envelope_state=stored_envelope_state, &
                 stored_vertex_index=stored_vertex_index, &
                 vtx_logp=vtx_logp,vtx_logr=vtx_logr,vtx_logt=vtx_logt)
! An integration failure (numerics_termination) or a table error comes
! back in jerr_atm; the caller (refresh_surface_bc) propagates it.
          if (jerr_atm /= 0) then
             ierr = jerr_atm
             return
          end if
          envelope_needs_recompute = .true.
       endif
      end do
end subroutine integrate_vertex_envelopes

! ---------------------------------------------------------------
! Refit the coefficients when any envelope is new: for log R, log L and
! log Teff, the plane in the fitting-point (log P, log T) through the
! three vertices.
subroutine refit_coefficients()
      if (envelope_needs_recompute) then
       call fit_triangle_coefficients(vtx_logp, vtx_logt, vtx_logr, &
            tri_logl, tri_teffl, envelope_coeffs)
       if (solver_diagnostics()) then
        write(run_log_unit,50)(i,tri_vertex_valid(i),tri_teffl(i),tri_logl(i),vtx_logp(i),vtx_logt(i),vtx_logr(i), &
              (envelope_coeffs(3*(i-1)+j),j=1,3), i=1,3)
       end if
 50      format(' ENVELOPE TRIANGLE  LOG(TE)  LOG(L)   LOG(P)   LOG(T)', &
              3X,'LOG(R)   COEFFICIENTS'/I14,L1,2X,5F9.5,8X,'LOG(R)', &
              3F10.5/ I14,L1,2X,5F9.5,5X,'LOG(L/L0)',3F10.5/ &
              I14,L1,2X,5F9.5,7X,'LOG(TE)',3F10.5)
      else
       if (solver_diagnostics()) write(run_log_unit,70)
 70      format(' ENVELOPE TRIANGLE NOT CHANGED')
      endif
end subroutine refit_coefficients

end subroutine compute_surface_boundary

end module surface_boundary
