!----------------------------------------------------------------------
! intmom
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original intmom.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Returns the moment of inertia per unit mass (moment_of_inertia_per_mass,
! originally HIM) and dI/d(omega) per unit mass (di_domega_per_mass,
! originally QIWM) at a single zone center, given the local rotation
! parameter (rotation_param, originally A -- essentially
! omega**2*r0**3/(GM*(2+eta**2)), a measure of the centrifugal-to-
! gravitational force ratio used to expand the moment of inertia in
! powers of the rotational distortion) and the oblateness parameter
! eta_squared (originally ETA2X). Called once per zone by momi.f90.
subroutine intmom(rotation_param, eta_squared, dlnr0_dlnr, r0_geom_factor, &
     moment_of_inertia_per_mass, di_domega_per_mass)
      use const_lib
      implicit none

      double precision, intent(in) :: rotation_param, eta_squared, &
           dlnr0_dlnr, r0_geom_factor
      double precision, intent(out) :: moment_of_inertia_per_mass, &
           di_domega_per_mass


! series_coeff: coefficients of the power series (in rotation_param) used
! to correct the spherical moment of inertia for rotational distortion.
      double precision :: series_coeff(5)
      data series_coeff/0.6D0,0.72857142857143D0,0.50884115884115D0, &
           0.49728653699242D0,0.42852079651326D0/
      integer :: term_idx
      double precision :: leading_factor, series_sum, series_sum_domega

! RETURN MOMENT OF INERTIA PER UNIT MASS AT THE SPECIFIED POINT(HI)
! AND DI/D OMEGA PER UNIT MASS(QIW)
      leading_factor = cc23*dlnr0_dlnr*r0_geom_factor
      series_sum = 0.0d0
      series_sum_domega = 0.0d0
      do term_idx = 1,5
         series_sum = series_sum + series_coeff(term_idx)* &
              (term_idx*eta_squared + 5.0d0)*rotation_param**term_idx
         series_sum_domega = series_sum_domega + &
              term_idx*series_coeff(term_idx)*(term_idx*eta_squared + &
              5.0d0)*rotation_param**term_idx
      end do
      moment_of_inertia_per_mass = leading_factor*(1.0d0 + 1.5d-1*series_sum)
      di_domega_per_mass = leading_factor*(1.5d-1*series_sum_domega)

      return
end subroutine intmom
