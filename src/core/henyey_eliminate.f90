!----------------------------------------------------------------------
! henyey_eliminate_lib (was reduce)
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original reduce.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! 2026 de-tramp (ROADMAP item 3): the 39-argument signature is gone.
! The eight star% arrays the single caller (core/henyey_coefficients)
! passed through (elim_coeff/elim_rhs/luminosity_lsun/max_residual/
! logP/logR/log_mass/logT) are read from star directly -- this is
! star-layer solver code, relocated from setup/ to core/ accordingly.
! The 30 per-shell equation terms are two henyey_shell_terms records:
! the equation values and their partials at the current shell (cur)
! and, from the previous call, at the shell below (prev). The caller
! keeps prev with a plain "prev = cur" after each shell.
!
! Classic Henyey-method forward elimination step: given the structure
! equation residuals and their partial derivatives (with respect to
! P, T, R, L) at the two mesh points bracketing shell I, eliminates
! the previous shell's (I-1) unknowns using the recursion coefficients
! carried in star%elim_coeff/star%elim_rhs from the prior call, then
! reduces the resulting 4x4 system for shell I down to expressions for
! P and T in terms of R and L (stored back into elim_coeff/elim_rhs at
! index I for use by the next shell). See Kippenhahn & Weigert-style
! Henyey solvers for the general method; the four rows correspond to
! the pressure, temperature, radius, and luminosity structure
! equations in that order, and the interpretation of the individual
! q*/q*_d* terms below (equation value and its partials at the lower/
! upper mesh point) is inferred from their use, not from surviving
! documentation.
module henyey_eliminate_lib
      implicit none

! Structure-equation terms at one mesh shell: for each of the four
! Henyey equations (pressure, temperature, radius, luminosity), the
! equation value q* and its nonzero partial derivatives. Note qp/qt/
! qt_dl double as pulse scratch at the caller (pt_scr%qp/qt/qtl).
      type :: henyey_shell_terms
         double precision :: qp, qp_dr, qp_dp
         double precision :: qt, qt_dr, qt_dl, qt_dp, qt_dt
         double precision :: qr, qr_dr, qr_dp, qr_dt
         double precision :: ql, ql_dp, ql_dt
      end type henyey_shell_terms

contains

subroutine henyey_eliminate(zone_index, prev, cur)
      use star_info_lib, only: star
      use phys_const_lib
      implicit none

      integer, intent(in) :: zone_index
      type(henyey_shell_terms), intent(in) :: prev, cur

      double precision :: q(4,4)
      double precision :: half_delta_log_mass, pivot, pivotb, div
      integer :: j

! DEFINE R.H.S
      half_delta_log_mass = 0.5D0*(star%log_mass(zone_index)-star%log_mass(zone_index-1))
      star%elim_rhs(1,zone_index) = half_delta_log_mass*(prev%qp + cur%qp) - &
           (star%logP(zone_index) - star%logP(zone_index-1))
      star%max_residual(1) = dmax1(star%max_residual(1),dabs(star%elim_rhs(1,zone_index)))
      star%elim_rhs(2,zone_index) = half_delta_log_mass*(prev%qt + cur%qt) - &
           (star%logT(zone_index) - star%logT(zone_index-1))
      star%max_residual(2) = dmax1(star%max_residual(2),dabs(star%elim_rhs(2,zone_index)))
      star%elim_rhs(3,zone_index) = half_delta_log_mass*(prev%qr + cur%qr) - &
           (star%logR(zone_index) - star%logR(zone_index-1))
      star%max_residual(3) = dmax1(star%max_residual(3),dabs(star%elim_rhs(3,zone_index)))
      star%elim_rhs(4,zone_index) = half_delta_log_mass*(prev%ql + cur%ql) - &
           (star%luminosity_lsun(zone_index) - star%luminosity_lsun(zone_index-1))
      star%max_residual(4) = dmax1(star%max_residual(4),dabs(star%elim_rhs(4,zone_index)))
      half_delta_log_mass = half_delta_log_mass*ln10
! ELIMINATE COLUMNS 1 AND 2
      pivot = -half_delta_log_mass*prev%qp_dr
      q(1,1) = (-1.0D0-half_delta_log_mass*prev%qp_dp) - &
           pivot*star%elim_coeff(3,1,zone_index-1)
      q(1,2) =                  - pivot*star%elim_coeff(3,2,zone_index-1)
      star%elim_rhs(1,zone_index) = star%elim_rhs(1,zone_index) - &
           pivot*star%elim_rhs(3,zone_index-1)
      pivot = -half_delta_log_mass*prev%qt_dr
      pivotb= -half_delta_log_mass*prev%qt_dl
      q(2,1) = -half_delta_log_mass*prev%qt_dp &
           -pivot*star%elim_coeff(3,1,zone_index-1)-pivotb*star%elim_coeff(4,1,zone_index-1)
      q(2,2) =(-1.0D0-half_delta_log_mass*prev%qt_dt) &
           -pivot*star%elim_coeff(3,2,zone_index-1)-pivotb*star%elim_coeff(4,2,zone_index-1)
      star%elim_rhs(2,zone_index) = star%elim_rhs(2,zone_index) - &
           pivot*star%elim_rhs(3,zone_index-1) - pivotb*star%elim_rhs(4,zone_index-1)
      pivot = -1.0D0 - half_delta_log_mass*prev%qr_dr
      q(3,1) = -half_delta_log_mass*prev%qr_dp - pivot*star%elim_coeff(3,1,zone_index-1)
      q(3,2) = -half_delta_log_mass*prev%qr_dt - pivot*star%elim_coeff(3,2,zone_index-1)
      star%elim_rhs(3,zone_index) = star%elim_rhs(3,zone_index) - &
           pivot*star%elim_rhs(3,zone_index-1)
      q(4,1) = -half_delta_log_mass*prev%ql_dp + star%elim_coeff(4,1,zone_index-1)
      q(4,2) = -half_delta_log_mass*prev%ql_dt + star%elim_coeff(4,2,zone_index-1)
      star%elim_rhs(4,zone_index) = star%elim_rhs(4,zone_index) + star%elim_rhs(4,zone_index-1)
! REDUCE 4*4 MATRIX
! PIVOT ON ROW-4 AND COLUMN-6
      star%elim_coeff(4,1,zone_index) = -half_delta_log_mass*cur%ql_dp
      star%elim_coeff(4,2,zone_index) = -half_delta_log_mass*cur%ql_dt
      pivot = -half_delta_log_mass*cur%qt_dl
      q(2,1) = q(2,1) - pivot*q(4,1)
      q(2,2) = q(2,2) - pivot*q(4,2)
      star%elim_rhs(2,zone_index) = star%elim_rhs(2,zone_index) - pivot*star%elim_rhs(4,zone_index)
      star%elim_coeff(2,1,zone_index) = -half_delta_log_mass*cur%qt_dp      - &
           pivot*star%elim_coeff(4,1,zone_index)
      star%elim_coeff(2,2,zone_index) = (1.0D0-half_delta_log_mass*cur%qt_dt) - &
           pivot*star%elim_coeff(4,2,zone_index)
! PIVOT ON ROW-3 AND COLUMN-5
      div = 1.0D0/(1.0D0 - half_delta_log_mass*cur%qr_dr)
      q(3,1) = q(3,1)*div
      q(3,2) = q(3,2)*div
      star%elim_rhs(3,zone_index) = star%elim_rhs(3,zone_index)*div
      div = -half_delta_log_mass*div
      star%elim_coeff(3,1,zone_index) = div*cur%qr_dp
      star%elim_coeff(3,2,zone_index) = div*cur%qr_dt
      pivot = -half_delta_log_mass*cur%qp_dr
      q(1,1) = q(1,1) - pivot*q(3,1)
      q(1,2) = q(1,2) - pivot*q(3,2)
      star%elim_rhs(1,zone_index) = star%elim_rhs(1,zone_index) - pivot*star%elim_rhs(3,zone_index)
      star%elim_coeff(1,1,zone_index) = (1.0D0-half_delta_log_mass*cur%qp_dp) - &
           pivot*star%elim_coeff(3,1,zone_index)
      star%elim_coeff(1,2,zone_index) =                - pivot*star%elim_coeff(3,2,zone_index)
      pivot = -half_delta_log_mass*cur%qt_dr
      q(2,1) = q(2,1) - pivot*q(3,1)
      q(2,2) = q(2,2) - pivot*q(3,2)
      star%elim_rhs(2,zone_index) = star%elim_rhs(2,zone_index) - pivot*star%elim_rhs(3,zone_index)
      star%elim_coeff(2,1,zone_index) = star%elim_coeff(2,1,zone_index) - &
           pivot*star%elim_coeff(3,1,zone_index)
      star%elim_coeff(2,2,zone_index) = star%elim_coeff(2,2,zone_index) - &
           pivot*star%elim_coeff(3,2,zone_index)
! PIVOT ON ROW-2 AND COLUMN-4
      div = 1.0D0/q(2,2)
      q(2,1) = q(2,1)*div
      star%elim_rhs(2,zone_index) = star%elim_rhs(2,zone_index)*div
      star%elim_coeff(2,1,zone_index) = star%elim_coeff(2,1,zone_index)*div
      star%elim_coeff(2,2,zone_index) = star%elim_coeff(2,2,zone_index)*div
      do j=1,4
      if (j.eq.2) cycle
      q(j,1) = q(j,1) - q(j,2)*q(2,1)
      star%elim_rhs(j,zone_index) = star%elim_rhs(j,zone_index) - q(j,2)*star%elim_rhs(2,zone_index)
      star%elim_coeff(j,1,zone_index) = star%elim_coeff(j,1,zone_index) - &
           q(j,2)*star%elim_coeff(2,1,zone_index)
      star%elim_coeff(j,2,zone_index) = star%elim_coeff(j,2,zone_index) - &
           q(j,2)*star%elim_coeff(2,2,zone_index)
      end do
! PIVOT ON ROW-1 AND COLUMN-3
      div = 1.0D0/q(1,1)
      star%elim_rhs(1,zone_index) = star%elim_rhs(1,zone_index)*div
      star%elim_coeff(1,1,zone_index) = star%elim_coeff(1,1,zone_index)*div
      star%elim_coeff(1,2,zone_index) = star%elim_coeff(1,2,zone_index)*div
      do j=2,4
      star%elim_rhs(j,zone_index) = star%elim_rhs(j,zone_index) - q(j,1)*star%elim_rhs(1,zone_index)
      star%elim_coeff(j,1,zone_index) = star%elim_coeff(j,1,zone_index) - &
           q(j,1)*star%elim_coeff(1,1,zone_index)
      star%elim_coeff(j,2,zone_index) = star%elim_coeff(j,2,zone_index) - &
           q(j,1)*star%elim_coeff(1,2,zone_index)
      end do
      return
end subroutine henyey_eliminate

end module henyey_eliminate_lib
