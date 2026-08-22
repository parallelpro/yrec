!----------------------------------------------------------------------
! reduce
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original reduce.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Classic Henyey-method forward elimination step: given the structure
! equation residuals and their partial derivatives (with respect to
! P, T, R, L) at the two mesh points bracketing shell I, eliminates
! the previous shell's (I-1) unknowns using the recursion coefficients
! carried in elim_coeff/elim_rhs from the prior call, then reduces the
! resulting 4x4 system for shell I down to expressions for P and T in
! terms of R and L (stored back into elim_coeff/elim_rhs at index I for
! use by the next shell). See Kippenhahn & Weigert-style Henyey solvers
! for the general method; the four rows correspond to the pressure,
! temperature, radius, and luminosity structure equations in that
! order, and the interpretation of the individual Q* arguments below
! (equation value and its partials at the lower/upper mesh point) is
! inferred from their use, not from surviving documentation.
subroutine reduce(zone_index,elim_coeff,elim_rhs,log_luminosity,max_residual, &
     log_pressure,log_radius,log_mass,log_temperature, &
     eq_p_val0,eq_p_val,dqp_dr0,dqp_dr, &
     dqp_dp0,dqp_dp,eq_t_val0,eq_t_val,dqt_dr0,dqt_dr,dqt_dl0,dqt_dl, &
     dqt_dp0,dqt_dp,dqt_dt0,dqt_dt,eq_r_val0,eq_r_val,dqr_dr0,dqr_dr, &
     dqr_dp0,dqr_dp,dqr_dt0,dqr_dt,eq_l_val0,eq_l_val,dql_dp0,dql_dp, &
     dql_dt0,dql_dt)

      use const_lib
      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: zone_index
      double precision, intent(inout) :: elim_coeff(4,2,json), elim_rhs(4,json)
      double precision, intent(in) :: log_luminosity(json)
      double precision, intent(inout) :: max_residual(4)
      double precision, intent(in) :: log_pressure(json), log_radius(json), &
           log_mass(json), log_temperature(json)
      double precision, intent(in) :: eq_p_val0,eq_p_val,dqp_dr0,dqp_dr
      double precision, intent(in) :: dqp_dp0,dqp_dp,eq_t_val0,eq_t_val, &
           dqt_dr0,dqt_dr,dqt_dl0,dqt_dl
      double precision, intent(in) :: dqt_dp0,dqt_dp,dqt_dt0,dqt_dt, &
           eq_r_val0,eq_r_val,dqr_dr0,dqr_dr
      double precision, intent(in) :: dqr_dp0,dqr_dp,dqr_dt0,dqr_dt, &
           eq_l_val0,eq_l_val,dql_dp0,dql_dp
      double precision, intent(in) :: dql_dt0,dql_dt


      double precision :: q(4,4)
      double precision :: half_delta_log_mass, pivot, pivotb, div
      integer :: j

! DEFINE R.H.S
      half_delta_log_mass = 0.5D0*(log_mass(zone_index)-log_mass(zone_index-1))
      elim_rhs(1,zone_index) = half_delta_log_mass*(eq_p_val0 + eq_p_val) - &
           (log_pressure(zone_index) - log_pressure(zone_index-1))
      max_residual(1) = dmax1(max_residual(1),dabs(elim_rhs(1,zone_index)))
      elim_rhs(2,zone_index) = half_delta_log_mass*(eq_t_val0 + eq_t_val) - &
           (log_temperature(zone_index) - log_temperature(zone_index-1))
      max_residual(2) = dmax1(max_residual(2),dabs(elim_rhs(2,zone_index)))
      elim_rhs(3,zone_index) = half_delta_log_mass*(eq_r_val0 + eq_r_val) - &
           (log_radius(zone_index) - log_radius(zone_index-1))
      max_residual(3) = dmax1(max_residual(3),dabs(elim_rhs(3,zone_index)))
      elim_rhs(4,zone_index) = half_delta_log_mass*(eq_l_val0 + eq_l_val) - &
           (log_luminosity(zone_index) - log_luminosity(zone_index-1))
      max_residual(4) = dmax1(max_residual(4),dabs(elim_rhs(4,zone_index)))
      half_delta_log_mass = half_delta_log_mass*ln10
! ELIMINATE COLUMNS 1 AND 2
      pivot = -half_delta_log_mass*dqp_dr0
      q(1,1) = (-1.0D0-half_delta_log_mass*dqp_dp0) - &
           pivot*elim_coeff(3,1,zone_index-1)
      q(1,2) =                  - pivot*elim_coeff(3,2,zone_index-1)
      elim_rhs(1,zone_index) = elim_rhs(1,zone_index) - &
           pivot*elim_rhs(3,zone_index-1)
      pivot = -half_delta_log_mass*dqt_dr0
      pivotb= -half_delta_log_mass*dqt_dl0
      q(2,1) = -half_delta_log_mass*dqt_dp0 &
           -pivot*elim_coeff(3,1,zone_index-1)-pivotb*elim_coeff(4,1,zone_index-1)
      q(2,2) =(-1.0D0-half_delta_log_mass*dqt_dt0) &
           -pivot*elim_coeff(3,2,zone_index-1)-pivotb*elim_coeff(4,2,zone_index-1)
      elim_rhs(2,zone_index) = elim_rhs(2,zone_index) - &
           pivot*elim_rhs(3,zone_index-1) - pivotb*elim_rhs(4,zone_index-1)
      pivot = -1.0D0 - half_delta_log_mass*dqr_dr0
      q(3,1) = -half_delta_log_mass*dqr_dp0 - pivot*elim_coeff(3,1,zone_index-1)
      q(3,2) = -half_delta_log_mass*dqr_dt0 - pivot*elim_coeff(3,2,zone_index-1)
      elim_rhs(3,zone_index) = elim_rhs(3,zone_index) - &
           pivot*elim_rhs(3,zone_index-1)
      q(4,1) = -half_delta_log_mass*dql_dp0 + elim_coeff(4,1,zone_index-1)
      q(4,2) = -half_delta_log_mass*dql_dt0 + elim_coeff(4,2,zone_index-1)
      elim_rhs(4,zone_index) = elim_rhs(4,zone_index) + elim_rhs(4,zone_index-1)
! REDUCE 4*4 MATRIX
! PIVOT ON ROW-4 AND COLUMN-6
      elim_coeff(4,1,zone_index) = -half_delta_log_mass*dql_dp
      elim_coeff(4,2,zone_index) = -half_delta_log_mass*dql_dt
      pivot = -half_delta_log_mass*dqt_dl
      q(2,1) = q(2,1) - pivot*q(4,1)
      q(2,2) = q(2,2) - pivot*q(4,2)
      elim_rhs(2,zone_index) = elim_rhs(2,zone_index) - pivot*elim_rhs(4,zone_index)
      elim_coeff(2,1,zone_index) = -half_delta_log_mass*dqt_dp      - &
           pivot*elim_coeff(4,1,zone_index)
      elim_coeff(2,2,zone_index) = (1.0D0-half_delta_log_mass*dqt_dt) - &
           pivot*elim_coeff(4,2,zone_index)
! PIVOT ON ROW-3 AND COLUMN-5
      div = 1.0D0/(1.0D0 - half_delta_log_mass*dqr_dr)
      q(3,1) = q(3,1)*div
      q(3,2) = q(3,2)*div
      elim_rhs(3,zone_index) = elim_rhs(3,zone_index)*div
      div = -half_delta_log_mass*div
      elim_coeff(3,1,zone_index) = div*dqr_dp
      elim_coeff(3,2,zone_index) = div*dqr_dt
      pivot = -half_delta_log_mass*dqp_dr
      q(1,1) = q(1,1) - pivot*q(3,1)
      q(1,2) = q(1,2) - pivot*q(3,2)
      elim_rhs(1,zone_index) = elim_rhs(1,zone_index) - pivot*elim_rhs(3,zone_index)
      elim_coeff(1,1,zone_index) = (1.0D0-half_delta_log_mass*dqp_dp) - &
           pivot*elim_coeff(3,1,zone_index)
      elim_coeff(1,2,zone_index) =                - pivot*elim_coeff(3,2,zone_index)
      pivot = -half_delta_log_mass*dqt_dr
      q(2,1) = q(2,1) - pivot*q(3,1)
      q(2,2) = q(2,2) - pivot*q(3,2)
      elim_rhs(2,zone_index) = elim_rhs(2,zone_index) - pivot*elim_rhs(3,zone_index)
      elim_coeff(2,1,zone_index) = elim_coeff(2,1,zone_index) - &
           pivot*elim_coeff(3,1,zone_index)
      elim_coeff(2,2,zone_index) = elim_coeff(2,2,zone_index) - &
           pivot*elim_coeff(3,2,zone_index)
! PIVOT ON ROW-2 AND COLUMN-4
      div = 1.0D0/q(2,2)
      q(2,1) = q(2,1)*div
      elim_rhs(2,zone_index) = elim_rhs(2,zone_index)*div
      elim_coeff(2,1,zone_index) = elim_coeff(2,1,zone_index)*div
      elim_coeff(2,2,zone_index) = elim_coeff(2,2,zone_index)*div
      do j=1,4
      if (j.eq.2) cycle
      q(j,1) = q(j,1) - q(j,2)*q(2,1)
      elim_rhs(j,zone_index) = elim_rhs(j,zone_index) - q(j,2)*elim_rhs(2,zone_index)
      elim_coeff(j,1,zone_index) = elim_coeff(j,1,zone_index) - &
           q(j,2)*elim_coeff(2,1,zone_index)
      elim_coeff(j,2,zone_index) = elim_coeff(j,2,zone_index) - &
           q(j,2)*elim_coeff(2,2,zone_index)
    2 continue
      end do
! PIVOT ON ROW-1 AND COLUMN-3
      div = 1.0D0/q(1,1)
      elim_rhs(1,zone_index) = elim_rhs(1,zone_index)*div
      elim_coeff(1,1,zone_index) = elim_coeff(1,1,zone_index)*div
      elim_coeff(1,2,zone_index) = elim_coeff(1,2,zone_index)*div
      do j=2,4
      elim_rhs(j,zone_index) = elim_rhs(j,zone_index) - q(j,1)*elim_rhs(1,zone_index)
      elim_coeff(j,1,zone_index) = elim_coeff(j,1,zone_index) - &
           q(j,1)*elim_coeff(1,1,zone_index)
      elim_coeff(j,2,zone_index) = elim_coeff(j,2,zone_index) - &
           q(j,1)*elim_coeff(1,2,zone_index)
    4 continue
      end do
      return
end subroutine reduce
