!----------------------------------------------------------------------
! henyey_elimination (legacy reduce and hsolve)
!----------------------------------------------------------------------
! The linear algebra of the Henyey solve. compute_henyey_coefficients
! (structure_equations.f90) forms the terms of the four structure
! equations at each shell, centre out, and calls henyey_eliminate for
! each pair of shells; henyey_solve then closes the sweep with the
! surface boundary condition and back-substitutes the corrections. For
! the general method see Kippenhahn & Weigert (1990), Stellar Structure
! and Evolution. The rows are the pressure, temperature, radius and
! luminosity equations, in that order.
!
! Conventions, the contract for every henyey_shell_terms record. (No
! original documentation of the terms survives; these were rederived
! from the arithmetic here and in compute_henyey_coefficients in the
! 2026-10 audit.)
!  - Variables y = (log P, log T, log R, L): log10 in cgs, and L linear
!    in Lsun (star%luminosity_lsun). Coordinate: log m, log10 of the
!    enclosed mass in g (star%log_mass).
!  - Each term is the derivative of its variable with respect to log m:
!      qp = dlog P/dlog m, qt = dlog T/dlog m, qr = dlog R/dlog m,
!      ql = dL/dlog m = (ln10 m/Lsun) eps,
!    with eps the total specific energy rate in erg/g/s.
!  - The difference equation for equation k across shell i (between
!    points i-1 and i) is the trapezoid rule, with residual
!      F_k(i) = 0.5 Delta log m (q_k(i-1) + q_k(i)) - (y_k(i) - y_k(i-1)),
!    Delta log m = log m(i) - log m(i-1).
!  - The partials q_dp, q_dt, q_dr are natural-log derivatives dq/dln P,
!    dq/dln T, dq/dln R at fixed values of the other three variables.
!    henyey_eliminate multiplies 0.5 Delta log m by ln10 once, which makes
!    them derivatives with respect to the log10 variables.
!  - The luminosity partial q_dl is dq/dL divided by ln10, so that the
!    same factor leaves it as dq/dL (qt_dl = qt/(L ln10) for a radiative
!    gradient).
!  - Hence ql_dp and ql_dt carry the absolute derivatives d eps/dln P
!    and d eps/dln T (times ln10 m/Lsun), not logarithmic ones; that is
!    what compute_zone_energy returns in burn%deps_dlnrho/deps_dlnt.
!  - Only the partials in henyey_shell_terms exist. The others (qp_dt,
!    qp_dl, qr_dl, ql_dr, ql_dl) are zero, so an energy term that depends
!    on R or L cannot be linearized without extending the type and the
!    pivots below.
!
! The linear system. Newton's method asks for corrections with
! F + J (corrections) = 0. Write dX for the correction of log X (of L
! itself for the luminosity) and h = 0.5 Delta log m ln10. For the pair
! (i-1, i), equation k reads
!   sum_X [-h q_kX(i-1) - delta_kX] dX(i-1)
!     + sum_X [-h q_kX(i) + delta_kX] dX(i) = F_k(i),
! with X over P, T, R, L and q_kX the partial q_k_dX. The elimination
! leaves, at each shell i, four relations that give corrections in terms
! of dlog P(i) and dlog T(i):
!   dlog P(i-1) = e(1,i) - c(1,1,i) dlog P(i) - c(1,2,i) dlog T(i)
!   dlog T(i-1) = e(2,i) - c(2,1,i) dlog P(i) - c(2,2,i) dlog T(i)
!   dlog R(i)   = e(3,i) - c(3,1,i) dlog P(i) - c(3,2,i) dlog T(i)
!   dL(i)       = e(4,i) - c(4,1,i) dlog P(i) - c(4,2,i) dlog T(i)
! with e = star%elim_rhs and c = star%elim_coeff. The central boundary
! conditions (compute_henyey_coefficients) give rows 3 and 4 at i = 1 in
! the same form, and the surface boundary condition (star%surface_bc,
! from henyey_iterate) gives dlog R and dL at the surface.
module henyey_elimination
! The structure-equation terms at one shell (the term q of each of the
! four equations and its nonzero partials) come from the kernel
! evaluate_shell_equations (shell_equations.f90), which defines their
! record.
      use shell_equations, only: henyey_shell_terms
      implicit none

contains

!----------------------------------------------------------------------
! henyey_eliminate (formerly reduce)
!----------------------------------------------------------------------
! One step of the forward elimination, for the pair of points
! i - 1 and i = zone_index, with prev and cur their terms. It forms the
! four residuals F_k(i) (and updates star%max_residual, the largest
! |F_k| over the shells), substitutes the previous shell's relations for
! dlog R(i-1) and dL(i-1), and solves the four equations for dL(i),
! dlog R(i), dlog T(i-1) and dlog P(i-1) in terms of dlog P(i) and
! dlog T(i). The result goes to index i of star%elim_coeff and
! star%elim_rhs (the relations in the module header).
! In the legacy column numbering the unknowns are 1 = dlog R(i-1),
! 2 = dL(i-1), 3 = dlog P(i-1), 4 = dlog T(i-1), 5 = dlog R(i),
! 6 = dL(i); dlog P(i) and dlog T(i) stay on the right as
! elim_coeff(:,1:2,i). q(k,1) and q(k,2) hold columns 3 and 4.
subroutine henyey_eliminate(zone_index, prev, cur)
      use star_info_lib, only: star
      use physical_constants
      implicit none

      integer, intent(in) :: zone_index
      type(henyey_shell_terms), intent(in) :: prev, cur

      double precision :: q(4,4)
      double precision :: half_delta_log_mass, pivot, pivotb, div
      integer :: j

! The residuals
!   F_k(i) = 0.5 Delta log m (q_k(i-1) + q_k(i)) - (y_k(i) - y_k(i-1))
! and their running maxima.
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
! From here half_delta_log_mass is h = 0.5 Delta log m ln10, which turns
! the dln partials into dlog10 ones.
      half_delta_log_mass = half_delta_log_mass*ln10
! Columns 1 and 2: the previous point's coefficients are
! -h q_kX(i-1) - delta_kX. Substitute its relations
!   dlog R(i-1) = e(3,i-1) - c(3,1,i-1) dlog P(i-1) - c(3,2,i-1) dlog T(i-1)
! and the same for dL(i-1) (row 4), moving the e terms to the right.
! Pressure row: dlog R(i-1) enters through qp_dr.
      pivot = -half_delta_log_mass*prev%qp_dr
      q(1,1) = (-1.0D0-half_delta_log_mass*prev%qp_dp) - &
           pivot*star%elim_coeff(3,1,zone_index-1)
      q(1,2) =                  - pivot*star%elim_coeff(3,2,zone_index-1)
      star%elim_rhs(1,zone_index) = star%elim_rhs(1,zone_index) - &
           pivot*star%elim_rhs(3,zone_index-1)
! Temperature row: dlog R(i-1) and dL(i-1) enter through qt_dr and qt_dl.
      pivot = -half_delta_log_mass*prev%qt_dr
      pivotb= -half_delta_log_mass*prev%qt_dl
      q(2,1) = -half_delta_log_mass*prev%qt_dp &
           -pivot*star%elim_coeff(3,1,zone_index-1)-pivotb*star%elim_coeff(4,1,zone_index-1)
      q(2,2) =(-1.0D0-half_delta_log_mass*prev%qt_dt) &
           -pivot*star%elim_coeff(3,2,zone_index-1)-pivotb*star%elim_coeff(4,2,zone_index-1)
      star%elim_rhs(2,zone_index) = star%elim_rhs(2,zone_index) - &
           pivot*star%elim_rhs(3,zone_index-1) - pivotb*star%elim_rhs(4,zone_index-1)
! Radius row: the coefficient of dlog R(i-1) is -1 - h qr_dr.
      pivot = -1.0D0 - half_delta_log_mass*prev%qr_dr
      q(3,1) = -half_delta_log_mass*prev%qr_dp - pivot*star%elim_coeff(3,1,zone_index-1)
      q(3,2) = -half_delta_log_mass*prev%qr_dt - pivot*star%elim_coeff(3,2,zone_index-1)
      star%elim_rhs(3,zone_index) = star%elim_rhs(3,zone_index) - &
           pivot*star%elim_rhs(3,zone_index-1)
! Luminosity row: the coefficient of dL(i-1) is -1 (there is no ql_dl).
      q(4,1) = -half_delta_log_mass*prev%ql_dp + star%elim_coeff(4,1,zone_index-1)
      q(4,2) = -half_delta_log_mass*prev%ql_dt + star%elim_coeff(4,2,zone_index-1)
      star%elim_rhs(4,zone_index) = star%elim_rhs(4,zone_index) + star%elim_rhs(4,zone_index-1)
! The current point's coefficients are -h q_kX(i) + delta_kX. Reduce the
! remaining system for columns 3 to 6.
! Pivot on row 4 (luminosity) and column 6: the coefficient of dL(i) is
! 1, so row 4 already gives dL(i). Remove dL(i) from row 2, the only
! other row it enters (through qt_dl).
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
! Pivot on row 3 (radius) and column 5: divide by the coefficient of
! dlog R(i), 1 - h qr_dr, and remove dlog R(i) from rows 1 and 2
! (through qp_dr and qt_dr).
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
! Pivot on row 2 (temperature) and column 4: divide by the coefficient
! of dlog T(i-1) and remove it from rows 1, 3 and 4.
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
! Pivot on row 1 (pressure) and column 3: divide by the coefficient of
! dlog P(i-1) and remove it from rows 2 to 4. What is left at index i
! is the four relations of the module header.
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

!----------------------------------------------------------------------
! henyey_solve (formerly hsolve)
!----------------------------------------------------------------------
! Closes the forward elimination at the surface and back-substitutes,
! surface to centre. On entry elim_coeff and elim_rhs (star%elim_coeff,
! star%elim_rhs) hold the relations henyey_eliminate and the central
! boundary conditions left at every shell (module header), and
! surface_bc (star%surface_bc, from henyey_iterate's form_surface_bc_rows)
! the outer boundary condition in the same form, with N = num_shells and
! s = surface_bc:
!   dlog R(N) = s(3) - s(1) dlog P(N) - s(2) dlog T(N)
!   dL(N)     = s(6) - s(4) dlog P(N) - s(5) dlog T(N)
! Equating these with rows 3 and 4 at N gives two equations for
! dlog P(N) and dlog T(N). Rows 1 and 2 then give dlog P and dlog T at
! N-1, N-2, ..., 1 in turn, and rows 3 and 4 dlog R and dL at each shell.
! On exit elim_rhs(1:4,i) holds the corrections (dlog P, dlog T, dlog R,
! dL) at shell i. These row roles were derived from the arithmetic here
! and in henyey_eliminate (2026-10), not from original documentation.
subroutine henyey_solve(num_shells, elim_coeff, elim_rhs, surface_bc)
      use star_info_lib, only: json

      implicit none

      integer, intent(in) :: num_shells
      double precision, intent(in) :: elim_coeff(4,2,json)
      double precision, intent(inout) :: elim_rhs(4,json)
      double precision, intent(in) :: surface_bc(6)
      double precision :: pivot_divisor, pivot_ratio, t_coupling_term
      double precision :: correction_p, correction_t
      double precision :: p_store, t_store
      integer :: shell_idx, im

! The surface system, with c = elim_coeff and e = elim_rhs:
!   (s(1) - c(3,1,N)) dlog P + (s(2) - c(3,2,N)) dlog T = s(3) - e(3,N)
!   (s(4) - c(4,1,N)) dlog P + (s(5) - c(4,2,N)) dlog T = s(6) - e(4,N)
! solved by elimination (pivot_ratio = the ratio of the dlog P
! coefficients). elim_rhs(1:2,1) serve as scratch: shell 1 has no rows
! 1 and 2.
      pivot_divisor = 1.0d0/(surface_bc(1) - elim_coeff(3,1,num_shells))
      pivot_ratio = (surface_bc(4) - elim_coeff(4,1,num_shells))*pivot_divisor
      elim_rhs(1,1) = surface_bc(3) - elim_rhs(3,num_shells)
      elim_rhs(2,1) = (surface_bc(6) - elim_rhs(4,num_shells)) - &
           pivot_ratio*elim_rhs(1,1)
      t_coupling_term = surface_bc(2) - elim_coeff(3,2,num_shells)
      correction_t = elim_rhs(2,1)/(surface_bc(5) - &
           elim_coeff(4,2,num_shells) - pivot_ratio*t_coupling_term)
      correction_p = pivot_divisor*(elim_rhs(1,1) - correction_t*t_coupling_term)
! Back-substitute inward, im = N, N-1, ..., 2: with dlog P(im) and
! dlog T(im) known, rows 1 and 2 give them at im-1 (p_store, t_store) and
! rows 3 and 4 give dlog R(im) and dL(im).
      do shell_idx = 2, num_shells
         im = num_shells + 2 - shell_idx
         p_store = elim_rhs(1,im) - correction_p*elim_coeff(1,1,im) - &
              correction_t*elim_coeff(1,2,im)
         t_store = elim_rhs(2,im) - correction_p*elim_coeff(2,1,im) - &
              correction_t*elim_coeff(2,2,im)
         elim_rhs(3,im) = elim_rhs(3,im) - correction_p*elim_coeff(3,1,im) - &
              correction_t*elim_coeff(3,2,im)
         elim_rhs(4,im) = elim_rhs(4,im) - correction_p*elim_coeff(4,1,im) - &
              correction_t*elim_coeff(4,2,im)
         elim_rhs(1,im) = correction_p
         elim_rhs(2,im) = correction_t
         correction_p = p_store
         correction_t = t_store
      end do
! The centre: rows 3 and 4 at shell 1 are the central boundary
! conditions.
      elim_rhs(1,1) = correction_p
      elim_rhs(2,1) = correction_t
      elim_rhs(3,1) = elim_rhs(3,1) - correction_p*elim_coeff(3,1,1) - &
           correction_t*elim_coeff(3,2,1)
      elim_rhs(4,1) = elim_rhs(4,1) - correction_p*elim_coeff(4,1,1) - &
           correction_t*elim_coeff(4,2,1)

      return
end subroutine henyey_solve

end module henyey_elimination
