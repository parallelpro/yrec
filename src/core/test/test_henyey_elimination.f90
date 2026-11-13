!----------------------------------------------------------------------
! test_henyey_elimination
!----------------------------------------------------------------------
! The linear algebra of the Henyey solve (core/structure/
! henyey_elimination.f90) against a dense solve; no star, no tables.
! Sections:
!   solve     systems of N shells (N = 2, 3 and 40) with pseudo-random
!             terms and partials (each shell's henyey_shell_terms),
!             structure variables, mass coordinates, central rows and
!             surface rows. henyey_eliminate over the pairs of shells,
!             then henyey_solve, must give the corrections of the
!             4N x 4N Newton system the module header writes out
!             (assembled here row by row and solved by LU decomposition,
!             numerics' ludcmp and lubksb): the largest difference,
!             relative to the largest correction, within 1e-10. The
!             residual maxima henyey_eliminate keeps must be the largest
!             |F_k| over the shells, exactly.
!   singular  a zero pivot in each of henyey_eliminate's three divisions
!             (the radius row, then the temperature and pressure rows
!             of the reduced system) and in each of the two of
!             henyey_solve's surface system returns ierr =
!             numerics_singular, with no division by zero.
! The numbers come from a fixed linear congruential sequence, so every
! platform builds the same systems. Each check prints its value (1pe9.1)
! and ok or FAIL; any FAIL stops the program with status 1.
program test_henyey_elimination
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use henyey_elimination, only: henyey_eliminate, henyey_solve
      use shell_equations, only: henyey_shell_terms
      use numerics_lib, only: ludcmp, lubksb, numerics_singular
      use physical_constants, only: ln10
      implicit none

      integer :: nfail
! the state of the pseudo-random sequence
      integer(8) :: seed

      nfail = 0
      seed = 20261112_8
      write(*,'(a)') 'test_henyey_elimination report'
      write(*,'(a)') '-- solve'
      call check_solve(2)
      call check_solve(3)
      call check_solve(40)
      write(*,'(a)') '-- singular'
      call check_singular_eliminate()
      call check_singular_solve()
      if (nfail > 0) then
         write(*,'(a,i0,a)') 'test_henyey_elimination: ', nfail, ' FAIL'
         stop 1
      end if
      write(*,'(a)') 'test_henyey_elimination: done'

contains

! ---------------------------------------------------------------
! A number in [lo, hi) from the sequence (Knuth's MMIX constants).
real(dp) function uniform(lo, hi)
      real(dp), intent(in) :: lo, hi
      seed = 6364136223846793005_8*seed + 1442695040888963407_8
      uniform = lo + (hi - lo)*real(ishft(seed, -11), dp)/2.0_dp**53
end function uniform

! ---------------------------------------------------------------
! One shell's terms and partials, each in [-1, 1).
type(henyey_shell_terms) function random_terms()
      random_terms%qp = uniform(-1.0_dp, 1.0_dp)
      random_terms%qp_dr = uniform(-1.0_dp, 1.0_dp)
      random_terms%qp_dp = uniform(-1.0_dp, 1.0_dp)
      random_terms%qt = uniform(-1.0_dp, 1.0_dp)
      random_terms%qt_dr = uniform(-1.0_dp, 1.0_dp)
      random_terms%qt_dl = uniform(-1.0_dp, 1.0_dp)
      random_terms%qt_dp = uniform(-1.0_dp, 1.0_dp)
      random_terms%qt_dt = uniform(-1.0_dp, 1.0_dp)
      random_terms%qr = uniform(-1.0_dp, 1.0_dp)
      random_terms%qr_dr = uniform(-1.0_dp, 1.0_dp)
      random_terms%qr_dp = uniform(-1.0_dp, 1.0_dp)
      random_terms%qr_dt = uniform(-1.0_dp, 1.0_dp)
      random_terms%ql = uniform(-1.0_dp, 1.0_dp)
      random_terms%ql_dp = uniform(-1.0_dp, 1.0_dp)
      random_terms%ql_dt = uniform(-1.0_dp, 1.0_dp)
end function random_terms

! ---------------------------------------------------------------
! The partial of term k (1 P, 2 T, 3 R, 4 L) with respect to variable v
! (the same order), zero where henyey_shell_terms has none.
real(dp) function partial(q, k, v)
      type(henyey_shell_terms), intent(in) :: q
      integer, intent(in) :: k, v
      partial = 0.0_dp
      select case (10*k + v)
      case (11); partial = q%qp_dp
      case (13); partial = q%qp_dr
      case (21); partial = q%qt_dp
      case (22); partial = q%qt_dt
      case (23); partial = q%qt_dr
      case (24); partial = q%qt_dl
      case (31); partial = q%qr_dp
      case (32); partial = q%qr_dt
      case (33); partial = q%qr_dr
      case (41); partial = q%ql_dp
      case (42); partial = q%ql_dt
      end select
end function partial

! ---------------------------------------------------------------
! Term k of q.
real(dp) function term(q, k)
      type(henyey_shell_terms), intent(in) :: q
      integer, intent(in) :: k
      select case (k)
      case (1); term = q%qp
      case (2); term = q%qt
      case (3); term = q%qr
      case default; term = q%ql
      end select
end function term

! ---------------------------------------------------------------
! An ok/FAIL line.
subroutine report(name, value, pass)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: value
      logical, intent(in) :: pass
      if (pass) then
         write(*,'(a,t58,1pe9.1,2x,a)') name, value, 'ok'
      else
         write(*,'(a,t58,1pe9.1,2x,a)') name, value, 'FAIL'
         nfail = nfail + 1
      end if
end subroutine report

! ---------------------------------------------------------------
! A random system of n shells: the elimination and the back
! substitution against the dense solve, and the residual maxima.
subroutine check_solve(n)
      integer, intent(in) :: n
      type(henyey_shell_terms) :: q(n)
      real(dp) :: log_mass(n), y(4,n), coeff(4,2,n), rhs(4,n), &
           surface_bc(6), max_residual(4), expected_max(4)
      real(dp) :: a(4*n,4*n), b(4*n), d, h, half_dlm, f, worst, scale
      integer :: indx(4*n), i, k, v, row, ierr, jerr
      character(len=60) :: name

      do i = 1, n
         q(i) = random_terms()
         log_mass(i) = 30.0_dp + 0.1_dp*i + uniform(0.0_dp, 0.05_dp)
         do v = 1, 4
            y(v,i) = uniform(-2.0_dp, 2.0_dp)
         end do
      end do
      coeff = 0.0_dp
      rhs = 0.0_dp
! The central rows 3 and 4 at shell 1, and the surface rows.
      do k = 3, 4
         coeff(k,1,1) = uniform(-1.0_dp, 1.0_dp)
         coeff(k,2,1) = uniform(-1.0_dp, 1.0_dp)
         rhs(k,1) = uniform(-1.0_dp, 1.0_dp)
      end do
      do k = 1, 6
         surface_bc(k) = uniform(-1.0_dp, 1.0_dp)
      end do

! The dense system, unknown (v, i) in column 4(i-1) + v. The central
! rows: dlog R(1) + c(3,1,1) dlog P(1) + c(3,2,1) dlog T(1) = e(3,1), and
! the same for dL(1) with row 4.
      a = 0.0_dp
      b = 0.0_dp
      do k = 3, 4
         row = k - 2
         a(row, k) = 1.0_dp
         a(row, 1) = coeff(k,1,1)
         a(row, 2) = coeff(k,2,1)
         b(row) = rhs(k,1)
      end do
! The pairs (i-1, i): sum_X [-h q_kX(i-1) - delta_kX] dX(i-1)
!   + sum_X [-h q_kX(i) + delta_kX] dX(i) = F_k(i), with
! h = 0.5 Delta log m ln10 and
! F_k(i) = 0.5 Delta log m (q_k(i-1) + q_k(i)) - (y_k(i) - y_k(i-1)).
      expected_max = 0.0_dp
      do i = 2, n
         half_dlm = 0.5_dp*(log_mass(i) - log_mass(i-1))
         h = half_dlm*ln10
         do k = 1, 4
            row = 2 + 4*(i-2) + k
            do v = 1, 4
               a(row, 4*(i-2)+v) = -h*partial(q(i-1), k, v)
               a(row, 4*(i-1)+v) = -h*partial(q(i), k, v)
            end do
            a(row, 4*(i-2)+k) = a(row, 4*(i-2)+k) - 1.0_dp
            a(row, 4*(i-1)+k) = a(row, 4*(i-1)+k) + 1.0_dp
            f = half_dlm*(term(q(i-1), k) + term(q(i), k)) - (y(k,i) - y(k,i-1))
            b(row) = f
            expected_max(k) = max(expected_max(k), abs(f))
         end do
      end do
! The surface rows: dlog R(N) + s(1) dlog P(N) + s(2) dlog T(N) = s(3),
! dL(N) + s(4) dlog P(N) + s(5) dlog T(N) = s(6).
      row = 4*n - 1
      a(row, 4*(n-1)+3) = 1.0_dp
      a(row, 4*(n-1)+1) = surface_bc(1)
      a(row, 4*(n-1)+2) = surface_bc(2)
      b(row) = surface_bc(3)
      row = 4*n
      a(row, 4*(n-1)+4) = 1.0_dp
      a(row, 4*(n-1)+1) = surface_bc(4)
      a(row, 4*(n-1)+2) = surface_bc(5)
      b(row) = surface_bc(6)
      call ludcmp(a, 4*n, 4*n, indx, d, jerr)
      if (jerr /= 0) then
         write(*,'(a)') 'test_henyey_elimination: the dense system is singular'
         stop 1
      end if
      call lubksb(a, 4*n, 4*n, indx, b)

! The Henyey solve.
      max_residual = 0.0_dp
      ierr = 0
      do i = 2, n
         call henyey_eliminate(i, log_mass, y(1,:), y(2,:), y(3,:), y(4,:), &
              q(i-1), q(i), coeff, rhs, max_residual, jerr)
         if (jerr /= 0) ierr = jerr
      end do
      if (ierr == 0) call henyey_solve(n, coeff, rhs, surface_bc, ierr)
      write(name,'(a,i0,a)') 'N = ', n, ': no error'
      call report(trim(name), real(ierr, dp), ierr == 0)
      worst = 0.0_dp
      scale = maxval(abs(b))
      do i = 1, n
         do v = 1, 4
            worst = max(worst, abs(rhs(v,i) - b(4*(i-1)+v)))
         end do
      end do
      write(name,'(a,i0,a)') 'N = ', n, ': corrections = dense solve'
      call report(trim(name), worst/scale, worst <= 1.0e-10_dp*scale)
      write(name,'(a,i0,a)') 'N = ', n, ': residual maxima = max |F_k|'
      call report(trim(name), maxval(abs(max_residual - expected_max)), &
           all(max_residual == expected_max))
end subroutine check_solve

! ---------------------------------------------------------------
! henyey_eliminate's three pivots, each made exactly zero in turn on a
! pair of shells whose partials are zero but for the one under test,
! with h = 0.5 Delta log m ln10 = 1/4 exactly (Delta log m chosen so):
! the radius row's 1 - h qr_dr(i) (qr_dr(i) = 4); the temperature row's
! coefficient of dlog T(i-1), -1 - h qt_dt(i-1) (qt_dt(i-1) = -4); the
! pressure row's, -1 - h qp_dp(i-1) (qp_dp(i-1) = -4).
subroutine check_singular_eliminate()
      type(henyey_shell_terms) :: prev, cur
      real(dp) :: log_mass(2), y(4,2), coeff(4,2,2), rhs(4,2), max_residual(4)
      integer :: ierr, which, k
      character(len=60), parameter :: names(3) = [character(len=60) :: &
           'eliminate: zero radius pivot', 'eliminate: zero temperature pivot', &
           'eliminate: zero pressure pivot']

! The Delta log m nearest 0.5/ln10 that makes h exactly 1/4.
      log_mass = [0.0_dp, 0.5_dp/ln10]
      do k = 1, 64
         if ((0.5_dp*(log_mass(2) - log_mass(1)))*ln10 == 0.25_dp) exit
         if ((0.5_dp*(log_mass(2) - log_mass(1)))*ln10 < 0.25_dp) then
            log_mass(2) = nearest(log_mass(2), 1.0_dp)
         else
            log_mass(2) = nearest(log_mass(2), -1.0_dp)
         end if
      end do
      if ((0.5_dp*(log_mass(2) - log_mass(1)))*ln10 /= 0.25_dp) then
         write(*,'(a)') 'test_henyey_elimination: no Delta log m gives h = 1/4'
         stop 1
      end if
      y = 1.0_dp
      do which = 1, 3
         prev = henyey_shell_terms(0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, &
              0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
              0.5_dp, 0.0_dp, 0.0_dp)
         cur = prev
         select case (which)
         case (1); cur%qr_dr = 4.0_dp
         case (2); prev%qt_dt = -4.0_dp
         case (3); prev%qp_dp = -4.0_dp
         end select
         coeff = 0.0_dp
         rhs = 0.0_dp
         max_residual = 0.0_dp
         call henyey_eliminate(2, log_mass, y(1,:), y(2,:), y(3,:), y(4,:), &
              prev, cur, coeff, rhs, max_residual, ierr)
         call report(trim(names(which)), real(ierr, dp), ierr == numerics_singular)
      end do
end subroutine check_singular_eliminate

! ---------------------------------------------------------------
! henyey_solve's surface system: its dlog P pivot s(1) - c(3,1,N) zero,
! then its dlog T pivot s(5) - c(4,2,N) - r (s(2) - c(3,2,N)) zero (r the
! ratio of the dlog P coefficients).
subroutine check_singular_solve()
      real(dp) :: coeff(4,2,3), rhs(4,3), surface_bc(6)
      integer :: ierr

      coeff = 0.5_dp
      rhs = 1.0_dp
      surface_bc = [0.5_dp, 0.25_dp, 1.0_dp, 0.75_dp, 0.5_dp, 1.0_dp]
      call henyey_solve(3, coeff, rhs, surface_bc, ierr)
      call report('solve: zero dlog P pivot', real(ierr, dp), ierr == numerics_singular)
! With every c = 0.5: s(1) - c = 0.25, r = (s(4) - c)/0.25 = 1, and
! s(5) - c - r (s(2) - c) = 0.25 - 0.5 - (0.25 - 0.5) = 0.
      coeff = 0.5_dp
      rhs = 1.0_dp
      surface_bc = [0.75_dp, 0.25_dp, 1.0_dp, 0.75_dp, 0.25_dp, 1.0_dp]
      call henyey_solve(3, coeff, rhs, surface_bc, ierr)
      call report('solve: zero dlog T pivot', real(ierr, dp), ierr == numerics_singular)
end subroutine check_singular_solve

end program test_henyey_elimination
