!----------------------------------------------------------------------
! test_math
!----------------------------------------------------------------------
! math_lib's contract under the backend it was built with (the build's
! own: test_math; crmath: test_math_crmath, when the MESA SDK is there).
! Each check prints ok or FAIL; the program stops with status 1 on any
! FAIL. Checks:
!   integer paths  under crmath, pow(x, n) for integer-valued n (real
!                  and integer exponents) is x*x*...*x, and 1/(x*...*x)
!                  for n < 0, bit for bit, and exp10(n) is exact for
!                  0 <= n <= 22 and 1/10**|n| for -22 <= n < 0; under
!                  the intrinsic backend they are x**y, x**n and 10**x;
!   real paths     pow(x, y) and exp10(x) against x**y and 10**x to a
!                  few ulp of y ln x, the composition's error (bit for
!                  bit under the intrinsic backend);
!                  a single-precision exponent is its exact binary value
!                  promoted (0.6666667 is 0.66666668653488159...), not
!                  its decimal 0.6666667d0: pow(x, y_sp) matches
!                  x**dble(y_sp) to the composition's error (bit for
!                  bit under the intrinsic backend) and differs from
!                  x**0.6666667d0 by about 2e-8;
!   negative bases pow(-2, 3) = -8, pow(-2, 101) = -2**101 and
!                  pow(-2, 100) = 2**100; pow(-2, 0.5) is NaN;
!   pow(0, y)      0 for y > 0 under both backends; for y = 0 and y < 0
!                  the intrinsic backend gives 1 and +Inf, crmath 0
!                  (MESA's definition), which the check expects per
!                  backend;
!   large y        an exponent beyond the default integer range (3e9)
!                  takes the composition and stays finite.
! The NaN and the infinity of the negative-base and zero-base checks are
! made on purpose, so those checks run with the invalid, divide-by-zero
! and overflow traps off (a -ffpe-trap build, make devel_tests, sets
! them) and the halting mode is restored after them.
program test_math
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
      use, intrinsic :: ieee_exceptions, only: ieee_usual, ieee_get_halting_mode, &
           ieee_set_halting_mode
      use math_lib
      implicit none

      integer :: failures, n, i
      real(dp) :: x, product, expected, y
      logical :: holds, crmath_backend, halting(size(ieee_usual))
      real(dp), parameter :: ulp_tolerance = 4.0d0*epsilon(1.0d0)

      failures = 0
      crmath_backend = math_backend == 'crmath'
      write(*,'(a,a)') "math_lib backend: ", math_backend

! integer paths: under crmath repeated multiplication (exact powers of
! ten up to 1e22), under the intrinsic backend x**y and x**n verbatim
      x = 1.0123456789d0
      if (crmath_backend) then
         holds = .true.
         do n = -12, 12
            product = 1.0d0
            do i = 1, abs(n)
               product = product*x
            end do
            if (n < 0) product = 1.0d0/product
            holds = holds .and. pow(x, dble(n)) == product .and. pow(x, n) == product
         end do
         call report_true("pow: integer exponents are repeated multiplication", holds)
         holds = .true.
         do n = 0, 22
            expected = 1.0d0
            do i = 1, n
               expected = expected*10.0d0
            end do
            holds = holds .and. exp10(dble(n)) == expected .and. &
                 exp10(-dble(n)) == 1.0d0/expected
         end do
         call report_true("exp10: exact for 0 <= n <= 22, 1/10**n below", holds)
      else
         holds = .true.
         do n = -12, 12
            holds = holds .and. pow(x, dble(n)) == x**dble(n) .and. pow(x, n) == x**n &
                 .and. exp10(dble(n)) == 10.0d0**dble(n)
         end do
         call report_true("pow, exp10: integer exponents are x**y, x**n, 10**x", holds)
      end if

! real paths
      holds = .true.
      do i = 1, 20
         x = 0.37d0*i
         y = -2.3d0 + 0.41d0*i
         holds = holds .and. abs(pow(x, y) - x**y) <= &
              ulp_tolerance*max(1.0d0, abs(y*log(x)))*abs(x**y)
         holds = holds .and. abs(exp10(y) - 10.0d0**y) <= &
              ulp_tolerance*max(1.0d0, abs(y*log(10.0d0)))*10.0d0**y
      end do
      call report_true("pow, exp10: x**y and 10**x to a few ulp of y ln x", holds)
      if (.not. crmath_backend) then
         holds = .true.
         do i = 1, 20
            x = 0.37d0*i
            y = -2.3d0 + 0.41d0*i
            holds = holds .and. pow(x, y) == x**y .and. exp10(y) == 10.0d0**y
         end do
         call report_true("pow, exp10: the intrinsic backend is x**y, bit for bit", holds)
      end if
      y = dble(0.6666667)
      if (crmath_backend) then
         holds = abs(pow(3.7d0, 0.6666667) - 3.7d0**y) <= ulp_tolerance*abs(3.7d0**y)
      else
         holds = pow(3.7d0, 0.6666667) == 3.7d0**y
      end if
      holds = holds .and. abs(pow(3.7d0, 0.6666667) - 3.7d0**0.6666667d0) > 1.0d-8*3.7d0**y
      call report_true("pow: a single-precision exponent is its binary value promoted", holds)

! negative bases
      call report_true("pow: (-2)**3 = -8", pow(-2.0d0, 3.0d0) == -8.0d0)
      call report_true("pow: (-2)**101 and (-2)**100 to a few ulp of y ln 2", &
           abs(pow(-2.0d0, 101.0d0) + 2.0d0**101) <= 101*ulp_tolerance*2.0d0**101 .and. &
           abs(pow(-2.0d0, 100.0d0) - 2.0d0**100) <= 100*ulp_tolerance*2.0d0**100)
      call ieee_get_halting_mode(ieee_usual, halting)
      call ieee_set_halting_mode(ieee_usual, .false.)
      call report_true("pow: (-2)**0.5 is NaN", ieee_is_nan(pow(-2.0d0, 0.5d0)))

! pow(0, y)
      call report_true("pow: 0**y = 0 for y > 0", pow(0.0d0, 2.5d0) == 0.0d0 .and. &
           pow(0.0d0, 2.0d0) == 0.0d0 .and. pow(0.0d0, 3) == 0.0d0)
      if (crmath_backend) then
         call report_true("pow: 0**0 and 0**(-1) are 0 under crmath (MESA)", &
              pow(0.0d0, 0.0d0) == 0.0d0 .and. pow(0.0d0, -1.0d0) == 0.0d0)
      else
         call report_true("pow: 0**0 = 1 and 0**(-1) = +Inf under the intrinsic backend", &
              pow(0.0d0, 0.0d0) == 1.0d0 .and. .not. ieee_is_finite(pow(0.0d0, -1.0d0)))
      end if
      call ieee_set_halting_mode(ieee_usual, halting)

! an exponent beyond the default integer range
      x = 1.0d0 + 1.0d-10
      call report_true("pow: y = 3e9 takes the composition, finite", &
           ieee_is_finite(pow(x, 3.0d9)) .and. &
           abs(pow(x, 3.0d9) - x**3.0d9) <= 1.0d-6*x**3.0d9)

      if (failures > 0) then
         write(*,'(i0,a)') failures, " check(s) failed"
         stop 1
      end if

contains

subroutine report_true(label, holds)
      character(len=*), intent(in) :: label
      logical, intent(in) :: holds
      if (holds) then
         write(*,'(a,t66,a)') label, "ok"
      else
         write(*,'(a,t66,a)') label, "FAIL"
         failures = failures + 1
      end if
end subroutine report_true

end program test_math
