!----------------------------------------------------------------------
! math_lib
!----------------------------------------------------------------------
! The one place YREC gets its elementary math (2026, reproducibility
! campaign; modeled on MESA's math_lib; src/math/README.md). Two
! backends, chosen at build time:
!
!   make USE_CRMATH=1   correctly-rounded backend: crmath (Townsend's
!                       Fortran wrapper over crlibm, shipped in the
!                       MESA SDK). Every elementary function returns
!                       the IEEE-nearest result, so results are
!                       bit-identical across platforms and architectures
!                       with -ffp-contract=off (set unconditionally in
!                       the Makefile) and the same gfortran major
!                       version.
!   make                intrinsic backend: pow and exp10 are ** and the
!                       compiler's functions, verbatim -- the historical
!                       behaviour bit for bit, so converting a call site
!                       leaves the pins unchanged. Fast, and
!                       platform-dependent in the last bits.
!
! Usage contract (enforced by tools/check_boundaries.py):
!   * every procedure that calls an elementary transcendental
!     (exp, log, log10, sin, cos, ...) has `use math_lib` in scope:
!     under the crmath backend the module procedures shadow the
!     intrinsics, so call sites are unchanged;
!   * no real exponent after **: x**y with real y is a hidden libm pow,
!     written pow(x, y); 10**x is written exp10(x). Integer exponents
!     (x**2, x**n with integer n) stay as **.
!
! Under crmath, pow and exp10 follow MESA's definitions: an
! integer-valued exponent below 100 in magnitude takes repeated
! multiplication, anything else pow(x,y) = exp(log(x)*y) and exp10(x) =
! exp(x*ln10), which crmath's exp and log determine completely. Under
! the intrinsic backend they are passthroughs (x**y, 10**x) with no
! fast path. Where the two differ (pow(0, y <= 0), integer exponents,
! negative bases) is listed in src/math/README.md and checked by
! test_math.
module math_lib
      use, intrinsic :: iso_fortran_env, only: dp => real64
#ifdef YREC_CRMATH
      use crmath
#endif
      implicit none
      private

      public :: pow, exp10

! pow is generic over real and integer exponents: pow(x, iy) with
! integer iy exists so ** conversions never need to know the
! exponent's type (integer-typed sites take the repeated
! multiplication; under the intrinsic backend it is x**iy verbatim).
      interface pow
         module procedure pow_r
         module procedure pow_i
         module procedure pow_r_sp
      end interface pow
#ifdef YREC_CRMATH
! re-export crmath's correctly-rounded shadows of the intrinsics for
! every elementary function YREC uses (plus the common siblings) so a
! consumer only ever needs `use math_lib`
      public :: exp, expm1, log, log1p, log10, log2
      public :: sin, cos, tan, asin, acos, atan, atan2
      public :: sinh, cosh, tanh, hypot
      character(len=*), parameter, public :: math_backend = 'crmath'

! ln(10) as the correctly-rounded double literal (cannot be written
! log(10.0d0) here: with crmath use-associated that names a module
! function, which is not a constant expression -- and gfortran 15.2
! ICEs on it rather than diagnosing it)
      real(dp), parameter :: ln10 = 2.302585092994045684d0
#else
      character(len=*), parameter, public :: math_backend = 'intrinsic'
#endif

contains

#ifdef YREC_CRMATH

! ---------------------------------------------------------------
! x**y for real y (MESA's definition). Integer-valued |y| < 100 uses
! repeated multiplication; otherwise the exp(log(x)*y) composition,
! fully determined by crmath, with the sign of (-1)**y for a negative x
! and an integer-valued y (where the composition alone would take the
! log of a negative number). x = 0 returns 0 for every y, as in MESA
! (the intrinsic x**y gives 1 for y = 0 and +Inf for y < 0).
elemental function pow_r(x, y) result(pow_x)
      real(dp), intent(in) :: x, y
      real(dp) :: pow_x
      integer :: iy, i

      if (x == 0.0d0) then
         pow_x = 0.0d0
         return
      end if
! the test on |y| comes first, so that int(y) cannot overflow
      if (abs(y) < 100.0d0 .and. y == aint(y)) then
         iy = int(y)
         pow_x = 1.0d0
         do i = 1, abs(iy)
            pow_x = pow_x*x
         end do
         if (iy < 0) pow_x = 1.0d0/pow_x
      else if (x < 0.0d0 .and. y == aint(y)) then
         pow_x = exp(log(-x)*y)
         if (mod(y, 2.0d0) /= 0.0d0) pow_x = -pow_x
      else
         pow_x = exp(log(x)*y)
      end if
end function pow_r

! integer exponent: repeated multiplication (MESA's pow_i); x = 0
! returns 0 for every iy, as pow_r does
elemental function pow_i(x, iy) result(pow_x)
      real(dp), intent(in) :: x
      integer, intent(in) :: iy
      real(dp) :: pow_x
      integer :: i

      if (x == 0.0d0) then
         pow_x = 0.0d0
         return
      end if
      pow_x = 1.0d0
      do i = 1, abs(iy)
         pow_x = pow_x*x
      end do
      if (iy < 0) pow_x = 1.0d0/pow_x
end function pow_i

! single-precision exponent (legacy literals like 0.6666667): Fortran's
! mixed-mode ** promotes the exponent to real(dp), so pow_r(x,
! real(y, dp)) matches x**y exactly
elemental function pow_r_sp(x, y) result(pow_x)
      real(dp), intent(in) :: x
      real, intent(in) :: y
      real(dp) :: pow_x
      pow_x = pow_r(x, real(y, dp))
end function pow_r_sp

! ---------------------------------------------------------------
! 10**x for real x (MESA's definition): for an integer-valued |x| < 100,
! repeated multiplication by 10, exact for 0 <= x <= 22 (10**22 is the
! largest power of ten a double holds exactly) and correctly rounded
! for -22 <= x < 0 (one division); exp(x*ln10) otherwise.
elemental function exp10(x) result(exp10_x)
      real(dp), intent(in) :: x
      real(dp) :: exp10_x
      integer :: ix, i

! the test on |x| comes first, so that int(x) cannot overflow
      if (abs(x) < 100.0d0 .and. x == aint(x)) then
         ix = int(x)
         exp10_x = 1.0d0
         do i = 1, abs(ix)
            exp10_x = exp10_x*10.0d0
         end do
         if (ix < 0) exp10_x = 1.0d0/exp10_x
      else
         exp10_x = exp(x*ln10)
      end if
end function exp10

#else

! ---------------------------------------------------------------
! Intrinsic backend: passthroughs, bit-identical to the expressions
! they replace at the call sites.
elemental function pow_r(x, y) result(pow_x)
      real(dp), intent(in) :: x, y
      real(dp) :: pow_x
      pow_x = x**y
end function pow_r

elemental function pow_i(x, iy) result(pow_x)
      real(dp), intent(in) :: x
      integer, intent(in) :: iy
      real(dp) :: pow_x
      pow_x = x**iy
end function pow_i

elemental function pow_r_sp(x, y) result(pow_x)
      real(dp), intent(in) :: x
      real, intent(in) :: y
      real(dp) :: pow_x
      pow_x = x**y
end function pow_r_sp

elemental function exp10(x) result(exp10_x)
      real(dp), intent(in) :: x
      real(dp) :: exp10_x
      exp10_x = 10.0d0**x
end function exp10

#endif

end module math_lib
