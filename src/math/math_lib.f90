!----------------------------------------------------------------------
! math_lib
!----------------------------------------------------------------------
! The one place YREC gets its elementary math (2026, reproducibility
! campaign; modeled on MESA's math_lib). Two backends, chosen at
! build time:
!
!   make USE_CRMATH=1   correctly-rounded backend: crmath (Townsend's
!                       Fortran wrapper over crlibm, shipped in the
!                       MESA SDK). Every elementary function returns
!                       the IEEE-nearest result, so results are
!                       bit-identical across platforms/architectures
!                       when combined with -ffp-contract=off (set
!                       unconditionally in the Makefile) and the same
!                       gfortran major version.
!   make                intrinsic backend: pow/exp10 reduce to
!                       EXACTLY ** and the compiler's functions -- the
!                       historical behavior, byte-for-byte, so the
!                       call-site migration is verifiable against the
!                       existing pins. Fast, platform-dependent in the
!                       last bits.
!
! Usage contract (enforced by tools/check_boundaries.py):
!   * every source file that calls an elementary transcendental
!     (exp/log/log10/sin/cos/...) has `use math_lib` -- under the
!     crmath backend the module procedures SHADOW the intrinsics, so
!     call sites are unchanged;
!   * no real-exponent ** anywhere: x**y with real y is a hidden libm
!     pow -- write pow(x,y); 10**x is written exp10(x). Integer
!     exponents (x**2, x**n with integer n) are exact repeated
!     multiplication and stay as **.
!
! Under crmath, pow and exp10 follow MESA's definitions exactly: an
! integer-valued exponent takes the repeated-multiplication fast path
! (bit-exact), otherwise pow(x,y) = exp(log(x)*y) and exp10(x) =
! exp(x*ln10) -- fully determined by crmath's exp/log, hence
! reproducible. Under the intrinsic backend they are PURE passthroughs
! (x**y, 10**x) with no fast path, so converting a call site is a
! no-op bit-for-bit.
module math_lib
#ifdef YREC_CRMATH
      use crmath
#endif
      implicit none
      private

      public :: pow, exp10

! pow is generic over real and integer exponents: pow(x, iy) with
! integer iy exists so ** conversions never need to know the
! exponent's type (integer-typed sites route through the exact
! integer path; under the intrinsic backend it is x**iy verbatim).
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
      double precision, parameter :: ln10 = 2.302585092994045684d0
#else
      character(len=*), parameter, public :: math_backend = 'intrinsic'
#endif

contains

#ifdef YREC_CRMATH

! ---------------------------------------------------------------
! x**y for real y (MESA's definition). Integer-valued |y| < 100 uses
! exact repeated multiplication; otherwise the exp(log(x)*y)
! composition, fully determined by crmath. x = 0 returns 0 (the
! composition would take log(0)).
elemental function pow_r(x, y) result(pow_x)
      double precision, intent(in) :: x, y
      double precision :: pow_x
      integer :: iy, i

      if (x == 0.0d0) then
         pow_x = 0.0d0
         return
      end if
      iy = floor(y)
      if (y == iy .and. abs(iy) < 100) then
         pow_x = 1.0d0
         do i = 1, abs(iy)
            pow_x = pow_x*x
         end do
         if (iy < 0) pow_x = 1.0d0/pow_x
      else
         pow_x = exp(log(x)*y)
      end if
end function pow_r

! integer exponent: exact repeated multiplication (MESA's pow_i)
elemental function pow_i(x, iy) result(pow_x)
      double precision, intent(in) :: x
      integer, intent(in) :: iy
      double precision :: pow_x
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
! mixed-mode ** promotes the exponent to dp, so pow_r(x, dble(y))
! matches x**y exactly
elemental function pow_r_sp(x, y) result(pow_x)
      double precision, intent(in) :: x
      real, intent(in) :: y
      double precision :: pow_x
      pow_x = pow_r(x, dble(y))
end function pow_r_sp

! ---------------------------------------------------------------
! 10**x for real x (MESA's definition): exact for integer x,
! exp(x*ln10) otherwise.
elemental function exp10(x) result(exp10_x)
      double precision, intent(in) :: x
      double precision :: exp10_x
      integer :: ix, i

      ix = floor(x)
      if (x == ix .and. abs(ix) < 100) then
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
! Intrinsic backend: pure passthroughs, bit-identical to the
! expressions they replace at the call sites.
elemental function pow_r(x, y) result(pow_x)
      double precision, intent(in) :: x, y
      double precision :: pow_x
      pow_x = x**y
end function pow_r

elemental function pow_i(x, iy) result(pow_x)
      double precision, intent(in) :: x
      integer, intent(in) :: iy
      double precision :: pow_x
      pow_x = x**iy
end function pow_i

elemental function pow_r_sp(x, y) result(pow_x)
      double precision, intent(in) :: x
      real, intent(in) :: y
      double precision :: pow_x
      pow_x = x**y
end function pow_r_sp

elemental function exp10(x) result(exp10_x)
      double precision, intent(in) :: x
      double precision :: exp10_x
      exp10_x = 10.0d0**x
end function exp10

#endif

end module math_lib
