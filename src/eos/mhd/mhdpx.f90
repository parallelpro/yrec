!-------------------------    GROUP: SR_PX    -------------------------------
!
!     MHD INTERPOLATION PACKAGE
!     PRESSURE AS ARGUMENT INTERPOLATION IN X, Z FIXED
!
!     INPUT   PGL = LOG10 (GAS PRESSURE [DYN/CM2])
!             TL  = LOG10 (TEMPERATURE [K])
!             XC  = HYDROGEN ABUNDANCE (BY MASS)
!
!     OUTPUT  RL  = LOG10 (DENSITY [G/CM3])        : ARGUMENT
!             OTHER (SEE SEPARATE INSTRUCTIONS)    : COMMON/MHDOUT/...
!
!     ERROR   IERR = 1 SIGNALS PGL,TL OUTSIDE THE DOMAIN OF TABLES
!             (OTHERWISE IERR = 0).
!
!----------------------------------------------------------------------
! mhdpx
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mhdpx.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model). common/luout/, common/mhdout/,
! and common/ccout2/ member names match those established in
! meqos.f90.
!
! Top-level entry point to the MHD (Mihalas, Hummer, Dappen)
! equation-of-state table interpolator. S/R MHDST (via mhdtbl) must
! be called in main to load the tables first. Dispatches to the
! variable-X interpolation in mhdpx1, then returns log10(density).
subroutine mhdpx(log10_pressure, log10_temperature, hydrogen_fraction, &
     log10_density)
      use mhd_eos_lib
      use const_lib
      use luout_lib
      implicit none
      integer, parameter :: ivarx = 25
      integer, parameter :: nchem0 = 6

      double precision, intent(in) :: log10_pressure, log10_temperature, &
           hydrogen_fraction
      double precision, intent(out) :: log10_density

      save

      call mhdpx1(log10_pressure, log10_temperature, hydrogen_fraction)
!     IERR = 0
      log10_density = mhd_eos%mhd_output(1)
      return
!   999 CONTINUE
      write(iowr,*) 'ERROR (MHD): OUT OF TABLE RANGE. RETURN'
      write(short_file_unit,*) 'ERROR (MHD): OUT OF TABLE RANGE. RETURN'
      stop
end subroutine mhdpx
