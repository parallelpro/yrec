!----------------------------------------------------------------------
! neutrino
!----------------------------------------------------------------------
! Modernized (free-form) 2026 as part of the YREC readability refactor.
! Logic and numerics are unchanged from the original neutrino.f; only
! the source form, explicit typing, and header/comment style were
! updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Added by Grant Newsham 9/06.
!
! Sets up the (currently 4-species: H, He, C, O) composition arrays
! expected by nulosses/sneut and calls through to them. Returns the
! neutrino energy-loss rate (snu) and its derivatives with respect to
! temperature and density (dsnudt, dsnudd).
subroutine neutrino(temp,den,x,y,zc,zo,snu,dsnudt,dsnudd)
      implicit none

      double precision, intent(in) :: temp, den, x, y, zc, zo
      double precision, intent(out) :: snu, dsnudt, dsnudd

      integer, parameter :: ionmax = 4
      double precision :: xmass(ionmax), ymass(ionmax), aion(ionmax), &
           zion(ionmax)

!..set the mass fractions, z's and a's of the composition
!..hydrogen

      aion(1)  = 1.0d0
      zion(1)  = 1.0d0

!..helium

      aion(2)  = 4.0d0
      zion(2)  = 2.0d0

!..carbon 12

      aion(3)  = 12.0d0
      zion(3)  = 6.0d0

!..oxygen 16

      aion(4)  = 16.0d0
      zion(4)  = 8.0d0

      xmass(1) = x
      xmass(2) = y
      xmass(3) = zc
      xmass(4) = zo

      call nulosses(temp,den,snu,xmass,ymass,aion,zion,dsnudt,dsnudd)

      return
end subroutine neutrino
