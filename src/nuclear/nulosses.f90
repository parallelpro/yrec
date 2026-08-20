!----------------------------------------------------------------------
! nulosses
!----------------------------------------------------------------------
! Modernized (free-form) 2026 as part of the YREC readability refactor.
! Logic and numerics are unchanged from the original nulosses.f; only
! the source form and header/comment style were updated -- this file's
! variable names were already lowercase and descriptive and were left
! as-is. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Neutrino loss subroutine modified by G. Newsham from codes of
! N. Itoh & F. Timmes (9/06).
!
! Currently uses just the carbon coefficients for the bremsstrahlung
! neutrino losses (the differences for other elements is slight.)
! We also currently set the fermi temperature cuttoff at 0.01 Tfermi.
! This is applicable for helium but will need to be amended for
! carbon and oxygen white dwarfs and cores of AGB stars.

subroutine nulosses(temp,den,snu,xmass,ymass,aion,zion, &
                     dsnudt,dsnudd)

      implicit none
      save

!..tests the neutrino loss rate routine

!..ionmax  = number of isotopes in the network
!..xmass   = mass fractions
!..ymass   = molar fractions
!..aion    = number of nucleons
!..zion    = number of protons

      integer          ionmax
      parameter        (ionmax=4)
      double precision xmass(ionmax),ymass(ionmax), &
                       aion(ionmax),zion(ionmax),temp,den,abar,zbar, &
                       snu,dsnudt,dsnudd,dsnuda,dsnudz



!..get abar and zbar
      call azbar(xmass,aion,zion,ionmax, &
                 ymass,abar,zbar)



!..get the neutrino losses
      call sneut(temp,den,abar,zbar, &
                  snu,dsnudt,dsnudd,dsnuda,dsnudz)



      return
end subroutine nulosses


subroutine azbar(xmass,aion,zion,ionmax, &
                 ymass,abar,zbar)
      implicit none
      save

!..input
!..mass fractions     = xmass(1:ionmax)
!..number of nucleons = aion(1:ionmax)
!..charge of nucleus  = zion(1:ionmax)
!..number of isotopes = ionmax

!..output:
!..molar abundances        = ymass(1:ionmax),
!..mean number of nucleons = abar
!..mean nucleon charge     = zbar



!..declare
      integer          i,ionmax
      double precision xmass(ionmax),aion(ionmax),zion(ionmax), &
                       ymass(ionmax),abar,zbar,zbarxx,ytot1


      zbarxx  = 0.0d0
      ytot1   = 0.0d0
      do i=1,ionmax
       ymass(i) = xmass(i)/aion(i)
       ytot1    = ytot1 + ymass(i)
       zbarxx   = zbarxx + zion(i) * ymass(i)
      enddo
      abar   = 1.0d0/ytot1
      zbar   = zbarxx * abar
      return
end subroutine azbar
