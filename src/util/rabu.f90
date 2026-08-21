!----------------------------------------------------------------------
! rabu
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original rabu.f; only variable names, source form, and comment style
! were updated.
!
! Reads one composition-table record (element count, atomic weights,
! number abundances, mass fractions, and mean molecular weight) from
! an unformatted file unit.
subroutine rabu(table_unit, nchem0, num_chem, atomic_weight, &
     number_abundance, mass_fraction, mean_molecular_weight)

      use const_lib
      implicit none

      integer, intent(in) :: table_unit, nchem0
      integer, intent(out) :: num_chem
      double precision, intent(out) :: atomic_weight(nchem0), &
           number_abundance(nchem0), mass_fraction(nchem0)
      double precision, intent(out) :: mean_molecular_weight


      save

      integer :: ic

!     NCHEM,ATWT,ABUN,ABFRCS ARE OUTPUT
!     READ(IR,99) NCHEM,(ATWT(IC),ABUN(IC),ABFRCS(IC),
!    1       IC=1,NCHEM),GASMU
      read(table_unit   ) num_chem,(atomic_weight(ic),number_abundance(ic), &
           mass_fraction(ic), ic=1,num_chem),mean_molecular_weight
      if (nchem0.lt.num_chem) then
         stop
      end if
      return
!  99   FORMAT(1X,I5,(/1X,3E15.7))
! 9009  FORMAT(' ERROR IN RABU. NCHEM READ FROM TABLE IS',
!      1 ' BIGGER THAN THE VALUE USED IN THE COMMONS.',
!      2 ' NCHEM0,NCHEM= ',/1X,2I8)
end subroutine rabu
