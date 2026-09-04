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
     number_abundance, mass_fraction, mean_molecular_weight, ierr)

      implicit none

      integer, intent(in) :: table_unit, nchem0
      integer, intent(out) :: num_chem
      double precision, intent(out) :: atomic_weight(nchem0), &
           number_abundance(nchem0), mass_fraction(nchem0)
      double precision, intent(out) :: mean_molecular_weight
      integer, intent(out) :: ierr

      integer :: ic

      ierr = 0

      read(table_unit   ) num_chem,(atomic_weight(ic),number_abundance(ic), &
           mass_fraction(ic), ic=1,num_chem),mean_molecular_weight
      if (nchem0.lt.num_chem) then
! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      return
end subroutine rabu
