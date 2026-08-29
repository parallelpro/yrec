!----------------------------------------------------------------------
! luout_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/luout/: I/O logical unit numbers,
! assigned once at run startup (core/read_input.f90) and read broadly --
! global configuration, not per-call data, so (like const_lib) this
! becomes a module of plain module-level variables rather than
! subroutine arguments. Lives alongside const_lib in this folder for
! the same reason: a codebase-wide foundational module, not belonging
! to any one physics domain.
!
! Canonical names below were chosen as the majority spelling already
! used across the converted files (short_file_unit/imilne/iowr); a
! handful of files used non-canonical local names for these three
! slots (ishort, milne_file_unit, main_output_unit) and had their
! body references renamed to the canonical names when converted.
module luout_lib
      implicit none

      integer :: ilast, idebug, short_file_unit, &
           iowr
! 2026 phase A (controls->star% campaign): the two opacity-table
! units from former common/zramp/, assigned by parmin like the rest
! of this module, moved here from controls_lib.
      integer :: iolaol2, ioopal2

end module luout_lib
