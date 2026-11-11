!----------------------------------------------------------------------
! eos_source_list
!----------------------------------------------------------------------
! The table sources eos consults, in order (eos_table_source says what
! an entry holds and what its stage means): eos_init reads the selected
! ones' tables in this order, after the Fermi-Dirac table; blend_eos
! applies the selected ones of each of its two stages in this order;
! eos_get takes a selected whole-EOS source's result alone. Adding a
! source is one entry here (and its use).
module eos_source_list
      use eos_table_source, only: eos_table_source_t
      use scv_eos, only: scv_table_source
      use opal_eos, only: opal_table_source
      use mhd_eos, only: mhd_table_source
      implicit none
      private
      public :: table_sources

      integer, parameter, public :: n_table_sources = 3

contains

function table_sources() result(list)
      type(eos_table_source_t) :: list(n_table_sources)
      list = [scv_table_source(), opal_table_source(), mhd_table_source()]
end function table_sources

end module eos_source_list
