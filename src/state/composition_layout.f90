!----------------------------------------------------------------------
! composition_layout
!----------------------------------------------------------------------
! The rows of the composition vector xa (one column per zone): the
! species each slot holds, and how many slots the two composition sets
! use. Constants only, shared by the model state (star_info_lib
! re-exports them) and the physics domains, which take the layout from
! here rather than from star_info_lib.
!
! Slot 3 is the TOTAL metal mass fraction Z (YREC evolves Z as one
! species); the heavy elements in slots 5..11 and 13..15 are diffusion
! and burning tracers carried inside Z, not independent contributions
! to it. The basic set is slots i_h1..i_o18; the extended set
! (job%use_extended_composition) adds i_h2..i_be9. Composition arrays
! are always dimensioned for the extended set.
module composition_layout
      implicit none
      private
      integer, parameter, public :: &
           i_h1  = 1,  i_he4 = 2,  i_metals = 3,  i_he3 = 4, &
           i_c12 = 5,  i_c13 = 6,  i_n14 = 7,     i_n15 = 8, &
           i_o16 = 9,  i_o17 = 10, i_o18 = 11,    i_h2  = 12, &
           i_li6 = 13, i_li7 = 14, i_be9 = 15
      integer, parameter, public :: &
           n_species_basic = 11, n_species_extended = 15
end module composition_layout
