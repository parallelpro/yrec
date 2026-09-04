!----------------------------------------------------------------------
! eos_mixture_lib
!----------------------------------------------------------------------
! Added 2026 (physics-purity pass -- see ROADMAP.md "Decoupling the
! physics domains from the model"). The envelope mixture the eos
! domain evaluates with: surface hydrogen/metal fractions, the mean
! atomic mass, and the 12-species number fractions consumed by the
! Yale (eqrelv/eqsaha) and SCV branches of eqstat2 and by eos/mu.f90.
!
! Previously the eos internals read this straight from the model
! (star%...), the one place the eos domain touched
! star_info. Now the STAR LAYER pushes the mixture through the
! eos_set_mixture facade entry (eos_lib) whenever it recomputes the
! envelope composition (core/read_starting_model.f90's mixture blocks,
! setup/rescale_model.f90's rescale), and the eos domain reads only this,
! its own state -- same pattern as every other physics-domain state
! module. MESA-equivalently: composition is an INPUT of the eos,
! owned on the physics side of the boundary.
module eos_mixture_lib
      implicit none
      private

      type, public :: eos_mixture_state
            double precision :: envelope_hydrogen_fraction, &
                 envelope_metal_fraction
            double precision :: amuenv
! fxenv(i): number fraction of species i in the envelope mixture,
! in the order of eqstat2's atomic_weights_full table:
! 1 Na, 2 Al, 3 Mg, 4 Fe, 5 Si, 6 C, 7 H, 8 O, 9 N, 10 Ar, 11 Ne, 12 He
! (slots 7 and 12 are set from X and Y directly, the metals from the
! mixture_weights_seed control scaled to Z).
            double precision :: fxenv(12)
      end type eos_mixture_state

      type(eos_mixture_state), public, save :: eos_mix

end module eos_mixture_lib
