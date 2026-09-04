!----------------------------------------------------------------------
! gmass
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original gmass.f; only variable names, source form, and comment
! style were updated.
!
! Called from esac.f90 (OPAL 1995 EOS reader) to get the total mass
! per particle (in amu) for a mix of fully-ionized H, He, C, N, O,
! Ne at composition (hydrogen_fraction, metal_fraction), assuming a
! fixed solar-relative C:N:O:Ne mix within the metals. Also returns
! the ground-state ionization energy and mole/mass fractions of each
! species, used by radsub.f90/esac.f90 for the radiation and specific
! heat corrections.
double precision function gmass(hydrogen_fraction, metal_fraction, &
     total_moles, ground_state_energy, metal_mole_fraction, &
     species_mass_fraction)

      implicit none

      double precision, intent(in) :: hydrogen_fraction, metal_fraction
      double precision, intent(out) :: total_moles, ground_state_energy, &
           metal_mole_fraction
      double precision, intent(out) :: species_mass_fraction(7)

! atomic_number(i): atomic number (= electrons per fully-ionized atom)
! of species i, for i=1,6 = Ne, O, N, C, He, H.
      double precision :: atomic_number(6)
! ionization_energy(i): total ground-state ionization energy (eV) of
! species i, same ordering as atomic_number.
      double precision :: ionization_energy(6)
! atomic_weight(i): atomic weight (amu), i=1 electron, 2 Ne, 3 O, 4 N,
! 5 C, 6 He, 7 H.
      double precision :: atomic_weight(7)
! solar-relative mole-fraction-per-unit-Z of C, N, O, Ne in the metal
! mix (fixed composition assumed for all Z).
      double precision :: carbon_mix_fraction, nitrogen_mix_fraction, &
           oxygen_mix_fraction, neon_mix_fraction
      double precision :: carbon_moles, nitrogen_moles, oxygen_moles, &
           neon_moles
      double precision :: hydrogen_moles, helium_moles, total_moles_raw
! electron_mole_excess: electrons per mole of ions. total_moles sums
! (1 + Z_i) over the species mole fractions (which sum to 1), i.e.
! ions plus their electrons, so total_moles - 1 is the electron count.
      double precision :: electron_mole_excess
      integer :: species_idx

      data (ionization_energy(species_idx), species_idx=1,6) &
           /-3388.637d0, -1970.918d0, -1431.717d0, -993.2303d0, &
           -76.2315d0, -15.29409d0/
      data (atomic_number(species_idx), species_idx=1,6) &
           /10.0d0, 8.0d0, 7.0d0, 6.0d0, 2.0d0, 1.0d0/
      carbon_mix_fraction = 0.247137766d0
      nitrogen_mix_fraction = 0.0620782d0
      oxygen_mix_fraction = 0.52837118d0
      neon_mix_fraction = 0.1624188d0
      atomic_weight(7) = 1.0079d0
      atomic_weight(6) = 4.0026d0
      atomic_weight(5) = 12.011d0
      atomic_weight(4) = 14.0067d0
      atomic_weight(3) = 15.9994d0
      atomic_weight(2) = 20.179d0
      atomic_weight(1) = 0.00054858d0
      metal_mole_fraction = metal_fraction/(carbon_mix_fraction* &
           atomic_weight(5) + nitrogen_mix_fraction*atomic_weight(4) + &
           oxygen_mix_fraction*atomic_weight(3) + neon_mix_fraction* &
           atomic_weight(2))
      carbon_moles = metal_mole_fraction*carbon_mix_fraction
      nitrogen_moles = metal_mole_fraction*nitrogen_mix_fraction
      oxygen_moles = metal_mole_fraction*oxygen_mix_fraction
      neon_moles = metal_mole_fraction*neon_mix_fraction
      hydrogen_moles = hydrogen_fraction/atomic_weight(7)
      helium_moles = (1.0d0 - hydrogen_fraction - metal_fraction)/ &
           atomic_weight(6)
      total_moles_raw = hydrogen_moles + helium_moles + carbon_moles + &
           nitrogen_moles + oxygen_moles + neon_moles
      species_mass_fraction(6) = hydrogen_moles/total_moles_raw
      species_mass_fraction(5) = helium_moles/total_moles_raw
      species_mass_fraction(4) = carbon_moles/total_moles_raw
      species_mass_fraction(3) = nitrogen_moles/total_moles_raw
      species_mass_fraction(2) = oxygen_moles/total_moles_raw
      species_mass_fraction(1) = neon_moles/total_moles_raw
      ground_state_energy = 0.0d0
      total_moles = 0.0d0
      do species_idx = 1, 6
         ground_state_energy = ground_state_energy + &
              ionization_energy(species_idx)*species_mass_fraction(species_idx)
         total_moles = total_moles + (1.0d0 + atomic_number(species_idx))* &
              species_mass_fraction(species_idx)
      end do
      electron_mole_excess = total_moles - 1.0d0
      gmass = electron_mole_excess*atomic_weight(1)
      do species_idx = 2, 7
         gmass = gmass + atomic_weight(species_idx)* &
              species_mass_fraction(species_idx - 1)
      end do

      return
end function gmass
