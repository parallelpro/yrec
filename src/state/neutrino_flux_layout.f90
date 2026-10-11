!----------------------------------------------------------------------
! neutrino_flux_layout
!----------------------------------------------------------------------
! The slots of the solar-neutrino flux vectors: net/'s
! nuclear_energy_result%neutrino_flux and the model's
! neutrino_flux_total(n_nu_fluxes) and neutrino_flux_zone(n_nu_fluxes,json)
! (the SNU tables). Slots 1-8 are the real neutrino sources; slots 9-10
! hold the He3+He3 and He3+He4 reaction rates in the same flux units, as
! if each reaction emitted one neutrino ("fictional" fluxes, tabulated
! with the others by neutrino_flux_table). Constants only, shared like
! composition_layout.
module neutrino_flux_layout
      implicit none
      private
      integer, parameter, public :: &
           i_nu_pp = 1, i_nu_pep = 2, i_nu_hep = 3, i_nu_be7 = 4, &
           i_nu_b8 = 5, i_nu_n13 = 6, i_nu_o15 = 7, i_nu_f17 = 8, &
           i_nu_he3he3 = 9, i_nu_he3he4 = 10, n_nu_fluxes = 10
end module neutrino_flux_layout
