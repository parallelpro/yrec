!----------------------------------------------------------------------
! wrthead
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original wrthead.f; only variable names, source form, and comment
! style were updated.
!
! Write the headers for all the appropriate output files.
subroutine write_output_headers(total_mass_msun)

      use luout_lib
      use star_info_lib, only: star
      implicit none

      double precision, intent(in) :: total_mass_msun
! --- locals ---
      double precision :: total_mass_grams

      if (star%job%rescale_kind(star%job%nk) .eq. 1) then
         write(iowr, 47) star%job%nk, &
              star%mixing_length_alpha, star%job%num_models(star%job%nk)
      else if (star%job%rescale_kind(star%job%nk) .eq. 2) then
         write(iowr, 48) star%job%nk, star%job%initial_envelope_x, star%job%initial_envelope_z, &
              star%mixing_length_alpha, star%job%num_models(star%job%nk)
      else if (star%job%rescale_kind(star%job%nk) .eq. 3) then
         write(iowr, 49) star%job%nk, star%job%initial_envelope_x, star%job%initial_envelope_z, &
              star%mixing_length_alpha, star%job%num_models(star%job%nk)
      end if
  47  format(/,1x,'card',i3,' (evolve):  mixing length =',f8.6, &
             '  models =',i6)
  48  format(/,1x,'card',i3,' (rescale):         X =',f8.6, &
             '  Z =',f8.6,'  mixing length =',f8.6,'  models =',i6)
  49  format(/,1x,'card',i3,' (rescale+evolve):  X =',f8.6, &
             '  Z =',f8.6,'  mixing length =',f8.6,'  models =',i6)

      if (star%ctrl%isochrone_output_active) then
! header stuff for isochrone output
         total_mass_grams = total_mass_msun*star%solar_mass_cgs
         write(star%ctrl%isochrone_file_unit, 1495) total_mass_grams, &
              star%job%initial_envelope_x,star%job%initial_envelope_z,star%mixing_length_alpha,star%solar_bolometric_magnitude
 1495    format(7X, 1P5E16.8)
      end if

! 2026 retire-legacy: the .track header block is deleted with the
! .track writer (history.data carries the run metadata).

      return
end subroutine write_output_headers
