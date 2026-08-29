!----------------------------------------------------------------------
! write_mod_model
!----------------------------------------------------------------------
! Added 2026 (retire-legacy campaign): the .mod model file, replacing
! the .last MODEL2 format. Contains exactly what a reload actually
! uses (established by audit of read_starting_model): the global
! scalars, the per-shell structure/composition/omega table, and the
! two deck-consistency values (mixing length, extended-composition
! flag). Everything the old formats stored beyond that -- the
! envelope-triangle ENV records, fit coefficients, trial values,
! TLUMX luminosity breakdown, CZ indices, and the config echo line --
! was read and then discarded or unconditionally recomputed by the
! loader ("ENVELOPE DATA (Now bypassed)"), so it is not written.
!
! MESA-style named layout: a version tag, named global values (one
! per line; readers ignore unknown names), a named column-header
! line, then one row per shell at full double precision. The old
! MODEL2 format truncated the species to ~9 significant digits;
! this format round-trips doubles exactly, which matters for
! multi-card runs that reload the model between cards.
!
! Everything is read from star%; the unit is the only argument
! (rewound first: the .mod file always holds ONE model, the latest).
subroutine write_mod_model(iwrite)

      use star_info_lib, only: star, i_h1, i_he4, i_metals, i_he3, &
           i_c12, i_c13, i_n14, i_n15, i_o16, i_o17, i_o18, i_h2, &
           i_li6, i_li7, i_be9
      use luout_lib
      use run_log_lib, only: solver_diagnostics
      implicit none

      integer, intent(in) :: iwrite

      integer :: k
      integer, parameter :: slots(15) = [i_h1, i_he4, i_metals, i_he3, &
           i_c12, i_c13, i_n14, i_n15, i_o16, i_o17, i_o18, i_h2, &
           i_li6, i_li7, i_be9]

      rewind iwrite
      write(iwrite,'(a)') 'YMOD 1'
      write(iwrite,'(a,1x,i12)')        'model_number        ', star%model_number
      write(iwrite,'(a,1x,i12)')        'num_zones           ', star%nz
      write(iwrite,'(a,1x,es26.16e3)')  'star_mass_msun      ', star%star_mass
      write(iwrite,'(a,1x,es26.16e3)')  'log_total_mass      ', star%log_total_mass
      write(iwrite,'(a,1x,es26.16e3)')  'star_age_gyr        ', star%dage
      write(iwrite,'(a,1x,es26.16e3)')  'timestep_yr         ', star%timestep_yr
      write(iwrite,'(a,1x,es26.16e3)')  'log_L               ', star%log_L
      write(iwrite,'(a,1x,es26.16e3)')  'log_Teff            ', star%log_Teff
      write(iwrite,'(a,1x,es26.16e3)')  'mixing_length_alpha ', star%mixing_length_alpha
      write(iwrite,'(a,1x,l1)')         'extended_composition', star%job%use_extended_composition
      write(iwrite,'(a,1x,l1)')         'rotation_active     ', star%job%rotation_active
      write(iwrite,'(a)') 'columns: zone log_mass logR luminosity logP logT' &
           // ' logRho conv omega h1 he4 z he3 c12 c13 n14 n15 o16 o17 o18' &
           // ' h2 li6 li7 be9'
      do k = 1, star%nz
         write(iwrite,'(i6,6(1x,es26.16e3),1x,l1,16(1x,es26.16e3))') k, &
              star%log_mass(k), star%logR(k), star%luminosity_lsun(k), &
              star%logP(k), star%logT(k), star%logRho(k), &
              star%convective_flag(k), star%omega(k), &
              star%xa(slots(1),k),  star%xa(slots(2),k),  star%xa(slots(3),k), &
              star%xa(slots(4),k),  star%xa(slots(5),k),  star%xa(slots(6),k), &
              star%xa(slots(7),k),  star%xa(slots(8),k),  star%xa(slots(9),k), &
              star%xa(slots(10),k), star%xa(slots(11),k), star%xa(slots(12),k), &
              star%xa(slots(13),k), star%xa(slots(14),k), star%xa(slots(15),k)
      end do
      flush(iwrite)

! 2026 log redesign: the per-model terminal echo is replaced by the
! run-log progress line (run_log_lib); the "DUMPED MODEL" bookkeeping
! is solver forensics behind the diagnostics flag.
      if (solver_diagnostics()) then
         if (iwrite .eq. last_model_unit) then
          write(run_log_unit,330) star%model_number, iwrite
         else
          write(run_log_unit,340) star%model_number, star%dage, iwrite
         endif
      end if
  330 format(' DUMPED MODEL',I5,'  FILE',I3)
  340 format(' DUMPED MODEL',I5,' AGE',F13.9,'  FILE',I3)

      return
end subroutine write_mod_model
