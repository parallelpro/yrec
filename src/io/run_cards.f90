!----------------------------------------------------------------------
! run_cards
!----------------------------------------------------------------------
! The run plan: the kind cards (evolve, rescale, rescale and evolve) read
! from the inlist turned into the per-card stop conditions, rescale
! targets and starting-model flags, echoed to the terminal and the run
! log (formerly part of parmin).
module run_cards
      use io_units
      use eos_lib, only: eos_set_mixture
      use legacy_controls
      implicit none
      private
      public :: interpret_kind_cards

contains

! ---------------------------------------------------------------
! Interpret the per-run "kind" cards: stop conditions, rescale
! parameters, and the rescaled envelope mixture for each run.
subroutine interpret_kind_cards(s)
      use star_info_lib, only: star_info
      type(star_info), intent(inout) :: s
      integer :: i, nkind

      write(run_log_unit,200)
      200 format(/1x,'run plan',/)

      lfirst(1) = .true.

!     RUN LOOP
      do nkind=1, numrun
! READ IN NMODLS AND MODEL SOURCE(MEMORY OR FIRST MODEL)-SAME FOR ALL
       rescale_kind(nkind) = kindrn(nkind)
! The cards that evolve (kinds 1 and 3) stop at their age and central
! abundance stops, and may fix the timestep.
       if (kindrn(nkind) == 1 .or. kindrn(nkind) == 3) then
          end_age_stop_active(nkind) = endage(nkind) > 0.0d0 .or. &
               end_dcen(nkind) > 0.0d0 .or. end_xcen(nkind) > 0.0d0 .or. &
               end_ycen(nkind) > 0.0d0
          timestep_override_active(nkind) = setdt(nkind) > 0.0d0
       end if
       if(kindrn(nkind).eq.1) then
! EVOLVE CARD
            if (nmodls(nkind).gt.0) then
          if (lfirst(nkind)) then
             write(terminal_unit,350) nkind,nmodls(nkind)
             write(run_log_unit,350) nkind,nmodls(nkind)
      350          format(1x,'card',i3,': evolve up to',i6, &
           &          ' models from the starting model')
          else
             write(terminal_unit,351) nkind, nmodls(nkind)
             write(run_log_unit,351) nkind, nmodls(nkind)
      351          format(1x,'card',i3,': evolve up to',i6, &
           &          ' models from the previous card''s model')
          end if
! GENERALIZE STOP CONDITIONS
          if(end_age_stop_active(nkind).or.timestep_override_active(nkind)) then
             write(terminal_unit,370) endage(nkind), setdt(nkind), &
           &          end_dcen(nkind), end_xcen(nkind), end_ycen(nkind)
             write(run_log_unit,370) endage(nkind), setdt(nkind), &
           &          end_dcen(nkind), end_xcen(nkind), end_ycen(nkind)
      370          format(9x,'stop conditions (0 = unused): age =',es9.2, &
           &          ' yr   fixed dt =',es9.2,' yr   central: D =',es10.3, &
           &          '  X =',es10.3,'  Y =',es10.3)
          endif
            end if
       else if(kindrn(nkind).eq.2) then
! RESCALE CARD:  RESCALE STARTING MODEL
! QUANTITIES TO BE RESCALED STORED IN ARRAY RESCALE(4,50)
! WHERE THE ELEMENTS MASS,X,Z,CORE MASS ARE STORED IN ORDER
          rescale_params(1,nkind) = rsclm(nkind)
          rescale_params(2,nkind) = rsclx(nkind)
          rescale_params(3,nkind) = rsclz(nkind)
          rescale_params(4,nkind) = rsclcm(nkind)
            if (nmodls(nkind) .gt. 0) then
          if (lfirst(nkind)) then
             write(terminal_unit,450) nkind
             write(run_log_unit,450) nkind
      450          format(1x,'card',i3,': rescale the starting model')
          else
             write(terminal_unit,451) nkind
             write(run_log_unit,451) nkind
      451          format(1x,'card',i3, &
           &          ': rescale the previous card''s model')
          end if
          write(terminal_unit,452) nmodls(nkind),(rescale_params(i,nkind),i = 1,4)
          write(run_log_unit,452) nmodls(nkind), &
           &       (rescale_params(i,nkind),i = 1,4)
      452       format(9x,'relax for',i3,' models; targets', &
           &       ' (0 or negative = keep current):', &
           &       '  M =',f9.6,'  X =',f9.6,'  Z =',f9.6, &
           &       '  core mass =',f10.6)
            end if
         else if(kindrn(nkind).eq.3) then
! RESCALE AND EVOLVE CARD:  RESCALE STARTING MODEL
! QUANTITIES TO BE RESCALED STORED IN ARRAY RESCALE(4,50)
! WHERE THE ELEMENTS MASS,X,Z,CORE MASS ARE STORED IN ORDER
            rescale_params(1,nkind) = rsclm(nkind)
            rescale_params(2,nkind) = rsclx(nkind)
            rescale_params(3,nkind) = rsclz(nkind)
            rescale_params(4,nkind) = rsclcm(nkind)
            if (lfirst(nkind)) then
               write(terminal_unit,550) nkind
               write(run_log_unit,550) nkind
      550          format(1x,'card',i3, &
           &          ': rescale and evolve the starting model')
            else
               write(terminal_unit,551) nkind
               write(run_log_unit,551) nkind
      551          format(1x,'card',i3, &
           &          ': rescale and evolve the previous card''s model')
            end if
            write(terminal_unit,452) nmodls(nkind),(rescale_params(i,nkind),i = 1,4)
            write(run_log_unit,452) nmodls(nkind), &
           &       (rescale_params(i,nkind),i = 1,4)
            if (end_age_stop_active(nkind) .or. timestep_override_active(nkind)) then
               write(terminal_unit,370) endage(nkind), setdt(nkind), &
                    end_dcen(nkind), end_xcen(nkind), end_ycen(nkind)
               write(run_log_unit,370) endage(nkind), setdt(nkind), &
                    end_dcen(nkind), end_xcen(nkind), end_ycen(nkind)
            end if
       end if
! a positive rescale target sets the envelope Z (0, the default, keeps it)
         if(rescale_params(3,nkind).gt.0.0d0)  s%envelope_metal_fraction=rescale_params(3,nkind)
! keep the eos-side mixture in step with the rescaled Z
         call eos_set_mixture(s%envelope_hydrogen_fraction, &
              s%envelope_metal_fraction, s%amuenv, &
              s%fxenv)
      end do
end subroutine interpret_kind_cards

end module run_cards
