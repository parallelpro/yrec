!----------------------------------------------------------------------
! atm_lib
!----------------------------------------------------------------------
! Added 2026 as part of the YREC readability refactor's phase two
! (disentangling the solver from the physics domains -- see
! GUIDELINES.md's "Physics domains still entangled with the solver").
! Like kap_lib's kap_get, this is not a new dispatch consolidation:
! envint (renamed atm_get) was already the domain's generic "solve one
! envelope for this (Teff, L)" primitive, called uniformly (no
! duplicated dispatch) from 7 sites -- 6 external, plus atm/surfbc.f90
! (the domain's other, more specialized entry point: a single-caller,
! solver-loop-internal wrapper around atm_get that adds a cached
! (Teff, L) triangle for cheap derivative interpolation and the
! hot-edge gray-atmosphere fallback -- surfbc.f90 itself is unchanged
! by this rename). This rename/module-wrap is purely to give atm/ the
! same public-facade shape as eos_lib/kap_lib.
!
! Computes one envelope solution (P, T, R at the model's fitting
! point) for a given (Teff, L) vertex: first integrates the gray (or
! Eddington/Krishna-Swamy) atmosphere in optical depth from a starting
! guess down to tau=2/3 via QATM (bypassed entirely for tabulated
! Kurucz/Allard atmospheres, which look up the boundary pressure
! directly via SURFP/KCSURFP/ALSURFP, in atm/tables/), then integrates
! the envelope structure equations in pressure from tau=2/3 to the
! fitting mass point via QENV. Both integrations use the Bulirsch-
! Stoer stepper BSSTEP. Also saves the full atmosphere/envelope
! structure (for profile output and pulsation) and, unless the caller
! has switched to the newer TAUINTNEW-based method, locates the
! surface convection zone within the envelope via TAUCAL (in
! atm/turnover/).
module atm_lib
      implicit none
contains

!----------------------------------------------------------------------
! atm_init
!----------------------------------------------------------------------
! Added 2026 (phase three, ROADMAP.md stage 1): the atm domain's
! startup-time table-load lifecycle entry, following MESA's
! <mod>_init convention. Covers, per the atm_choice option in force:
! the Kurucz surface-pressure table (atm_choice 3 and 4, into
! const_lib's kurucz_* tables plus atm_table's gmin/gmax edge
! indices), the Allard NextGen atmosphere tables (atm_choice 4, via
! tables/alfilein.f90), and the Kurucz/Castelli surface-pressure
! table (atm_choice 5). All of this previously lived as an inline
! open/read block in setup/setups.f90 -- moved here verbatim,
! including the G Somers 5/15 invalid-pressure-edge catches and the
! preserved use of `ng` (not `ngc`) in the Castelli branch's final
! failsafe (harmless, both are 11; not "fixed").
subroutine atm_init(atm_table_path, allard_table_path, ierr)

      use atm_table_lib
      use const_lib
      implicit none
! PARAMETERS NT AND NG FOR TABULATED SURFACE PRESSURES OF KURUCZ,
! NTC/NGC FOR KURUCZ/CASTELLI (matching setups.f90's historical
! declarations).
      integer, parameter :: nt = 57, ng = 11
      integer, parameter :: ntc = 76, ngc = 11

      character(len=256), intent(in) :: atm_table_path
      character(len=256), intent(in) :: allard_table_path
! 2026 (ROADMAP.md stage 3): OPTIONAL ierr, same contract as atm_get's.
      integer, intent(out), optional :: ierr

      integer :: jerr
      integer :: teff_idx, logg_idx
      logical :: found_valid_pressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INPUT PRESSURE TABLE FOR SURFACE BOUNDARY CONDITIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (present(ierr)) ierr = 0
      jerr = 0

      if ((atm_choice .eq. 3) .or. (atm_choice .eq. 4)) then
! OPEN SURFACE PRESSURE TABLE
        open(atm_table_file_unit,file=atm_table_path, status='OLD')
! GET ABUNDANCE:
        read(atm_table_file_unit,200) kurucz_table_z
  200   format(1x,10x,1pe16.8,/)
! GET VALUES OF LOG Teff:
        read(atm_table_file_unit,201) (kurucz_teff_table(teff_idx),teff_idx=1,nt)
  201   format(1p4e16.8,13(/1p4e16.8),/1pe16.8,/)
! GET VALUES OF LOG G:
        read(atm_table_file_unit,202) (kurucz_logg_table(logg_idx),logg_idx=1,ng)
  202   format(1p4e16.8,/1p4e16.8,/1p3e16.8,/)
! GET GRID OF LOG PRESSURE VALUES:
        do teff_idx=1,nt
          read(atm_table_file_unit,203) (kurucz_log10_pressure_table(teff_idx,logg_idx),logg_idx=1,ng)
  203     format(1p4e16.8,/1p4e16.8,/1p3e16.8)
  205   continue
        end do
        rewind atm_table_file_unit
        close(atm_table_file_unit)
!       G Somers 5/15; added a catch that allows the code to work
!                      if the highest gravity P value for a specific
!                      temperature is -999. LCATCH is set to true
!                      once the first valid P is read, and it can
!                      set the minimum P thereafter.
         do teff_idx = 1,nt
            found_valid_pressure = .false.
            do logg_idx = ng,1,-1
               if(kurucz_log10_pressure_table(teff_idx,logg_idx).le.0.0d0)then
!                 check if the first non-999 value has been reached.
!                 if so, set IMIN.
                  if(found_valid_pressure)then
                     atm_table%kurucz_gmin_index(teff_idx) = logg_idx + 1
                     exit
                  endif
               else
!                 if we have reached a non-negative pressure value,
!                 turn off the catch so IMIN can be set. also record
!                 the highest gravity with a pressure for later int.
                  if(.not.found_valid_pressure) atm_table%kurucz_gmax_index(teff_idx) = logg_idx
                  found_valid_pressure = .true.
               endif
            end do
            if (logg_idx .lt. 1) atm_table%kurucz_gmin_index(teff_idx) = 1
!           if all of the P values at a given T are -999, set IMIN
!           to the number of gravity terms. in responce, the code
!           should break when trying to find surface P.
            if(.not.found_valid_pressure) atm_table%kurucz_gmin_index(teff_idx) = ng
         end do
!        G Somers 5/15 END
! MHP 6/97 ADDED OPTION FOR ALLARD MODEL ATMOSPHERES; USED INSTEAD OF
! KURUCZ FOR TEFF < 10,000 K.
         if(atm_choice .eq. 4)then
!            ATMZA = 0.02D0
          call alfilein(allard_table_path, jerr)      ! Get Allard Atmospheres files and
          if (jerr /= 0) then
             if (present(ierr)) then
                ierr = jerr
                return
             end if
             stop
          end if
         endif                  ! initialize tables. 9/23/08 LLP


! JNT 6/14 ADD OPTION FOR NEW KURUCZ/CASTELLI ATMOSHPERES
      else if (atm_choice .eq. 5) then
! OPEN SURFACE PRESSURE TABLE
        open(atm_table_file_unit,file=atm_table_path, status='OLD')
! GET ABUNDANCE:
        read(atm_table_file_unit,200) kurucz_table_z
! GET VALUES OF LOG Teff:
        read(atm_table_file_unit,206) (kurucz_castelli_teff_table(teff_idx),teff_idx=1,ntc)
  206   format(1p4e16.8,18(/1p4e16.8),/)
! GET VALUES OF LOG G:
        read(atm_table_file_unit,208) (kurucz_castelli_logg_table(logg_idx),logg_idx=1,ngc)
  208   format(1p4e16.8,/1p4e16.8,/1p3e16.8,/)
! GET GRID OF LOG PRESSURE VALUES:
        do teff_idx=1,ntc
          read(atm_table_file_unit,203) (kurucz_castelli_log10_pressure_table(teff_idx,logg_idx),logg_idx=1,ngc)
  207   continue
        end do
        rewind atm_table_file_unit
        close(atm_table_file_unit)
!       G Somers 5/15; added a catch that allows the code to work
!                      if the highest gravity P value for a specific
!                      temperature is -999. LCATCH is set to true
!                      once the first valid P is read, and it can
!                      set the minimum P thereafter.
        do teff_idx = 1,ntc
           found_valid_pressure = .false.
           do logg_idx = ngc,1,-1
              if(kurucz_castelli_log10_pressure_table(teff_idx,logg_idx).le.0.0d0)then
!                check if the first non-999 value has been reached.
!                if so, set IMIN2.
                 if(found_valid_pressure)then
                    atm_table%castelli_gmin_index(teff_idx) = logg_idx + 1
                    exit
                 endif
              else
!                if we have reached a non-negative pressure value,
!                turn off the catch so IMIN2 can be set. also record
!                the highest gravity with a pressure for later int.
                 if(.not.found_valid_pressure) atm_table%castelli_gmax_index(teff_idx) = logg_idx
                 found_valid_pressure = .true.
              endif
           end do
           if (logg_idx .lt. 1) atm_table%castelli_gmin_index(teff_idx) = 1
!          if all of the P values at a given T are -999, set IMIN
!          to the number of gravity terms. in responce, the code
!          should break when trying to find surface P.
           if(.not.found_valid_pressure) atm_table%castelli_gmin_index(teff_idx) = ng
        end do
! END JNT 6/14
      endif

      return

! error funnel: same contract as atm_get's.
end subroutine atm_init

!----------------------------------------------------------------------
! atm_get_surface_pt
!----------------------------------------------------------------------
! Added 2026 (phase three, ROADMAP.md stage 1): public accessor for
! the Allard-atmosphere surface lookup, created so wind/massloss.f90
! (previously the last file calling tables/alsurfp.f90 directly) can
! go through the facade. Given log10(Teff) and log10(g), performs the
! Allard table interpolation; results land in atm_table state
! (atm_table%atm_log10_pressure/atm_log10_temperature), exactly as
! alsurfp has always delivered them. lookup_failed is set when the
! requested point falls outside the Allard tables (the caller then
! falls back per its own policy -- see alsurfp.f90's header).
subroutine atm_get_surface_pt(log_teff, log_g, print_to_files, &
     lookup_failed, ierr)

      implicit none

      double precision, intent(in) :: log_teff, log_g
      logical, intent(in) :: print_to_files
      logical, intent(out) :: lookup_failed
! 2026 (ROADMAP.md stage 3): OPTIONAL ierr, same contract as atm_get's.
      integer, intent(out), optional :: ierr

      integer :: jerr

      if (present(ierr)) ierr = 0
      jerr = 0

      call alsurfp(log_teff, log_g, print_to_files, lookup_failed, jerr)
      if (.not. (jerr /= 0)) then

      return

! error funnel: same contract as atm_get's.
      end if
end subroutine atm_get_surface_pt

end module atm_lib
