!----------------------------------------------------------------------
! atm_lib
!----------------------------------------------------------------------
! Added 2026 as part of the YREC readability refactor's phase two
! (disentangling the solver from the physics domains -- see
! GUIDELINES.md's "Physics domains still entangled with the solver").
! Public facade of the atm domain (same shape as eos_lib/kap_lib):
! atm_init loads the tabulated-atmosphere surface-pressure tables at
! startup, and atm_get_surface_pt is the Allard-table surface lookup
! used from outside the domain. The envelope/atmosphere integrator
! itself (formerly envint, "atm_get") lives in core/envint_lib.f90;
! the table interpolators are in atm/tables/ (surfp, kcsurfp,
! alsurfp) and the T-tau relations in ttau_lib.f90.
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
! table sizes atm_table_nt/ntc x atm_table_ng from atm_table_lib.
!
! 2026 wave 3 (R5): the Kurucz and Kurucz/Castelli tables are the two
! surface_p_table instances atm_table%kurucz_surface_p (atm_table_nt
! rows) and castelli_surface_p (atm_table_ntc rows). Each is
! allocated on first use with exactly the extents the fixed arrays
! had and zero-filled as they were (allocate_surface_p); the READ
! statements are unchanged apart from the names; the identical
! G Somers 5/15 edge scans became the internal find_gravity_edges,
! its working scalars host-associated.
subroutine atm_init(atm_table_path, allard_table_path, ierr)
      use star_info_lib, only: star

      use atm_table_lib
      implicit none
      character(len=256), intent(in) :: atm_table_path
      character(len=256), intent(in) :: allard_table_path
! 2026 (ROADMAP.md stage 3): ierr-not-stop; table-load failures print
! their diagnostic at the point of failure and come back as ierr /= 0.
      integer, intent(out) :: ierr

      integer :: jerr
      integer :: teff_idx, logg_idx
      logical :: found_valid_pressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INPUT PRESSURE TABLE FOR SURFACE BOUNDARY CONDITIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ierr = 0
      jerr = 0

      if ((star%job%atm_choice .eq. 3) .or. (star%job%atm_choice .eq. 4)) then
        call allocate_surface_p(atm_table%kurucz_surface_p, atm_table_nt)
! OPEN SURFACE PRESSURE TABLE
        open(atm_table_file_unit,file=atm_table_path, status='OLD')
! GET ABUNDANCE:
        read(atm_table_file_unit,200) kurucz_table_z
  200   format(1x,10x,1pe16.8,/)
! GET VALUES OF LOG Teff:
        read(atm_table_file_unit,201) (atm_table%kurucz_surface_p%teff(teff_idx),teff_idx=1,atm_table_nt)
  201   format(1p4e16.8,13(/1p4e16.8),/1pe16.8,/)
! GET VALUES OF LOG G:
        read(atm_table_file_unit,202) (atm_table%kurucz_surface_p%logg(logg_idx),logg_idx=1,atm_table_ng)
  202   format(1p4e16.8,/1p4e16.8,/1p3e16.8,/)
! GET GRID OF LOG PRESSURE VALUES:
        do teff_idx=1,atm_table_nt
          read(atm_table_file_unit,203) (atm_table%kurucz_surface_p%log10_pressure(teff_idx,logg_idx),logg_idx=1,atm_table_ng)
  203     format(1p4e16.8,/1p4e16.8,/1p3e16.8)
        end do
        rewind atm_table_file_unit
        close(atm_table_file_unit)
!       G Somers 5/15; added a catch that allows the code to work
!                      if the highest gravity P value for a specific
!                      temperature is -999 (find_gravity_edges).
         call find_gravity_edges(atm_table%kurucz_surface_p)
!        G Somers 5/15 END
! MHP 6/97 ADDED OPTION FOR ALLARD MODEL ATMOSPHERES; USED INSTEAD OF
! KURUCZ FOR TEFF < 10,000 K.
         if(star%job%atm_choice .eq. 4)then
!            ATMZA = 0.02D0
          call alfilein(allard_table_path, jerr)      ! Get Allard Atmospheres files and
          if (jerr /= 0) then
             ierr = jerr
             return
          end if
         endif                  ! initialize tables. 9/23/08 LLP


! JNT 6/14 ADD OPTION FOR NEW KURUCZ/CASTELLI ATMOSHPERES
      else if (star%job%atm_choice .eq. 5) then
        call allocate_surface_p(atm_table%castelli_surface_p, atm_table_ntc)
! OPEN SURFACE PRESSURE TABLE
        open(atm_table_file_unit,file=atm_table_path, status='OLD')
! GET ABUNDANCE:
        read(atm_table_file_unit,200) kurucz_table_z
! GET VALUES OF LOG Teff:
        read(atm_table_file_unit,206) (atm_table%castelli_surface_p%teff(teff_idx),teff_idx=1,atm_table_ntc)
  206   format(1p4e16.8,18(/1p4e16.8),/)
! GET VALUES OF LOG G:
        read(atm_table_file_unit,208) (atm_table%castelli_surface_p%logg(logg_idx),logg_idx=1,atm_table_ng)
  208   format(1p4e16.8,/1p4e16.8,/1p3e16.8,/)
! GET GRID OF LOG PRESSURE VALUES:
        do teff_idx=1,atm_table_ntc
          read(atm_table_file_unit,203) (atm_table%castelli_surface_p%log10_pressure(teff_idx,logg_idx),logg_idx=1,atm_table_ng)
        end do
        rewind atm_table_file_unit
        close(atm_table_file_unit)
!       G Somers 5/15; added a catch that allows the code to work
!                      if the highest gravity P value for a specific
!                      temperature is -999 (find_gravity_edges).
        call find_gravity_edges(atm_table%castelli_surface_p)
! END JNT 6/14
      endif

      return

contains

! allocate_surface_p: on first use give tbl num_teff rows x
! atm_table_ng columns, zero-filled like the static arrays it
! replaces (an all -999 row leaves gmax_index untouched, so it must
! read 0 as before). Later calls (a re-read of the same table) reuse
! the storage as the static arrays did.
subroutine allocate_surface_p(tbl, num_teff)
      type(surface_p_table), intent(inout) :: tbl
      integer, intent(in) :: num_teff
      if (allocated(tbl%teff)) return
      allocate(tbl%teff(num_teff), tbl%logg(atm_table_ng), &
           tbl%log10_pressure(num_teff,atm_table_ng), &
           tbl%gmin_index(num_teff), tbl%gmax_index(num_teff))
      tbl%teff = 0.0d0
      tbl%logg = 0.0d0
      tbl%log10_pressure = 0.0d0
      tbl%gmin_index = 0
      tbl%gmax_index = 0
      tbl%num_teff = num_teff
end subroutine allocate_surface_p

! find_gravity_edges: G Somers 5/15; per log Teff row, the lowest and
! highest log g columns with a valid (positive) tabulated pressure.
! LCATCH (found_valid_pressure) is set to true once the first valid P
! is read, and it can set the minimum P thereafter.
subroutine find_gravity_edges(tbl)
      type(surface_p_table), intent(inout) :: tbl
         do teff_idx = 1,tbl%num_teff
            found_valid_pressure = .false.
            do logg_idx = atm_table_ng,1,-1
               if(tbl%log10_pressure(teff_idx,logg_idx).le.0.0d0)then
!                 check if the first non-999 value has been reached.
!                 if so, set IMIN.
                  if(found_valid_pressure)then
                     tbl%gmin_index(teff_idx) = logg_idx + 1
                     exit
                  endif
               else
!                 if we have reached a non-negative pressure value,
!                 turn off the catch so IMIN can be set. also record
!                 the highest gravity with a pressure for later int.
                  if(.not.found_valid_pressure) tbl%gmax_index(teff_idx) = logg_idx
                  found_valid_pressure = .true.
               endif
            end do
            if (logg_idx .lt. 1) tbl%gmin_index(teff_idx) = 1
!           if all of the P values at a given T are -999, set IMIN
!           to the number of gravity terms. in responce, the code
!           should break when trying to find surface P.
            if(.not.found_valid_pressure) tbl%gmin_index(teff_idx) = atm_table_ng
         end do
end subroutine find_gravity_edges

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
! 2026 (ROADMAP.md stage 3): ierr-not-stop. alsurfp's own ierr is
! deliberately not propagated here (it has always been ignored on
! this path); ierr is returned as 0 and lookup_failed carries the
! table-range result.
      integer, intent(out) :: ierr

      integer :: jerr

      ierr = 0
      jerr = 0

      call alsurfp(log_teff, log_g, print_to_files, lookup_failed, jerr)
      return
end subroutine atm_get_surface_pt

end module atm_lib
