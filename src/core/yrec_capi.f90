!----------------------------------------------------------------------
! yrec_capi
!----------------------------------------------------------------------
! Added 2026 (libyrec milestone -- ROADMAP.md). The C-interoperable
! surface of the embeddable engine, consumed by pyyrec (the ctypes
! binding under /Applications/YREC/pyyrec). One entry:
!
!   int yrec_run(const char *control_nml, const char *physics_nml);
!
! Runs one full YREC job (equivalent to `yrec <nml1> <nml2>` from the
! current working directory) and returns 0 on success, nonzero on the
! errors that run_yrec reports via ierr. The namelist paths are
! injected through star_job_lib's overrides rather than argv, since
! an embedded engine must not read the host process's command line.
! Empty strings fall back to the CLI defaults (yrec8.nml1/.nml2).
!
! Re-entrancy contract: repeated yrec_run calls in one process are
! supported and byte-identical to fresh processes (phase-five step C;
! enforced by test_reentry and test_pyyrec). Known residual, same as
! the CLI: legacy `stop` endings behind the numerics gates (e.g. the
! m0030 case's historical BSSTEP termination) end the *process*, not
! just the call -- avoid such configurations when embedding until the
! numerics-gate ierr opt-in lands. 2026: it has landed -- a NEGATIVE
! status is numerics_termination (numerics_lib): the run ended in the
! historical "solution diverged" mode (outputs up to that model are
! valid); positive is a hard error; zero is success.
module yrec_capi
      use iso_c_binding, only: c_int, c_char, c_null_char
      use star_info_lib, only: control_nml_override, physics_nml_override
      implicit none
      private
      public :: yrec_run

contains

      integer(c_int) function yrec_run(control_nml, physics_nml) &
           bind(c, name='yrec_run')
            character(kind=c_char), intent(in) :: control_nml(*)
            character(kind=c_char), intent(in) :: physics_nml(*)
            integer :: ierr

            control_nml_override = c_string_value(control_nml)
            physics_nml_override = c_string_value(physics_nml)
            ierr = 0
            call run_yrec(ierr)
            ! A library call must leave no dangling open output units:
            ! the CLI relies on process exit to flush/close, and the
            ! re-entrancy prologue only closes them on the NEXT entry
            ! -- an embedder reading the output files right after this
            ! return would otherwise see unflushed tails.
            call close_open_units
            ! Leave the overrides clean so a later CLI-style use of the
            ! same process (or a blank-path call) behaves historically.
            control_nml_override = ' '
            physics_nml_override = ' '
            yrec_run = int(ierr, c_int)
      end function yrec_run

      subroutine close_open_units
            ! Same sweep as yrec_run_prologue's: every unit the engine
            ! may have written (7-99), skipping stdin/stdout/stderr.
            integer :: u
            logical :: is_open

            do u = 7, 99
               inquire(unit=u, opened=is_open)
               if (is_open) then
                  close(u)
               end if
            end do
      end subroutine close_open_units

      function c_string_value(cstr) result(fstr)
            character(kind=c_char), intent(in) :: cstr(*)
            character(len=256) :: fstr
            integer :: i

            fstr = ' '
            do i = 1, len(fstr)
               if (cstr(i) == c_null_char) exit
               fstr(i:i) = cstr(i)
            end do
      end function c_string_value

end module yrec_capi
