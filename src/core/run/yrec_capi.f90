!----------------------------------------------------------------------
! yrec_capi
!----------------------------------------------------------------------
! Added 2026 (libyrec milestone -- docs/ROADMAP.md). The C-interoperable
! surface of the embeddable engine, consumed by pyyrec (the ctypes
! binding in pyyrec/ at the repository root). One entry:
!
!   int yrec_run(const char *control_nml, const char *physics_nml);
!
! Runs one full YREC job (equivalent to `yrec <nml1> <nml2>` from the
! current working directory) and returns 0 on success, nonzero on the
! errors that run_yrec reports via ierr. The namelist paths are
! injected through star_info_lib's overrides rather than argv, since
! an embedded engine must not read the host process's command line.
! Empty strings fall back to the CLI defaults (yrec8.nml1/.nml2).
!
! Re-entrancy contract: repeated yrec_run calls in one process are
! supported and byte-identical to fresh processes (phase-five step C;
! enforced by test_reentry and test_pyyrec). The numerics gates
! report through ierr rather than stopping the process: a NEGATIVE
! status is envelope_termination (envelope_integration_kernel), i.e. the run ended
! in the historical "solution diverged" mode (outputs up to that
! model are valid); positive is a hard error; zero is success.
! A namelist path longer than 256 characters (the override slots) is
! refused with status 1 before anything runs: cut, it would name another
! file. After the run the C entry closes the run's own units
! (yrec_reset's close_run_units), not the host's. The two arguments are
! the C strings' addresses (type(c_ptr), passed by value, as a char *
! is); their lengths are C's strlen.
module yrec_capi
      use run_driver, only: run_yrec
      use, intrinsic :: iso_c_binding, only: c_char, c_f_pointer, c_int, c_ptr, c_size_t
      use star_info_lib, only: control_nml_override, physics_nml_override
      use yrec_reset, only: close_run_units
      use core_messages, only: core_terminal, message_length
      implicit none
      private
      public :: yrec_run

      interface
         ! C's strlen: the length of a NUL-terminated string.
         integer(c_size_t) function c_strlen(cstr) bind(c, name='strlen')
            import :: c_ptr, c_size_t
            type(c_ptr), value :: cstr
         end function c_strlen
      end interface

contains

      integer(c_int) function yrec_run(control_nml, physics_nml) &
           bind(c, name='yrec_run')
            type(c_ptr), value :: control_nml, physics_nml
            integer :: ierr
            character(len=message_length) :: line

            if (c_strlen(control_nml) > len(control_nml_override) .or. &
                 c_strlen(physics_nml) > len(physics_nml_override)) then
               write(line,'(a,i0,a)') ' ERROR: yrec_run: a namelist path is longer than ', &
                    len(control_nml_override), ' characters; nothing run'
               call core_terminal([line])
               yrec_run = 1_c_int
               return
            end if
            control_nml_override = c_string_value(control_nml)
            physics_nml_override = c_string_value(physics_nml)
            ierr = 0
            call run_yrec(ierr)
            ! A library call must leave no dangling open output units:
            ! the CLI relies on process exit to flush/close, and the
            ! re-entrancy prologue only closes them on the NEXT entry
            ! -- an embedder reading the output files right after this
            ! return would otherwise see unflushed tails. Only the
            ! run's own units are closed; the host's stay open.
            call close_run_units
            ! Leave the overrides clean so a later CLI-style use of the
            ! same process (or a blank-path call) behaves historically.
            control_nml_override = ' '
            physics_nml_override = ' '
            yrec_run = int(ierr, c_int)
      end function yrec_run

      ! The C string at cstr (NUL-terminated, at most 256 characters,
      ! which yrec_run has checked), blank-padded.
      function c_string_value(cstr) result(fstr)
            type(c_ptr), intent(in) :: cstr
            character(len=256) :: fstr
            character(kind=c_char), pointer :: chars(:)
            integer :: i, n

            n = int(c_strlen(cstr))
            fstr = ' '
            call c_f_pointer(cstr, chars, [n])
            do i = 1, n
               fstr(i:i) = chars(i)
            end do
      end function c_string_value

end module yrec_capi
