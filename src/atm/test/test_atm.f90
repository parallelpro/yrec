!----------------------------------------------------------------------
! test_atm
!----------------------------------------------------------------------
! Standalone test of atm/, through atm_lib. Loads each tabulated
! atmosphere with atm_init and looks the surface boundary condition
! up with atm_get_surface_pt:
!   atm_choice=3 -- the Kurucz (1993) solar table;
!   atm_choice=5 -- the Castelli & Kurucz (2003) solar table;
!   atm_choice=4 -- the Allard NextGen solar tables (loaded together
!     with the Kurucz table, as in production): NextGen 2, NextGen 2
!     at tau = 100, the old NextGen 1 format, and a lookup outside the
!     tables;
!   atm_choice=0/1/2 -- the analytic T-tau relations (Eddington,
!     Krishna-Swamy, HSRA) on a grid of optical depths;
! and the error paths. The envelope integration that uses these lives
! in core and is covered by the full-model tests. Results print for
! byte-comparison against expected_test_atm.out.
program test_atm
      use io_units, only: run_log_unit, terminal_unit
      use atm_lib
      implicit none
      integer :: gerr
      character(len=256) :: yrec_input
      character(len=256) :: kurucz_atm_path, castelli_path, allard_path, &
           allard_old_path
      integer :: ipt, atm_ierr
      double precision :: teffl, gl, hras, tau, log10_p, log10_t
      integer :: choice, itau
      double precision, parameter :: grid_tau(6) = [1.0d-4, 1.0d-2, 0.312d0, &
           2.0d0/3.0d0, 1.0d0, 10.0d0]
      logical :: failed
      type(atm_settings_t) :: settings

      integer, parameter :: npts = 4
      double precision, parameter :: grid_teffl(npts) = [3.60d0, 3.70d0, 3.762d0, 3.90d0]
      double precision, parameter :: grid_gl(npts) = [4.9d0, 4.6d0, 4.438d0, 4.2d0]

      call get_environment_variable("YREC_INPUT", yrec_input)
      if (len_trim(yrec_input) == 0) then
         write(*,'(a)') "test_atm: FAIL (YREC_INPUT not set)"
         stop 1
      end if
      kurucz_atm_path = trim(yrec_input)//"/atmos/kurucz/atmk1990p00.tab"
      castelli_path   = trim(yrec_input)//"/atmos/CastelliKurucz/atmk2004p00.tab"
      allard_path     = trim(yrec_input)//"/atmos/allard/Nextgen2.all"
      allard_old_path = trim(yrec_input)//"/atmos/allard/Nextgen.p00"

! The diagnostics (the table readers' notes, the out-of-table messages)
! go to the run log and terminal_unit; both point at a scratch file so
! that stdout stays deterministic.
      open(newunit=run_log_unit, file="test_atm.short", status="replace")
      terminal_unit = run_log_unit

! --- Kurucz (1993), atm_choice=3 ---
      settings%atm_choice = 3
      call atm_init(settings, kurucz_atm_path, "", gerr)
      if (gerr /= 0) stop 1
      write(*,'(a)') "# test_atm: Kurucz (choice 3), surfp lookups"
      do ipt = 1, npts
         call lookup(3, grid_teffl(ipt), grid_gl(ipt))
         write(*,'(a,i2,2(1pe24.15))') "krz ", ipt, log10_p, log10_t
      end do

! --- Castelli & Kurucz (2003), atm_choice=5 ---
      settings%atm_choice = 5
      call atm_init(settings, castelli_path, "", gerr)
      if (gerr /= 0) stop 1
      write(*,'(a)') "# test_atm: Castelli/Kurucz (choice 5), " // &
           "kcsurfp lookups"
      do ipt = 1, npts
         call lookup(5, grid_teffl(ipt), grid_gl(ipt))
         write(*,'(a,i2,2(1pe24.15))') "ck  ", ipt, log10_p, log10_t
      end do

! --- Allard NextGen 2, atm_choice=4 ---
      settings%atm_choice = 4
      call atm_init(settings, kurucz_atm_path, allard_path, gerr)
      if (gerr /= 0) stop 1
      write(*,'(a)') "# test_atm: Allard (choice 4), " // &
           "atm_get_surface_pt lookups"
      do ipt = 1, npts
         call lookup(4, grid_teffl(ipt), grid_gl(ipt))
         write(*,'(a,i2,l2,2(1pe24.15))') "al  ", ipt, failed, log10_p, log10_t
      end do

! --- Allard NextGen 2 with P and T at tau = 100 (allard_use_tau100) ---
      settings%allard_use_tau100 = .true.
      call atm_init(settings, kurucz_atm_path, allard_path, gerr)
      if (gerr /= 0) stop 1
      write(*,'(a)') "# test_atm: Allard at tau = 100, atm_get_surface_pt lookups"
      do ipt = 1, npts
         call lookup(4, grid_teffl(ipt), grid_gl(ipt))
         write(*,'(a,i2,l2,2(1pe24.15))') "a100", ipt, failed, log10_p, log10_t
      end do
      settings%allard_use_tau100 = .false.

! --- Allard NextGen 1, the old table format ---
      call atm_init(settings, kurucz_atm_path, allard_old_path, gerr)
      if (gerr /= 0) stop 1
      write(*,'(a)') "# test_atm: Allard NextGen 1 (old format), atm_get_surface_pt lookups"
      do ipt = 1, npts
         call lookup(4, grid_teffl(ipt), grid_gl(ipt))
         write(*,'(a,i2,l2,2(1pe24.15))') "al1 ", ipt, failed, log10_p, log10_t
      end do
! a hot star is outside the Allard tables: the lookup fails (the caller
! falls back to a gray atmosphere)
      call atm_get_surface_pt(4, 4.5d0, 4.0d0, log10_p, log10_t, failed, gerr)
      write(*,'(a,l2,i4)') "al1 out of range: failed, ierr =", failed, gerr

! --- the analytic T-tau relations (atm_choice 0 Eddington, 1
! Krishna-Swamy, 2 HSRA): log T at a grid of optical depths, the
! starting temperature and the photospheric limit, at log Teff 3.762
      write(*,'(a)') "# test_atm: T-tau relations at log Teff = 3.762"
      hras = hsra_t_tau_offset()
      write(*,'(a,1pe24.15)') "hsra offset ", hras
      do choice = 0, 2
         write(*,'(a,i2,2(1pe24.15))') "ttau start, x_limit", choice, &
              ttau_start_log10_temperature(3.762d0, choice, hras), &
              ttau_photosphere_x_limit(choice)
         do itau = 1, size(grid_tau)
            tau = grid_tau(itau)
            write(*,'(a,i2,2(1pe24.15))') "ttau", choice, tau, &
                 ttau_log10_temperature(tau, 3.762d0, choice, hras)
         end do
      end do

! --- error paths: a Kurucz lookup outside the table (log Teff < 3.5 or
! log g < -0.5) returns ierr = 1 with its diagnostic in the run log; a
! lookup inside it returns 0
      write(*,'(a)') "# test_atm: error paths via ierr"
      settings%atm_choice = 3
      call atm_init(settings, kurucz_atm_path, "", gerr)
      if (gerr /= 0) stop 1
      teffl = 3.40d0
      gl = 4.5d0
      call atm_get_surface_pt(3, teffl, gl, log10_p, log10_t, failed, atm_ierr)
      write(*,'(a,i4)') "err out-of-table   ierr = ", atm_ierr
      call atm_get_surface_pt(3, grid_teffl(1), grid_gl(1), log10_p, log10_t, &
           failed, atm_ierr)
      write(*,'(a,i4)') "err facade-success ierr = ", atm_ierr

      close(run_log_unit)
      write(*,'(a)') "test_atm: done"
contains

subroutine lookup(atm_choice, log_teff, log_g)
      integer, intent(in) :: atm_choice
      double precision, intent(in) :: log_teff, log_g
      call atm_get_surface_pt(atm_choice, log_teff, log_g, log10_p, log10_t, &
           failed, gerr)
      if (gerr /= 0) stop 1
end subroutine lookup

end program test_atm
