!----------------------------------------------------------------------
! test_timestep_limits
!----------------------------------------------------------------------
! The timestep limits (core/step/timestep_limits.f90) on synthetic
! models of a few points; no star, no tables (the hydrogen limit's shell
! branch takes its point below the network's first temperature cutoff,
! where the network burns nothing). Each limit against its formula, to
! a relative 1e-14 (round-off of a few operations):
!   structure  off -> unlimited; nothing changed -> the previous dt; the
!              largest change of log P, log T, log R (each over its
!              factor) and of L relative to the mean of the point and
!              the one below (over ln 10 times its factor; the centre
!              against L(2)) sets f and dt = dt_prev/f; f clipped to
!              [1/a, a]; a limit below 1 s -> unlimited;
!   omega      off -> unlimited, max_domega_frac kept; the largest
!              2|w - w0|/(w + w0) and dt = dt_prev max_domega_global/that;
!              no change -> unlimited;
!   hydrogen   L_H at or below 1e-34 -> unlimited; the core branch
!              E_H Delta X m_cc/(L_cc Lsun) with the smaller Delta X; the
!              shell branch E_H atime(dx_total) X_surf Msun/(L_H Lsun);
!   helium     off -> unlimited; the flash's 1e15/eps_3alpha, the rate
!              floored at 1e-22; the core branch; the shell branch's
!              smaller limit, and none where L <= 0; the hottest zone;
!   hr         off unless asked for with the structure limits off; f from
!              log Teff or log L over the triangle's size; nothing
!              changed -> the previous dt;
!   tightest   the growth limit a dt_prev for dt_prev > 1 s, the smallest
!              limit otherwise.
! Each check prints its relative error (or its value) and ok or FAIL;
! any FAIL stops the program with status 1.
program test_timestep_limits
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use timestep_limits, only: timestep_settings, structure_change_limit, &
           omega_change_limit, hydrogen_burning_limit, helium_burning_limit, &
           hr_motion_limit, tightest_timestep, hottest_zone, dt_unlimited
      use control_slots, only: itime_core_min, itime_dx_core_tot, &
           itime_dx_core_frac, itime_dy_core_tot, itime_dy_core_frac, &
           itime_dx_total, itime_dx_shell, itime_dt, itime_dp, itime_dr, &
           itime_dl, itime_dy_shell, itime_max_dt_frac, itime_dy_total
      use composition_layout, only: n_species_extended, i_h1, i_he4, i_metals
      implicit none

      integer, parameter :: n = 5
      real(dp), parameter :: dt_prev = 1.0e13_dp, a = 1.2_dp
      real(dp), parameter :: lsun = 3.828e33_dp, msun = 1.989e33_dp
      type(timestep_settings) :: set
      real(dp) :: lp0(n), lt0(n), lr0(n), l0(n), lp(n), lt(n), lr(n), l(n)
      real(dp) :: xa(n_species_extended,n), m(n), lrho(n)
      real(dp) :: w0(n), w(n), frac, dt, expected
      integer :: nfail

      nfail = 0
      call default_settings()
      write(*,'(a)') 'test_timestep_limits report'
      call check_structure()
      call check_omega()
      call check_hydrogen()
      call check_helium()
      call check_hr()
      call check_tightest()
      if (nfail > 0) then
         write(*,'(a,i0,a)') 'test_timestep_limits: ', nfail, ' FAIL'
         stop 1
      end if
      write(*,'(a)') 'test_timestep_limits: done'

contains

subroutine default_settings()
      set%timestep_factors = 0.0_dp
      set%timestep_factors(itime_core_min) = 1.0e-3_dp
      set%timestep_factors(itime_dx_core_tot) = 2.0e-2_dp
      set%timestep_factors(itime_dx_core_frac) = 5.0e-1_dp
      set%timestep_factors(itime_dy_core_tot) = 2.0e-2_dp
      set%timestep_factors(itime_dy_core_frac) = 3.0e-1_dp
      set%timestep_factors(itime_dx_total) = 1.0e-4_dp
      set%timestep_factors(itime_dx_shell) = 1.0e-2_dp
      set%timestep_factors(itime_dt) = 2.0e-2_dp
      set%timestep_factors(itime_dp) = 5.0e-2_dp
      set%timestep_factors(itime_dr) = 3.0e-2_dp
      set%timestep_factors(itime_dl) = 1.0e-1_dp
      set%timestep_factors(itime_dy_shell) = 4.0e-2_dp
      set%timestep_factors(itime_max_dt_frac) = a
      set%timestep_factors(itime_dy_total) = 2.0e-4_dp
      set%max_domega_global = 5.0e-2_dp
      set%tri_delta_teffl = 1.0e-2_dp
      set%tri_delta_logl = 8.0e-2_dp
      set%solar_luminosity_cgs = lsun
      set%solar_mass_cgs = msun
      set%structure_limits_active = .true.
      set%rotation_active = .true.
      set%use_envelope_triangle_dt = .true.
      lp0 = [17.0_dp, 16.5_dp, 15.0_dp, 12.0_dp, 5.0_dp]
      lt0 = [7.2_dp, 7.1_dp, 6.9_dp, 6.3_dp, 4.0_dp]
      lr0 = [8.5_dp, 9.6_dp, 10.2_dp, 10.6_dp, 10.84_dp]
      l0 = [0.1_dp, 0.6_dp, 0.95_dp, 1.0_dp, 1.0_dp]
      m = msun*[0.01_dp, 0.2_dp, 0.6_dp, 0.95_dp, 1.0_dp]
      lrho = [2.2_dp, 1.8_dp, 0.5_dp, -1.0_dp, -7.0_dp]
      xa = 0.0_dp
      xa(i_h1,:) = 0.7_dp
      xa(i_metals,:) = 0.02_dp
      xa(i_he4,:) = 0.28_dp
end subroutine default_settings

subroutine report(name, value, pass)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: value
      logical, intent(in) :: pass
      if (pass) then
         write(*,'(a,t62,1pe9.1,2x,a)') name, value, 'ok'
      else
         write(*,'(a,t62,1pe9.1,2x,a)') name, value, 'FAIL'
         nfail = nfail + 1
      end if
end subroutine report

! The relative error of got against want.
real(dp) function relerr(got, want)
      real(dp), intent(in) :: got, want
      relerr = abs(got - want)/abs(want)
end function relerr

! The structure limit of the current arrays against the start ones.
real(dp) function structure_dt(previous)
      real(dp), intent(in) :: previous
      structure_dt = structure_change_limit(set, previous, lp0, lt0, lr0, l0, &
           lp, lt, lr, l, n)
end function structure_dt

subroutine reset_structure()
      lp = lp0
      lt = lt0
      lr = lr0
      l = l0
end subroutine reset_structure

subroutine check_structure()
      write(*,'(a)') '-- structure'
      call reset_structure()
      set%structure_limits_active = .false.
      lp(3) = lp0(3) + 1.0_dp
      dt = structure_dt(dt_prev)
      call report('structure: off -> unlimited', dt, dt == dt_unlimited)
      set%structure_limits_active = .true.
      call reset_structure()
      dt = structure_dt(dt_prev)
      call report('structure: nothing changed -> the previous dt', relerr(dt, dt_prev), &
           dt == dt_prev)
      lp(3) = lp0(3) + 1.1_dp*set%timestep_factors(itime_dp)
      dt = structure_dt(dt_prev)
      expected = dt_prev/(abs(lp0(3) - lp(3))/set%timestep_factors(itime_dp))
      call report('structure: f = dlog P/atime(dp)', relerr(dt, expected), &
           relerr(dt, expected) <= 1.0e-14_dp)
      call reset_structure()
      lt(2) = lt0(2) - 1.05_dp*set%timestep_factors(itime_dt)
      lp(4) = lp0(4) + 0.5_dp*set%timestep_factors(itime_dp)
      dt = structure_dt(dt_prev)
      expected = dt_prev/(abs(lt0(2) - lt(2))/set%timestep_factors(itime_dt))
      call report('structure: f = dlog T/atime(dt), the largest change', relerr(dt, expected), &
           relerr(dt, expected) <= 1.0e-14_dp)
      call reset_structure()
      lr(5) = lr0(5) + 1.15_dp*set%timestep_factors(itime_dr)
      dt = structure_dt(dt_prev)
      expected = dt_prev/(abs(lr0(5) - lr(5))/set%timestep_factors(itime_dr))
      call report('structure: f = dlog R/atime(dr)', relerr(dt, expected), &
           relerr(dt, expected) <= 1.0e-14_dp)
      call reset_structure()
      l(3) = l0(3)*1.25_dp
      dt = structure_dt(dt_prev)
      expected = dt_prev/(abs((l0(3) - l(3))*2.0_dp/(l(3) + l(2)))/ &
           (set%timestep_factors(itime_dl)*2.3026_dp))
      call report('structure: f = (dL/mean L)/(2.3026 atime(dl))', relerr(dt, expected), &
           relerr(dt, expected) <= 1.0e-14_dp)
      call reset_structure()
      l(1) = l0(1)*2.0_dp
      dt = structure_dt(dt_prev)
      expected = dt_prev/(abs((l0(1) - l(1))*2.0_dp/(l(2) + l(1)))/ &
           (set%timestep_factors(itime_dl)*2.3026_dp))
      call report('structure: the centre''s L against L(2)', relerr(dt, expected), &
           relerr(dt, expected) <= 1.0e-14_dp)
      call reset_structure()
      lp(2) = lp0(2) + 10.0_dp
      dt = structure_dt(dt_prev)
      call report('structure: f clipped to a', relerr(dt, dt_prev/a), &
           relerr(dt, dt_prev/a) <= 1.0e-15_dp)
      call reset_structure()
      lp(2) = lp0(2) + 1.0e-6_dp
      dt = structure_dt(dt_prev)
      call report('structure: f clipped to 1/a', relerr(dt, dt_prev*a), &
           relerr(dt, dt_prev*a) <= 1.0e-15_dp)
      dt = structure_dt(0.5_dp)
      call report('structure: a limit below 1 s -> unlimited', dt, dt == dt_unlimited)
end subroutine check_structure

subroutine check_omega()
      write(*,'(a)') '-- omega'
      w0 = [3.0e-6_dp, 2.9e-6_dp, 2.8e-6_dp, 2.7e-6_dp, 2.6e-6_dp]
      w = w0
      w(4) = w0(4)*1.05_dp
      w(2) = w0(2)*1.01_dp
      set%rotation_active = .false.
      frac = -1.0_dp
      call omega_change_limit(set, dt_prev, w0, w, n, dt, frac)
      call report('omega: off -> unlimited, max_domega_frac kept', dt, &
           dt == dt_unlimited .and. frac == -1.0_dp)
      set%rotation_active = .true.
      call omega_change_limit(set, dt_prev, w0, w, n, dt, frac)
      expected = 2.0_dp*abs(w(4) - w0(4))/(w(4) + w0(4))
      call report('omega: max_domega_frac = max 2|w - w0|/(w + w0)', relerr(frac, expected), &
           relerr(frac, expected) <= 1.0e-15_dp)
      expected = dt_prev/(expected/set%max_domega_global)
      call report('omega: dt = dt_prev max_domega_global/max_domega_frac', &
           relerr(dt, expected), relerr(dt, expected) <= 1.0e-14_dp)
      call omega_change_limit(set, dt_prev, w0, w0, n, dt, frac)
      call report('omega: no change -> unlimited', dt, dt == dt_unlimited)
end subroutine check_omega

subroutine check_hydrogen()
      real(dp) :: delta_x
      write(*,'(a)') '-- hydrogen'
      dt = hydrogen_burning_limit(set, 1.0e-35_dp, xa, xa, lrho, lt0, m, l0, 2, 3, n)
      call report('hydrogen: L_H <= 1e-34 -> unlimited', dt, dt == dt_unlimited)
! Core burning: Delta X = min(0.02, 0.5 X_cc) = 0.02.
      dt = hydrogen_burning_limit(set, 1.0_dp, xa, xa, lrho, lt0, m, l0, 2, 3, n)
      delta_x = set%timestep_factors(itime_dx_core_tot)
      expected = 6.0e18_dp/lsun*(m(2)/l0(2))*delta_x
      call report('hydrogen: core E_H dX m_cc/(L_cc Lsun)', &
           relerr(dt, expected), relerr(dt, expected) <= 1.0e-14_dp)
      xa(i_h1,2) = 0.01_dp
      dt = hydrogen_burning_limit(set, 1.0_dp, xa, xa, lrho, lt0, m, l0, 2, 3, n)
      delta_x = set%timestep_factors(itime_dx_core_frac)*0.01_dp
      expected = 6.0e18_dp/lsun*(m(2)/l0(2))*delta_x
      call report('hydrogen: core dX = atime(dx_core_frac) X_cc', &
           relerr(dt, expected), relerr(dt, expected) <= 1.0e-14_dp)
! Shell burning (central X below atime(core_min)); the shell midpoint,
! point 4, lies below the network's first cutoff (log T 6.3 < 6.5).
      xa(i_h1,1) = 1.0e-4_dp
      dt = hydrogen_burning_limit(set, 0.5_dp, xa, xa, lrho, lt0, m, l0, 2, 4, n)
      expected = 6.0e18_dp*set%timestep_factors(itime_dx_total)*xa(i_h1,n)*(msun/lsun)/0.5_dp
      call report('hydrogen: shell E_H atime(dx_total) X_surf Msun/(L_H Lsun)', &
           relerr(dt, expected), relerr(dt, expected) <= 1.0e-14_dp)
      xa(i_h1,:) = 0.7_dp
end subroutine check_hydrogen

subroutine check_helium()
      write(*,'(a)') '-- helium'
      dt = helium_burning_limit(set, 1.0e-35_dp, 1.0_dp, xa, lrho, m, l0, 2, 3)
      call report('helium: L_3alpha <= 1e-34 -> unlimited', dt, dt == dt_unlimited)
      lrho(1) = 5.5_dp
      dt = helium_burning_limit(set, 1.0_dp, 4.0e3_dp, xa, lrho, m, l0, 2, 3)
      call report('helium: flash 1e15/eps_3alpha', relerr(dt, 1.0e15_dp/4.0e3_dp), &
           relerr(dt, 1.0e15_dp/4.0e3_dp) <= 1.0e-15_dp)
      dt = helium_burning_limit(set, 1.0_dp, 1.0e-30_dp, xa, lrho, m, l0, 2, 3)
      call report('helium: flash rate floored at 1e-22', relerr(dt, 1.0e37_dp), &
           relerr(dt, 1.0e37_dp) <= 1.0e-15_dp)
      lrho(1) = 2.2_dp
! Core burning: Y_core = 1 - Z(2) = 0.98, Delta Y = min(0.02, 0.3 0.98).
      dt = helium_burning_limit(set, 1.0_dp, 0.0_dp, xa, lrho, m, l0, 2, 3)
      expected = (5.85e17_dp/lsun)*set%timestep_factors(itime_dy_core_tot)*(m(2)/l0(2))
      call report('helium: core E_He dY m_cc/(L_cc Lsun)', relerr(dt, expected), &
           relerr(dt, expected) <= 1.0e-14_dp)
! Shell burning: core helium gone (1 - Z below atime(core_min)).
      xa(i_metals,2) = 1.0_dp
      dt = helium_burning_limit(set, 1.0_dp, 0.0_dp, xa, lrho, m, l0, 2, 4)
      expected = min((5.85e17_dp/lsun)*set%timestep_factors(itime_dy_total)*(msun/l0(3)), &
           (5.85e17_dp/lsun)*set%timestep_factors(itime_dy_shell)*xa(i_he4,3)*(m(3)/l0(3)))
      call report('helium: shell the smaller of the two limits', relerr(dt, expected), &
           relerr(dt, expected) <= 1.0e-14_dp)
      l0(3) = -0.5_dp
      dt = helium_burning_limit(set, 1.0_dp, 0.0_dp, xa, lrho, m, l0, 2, 4)
      call report('helium: shell with L <= 0 -> unlimited', dt, dt == dt_unlimited)
      l0(3) = 0.95_dp
      xa(i_metals,2) = 0.02_dp
      call report('hottest_zone: the first largest log T', &
           real(hottest_zone([6.0_dp, 7.0_dp, 7.0_dp, 6.5_dp], 4), dp), &
           hottest_zone([6.0_dp, 7.0_dp, 7.0_dp, 6.5_dp], 4) == 2)
end subroutine check_helium

subroutine check_hr()
      write(*,'(a)') '-- hr'
      set%structure_limits_active = .true.
      dt = hr_motion_limit(set, dt_prev, 3.76_dp, 1.0_dp, 3.77_dp, 1.1_dp)
      call report('hr: off while the structure limits are on', dt, dt == dt_unlimited)
      set%structure_limits_active = .false.
      set%use_envelope_triangle_dt = .false.
      dt = hr_motion_limit(set, dt_prev, 3.76_dp, 1.0_dp, 3.77_dp, 1.1_dp)
      call report('hr: off unless use_envelope_triangle_dt', dt, dt == dt_unlimited)
      set%use_envelope_triangle_dt = .true.
      dt = hr_motion_limit(set, dt_prev, 3.76_dp, 1.0_dp, 3.76_dp + 1.1_dp*set%tri_delta_teffl, 1.0_dp)
      expected = dt_prev/(abs(-1.1_dp*set%tri_delta_teffl)/set%tri_delta_teffl)
      call report('hr: f = dlog Teff/tri_delta_teffl', relerr(dt, expected), &
           relerr(dt, expected) <= 1.0e-12_dp)
      dt = hr_motion_limit(set, dt_prev, 3.76_dp, 1.0_dp, 3.76_dp, 10.0_dp**(0.09_dp))
      expected = dt_prev/(abs(log10(10.0_dp**0.09_dp))/set%tri_delta_logl)
      call report('hr: f = dlog L/tri_delta_logl', relerr(dt, expected), &
           relerr(dt, expected) <= 1.0e-12_dp)
      dt = hr_motion_limit(set, dt_prev, 3.76_dp, 1.0_dp, 3.76_dp, 1.0_dp)
      call report('hr: nothing changed -> the previous dt', relerr(dt, dt_prev), dt == dt_prev)
      set%structure_limits_active = .true.
end subroutine check_hr

subroutine check_tightest()
      write(*,'(a)') '-- tightest'
      dt = tightest_timestep(set, dt_prev, [2.0_dp*dt_prev, 3.0_dp*dt_prev])
      call report('tightest: the growth limit a dt_prev', relerr(dt, a*dt_prev), &
           dt == a*dt_prev)
      dt = tightest_timestep(set, dt_prev, [2.0_dp*dt_prev, 0.5_dp*dt_prev, dt_unlimited])
      call report('tightest: the smallest limit', relerr(dt, 0.5_dp*dt_prev), &
           dt == 0.5_dp*dt_prev)
      dt = tightest_timestep(set, 0.5_dp, [7.0_dp, 3.0_dp])
      call report('tightest: dt_prev <= 1 s: the smallest limit alone', dt, dt == 3.0_dp)
end subroutine check_tightest

end program test_timestep_limits
