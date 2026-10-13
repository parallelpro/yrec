!----------------------------------------------------------------------
! test_mass_change
!----------------------------------------------------------------------
! The mass-change kernels (core/step/mass_change_kernel.f90) against
! the conservation laws they must keep, on a synthetic model of 40
! shells: m(k) = 0.9999 M (k/40)^2 with M = 2e33 g, the shell masses dm
! as rezoning sets them (dm(1) = (m(1) + m(2))/2, the last reaching M),
! r(k) = R (0.02 + 0.97 k/40) with R = 7e10 cm, and a surface
! convection zone from shell jenv to 40 (jenv = 25, or 1 for a fully
! convective star) that rotates as omega ~ r^walpcz (omega(40) =
! 3e-6 rad/s) with j_rot = (2/3) r^2 omega, and is uniform in
! composition. m, dm, r, omega and j_rot are built with arithmetic only
! (the kernels are also given log10 m and log10 r), and each expected
! value is formed here from the model, not from the kernels. Checks:
!   mass         after a loss or gain dM = -+1e-3 M (change_shell_masses):
!                the shell masses sum to M + dM; with a loss the
!                zone-center masses above the base scale about m(jenv)
!                by (M_cz + dM)/M_cz (M_cz = m(40) - m(jenv), or m(40)
!                for jenv = 1), with a gain they stay; the masses at and
!                below the base stay; for jenv = 1, dm(1) = (m(1) +
!                m(2))/2;
!   angular momentum  after a loss of 1e-3 M from a rotating star, for
!                walpcz = 0, -1 and -2 (wind_specific_angular_momentum,
!                change_shell_masses, remove_wind_angular_momentum):
!                dJ = j_s dM with j_s = (2/3) R^2 omega(40) (R/r(40))^walpcz,
!                J_cz = sum of j_rot dm over the zone before the change,
!                and the zone holds J_cz + dJ after it, sum of j_rot dm
!                on its new masses;
!   composition  after accreting 1e-3 M (change_shell_masses,
!                mix_accreted_composition): for each basic species, the
!                zone's sum of X dm is the old sum plus X_acc dM.
! Each check prints its relative error and ok (below 1e-13, round-off)
! or FAIL; the program stops with status 1 on any FAIL.
program test_mass_change
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: n_species_basic, n_species_extended
      use mass_change_kernel, only: wind_angular_momentum, &
           convection_zone_mass, change_shell_masses, &
           wind_specific_angular_momentum, remove_wind_angular_momentum, &
           mix_accreted_composition
      implicit none

      integer, parameter :: nz = 40
      real(dp), parameter :: tolerance = 1.0d-13
      real(dp), parameter :: total_mass = 2.0d33, fit_mass_fraction = 0.9999d0
      real(dp), parameter :: total_radius = 7.0d10
      real(dp), parameter :: surface_omega = 3.0d-6
      real(dp), parameter :: mass_change_fraction = 1.0d-3
! The zone's composition and the accreted gas's (the default
! accreted_composition), slots of composition_layout.
      real(dp), parameter :: x_zone(n_species_extended) = [0.70d0, 0.28d0, &
           0.02d0, 3.0d-5, 3.0d-3, 4.0d-5, 1.0d-3, 4.0d-6, 9.0d-3, 4.0d-6, &
           2.0d-5, 2.0d-5, 1.0d-9, 1.0d-8, 1.0d-10]
      real(dp), parameter :: x_accreted(n_species_extended) = [0.71668d0, &
           0.265721d0, 0.01757d0, 2.9d-5, 3.013d-3, 3.385d-5, 9.346d-4, &
           0.0d0, 8.462d-3, 0.0d0, 1.696d-5, 0.0d0, 2.0d-9, 2.0d-9, 3.0d-11]
      character(len=3), parameter :: species_name(n_species_basic) = [ &
           'h1 ', 'he4', 'Z  ', 'he3', 'c12', 'c13', 'n14', 'n15', 'o16', &
           'o17', 'o18']

      real(dp) :: m(nz), log_mass(nz), dm(nz), radius(nz), log_radius(nz), &
           omega(nz), j_rot(nz), xa(n_species_extended,nz)
      integer :: nfail

      nfail = 0
      call check_shell_masses(25, -mass_change_fraction)
      call check_shell_masses(25, mass_change_fraction)
      call check_shell_masses(1, -mass_change_fraction)
      call check_shell_masses(1, mass_change_fraction)
      call check_wind_angular_momentum(25, 0.0d0)
      call check_wind_angular_momentum(25, -1.0d0)
      call check_wind_angular_momentum(25, -2.0d0)
      call check_wind_angular_momentum(1, 0.0d0)
      call check_accreted_composition(25)
      call check_accreted_composition(1)
      if (nfail > 0) then
         write(*,'(a,i4,a)') 'test_mass_change: ', nfail, ' FAIL'
         stop 1
      end if
      write(*,'(a)') 'test_mass_change: done'

contains

! ---------------------------------------------------------------
! The synthetic model with its convection zone from shell jenv up,
! rotating as omega ~ r^walpcz (walpcz a whole number here).
subroutine build_model(jenv, walpcz)
      integer, intent(in) :: jenv
      real(dp), intent(in) :: walpcz
      integer :: k

      do k = 1, nz
         m(k) = fit_mass_fraction*total_mass*real(k*k, dp)/real(nz*nz, dp)
         log_mass(k) = log10(m(k))
         radius(k) = total_radius*(0.02d0 + 0.97d0*real(k, dp)/real(nz, dp))
         log_radius(k) = log10(radius(k))
      end do
      dm(1) = 0.5d0*(m(1) + m(2))
      do k = 2, nz-1
         dm(k) = 0.5d0*(m(k+1) - m(k-1))
      end do
      dm(nz) = total_mass - 0.5d0*(m(nz) + m(nz-1))
      do k = 1, nz
         omega(k) = surface_omega*(radius(max(k,jenv))/radius(nz))**nint(walpcz)
         j_rot(k) = 2.0d0/3.0d0*radius(k)**2*omega(k)
      end do
      do k = 1, nz
         if (k >= jenv) then
            xa(:,k) = x_zone
         else
            xa(:,k) = 0.5d0*x_zone
            xa(1,k) = 0.35d0
            xa(2,k) = 1.0d0 - 0.35d0 - 0.5d0*x_zone(3)
         end if
      end do
end subroutine build_model

! ---------------------------------------------------------------
! One check: its relative error against the tolerance.
subroutine report(label, relative_error)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: relative_error
      character(len=4) :: status

      if (relative_error <= tolerance) then
         status = 'ok'
      else
         status = 'FAIL'
         nfail = nfail + 1
      end if
      write(*,'(2x,a,t52,1pe9.1,2x,a)') label, relative_error, trim(status)
end subroutine report

! ---------------------------------------------------------------
! The sum of x(k) y(k) over k, in order.
pure function weighted_sum(x, y) result(total)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: total
      integer :: k

      total = 0.0d0
      do k = 1, size(x)
         total = total + x(k)*y(k)
      end do
end function weighted_sum

! ---------------------------------------------------------------
! The masses after a loss or gain dM = fraction M.
subroutine check_shell_masses(jenv, fraction)
      integer, intent(in) :: jenv
      real(dp), intent(in) :: fraction
      real(dp) :: m_old(nz), delta_mass, cz_mass, scale, err, new_total
      integer :: k

      call build_model(jenv, 0.0d0)
      write(*,'(a,i0,a,i0,a,1pe9.1)') '# mass: zone ', jenv, '..', nz, &
           ', dM/M =', fraction
      m_old = m
      delta_mass = fraction*total_mass
      new_total = total_mass + delta_mass
      if (jenv > 1) then
         cz_mass = m_old(nz) - m_old(jenv)
      else
         cz_mass = m_old(nz)
      end if
      if (delta_mass < 0.0d0) then
         scale = (cz_mass + delta_mass)/cz_mass
      else
         scale = 1.0d0
      end if
      call report('M_cz = convection_zone_mass', &
           abs(convection_zone_mass(jenv, m) - cz_mass)/cz_mass)
      call change_shell_masses(jenv, delta_mass, new_total, m, log_mass, dm)
      call report('sum of dm = M + dM', &
           abs(sum(dm) - new_total)/new_total)
      err = 0.0d0
      do k = jenv+1, nz
         err = max(err, abs((m(k) - m_old(jenv)) - scale*(m_old(k) - m_old(jenv)))/ &
              (scale*(m_old(k) - m_old(jenv))))
      end do
      if (delta_mass < 0.0d0) then
         call report('m(k > jenv) - m(jenv) scale by (M_cz+dM)/M_cz', err)
      else
         call report('m(k > jenv) unchanged (dM lands in dm(nz))', err)
      end if
      err = 0.0d0
      do k = 1, jenv
         err = max(err, abs(m(k) - m_old(k))/m_old(k))
      end do
      call report('m(k <= jenv) unchanged', err)
      if (jenv == 1) call report('dm(1) = (m(1) + m(2))/2', &
           abs(dm(1) - 0.5d0*(m(1) + m(2)))/(0.5d0*(m(1) + m(2))))
end subroutine check_shell_masses

! ---------------------------------------------------------------
! The angular momentum after a loss of 1e-3 M from the zone rotating as
! omega ~ r^walpcz.
subroutine check_wind_angular_momentum(jenv, walpcz)
      integer, intent(in) :: jenv
      real(dp), intent(in) :: walpcz
      real(dp) :: dm_before(nz-jenv+1), delta_mass, j_surface, &
           j_surface_expected, cz_am, delta_am, cz_am_after
      type(wind_angular_momentum) :: am

      call build_model(jenv, walpcz)
      write(*,'(a,i0,a,i0,a,f5.1,a,1pe9.1)') '# angular momentum: zone ', &
           jenv, '..', nz, ', walpcz =', walpcz, ', dM/M =', -mass_change_fraction
      delta_mass = -mass_change_fraction*total_mass
      j_surface_expected = 2.0d0/3.0d0*total_radius**2*surface_omega* &
           (total_radius/radius(nz))**nint(walpcz)
      delta_am = j_surface_expected*delta_mass
      dm_before = dm(jenv:nz)
      cz_am = weighted_sum(j_rot(jenv:nz), dm_before)
      j_surface = wind_specific_angular_momentum(walpcz, total_radius, &
           omega(nz), log_radius(nz), j_rot(nz))
      call change_shell_masses(jenv, delta_mass, total_mass + delta_mass, &
           m, log_mass, dm)
      call remove_wind_angular_momentum(j_surface, delta_mass, dm_before, &
           dm(jenv:nz), j_rot(jenv:nz), am)
      call report('j_s = (2/3) R^2 omega(nz) (R/r(nz))^walpcz', &
           abs(j_surface - j_surface_expected)/j_surface_expected)
      call report('dJ = j_s dM', abs(am%lost - delta_am)/abs(delta_am))
      call report('J_cz = sum of j dm before', abs(am%zone_before - cz_am)/cz_am)
      cz_am_after = weighted_sum(j_rot(jenv:nz), dm(jenv:nz))
      call report('sum of j dm after = J_cz + dJ', &
           abs(cz_am_after - (cz_am + delta_am))/(cz_am + delta_am))
end subroutine check_wind_angular_momentum

! ---------------------------------------------------------------
! The composition after accreting 1e-3 M into the zone.
subroutine check_accreted_composition(jenv)
      integer, intent(in) :: jenv
      real(dp) :: dm_before(nz-jenv+1), xa_before(n_species_extended,nz-jenv+1), &
           delta_mass, cz_mass, expected
      integer :: s

      call build_model(jenv, 0.0d0)
      write(*,'(a,i0,a,i0,a,1pe9.1)') '# composition: zone ', jenv, '..', nz, &
           ', dM/M =', mass_change_fraction
      delta_mass = mass_change_fraction*total_mass
      dm_before = dm(jenv:nz)
      xa_before = xa(:,jenv:nz)
      cz_mass = sum(dm_before)
      call change_shell_masses(jenv, delta_mass, total_mass + delta_mass, &
           m, log_mass, dm)
      call mix_accreted_composition(cz_mass, delta_mass, x_accreted, &
           xa(:,jenv:nz))
      do s = 1, n_species_basic
         expected = weighted_sum(xa_before(s,:), dm_before) + x_accreted(s)*delta_mass
         call report('sum of X dm = old + X_acc dM, '//species_name(s), &
              abs(weighted_sum(xa(s,jenv:nz), dm(jenv:nz)) - expected)/expected)
      end do
end subroutine check_accreted_composition

end program test_mass_change
