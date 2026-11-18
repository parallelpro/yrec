!----------------------------------------------------------------------
! test_composition_kernels
!----------------------------------------------------------------------
! The star-blind kernels of the composition step
! (core/step/composition_kernels.f90) on synthetic models of a few
! shells; no star, no tables.
!   deep mixing  the shells mixed out to the mass limit, the convective
!                flag beyond it;
!   homogenize   each species keeps its mass in a mixed region, which
!                ends uniform; the shells and slots outside are
!                untouched; a one-shell region away from the centre is
!                skipped;
!   closure      every abundance in [0, 1], Z <= 1 - X - He3, and
!                X + Y + Z + He3 = 1.
! Each check prints its error (or its value) and ok or FAIL; any FAIL
! stops the program with status 1.
program test_composition_kernels
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: n_species_extended, n_species_basic, &
           i_h1, i_he3, i_he4, i_metals
      use composition_kernels, only: deep_mixing_flags, homogenize_regions, &
           clip_and_close_abundances
      implicit none

      integer :: nfail

      nfail = 0
      write(*,'(a)') 'test_composition_kernels report'
      call check_deep_mixing()
      call check_homogenize()
      call check_closure()
      if (nfail > 0) then
         write(*,'(a,i0,a)') 'test_composition_kernels: ', nfail, ' FAIL'
         stop 1
      end if
      write(*,'(a)') 'test_composition_kernels: done'

contains

subroutine report(name, value, pass)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: value
      logical, intent(in) :: pass
      if (pass) then
         write(*,'(a,t68,1pe9.1,2x,a)') name, value, 'ok'
      else
         write(*,'(a,t68,1pe9.1,2x,a)') name, value, 'FAIL'
         nfail = nfail + 1
      end if
end subroutine report

! The relative error of got against want.
real(dp) function relerr(got, want)
      real(dp), intent(in) :: got, want
      relerr = abs(got - want)/abs(want)
end function relerr

subroutine check_deep_mixing()
      logical :: conv(6), flag(6)
      real(dp) :: mass(6)
      write(*,'(a)') '-- deep mixing'
      mass = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
      conv = [.false., .false., .false., .false., .true., .false.]
      flag = conv
      call deep_mixing_flags(conv, mass, 3.5_dp, 6, flag)
      call report('deep mixing: mixed out to the mass limit, the flag beyond', 0.0_dp, &
           all(flag .eqv. [.true., .true., .true., .false., .true., .false.]))
      flag = conv
      call deep_mixing_flags(conv, mass, 6.0_dp, 6, flag)
      call report('deep mixing: a shell at the limit is mixed', 0.0_dp, all(flag))
end subroutine check_deep_mixing

subroutine check_homogenize()
      integer, parameter :: nh = 8, ns = 5
      real(dp) :: x(ns+1,nh), xb(ns+1,nh), w(nh), before, after, err, lone
      integer :: bounds(3,2), k, i
      write(*,'(a)') '-- homogenize'
      do k = 1, nh
         w(k) = 1.0_dp + 0.37_dp*k
         do i = 1, ns+1
            x(i,k) = 0.01_dp*i + 0.003_dp*k*k
         end do
      end do
! a value that the mean of one shell, (x w)/w, does not give back
      do i = 1, 200
         lone = 0.1_dp + 0.001_dp*i
         if ((lone*w(6))/w(6) /= lone) exit
      end do
      x(1,6) = lone
      xb = x
      bounds(1,:) = [1, 1]
      bounds(2,:) = [2, 4]
      bounds(3,:) = [6, 6]
      call homogenize_regions(x, w, bounds, 3, ns)
      err = 0.0_dp
      do i = 1, ns
         before = sum(xb(i,2:4)*w(2:4))
         after = sum(x(i,2:4)*w(2:4))
         err = max(err, relerr(after, before))
      end do
      call report('homogenize: each species keeps its mass in a region', err, err <= 1.0e-14_dp)
      call report('homogenize: the region is uniform', maxval(abs(x(1:ns,2:4) - spread(x(1:ns,2),2,3))), &
           all(x(1:ns,2:4) == spread(x(1:ns,2),2,3)) .and. x(1,2) /= xb(1,2))
      call report('homogenize: the shells and slots outside are untouched', 0.0_dp, &
           all(x(:,5) == xb(:,5)) .and. all(x(:,7:8) == xb(:,7:8)) .and. &
           all(x(ns+1,:) == xb(ns+1,:)))
      call report('homogenize: a one-shell region away from the centre is skipped', &
           abs(x(1,6) - lone), x(1,6) == lone .and. (lone*w(6))/w(6) /= lone)
end subroutine check_homogenize

subroutine check_closure()
      real(dp) :: x(n_species_extended,3), total(3)
      write(*,'(a)') '-- closure'
      x = 0.2_dp
! shell 1: values outside [0, 1]; shell 2: Z above 1 - X - He3; shell 3:
! an ordinary mixture whose Y is not the closure's
      x(:,1) = [1.3_dp, -0.2_dp, -0.1_dp, -0.3_dp, 1.7_dp, x(6:,1)]
      x(i_h1,2) = 0.7_dp
      x(i_he3,2) = 0.1_dp
      x(i_metals,2) = 0.25_dp
      x(i_h1,3) = 0.71_dp
      x(i_he3,3) = 1.0e-4_dp
      x(i_metals,3) = 0.02_dp
      x(i_he4,3) = 0.5_dp
      call clip_and_close_abundances(x, 3, n_species_basic)
      call report('closure: every abundance in [0, 1]', maxval(x(1:n_species_basic,:)), &
           all(x(1:n_species_basic,:) >= 0.0_dp) .and. all(x(1:n_species_basic,:) <= 1.0_dp) &
           .and. x(i_h1,1) == 1.0_dp .and. x(5,1) == 1.0_dp)
      call report('closure: Z <= 1 - X - He3', x(i_metals,2), &
           x(i_metals,2) == 1.0_dp - x(i_h1,2) - x(i_he3,2))
      total = x(i_h1,:) + x(i_he4,:) + x(i_metals,:) + x(i_he3,:)
      call report('closure: X + Y + Z + He3 = 1', maxval(abs(total - 1.0_dp)), &
           maxval(abs(total - 1.0_dp)) <= 4.0e-16_dp)
end subroutine check_closure

end program test_composition_kernels
