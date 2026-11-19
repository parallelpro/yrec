!----------------------------------------------------------------------
! test_li_be_burn
!----------------------------------------------------------------------
! The star-blind kernels of the Li6, Li7 and Be9 burn
! (core/step/li_be_burn.f90) on a synthetic model of twelve shells; no
! star, no tables.
!   base         the base of the surface mixed region
!                (convection_geometry's surface_mixing_base):
!                without overshoot the convection zone's base zones and
!                the edge where del_ad = del_rad, clipped to half a
!                zone; a base at the outermost shell passed through;
!                with overshoot, a step's first call from the
!                start-of-step model and a later call from the last
!                call's zone, the end base the overshoot depth below
!                the interpolated edge, and what the call leaves for
!                the next;
!   average      the mixed region's mean abundances, its mass with the
!                base zone's edge, and the log of its mean rate;
!   radiative    shells below the mixed region against the closed-form
!                survival X exp(-Lambda), for a constant rate and for a
!                rate log-linear in time; the cleared, too-cool, depleted
!                and floored shells;
!   mixed        the mixed region at its mean times exp(-Lambda),
!                conserving each species when nothing burns, and left
!                alone when its base is too cool;
!   accretion    the accreted share: conserved without burning, burning
!                from its arrival on, none without inflow;
!   moving base  the shells a retreating region leaves, and a region
!                that deepened.
! Each check prints its error (or its value) and ok or FAIL; any FAIL
! stops the program with status 1.
program test_li_be_burn
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: n_species_extended
      use convection_geometry, only: surface_mixing_base
      use li_be_burn, only: mixed_region_bases, li_be_accretion, &
           average_over_mixed_region, burn_li_be_zones, li_be_species, &
           num_li_be
      implicit none

      integer, parameter :: n = 12
      real(dp), parameter :: dt = 1.0e10_dp
! the rates (per s) at the start and end of the step, Li6, Li7, Be9
      real(dp), parameter :: lam0(num_li_be) = [2.0e-10_dp, 1.0e-10_dp, 0.5e-10_dp], &
           lam1(num_li_be) = [4.0e-10_dp, 3.0e-10_dp, 1.0e-10_dp]
      real(dp) :: xa(n_species_extended,n), x0(n_species_extended,n)
      real(dp) :: m(n), dm(n), logt(n), rate0(num_li_be,n), rate1(num_li_be,n)
      type(mixed_region_bases) :: base
      type(li_be_accretion) :: no_accretion, accretion
      integer :: nfail

      nfail = 0
      write(*,'(a)') 'test_li_be_burn report'
      call check_base()
      call check_average()
      call check_radiative()
      call check_mixed_region()
      call check_accretion()
      call check_moving_base()
      if (nfail > 0) then
         write(*,'(a,i0,a)') 'test_li_be_burn: ', nfail, ' FAIL'
         stop 1
      end if
      write(*,'(a)') 'test_li_be_burn: done'

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

! exp(-Lambda) for a rate that goes from r0 to r1 linearly in its log
! over the step: Lambda = dt (r1 - r0)/ln(r1/r0).
real(dp) function log_linear_survival(r0, r1)
      real(dp), intent(in) :: r0, r1
      log_linear_survival = exp(-dt*(r1 - r0)/log(r1/r0))
end function log_linear_survival

! The synthetic model: shell masses that grow outward, the zone centres'
! masses, a graded light-element composition, uniform rates and cool
! shells.
subroutine reset_model()
      integer :: k, isp
      do k = 1, n
         dm(k) = (1.0_dp + 0.1_dp*k)*1.0e30_dp
      end do
      m(1) = 0.5_dp*dm(1)
      do k = 2, n
         m(k) = m(k-1) + 0.5_dp*(dm(k-1) + dm(k))
      end do
      xa = 0.0_dp
      do k = 1, n
         do isp = 1, num_li_be
            xa(li_be_species(isp),k) = 1.0e-9_dp*(isp + 0.1_dp*k)
            rate0(isp,k) = lam0(isp)
            rate1(isp,k) = lam1(isp)
         end do
      end do
      x0 = xa
      logt = 6.5_dp
      base = mixed_region_bases(zone_start=9, zone_end=9)
end subroutine reset_model

! The mass of li_be species isp in shells first..n.
real(dp) function species_mass(x, isp, first)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: isp, first
      species_mass = sum(x(li_be_species(isp),first:n)*dm(first:n))
end function species_mass

! The base of the surface mixed region on a model whose shells sit at
! r(k) = k 1e10 cm at the end of the step and 2% further out at its
! start, with del_ad - del_rad = -0.1 in the convection zone's base
! zone, 9, and +0.3 in the zone below: the edge is 3/4 of the way up
! from zone 8, an offset of -0.25 from the zones' midpoint.
subroutine check_base()
      real(dp) :: logr(n), logr0(n), grada(n), gradr(n), frac, want, last_radius
      integer :: k, zone, zone_old, last_zone
      write(*,'(a)') '-- the base of the mixed region'
      do k = 1, n
         logr(k) = log10(k*1.0e10_dp)
         logr0(k) = log10(k*1.02e10_dp)
      end do
      grada = 0.4_dp
      gradr = 0.1_dp
      gradr(9:n) = 0.5_dp
      want = 0.5_dp - 0.3_dp/(0.3_dp - (-0.1_dp))
      frac = 0.0_dp; last_radius = 0.0_dp; last_zone = 0
      call surface_mixing_base(logr, logr0, grada, gradr, 9, 8, n, .false., &
           1.5e10_dp, 2.2e10_dp, last_radius, last_zone, zone, zone_old, frac)
      call report('base: the zones and the edge where del_ad = del_rad', abs(frac - want), &
           zone == 9 .and. zone_old == 8 .and. abs(frac - want) <= 1.0e-15_dp .and. last_zone == 9)
      gradr(8) = 0.5_dp
      gradr(9) = 0.7_dp
      call surface_mixing_base(logr, logr0, grada, gradr, 9, 8, n, .false., &
           1.5e10_dp, 2.2e10_dp, last_radius, last_zone, zone, zone_old, frac)
      call report('base: the edge clipped to half a zone', frac, frac == 0.5_dp)
      gradr(8) = 0.1_dp
      gradr(9) = 0.5_dp
      frac = 0.123_dp
      call surface_mixing_base(logr, logr0, grada, gradr, n, 8, n, .false., &
           1.5e10_dp, 2.2e10_dp, last_radius, last_zone, zone, zone_old, frac)
      call report('base: at the outermost shell, the zones passed through', frac, &
           zone == n .and. zone_old == 8 .and. frac == 0.123_dp)
! With overshoot, the first call of a step: the start base 1.5e10 cm
! below the midpoint of the start model's zones 7 and 8 (7.65e10 cm),
! the end base 2.2e10 cm below the edge at 8.75e10 cm.
      last_radius = 0.0_dp; last_zone = 0
      call surface_mixing_base(logr, logr0, grada, gradr, 9, 8, n, .true., &
           1.5e10_dp, 2.2e10_dp, last_radius, last_zone, zone, zone_old, frac)
      call report('base, overshoot: the start base from the start-of-step model', real(zone_old, dp), &
           zone_old == 7)
      call report('base, overshoot: the end base below the interpolated edge', abs(frac + 0.05_dp), &
           zone == 7 .and. abs(frac + 0.05_dp) <= 1.0e-12_dp)
      call report('base, overshoot: the edge and the zone kept for the next call', &
           abs(last_radius/8.75e10_dp - 1.0_dp), abs(last_radius/8.75e10_dp - 1.0_dp) <= 1.0e-14_dp &
           .and. last_zone == 7)
      last_zone = 5
      call surface_mixing_base(logr, logr0, grada, gradr, 9, 8, n, .true., &
           1.5e10_dp, 2.2e10_dp, last_radius, last_zone, zone, zone_old, frac)
      call report('base, overshoot: a later call starts from the last call''s zone', real(zone_old, dp), &
           zone_old == 5)
      last_radius = 0.0_dp
      call surface_mixing_base(logr, logr0, grada, gradr, 9, 1, n, .true., &
           1.5e10_dp, 2.2e10_dp, last_radius, last_zone, zone, zone_old, frac)
      call report('base, overshoot: fully convective at the start, zone 1', real(zone_old, dp), &
           zone_old == 1)
end subroutine check_base

subroutine check_average()
      real(dp) :: x_mean(num_li_be), log_rate(num_li_be), region_mass, want, err, &
           edge_mass, edge_rate
      integer :: isp, k
      write(*,'(a)') '-- average'
      call reset_model()
      do k = 1, n
         rate0(:,k) = lam0*(1.0_dp + 0.1_dp*k)
      end do
      call average_over_mixed_region(xa, m, dm, rate0, 9, 0.0_dp, n, x_mean, &
           log_rate, region_mass)
      err = 0.0_dp
      do isp = 1, num_li_be
         err = max(err, relerr(x_mean(isp), species_mass(xa, isp, 9)/sum(dm(9:n))))
      end do
      call report('average: the mass-weighted mean abundance', err, err <= 1.0e-15_dp)
      call average_over_mixed_region(xa, m, dm, rate0, 9, 0.25_dp, n, x_mean, &
           log_rate, region_mass)
      edge_mass = dm(9) + 0.25_dp*(m(9) - m(8))
      want = edge_mass + sum(dm(10:n))
      call report('average: the region''s mass with the base zone''s edge', &
           relerr(region_mass, want), relerr(region_mass, want) <= 1.0e-15_dp)
      err = 0.0_dp
      do isp = 1, num_li_be
         edge_rate = rate0(isp,9) + 0.5_dp*0.25_dp*(rate0(isp,8) - rate0(isp,9))
         want = log((edge_rate*edge_mass + sum(rate0(isp,10:n)*dm(10:n)))/region_mass)
         err = max(err, relerr(log_rate(isp), want))
      end do
      call report('average: ln of the mean rate, the base zone''s moved to the edge', &
           err, err <= 1.0e-14_dp)
      call average_over_mixed_region(xa, m, dm, rate0, n, 0.3_dp, n, x_mean, &
           log_rate, region_mass)
      call report('average: no edge for a base at the outermost shell', &
           relerr(region_mass, dm(n)), region_mass == dm(n))
end subroutine check_average

subroutine check_radiative()
      real(dp) :: err
      integer :: isp, k
      write(*,'(a)') '-- radiative shells'
      call reset_model()
      call burn_li_be_zones(dt, xa, m, dm, logt, rate0, rate1, base, no_accretion, n)
      err = 0.0_dp
      do k = 1, 8
         do isp = 1, num_li_be
            err = max(err, relerr(xa(li_be_species(isp),k), &
                 x0(li_be_species(isp),k)*log_linear_survival(lam0(isp), lam1(isp))))
         end do
      end do
      call report('li/be radiative: X exp(-Lambda), the rate log-linear in time', err, &
           err <= 1.0e-13_dp)
      call reset_model()
      call burn_li_be_zones(dt, xa, m, dm, logt, rate0, rate0, base, no_accretion, n)
      err = 0.0_dp
      do k = 1, 8
         do isp = 1, num_li_be
            err = max(err, relerr(xa(li_be_species(isp),k), &
                 x0(li_be_species(isp),k)*exp(-lam0(isp)*dt)))
         end do
      end do
      call report('li/be radiative: X exp(-lambda dt) at a constant rate', err, &
           err <= 1.0e-14_dp)
      call reset_model()
      logt(3) = 7.1_dp
      call burn_li_be_zones(dt, xa, m, dm, logt, rate0, rate1, base, no_accretion, n)
      call report('li/be radiative: a shell hotter than log T = 7 is cleared', &
           maxval(xa(li_be_species,3)), all(xa(li_be_species,3) == 0.0_dp) .and. &
           all(xa(li_be_species,2) > 0.0_dp) .and. all(xa(li_be_species,4) > 0.0_dp))
      call reset_model()
      rate1(3,5) = 0.0_dp
      call burn_li_be_zones(dt, xa, m, dm, logt, rate0, rate1, base, no_accretion, n)
      call report('li/be radiative: the burn stops at the first shell too cool', 0.0_dp, &
           all(xa(li_be_species,5:8) == x0(li_be_species,5:8)) .and. &
           all(xa(li_be_species,1:4) < x0(li_be_species,1:4)))
      call reset_model()
      xa(li_be_species(1),2) = 1.0e-23_dp
      call burn_li_be_zones(dt, xa, m, dm, logt, rate0, rate1, base, no_accretion, n)
      call report('li/be radiative: an abundance below 1e-24 becomes zero', &
           xa(li_be_species(1),2), xa(li_be_species(1),2) == 0.0_dp .and. &
           xa(li_be_species(2),2) > 0.0_dp)
      call reset_model()
      xa(li_be_species,4) = 5.0e-25_dp
      call burn_li_be_zones(dt, xa, m, dm, logt, rate0, rate1, base, no_accretion, n)
      call report('li/be radiative: a depleted shell is skipped', &
           xa(li_be_species(1),4), all(xa(li_be_species,4) == 5.0e-25_dp))
end subroutine check_radiative

subroutine check_mixed_region()
      real(dp) :: err, want
      integer :: isp
      write(*,'(a)') '-- the mixed region'
      call reset_model()
      call burn_li_be_zones(dt, xa, m, dm, logt, rate0, rate1, base, no_accretion, n)
      err = 0.0_dp
      do isp = 1, num_li_be
         want = species_mass(x0, isp, 9)/sum(dm(9:n))*log_linear_survival(lam0(isp), lam1(isp))
         err = max(err, maxval(abs(xa(li_be_species(isp),9:n) - want))/want)
      end do
      call report('li/be mixed: uniform at the mean times exp(-Lambda)', err, &
           err <= 1.0e-13_dp)
! rates too small to burn anything over the step: the mixing alone
      call reset_model()
      rate0 = 1.0e-30_dp
      rate1 = 1.0e-30_dp
      call burn_li_be_zones(dt, xa, m, dm, logt, rate0, rate1, base, no_accretion, n)
      err = 0.0_dp
      do isp = 1, num_li_be
         err = max(err, relerr(species_mass(xa, isp, 9), species_mass(x0, isp, 9)))
      end do
      call report('li/be mixed: the mixing conserves each species', err, &
           err <= 1.0e-14_dp .and. xa(li_be_species(1),9) /= x0(li_be_species(1),9))
      call reset_model()
      rate0(3,9) = 0.0_dp
      call burn_li_be_zones(dt, xa, m, dm, logt, rate0, rate1, base, no_accretion, n)
      call report('li/be mixed: nothing changes when its base is too cool', 0.0_dp, &
           all(xa(li_be_species,9:n) == x0(li_be_species,9:n)))
end subroutine check_mixed_region

subroutine check_accretion()
      real(dp) :: err, cz_mass, want, survive, xb(n_species_extended,n)
      integer :: isp
      write(*,'(a)') '-- accretion'
      call reset_model()
      cz_mass = sum(dm(9:n))
      accretion%inflow = .true.
      accretion%mass = 0.3_dp*cz_mass
      accretion%abundance = [3.0e-9_dp, 2.0e-9_dp, 1.0e-9_dp]
      rate0 = 1.0e-30_dp
      rate1 = 1.0e-30_dp
      call burn_li_be_zones(dt, xa, m, dm, logt, rate0, rate1, base, accretion, n)
      err = 0.0_dp
      do isp = 1, num_li_be
         want = species_mass(x0, isp, 9) + accretion%abundance(isp)*accretion%mass
         err = max(err, relerr(xa(li_be_species(isp),n)*(cz_mass + accretion%mass), want))
      end do
      call report('li/be accretion: the zone gains the accreted share, x_acc M_acc', err, &
           err <= 1.0e-14_dp)
! bug_sweep M383: a region too cool to burn (rates at zero) gains it too.
      call reset_model()
      rate0 = 0.0_dp
      rate1 = 0.0_dp
      call burn_li_be_zones(dt, xa, m, dm, logt, rate0, rate1, base, accretion, n)
      err = 0.0_dp
      do isp = 1, num_li_be
         want = species_mass(x0, isp, 9) + accretion%abundance(isp)*accretion%mass
         err = max(err, relerr(xa(li_be_species(isp),n)*(cz_mass + accretion%mass), want))
      end do
      call report('li/be accretion: too cool to burn, the zone gains it (M383)', err, &
           err <= 1.0e-14_dp)
! a constant rate with lambda dt = 2: the gas that arrives uniformly
! over the step keeps (1 - exp(-2))/2 of its Li6 on average
      call reset_model()
      rate0 = 2.0e-10_dp
      call burn_li_be_zones(dt, xa, m, dm, logt, rate0, rate0, base, accretion, n)
      survive = (1.0_dp - exp(-2.0_dp))/2.0_dp
      want = (species_mass(x0, 1, 9)*exp(-2.0_dp) + &
           accretion%abundance(1)*survive*accretion%mass)/(cz_mass + accretion%mass)
      err = relerr(xa(li_be_species(1),n), want)
      call report('li/be accretion: the accreted gas burns from its arrival on', err, &
           err <= 1.0e-6_dp)
      call reset_model()
      rate0 = 2.0e-10_dp
      call burn_li_be_zones(dt, xa, m, dm, logt, rate0, rate0, base, no_accretion, n)
      accretion%inflow = .false.
      xb = xa
      call reset_model()
      rate0 = 2.0e-10_dp
      call burn_li_be_zones(dt, xa, m, dm, logt, rate0, rate0, base, accretion, n)
      call report('li/be accretion: none without inflow', 0.0_dp, &
           all(xa(li_be_species,:) == xb(li_be_species,:)))
end subroutine check_accretion

subroutine check_moving_base()
      real(dp) :: err, want, m_beg, m_end, frac, rate
      integer :: isp
      write(*,'(a)') '-- a base that moves'
! The region retreats: its base moves out from shell 7 to shell 10.
      call reset_model()
      base = mixed_region_bases(zone_start=7, zone_end=10)
      call burn_li_be_zones(dt, xa, m, dm, logt, rate0, rate1, base, no_accretion, n)
      m_beg = 0.5_dp*(m(7) + m(6))
      m_end = 0.5_dp*(m(10) + m(9))
      frac = (m_end - m(8))/(m_end - m_beg)
      err = 0.0_dp
      do isp = 1, num_li_be
         rate = exp(frac*log(lam1(isp)) + (1.0_dp - frac)*log(lam0(isp)))
         err = max(err, relerr(xa(li_be_species(isp),8), x0(li_be_species(isp),8)*exp(-rate*dt)))
      end do
      call report('li/be retreat: a shell left burns on its time-weighted rate', err, &
           err <= 1.0e-13_dp)
      err = 0.0_dp
      do isp = 1, num_li_be
         want = species_mass(x0, isp, 7)/sum(dm(7:n))*log_linear_survival(lam0(isp), lam1(isp))
         err = max(err, maxval(abs(xa(li_be_species(isp),10:n) - want))/want)
      end do
      call report('li/be retreat: the region keeps its start mean times exp(-Lambda)', err, &
           err <= 1.0e-13_dp)
! The region deepens: its base moves in from shell 10 to shell 7.
      call reset_model()
      base = mixed_region_bases(zone_start=10, zone_end=7)
      call burn_li_be_zones(dt, xa, m, dm, logt, rate0, rate1, base, no_accretion, n)
      err = 0.0_dp
      do isp = 1, num_li_be
         want = species_mass(x0, isp, 7)/sum(dm(7:n))*log_linear_survival(lam0(isp), lam1(isp))
         err = max(err, maxval(abs(xa(li_be_species(isp),10:n) - want))/want)
      end do
      call report('li/be deepening: starts from the mean over the deeper region', err, &
           err <= 1.0e-13_dp)
      call report('li/be deepening: the shells between the bases are left', 0.0_dp, &
           all(xa(li_be_species,7:9) == x0(li_be_species,7:9)))
end subroutine check_moving_base

end program test_li_be_burn
