!----------------------------------------------------------------------
! test_element_diffusion
!----------------------------------------------------------------------
! Standalone test of element diffusion, printed for byte-comparison
! against expected_test_element_diffusion.out. Both settling routines run one
! step on a synthetic 200-zone solar-like model (radiative interior,
! convective envelope above 0.71 R):
!   gravitational_settling  the original routines (use_new_diffusion_
!                           routines false): helium only, helium and
!                           metals, each Coulomb-logarithm choice, the
!                           constant mixing and settling factors, the
!                           helium and metal scale factors, the
!                           Thoul et al. fit and the older analytic
!                           coefficients instead of the Thoul et al.
!                           ones, and a solve stopped short of
!                           convergence
!   microdiff               the newer routines (Bahcall & Loeb 1990
!                           units): helium only, helium and metals, the
!                           light elements (they ignore use_thoul_fit and
!                           use_thoul_diffusion)
! and, for both, the abundance floors: the outer edge of settling where
! Y first exceeds helium_diffusion_min, and no settling when Y or X is
! below its floor everywhere.
! Each case prints the change of H, He4, Z, Li7 and Be9 at nine zones
! and the mass-weighted change of each species over the star, all at
! full precision. No deck pins the original routines.
program test_element_diffusion
      use star_info_lib, only: star, json, i_h1, i_he4, i_metals, i_he3, &
           i_c12, i_c13, i_n14, i_n15, i_o16, i_o17, i_o18, i_h2, i_li6, &
           i_li7, i_be9
      use io_units, only: run_log_unit, terminal_unit
      use element_diffusion_lib, only: diffuse_elements
      use element_diffusion_scratch, only: diff_scr
      implicit none
      integer, parameter :: nz = 200
      double precision, parameter :: msun = 1.989d33, rsun = 6.9598d10, &
           year = 3.1558d7, cgs_g = 6.6726d-8, k_over_mu_mh = 1.3806d-16/(0.62d0*1.6605d-24)
      double precision :: xa0(15,json), xa(15,json), dlnp_dr(json), log_r(json), &
           log_rho(json), m(json), log_t(json), del(json), dm(json)
      logical :: convective(json)
      double precision :: dt, total_mass
      integer :: i

      run_log_unit = 6
      terminal_unit = 6
      call build_model()

      call reset_controls()
      star%job%diffuse_helium_active = .true.
      call settle_case('old: helium', .false.)
      star%job%use_diffusion_z = .true.
      call settle_case('old: helium and metals', .false.)
      star%ctrl%coulomb_log_choice = 2
      call settle_case('old: Coulomb logarithm choice 2', .false.)
      star%ctrl%coulomb_log_choice = 3
      call settle_case('old: Coulomb logarithm choice 3', .false.)
      star%ctrl%coulomb_log_choice = 4
      call settle_case('old: Coulomb logarithm choice 4', .false.)
      star%ctrl%coulomb_log_choice = 1
      star%ctrl%constant_mixing_coeff = 0.5d0
      star%ctrl%constant_settling_reduction = 0.8d0
      call settle_case('old: constant mixing 0.5, settling reduction 0.8', .false.)
      star%ctrl%constant_mixing_coeff = 1.0d0
      star%ctrl%constant_settling_reduction = 1.0d0
      star%ctrl%fgry = 0.7d0
      star%job%fgrz = 0.9d0
      call settle_case('old: helium factor 0.7, metal factor 0.9', .false.)
      star%ctrl%fgry = 1.0d0
      star%job%fgrz = 1.0d0
      star%ctrl%use_thoul_fit = .true.
      call settle_case('old: the Thoul et al. fit for metals', .false.)
      star%ctrl%use_thoul_fit = .false.
      star%ctrl%use_thoul_diffusion = .false.
      call settle_case('old: without the Thoul et al. coefficients', .false.)
      star%ctrl%use_thoul_diffusion = .true.
      star%ctrl%settling_num_iterations = 1
      call settle_case('old: one iteration, short of convergence', .false.)
      star%ctrl%settling_num_iterations = 10
      call floors(.false.)

      call reset_controls()
      star%job%diffuse_helium_active = .true.
      call settle_case('new: helium', .true.)
      star%job%use_diffusion_z = .true.
      call settle_case('new: helium and metals', .true.)
      star%ctrl%diffuse_lithium = .true.
      star%ctrl%coulomb_log_choice = 4
      call settle_case('new: helium, metals and light elements', .true.)
      call floors(.true.)

contains

! A solar-like model: center (zone 1) to surface (zone nz), the mass
! and pressure gradient consistent with an ideal gas in hydrostatic
! equilibrium, the composition with a burned core and the light
! elements present only outside it.
subroutine build_model()
      double precision :: f, r(json), rho(json), p
      xa0 = 0.0d0
      do i = 1, nz
         f = dble(i - 1)/dble(nz - 1)
         r(i) = rsun*(0.002d0 + 0.998d0*f)
         log_r(i) = log10(r(i))
         log_t(i) = 7.196d0 - 0.85d0*f**1.3d0 - 2.6d0*f**25
         log_rho(i) = 2.18d0 - 3.8d0*f**1.1d0 - 5.0d0*f**30
         rho(i) = 10.0d0**log_rho(i)
         convective(i) = f > 0.71d0
         del(i) = merge(0.40d0, 0.25d0 + 0.1d0*f, convective(i))
         xa0(i_h1,i) = 0.71d0 - 0.35d0*exp(-(f/0.12d0)**2)
         xa0(i_metals,i) = 0.018d0
         xa0(i_he3,i) = 1.0d-5*(1.0d0 + 4.0d0*f)
         xa0(i_c12,i) = 3.0d-3
         xa0(i_c13,i) = 4.0d-5
         xa0(i_n14,i) = 1.0d-3
         xa0(i_n15,i) = 4.0d-6
         xa0(i_o16,i) = 9.0d-3
         xa0(i_o17,i) = 4.0d-6
         xa0(i_o18,i) = 2.0d-5
         if (f > 0.6d0) then
            xa0(i_h2,i) = 2.8d-5
            xa0(i_li6,i) = 7.0d-10
            xa0(i_li7,i) = 1.0d-8
            xa0(i_be9,i) = 1.6d-10
         end if
         xa0(i_he4,i) = 1.0d0 - xa0(i_h1,i) - xa0(i_metals,i) - xa0(i_he3,i)
      end do
! the mass from the density (trapezoid rule), then the pressure
! gradient of an ideal gas in hydrostatic equilibrium
      m(1) = 4.0d0/3.0d0*acos(-1.0d0)*r(1)**3*rho(1)
      do i = 2, nz
         m(i) = m(i-1) + 2.0d0*acos(-1.0d0)*(r(i)**2*rho(i) + r(i-1)**2*rho(i-1))*(r(i) - r(i-1))
      end do
      total_mass = m(nz)
      do i = 1, nz
         p = rho(i)*k_over_mu_mh*10.0d0**log_t(i)
         dlnp_dr(i) = -cgs_g*m(i)*rho(i)/(p*r(i)*r(i))
      end do
      dm(1) = m(1)
      do i = 2, nz
         dm(i) = m(i) - m(i-1)
      end do
      star%gradT(1:nz) = del(1:nz)
      dt = 1.0d8*year
end subroutine build_model

! The settling controls at their defaults (the tolerance 1e-7: the original
! solve stalls near 4e-8 at this model's centre), the switches off.
subroutine reset_controls()
      star%job%diffuse_helium_active = .false.
      star%job%use_diffusion_z = .false.
      star%ctrl%diffuse_lithium = .false.
      star%ctrl%use_thoul_fit = .false.
      star%ctrl%use_thoul_diffusion = .true.
      star%ctrl%coulomb_log_choice = 1
      star%ctrl%settling_num_iterations = 10
      star%ctrl%settling_tolerance = 1.0d-7
      star%ctrl%hydrogen_diffusion_floor = 1.0d-3
      star%ctrl%helium_diffusion_min = 1.0d-3
      star%ctrl%constant_mixing_coeff = 1.0d0
      star%ctrl%constant_settling_reduction = 1.0d0
      star%ctrl%fgry = 1.0d0
      star%job%fgrz = 1.0d0
      diff_scr%settling_suspended_reported = .false.
end subroutine reset_controls

! The abundance floors: settling only below the outermost zone with
! Y > helium_diffusion_min, and none at all (a run-log note) when Y or X
! is below its floor everywhere.
subroutine floors(new_routines)
      logical, intent(in) :: new_routines
      character(len=4) :: tag
      tag = merge('new:', 'old:', new_routines)
      star%ctrl%helium_diffusion_min = 0.5d0
      call settle_case(tag//' outer edge where Y > 0.5', new_routines)
      star%ctrl%helium_diffusion_min = 0.99d0
      diff_scr%settling_suspended_reported = .false.
      call settle_case(tag//' Y below its floor everywhere', new_routines)
      star%ctrl%helium_diffusion_min = 1.0d-3
      star%ctrl%hydrogen_diffusion_floor = 0.9d0
      diff_scr%settling_suspended_reported = .false.
      call settle_case(tag//' X below its floor everywhere', new_routines)
      star%ctrl%hydrogen_diffusion_floor = 1.0d-3
end subroutine floors

! One settling step from the pristine model; prints the changes.
subroutine settle_case(label, new_routines)
      character(len=*), intent(in) :: label
      logical, intent(in) :: new_routines
      double precision :: step, mass, local_dlnp_dr(json), local_m(json)
      integer :: ierr, k, s
      integer, parameter :: shown(9) = [1, 40, 80, 120, 140, 143, 160, 190, 200]
      integer, parameter :: species(5) = [i_h1, i_he4, i_metals, i_li7, i_be9]

      xa = xa0
      step = dt
      mass = total_mass
      local_dlnp_dr = dlnp_dr
      local_m = m
      ierr = 0
      star%ctrl%use_new_diffusion_routines = new_routines
      call diffuse_elements(step, xa, local_dlnp_dr, log_r, log_rho, local_m, &
           log_t, del, convective, nz, mass, ierr)
      write(*,'(/,a)') '# '//label
      write(*,'(a,i0)') 'ierr = ', ierr
      write(*,'(a)') ' zone   dX  dY  dZ  dLi7  dBe9'
      do k = 1, size(shown)
         write(*,'(i5,5es25.16e3)') shown(k), (xa(species(s),shown(k)) - &
              xa0(species(s),shown(k)), s = 1, size(species))
      end do
      write(*,'(a)') ' species  mass-weighted change'
      do s = 1, 15
         write(*,'(i5,es25.16e3)') s, sum((xa(s,1:nz) - xa0(s,1:nz))*dm(1:nz))/total_mass
      end do
end subroutine settle_case

end program test_element_diffusion
