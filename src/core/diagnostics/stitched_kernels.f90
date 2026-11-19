!----------------------------------------------------------------------
! stitched_kernels
!----------------------------------------------------------------------
! The arithmetic of the stitched model on plain arrays: a star-blind
! kernel (docs/GUIDELINES.md, star rule 6). build_stitched_model
! (stitched_model.f90) and compute_seismic_observables (observables.f90)
! gather from the star, call these and store the results; nothing here
! reads the star, prints or keeps state between calls.
!   accumulate_atmosphere_heights  the atmosphere points' heights above
!                                  the photosphere;
!   map_extended_grid              the index map of the extended grid:
!                                  interior, envelope, atmosphere;
!   thermal_brunt_n2               N^2 of one point, the thermal form;
!   difference_brunt_n2            N^2 by centred differences where the
!                                  thermal form is not negative;
!   acoustic_radius                int dr/c, for the large separation;
!   buoyancy_radius                int N/r dr over the innermost g-mode
!                                  cavity, for the period spacing.
! test_stitched checks each on synthetic arrays.
module stitched_kernels
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use math_lib
      implicit none
      private
      public :: accumulate_atmosphere_heights, map_extended_grid, thermal_brunt_n2, &
           difference_brunt_n2, acoustic_radius, buoyancy_radius
      public :: region_interior, region_envelope, region_atmosphere

! The regions of the extended grid (the star's stitched%region).
      integer, parameter :: region_interior = 1, region_envelope = 2, &
           region_atmosphere = 3

contains

! ---------------------------------------------------------------
! The heights [cm] of the atmosphere's points above its deepest one, the
! photosphere (tau = 2/3). The points run from the outermost (1) inward,
! and step_length(k) is the geometric length between points k-1 and k
! (the envelope kernel's Cox p. 590 quadrature), so point i sits
! sum(step_length(i+1:)) above the last. step_length(1) enters no
! height. height has the size of step_length; nothing is set when there
! is no point.
pure subroutine accumulate_atmosphere_heights(step_length, height)
      real(dp), intent(in) :: step_length(:)
      real(dp), intent(inout) :: height(:)
      integer :: i, num_points

      num_points = size(step_length)
      if (num_points > 0) then
         height(num_points) = 0.0d0
         do i = num_points - 1, 1, -1
            height(i) = height(i+1) + step_length(i+1)
         end do
      end if
end subroutine accumulate_atmosphere_heights

! ---------------------------------------------------------------
! The index map of the extended grid, from the centre outward: point j
! is point index(j) of region(j). First the num_interior interior points.
! Then the envelope's, which run from the fitting point to the
! photosphere: every one whose log r is above the interior's outermost,
! interior_top_log_radius, but the first, which repeats the fitting point
! (to the integration's round-off, so its radius may land just above the
! interior's, with a density that rises outward: bug_sweep M359) and is
! skipped by its index. Then the
! atmosphere's, which run from the outside in and are taken in reverse:
! every one with a positive height above the photosphere, which skips the
! deepest (the envelope's last point again). num_interior_envelope is the
! last point that is not in the atmosphere, and num_points the last point.
! The map stops, without an error, at the size of region.
pure subroutine map_extended_grid(num_interior, interior_top_log_radius, &
     envelope_log_radius, atmosphere_height, region, index, num_interior_envelope, &
     num_points)
      integer, intent(in) :: num_interior
      real(dp), intent(in) :: interior_top_log_radius, envelope_log_radius(:), &
           atmosphere_height(:)
      integer, intent(inout) :: region(:), index(:)
      integer, intent(out) :: num_interior_envelope, num_points
      integer :: i

      num_points = 0
      do i = 1, num_interior
         num_points = num_points + 1
         region(num_points) = region_interior
         index(num_points) = i
      end do
      do i = 1, size(envelope_log_radius)
         if (num_points >= size(region)) exit
         if (i == 1 .or. envelope_log_radius(i) <= interior_top_log_radius) cycle
         num_points = num_points + 1
         region(num_points) = region_envelope
         index(num_points) = i
      end do
      num_interior_envelope = num_points
      do i = size(atmosphere_height), 1, -1
         if (num_points >= size(region)) exit
         if (atmosphere_height(i) <= 0.0d0) cycle
         num_points = num_points + 1
         region(num_points) = region_atmosphere
         index(num_points) = i
      end do
end subroutine map_extended_grid

! ---------------------------------------------------------------
! The Brunt-Vaisala frequency squared of one point in its thermal form
! [s^-2],
!   N^2 = g^2 (rho/P) delta (nabla_ad - nabla),
! exact for a homogeneous composition (it has no mu-gradient term), from
! the gravity [cm/s^2], the density [g/cm^3], the pressure [dyn/cm^2],
! delta = -(dln rho/dln T)_P and the adiabatic and actual gradients.
pure real(dp) function thermal_brunt_n2(gravity, density, pressure, delta, &
     nabla_ad, nabla)
      real(dp), intent(in) :: gravity, density, pressure, delta, nabla_ad, nabla

      thermal_brunt_n2 = gravity*gravity*(density/pressure)*delta*(nabla_ad - nabla)
end function thermal_brunt_n2

! ---------------------------------------------------------------
! The second pass over N^2 (brunt_n2, on entry the thermal form at every
! point): where it is positive, the gradient form by centred
! differences in r,
!   N^2 = g [ (1/Gamma_1) dln P/dr - dln rho/dr ],  g = G m/r^2,
! which carries the mu-gradient term the thermal form leaves out: in the
! composition gradient above a retreating core or around the hydrogen
! shell that term dominates (on a 1.2 Msun subgiant the thermal form
! understated the g cavity's buoyancy by about 30%, 12% in the period
! spacing, and moved the l = 1 mixed modes by up to 16 uHz). A point
! whose thermal N^2 is negative (Schwarzschild-unstable) keeps it: the
! mixture there is homogeneous, so the thermal form is exact, while the
! centred difference of a near-adiabatic stratification is cancellation
! noise. A point where nabla equals nabla_ad exactly, as in an efficient
! convection zone, keeps its thermal 0 for the same reason (the test is
! "> 0"; bug_sweep M361). The positive N^2 the difference gives at the base of
! a convection zone (the buoyancy of the mu discontinuity) lies on the
! radiative side and is kept. The innermost point takes its neighbour's
! N^2 scaled as r^2, N^2(2) (r_1/r_2)^2 (N^2 is proportional to r^2 about
! the centre; bug_sweep M363), and the outermost copies its neighbour's.
! The arrays hold the points from
! the centre outward: r [cm], m [g], P, rho, Gamma_1; newton_g is G.
pure subroutine difference_brunt_n2(radius, mass, pressure, density, gamma1, &
     newton_g, brunt_n2)
      real(dp), intent(in) :: radius(:), mass(:), pressure(:), density(:), gamma1(:), &
           newton_g
      real(dp), intent(inout) :: brunt_n2(:)
      real(dp) :: gravity
      integer :: j, num_points

      num_points = size(radius)
      do j = 2, num_points - 1
         if (brunt_n2(j) > 0.0d0 .and. radius(j) > 0.0d0 .and. &
              radius(j+1) > radius(j-1) .and. gamma1(j) > 0.0d0) then
            gravity = newton_g*mass(j)/(radius(j)*radius(j))
            brunt_n2(j) = gravity*( &
                 (log(pressure(j+1)) - log(pressure(j-1)))/gamma1(j) - &
                 (log(density(j+1)) - log(density(j-1))) ) / &
                 (radius(j+1) - radius(j-1))
         end if
      end do
      if (num_points >= 2) then
         if (radius(2) > 0.0d0) then
            brunt_n2(1) = brunt_n2(2)*(radius(1)/radius(2))**2
         else
            brunt_n2(1) = brunt_n2(2)
         end if
         brunt_n2(num_points) = brunt_n2(num_points-1)
      end if
end subroutine difference_brunt_n2

! ---------------------------------------------------------------
! The sound travel time [s] from the centre to the last point,
!   int_0 dr/c,
! by the trapezoidal rule over the points (log10 r in cm, c in cm/s), and
! inside the first point r_1/c_1 (c is constant to first order there;
! bug_sweep M362); an interval with a sound speed that is not positive at
! either end is left out, and so is the centre's part when c_1 is not
! positive.
pure real(dp) function acoustic_radius(log_radius, sound_speed)
      real(dp), intent(in) :: log_radius(:), sound_speed(:)
      real(dp) :: dr
      integer :: k

      acoustic_radius = 0.0d0
      if (size(log_radius) < 1) return
      if (sound_speed(1) > 0.0d0) acoustic_radius = exp10(log_radius(1))/sound_speed(1)
      do k = 2, size(log_radius)
         if (sound_speed(k) <= 0.0d0 .or. sound_speed(k-1) <= 0.0d0) cycle
         dr = exp10(log_radius(k)) - exp10(log_radius(k-1))
         acoustic_radius = acoustic_radius + dr*0.5d0 &
              *(1.0d0/sound_speed(k) + 1.0d0/sound_speed(k-1))
      end do
end function acoustic_radius

! ---------------------------------------------------------------
! The buoyancy radius [rad/s]
!   int N/r dr
! over the innermost g-mode cavity: walking outward from the second
! point, each point with N^2 > 0 adds N/r times half the span to its
! neighbours, and the walk stops at the first point with N^2 <= 0 after
! the cavity (after WB's MESA treatment). Zero without a cavity. The
! first point is not walked (N/r is 0/0 at the centre); when the second
! is in the cavity, the part of the integral inside (r_1 + r_2)/2 is its
! N/r times (r_1 + r_2)/2, N being proportional to r near the centre
! (bug_sweep M362). The last point has no outer neighbour. log_radius is
! log10 r in cm, brunt_n2 in s^-2.
pure real(dp) function buoyancy_radius(log_radius, brunt_n2)
      real(dp), intent(in) :: log_radius(:), brunt_n2(:)
      real(dp) :: r_inner, r_here, r_outer
      logical :: entered_g_mode_cavity
      integer :: k

      buoyancy_radius = 0.0d0
      entered_g_mode_cavity = .false.
      if (size(log_radius) >= 3) then
         if (brunt_n2(2) > 0.0d0) buoyancy_radius = sqrt(brunt_n2(2))/exp10(log_radius(2)) &
              *0.5d0*(exp10(log_radius(1)) + exp10(log_radius(2)))
      end if
      do k = 2, size(log_radius) - 1
         if (brunt_n2(k) > 0.0d0) then
            entered_g_mode_cavity = .true.
            r_inner = exp10(log_radius(k-1))
            r_here = exp10(log_radius(k))
            r_outer = exp10(log_radius(k+1))
            buoyancy_radius = buoyancy_radius &
                 + sqrt(brunt_n2(k))/r_here*0.5d0*(r_outer - r_inner)
         else if (entered_g_mode_cavity) then
            exit
         end if
      end do
end function buoyancy_radius

end module stitched_kernels
