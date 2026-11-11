!----------------------------------------------------------------------
! convection_zones
!----------------------------------------------------------------------
! Where the model convects (formerly convec and mixcz):
!   find_convection_zones        the mixed regions (the convection zones,
!                                extended by overshoot) and the radiative
!                                regions between them
!   homogenize_convection_zones  mixes each convection zone's composition
!   merge_overlapping_zones      merges mixed regions that overshoot or
!                                semiconvection made meet
!   homogenize_region            mixes one run of shells
! A zone list holds, in row k, the first and the last shell of the k-th
! region from the centre: bounds(k,1) and bounds(k,2). The overshoot
! switches and extents come from convection_settings.
module convection_zones
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use convection_configuration, only: convection_settings, convection_too_many_zones
      use composition_layout, only: i_h1, i_metals, i_o16, n_species_basic, &
           n_species_extended
      use overshoot, only: extend_by_overshoot
      use convection_messages, only: log_too_many_zones, log_zones_merged, &
           merged_by_overshoot
      implicit none
      private
      public :: find_convection_zones, homogenize_convection_zones, &
           merge_overlapping_zones, homogenize_region

! After homogenization O16 is kept at least this fraction of the shell's
! C, N and O, Z minus the envelope's metals outside C, N and O (the F77
! MIXCZ's floor; its source is not recorded).
      real(dp), parameter :: o16_floor_fraction_of_cno = 0.99d-3

contains

!----------------------------------------------------------------------
! find_convection_zones (formerly convec)
!----------------------------------------------------------------------
! The mixed regions of the model: the convection zones from the shells'
! convective flags (runs of more than one shell), extended by overshoot
! (extend_by_overshoot) when an overshoot switch is on, with the regions
! that overshoot makes meet merged; and the radiative regions between
! them. Overshoot is overmixing: the shells it reaches are mixed with the
! zone but keep the radiative gradient (unless adiabatic_overshoot puts
! the adiabatic one below the envelope, compute_temperature_gradients).
!
! In: the composition and log10 rho, P, r, m and T by shell (for the
! overshoot's scale heights); convective_flag(1:num_zones). Its element
! num_zones + 1 is set .false., a sentinel that callers' own scans read.
! Out:
!   mixed_zone_bounds, num_mixed_zones  the mixed regions with overshoot
!       (a zone list, above);
!   mixed_zone_bounds_no_overshoot, num_mixed_zones_no_overshoot  the
!       convection zones themselves (the counts differ after a merge);
!   core_cz_edge  the top shell of the convective core (1 without one);
!   envelope_cz_edge  the bottom shell of the surface convection zone
!       (num_zones without one); both edges are 1 for a fully convective
!       model, which has no overshoot;
!   radiative_zone_bounds, num_radiative_zones  the radiative regions;
!   overshoot_bottom, overshoot_top  the overshoot region below the
!       envelope (extend_by_overshoot), which the caller keeps; left as
!       they were when no overshoot is computed;
!   ierr  convection_too_many_zones (logged) when the model has more
!       convection zones than mixed_zone_bounds has rows; an EOS
!       failure's code from the scale heights.
! Without a convection zone it writes no row of mixed_zone_bounds and
! lists one radiative region, the model; a fully convective model has
! no overshoot and no radiative region (num_radiative_zones = 0).
subroutine find_convection_zones(composition, log_density, log_pressure, log_radius, &
     log_mass, log_temperature, convective_flag, num_zones, &
     radiative_zone_bounds, mixed_zone_bounds, &
     mixed_zone_bounds_no_overshoot, core_cz_edge, envelope_cz_edge, &
     num_radiative_zones, num_mixed_zones, num_mixed_zones_no_overshoot, &
     overshoot_bottom, overshoot_top, ierr)
      real(dp), intent(in) :: composition(:,:), log_density(:), log_pressure(:), &
           log_radius(:), log_mass(:), log_temperature(:)
      logical, intent(inout) :: convective_flag(:)
      integer, intent(in) :: num_zones
      integer, intent(out) :: radiative_zone_bounds(:,:)
      integer, intent(inout) :: mixed_zone_bounds(:,:)
      integer, intent(out) :: mixed_zone_bounds_no_overshoot(:,:)
      integer, intent(out) :: core_cz_edge, envelope_cz_edge
      integer, intent(out) :: num_radiative_zones, num_mixed_zones, &
           num_mixed_zones_no_overshoot
      integer, intent(inout) :: overshoot_bottom, overshoot_top
      integer, intent(out) :: ierr

      ierr = 0
      convective_flag(num_zones+1) = .false.
      call list_convection_zones(convective_flag, num_zones, mixed_zone_bounds, &
           num_mixed_zones, ierr)
      if (ierr /= 0) return
! No convection zone: one radiative region, the whole model.
      if (num_mixed_zones == 0) then
         core_cz_edge = 1
         envelope_cz_edge = num_zones
         num_mixed_zones_no_overshoot = 0
         num_radiative_zones = 1
         radiative_zone_bounds(1,1) = 1
         radiative_zone_bounds(1,2) = num_zones
         return
      end if
      mixed_zone_bounds_no_overshoot(1:num_mixed_zones,:) = &
           mixed_zone_bounds(1:num_mixed_zones,:)
      num_mixed_zones_no_overshoot = num_mixed_zones
! The edges of the convective core and of the surface convection zone,
! before overshoot.
      if (mixed_zone_bounds(1,1) == 1) then
         if (mixed_zone_bounds(1,2) == num_zones) then
! A fully convective model: no overshoot and no radiative region.
            core_cz_edge = 1
            envelope_cz_edge = 1
            num_radiative_zones = 0
            return
         end if
         core_cz_edge = mixed_zone_bounds(1,2)
      else
         core_cz_edge = 1
      end if
      if (mixed_zone_bounds(num_mixed_zones,2) == num_zones) then
         envelope_cz_edge = mixed_zone_bounds(num_mixed_zones,1)
      else
         envelope_cz_edge = num_zones
      end if
! Overshoot moves the edges of the zones it applies to; regions that
! then meet, or share a shell, become one.
      if (convection_settings%core_overshoot .or. convection_settings%envelope_overshoot .or. &
           convection_settings%intermediate_overshoot) then
         call extend_by_overshoot(composition, log_density, log_pressure, log_radius, &
              log_mass, log_temperature, num_zones, mixed_zone_bounds, &
              mixed_zone_bounds_no_overshoot, num_mixed_zones, overshoot_bottom, overshoot_top, ierr)
         if (ierr /= 0) return
         call merge_overlapping_zones(mixed_zone_bounds, num_mixed_zones, &
              merge_shared_shell=.true., cause=merged_by_overshoot)
      end if
      call list_radiative_regions(mixed_zone_bounds, num_mixed_zones, num_zones, &
           radiative_zone_bounds, num_radiative_zones)
end subroutine find_convection_zones

!----------------------------------------------------------------------
! next_convective_run: the first run of convective shells at or above
! shell `from` (at most num_zones): shells first..last, found .false. when
! there is none.
pure subroutine next_convective_run(convective_flag, num_zones, from, first, last, found)
      logical, intent(in) :: convective_flag(:)
      integer, intent(in) :: num_zones, from
      integer, intent(out) :: first, last
      logical, intent(out) :: found
      integer :: k

      first = 0
      last = 0
      found = .false.
      do k = from, num_zones
         if (convective_flag(k)) then
            first = k
            found = .true.
            exit
         end if
      end do
      if (.not. found) return
      last = num_zones
      do k = first + 1, num_zones
         if (.not. convective_flag(k)) then
            last = k - 1
            exit
         end if
      end do
end subroutine next_convective_run

!----------------------------------------------------------------------
! list_convection_zones: the convection zones, the runs of more than one
! convective shell from the centre out, into the zone list bounds. A
! model with more zones than bounds has rows is an error
! (convection_too_many_zones, logged).
subroutine list_convection_zones(convective_flag, num_zones, bounds, num_found, ierr)
      logical, intent(in) :: convective_flag(:)
      integer, intent(in) :: num_zones
      integer, intent(inout) :: bounds(:,:)
      integer, intent(out) :: num_found, ierr
      integer :: from, first, last
      logical :: found

      ierr = 0
      num_found = 0
      from = 1
      do
         call next_convective_run(convective_flag, num_zones, from, first, last, found)
         if (.not. found) return
         from = last + 1
         if (last == first) cycle
         if (num_found == size(bounds, 1)) then
            call log_too_many_zones()
            ierr = convection_too_many_zones
            return
         end if
         num_found = num_found + 1
         bounds(num_found,1) = first
         bounds(num_found,2) = last
      end do
end subroutine list_convection_zones

!----------------------------------------------------------------------
! merge_overlapping_zones: neighbouring regions of the zone list that
! overlap become one. Regions k and k+1 overlap when region k+1 starts
! below the top of region k, or at it when merge_shared_shell is set
! (overshoot sets it; semiconvection does not, as in the F77, so its
! regions may share an edge shell). The merged region runs from region
! k's first shell to region k+1's last; the regions above move down a
! row, and each merge is logged with its cause (merged_by_overshoot or
! merged_by_semiconvection).
subroutine merge_overlapping_zones(bounds, num_regions, merge_shared_shell, cause)
      integer, intent(inout) :: bounds(:,:), num_regions
      logical, intent(in) :: merge_shared_shell
      integer, intent(in) :: cause
      integer :: k, j
      logical :: overlap

      if (num_regions <= 1) return
      k = 1
      do
         if (merge_shared_shell) then
            overlap = bounds(k,2) >= bounds(k+1,1)
         else
            overlap = bounds(k,2) > bounds(k+1,1)
         end if
         if (overlap) then
            call log_zones_merged(cause, bounds(k,:), bounds(k+1,:))
            bounds(k+1,1) = bounds(k,1)
            do j = k, num_regions - 1
               bounds(j,:) = bounds(j+1,:)
            end do
            num_regions = num_regions - 1
            if (k > num_regions - 1) exit
            cycle
         end if
         k = k + 1
         if (k > num_regions - 1) exit
      end do
end subroutine merge_overlapping_zones

!----------------------------------------------------------------------
! list_radiative_regions: the radiative regions of the model, the shells
! outside the mixed regions of the zone list mixed: below the first, between
! each pair, and above the last.
pure subroutine list_radiative_regions(mixed, num_mixed, num_zones, radiative, &
     num_radiative)
      integer, intent(in) :: mixed(:,:), num_mixed, num_zones
      integer, intent(inout) :: radiative(:,:)
      integer, intent(out) :: num_radiative
      integer :: k

      num_radiative = 0
      if (mixed(1,1) > 1) then
         num_radiative = 1
         radiative(1,1) = 1
         radiative(1,2) = mixed(1,1) - 1
      end if
      do k = 1, num_mixed - 1
         num_radiative = num_radiative + 1
         radiative(num_radiative,1) = mixed(k,2) + 1
         radiative(num_radiative,2) = mixed(k+1,1) - 1
      end do
      if (mixed(num_mixed,2) < num_zones) then
         num_radiative = num_radiative + 1
         radiative(num_radiative,1) = mixed(num_mixed,2) + 1
         radiative(num_radiative,2) = num_zones
      end if
end subroutine list_radiative_regions

!----------------------------------------------------------------------
! homogenize_convection_zones (formerly mixcz)
!----------------------------------------------------------------------
! Mixes the composition of each convection zone, the runs of convective
! shells of convective_flag(1:num_zones) (homogenize_region), then floors
! the composition everywhere: X >= 0, Z <= 1 - X, and O16 at least
! o16_floor_fraction_of_cno of Z - Z_noCNO (convection_settings'
! envelope_metals_no_cno). A zone of one shell is not mixed, except at the
! centre (as MIXCZ did). The flags are the caller's: the convection zones
! without overshoot (core's), or with it (rotation's transport flags).
! The basic composition's slots are mixed, or the extended one's with
! use_extended_composition. convective_flag(num_zones + 1) is set
! .false., as find_convection_zones does.
subroutine homogenize_convection_zones(composition, shell_mass, convective_flag, num_zones)
      real(dp), intent(inout) :: composition(:,:)
      real(dp), intent(in) :: shell_mass(:)
      logical, intent(inout) :: convective_flag(:)
      integer, intent(in) :: num_zones
      integer :: num_species, from, first, last
      logical :: found

      convective_flag(num_zones+1) = .false.
      num_species = n_species_basic
      if (convection_settings%use_extended_composition) num_species = n_species_extended
      from = 1
      do
         call next_convective_run(convective_flag, num_zones, from, first, last, found)
         if (.not. found) exit
         from = last + 1
         if (first /= 1 .and. first >= last) cycle
         call homogenize_region(composition, shell_mass, first, last, num_species)
      end do
      call apply_composition_floors(composition, num_zones)
end subroutine homogenize_convection_zones

!----------------------------------------------------------------------
! homogenize_region: shells first..last take, in each of the first
! num_species slots, the region's mass-weighted mean
!   x_i = sum_k x_i(k) dm(k) / sum_k dm(k),
! which keeps each species' mass in the region.
pure subroutine homogenize_region(composition, shell_mass, first, last, num_species)
      real(dp), intent(inout) :: composition(:,:)
      real(dp), intent(in) :: shell_mass(:)
      integer, intent(in) :: first, last, num_species
      real(dp) :: species_sum(num_species), weight_sum
      integer :: k, i

      weight_sum = 0.0d0
      species_sum = 0.0d0
      do k = first, last
         weight_sum = weight_sum + shell_mass(k)
         do i = 1, num_species
            species_sum(i) = species_sum(i) + composition(i,k)*shell_mass(k)
         end do
      end do
      species_sum = species_sum/weight_sum
      do k = first, last
         composition(1:num_species,k) = species_sum
      end do
end subroutine homogenize_region

!----------------------------------------------------------------------
! apply_composition_floors: X >= 0, Z <= 1 - X and the O16 floor in
! shells 1..num_zones (homogenize_convection_zones).
pure subroutine apply_composition_floors(composition, num_zones)
      real(dp), intent(inout) :: composition(:,:)
      integer, intent(in) :: num_zones
      integer :: k

      do k = 1, num_zones
         composition(i_h1,k) = max(composition(i_h1,k), 0.0d0)
         composition(i_metals,k) = min(composition(i_metals,k), 1.0d0 - composition(i_h1,k))
         composition(i_o16,k) = max(composition(i_o16,k), o16_floor_fraction_of_cno* &
              (composition(i_metals,k) - convection_settings%envelope_metals_no_cno))
      end do
end subroutine apply_composition_floors

! Sources: CONVEC and MIXCZ of the Yale code. The early return without
! a convection zone is MHP's (5/91). The F77 wrote the radiative region
! above the last zone into the first region's row (MRZONE(1,2) = M),
! which the 2026 audit fixed; MIXCZ's convective turnover time moved to
! rotation (KC, 2025-05-31). The scans' limits, 11 zones in CONVEC and
! 12 in MIXCZ, became the bounds array's rows (an error beyond them) and
! none in 2026-11 (bug_sweep M248).

end module convection_zones
