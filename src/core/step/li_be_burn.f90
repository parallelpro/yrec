!----------------------------------------------------------------------
! li_be_burn
!----------------------------------------------------------------------
! The arithmetic of the Li6, Li7 and Be9 burn over a step, as star-blind
! kernels: the model's arrays, the burning rates of every shell at the
! start and end of the step, the base of the surface mixed region and the
! step's accretion come in as arguments, and only the composition
! changes. burn_step's burn_li_be locates the base, gathers the rates and
! the accretion from the star and calls burn_li_be_zones;
! core/test/test_li_be_burn checks the kernels (the closed-form
! survival, the mixed average, the conservation of the accreted share).
!
! Each species decays at its proton-capture rate lambda (per s),
! interpolated in log between the start- and end-of-step rates,
!   lambda(t) = lambda_0 (lambda_1/lambda_0)^(t/dt),
! and survives as X_end = X_start exp(-integral of lambda dt), in closed
! form (net's surviving_fraction).
!
! composition is the star's xa(species, shell); mass_coordinate is m, the
! mass inside each zone centre (g), shell_mass dm (g); rate_start and
! rate_end are (li_be species, shell), in the order of li_be_species.
module li_be_burn
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: i_li6, i_li7, i_be9
      use net_lib, only: surviving_fraction, accreted_survival, &
           linear_in_time, linear_in_log
! exp and log: crmath's under USE_CRMATH
      use math_lib
      implicit none
      private
      public :: mixed_region_bases, li_be_accretion, &
           average_over_mixed_region, burn_li_be_zones

! The three species, by their composition slots, and Be9's place among
! them: its rate, the slowest, decides whether a shell burns at all.
      integer, parameter, public :: num_li_be = 3
      integer, parameter, public :: li_be_species(num_li_be) = [i_li6, i_li7, i_be9]
      integer, parameter :: k_be9 = 3

! A shell whose Be9 rate (per s) is at or below this is too cool to burn.
      real(dp), parameter, public :: min_li_be_rate = 1.0e-32_dp
! An abundance (mass fraction) below this is set to zero.
      real(dp), parameter, public :: min_li_be_abundance = 1.0e-24_dp
! A radiative shell hotter than this log10 T (K) has burned all three.
      real(dp), parameter, public :: li_be_cleared_log_temperature = 7.0_dp

! The base of the surface mixed region (the convection zone plus the
! envelope overshoot) at the start and at the end of the step: the base
! zone, and where the region's true edge sits within it, as a fraction
! of m(k) - m(k-1) below the zone's lower boundary, in [-0.5, 0.5]
! (convection_geometry's mixed_base_mass convention; used only for an
! interior base zone, 1 < k < num_zones).
      type :: mixed_region_bases
         integer :: zone_start = 0, zone_end = 0
         real(dp) :: frac_start = 0.0_dp, frac_end = 0.0_dp
      end type mixed_region_bases

! The step's accretion, for the mixed region's burn: whether gas flows
! in, its mass over the step (g) and its Li6, Li7 and Be9 mass fractions.
      type :: li_be_accretion
         logical :: inflow = .false.
         real(dp) :: mass = 0.0_dp
         real(dp) :: abundance(num_li_be) = 0.0_dp
      end type li_be_accretion

contains

!----------------------------------------------------------------------
! average_over_mixed_region
!----------------------------------------------------------------------
! The Li6, Li7 and Be9 abundances (x_mean) and the natural log of their
! burning rates (log_rate_mean), mass-averaged over the mixed region
! from base_zone = k up to num_zones, and the region's mass
! (region_mass, g). An interior base zone (1 < k < num_zones) is
! adjusted to the region's true edge: its mass changed by
! base_frac (m(k) - m(k-1)), and its rate moved toward the zone below's,
!   rate(k) + (base_frac/2) (rate(k-1) - rate(k)).
pure subroutine average_over_mixed_region(composition, mass_coordinate, &
     shell_mass, rate, base_zone, base_frac, num_zones, x_mean, &
     log_rate_mean, region_mass)
      real(dp), intent(in) :: composition(:,:), mass_coordinate(:), &
           shell_mass(:), rate(:,:)
      integer, intent(in) :: base_zone, num_zones
      real(dp), intent(in) :: base_frac
      real(dp), intent(out) :: x_mean(:), log_rate_mean(:), region_mass
      real(dp) :: rate_sum(num_li_be), zone_rate(num_li_be), mass
      integer :: zone_idx

      x_mean = 0.0d0
      rate_sum = 0.0d0
      region_mass = 0.0d0
      do zone_idx = base_zone,num_zones
         mass = shell_mass(zone_idx)
         zone_rate = rate(:,zone_idx)
         if (zone_idx == base_zone .and. base_zone > 1 .and. &
              base_zone < num_zones) then
            mass = mass + base_frac* &
                 (mass_coordinate(base_zone)-mass_coordinate(base_zone-1))
            zone_rate = zone_rate + 0.5d0*base_frac*(rate(:,base_zone-1) - zone_rate)
         end if
         x_mean = x_mean + composition(li_be_species,zone_idx)*mass
         rate_sum = rate_sum + zone_rate*mass
         region_mass = region_mass + mass
      end do
      x_mean = x_mean/region_mass
      log_rate_mean = log(rate_sum/region_mass)
end subroutine average_over_mixed_region

! ---------------------------------------------------------------
! The mass-averaged Li6, Li7 and Be9 of the mixed region from base_zone
! out, x_mean, and its mass, region_mass [g], as
! average_over_mixed_region takes them, without the rates (a region too
! cool to burn).
pure subroutine mixed_region_mean(composition, mass_coordinate, shell_mass, &
     base_zone, base_frac, num_zones, x_mean, region_mass)
      real(dp), intent(in) :: composition(:,:), mass_coordinate(:), shell_mass(:)
      integer, intent(in) :: base_zone, num_zones
      real(dp), intent(in) :: base_frac
      real(dp), intent(out) :: x_mean(:), region_mass
      real(dp) :: mass
      integer :: zone_idx

      x_mean = 0.0d0
      region_mass = 0.0d0
      do zone_idx = base_zone,num_zones
         mass = shell_mass(zone_idx)
         if (zone_idx == base_zone .and. base_zone > 1 .and. &
              base_zone < num_zones) mass = mass + base_frac* &
              (mass_coordinate(base_zone)-mass_coordinate(base_zone-1))
         x_mean = x_mean + composition(li_be_species,zone_idx)*mass
         region_mass = region_mass + mass
      end do
      x_mean = x_mean/region_mass
end subroutine mixed_region_mean

!----------------------------------------------------------------------
! burn_li_be_zones
!----------------------------------------------------------------------
! Li6, Li7 and Be9 over the step timestep (s), in every shell, in order:
!   1. radiative shells below both bases, each on its own rates; the
!      loop stops at the first shell too cool to burn (Be9's rate at or
!      below min_li_be_rate at the start or end), skips shells already
!      depleted, and clears shells hotter than
!      li_be_cleared_log_temperature;
!   2. the mixed region as one unit, on its mass-averaged rates
!      (average_over_mixed_region), starting from the abundance mixed
!      over the deeper region when it deepened; with accretion (MHP
!      05/02), the accreted gas burns from its arrival on
!      (accreted_survival, S_acc) and is mixed in by mass,
!        x = (x M_cz + x_acc S_acc M_acc)/(M_cz + M_acc),
!      which assumes that all of it lands in a region of the end-of-step
!      depth; the new abundance is set from the shallower base out (when
!      the region deepened, the shells between the two bases are left as
!      they are). Nothing more burns when the region's base is too cool
!      at the start or at the end, but the accreted gas is still mixed in
!      by mass, unburnt, over the end-of-step region (bug_sweep M383);
!   3. shells the region left during the step (a retreating convection
!      zone), each on one rate held over the step, weighted in the log
!      between its own end-of-step rate and the region's start-of-step
!      rate by the time spent in each: the base moved out in mass from
!      m_beg (the midpoint of the old base zone and the one below) to
!      m_end (likewise for the new), and a shell's radiative fraction of
!      the step is (m_end - m)/(m_end - m_beg), 1 at the old base and 0
!      at the new one.
! Abundances below min_li_be_abundance are set to zero.
pure subroutine burn_li_be_zones(timestep, composition, mass_coordinate, &
     shell_mass, log_temperature, rate_start, rate_end, base, accretion, &
     num_zones)
      real(dp), intent(in) :: timestep
      real(dp), intent(inout) :: composition(:,:)
      real(dp), intent(in) :: mass_coordinate(:), shell_mass(:), &
           log_temperature(:), rate_start(:,:), rate_end(:,:)
      type(mixed_region_bases), intent(in) :: base
      type(li_be_accretion), intent(in) :: accretion
      integer, intent(in) :: num_zones

      integer :: zone_idx, isp, min_zone, max_zone
      real(dp) :: rate
      real(dp) :: x_cz_start(num_li_be), x_cz_end(num_li_be), &
           log_rate_cz_start(num_li_be), log_rate_cz_end(num_li_be), &
           cz_mass_start, cz_mass_end
      real(dp) :: mass_coord_beg, mass_coord_end, delta_mass, radiative_frac

! min_zone is the deeper of the two bases, max_zone the shallower.
      min_zone = min(base%zone_end,base%zone_start)
      max_zone = max(base%zone_end,base%zone_start)

! 1. The radiative interior, up to the deeper of the two bases.
      do zone_idx = 1,min_zone-1
         if(rate_end(k_be9,zone_idx) <= min_li_be_rate .or. &
              rate_start(k_be9,zone_idx) <= min_li_be_rate)exit
         if(all(composition(li_be_species,zone_idx) < min_li_be_abundance))cycle
         if(log_temperature(zone_idx) > li_be_cleared_log_temperature)then
            composition(li_be_species,zone_idx) = 0.0d0
            cycle
         end if
         do isp = 1,num_li_be
            composition(li_be_species(isp),zone_idx) = composition(li_be_species(isp),zone_idx)* &
                 surviving_fraction(rate_start(isp,zone_idx), rate_end(isp,zone_idx), &
                 timestep, linear_in_log)
            if(composition(li_be_species(isp),zone_idx) < min_li_be_abundance) &
                 composition(li_be_species(isp),zone_idx) = 0.0d0
         end do
      end do

! 2. The mixed region as one unit.
      if (rate_start(k_be9,base%zone_start) <= min_li_be_rate .or. &
           rate_end(k_be9,base%zone_end) <= min_li_be_rate) then
         if (accretion%inflow) then
            call mixed_region_mean(composition, mass_coordinate, shell_mass, &
                 base%zone_end, base%frac_end, num_zones, x_cz_end, cz_mass_end)
            x_cz_end = (x_cz_end*cz_mass_end + accretion%abundance*accretion%mass)/ &
                 (accretion%mass + cz_mass_end)
            where (x_cz_end < min_li_be_abundance) x_cz_end = 0.0d0
            do zone_idx = max_zone,num_zones
               composition(li_be_species,zone_idx) = x_cz_end
            end do
         end if
         return
      end if
      call average_over_mixed_region(composition, mass_coordinate, shell_mass, &
           rate_start, base%zone_start, base%frac_start, num_zones, &
           x_cz_start, log_rate_cz_start, cz_mass_start)
      call average_over_mixed_region(composition, mass_coordinate, shell_mass, &
           rate_end, base%zone_end, base%frac_end, num_zones, &
           x_cz_end, log_rate_cz_end, cz_mass_end)
! When the region deepened, start from the abundance mixed over the
! deeper region.
      if(base%zone_end < base%zone_start) x_cz_start = x_cz_end
      do isp = 1,num_li_be
         x_cz_end(isp) = x_cz_start(isp)*surviving_fraction(exp(log_rate_cz_start(isp)), &
              exp(log_rate_cz_end(isp)), timestep, linear_in_log)
         if (accretion%inflow) x_cz_end(isp) = (x_cz_end(isp)*cz_mass_end + &
              accretion%abundance(isp)* &
              accreted_survival(exp(log_rate_cz_start(isp)), exp(log_rate_cz_end(isp)), &
              timestep, linear_in_log)*accretion%mass)/ &
              (cz_mass_end + accretion%mass)
         if(x_cz_end(isp) < min_li_be_abundance) x_cz_end(isp) = 0.0d0
      end do
      do zone_idx = max_zone,num_zones
         composition(li_be_species,zone_idx) = x_cz_end
      end do

! 3. Shells that began in the mixed region and ended radiative.
      if (base%zone_end <= base%zone_start) return
      if(base%zone_start > 1)then
         mass_coord_beg = 0.5d0*(mass_coordinate(base%zone_start)+mass_coordinate(base%zone_start-1))
      else
         mass_coord_beg = 0.0d0
      end if
      mass_coord_end = 0.5d0*(mass_coordinate(base%zone_end)+mass_coordinate(base%zone_end-1))
      delta_mass = mass_coord_end - mass_coord_beg
      do zone_idx = base%zone_start,base%zone_end-1
! Stop at the first shell too cool to burn.
         if(rate_end(k_be9,zone_idx) <= min_li_be_rate)exit
         radiative_frac = (mass_coord_end - mass_coordinate(zone_idx))/delta_mass
         do isp = 1,num_li_be
!   ln rate = radiative_frac ln rate_end + (1 - radiative_frac) ln rate_cz.
            rate = exp(radiative_frac*log(rate_end(isp,zone_idx)) + &
                 (1.0d0-radiative_frac)*log_rate_cz_start(isp))
            composition(li_be_species(isp),zone_idx) = composition(li_be_species(isp),zone_idx)* &
                 surviving_fraction(rate, rate, timestep, linear_in_time)
            if(composition(li_be_species(isp),zone_idx) < min_li_be_abundance) &
                 composition(li_be_species(isp),zone_idx) = 0.0d0
         end do
      end do
end subroutine burn_li_be_zones

end module li_be_burn
