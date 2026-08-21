!----------------------------------------------------------------------
! getgrid
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original getgrid.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! This routine forms the sum of chi, defined as
! chi = log(M)/dm + L/(Ltot*dl) - log(P)/dp, where dm, dp, and dl are
! the point spacings for log mass, luminosity, and pressure
! respectively. This routine transforms these variables to an equally
! spaced grid from the first point to the last point. It returns a
! set of equally spaced chi values and their location in mass.
subroutine getgrid(log_luminosity, log_pressure, log_mass, zone_begin, &
     zone_end, num_zones)
      use osplin_mod
      implicit none
      integer, parameter :: json = 5000

! log_luminosity (originally HL): this file's own comment describes
! this term as "L/(LTOT*DL)" (not "LOG(L)/DL", unlike the LOG(M) and
! LOG(P) terms), but both already-converted callers (mixgrid.f90,
! rotgrid.f90) pass their log_luminosity array into this argument
! position, so it is named accordingly here.
      double precision, intent(in) :: log_luminosity(json), log_pressure(json), &
           log_mass(json)
      integer, intent(in) :: zone_begin, zone_end, num_zones

! common/const1/: only ln10 is used here. Naming matches eqburn.f90.
      double precision :: ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi
      common/const1/ ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi

! common/ctol/: only chi_grid_scale (originally HPTTOL) is used here.
! Naming matches mixgrid.f90.
      double precision :: htoler(5,2), fcorr0, fcorri, fcorr, &
           chi_grid_scale(12)
      integer :: niter1, niter2, niter3
      common/ctol/ htoler, fcorr0, fcorri, fcorr, chi_grid_scale, niter1, &
           niter2, niter3

! common/egrid/: the equally spaced coordinate grid (chi/echi/es1) and
! its spacing/point count (dchi/ntot), all set here. Naming matches
! mixgrid.f90.
      double precision :: chi(json), echi(json), es1(json), dchi
      integer :: ntot
      common/egrid/ chi, echi, es1, dchi, ntot

      double precision :: log_mass_table(json)
      save

      integer :: num_points_in_range, zone_index, range_index
      double precision :: mass_scale, luminosity_scale, pressure_scale

! USE THE MODEL CRITERIA FOR ASSIGNING THE SPACING BETWEEN GRID POINTS.
      mass_scale = chi_grid_scale(2)
      luminosity_scale = chi_grid_scale(9)*log_luminosity(num_zones)
      pressure_scale = chi_grid_scale(11)
      num_points_in_range = zone_end - zone_begin + 1
      do range_index = 1, num_points_in_range
         zone_index = zone_begin + range_index - 1
         chi(range_index) = log_mass(zone_index)/mass_scale + &
              log_luminosity(zone_index)/luminosity_scale - &
              log_pressure(zone_index)/pressure_scale
      end do
! DEFINE NTOT EQUAL TO NTAB
      ntot = num_points_in_range
! TOTAL NUMBER OF ZONES IS MODULUS OF FINAL CHI PLUS ONE.
!      NTOT = INT(CHI(NTAB)-CHI(1))+1
! FOR ROTATION PURPOSES, DEFINE THE MINIMUM NUMBER OF
! EQUALLY SPACED SHELLS AS 3.
      ntot = max(ntot,3)
! EQUALLY SPACED INCREMENT IN CHI
!      CHT = 0.5D0*(CHI(NTAB)+CHI(NTAB-1))
!      CHB = 0.5D0*(CHI(2)+CHI(1))
!      DCHI = (CHT-CHB)/FLOAT(NTOT-2)
      dchi = (chi(num_points_in_range)-chi(1))/dfloat(ntot-1)
! ASSIGN VECTOR OF EQUALLY SPACED CHI
!       ECHI(1) = CHB - 0.5D0*DCHI
!       DO I = 2,NTOT
!          ECHI(I) = ECHI(I-1)+DCHI
!       END DO
      echi(1) = chi(1)
      do zone_index = 2, ntot
         echi(zone_index) = echi(zone_index-1)+dchi
      end do
! NOW ASSIGN MASSES TO THE NEW EQUALLY SPACED GRID POINTS.
! PERFORM INTERPOLATION IN LOG MASS - CONSISTENT WITH CHI DEFINED
! IN LOG M, LOG P, ETC.
      do zone_index = 1, num_points_in_range
         log_mass_table(zone_index) = log_mass(zone_begin+zone_index-1)
      end do
      call osplin(echi,es1,chi,log_mass_table,num_points_in_range,ntot)
! TRANSFORM TO PHYSICAL MASS (GM)
      do zone_index = 1, ntot
         es1(zone_index) = exp(ln10*es1(zone_index))
      end do
      return
end subroutine getgrid
