!----------------------------------------------------------------------
! mixgrid
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mixgrid.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Builds the equally-spaced-grid quantities needed by mixcom.f90's
! diffusive composition solve: given the first/last unstable zones of
! a region (zone_begin/zone_end), calls getgrid to lay down the
! equally spaced coordinate chi, then computes the equally spaced grid
! masses (equally_spaced_mass) and the geometrically weighted
! diffusion coefficients (equally_spaced_diffusion_coeff) at the zone
! edges, including the Jacobian factor for the transformation from
! radius to the chi coordinate.
subroutine mixgrid(diffusion_coeff, log_density, log_luminosity, &
     log_pressure, log_radius, log_mass, enclosed_mass, shell_mass, &
     log_total_mass, zone_begin, zone_end, convective_flag, num_zones, &
     equally_spaced_diffusion_coeff, equally_spaced_mass, &
     single_interface_flag)
      use numerics_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: diffusion_coeff(json), &
           log_density(json), log_luminosity(json), log_pressure(json), &
           log_radius(json), log_mass(json), enclosed_mass(json), &
           shell_mass(json)
      double precision, intent(in) :: log_total_mass
      integer, intent(in) :: zone_begin, zone_end
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(out) :: equally_spaced_diffusion_coeff(json), &
           equally_spaced_mass(json)
      logical, intent(out) :: single_interface_flag

! common/const/: only solar_luminosity_cgs is used here (to convert
! log_luminosity(M), which is in solar units, to cgs). Remaining
! members are placeholders (uncertain interpretation) preserving
! layout; not referenced in any already-converted file.
      double precision :: solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude
      common/const/ solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude

! common/const1/: ln10 and c4pi are used here. Naming matches
! eqburn.f90/dburn.f90.
      double precision :: ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi
      common/const1/ ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi

! common/const2/: only cgl (used in the log(g) term of the chi
! Jacobian) is used here. Naming matches eqstat2.f90.
      double precision :: gas_constant, radiation_constant_over_3, ca3l, &
           csig, csigl, cgl, cmkh, cmkhn
      common/const2/ gas_constant, radiation_constant_over_3, ca3l, csig, &
           csigl, cgl, cmkh, cmkhn

! common/ctol/: only chi_grid_scale (originally HPTTOL) is used here,
! for the mass/luminosity/pressure scale factors (indices 2/9/11) that
! normalize the terms of the chi-coordinate Jacobian. Remaining
! members are unused layout placeholders; not referenced in any
! already-converted file.
      double precision :: htoler(5,2), fcorr0, fcorri, fcorr, &
           chi_grid_scale(12)
      integer :: niter1, niter2, niter3
      common/ctol/ htoler, fcorr0, fcorri, fcorr, chi_grid_scale, niter1, &
           niter2, niter3

! common/egrid/: the equally spaced coordinate grid (chi/echi/es1) and
! its spacing/point count (dchi/ntot), all used here. Not referenced
! in any already-converted file.
      double precision :: chi(json), echi(json), es1(json), dchi
      integer :: ntot
      common/egrid/ chi, echi, es1, dchi, ntot

! common/mdphy/: only epsm (nuclear energy generation rate per unit
! mass, used in the luminosity term of the chi Jacobian) is used here.
! Naming matches liburn.f90.
      double precision :: amum(json), cpm(json), delm(json), &
           del_adiabatic_mix(json), del_radiative_mix(json), esumm(json), &
           om(json), qdtm(json), thdifm(json), velm(json), viscm(json), &
           epsm(json)
      common/mdphy/ amum, cpm, delm, del_adiabatic_mix, del_radiative_mix, &
           esumm, om, qdtm, thdifm, velm, viscm, epsm

! common/splin/: scratch spline in/out arrays, used as work space here
! by the calls to osplin. Not referenced in any already-converted
! file.
      double precision :: xval(json), yval(json), xtab(json), ytab(json)
      common/splin/ xval, yval, xtab, ytab

      save

      integer :: idx, search_idx, i0, i1, ntab, ntabb
      double precision :: em_top, em_bot
      double precision :: mass_scale, luminosity_scale, pressure_scale
      double precision :: four_pi_rho_r2, dchidr

! FLAG THE SPECIAL CASE OF A SINGLE UNSTABLE INTERFACE AND EXIT
      if (zone_end - zone_begin.le.1) then
         single_interface_flag = .true.
         goto 9999
      else
         single_interface_flag = .false.
      end if
! DEFINE A GRID OF EQUALLY SPACED POINTS.
      call getgrid(log_luminosity, log_pressure, log_mass, zone_begin, &
           zone_end, num_zones)
! GETGRID HAS DEFINED A SET OF CO-ORDINATES (CHI) AND EQUALLY SPACED
! MASS POINTS.  NOW FIND THE OTHER QUANTITIES OF INTEREST AT ZONE
! CENTERS:
! TOTAL ZONE MASSES
! INTERMEDIATE POINTS
      do idx = 2, ntot-1
         equally_spaced_mass(idx) = 0.5d0*(es1(idx+1) - es1(idx-1))
      end do
! SPECIAL TREATMENT OF THE BOUNDARIES; CAN BE CONVECTIVE.
! IF CONVECTIVE SUM OVER ALL SHELLS.  CARE IS NEEDED TO DO BOOK-KEEPING
! PROPERLY AT THE EDGES - TOP IS HALFWAY TO EQUALLY SPACED POINT, NOT
! HALFWAY TO EDGE OF UNEQUALLY SPACED ORIGINAL SET OF POINTS.
!
! CENTER
      em_top = 0.5d0*(es1(2) + es1(1))
      if (zone_begin.gt.1) then
         em_bot = 0.5d0*(enclosed_mass(zone_begin) + &
              enclosed_mass(zone_begin-1))
      else
         em_bot = 0.0d0
      end if
      equally_spaced_mass(1) = em_top - em_bot
      if (zone_begin.gt.1) then
         do search_idx = zone_begin-1, 1, -1
            if (.not.convective_flag(search_idx)) then
               i0 = idx + 1
               goto 10
            end if
            equally_spaced_mass(1) = equally_spaced_mass(1) + &
                 shell_mass(search_idx)
         end do
         i0 = 1
   10    continue
      else
         i0 = 1
      end if
! SURFACE
      em_bot = 0.5d0*(es1(ntot) + es1(ntot-1))
      if (zone_end.lt.num_zones) then
         em_top = 0.5d0*(enclosed_mass(zone_end) + &
              enclosed_mass(zone_end+1))
      else
         em_top = exp(ln10*log_total_mass)
      end if
      equally_spaced_mass(ntot) = em_top - em_bot
      if (zone_end.lt.num_zones) then
         do search_idx = zone_end+1, num_zones
            if (.not.convective_flag(search_idx)) then
               i1 = idx - 1
               goto 20
            end if
            equally_spaced_mass(ntot) = equally_spaced_mass(ntot) + &
                 shell_mass(search_idx)
         end do
         i1 = num_zones
   20    continue
      else
         i1 = num_zones
      end if
! NOW SOLVE FOR QUANTITIES NEEDED AT THE ZONE EDGES.  THESE ARE
! RELATED TO THE DIFFUSION COEFFICIENTS.  UNLIKE THE EQUALLY SPACED
! GRID IN R, WE NEED TO INCLUDE A JACOBIAN TERM FOR THE TRANSFORMATION
! OF VARIABLES.
      ntab = zone_end - zone_begin + 1
      xtab(1) = chi(1)
      do idx = 2, ntab
         xtab(idx) = 0.5d0*(chi(idx) + chi(idx-1))
      end do
      ntabb = ntab + 1
      xtab(ntabb) = chi(ntab)
! DIFFUSION COEFFICIENT FOR MIXING - ASSUME CONSTANT BELOW
! BOTTOM INTERFACE OR ABOVE TOP INTERFACE
      ytab(1) = diffusion_coeff(zone_begin+1)
      do idx = 2, ntab
         search_idx = zone_begin + idx - 1
         ytab(idx) = diffusion_coeff(search_idx)
      end do
      ytab(ntabb) = diffusion_coeff(zone_end)
      xval(1) = chi(1)
      do idx = 2, ntot
         xval(idx) = echi(idx) - 0.5d0*dchi
      end do
      call osplin(xval, equally_spaced_diffusion_coeff, xtab, ytab, ntabb, &
           ntot)
! PRODUCT OF RHO R^2 BY D CHI/DR
      mass_scale = chi_grid_scale(2)
      luminosity_scale = chi_grid_scale(9)*log_luminosity(num_zones)* &
           solar_luminosity_cgs
      pressure_scale = chi_grid_scale(11)
      do idx = 1, ntab
         search_idx = zone_begin + idx - 1
         xtab(idx) = chi(idx)
! D CHI/DR = 1/DM*( D LOG M/DR) + 1/DL*(DL/DR) - 1/DP*(D LOG P/DR)
! OR, USING FAC = 4*PI*RHO*R**2
! D CHI/DR = FAC/(LN 10 * DM * M) + FAC*EPSILON/DL + RHO*GM/(LN10*DP*R**2)
! STORED IN YVAL
         four_pi_rho_r2 = c4pi*exp(ln10*(log_density(search_idx) + &
              2.0d0*log_radius(search_idx)))
         dchidr = four_pi_rho_r2/(ln10*mass_scale*enclosed_mass(search_idx)) &
              + four_pi_rho_r2*epsm(search_idx)/luminosity_scale + &
              exp(ln10*(cgl + log_density(search_idx) + log_mass(search_idx) &
              - log_pressure(search_idx) - &
              2.0d0*log_radius(search_idx)))/(ln10*pressure_scale)
         ytab(idx) = log_density(search_idx) + log10(dchidr) + &
              2.0d0*log_radius(search_idx)
      end do
      call osplin(xval, yval, xtab, ytab, ntab, ntot)
! NOW ADD MULTIPLICATIVE FACTORS TO DIFFUSION COEFFICIENTS
      do idx = 1, ntot
         equally_spaced_diffusion_coeff(idx) = &
              equally_spaced_diffusion_coeff(idx)*exp(ln10*yval(idx))
      end do
 9999 continue
      return
end subroutine mixgrid
