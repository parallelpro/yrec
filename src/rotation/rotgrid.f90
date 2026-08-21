!----------------------------------------------------------------------
! rotgrid
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original rotgrid.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Builds the equally-spaced-grid quantities needed by dadcoeft.f90's
! angular-momentum-transport diffusion solve: given the first/last
! unstable zones of a region (zone_begin/zone_end), calls getgrid to
! lay down the equally spaced coordinate chi, then computes the
! equally spaced grid moment of inertia (eq_moment_of_inertia),
! specific angular momentum (eq_angular_momentum), mass (eq_mass),
! angular velocity (eq_omega), and the geometrically weighted
! angular-momentum/mixing diffusion coefficients
! (eq_am_diffusion_coeff/eq_mixing_diffusion_coeff) at the zone edges,
! including the Jacobian factor for the transformation from radius to
! the chi coordinate. Companion routine to mixgrid.f90 (composition
! diffusion) for the rotation/angular-momentum-transport case.
subroutine rotgrid(am_diffusion_coeff, mixing_diffusion_coeff, log_density, &
     moment_of_inertia, specific_angular_momentum, log_luminosity, &
     log_pressure, log_radius, radius, log_mass, enclosed_mass, shell_mass, &
     log_total_mass, zone_begin, zone_end, am_transport_convective_flag, &
     num_zones, omega, grid_spacing, eq_am_diffusion_coeff, &
     eq_mixing_diffusion_coeff, eq_moment_of_inertia, eq_angular_momentum, &
     eq_mass, eq_omega, single_interface_flag)
      use mdphy_lib
      use const_lib
      use numerics_lib
      implicit none
      integer, parameter :: json = 5000

! INPUT VARIABLES
      double precision, intent(in) :: am_diffusion_coeff(json), &
           mixing_diffusion_coeff(json), log_density(json), &
           moment_of_inertia(json), specific_angular_momentum(json), &
           log_luminosity(json), log_pressure(json), log_radius(json), &
           radius(json), log_mass(json), enclosed_mass(json), &
           shell_mass(json), log_total_mass
      integer, intent(in) :: zone_begin, zone_end
      logical, intent(in) :: am_transport_convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(in) :: omega(json)
! OUTPUT VARIABLES
      double precision, intent(out) :: grid_spacing, &
           eq_am_diffusion_coeff(json), eq_mixing_diffusion_coeff(json), &
           eq_moment_of_inertia(json), eq_angular_momentum(json), &
           eq_mass(json), eq_omega(json)
      logical, intent(out) :: single_interface_flag





! common/egrid/: the equally spaced coordinate grid (chi/echi/es1) and
! its spacing/point count (dchi/ntot), all used here. Naming matches
! mixgrid.f90.
      double precision :: chi(json), echi(json), es1(json), dchi
      integer :: ntot
      common/egrid/ chi, echi, es1, dchi, ntot


! common/splin/: scratch spline in/out arrays, used as work space here
! by the calls to osplin. Naming matches mixgrid.f90.
      double precision :: xval(json), yval(json), xtab(json), ytab(json)
      common/splin/ xval, yval, xtab, ytab

! MHP 05/02
! common/difad/: am_advective_coeff/am_diffusive_coeff (originally
! ECOD3/ECOD4), the diffusion coefficients for the terms involving
! d(omega)/dr (advective) and d/dr(d omega/dr) (diffusive) in the
! angular-momentum transport equation. Both are written here. Naming
! is local to this batch (rotation-transport files); not referenced in
! any already-converted file outside this batch.
      double precision :: am_advective_coeff(json), am_diffusive_coeff(json)
      common/difad/ am_advective_coeff, am_diffusive_coeff

! common/difad2/: es_advective_velocity/es_diffusive_velocity
! (originally VESA/VESD, MHP 11/94), the advective and diffusive
! components of the Eddington-Sweet meridional-circulation velocity
! when the combined diffusion+advection treatment is used. Only these
! two (not the "_prev" pair) are used here; naming is local to this
! batch.
      double precision :: es_advective_velocity(json), &
           es_advective_velocity_prev(json), es_diffusive_velocity(json), &
           es_diffusive_velocity_prev(json)
      common/difad2/ es_advective_velocity, es_advective_velocity_prev, &
           es_diffusive_velocity, es_diffusive_velocity_prev



      double precision :: eq_reduced_moment_of_inertia(json)
      save

      integer :: ntab, i, ii, i0, i1, ntabb
      double precision :: emtop, embot, mass_scale_factor, &
           luminosity_scale_factor, pressure_scale_factor, scale_factor, &
           dchi_dr

! FLAG THE SPECIAL CASE OF A SINGLE UNSTABLE INTERFACE AND EXIT
      if (zone_end-zone_begin.le.1) then
         single_interface_flag = .true.
         goto 9999
      else
         single_interface_flag = .false.
      end if
! DEFINE A GRID OF EQUALLY SPACED POINTS.
      call getgrid(log_luminosity,log_pressure,log_mass,zone_begin, &
           zone_end,num_zones)
! GETGRID HAS DEFINED A SET OF CO-ORDINATES (CHI) AND EQUALLY SPACED
! MASS POINTS.  NOW FIND THE OTHER QUANTITIES OF INTEREST AT ZONE
! CENTERS:
! J/M (TO GET JTOT)
! I/MR^2 (TO GET ITOT)
! OMEGA
! R (NEEDED FOR I CALCULATION) - STORED IN YVAL
      ntab = zone_end - zone_begin + 1
      do i = 1,ntab
         ii = zone_begin + i - 1
         xtab(i) = chi(i)
         ytab(i) = log_radius(ii)
      end do
      do i = 1,ntot
         xval(i) = echi(i)
      end do
      call osplin(xval,yval,xtab,ytab,ntab,ntot)
! I/MR^2
      do i = 1,ntab
         ii = zone_begin + i - 1
         ytab(i) = moment_of_inertia(ii)/(shell_mass(ii)*radius(ii)**2)
      end do
! MHP 05/02 STORE I/MR^2 FOR LATER USE IN
! DIFFUSION COEFFICIENTS FOR ANGULAR MOMENTUM
! TRANSPORT
      call osplin(xval,eq_reduced_moment_of_inertia,xtab,ytab,ntab,ntot)
! J/M
      do i = 1,ntab
         ii = zone_begin + i - 1
         ytab(i) = specific_angular_momentum(ii)
      end do
      call osplin(xval,eq_angular_momentum,xtab,ytab,ntab,ntot)
! OMEGA
      do i = 1,ntab
         ii = zone_begin + i - 1
         ytab(i) = omega(ii)
      end do
      call osplin(xval,eq_omega,xtab,ytab,ntab,ntot)
! CONVERT TO TOTAL I AND TOTAL J OF THE ZONES.
! INTERMEDIATE POINTS
      do i = 2, ntot-1
         eq_mass(i) = 0.5d0*(es1(i+1)-es1(i-1))
! MHP 05/02 CHANGED TO REFLECT THE FACT THAT
! THE INFO PREVIOUSLY STORED IN EI IS NOW IN EI0
         eq_moment_of_inertia(i) = eq_reduced_moment_of_inertia(i)* &
              eq_mass(i)*exp(ln10*2.0d0*yval(i))
         eq_angular_momentum(i) = eq_angular_momentum(i)*eq_mass(i)
      end do
! SPECIAL TREATMENT OF THE BOUNDARIES; CAN BE CONVECTIVE.
! IF CONVECTIVE SUM OVER ALL SHELLS.  CARE IS NEEDED TO DO BOOK-KEEPING
! PROPERLY AT THE EDGES - TOP IS HALFWAY TO EQUALLY SPACED POINT, NOT
! HALFWAY TO EDGE OF UNEQUALLY SPACED ORIGINAL SET OF POINTS.
!
! CENTER
      emtop = 0.5d0*(es1(2)+es1(1))
      if (zone_begin.gt.1) then
         embot = 0.5d0*(enclosed_mass(zone_begin)+enclosed_mass(zone_begin-1))
      else
         embot = 0.0d0
      end if
      eq_mass(1) = emtop - embot
! MHP 05/02 CHANGED TO REFLECT THE FACT THAT
! THE INFO PREVIOUSLY STORED IN EI IS NOW IN EI0
      eq_moment_of_inertia(1) = eq_reduced_moment_of_inertia(1)*eq_mass(1)* &
           radius(zone_begin)**2
      eq_angular_momentum(1) = eq_angular_momentum(1)*eq_mass(1)
      if (zone_begin.gt.1) then
         do ii = zone_begin-1,1,-1
            if (.not.am_transport_convective_flag(ii)) then
               i0 = i + 1
               goto 10
            end if
            eq_mass(1) = eq_mass(1)+shell_mass(ii)
            eq_moment_of_inertia(1) = eq_moment_of_inertia(1)+ &
                 moment_of_inertia(ii)
            eq_angular_momentum(1) = eq_angular_momentum(1)+ &
                 specific_angular_momentum(ii)*shell_mass(ii)
         end do
         i0 = 1
 10      continue
      else
         i0 = 1
      end if
! SURFACE
      embot = 0.5d0*(es1(ntot)+es1(ntot-1))
      if (zone_end.lt.num_zones) then
         emtop = 0.5d0*(enclosed_mass(zone_end)+enclosed_mass(zone_end+1))
      else
         emtop = exp(ln10*log_total_mass)
      end if
      eq_mass(ntot) = emtop - embot
! MHP 05/02 CHANGED TO REFLECT THE FACT THAT
! THE INFO PREVIOUSLY STORED IN EI IS NOW IN EI0
      eq_moment_of_inertia(ntot) = eq_reduced_moment_of_inertia(ntot)* &
           eq_mass(ntot)*radius(zone_end)**2
      eq_angular_momentum(ntot) = eq_angular_momentum(ntot)*eq_mass(ntot)
      if (zone_end.lt.num_zones) then
         do ii = zone_end+1,num_zones
            if (.not.am_transport_convective_flag(ii)) then
               i1 = i -1
               goto 20
            end if
            eq_mass(ntot) = eq_mass(ntot)+shell_mass(ii)
            eq_moment_of_inertia(ntot) = eq_moment_of_inertia(ntot)+ &
                 moment_of_inertia(ii)
            eq_angular_momentum(ntot) = eq_angular_momentum(ntot)+ &
                 specific_angular_momentum(ii)*shell_mass(ii)
         end do
         i1 = num_zones
 20      continue
      else
         i1 = num_zones
      end if
! NOW SOLVE FOR QUANTITIES NEEDED AT THE ZONE EDGES.  THESE ARE
! RELATED TO THE DIFFUSION COEFFICIENTS.  UNLIKE THE EQUALLY SPACED
! GRID IN R, WE NEED TO INCLUDE A JACOBIAN TERM FOR THE TRANSFORMATION
! OF VARIABLES.
! DIFFUSION COEFFICIENT FOR ANGULAR MOMENTUM - ASSUME CONSTANT BELOW
! BOTTOM INTERFACE OR ABOVE TOP INTERFACE
      xtab(1) = chi(1)
      ytab(1) = am_diffusion_coeff(zone_begin+1)
      do i = 2,ntab
         ii = zone_begin + i - 1
         xtab(i) = 0.5d0*(chi(i)+chi(i-1))
         ytab(i) = am_diffusion_coeff(ii)
      end do
      ntabb = ntab + 1
      xtab(ntabb) = chi(ntab)
      ytab(ntabb) = am_diffusion_coeff(zone_end)
      xval(1) = chi(1)
      do i = 2, ntot
         xval(i) = echi(i)-0.5d0*dchi
      end do
      call osplin(xval,eq_am_diffusion_coeff,xtab,ytab,ntabb,ntot)
! DIFFUSION COEFFICIENT FOR MIXING - ASSUME CONSTANT BELOW
! BOTTOM INTERFACE OR ABOVE TOP INTERFACE
      ytab(1) = mixing_diffusion_coeff(zone_begin+1)
      do i = 2,ntab
         ii = zone_begin + i - 1
         ytab(i) = mixing_diffusion_coeff(ii)
      end do
      ytab(ntabb) = mixing_diffusion_coeff(zone_end)
      call osplin(xval,eq_mixing_diffusion_coeff,xtab,ytab,ntabb,ntot)
! ADD DIFFUSION PLUS ADVECTION TREATMENT IF DESIRED
      if (use_diffusion_advection_transport) then
         scale_factor = 0.2d0*c4pi*difad_velocity_scale
         ytab(1) = scale_factor*es_advective_velocity(zone_begin + 1)
         do i = 2,ntab
            ii = zone_begin + i - 1
            ytab(i) = scale_factor*es_advective_velocity(ii)
         end do
         ytab(ntabb) = scale_factor*es_advective_velocity(zone_end)
         call osplin(xval,am_advective_coeff,xtab,ytab,ntabb,ntot)
         ytab(1) = scale_factor*es_diffusive_velocity(zone_begin + 1)
         do i = 2,ntab
            ii = zone_begin + i - 1
            ytab(i) = scale_factor*es_diffusive_velocity(ii)
         end do
         ytab(ntabb) = scale_factor*es_diffusive_velocity(zone_end)
         call osplin(xval,am_diffusive_coeff,xtab,ytab,ntabb,ntot)
! MHP 05/02 NOTE THAT THE COMPOSITION DIFFUSION COEFFICIENTS
! SHOULD BE ADDED TO THE DIFFUSIVE PART OF THE DIFFUSION + ADVECTION
! ANGULAR MOMENTUM TRANSPORT; C.F. ZAHN 1992.  ESSENTIALLY THE
! ORIGINAL VESD TERM IS THE DIFFUSIVE COMPONENT FROM HORIZONTAL
! TRANSPORT WHILE THE COD2 TERM REPRESENTS THE DIFFUSIVE COMPONENT
! FROM VERTICAL TRANSPORT.
         do i = 1, ntot
            am_diffusive_coeff(i) = am_diffusive_coeff(i)+ &
                 eq_mixing_diffusion_coeff(i)
         end do
      end if
! PRODUCT OF RHO R^2 BY D CHI/DR
      mass_scale_factor = chi_grid_scale(2)
      luminosity_scale_factor = chi_grid_scale(9)*log_luminosity(num_zones)* &
           solar_luminosity_cgs
      pressure_scale_factor = chi_grid_scale(11)
      do i = 1, ntab
         ii = zone_begin + i - 1
         xtab(i) = chi(i)
! D CHI/DR = 1/DM*( D LOG M/DR) + 1/DL*(DL/DR) - 1/DP*(D LOG P/DR)
! OR, USING FAC = 4*PI*RHO*R**2
! D CHI/DR = FAC/(LN 10 * DM * M) + FAC*EPSILON/DL + RHO*GM/(LN10*DP*R**2)
! STORED IN YVAL
         scale_factor = c4pi*exp(ln10*(log_density(ii)+2.0d0*log_radius(ii)))
         dchi_dr = scale_factor/(ln10*mass_scale_factor*enclosed_mass(ii))+ &
              scale_factor*mix_phys%epsm(ii)/luminosity_scale_factor+ &
              exp(ln10*(cgl+log_density(ii)+log_mass(ii)-log_pressure(ii)- &
              2.0d0*log_radius(ii)))/(ln10*pressure_scale_factor)
         ytab(i) = log_density(ii) + log10(dchi_dr) + 2.0d0*log_radius(ii)
      end do
      call osplin(xval,yval,xtab,ytab,ntab,ntot)
! NOW ADD MULTIPLICATIVE FACTORS TO DIFFUSION COEFFICIENTS
! NOTE THAT A FACTOR OF 4PI HAS ALREADY BEEN INCLUDED IN
! CODIFF.
      do i = 1, ntot
         eq_mixing_diffusion_coeff(i) = eq_mixing_diffusion_coeff(i)* &
              exp(ln10*yval(i))
      end do
! PRODUCT OF RHO R^4 BY D CHI/DR - STORED IN YVAL
! MHP 05/02 ADDED FACTOR OF I/MR^2 - 2/3 FOR A SPHERICAL SHELL
      do i = 1, ntab
         ii = zone_begin + i - 1
!         YTAB(I) = YTAB(I) + 2.0D0*HR(II) + LOG10(EI0(I))
         ytab(i) = ytab(i)  + log10(moment_of_inertia(ii)/shell_mass(ii))
      end do
      call osplin(xval,yval,xtab,ytab,ntab,ntot)
! NOW ADD MULTIPLICATIVE FACTORS TO DIFFUSION COEFFICIENTS
      do i = 1, ntot
         eq_am_diffusion_coeff(i) = eq_am_diffusion_coeff(i)*exp(ln10*yval(i))
      end do
! MHP 05/02
      if (use_diffusion_advection_transport) then
         do i = 1,ntot
            am_advective_coeff(i) = am_advective_coeff(i)*exp(ln10*yval(i))
            am_diffusive_coeff(i) = am_diffusive_coeff(i)*exp(ln10*yval(i))
         end do
      end if
! REDEFINE DR AS DCHI
      grid_spacing = dchi
 9999 continue
      return
end subroutine rotgrid
