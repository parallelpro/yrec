!----------------------------------------------------------------------
! am_grid
!----------------------------------------------------------------------
! The equally spaced grid of angular momentum transport
! (am_transport_grid), and the map of its results back to the model
! (equal_grid_to_model).
module am_grid
      implicit none
      private
      public :: am_transport_grid, equal_grid_to_model

contains

! ---------------------------------------------------------------
! Builds the equally-spaced-grid quantities needed by am_advection_diffusion_coeffs.f90's
! angular-momentum-transport diffusion solve: given the first/last
! unstable zones of a region (zone_begin/zone_end), calls equal_spaced_grid to
! lay down the equally spaced coordinate rot_scr%chi, then computes the
! equally spaced grid moment of inertia (eq_moment_of_inertia),
! specific angular momentum (eq_angular_momentum), mass (eq_mass),
! angular velocity (eq_omega), and the geometrically weighted
! angular-momentum/mixing diffusion coefficients
! (eq_am_diffusion_coeff/eq_mixing_diffusion_coeff) at the zone edges,
! including the Jacobian factor for the transformation from radius to
! the rot_scr%chi coordinate. Companion routine to composition_grid.f90 (composition
! diffusion) for the rotation/angular-momentum-transport case.
subroutine am_transport_grid(am_diffusion_coeff, mixing_diffusion_coeff, log_density, &
     moment_of_inertia, specific_angular_momentum, luminosity_lsun, &
     log_pressure, log_radius, radius, log_mass, enclosed_mass, shell_mass, &
     log_total_mass, zone_begin, zone_end, am_transport_convective_flag, &
     num_zones, omega, grid_spacing, eq_am_diffusion_coeff, &
     eq_mixing_diffusion_coeff, eq_moment_of_inertia, eq_angular_momentum, &
     eq_mass, eq_omega, single_interface_flag)
      use rotation_scratch_lib
      use equal_grid
      use star_info_lib, only: star, json
      use physical_constants
      use numerics_lib
      use math_lib
      implicit none

! INPUT VARIABLES
      double precision, intent(in) :: am_diffusion_coeff(json), &
           mixing_diffusion_coeff(json), log_density(json), &
           moment_of_inertia(json), specific_angular_momentum(json), &
           luminosity_lsun(json), log_pressure(json), log_radius(json), &
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

      double precision :: eq_reduced_moment_of_inertia(json)
      integer :: ntab, i, ii, ntabb
      double precision :: emtop, embot, scale_factor
! the spline's table (the model points) and evaluation points (the grid)
      double precision :: table_x(json), table_y(json), eval_x(json), eval_y(json)

! FLAG THE SPECIAL CASE OF A SINGLE UNSTABLE INTERFACE AND EXIT
      if (zone_end-zone_begin.le.1) then
         single_interface_flag = .true.
         return
      else
         single_interface_flag = .false.
      end if
! DEFINE A GRID OF EQUALLY SPACED POINTS.
      call equal_spaced_grid(luminosity_lsun,log_pressure,log_mass,zone_begin, &
           zone_end,num_zones)
      call interpolate_to_grid_centers()
      call sum_zone_contents()
      call edge_coefficients()
      call apply_jacobians()
! REDEFINE DR AS DCHI
      grid_spacing = rot_scr%eq_chi_spacing

contains

! ---------------------------------------------------------------
subroutine interpolate_to_grid_centers()
! EQUAL_SPACED_GRID HAS DEFINED A SET OF CO-ORDINATES (CHI) AND EQUALLY
! SPACED MASS POINTS.  NOW FIND THE OTHER QUANTITIES OF INTEREST AT ZONE
! CENTERS:
! J/M (TO GET JTOT)
! I/MR^2 (TO GET ITOT)
! OMEGA
! R (NEEDED FOR I CALCULATION) - STORED IN YVAL
      ntab = zone_end - zone_begin + 1
      do i = 1,ntab
         ii = zone_begin + i - 1
         table_x(i) = rot_scr%chi(i)
         table_y(i) = log_radius(ii)
      end do
      do i = 1,rot_scr%num_eq_points
         eval_x(i) = rot_scr%eq_chi(i)
      end do
      call osplin(eval_x,eval_y,table_x,table_y,ntab,rot_scr%num_eq_points)
! I/MR^2
      do i = 1,ntab
         ii = zone_begin + i - 1
         table_y(i) = moment_of_inertia(ii)/(shell_mass(ii)*radius(ii)**2)
      end do
! MHP 05/02 STORE I/MR^2 FOR LATER USE IN
! DIFFUSION COEFFICIENTS FOR ANGULAR MOMENTUM
! TRANSPORT
      call osplin(eval_x,eq_reduced_moment_of_inertia,table_x,table_y,ntab,rot_scr%num_eq_points)
! J/M
      do i = 1,ntab
         ii = zone_begin + i - 1
         table_y(i) = specific_angular_momentum(ii)
      end do
      call osplin(eval_x,eq_angular_momentum,table_x,table_y,ntab,rot_scr%num_eq_points)
! OMEGA
      do i = 1,ntab
         ii = zone_begin + i - 1
         table_y(i) = omega(ii)
      end do
      call osplin(eval_x,eq_omega,table_x,table_y,ntab,rot_scr%num_eq_points)
end subroutine interpolate_to_grid_centers

! ---------------------------------------------------------------
subroutine sum_zone_contents()
! CONVERT TO TOTAL I AND TOTAL J OF THE ZONES.
! INTERMEDIATE POINTS
      do i = 2, rot_scr%num_eq_points-1
         eq_mass(i) = 0.5d0*(rot_scr%eq_enclosed_mass(i+1)-rot_scr%eq_enclosed_mass(i-1))
! MHP 05/02 CHANGED TO REFLECT THE FACT THAT
! THE INFO PREVIOUSLY STORED IN EI IS NOW IN EI0
         eq_moment_of_inertia(i) = eq_reduced_moment_of_inertia(i)* &
              eq_mass(i)*exp(ln10*2.0d0*eval_y(i))
         eq_angular_momentum(i) = eq_angular_momentum(i)*eq_mass(i)
      end do
! SPECIAL TREATMENT OF THE BOUNDARIES; CAN BE CONVECTIVE.
! IF CONVECTIVE SUM OVER ALL SHELLS.  CARE IS NEEDED TO DO BOOK-KEEPING
! PROPERLY AT THE EDGES - TOP IS HALFWAY TO EQUALLY SPACED POINT, NOT
! HALFWAY TO EDGE OF UNEQUALLY SPACED ORIGINAL SET OF POINTS.
!
! CENTER
      emtop = 0.5d0*(rot_scr%eq_enclosed_mass(2)+rot_scr%eq_enclosed_mass(1))
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
! (INCLUDE EVERY CONVECTIVE SHELL BELOW ZONE_BEGIN.)
      if (zone_begin.gt.1) then
         do ii = zone_begin-1,1,-1
            if (.not.am_transport_convective_flag(ii)) exit
            eq_mass(1) = eq_mass(1)+shell_mass(ii)
            eq_moment_of_inertia(1) = eq_moment_of_inertia(1)+ &
                 moment_of_inertia(ii)
            eq_angular_momentum(1) = eq_angular_momentum(1)+ &
                 specific_angular_momentum(ii)*shell_mass(ii)
         end do
      end if
! SURFACE
      embot = 0.5d0*(rot_scr%eq_enclosed_mass(rot_scr%num_eq_points)+rot_scr%eq_enclosed_mass(rot_scr%num_eq_points-1))
      if (zone_end.lt.num_zones) then
         emtop = 0.5d0*(enclosed_mass(zone_end)+enclosed_mass(zone_end+1))
      else
         emtop = exp(ln10*log_total_mass)
      end if
      eq_mass(rot_scr%num_eq_points) = emtop - embot
! MHP 05/02 CHANGED TO REFLECT THE FACT THAT
! THE INFO PREVIOUSLY STORED IN EI IS NOW IN EI0
      eq_moment_of_inertia(rot_scr%num_eq_points) = eq_reduced_moment_of_inertia(rot_scr%num_eq_points)* &
           eq_mass(rot_scr%num_eq_points)*radius(zone_end)**2
      eq_angular_momentum(rot_scr%num_eq_points) = eq_angular_momentum(rot_scr%num_eq_points)*eq_mass(rot_scr%num_eq_points)
! (INCLUDE EVERY CONVECTIVE SHELL ABOVE ZONE_END.)
      if (zone_end.lt.num_zones) then
         do ii = zone_end+1,num_zones
            if (.not.am_transport_convective_flag(ii)) exit
            eq_mass(rot_scr%num_eq_points) = eq_mass(rot_scr%num_eq_points)+shell_mass(ii)
            eq_moment_of_inertia(rot_scr%num_eq_points) = eq_moment_of_inertia(rot_scr%num_eq_points)+ &
                 moment_of_inertia(ii)
            eq_angular_momentum(rot_scr%num_eq_points) = eq_angular_momentum(rot_scr%num_eq_points)+ &
                 specific_angular_momentum(ii)*shell_mass(ii)
         end do
      end if
end subroutine sum_zone_contents

! ---------------------------------------------------------------
subroutine edge_coefficients()
! NOW SOLVE FOR QUANTITIES NEEDED AT THE ZONE EDGES.  THESE ARE
! RELATED TO THE DIFFUSION COEFFICIENTS.  UNLIKE THE EQUALLY SPACED
! GRID IN R, WE NEED TO INCLUDE A JACOBIAN TERM FOR THE TRANSFORMATION
! OF VARIABLES.
      call edge_grid_abscissae(rot_scr%chi, rot_scr%eq_chi, rot_scr%eq_chi_spacing, ntab, &
           rot_scr%num_eq_points, table_x, eval_x, ntabb)
! DIFFUSION COEFFICIENT FOR ANGULAR MOMENTUM - ASSUME CONSTANT BELOW
! BOTTOM INTERFACE OR ABOVE TOP INTERFACE
      call interp_edge_coeff(am_diffusion_coeff, zone_begin, zone_end, ntab, ntabb, &
           rot_scr%num_eq_points, eval_x, table_x, table_y, eq_am_diffusion_coeff)
! DIFFUSION COEFFICIENT FOR MIXING - ASSUME CONSTANT BELOW
! BOTTOM INTERFACE OR ABOVE TOP INTERFACE
      call interp_edge_coeff(mixing_diffusion_coeff, zone_begin, zone_end, ntab, ntabb, &
           rot_scr%num_eq_points, eval_x, table_x, table_y, eq_mixing_diffusion_coeff)
! ADD DIFFUSION PLUS ADVECTION TREATMENT IF DESIRED
      if (star%ctrl%use_diffusion_advection_transport) then
         scale_factor = 0.2d0*c4pi*star%ctrl%difad_velocity_scale
         table_y(1) = scale_factor*rot_scr%es_advective_velocity(zone_begin + 1)
         do i = 2,ntab
            ii = zone_begin + i - 1
            table_y(i) = scale_factor*rot_scr%es_advective_velocity(ii)
         end do
         table_y(ntabb) = scale_factor*rot_scr%es_advective_velocity(zone_end)
         call osplin(eval_x,rot_scr%am_advective_coeff,table_x,table_y,ntabb,rot_scr%num_eq_points)
         table_y(1) = scale_factor*rot_scr%es_diffusive_velocity(zone_begin + 1)
         do i = 2,ntab
            ii = zone_begin + i - 1
            table_y(i) = scale_factor*rot_scr%es_diffusive_velocity(ii)
         end do
         table_y(ntabb) = scale_factor*rot_scr%es_diffusive_velocity(zone_end)
         call osplin(eval_x,rot_scr%am_diffusive_coeff,table_x,table_y,ntabb,rot_scr%num_eq_points)
! MHP 05/02 NOTE THAT THE COMPOSITION DIFFUSION COEFFICIENTS
! SHOULD BE ADDED TO THE DIFFUSIVE PART OF THE DIFFUSION + ADVECTION
! ANGULAR MOMENTUM TRANSPORT; C.F. ZAHN 1992.  ESSENTIALLY THE
! ORIGINAL VESD TERM IS THE DIFFUSIVE COMPONENT FROM HORIZONTAL
! TRANSPORT WHILE THE COD2 TERM REPRESENTS THE DIFFUSIVE COMPONENT
! FROM VERTICAL TRANSPORT.
         do i = 1, rot_scr%num_eq_points
            rot_scr%am_diffusive_coeff(i) = rot_scr%am_diffusive_coeff(i)+ &
                 eq_mixing_diffusion_coeff(i)
         end do
      end if
end subroutine edge_coefficients

! ---------------------------------------------------------------
subroutine apply_jacobians()
! PRODUCT OF RHO R^2 BY D CHI/DR (LOG10, ON THE EQUAL GRID IN YVAL;
! THE MODEL-POINT TABLE STAYS IN YTAB FOR THE RHO R^4 STEP BELOW)
      call dchi_dr_jacobian(log_density, log_radius, log_mass, log_pressure, &
           enclosed_mass, mid_scr%eps_plus_neu_grav, luminosity_lsun(num_zones), zone_begin, &
           ntab, rot_scr%num_eq_points, rot_scr%chi, eval_x, table_x, &
           table_y, eval_y)
! NOW ADD MULTIPLICATIVE FACTORS TO DIFFUSION COEFFICIENTS
! NOTE THAT A FACTOR OF 4PI HAS ALREADY BEEN INCLUDED IN
! diffusion_velocity_scales.
      call multiply_by_exp10(eq_mixing_diffusion_coeff, eval_y, rot_scr%num_eq_points)
! PRODUCT OF RHO R^4 BY D CHI/DR - STORED IN YVAL
! MHP 05/02 ADDED FACTOR OF I/MR^2 - 2/3 FOR A SPHERICAL SHELL
      do i = 1, ntab
         ii = zone_begin + i - 1
         table_y(i) = table_y(i)  + log10(moment_of_inertia(ii)/shell_mass(ii))
      end do
      call osplin(eval_x,eval_y,table_x,table_y,ntab,rot_scr%num_eq_points)
! NOW ADD MULTIPLICATIVE FACTORS TO DIFFUSION COEFFICIENTS
      call multiply_by_exp10(eq_am_diffusion_coeff, eval_y, rot_scr%num_eq_points)
! MHP 05/02
      if (star%ctrl%use_diffusion_advection_transport) then
         call multiply_by_exp10(rot_scr%am_advective_coeff, eval_y, rot_scr%num_eq_points)
         call multiply_by_exp10(rot_scr%am_diffusive_coeff, eval_y, rot_scr%num_eq_points)
      end if
end subroutine apply_jacobians

end subroutine am_transport_grid

! ---------------------------------------------------------------
! The diffusion equation for angular momentum transport is solved using
! an equally spaced grid for each unstable region. Once the equation
! has been solved, we need the new run of angular momentum at the
! original model points.
!
! eq2mod transforms back to the original grid and calculates the
! change in specific angular momentum in each shell.
!
! Input variables:
!
!   Variables pertaining to the equally spaced grid points used to
!   solve the diffusion equation:
!
!   delta_angular_momentum - change in J at each point due to
!       rotationally induced instabilities.
!   angular_momentum - contains the run of angular momentum prior to
!       angular momentum redistribution.
!   rot_scr%num_eq_points (common/egrid/) - number of equally spaced grid points in the
!       unstable region.
!   total_delta_angular_momentum - total change in angular momentum
!       across grid.
!
!   Variables pertaining to the original model points:
!
!   shell_mass - the unlogged mass content of the model points.
!   zone_begin, zone_end - the first and last model points contained
!       in the unstable region.
!   convective_flag - a logical array which is T if the model point is
!       convective for angular momentum purposes and F otherwise.
!   num_points - the total number of model points.
!
! Scratch (rotation_scratch_lib): rot_scr%chi, eq_chi and num_eq_points
! hold the equally spaced grid.
!
! Output variable:
!   specific_angular_momentum - the new run of specific angular
!       momentum at the model points.
!-------------------------------------------------------------------------
subroutine equal_grid_to_model(delta_angular_momentum, angular_momentum, shell_mass, &
     zone_begin, zone_end, convective_flag, num_points, &
     total_delta_angular_momentum, specific_angular_momentum)
      use rotation_scratch_lib
      use star_info_lib, only: json
      use numerics_lib
      implicit none

      double precision, intent(in) :: delta_angular_momentum(json), &
           angular_momentum(json), shell_mass(json)
      integer, intent(in) :: zone_begin, zone_end
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_points
      double precision, intent(in) :: total_delta_angular_momentum
      double precision, intent(inout) :: specific_angular_momentum(json)
      integer :: i0, i1, i, nmod, ii
      double precision :: sumjmod, sumdjmod, test, val, ratio
! the spline's table (the grid) and evaluation points (the model points)
      double precision :: table_x(json), table_y(json), eval_x(json), eval_y(json)

! CHECK FOR LOWER AND UPPER CONVECTION ZONES.
      if (.not.convective_flag(zone_begin) .or. zone_begin.eq.1) then
         i0 = zone_begin
      else
         do i = zone_begin - 1, 1, -1
            if (.not.convective_flag(i)) then
               i0 = i + 1
               exit
            end if
         end do
         if (i < (1)) then
         i0 = 1
         end if
      end if
      if (.not.convective_flag(zone_end) .or. zone_end.eq.num_points) then
         i1 = zone_end
      else
         do i = zone_end+1, num_points
            if (.not.convective_flag(i)) then
               i1 = i - 1
               exit
            end if
         end do
         if (i > num_points) then
         i1 = num_points
         end if
      end if
! INTERPOLATE IN DJ/J AS A FUNCTION OF ECHI
      nmod = zone_end - zone_begin + 1
      do i = 1, nmod
         eval_x(i) = rot_scr%chi(i)
      end do
      do i = 1, rot_scr%num_eq_points
         table_x(i) = rot_scr%eq_chi(i)
         table_y(i) = delta_angular_momentum(i)/angular_momentum(i)
      end do
      call osplin(eval_x, eval_y, table_x, table_y, rot_scr%num_eq_points, nmod)
! APPLY THE FRACTIONAL CHANGE IN J/M ACROSS THE MODEL.
      sumjmod = 0.0d0
      sumdjmod = 0.0d0
      if (i0.lt.zone_begin) then
         do i = i0, zone_begin-1
            sumjmod = sumjmod + specific_angular_momentum(i)*shell_mass(i)
            sumdjmod = sumdjmod + eval_y(1)*specific_angular_momentum(i)* &
                 shell_mass(i)
            specific_angular_momentum(i) = specific_angular_momentum(i)* &
                 (1.0d0+eval_y(1))
         end do
      end if
      do i = zone_begin, zone_end
         ii = i - zone_begin + 1
         sumjmod = sumjmod + specific_angular_momentum(i)*shell_mass(i)
         sumdjmod = sumdjmod + eval_y(ii)*specific_angular_momentum(i)* &
              shell_mass(i)
         specific_angular_momentum(i) = specific_angular_momentum(i)* &
              (1.0d0+eval_y(ii))
      end do
      if (i1.gt.zone_end) then
         do i = zone_end+1, i1
            sumjmod = sumjmod + specific_angular_momentum(i)*shell_mass(i)
            sumdjmod = sumdjmod + eval_y(nmod)*specific_angular_momentum(i)* &
                 shell_mass(i)
            specific_angular_momentum(i) = specific_angular_momentum(i)* &
                 (1.0d0+eval_y(nmod))
         end do
      end if
! CHECK THAT ANGULAR MOMENTUM HAS BEEN CONSERVED.  THE TOTAL ANGULAR
! MOMENTUM CAN BE REDUCED FROM SURFACE ANGULAR MOMENTUM LOSS; CARE IS
! USED HERE BECAUSE DUMDJ SHOULD BE ZERO IN MODELS WITHOUT LOSS.
      test = abs(sumdjmod - total_delta_angular_momentum)
      val = 1.0d-5*abs(total_delta_angular_momentum)
! IF THE TEST SUM IS NOT THE SAME WITHIN THE RELATIVE TOLERANCE VAL,
! ADJUST THE TOTAL ANGULAR MOMENTUM OF THE ENTIRE REGION BY A CONSTANT
! FACTOR TO ENFORCE ANGULAR MOMENTUM CONSERVATION.
      if (test.gt.val) then
         ratio = (sumjmod + total_delta_angular_momentum)/ &
              (sumjmod + sumdjmod)
         do i = i0, i1
            specific_angular_momentum(i) = specific_angular_momentum(i)*ratio
         end do
      end if
      return
end subroutine equal_grid_to_model

end module am_grid
