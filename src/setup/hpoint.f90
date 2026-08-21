!----------------------------------------------------------------------
! hpoint
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original hpoint.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Rezones the model, spacing the shells such that the maximum distance
! between two successive points in P, L, X, and Z specified by the
! user is not exceeded. Also flags real discontinuities to prevent
! artificial smoothing, relocates the convective core/envelope
! boundaries in the new point distribution, and (for rotating models)
! recomputes the moment of inertia and angular velocity distribution
! consistent with the new mesh.
subroutine hpoint(num_zones,log_total_mass,log_mass,enclosed_mass, &
     shell_mass,log_temperature,log_pressure,log_radius,convective_flag, &
     log_luminosity,log_density,composition,stored_log_pressure, &
     stored_log_radius,stored_envelope_values,envelope_store_index, &
     point_reset_flag,h_shell_zone_begin,model_number, &
     h_shell_active,convective_core_edge_zone,convective_envelope_edge_zone, &
     omega,eta_squared,mean_radius,moment_of_inertia, &
     specific_angular_momentum,rotational_kinetic_energy, &
     total_angular_momentum,total_rotational_ke, &
! BL,DELTS,FP,FT,HG,QIW,SMASS,TEFFL)  ! KC 2025-05-31
     fp,ft,hg,qiw,log_teff)

      use light_burn_lib
      use scrtch_lib
      use oldmod_lib
      use luout_lib
      use const_lib
      use numerics_lib
      implicit none
      integer, parameter :: json = 5000

      integer, intent(inout) :: num_zones
      double precision, intent(in) :: log_total_mass
      double precision, intent(inout) :: log_mass(json)
      double precision, intent(out) :: enclosed_mass(json), shell_mass(json)
      double precision, intent(inout) :: log_temperature(json), &
           log_pressure(json), log_radius(json)
      logical, intent(inout) :: convective_flag(json)
      double precision, intent(inout) :: log_luminosity(json), &
           log_density(json)
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: stored_log_pressure(3), &
           stored_log_radius(3), stored_envelope_values(4)
      integer, intent(in) :: envelope_store_index
      logical, intent(out) :: point_reset_flag
      integer, intent(in) :: h_shell_zone_begin, model_number
      logical, intent(in) :: h_shell_active
      integer, intent(inout) :: convective_core_edge_zone, &
           convective_envelope_edge_zone
      double precision, intent(inout) :: omega(json), eta_squared(json), &
           mean_radius(json), moment_of_inertia(json), &
           specific_angular_momentum(json)
      double precision, intent(out) :: rotational_kinetic_energy(json)
      double precision, intent(inout) :: total_angular_momentum, &
           total_rotational_ke
      double precision, intent(inout) :: fp(json), ft(json), hg(json), &
           qiw(json)
      double precision, intent(in) :: log_teff

! MHP 10/02 added MRZONE,MXZONE to dimension statements
      double precision :: log10_omega(json), point_spacing_max(4)
      integer :: flag_point(100)
      logical :: am_transport_convective_flag(json)
      double precision :: ft_old(json), fp_old(json)
      integer :: radiative_zone_bounds(13,2), convective_zone_bounds(12,2)


! common/burn/: reaction_rate_by_zone (originally HCOMPM) is read here
! (interpolated to the new mesh alongside HCOMP). Naming matches
! mix.f90.
      double precision :: reaction_rate_by_zone(15,json)
      common/burn/ reaction_rate_by_zone

! common/ccout/: not used in this file. Naming matches ccoeft.f90.
      logical :: lstore, lstatm, lstenv, lstmod, lstphys, lstrot, lscrib, &
           lstch, lphhd
      common/ccout/ lstore, lstatm, lstenv, lstmod, lstphys, lstrot, &
           lscrib, lstch, lphhd

! common/ccout1/: only npoint is used here. Naming matches wrtmil.f90.
      integer :: npenv, nprtmod, print_point_interval, npoint
      common/ccout1/ npenv, nprtmod, print_point_interval, npoint

! common/ccout2/: only ldebug is used here. Naming matches meqos.f90.
      logical :: ldebug, lcorr, lmilne, ltrack, lstpch
      common/ccout2/ ldebug, lcorr, lmilne, ltrack, lstpch

! common/comp/: senv/xnew/znew are used here. Naming matches
! getopac.f90.
      double precision :: envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, amuenv, fxenv(12), xnew, znew, stotal, senv
      common/comp/ envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, amuenv, fxenv, xnew, znew, stotal, senv




! common/ctol/: chi_grid_scale (originally HPTTOL) is used here for the
! mesh-spacing tolerances. Naming matches mixgrid.f90.
      double precision :: htoler(5,2), fcorr0, fcorri, fcorr, &
           chi_grid_scale(12)
      integer :: niter1, niter2, niter3
      common/ctol/ htoler, fcorr0, fcorri, fcorr, chi_grid_scale, niter1, &
           niter2, niter3


! common/dwmax/: max shear (domega/dr) magnitude, current and from the
! previous timestep; only max_domega_dr is set here. Not referenced in
! any already-converted file.
      double precision :: max_domega_dr(json), max_domega_dr_old(json)
      common/dwmax/ max_domega_dr, max_domega_dr_old

! common/flag/: use_extended_composition (originally LEXCOM). Naming
! matches mixcz.f90.
      logical :: use_extended_composition
      common/flag/ use_extended_composition


! common/oldrot/: previous-timestep rotation-state snapshot, all used
! here. Naming is local to this batch (matches midmod.f90).
      double precision :: old_omega(json), old_specific_angular_momentum(json), &
           old_moment_of_inertia(json), old_hg(json), old_mean_radius(json), &
           old_eta_squared(json)
      common/oldrot/ old_omega, old_specific_angular_momentum, &
           old_moment_of_inertia, old_hg, old_mean_radius, old_eta_squared

! common/oldphy/: previous-timestep auxiliary physics quantities,
! stored here for use by the diffusion/mixing routines at the start of
! the next timestep. Naming is local to this batch (matches
! midmod.f90); member names parallel common/mdphy/'s amum/cpm/delm/
! delam/delrm/esumm/om/qdtm/thdifm/velm/viscm/epsm.
      double precision :: old_delm(json), old_del_adiabatic_mix(json), &
           old_amu(json), old_om(json), old_cp(json), old_qdt(json), &
           old_vel(json), old_visc(json), old_thdif(json), old_esum(json), &
           old_del_radiative_mix(json), old_eps(json)
      common/oldphy/ old_delm, old_del_adiabatic_mix, old_amu, old_om, &
           old_cp, old_qdt, old_vel, old_visc, old_thdif, old_esum, &
           old_del_radiative_mix, old_eps

! MHP 06/02
! Time change of theta
! common/oldrot2/: tho/qwrst (domega/dr and its exponentially-weighted
! form) are used here; the remaining members are unused placeholders.
! Naming matches vcirc.f90.
      double precision :: tho(json), theta_new(json), theta_mean(json), &
           del_grad_diff_interface(json), es_relaxation_factor(json), &
           theta_prev(json), qwrst(json), wmst(json), qwrmst(json)
      common/oldrot2/ tho, theta_new, theta_mean, del_grad_diff_interface, &
           es_relaxation_factor, theta_prev, qwrst, wmst, qwrmst


! common/optab/: not used in this file. Naming matches getopac.f90.
      double precision :: metal_fraction_match_tolerance, zsi
      integer :: idt, idd(4)
      common/optab/ metal_fraction_match_tolerance, zsi, idt, idd



! common/temp/: current-timestep counterparts of common/oldphy/'s
! old_cp/old_amu/old_qdt/old_thdif/old_visc, all used here. Not
! referenced in any already-converted file.
      double precision :: cp(json), mean_molecular_weight(json), qdt(json), &
           thdif(json), visc(json)
      common/temp/ cp, mean_molecular_weight, qdt, thdif, visc

! 7/91 entropy term common block added.
! common/entrop/: entropy-correction terms for T/P/L/R, all used here
! (interpolated onto the new mesh like the other physical variables).
! Not referenced in any already-converted file.
      double precision :: temperature_entropy_term(json), &
           pressure_entropy_term(json), luminosity_entropy_term(json), &
           radius_entropy_term(json)
      common/entrop/ temperature_entropy_term, pressure_entropy_term, &
           luminosity_entropy_term, radius_entropy_term


      integer :: reaction_rate_species_index(7)
      double precision :: z_new(json), x_new(json)
! MHP 6/00 added dummy vector
! MHP 7/02 added chi vector - a running total of the normalized
! differences in M, L, P between shells that is used to set the
! point spacing.
      double precision :: spline_x(json), spline_y(json), chi(json), &
           spline_second_deriv(json)
      integer :: gradient_flag_index(json)
!       DATA IDMAX/JSON/
      data reaction_rate_species_index/1,2,4,5,6,7,9/
      save

      integer :: num_species_tracked, i, j, k
      integer :: flag_count
      double precision :: pressure_test, dp_scale
      integer :: new_num_zones, working_num_zones, num_new_points
      double precision :: delta_x_over_max, delta_z_over_max, &
           point_insert_spacing
      double precision :: mass_scale, luminosity_scale
      double precision :: dchi
      double precision :: chi_prev
      integer :: point_prev_index
      double precision :: delta_chi
      integer :: segment_point_count
      double precision :: spline_eval_x, spline_eval_y
      integer :: spline_klo, spline_khi
      integer :: gradient_flag_count
      double precision :: luminosity_max
      double precision :: log_omega_top, log_omega_bot
      integer :: overshoot_base_zone, fine_zone_base
      double precision :: delta_log_pressure
      double precision :: pmax1, pmax2, pmax3, pmax4, pmax5
      integer :: overshoot_point_count
      logical :: sort_done
      integer :: temp_swap
      integer :: point_insert_flag, chi_start_index
      integer :: old_point_count, new_point_count
      integer :: zone_index
      integer :: min_common_count
      double precision :: sum_angular_momentum, sum_rotational_ke
      double precision :: angular_momentum_shell
      double precision :: delta_radius, delta_omega, omega_mid, log_factor
      double precision :: mass_two_back, mass_prev, mass_curr
      integer :: radiative_zone_begin, radiative_zone_end
      integer :: num_radiative_zones, num_convective_zones

!  THIS SUBROUTINE REZONES THE MODEL, SPACING THE SHELLS SUCH THAT
!  THE MAXIMUM DISTANCE BETWEEN 2 SUCCESSIVE POINTS IN P,L,X, AND Z
!  SPECIFIED BY THE USER IS NOT EXCEEDED.  IT ALSO FLAGS REAL
!  DISCONTINUITIES TO PREVENT ARTIFICIAL SMOOTHING.
      point_reset_flag = .false.
!  IEND IS THE NUMBER OF SPECIES THE PROGRAM IS KEEPING TRACK OF
      num_species_tracked = 11
      if (use_extended_composition) num_species_tracked = 15
! CHECK IF TEMPERATURE OF OUTERMOST HENYEY POINT > MIMIMUM ENVELOPE T
      if (log_temperature(num_zones).lt.tenv0) then
       do 20 i = num_zones-1,1,-1
          if (log_temperature(i).gt.tenv0) then
             write(short_file_unit,10) num_zones,i
   10     format(' OUTER POINTS DELETED OLD M =',I5,'  NEW M =',I5)
             num_zones = i
             senv = log_mass(num_zones) - log_total_mass
             point_reset_flag = .true.
             goto 40
          endif
   20    continue
!  ENTIRE MODEL HAS T<TENV0 - UNLIKELY - BUT STOP IF TRUE
       write(short_file_unit,30)
       write(iowr,30)
   30    format(1X,39('>'),40('<')/1X,'ERROR IN HPOINT'/ &
     1X,'MAX. STAR T LESS THAN MINIMUM ENVELOPE T.RUN STOPPED')
       stop
   40    continue
!  CHECK IF OUTER POINT T < MAXIMUM ENVELOPE T
      else if (log_temperature(num_zones).gt.tenv1.and. &
           envelope_store_index.ne.0) then
       num_zones = num_zones + 1
       senv = stored_envelope_values(4)
       log_mass(num_zones) = log_total_mass + senv
       log_pressure(num_zones) = log_pressure(num_zones-1) + &
            (stored_envelope_values(1) - stored_log_pressure(envelope_store_index))
       log_temperature(num_zones) = log_temperature(num_zones-1) + &
            (stored_envelope_values(1) - &
            stored_log_pressure(envelope_store_index))*0.250D0
       log_radius(num_zones) = log_radius(num_zones-1) + &
            (stored_envelope_values(3) - stored_log_radius(envelope_store_index))
       log_luminosity(num_zones) = log_luminosity(num_zones-1)
       do 50 i = 1,num_species_tracked
          composition(i,num_zones) = composition(i,num_zones-1)
   50    continue
       log_density(num_zones) = log_pressure(num_zones) - &
            log_temperature(num_zones) - 8.0D0
       j = num_zones - 1
       write(short_file_unit,60) j,log_mass(j),log_pressure(j), &
            log_temperature(j),log_radius(j),num_zones,log_mass(num_zones), &
            log_pressure(num_zones),log_temperature(num_zones), &
            log_radius(num_zones)
   60    format(' OUTER POINT ADDED',I5,F15.10,'  PTR',3F10.6)
       point_reset_flag = .true.
      endif


! SET UP FLAGGED POINTS - PROGRAM WILL NOT REZONE ACROSS FLAGGED POINTS
      flag_count = 1
! FLAG EDGES OF CENTRAL AND SURFACE CONVECTION ZONES
      if (convective_core_edge_zone.gt.1) then
       flag_point(flag_count) = convective_core_edge_zone
       flag_count = flag_count + 1
      endif
      if (convective_envelope_edge_zone.lt.num_zones .and. &
           convective_envelope_edge_zone.gt.1) then
       flag_point(flag_count) = convective_envelope_edge_zone
       flag_count = flag_count + 1
      endif
! FLAG EDGE OF H SHELL
      if (h_shell_active.and.h_shell_zone_begin.gt.1) then
       flag_point(flag_count) = h_shell_zone_begin - 1
       flag_count = flag_count + 1
      endif
      do 120 i = 2,num_zones
! TEST FOR FLAGGING DUE TO X GRADIENT
       if (dabs(composition(1,i)-composition(1,i-1)).gt.chi_grid_scale(3)) then
          flag_point(flag_count) = i
          flag_count = flag_count + 1
! TEST FOR FLAGGING DUE TO Z GRADIENT
       else if (dabs(composition(3,i)-composition(3,i-1)).gt. &
            chi_grid_scale(4)) then
          flag_point(flag_count) = i
          flag_count = flag_count + 1
! TEST FOR FLAGGING DUE TO GRADIENT IN LOG OMEGA.
       else if (rotation_active) then
          log_omega_top = dlog10(omega(i))
          log_omega_bot = dlog10(omega(i-1))
          if (dabs(log_omega_top-log_omega_bot).gt.chi_grid_scale(12)) then
             flag_point(flag_count) = i
             flag_count = flag_count + 1
          endif
       endif
       if (flag_count.ge.100) then
          write(short_file_unit,110)
  110       format(1X,'MORE THAN 100 FLAG POINTS-FIRST 100 RETAINED')
          goto 130
       endif
  120 continue
!  PMAX1 = MAX DEL LOG P BELOW SURFACE C.Z. AND BELOW FINELY ZONED
!  REGION AROUND IT.
!  PMAX2 = MAX DEL LOG P BETWEEN LOWER EDGE OF FINELY ZONED REGION
!  AROUND SURFACE C.Z. AND BASE OF OVERSHOOT REGION.
!  PMAX3 = SAME FOR OVERSHOOT REGION.
!  PMAX4 = MAX DEL LOG P ABOVE BASE OF SURFACE C.Z. IN FINELY ZONED
!  REGION AROUND IT.
!  PMAX5 = SAME FOR THE OUTER POINTS IN THE STAR.
      pmax1 = chi_grid_scale(11)
      pmax4 = chi_grid_scale(10)
      pmax5 = chi_grid_scale(8)
      if (.not.convective_flag(num_zones)) then
       overshoot_base_zone = num_zones
       fine_zone_base = num_zones
      else if (convective_envelope_edge_zone.eq.1) then
       overshoot_base_zone = 1
       fine_zone_base = 1
      else
! LOCATE BASE OF OVERSHOOT REGION IF APPLICABLE.
       if (.not.envelope_overshoot_active) then
          i = convective_envelope_edge_zone
          overshoot_base_zone = convective_envelope_edge_zone
       else
          do 191 overshoot_base_zone = convective_envelope_edge_zone-1,1,-1
             if (log_pressure(overshoot_base_zone)- &
                  log_pressure(convective_envelope_edge_zone).gt.alphae) goto 193
  191       continue
  193       overshoot_base_zone = overshoot_base_zone + 1
          delta_log_pressure = log_pressure(overshoot_base_zone)- &
               log_pressure(convective_envelope_edge_zone)
          if (delta_log_pressure.gt.0.0D0) then
             overshoot_point_count = int(delta_log_pressure/chi_grid_scale(10))
             if (mod(delta_log_pressure,chi_grid_scale(10)).ne.0D0) &
                  overshoot_point_count = overshoot_point_count+1
             pmax3 = delta_log_pressure/dfloat(overshoot_point_count)
          else
             pmax3 = chi_grid_scale(10)
          endif
            if (overshoot_base_zone.gt.1) then
               flag_point(flag_count) = overshoot_base_zone
               flag_count = flag_count + 1
            endif
       endif
       if (chi_grid_scale(7).eq.0.0D0) then
          fine_zone_base = overshoot_base_zone
       else
! NOW LOCATE BASE OF FINELY ZONED REGION.
          if (overshoot_base_zone.eq.1) then
             fine_zone_base = 1
             goto 198
          endif
          do 195 fine_zone_base = overshoot_base_zone-1,1,-1
             if (log_pressure(fine_zone_base) - &
                  log_pressure(convective_envelope_edge_zone).gt.chi_grid_scale(7)) &
                  goto 197
  195       continue
  197       fine_zone_base = fine_zone_base + 1
          if (fine_zone_base.eq.overshoot_base_zone) goto 198
          delta_log_pressure = log_pressure(fine_zone_base) - &
               log_pressure(overshoot_base_zone)
          overshoot_point_count = int(delta_log_pressure/chi_grid_scale(10))
          if (mod(delta_log_pressure,chi_grid_scale(10)).ne.0D0) &
               overshoot_point_count = overshoot_point_count+1
          pmax2 = delta_log_pressure/dfloat(overshoot_point_count)
          if (pmax2.eq.0.0D0) pmax2 = chi_grid_scale(10)
            if (fine_zone_base.gt.1) then
               flag_point(flag_count) = fine_zone_base
               flag_count = flag_count + 1
            endif
  198       continue
       endif
      endif
  130 flag_point(flag_count) = num_zones
! ARRANGE THE FLAG POINTS IN ASCENDING ORDER
      if (flag_count.eq.1) goto 180
  140 continue
      sort_done = .true.
      do 150 i = 1,flag_count-1
       if (flag_point(i+1).lt.flag_point(i)) then
          temp_swap = flag_point(i)
          flag_point(i) = flag_point(i+1)
          flag_point(i+1) = temp_swap
          sort_done = .false.
       endif
  150 continue
      if (.not.sort_done) goto 140
! ENSURE THAT POINTS ARENT FLAGGED MORE THAN ONCE.
      i = 2
  160 continue
      if (flag_point(i).eq.flag_point(i-1)) then
       if (i.lt.flag_count) then
          do 170 j = i,flag_count-1
             flag_point(j) = flag_point(j+1)
  170       continue
       endif
       flag_count = flag_count - 1
      endif
      i = i + 1
      if (i.le.flag_count) goto 160
  180 continue
      write(short_file_unit,185) (flag_point(j),j=1,flag_count)
  185 format(1X,'FLAG-POINTS',20I4)
! BEGIN REFLOATING OF POINTS
      if (rotation_active) then
       do 190 i = 1,num_zones
          if (omega(i).gt.0.0D0) then
             log10_omega(i) = dlog10(omega(i))
          else
             log10_omega(i) = 0.0D0
          endif
  190    continue
      endif
      prev_model%old_shell_mass(1) = log_mass(1)
      prev_model%old_pressure(1) = log_pressure(1)
      prev_model%old_luminosity(1) = log_luminosity(1)
      x_new(1) = composition(1,1)
      z_new(1) = composition(3,1)
      luminosity_max = log_luminosity(num_zones)
!       JVS 04/14 added Teff to saved variables
        prev_model%old_teff = log_teff
!  JVS 05/25 Added model number to list of saved values
      prev_model%old_num_zones = num_zones
      do i = num_zones-1,1,-1
         if (log_luminosity(i).gt.luminosity_max) then
            luminosity_max = log_luminosity(i)
         endif
      end do
      point_spacing_max(1) = chi_grid_scale(8)
      point_spacing_max(2) = chi_grid_scale(9)*luminosity_max
      point_spacing_max(3) = chi_grid_scale(5)
      point_spacing_max(4) = chi_grid_scale(6)
!      KFACT = 0
!  200 CONTINUE
      point_insert_flag = 1
      j = 2
      chi_start_index = 2
! CHI IS THE NORMALIZED VECTOR OF DIFFERENCES IN M,L,P:
! CHI = HS/DELTA M + HL/DELTA L - HP/DELTA P
      mass_scale = chi_grid_scale(2)
      luminosity_scale = point_spacing_max(2)
      chi(1) = 1.0D0
      do j = 2, num_zones
         pressure_test = log_pressure(j) - log_pressure(convective_envelope_edge_zone)
         if (abs(pressure_test).lt.chi_grid_scale(7)) then
! FINELY ZONED REGION
            dp_scale = chi_grid_scale(10)
         else if (pressure_test.gt.chi_grid_scale(7)) then
! BELOW SURFACE CZ
            dp_scale = chi_grid_scale(11)
         else
! IN SURFACE CZ
            dp_scale = chi_grid_scale(8)
         endif
         if (log_luminosity(j).gt.log_luminosity(j-1)) then
            dchi = (log_mass(j)-log_mass(j-1))/mass_scale + &
                 (log_luminosity(j)-log_luminosity(j-1))/luminosity_scale - &
                 (log_pressure(j)-log_pressure(j-1))/dp_scale
         else
            dchi = (log_mass(j)-log_mass(j-1))/mass_scale - &
                 (log_pressure(j)-log_pressure(j-1))/dp_scale
         endif
         chi(j) = chi(j-1)+dchi
      end do
      do j = 1,num_zones
         spline_x(j) = chi(j)
         spline_y(j) = log_mass(j)
      end do
! GET SPLINE COEFFICIENTS
      call splinc(spline_x,spline_y,spline_second_deriv,num_zones)
      chi_prev = chi(1)
      point_prev_index = 1
      new_num_zones = 1
      do i = 1,flag_count
         delta_chi = chi(flag_point(i))-chi_prev
         segment_point_count = int(delta_chi)+1
         if (new_num_zones+segment_point_count.gt.json) then
! EXCEEDED ARRAY DIMENSIONS; COMPLAIN AND STOP
            write(*,101) new_num_zones+segment_point_count,json
 101        format(' DESIRED NUMBER OF POINTS ',I6,'EXCEEEDS JSON ', &
     I6/' RUN STOPPED')
            stop 9999
         endif
         dchi = delta_chi/dfloat(segment_point_count)
! ASSIGN NEW POINTS
         do j = point_prev_index+1,point_prev_index+segment_point_count
            spline_eval_x = chi_prev + dchi
            call splintd2(spline_x, spline_y, num_zones, &
                 spline_second_deriv, spline_eval_x, spline_eval_y, &
                 spline_klo, spline_khi)
            prev_model%old_shell_mass(j) = spline_eval_y
            chi_prev = spline_eval_x
!

         end do
         new_num_zones = new_num_zones + segment_point_count
         point_prev_index = new_num_zones
!
      end do








! TEST FOR ASSIGNING POINTS BASED ON THE GRADIENT IN X.
      do j = 1,num_zones
         spline_x(j) = log_mass(j)
         spline_y(j) = composition(1,j)
      end do
! GET SPLINE COEFFICIENTS
      call splinc(spline_x,spline_y,spline_second_deriv,num_zones)
! ASSIGN INTERPOLATED VECTOR OF X VALUES TO HIO
      do i = 2,new_num_zones
         spline_eval_x = prev_model%old_shell_mass(i)
         call splintd2(spline_x, spline_y, num_zones, spline_second_deriv, &
              spline_eval_x, spline_eval_y, spline_klo, spline_khi)
         x_new(i) = spline_eval_y
      end do




! SKIP IF HPMAX(3) IS ZEROED OUT
      if (point_spacing_max(3).le.1.0D-15) goto 102
! TEST ON X-CHANGE (ONLY FOR INCREASING X) USING HIO AS DUMMY ARRAY
      gradient_flag_count = 0
      do j = new_num_zones,2,-1
         if (x_new(j)-x_new(j-1).gt.point_spacing_max(3)) then
            gradient_flag_count = gradient_flag_count + 1
            gradient_flag_index(gradient_flag_count) = j
         endif
      end do
      if (gradient_flag_count.eq.0) goto 102
      working_num_zones = new_num_zones
      do i = 1,gradient_flag_count
         j = gradient_flag_index(i)
         delta_x_over_max = (x_new(j) - x_new(j-1))/point_spacing_max(3)
! NUMBER OF NEW POINTS NEEDED
         num_new_points = int(delta_x_over_max)
!
         point_insert_spacing = (prev_model%old_shell_mass(j)-prev_model%old_shell_mass(j-1))/ &
              dfloat(num_new_points+1)
         do k = working_num_zones+num_new_points,j+num_new_points,-1
            prev_model%old_shell_mass(k) = prev_model%old_shell_mass(k-num_new_points)
!
         end do
         do k = j + num_new_points -1, j -1, -1
            prev_model%old_shell_mass(k) = prev_model%old_shell_mass(k+1) - point_insert_spacing
         end do
         working_num_zones = working_num_zones + num_new_points
      end do
      new_num_zones = working_num_zones
!
 102  continue
! TEST FOR ASSIGNING POINTS BASED ON THE GRADIENT IN Z.
      do j = 1,num_zones
         spline_x(j) = log_mass(j)
         spline_y(j) = composition(3,j)
      end do
! GET SPLINE COEFFICIENTS
      call splinc(spline_x,spline_y,spline_second_deriv,num_zones)
! ASSIGN INTERPOLATED VECTOR OF Z VALUES TO HGO
      do i = 2,new_num_zones
         spline_eval_x = prev_model%old_shell_mass(i)
         call splintd2(spline_x, spline_y, num_zones, spline_second_deriv, &
              spline_eval_x, spline_eval_y, spline_klo, spline_khi)
         z_new(i) = spline_eval_y
!
      end do
! TEST ON Z-CHANGE (ONLY FOR DECREASING Z) USING HIO AS DUMMY ARRAY
! SKIP IF HPMAX(4) IS ZEROED OUT
!

      if (point_spacing_max(4).le.1.0D-15) goto 103
      gradient_flag_count = 0
      do j = new_num_zones,2,-1
         if (z_new(j-1)-z_new(j).gt.point_spacing_max(4)) then
            gradient_flag_count = gradient_flag_count + 1
            gradient_flag_index(gradient_flag_count) = j
        endif
      end do
!
      if (gradient_flag_count.eq.0) goto 103
      working_num_zones = new_num_zones
!
      do i = 1,gradient_flag_count
         j = gradient_flag_index(i)
         delta_z_over_max = (z_new(j-1) - z_new(j))/point_spacing_max(4)
! NUMBER OF NEW POINTS NEEDED
         num_new_points = int(delta_z_over_max)
         point_insert_spacing = (prev_model%old_shell_mass(j)-prev_model%old_shell_mass(j-1))/ &
              dfloat(num_new_points+1)
         do k = working_num_zones+num_new_points,j+num_new_points,-1
            prev_model%old_shell_mass(k) = prev_model%old_shell_mass(k-num_new_points)
         end do
         do k = j + num_new_points -1, j -1, -1
            prev_model%old_shell_mass(k) = prev_model%old_shell_mass(k+1) - point_insert_spacing
         end do
         working_num_zones = working_num_zones + num_new_points
      end do
      new_num_zones = working_num_zones
!
 103  continue
! DELETE NEW POINTS THAT ARE TOO CLOSE TOGETHER.
! (NOTE HDO IS BEING USED AS A DUMMY ARRAY HERE).
      j = 1
      prev_model%old_density(j) = prev_model%old_shell_mass(j)
      do 810 k = 2,new_num_zones-1
!
       if (prev_model%old_shell_mass(k) - prev_model%old_density(j).gt.chi_grid_scale(1)) then
          j = j + 1
          prev_model%old_density(j) = prev_model%old_shell_mass(k)
       endif
  810 continue
      j = j + 1
      prev_model%old_density(j) = prev_model%old_shell_mass(new_num_zones)
      new_num_zones = j
      do 820 j = 2,new_num_zones
       prev_model%old_shell_mass(j) = prev_model%old_density(j)
  820 continue
!

!  NOW LOCATE OUTER EDGE OF CONVECTIVE CORE AND INNER EDGE OF CONVECTIVE
!  ENVELOPE IN THE NEW POINT DISTRIBUTION.
      if (convective_core_edge_zone.gt.1) then
       do 823 j = 2,new_num_zones
          if (prev_model%old_shell_mass(j).gt.log_mass(convective_core_edge_zone)) goto 824
  823    continue
  824    convective_core_edge_zone = j - 1
      else
       convective_core_edge_zone = 1
      endif
      if (convective_envelope_edge_zone.lt.num_zones) then
       do 825 j = new_num_zones-1,1,-1
          if (prev_model%old_shell_mass(j).lt.log_mass(convective_envelope_edge_zone)) &
               goto 826
  825    continue
  826    convective_envelope_edge_zone = j + 1
      else
       convective_envelope_edge_zone = new_num_zones
      endif
      if (convective_core_edge_zone.gt.1) then
       do 827 i = 1,convective_core_edge_zone
          convective_flag(i) = .true.
  827    continue
       radiative_zone_begin = convective_core_edge_zone + 1
      else
       radiative_zone_begin = 1
      endif
      if (convective_envelope_edge_zone.lt.new_num_zones) then
       do 828 i = convective_envelope_edge_zone,new_num_zones
          convective_flag(i) = .true.
  828    continue
       radiative_zone_end = convective_envelope_edge_zone - 1
      else
       radiative_zone_end = new_num_zones
      endif
      if (radiative_zone_end.lt.1.or.radiative_zone_begin.gt.new_num_zones) &
           goto 830
      do 829 j = radiative_zone_begin,radiative_zone_end
       convective_flag(j) = .false.
  829 continue
  830 continue


!  NOW USE AN OSCILLATORY SPLINE TO FIT THE OLD RUN OF PHYSICAL VARIABLES
!  AT THE NEW RUN OF MASS POINTS.
      old_point_count = num_zones
      new_point_count = new_num_zones
!  XVAL=RUN OF NEW MODEL MASS COORDINATES(HSS)
!  XTAB = RUN OF OLD MODEL MASS CO-ORDINATES(HS)
!  YTAB = RUN OF VARIABLE WHOSE VALUE AT THE NEW MASS POINTS
!         IS DESIRED(HP,HT,HL,HR,HD,ETC.)
!  NTAB = NUMBER OF DATA POINTS(M)
!  NTOT=NUMBER OF POINTS AT WHICH SPLINE IS TO BE EVALUTED(MNEW)
!  YVAL = OUTPUT RUN OF VARIABLE VALUES AT THE NEW RUN OF MASS POINTS.
!  FORM OF CALL IS CALL OSPLIN(XVAL,YVAL,XTAB,YTAB,NTAB,NTOT)
!  DO EACH COMPOSITION IN ORDER USING HPO AND HTO AS DUMMY ARRAYS.
! 7/91 ADD ENTROPY TERM INTERPOLATION.
      do 904 j = 1,num_zones
         prev_model%old_pressure(j) = temperature_entropy_term(j)
  904 continue
      call osplin(prev_model%old_shell_mass,prev_model%old_temperature,log_mass,prev_model%old_pressure, &
           old_point_count,new_point_count)
      do 905 j = 1,new_num_zones
         temperature_entropy_term(j) = prev_model%old_temperature(j)
  905 continue

!


      do 906 j = 1,num_zones
         prev_model%old_pressure(j) = pressure_entropy_term(j)
  906 continue
      call osplin(prev_model%old_shell_mass,prev_model%old_temperature,log_mass,prev_model%old_pressure, &
           old_point_count,new_point_count)
      do 907 j = 1,new_num_zones
         pressure_entropy_term(j) = prev_model%old_temperature(j)
  907 continue
      do 911 j = 1,num_zones
         prev_model%old_pressure(j) = luminosity_entropy_term(j)
  911 continue
      call osplin(prev_model%old_shell_mass,prev_model%old_temperature,log_mass,prev_model%old_pressure, &
           old_point_count,new_point_count)
      do 912 j = 1,new_num_zones
         luminosity_entropy_term(j) = prev_model%old_temperature(j)
  912 continue
      do 913 j = 1,num_zones
         prev_model%old_pressure(j) = radius_entropy_term(j)
  913 continue
      call osplin(prev_model%old_shell_mass,prev_model%old_temperature,log_mass,prev_model%old_pressure, &
           old_point_count,new_point_count)
      do 914 j = 1,new_num_zones
         radius_entropy_term(j) = prev_model%old_temperature(j)
  914 continue


      do 850 i = 1,num_species_tracked
       do 833 j = 1,num_zones
          prev_model%old_pressure(j) = composition(i,j)
  833    continue
         call osplin(prev_model%old_shell_mass,prev_model%old_temperature,log_mass,prev_model%old_pressure, &
              old_point_count,new_point_count)
       do 835 j = 1,new_num_zones
          composition(i,j) = prev_model%old_temperature(j)
  835    continue
!  HCOMPP IS THE ARRAY OF COMPOSITION AT THE BEGINNING OF THE TIMESTEP.
!  THIS IS NEEDED FOR COMPOSITION DIFFUSION IN ROTATING MODELS.
       do 840 j = 1,num_zones
          prev_model%old_pressure(j) = prev_model%old_composition(i,j)
  840    continue
         call osplin(prev_model%old_shell_mass,prev_model%old_temperature,log_mass,prev_model%old_pressure, &
              old_point_count,new_point_count)
       do 845 j = 1,new_num_zones
          prev_model%old_composition(i,j) = prev_model%old_temperature(j)
  845    continue
  850 continue


!  HCOMPM IS THE ARRAY OF CHANGES IN COMPOSITION DUE TO NUCLEAR BURNING.
!  THIS IS NEEDED FOR COMPOSITION DIFFUSION IN ROTATING MODELS.
      do 849 i = 1,7
       do 847 j = 1,num_zones
          prev_model%old_pressure(j) = reaction_rate_by_zone(reaction_rate_species_index(i),j)
  847    continue
         call osplin(prev_model%old_shell_mass,prev_model%old_temperature,log_mass,prev_model%old_pressure, &
              old_point_count,new_point_count)
       do 848 j = 1,new_num_zones
          reaction_rate_by_zone(reaction_rate_species_index(i),j) = &
               prev_model%old_temperature(j)
  848    continue
  849 continue
! MHP 05/02 IF THE SURFACE DEUTERIUM IS ABOVE
! THRESHOLD (1.0D-14) FIND THE NEW RUN OF
! DEUTERIUM BURNING RATES
      if (use_extended_composition .and. &
           composition(12,num_zones).ge.1.0D-14) then
         do j = 1,num_zones
            prev_model%old_pressure(j) = light_burn%deuterium_burning_rate_start(j)
         end do
         call osplin(prev_model%old_shell_mass,prev_model%old_temperature,log_mass,prev_model%old_pressure, &
              old_point_count,new_point_count)
         do j = 1,new_num_zones
            light_burn%deuterium_burning_rate_start(j) = prev_model%old_temperature(j)
         end do
      endif
! NOW FIND RUN OF P,R,L,T,AND RHO IN THAT ORDER FOR THE NEW POINTS.

      call osplin(prev_model%old_shell_mass,prev_model%old_pressure,log_mass,log_pressure, &
           old_point_count,new_point_count)
      call osplin(prev_model%old_shell_mass,prev_model%old_radius,log_mass,log_radius, &
           old_point_count,new_point_count)
      call osplin(prev_model%old_shell_mass,prev_model%old_luminosity,log_mass,log_luminosity, &
           old_point_count,new_point_count)
      call osplin(prev_model%old_shell_mass,prev_model%old_temperature,log_mass,log_temperature, &
           old_point_count,new_point_count)
      call osplin(prev_model%old_shell_mass,prev_model%old_density,log_mass,log_density, &
           old_point_count,new_point_count)

! FOR ROTATING MODELS FIND THE NEW RUN OF OMEGA,J/M,FP,FT,R0,AND ETA2.
      if (rotation_active) then
         call osplin(prev_model%old_shell_mass,old_omega,log_mass,omega, &
              old_point_count,new_point_count)
         call osplin(prev_model%old_shell_mass,old_specific_angular_momentum,log_mass, &
              specific_angular_momentum,old_point_count,new_point_count)
         call osplin(prev_model%old_shell_mass,fp_old,log_mass,fp, &
              old_point_count,new_point_count)
         call osplin(prev_model%old_shell_mass,ft_old,log_mass,ft, &
              old_point_count,new_point_count)
         call osplin(prev_model%old_shell_mass,old_mean_radius,log_mass,mean_radius, &
              old_point_count,new_point_count)
         call osplin(prev_model%old_shell_mass,old_eta_squared,log_mass,eta_squared, &
              old_point_count,new_point_count)
      endif
!

!     SPIT OUT POINT DISTRIBUTION DETAILS IF REQUESTED
      if (ldebug .and.  npoint.lt.9999) then
      if (mod(model_number,npoint).eq.0) then
         min_common_count = min0(num_zones,new_num_zones)
!
         write(idebug,910)
  910    format('1',20X,'OLD POINTS',54X,'NEW POINTS'/2(3X,'N',5X,'S', &
     8X,'P',7X,'T',7X,'R',8X,'L',7X,'X',4X,'Z',3X,'O16',1X) )
         write(idebug,920) (i,log_mass(i),log_pressure(i), &
              log_temperature(i),log_radius(i),log_luminosity(i), &
              x_new(i),z_new(i),composition(9,i),i,prev_model%old_shell_mass(i), &
              prev_model%old_pressure(i),prev_model%old_temperature(i),prev_model%old_radius(i), &
              prev_model%old_luminosity(i),composition(1,i),composition(3,i), &
              composition(9,i), i = 1,min_common_count)
  920    format( 2(1X,I3,F11.7,F8.4,F8.5,F8.4,1PE9.2,0PF6.3,2F5.3) )
         if (num_zones.gt.min_common_count) then
            min_common_count = min_common_count + 1
            write(idebug,930) (i,log_mass(i),log_pressure(i), &
                 log_temperature(i),log_radius(i),log_luminosity(i), &
                 x_new(i),z_new(i),composition(9,i),i=min_common_count, &
                 num_zones)
  930       format( 1X,I3,F11.7,F8.4,F8.5,F8.4,1PE9.2,0PF6.3,2F5.3)
         else if (new_num_zones.gt.min_common_count) then
            min_common_count = min_common_count + 1
            write(idebug,940)(i,prev_model%old_shell_mass(i),prev_model%old_pressure(i), &
                 prev_model%old_temperature(i),prev_model%old_radius(i),prev_model%old_luminosity(i), &
                 composition(1,i),composition(3,i),composition(9,i), &
                 i=min_common_count,new_num_zones)
  940       format(65X,I3,F11.7,F8.4,F8.5,F8.4,1PE9.2,0PF6.3,2F5.3)
         endif
      endif
      endif


! TRANSFER NEW POINTS.
      do 1000 j = 1,new_num_zones
       log_mass(j) = prev_model%old_shell_mass(j)
       log_pressure(j) = prev_model%old_pressure(j)
       log_temperature(j) = prev_model%old_temperature(j)
       log_radius(j) = prev_model%old_radius(j)
       log_luminosity(j) = prev_model%old_luminosity(j)
       log_density(j) = prev_model%old_density(j)
 1000 continue
      if (rotation_active) then
       do 1005 j = 1, new_num_zones
          specific_angular_momentum(j) = old_specific_angular_momentum(j)
          omega(j) = old_omega(j)
          fp(j) = fp_old(j)
          ft(j) = ft_old(j)
          eta_squared(j) = old_eta_squared(j)
          mean_radius(j) = old_mean_radius(j)
 1005    continue
      endif
! MHP 6/00 INTERPOLATED IN ENERGY GENERATION AT START OF TIMESTEP
      if (rotation_active .or. (use_extended_composition .and. &
           envelope_overshoot_active)) then
         call osplin(prev_model%old_shell_mass,old_esum,log_mass,shell_diag%sesum, &
              old_point_count,new_point_count)
         do zone_index = 1,num_zones
            spline_y(zone_index) = shell_diag%sesum(zone_index)+shell_diag%seg(6,zone_index)+ &
                 shell_diag%seg(7,zone_index)
         end do
         call osplin(prev_model%old_shell_mass,old_eps,log_mass,spline_y, &
              old_point_count,new_point_count)
      endif
      write(short_file_unit,1020) num_zones,new_num_zones
 1020 format(' POINTS  OLD',I5,'   NEW',I5)
      num_zones = new_num_zones
! SET UP WEIGHTS AND MASSES
      mass_curr = dexp(clndp*log_mass(1))
      mass_prev = - mass_curr
      do 1030 i = 2,num_zones
       mass_two_back = mass_prev
       mass_prev = mass_curr
       mass_curr = dexp(clndp*log_mass(i))
       enclosed_mass(i-1) = mass_prev
       shell_mass(i-1) = 0.5D0*(mass_curr-mass_two_back)
 1030 continue
      enclosed_mass(num_zones) = mass_curr
      shell_mass(num_zones) = dexp(ln10*log_total_mass) - &
           0.5D0*(mass_prev+mass_curr)
      if (rotation_active) then
!  FIRST GUESS AT MOMENT OF INERTIA(HI)
       do 1070 i=1,num_zones
          moment_of_inertia(i) = cc23*shell_mass(i)* &
               dexp(ln10*2.0D0*log_radius(i))
 1070    continue
!   CALCULATE OVERSHOOT
       call ovrot(composition,log_density,log_pressure,log_radius, &
            log_mass,log_temperature,convective_flag,num_zones, &
            am_transport_convective_flag,radiative_zone_bounds, &
            convective_zone_bounds,num_radiative_zones,num_convective_zones)
! JNT 2025/09/03 duplicating 2015/04/06 recompute moment of interia
! before recomputing the rotation I am less confident that this is
! necessary since WALPCZ does run in this version but I don't think
! it can hurt.
       call momi(eta_squared,log_radius,log_mass,shell_mass,1,num_zones, &
            omega,mean_radius,moment_of_inertia,qiw)
! END JNT

!   FIND THE ANGULAR VELOCITY OMEGA THAT CORRESPONDS TO THE GIVEN
!   SPECIFIC ANGULAR MOMENTUM HJM.
       call getrot(log_density,specific_angular_momentum,log_radius, &
            log_mass,shell_mass,am_transport_convective_flag,num_zones, &
            eta_squared,moment_of_inertia,omega,qiw,mean_radius)
!  CALCULATE FP,FT,R0 AND ETA2 GIVEN OMEGA
       call fpft(log_density,log_radius,log_mass,num_zones,omega, &
            eta_squared,fp,ft,hg,mean_radius)
!  FIND CORRECT MOMENT OF INERTIA(HI)
!        CALL MOMI(ETA2,HD,HR,HS,HS2,1,M,OMEGA,R0,HI,QIW,M)  ! KC 2025-05-31
       call momi(eta_squared,log_radius,log_mass,shell_mass,1,num_zones, &
            omega,mean_radius,moment_of_inertia,qiw)
!  FIND NEW TOTAL ANGULAR MOMENTUM
       sum_angular_momentum = 0.0D0
       sum_rotational_ke = 0.0D0
       do 1110 i = 1,num_zones
          angular_momentum_shell = specific_angular_momentum(i)*shell_mass(i)
          rotational_kinetic_energy(i) = 0.5D0*omega(i)*angular_momentum_shell
          sum_angular_momentum = sum_angular_momentum + angular_momentum_shell
          sum_rotational_ke = sum_rotational_ke + rotational_kinetic_energy(i)
 1110    continue
       write(short_file_unit,1120)total_angular_momentum, &
            sum_angular_momentum,total_rotational_ke,sum_rotational_ke
 1120    format(1X,'TOTAL J OF STAR - PREVIOUS ',1PE21.13,' NEW ', &
     1PE21.13/' TOTAL ROTATIONAL K.E. OF STAR-PREVIOUS ',1PE21.13, &
     ' NEW ',1PE21.13)
       total_angular_momentum = sum_angular_momentum
       total_rotational_ke = sum_rotational_ke
!  STORE THE OLD MODEL STRUCTURE FOR USE IN DIFFUSION.
       do 1130 i = 1,num_zones
          old_omega(i) = omega(i)
          old_hg(i) = hg(i)
          old_moment_of_inertia(i) = moment_of_inertia(i)
          old_eta_squared(i) = eta_squared(i)
          old_mean_radius(i) = mean_radius(i)
          prev_model%old_convective_flag(i) = convective_flag(i)
          prev_model%old_cz_flag(i) = am_transport_convective_flag(i)
! MHP 10/91 J/M STORED IN HJX FOR I/O USE.
            old_specific_angular_momentum(i) = specific_angular_momentum(i)
 1130    continue
! MHP 9/91 CHANGE : T GRADIENTS STORED IF LEXCOM=T AND LOVSTE=T; OR FOR
! ROTATION; THIS IS NEEDED SO THAT THE BASE OF THE OVERSHOOT REGION FOR
! PRE-MS MODELS CAN BE ACCURATELY LOCATED.
      endif
      if (rotation_active .or. (use_extended_composition .and. &
           envelope_overshoot_active)) then
! END OF 9/91 CHANGE
!   FIND THE NEW RUN OF PHYSICAL VARIABLES AT THE NEW SET OF POINTS;
!   THIS IS NEEDED EVEN IN THE ABSENCE OF DIFFUSION TO ACCURATELY LOCATE
!   THE EDGES OF CONVECTION ZONES.
!        CALL PHYSIC(FP,FT,HCOMP,HD,HG,HL,HP,HR,HS,HT,LC,LCZ,M,TEFFL)  ! KC 2025-05-31
       call physic(fp,ft,composition,log_density,hg,log_luminosity, &
            log_pressure,log_radius,log_mass,log_temperature, &
            convective_flag,num_zones,log_teff)
!   FOR DIFFUSION STORE THE AUXILLARY QUANTITIES NEEDED TO CALCULATE
!   VELOCITIES AT THE START OF THE TIMESTEP WITH THE NEW POINT DISTRIBUTION
!   SO THAT A SERIES OF SMALL DIFFUSION TIMESTEPS CAN BE TAKEN WITHIN
!   ONE LARGE EVOLUTIONARY TIMESTEP.
         do 1040 zone_index = 1,num_zones
            old_del_radiative_mix(zone_index) = shell_diag%del_grad(1,zone_index)
            old_delm(zone_index) = shell_diag%del_grad(2,zone_index)
            old_del_adiabatic_mix(zone_index) = shell_diag%del_grad(3,zone_index)
            old_amu(zone_index) = mean_molecular_weight(zone_index)
            old_om(zone_index) = shell_diag%so(zone_index)
            old_cp(zone_index) = cp(zone_index)
            old_qdt(zone_index) = qdt(zone_index)
            old_vel(zone_index) = shell_diag%svel(zone_index)
            old_visc(zone_index) = visc(zone_index)
            old_thdif(zone_index) = thdif(zone_index)
! MHP 06/02
            del_grad_diff_interface(zone_index) = &
                 old_del_adiabatic_mix(zone_index) - old_delm(zone_index)
! MHP 6/00 CALCULATED EARLIER
!            ESUMO(IM) = SESUM(IM)
            max_domega_dr_old(zone_index) = max_domega_dr(zone_index)
 1040    continue
! MHP 06/02 ADDED TERM FOR THE TIME EVOLUTION
! OF THE ANGULAR VELOCITY DISTRIBUTION
         do i = 2,num_zones
            delta_radius = exp(ln10*log_radius(i))-exp(ln10*log_radius(i-1))
            delta_omega = omega(i) - omega(i-1)
            omega_mid = 0.5D0*(omega(i)+omega(i-1))
            log_factor = 2.0D0*(log_radius(i)+log_radius(i-1))-0.5D0* &
     (log_mass(i)+log_mass(i-1))-cgl
            tho(i) = exp(ln10*log_factor)*omega_mid*delta_omega/delta_radius
            qwrst(i) = delta_omega/delta_radius
         end do
      endif
!  CALCULATE NEW SURFACE OPACITY TABLE IF NEEDED.
      if (dabs(xnew-composition(1,num_zones)).gt.1.0D-8) then
               xnew = composition(1,num_zones)
               znew = composition(3,num_zones)
               call surfopac(xnew)

      end if

      return
end subroutine hpoint
