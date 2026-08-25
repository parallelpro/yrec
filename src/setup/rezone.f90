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
subroutine rezone(envelope_store_index, point_reset_flag, &
     h_shell_zone_begin, h_shell_active, total_angular_momentum, &
     total_rotational_ke, ierr)
      use star_info_lib, only: star, i_eps_grav, i_eps_neu, i_grad_actual, i_grad_ad, i_grad_rad, i_h1, i_h2, i_metals, i_o16, json
      use kap_lib
      use luout_lib
      use phys_const_lib
      use numerics_lib

      implicit none

      integer, intent(in) :: envelope_store_index
      logical, intent(out) :: point_reset_flag
      integer, intent(in) :: h_shell_zone_begin
      logical, intent(in) :: h_shell_active
      double precision, intent(inout) :: total_angular_momentum, &
           total_rotational_ke

! MHP 10/02 added MRZONE,MXZONE to dimension statements
      double precision :: log10_omega(json), point_spacing_max(4)
      integer :: flag_point(100)
      logical :: am_transport_convective_flag(json)
      double precision :: ft_old(json), fp_old(json)
      integer :: radiative_zone_bounds(13,2), convective_zone_bounds(12,2)

























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
      ! 2026 (ROADMAP.md stage 3): library errors return here via ierr;
      ! this driver-side call site preserves the historical stop.
      integer :: jerr

      integer, intent(out) :: ierr

      ierr = 0

      point_reset_flag = .false.
!  IEND IS THE NUMBER OF SPECIES THE PROGRAM IS KEEPING TRACK OF
      num_species_tracked = 11
      if (star%job%use_extended_composition) num_species_tracked = 15
      call check_envelope_temperature_range
      if (ierr /= 0) return
      call flag_fixed_points
      if (ierr /= 0) return
      call assign_new_points
      if (ierr /= 0) return
      call locate_new_cz_edges
      if (ierr /= 0) return
      call interpolate_onto_new_grid
      if (ierr /= 0) return

      return

contains

! ---------------------------------------------------------------
! Verify the outermost Henyey point sits inside the allowed
! envelope temperature window (tenv0/tenv), deleting or flagging
! points as needed; ierr if the whole model is cooler than tenv0.
subroutine check_envelope_temperature_range
! CHECK IF TEMPERATURE OF OUTERMOST HENYEY POINT > MIMIMUM ENVELOPE T
      if (star%logT(star%nz).lt.star%ctrl%tenv0) then
       do i = star%nz-1,1,-1
          if (star%logT(i).gt.star%ctrl%tenv0) then
             write(short_file_unit,10) star%nz,i
   10     format(' OUTER POINTS DELETED OLD M =',I5,'  NEW M =',I5)
             star%nz = i
             star%senv = star%log_mass(star%nz) - star%log_total_mass
             point_reset_flag = .true.
             exit
          endif
       end do
       if (i < (1)) then
!  ENTIRE MODEL HAS T<TENV0 - UNLIKELY - BUT STOP IF TRUE
       write(short_file_unit,30)
       write(iowr,30)
   30    format(1X,39('>'),40('<')/1X,'ERROR IN HPOINT'/ &
     1X,'MAX. STAR T LESS THAN MINIMUM ENVELOPE T.RUN STOPPED')
       ! 2026 (phase five, step B): stop converted to ierr; run_yrec
       ! returns the error and the CLI wrapper (main) stops.
       ierr = 1
       return
       end if
!  CHECK IF OUTER POINT T < MAXIMUM ENVELOPE T
      else if (star%logT(star%nz).gt.star%ctrl%tenv1.and. &
           envelope_store_index.ne.0) then
       star%nz = star%nz + 1
       star%senv = star%stored_envelope_state(4)
       star%log_mass(star%nz) = star%log_total_mass + star%senv
       star%logP(star%nz) = star%logP(star%nz-1) + &
            (star%stored_envelope_state(1) - star%fit_point_pressure(envelope_store_index))
       star%logT(star%nz) = star%logT(star%nz-1) + &
            (star%stored_envelope_state(1) - &
            star%fit_point_pressure(envelope_store_index))*0.250D0
       star%logR(star%nz) = star%logR(star%nz-1) + &
            (star%stored_envelope_state(3) - star%fit_point_radius(envelope_store_index))
       star%luminosity_lsun(star%nz) = star%luminosity_lsun(star%nz-1)
       do i = 1,num_species_tracked
          star%xa(i,star%nz) = star%xa(i,star%nz-1)
       end do
       star%logRho(star%nz) = star%logP(star%nz) - &
            star%logT(star%nz) - 8.0D0
       j = star%nz - 1
       write(short_file_unit,60) j,star%log_mass(j),star%logP(j), &
            star%logT(j),star%logR(j),star%nz,star%log_mass(star%nz), &
            star%logP(star%nz),star%logT(star%nz), &
            star%logR(star%nz)
   60    format(' OUTER POINT ADDED',I5,F15.10,'  PTR',3F10.6)
       point_reset_flag = .true.
      endif


end subroutine check_envelope_temperature_range

! ---------------------------------------------------------------
! Set up the flagged points the rezoner must not smooth across:
! convection-zone edges, the H-burning shell edge, X/Z/omega
! gradient jumps, the overshoot base and the finely-zoned region
! around the surface CZ; sort and de-duplicate them.
subroutine flag_fixed_points
! SET UP FLAGGED POINTS - PROGRAM WILL NOT REZONE ACROSS FLAGGED POINTS
      flag_count = 1
! FLAG EDGES OF CENTRAL AND SURFACE CONVECTION ZONES
      if (star%core_cz_top_index.gt.1) then
       flag_point(flag_count) = star%core_cz_top_index
       flag_count = flag_count + 1
      endif
      if (star%envelope_cz_bottom_index.lt.star%nz .and. &
           star%envelope_cz_bottom_index.gt.1) then
       flag_point(flag_count) = star%envelope_cz_bottom_index
       flag_count = flag_count + 1
      endif
! FLAG EDGE OF H SHELL
      if (h_shell_active.and.h_shell_zone_begin.gt.1) then
       flag_point(flag_count) = h_shell_zone_begin - 1
       flag_count = flag_count + 1
      endif
      do i = 2,star%nz
! TEST FOR FLAGGING DUE TO X GRADIENT
       if (dabs(star%xa(i_h1,i)-star%xa(i_h1,i-1)).gt.star%ctrl%chi_grid_scale(3)) then
          flag_point(flag_count) = i
          flag_count = flag_count + 1
! TEST FOR FLAGGING DUE TO Z GRADIENT
       else if (dabs(star%xa(i_metals,i)-star%xa(i_metals,i-1)).gt. &
            star%ctrl%chi_grid_scale(4)) then
          flag_point(flag_count) = i
          flag_count = flag_count + 1
! TEST FOR FLAGGING DUE TO GRADIENT IN LOG OMEGA.
       else if (star%job%rotation_active) then
          log_omega_top = dlog10(star%omega(i))
          log_omega_bot = dlog10(star%omega(i-1))
          if (dabs(log_omega_top-log_omega_bot).gt.star%ctrl%chi_grid_scale(12)) then
             flag_point(flag_count) = i
             flag_count = flag_count + 1
          endif
       endif
       if (flag_count.ge.100) then
          write(short_file_unit,110)
  110       format(1X,'MORE THAN 100 FLAG POINTS-FIRST 100 RETAINED')
          exit
       endif
      end do
      if (i > (star%nz)) then
!  PMAX1 = MAX DEL LOG P BELOW SURFACE C.Z. AND BELOW FINELY ZONED
!  REGION AROUND IT.
!  PMAX2 = MAX DEL LOG P BETWEEN LOWER EDGE OF FINELY ZONED REGION
!  AROUND SURFACE C.Z. AND BASE OF OVERSHOOT REGION.
!  PMAX3 = SAME FOR OVERSHOOT REGION.
!  PMAX4 = MAX DEL LOG P ABOVE BASE OF SURFACE C.Z. IN FINELY ZONED
!  REGION AROUND IT.
!  PMAX5 = SAME FOR THE OUTER POINTS IN THE STAR.
      pmax1 = star%ctrl%chi_grid_scale(11)
      pmax4 = star%ctrl%chi_grid_scale(10)
      pmax5 = star%ctrl%chi_grid_scale(8)
      if (.not.star%convective_flag(star%nz)) then
       overshoot_base_zone = star%nz
       fine_zone_base = star%nz
      else if (star%envelope_cz_bottom_index.eq.1) then
       overshoot_base_zone = 1
       fine_zone_base = 1
      else
! LOCATE BASE OF OVERSHOOT REGION IF APPLICABLE.
       if (.not.star%job%envelope_overshoot_active) then
          i = star%envelope_cz_bottom_index
          overshoot_base_zone = star%envelope_cz_bottom_index
       else
          do overshoot_base_zone = star%envelope_cz_bottom_index-1,1,-1
             if (star%logP(overshoot_base_zone)- &
                  star%logP(star%envelope_cz_bottom_index).gt.star%ctrl%alphae) exit
          end do
            overshoot_base_zone = overshoot_base_zone + 1
          delta_log_pressure = star%logP(overshoot_base_zone)- &
               star%logP(star%envelope_cz_bottom_index)
          if (delta_log_pressure.gt.0.0D0) then
             overshoot_point_count = int(delta_log_pressure/star%ctrl%chi_grid_scale(10))
             if (mod(delta_log_pressure,star%ctrl%chi_grid_scale(10)).ne.0D0) &
                  overshoot_point_count = overshoot_point_count+1
             pmax3 = delta_log_pressure/dfloat(overshoot_point_count)
          else
             pmax3 = star%ctrl%chi_grid_scale(10)
          endif
            if (overshoot_base_zone.gt.1) then
               flag_point(flag_count) = overshoot_base_zone
               flag_count = flag_count + 1
            endif
       endif
       if (star%ctrl%chi_grid_scale(7).eq.0.0D0) then
          fine_zone_base = overshoot_base_zone
       else
! NOW LOCATE BASE OF FINELY ZONED REGION.
          if (overshoot_base_zone.eq.1) then
             fine_zone_base = 1
          else
          do fine_zone_base = overshoot_base_zone-1,1,-1
             if (star%logP(fine_zone_base) - &
                  star%logP(star%envelope_cz_bottom_index).gt.star%ctrl%chi_grid_scale(7)) &
                  exit
          end do
            fine_zone_base = fine_zone_base + 1
          if (fine_zone_base.ne.overshoot_base_zone) then
          delta_log_pressure = star%logP(fine_zone_base) - &
               star%logP(overshoot_base_zone)
          overshoot_point_count = int(delta_log_pressure/star%ctrl%chi_grid_scale(10))
          if (mod(delta_log_pressure,star%ctrl%chi_grid_scale(10)).ne.0D0) &
               overshoot_point_count = overshoot_point_count+1
          pmax2 = delta_log_pressure/dfloat(overshoot_point_count)
          if (pmax2.eq.0.0D0) pmax2 = star%ctrl%chi_grid_scale(10)
            if (fine_zone_base.gt.1) then
               flag_point(flag_count) = fine_zone_base
               flag_count = flag_count + 1
            endif
          end if
          end if
       endif
      endif
      end if
      flag_point(flag_count) = star%nz
! ARRANGE THE FLAG POINTS IN ASCENDING ORDER
      if (flag_count.ne.1) then
      do
         continue
      sort_done = .true.
      do i = 1,flag_count-1
       if (flag_point(i+1).lt.flag_point(i)) then
          temp_swap = flag_point(i)
          flag_point(i) = flag_point(i+1)
          flag_point(i+1) = temp_swap
          sort_done = .false.
       endif
      end do
      if (sort_done) exit
      end do
! ENSURE THAT POINTS ARENT FLAGGED MORE THAN ONCE.
      i = 2
      do
         continue
      if (flag_point(i).eq.flag_point(i-1)) then
       if (i.lt.flag_count) then
          do j = i,flag_count-1
             flag_point(j) = flag_point(j+1)
          end do
       endif
       flag_count = flag_count - 1
      endif
      i = i + 1
      if (i.gt.flag_count) exit
      end do
      end if
      write(short_file_unit,185) (flag_point(j),j=1,flag_count)
  185 format(1X,'FLAG-POINTS',20I4)
end subroutine flag_fixed_points

! ---------------------------------------------------------------
! Refloat the point distribution between flagged points: accumulate
! the normalized chi = dM + dL + dP spacing measure per region
! (pmax1..5), add points where X or Z gradients demand them, and
! delete new points that land too close together. ierr if the
! requested point count exceeds json.
subroutine assign_new_points
! BEGIN REFLOATING OF POINTS
      if (star%job%rotation_active) then
       do i = 1,star%nz
          if (star%omega(i).gt.0.0D0) then
             log10_omega(i) = dlog10(star%omega(i))
          else
             log10_omega(i) = 0.0D0
          endif
       end do
      endif
      star%old_shell_mass(1) = star%log_mass(1)
      star%logP_start(1) = star%logP(1)
      star%luminosity_lsun_start(1) = star%luminosity_lsun(1)
      x_new(1) = star%xa(i_h1,1)
      z_new(1) = star%xa(i_metals,1)
      luminosity_max = star%luminosity_lsun(star%nz)
!       JVS 04/14 added Teff to saved variables
        star%log_Teff_start = star%log_Teff
!  JVS 05/25 Added model number to list of saved values
      star%nz_start = star%nz
      do i = star%nz-1,1,-1
         if (star%luminosity_lsun(i).gt.luminosity_max) then
            luminosity_max = star%luminosity_lsun(i)
         endif
      end do
      point_spacing_max(1) = star%ctrl%chi_grid_scale(8)
      point_spacing_max(2) = star%ctrl%chi_grid_scale(9)*luminosity_max
      point_spacing_max(3) = star%ctrl%chi_grid_scale(5)
      point_spacing_max(4) = star%ctrl%chi_grid_scale(6)
!      KFACT = 0
!  200 CONTINUE
      point_insert_flag = 1
      j = 2
      chi_start_index = 2
! CHI IS THE NORMALIZED VECTOR OF DIFFERENCES IN M,L,P:
! CHI = HS/DELTA M + HL/DELTA L - HP/DELTA P
      mass_scale = star%ctrl%chi_grid_scale(2)
      luminosity_scale = point_spacing_max(2)
      chi(1) = 1.0D0
      do j = 2, star%nz
         pressure_test = star%logP(j) - star%logP(star%envelope_cz_bottom_index)
         if (abs(pressure_test).lt.star%ctrl%chi_grid_scale(7)) then
! FINELY ZONED REGION
            dp_scale = star%ctrl%chi_grid_scale(10)
         else if (pressure_test.gt.star%ctrl%chi_grid_scale(7)) then
! BELOW SURFACE CZ
            dp_scale = star%ctrl%chi_grid_scale(11)
         else
! IN SURFACE CZ
            dp_scale = star%ctrl%chi_grid_scale(8)
         endif
         if (star%luminosity_lsun(j).gt.star%luminosity_lsun(j-1)) then
            dchi = (star%log_mass(j)-star%log_mass(j-1))/mass_scale + &
                 (star%luminosity_lsun(j)-star%luminosity_lsun(j-1))/luminosity_scale - &
                 (star%logP(j)-star%logP(j-1))/dp_scale
         else
            dchi = (star%log_mass(j)-star%log_mass(j-1))/mass_scale - &
                 (star%logP(j)-star%logP(j-1))/dp_scale
         endif
         chi(j) = chi(j-1)+dchi
      end do
      do j = 1,star%nz
         spline_x(j) = chi(j)
         spline_y(j) = star%log_mass(j)
      end do
! GET SPLINE COEFFICIENTS
      call splinc(spline_x,spline_y,spline_second_deriv,star%nz)
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
            ! 2026 (phase five, step B): stop converted to ierr; run_yrec
            ! returns the error and the CLI wrapper (main) stops.
            ierr = 1
            return
         endif
         dchi = delta_chi/dfloat(segment_point_count)
! ASSIGN NEW POINTS
         do j = point_prev_index+1,point_prev_index+segment_point_count
            spline_eval_x = chi_prev + dchi
            call splintd2(spline_x, spline_y, star%nz, &
                 spline_second_deriv, spline_eval_x, spline_eval_y, &
                 spline_klo, spline_khi)
            star%old_shell_mass(j) = spline_eval_y
            chi_prev = spline_eval_x
!

         end do
         new_num_zones = new_num_zones + segment_point_count
         point_prev_index = new_num_zones
!
      end do








! TEST FOR ASSIGNING POINTS BASED ON THE GRADIENT IN X.
      do j = 1,star%nz
         spline_x(j) = star%log_mass(j)
         spline_y(j) = star%xa(i_h1,j)
      end do
! GET SPLINE COEFFICIENTS
      call splinc(spline_x,spline_y,spline_second_deriv,star%nz)
! ASSIGN INTERPOLATED VECTOR OF X VALUES TO HIO
      do i = 2,new_num_zones
         spline_eval_x = star%old_shell_mass(i)
         call splintd2(spline_x, spline_y, star%nz, spline_second_deriv, &
              spline_eval_x, spline_eval_y, spline_klo, spline_khi)
         x_new(i) = spline_eval_y
      end do




! SKIP IF HPMAX(3) IS ZEROED OUT
      if (point_spacing_max(3).gt.1.0D-15) then
! TEST ON X-CHANGE (ONLY FOR INCREASING X) USING HIO AS DUMMY ARRAY
      gradient_flag_count = 0
      do j = new_num_zones,2,-1
         if (x_new(j)-x_new(j-1).gt.point_spacing_max(3)) then
            gradient_flag_count = gradient_flag_count + 1
            gradient_flag_index(gradient_flag_count) = j
         endif
      end do
      if (gradient_flag_count.ne.0) then
      working_num_zones = new_num_zones
      do i = 1,gradient_flag_count
         j = gradient_flag_index(i)
         delta_x_over_max = (x_new(j) - x_new(j-1))/point_spacing_max(3)
! NUMBER OF NEW POINTS NEEDED
         num_new_points = int(delta_x_over_max)
!
         point_insert_spacing = (star%old_shell_mass(j)-star%old_shell_mass(j-1))/ &
              dfloat(num_new_points+1)
         do k = working_num_zones+num_new_points,j+num_new_points,-1
            star%old_shell_mass(k) = star%old_shell_mass(k-num_new_points)
!
         end do
         do k = j + num_new_points -1, j -1, -1
            star%old_shell_mass(k) = star%old_shell_mass(k+1) - point_insert_spacing
         end do
         working_num_zones = working_num_zones + num_new_points
      end do
      new_num_zones = working_num_zones
!
      end if
      end if
! TEST FOR ASSIGNING POINTS BASED ON THE GRADIENT IN Z.
      do j = 1,star%nz
         spline_x(j) = star%log_mass(j)
         spline_y(j) = star%xa(i_metals,j)
      end do
! GET SPLINE COEFFICIENTS
      call splinc(spline_x,spline_y,spline_second_deriv,star%nz)
! ASSIGN INTERPOLATED VECTOR OF Z VALUES TO HGO
      do i = 2,new_num_zones
         spline_eval_x = star%old_shell_mass(i)
         call splintd2(spline_x, spline_y, star%nz, spline_second_deriv, &
              spline_eval_x, spline_eval_y, spline_klo, spline_khi)
         z_new(i) = spline_eval_y
!
      end do
! TEST ON Z-CHANGE (ONLY FOR DECREASING Z) USING HIO AS DUMMY ARRAY
! SKIP IF HPMAX(4) IS ZEROED OUT
!

      if (point_spacing_max(4).gt.1.0D-15) then
      gradient_flag_count = 0
      do j = new_num_zones,2,-1
         if (z_new(j-1)-z_new(j).gt.point_spacing_max(4)) then
            gradient_flag_count = gradient_flag_count + 1
            gradient_flag_index(gradient_flag_count) = j
        endif
      end do
!
      if (gradient_flag_count.ne.0) then
      working_num_zones = new_num_zones
!
      do i = 1,gradient_flag_count
         j = gradient_flag_index(i)
         delta_z_over_max = (z_new(j-1) - z_new(j))/point_spacing_max(4)
! NUMBER OF NEW POINTS NEEDED
         num_new_points = int(delta_z_over_max)
         point_insert_spacing = (star%old_shell_mass(j)-star%old_shell_mass(j-1))/ &
              dfloat(num_new_points+1)
         do k = working_num_zones+num_new_points,j+num_new_points,-1
            star%old_shell_mass(k) = star%old_shell_mass(k-num_new_points)
         end do
         do k = j + num_new_points -1, j -1, -1
            star%old_shell_mass(k) = star%old_shell_mass(k+1) - point_insert_spacing
         end do
         working_num_zones = working_num_zones + num_new_points
      end do
      new_num_zones = working_num_zones
!
      end if
      end if
! DELETE NEW POINTS THAT ARE TOO CLOSE TOGETHER.
! (NOTE HDO IS BEING USED AS A DUMMY ARRAY HERE).
      j = 1
      star%logRho_start(j) = star%old_shell_mass(j)
      do k = 2,new_num_zones-1
!
       if (star%old_shell_mass(k) - star%logRho_start(j).gt.star%ctrl%chi_grid_scale(1)) then
          j = j + 1
          star%logRho_start(j) = star%old_shell_mass(k)
       endif
      end do
      j = j + 1
      star%logRho_start(j) = star%old_shell_mass(new_num_zones)
      new_num_zones = j
      do j = 2,new_num_zones
       star%old_shell_mass(j) = star%logRho_start(j)
      end do
!

end subroutine assign_new_points

! ---------------------------------------------------------------
! Locate the convective-core outer edge and convective-envelope
! inner edge in the new point distribution.
subroutine locate_new_cz_edges
!  NOW LOCATE OUTER EDGE OF CONVECTIVE CORE AND INNER EDGE OF CONVECTIVE
!  ENVELOPE IN THE NEW POINT DISTRIBUTION.
      if (star%core_cz_top_index.gt.1) then
       do j = 2,new_num_zones
          if (star%old_shell_mass(j).gt.star%log_mass(star%core_cz_top_index)) exit
       end do
         star%core_cz_top_index = j - 1
      else
       star%core_cz_top_index = 1
      endif
      if (star%envelope_cz_bottom_index.lt.star%nz) then
       do j = new_num_zones-1,1,-1
          if (star%old_shell_mass(j).lt.star%log_mass(star%envelope_cz_bottom_index)) &
               exit
       end do
         star%envelope_cz_bottom_index = j + 1
      else
       star%envelope_cz_bottom_index = new_num_zones
      endif
      if (star%core_cz_top_index.gt.1) then
       do i = 1,star%core_cz_top_index
          star%convective_flag(i) = .true.
       end do
       radiative_zone_begin = star%core_cz_top_index + 1
      else
       radiative_zone_begin = 1
      endif
      if (star%envelope_cz_bottom_index.lt.new_num_zones) then
       do i = star%envelope_cz_bottom_index,new_num_zones
          star%convective_flag(i) = .true.
       end do
       radiative_zone_end = star%envelope_cz_bottom_index - 1
      else
       radiative_zone_end = new_num_zones
      endif
      if (radiative_zone_end.ge.1 .and. radiative_zone_begin.le.new_num_zones) then
      do j = radiative_zone_begin,radiative_zone_end
       star%convective_flag(j) = .false.
      end do
      end if


end subroutine locate_new_cz_edges

! ---------------------------------------------------------------
! Osculatory-spline fit of every physical variable and species onto
! the new mass grid, refresh the rotation profile and diffusion
! auxiliaries, re-run physic, and update the surface opacity
! tables if the surface composition moved.
subroutine interpolate_onto_new_grid
!  NOW USE AN OSCILLATORY SPLINE TO FIT THE OLD RUN OF PHYSICAL VARIABLES
!  AT THE NEW RUN OF MASS POINTS.
      old_point_count = star%nz
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
      do j = 1,star%nz
         star%logP_start(j) = star%temperature_entropy_term(j)
      end do
      call osplin(star%old_shell_mass,star%logT_start,star%log_mass,star%logP_start, &
           old_point_count,new_point_count)
      do j = 1,new_num_zones
         star%temperature_entropy_term(j) = star%logT_start(j)
      end do

!


      do j = 1,star%nz
         star%logP_start(j) = star%pressure_entropy_term(j)
      end do
      call osplin(star%old_shell_mass,star%logT_start,star%log_mass,star%logP_start, &
           old_point_count,new_point_count)
      do j = 1,new_num_zones
         star%pressure_entropy_term(j) = star%logT_start(j)
      end do
      do j = 1,star%nz
         star%logP_start(j) = star%luminosity_entropy_term(j)
      end do
      call osplin(star%old_shell_mass,star%logT_start,star%log_mass,star%logP_start, &
           old_point_count,new_point_count)
      do j = 1,new_num_zones
         star%luminosity_entropy_term(j) = star%logT_start(j)
      end do
      do j = 1,star%nz
         star%logP_start(j) = star%radius_entropy_term(j)
      end do
      call osplin(star%old_shell_mass,star%logT_start,star%log_mass,star%logP_start, &
           old_point_count,new_point_count)
      do j = 1,new_num_zones
         star%radius_entropy_term(j) = star%logT_start(j)
      end do


      do i = 1,num_species_tracked
       do j = 1,star%nz
          star%logP_start(j) = star%xa(i,j)
       end do
         call osplin(star%old_shell_mass,star%logT_start,star%log_mass,star%logP_start, &
              old_point_count,new_point_count)
       do j = 1,new_num_zones
          star%xa(i,j) = star%logT_start(j)
       end do
!  HCOMPP IS THE ARRAY OF COMPOSITION AT THE BEGINNING OF THE TIMESTEP.
!  THIS IS NEEDED FOR COMPOSITION DIFFUSION IN ROTATING MODELS.
       do j = 1,star%nz
          star%logP_start(j) = star%xa_start(i,j)
       end do
         call osplin(star%old_shell_mass,star%logT_start,star%log_mass,star%logP_start, &
              old_point_count,new_point_count)
       do j = 1,new_num_zones
          star%xa_start(i,j) = star%logT_start(j)
       end do
      end do


!  HCOMPM IS THE ARRAY OF CHANGES IN COMPOSITION DUE TO NUCLEAR BURNING.
!  THIS IS NEEDED FOR COMPOSITION DIFFUSION IN ROTATING MODELS.
      do i = 1,7
       do j = 1,star%nz
          star%logP_start(j) = star%rot%reaction_rate_by_zone(reaction_rate_species_index(i),j)
       end do
         call osplin(star%old_shell_mass,star%logT_start,star%log_mass,star%logP_start, &
              old_point_count,new_point_count)
       do j = 1,new_num_zones
          star%rot%reaction_rate_by_zone(reaction_rate_species_index(i),j) = &
               star%logT_start(j)
       end do
      end do
! MHP 05/02 IF THE SURFACE DEUTERIUM IS ABOVE
! THRESHOLD (1.0D-14) FIND THE NEW RUN OF
! DEUTERIUM BURNING RATES
      if (star%job%use_extended_composition .and. &
           star%xa(i_h2,star%nz).ge.1.0D-14) then
         do j = 1,star%nz
            star%logP_start(j) = star%deuterium_burning_rate_start(j)
         end do
         call osplin(star%old_shell_mass,star%logT_start,star%log_mass,star%logP_start, &
              old_point_count,new_point_count)
         do j = 1,new_num_zones
            star%deuterium_burning_rate_start(j) = star%logT_start(j)
         end do
      endif
! NOW FIND RUN OF P,R,L,T,AND RHO IN THAT ORDER FOR THE NEW POINTS.

      call osplin(star%old_shell_mass,star%logP_start,star%log_mass,star%logP, &
           old_point_count,new_point_count)
      call osplin(star%old_shell_mass,star%logR_start,star%log_mass,star%logR, &
           old_point_count,new_point_count)
      call osplin(star%old_shell_mass,star%luminosity_lsun_start,star%log_mass,star%luminosity_lsun, &
           old_point_count,new_point_count)
      call osplin(star%old_shell_mass,star%logT_start,star%log_mass,star%logT, &
           old_point_count,new_point_count)
      call osplin(star%old_shell_mass,star%logRho_start,star%log_mass,star%logRho, &
           old_point_count,new_point_count)

! FOR ROTATING MODELS FIND THE NEW RUN OF OMEGA,J/M,FP,FT,R0,AND ETA2.
      if (star%job%rotation_active) then
         call osplin(star%old_shell_mass,star%old_omega,star%log_mass,star%omega, &
              old_point_count,new_point_count)
         call osplin(star%old_shell_mass,star%old_specific_angular_momentum,star%log_mass, &
              star%j_rot,old_point_count,new_point_count)
         call osplin(star%old_shell_mass,fp_old,star%log_mass,star%fp_rot, &
              old_point_count,new_point_count)
         call osplin(star%old_shell_mass,ft_old,star%log_mass,star%ft_rot, &
              old_point_count,new_point_count)
         call osplin(star%old_shell_mass,star%old_mean_radius,star%log_mass,star%mean_radius, &
              old_point_count,new_point_count)
         call osplin(star%old_shell_mass,star%old_eta_squared,star%log_mass,star%eta_squared, &
              old_point_count,new_point_count)
      endif
!

!     SPIT OUT POINT DISTRIBUTION DETAILS IF REQUESTED
      if (star%ctrl%ldebug .and.  star%ctrl%npoint.lt.9999) then
      if (mod(star%model_number,star%ctrl%npoint).eq.0) then
         min_common_count = min0(star%nz,new_num_zones)
!
         write(idebug,910)
  910    format('1',20X,'OLD POINTS',54X,'NEW POINTS'/2(3X,'N',5X,'S', &
     8X,'P',7X,'T',7X,'R',8X,'L',7X,'X',4X,'Z',3X,'O16',1X) )
         write(idebug,920) (i,star%log_mass(i),star%logP(i), &
              star%logT(i),star%logR(i),star%luminosity_lsun(i), &
              x_new(i),z_new(i),star%xa(i_o16,i),i,star%old_shell_mass(i), &
              star%logP_start(i),star%logT_start(i),star%logR_start(i), &
              star%luminosity_lsun_start(i),star%xa(i_h1,i),star%xa(i_metals,i), &
              star%xa(i_o16,i), i = 1,min_common_count)
  920    format( 2(1X,I3,F11.7,F8.4,F8.5,F8.4,1PE9.2,0PF6.3,2F5.3) )
         if (star%nz.gt.min_common_count) then
            min_common_count = min_common_count + 1
            write(idebug,930) (i,star%log_mass(i),star%logP(i), &
                 star%logT(i),star%logR(i),star%luminosity_lsun(i), &
                 x_new(i),z_new(i),star%xa(i_o16,i),i=min_common_count, &
                 star%nz)
  930       format( 1X,I3,F11.7,F8.4,F8.5,F8.4,1PE9.2,0PF6.3,2F5.3)
         else if (new_num_zones.gt.min_common_count) then
            min_common_count = min_common_count + 1
            write(idebug,940)(i,star%old_shell_mass(i),star%logP_start(i), &
                 star%logT_start(i),star%logR_start(i),star%luminosity_lsun_start(i), &
                 star%xa(i_h1,i),star%xa(i_metals,i),star%xa(i_o16,i), &
                 i=min_common_count,new_num_zones)
  940       format(65X,I3,F11.7,F8.4,F8.5,F8.4,1PE9.2,0PF6.3,2F5.3)
         endif
      endif
      endif


! TRANSFER NEW POINTS.
      do j = 1,new_num_zones
       star%log_mass(j) = star%old_shell_mass(j)
       star%logP(j) = star%logP_start(j)
       star%logT(j) = star%logT_start(j)
       star%logR(j) = star%logR_start(j)
       star%luminosity_lsun(j) = star%luminosity_lsun_start(j)
       star%logRho(j) = star%logRho_start(j)
      end do
      if (star%job%rotation_active) then
       do j = 1, new_num_zones
          star%j_rot(j) = star%old_specific_angular_momentum(j)
          star%omega(j) = star%old_omega(j)
          star%fp_rot(j) = fp_old(j)
          star%ft_rot(j) = ft_old(j)
          star%eta_squared(j) = star%old_eta_squared(j)
          star%mean_radius(j) = star%old_mean_radius(j)
       end do
      endif
! MHP 6/00 INTERPOLATED IN ENERGY GENERATION AT START OF TIMESTEP
      if (star%job%rotation_active .or. (star%job%use_extended_composition .and. &
           star%job%envelope_overshoot_active)) then
         call osplin(star%old_shell_mass,star%rot%old_esum,star%log_mass,star%sesum, &
              old_point_count,new_point_count)
         do zone_index = 1,star%nz
            spline_y(zone_index) = star%sesum(zone_index)+star%seg(i_eps_neu,zone_index)+ &
                 star%seg(i_eps_grav,zone_index)
         end do
         call osplin(star%old_shell_mass,star%rot%old_eps,star%log_mass,spline_y, &
              old_point_count,new_point_count)
      endif
      write(short_file_unit,1020) star%nz,new_num_zones
 1020 format(' POINTS  OLD',I5,'   NEW',I5)
      star%nz = new_num_zones
! SET UP WEIGHTS AND MASSES
      mass_curr = dexp(clndp*star%log_mass(1))
      mass_prev = - mass_curr
      do i = 2,star%nz
       mass_two_back = mass_prev
       mass_prev = mass_curr
       mass_curr = dexp(clndp*star%log_mass(i))
       star%m(i-1) = mass_prev
       star%dm(i-1) = 0.5D0*(mass_curr-mass_two_back)
      end do
      star%m(star%nz) = mass_curr
      star%dm(star%nz) = dexp(ln10*star%log_total_mass) - &
           0.5D0*(mass_prev+mass_curr)
      if (star%job%rotation_active) then
!  FIRST GUESS AT MOMENT OF INERTIA(HI)
       do i=1,star%nz
          star%i_rot(i) = cc23*star%dm(i)* &
               dexp(ln10*2.0D0*star%logR(i))
       end do
!   CALCULATE OVERSHOOT
       call am_convective_regions(star%xa,star%logRho,star%logP,star%logR, &
            star%log_mass,star%logT,star%convective_flag,star%nz, &
            am_transport_convective_flag,radiative_zone_bounds, &
            convective_zone_bounds,num_radiative_zones,num_convective_zones)
! JNT 2025/09/03 duplicating 2015/04/06 recompute moment of interia
! before recomputing the rotation I am less confident that this is
! necessary since WALPCZ does run in this version but I don't think
! it can hurt.
       call zone_moments_of_inertia(star%eta_squared,star%logR,star%log_mass,star%dm,1,star%nz, &
            star%omega,star%mean_radius,star%i_rot,star%qiw)
! END JNT

!   FIND THE ANGULAR VELOCITY OMEGA THAT CORRESPONDS TO THE GIVEN
!   SPECIFIC ANGULAR MOMENTUM HJM.
       call omega_from_j(star%logRho,star%j_rot,star%logR, &
            star%log_mass,star%dm,am_transport_convective_flag,star%nz, &
            star%eta_squared,star%i_rot,star%omega,star%qiw,star%mean_radius)
!  CALCULATE FP,FT,R0 AND ETA2 GIVEN OMEGA
       call rotation_shape_factors(star%logRho,star%logR,star%log_mass,star%nz,star%omega, &
            star%eta_squared,star%fp_rot,star%ft_rot,star%mean_gravity,star%mean_radius)
!  FIND CORRECT MOMENT OF INERTIA(HI)
!        CALL MOMI(ETA2,HD,HR,HS,HS2,1,M,OMEGA,R0,HI,QIW,M)  ! KC 2025-05-31
       call zone_moments_of_inertia(star%eta_squared,star%logR,star%log_mass,star%dm,1,star%nz, &
            star%omega,star%mean_radius,star%i_rot,star%qiw)
!  FIND NEW TOTAL ANGULAR MOMENTUM
       sum_angular_momentum = 0.0D0
       sum_rotational_ke = 0.0D0
       do i = 1,star%nz
          angular_momentum_shell = star%j_rot(i)*star%dm(i)
          star%kinetic_energy_rot(i) = 0.5D0*star%omega(i)*angular_momentum_shell
          sum_angular_momentum = sum_angular_momentum + angular_momentum_shell
          sum_rotational_ke = sum_rotational_ke + star%kinetic_energy_rot(i)
       end do
       write(short_file_unit,1120)total_angular_momentum, &
            sum_angular_momentum,total_rotational_ke,sum_rotational_ke
 1120    format(1X,'TOTAL J OF STAR - PREVIOUS ',1PE21.13,' NEW ', &
     1PE21.13/' TOTAL ROTATIONAL K.E. OF STAR-PREVIOUS ',1PE21.13, &
     ' NEW ',1PE21.13)
       total_angular_momentum = sum_angular_momentum
       total_rotational_ke = sum_rotational_ke
!  STORE THE OLD MODEL STRUCTURE FOR USE IN DIFFUSION.
       do i = 1,star%nz
          star%old_omega(i) = star%omega(i)
          star%old_hg(i) = star%mean_gravity(i)
          star%old_moment_of_inertia(i) = star%i_rot(i)
          star%old_eta_squared(i) = star%eta_squared(i)
          star%old_mean_radius(i) = star%mean_radius(i)
          star%convective_flag_start(i) = star%convective_flag(i)
          star%cz_flag_start(i) = am_transport_convective_flag(i)
! MHP 10/91 J/M STORED IN HJX FOR I/O USE.
            star%old_specific_angular_momentum(i) = star%j_rot(i)
       end do
! MHP 9/91 CHANGE : T GRADIENTS STORED IF LEXCOM=T AND LOVSTE=T; OR FOR
! ROTATION; THIS IS NEEDED SO THAT THE BASE OF THE OVERSHOOT REGION FOR
! PRE-MS MODELS CAN BE ACCURATELY LOCATED.
      endif
      if (star%job%rotation_active .or. (star%job%use_extended_composition .and. &
           star%job%envelope_overshoot_active)) then
! END OF 9/91 CHANGE
!   FIND THE NEW RUN OF PHYSICAL VARIABLES AT THE NEW SET OF POINTS;
!   THIS IS NEEDED EVEN IN THE ABSENCE OF DIFFUSION TO ACCURATELY LOCATE
!   THE EDGES OF CONVECTION ZONES.
!        CALL PHYSIC(FP,FT,HCOMP,HD,HG,HL,HP,HR,HS,HT,LC,LCZ,M,TEFFL)  ! KC 2025-05-31
       call shell_physics(star%fp_rot,star%ft_rot,star%xa,star%logRho,star%mean_gravity,star%luminosity_lsun, &
            star%logP,star%logR,star%log_mass,star%logT, &
            star%convective_flag,star%nz,star%log_Teff, jerr)
       if (jerr /= 0) then
       ! 2026 (phase five, step B): propagate instead of stopping
          ierr = jerr
          return
       end if
!   FOR DIFFUSION STORE THE AUXILLARY QUANTITIES NEEDED TO CALCULATE
!   VELOCITIES AT THE START OF THE TIMESTEP WITH THE NEW POINT DISTRIBUTION
!   SO THAT A SERIES OF SMALL DIFFUSION TIMESTEPS CAN BE TAKEN WITHIN
!   ONE LARGE EVOLUTIONARY TIMESTEP.
         do zone_index = 1,star%nz
            star%rot%old_del_radiative_mix(zone_index) = star%del_grad(i_grad_rad,zone_index)
            star%rot%old_delm(zone_index) = star%del_grad(i_grad_actual,zone_index)
            star%rot%old_del_adiabatic_mix(zone_index) = star%del_grad(i_grad_ad,zone_index)
            star%rot%old_amu(zone_index) = star%mean_molecular_weight(zone_index)
            star%rot%old_om(zone_index) = star%so(zone_index)
            star%rot%old_cp(zone_index) = star%cp(zone_index)
            star%rot%old_qdt(zone_index) = star%qdt(zone_index)
            star%rot%old_vel(zone_index) = star%svel(zone_index)
            star%rot%old_visc(zone_index) = star%visc(zone_index)
            star%rot%old_thdif(zone_index) = star%thdif(zone_index)
! MHP 06/02
            star%rot%del_grad_diff_interface(zone_index) = &
                 star%rot%old_del_adiabatic_mix(zone_index) - star%rot%old_delm(zone_index)
! MHP 6/00 CALCULATED EARLIER
!            ESUMO(IM) = SESUM(IM)
            star%rot%max_domega_dr_old(zone_index) = star%rot%max_domega_dr(zone_index)
         end do
! MHP 06/02 ADDED TERM FOR THE TIME EVOLUTION
! OF THE ANGULAR VELOCITY DISTRIBUTION
         do i = 2,star%nz
            delta_radius = exp(ln10*star%logR(i))-exp(ln10*star%logR(i-1))
            delta_omega = star%omega(i) - star%omega(i-1)
            omega_mid = 0.5D0*(star%omega(i)+star%omega(i-1))
            log_factor = 2.0D0*(star%logR(i)+star%logR(i-1))-0.5D0* &
     (star%log_mass(i)+star%log_mass(i-1))-cgl
            star%rot%tho(i) = exp(ln10*log_factor)*omega_mid*delta_omega/delta_radius
            star%rot%qwrst(i) = delta_omega/delta_radius
         end do
      endif
!  CALCULATE NEW SURFACE OPACITY TABLE IF NEEDED.
      if (dabs(star%xnew-star%xa(i_h1,star%nz)).gt.1.0D-8) then
               star%xnew = star%xa(i_h1,star%nz)
               star%znew = star%xa(i_metals,star%nz)
               call kap_update_surface_tables(star%xnew)

      end if
end subroutine interpolate_onto_new_grid

end subroutine rezone
