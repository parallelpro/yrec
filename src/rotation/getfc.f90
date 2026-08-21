!----------------------------------------------------------------------
! getfc
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original getfc.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! This subroutine solves for the variable FC using the prescription
! of Zahn (1992). FC = MIN(1, C/30*U/(2V-ALPHA*U)), where
! ALPHA = 1/2 D(LN J/M)/D LN R and
! V = 1/(6*RHO*R) *D/DR(RHO* R**2 * U) (U = diffusion velocity).
! If FC = 0, skip calculation (no mixing).
subroutine getfc(log_density, radius, diffusion_velocity, zone_min, &
     zone_max, angular_velocity)

      use const_lib
      implicit none
      integer, parameter :: json = 5000

! INPUT VARIABLES :
! log_density - LOG10 OF DENSITY(GM/CM**3) AT ZONE CENTERS.
! radius - UNLOGGED RUN OF RADII AT ZONE CENTERS.
! diffusion_velocity - UNLOGGED RUN OF DIFFUSION VELOCITIES AT ZONE EDGES.
! zone_min,zone_max - STARTING AND ENDING SHELLS.
! angular_velocity - RUN OF ANGULAR VELOCITIES AT ZONE CENTERS.
      double precision, intent(in) :: log_density(json), radius(json), &
           diffusion_velocity(json)
      integer, intent(in) :: zone_min, zone_max
      double precision, intent(in) :: angular_velocity(json)


! common/intvar/: only dm (RHO at zone edges) and interface_radius (RM,
! the run of radii at zone edges) are used here; the rest are unused
! placeholders. Naming matches vcirc.f90.
      double precision :: interface_luminosity(json), delami(json), &
           delmi(json), dm(json), epsilm(json), interface_gravity_factor(json), &
           hs3(json), pm(json), qdtmi(json), interface_radius(json), tm(json)
      common/intvar/ interface_luminosity, delami, delmi, dm, epsilm, &
           interface_gravity_factor, hs3, pm, qdtmi, interface_radius, tm



      double precision :: half_dlnj_dlnr(json), circ_velocity(json)
      save

      integer :: zone_index
      double precision :: omega_mid
      double precision :: interp_weight
      double precision :: v_minus, v_plus
      double precision :: density_r2_v_minus, density_r2_v_plus
      double precision :: density_radius, d_density_r2, du_dr
      double precision :: denom_test

! IF FC = 0, SKIP CALCULATION (NO MIXING).

      if (mixing_velocity_scale.le.0.0D0) then
         do zone_index = zone_min, zone_max
            vfc(zone_index) = 0.0D0
         end do
         goto 9999
      endif
! DETERMINE ALPHA.
      do zone_index = zone_min, zone_max
         omega_mid = 0.5D0*(angular_velocity(zone_index)+ &
              angular_velocity(zone_index-1))
         half_dlnj_dlnr(zone_index) = (angular_velocity(zone_index)* &
              radius(zone_index)**2-angular_velocity(zone_index-1)* &
              radius(zone_index-1)**2)/ &
              (radius(zone_index)-radius(zone_index-1))/ &
              (omega_mid*interface_radius(zone_index))
      end do
! DETERMINE 2V.
      interp_weight = (radius(zone_min)-interface_radius(zone_min))/ &
           (interface_radius(zone_min+1)-interface_radius(zone_min))
      v_minus = diffusion_velocity(zone_min)+interp_weight* &
           (diffusion_velocity(zone_min+1)-diffusion_velocity(zone_min))
      density_r2_v_minus = exp(ln10*log_density(zone_min))* &
           radius(zone_min)**2*v_minus
      do zone_index = zone_min + 1, zone_max - 1
         interp_weight = (radius(zone_index)-interface_radius(zone_index))/ &
              (interface_radius(zone_index+1)-interface_radius(zone_index))
         v_plus = diffusion_velocity(zone_index)+interp_weight* &
              (diffusion_velocity(zone_index+1)-diffusion_velocity(zone_index))
         density_r2_v_plus = exp(ln10*log_density(zone_index))* &
              radius(zone_index)**2*v_plus
         density_radius = dm(zone_index)*interface_radius(zone_index)
         circ_velocity(zone_index) = cc13*(density_r2_v_plus-density_r2_v_minus)/ &
              density_radius/(radius(zone_index)-radius(zone_index-1))
         density_r2_v_minus = density_r2_v_plus
      end do
! FOR FIRST AND LAST SHELLS, CALCULATE U*D/DR(RHO*R**2) AND ASSUME
! DU/DR = VALUE DERIVED FROM THE LAST NEIGHBOR SHELL WITH NONZERO V.
      density_radius = dm(zone_min)*interface_radius(zone_min)
      d_density_r2 = (exp(ln10*log_density(zone_min))*radius(zone_min)**2 - &
           exp(ln10*log_density(zone_min-1))*radius(zone_min-1)**2)/ &
           (radius(zone_min) - radius(zone_min-1))
      du_dr = (diffusion_velocity(zone_min+1)-diffusion_velocity(zone_min))/ &
           (interface_radius(zone_min+1)-interface_radius(zone_min))
      circ_velocity(zone_min) = cc13*(diffusion_velocity(zone_min)*d_density_r2/ &
           density_radius + interface_radius(zone_min)*du_dr)
      density_radius = dm(zone_max)*interface_radius(zone_max)
      d_density_r2 = (exp(ln10*log_density(zone_max))*radius(zone_max)**2 - &
           exp(ln10*log_density(zone_max-1))*radius(zone_max-1)**2)/ &
           (radius(zone_max) - radius(zone_max-1))
      du_dr = (diffusion_velocity(zone_max)-diffusion_velocity(zone_max-1))/ &
           (interface_radius(zone_max)-interface_radius(zone_max-1))
      circ_velocity(zone_max) = cc13*(diffusion_velocity(zone_max)*d_density_r2/ &
           density_radius + interface_radius(zone_max)*du_dr)
! NOW COMPUTE RUN OF FC; THIS ASSUMES THAT FC = SAME FOR ALL
! MECHANISMS AND IS LIMITED TO A MAXIMUM OF 1.
      do zone_index = zone_min, zone_max
         if (diffusion_velocity(zone_index).le.0.0D0) then
            vfc(zone_index) = 0.0D0
         else
            denom_test = max(abs(circ_velocity(zone_index)- &
                 half_dlnj_dlnr(zone_index)*diffusion_velocity(zone_index))/ &
                 diffusion_velocity(zone_index),mixing_velocity_scale)
            vfc(zone_index) = mixing_velocity_scale/denom_test
         endif
      end do
 9999 continue

      return
end subroutine getfc
