!----------------------------------------------------------------------
! microdiff_etm
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original microdiff_etm.f; only variable names, source form, and
! comment style were updated. Validated against the Stage 0 regression
! suite (examples/run_standard_solar_model).
!
! "ETM" = equally-spaced-to-model-grid transform, the reverse of
! microdiff_mte.f90: interpolates the changes in hydrogen/metal/light-
! element abundance computed by microdiff_run.f90 (on the equally
! spaced grid) back onto the original model grid, applies the standard
! floor/ceiling failsafes, renormalizes helium via X+Y+Z=1, and
! converts DT/STOT/HRU/HTU/HS1/HQPR back out of Bahcall & Loeb units.
! Last stage of the microdiff.f90 pipeline (see also
! microdiff_setup.f90, microdiff_mte.f90, microdiff_coefficients.f90,
! microdiff_run.f90).
!  TRANSFORM BACK TO ORIGINAL GRID OF MODEL POINTS FROM EQUALLY
!  SPACED GRID.
subroutine microdiff_etm(timestep, eq_radius, eq_delta_hydrogen, &
     eq_delta_metal, eq_delta_light, zone_begin, zone_end, num_eq_points, &
     composition, dlnp_dr, radius_bl, enclosed_mass, temperature_bl, &
     num_zones, total_mass, num_light, light_element_id)
      use star_info_lib, only: star, json
      use numerics_lib, only: intrp2
      use microdiff_mte_lib, only: lagrange4
      implicit none

      double precision, intent(inout) :: timestep
      double precision, intent(in) :: eq_radius(json), &
           eq_delta_hydrogen(json), eq_delta_metal(json)
      integer, intent(in) :: num_light
      double precision, intent(in) :: eq_delta_light(num_light,json)
      integer, intent(in) :: zone_begin, zone_end, num_eq_points
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(inout) :: dlnp_dr(json), radius_bl(json), &
           enclosed_mass(json), temperature_bl(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: total_mass
      integer, intent(in) :: light_element_id(num_light)
! --- locals ---
      double precision :: tabler(4), facinterp(4)
      integer :: i, ii, j, jmin, k, k0, kk
      double precision :: x_min_floor, z_max, zz, zz2, radmod, dxmod, &
           x_max, dzmod

!  TRANSFORM BACK TO ORIGINAL GRID OF MODEL POINTS FROM EQUALLY
!  SPACED GRID.
      x_min_floor = 0.0d0
      do i = zone_begin,1,-1
         composition(1,i)=max(composition(1,i) + eq_delta_hydrogen(1),x_min_floor)
      end do
! MHP 3/94 ADDED METAL DIFFUSION
! NOTE THAT BECAUSE METALS SINK, AND HYDROGEN RISES, THE FAILSAFES
! ARE OPPOSITE (GUARDING AGAINST NEGATIVE X AND Z>1 RESPECTIVELY).
      if(star%job%use_diffusion_z)then
         do i = zone_begin,1,-1
            z_max = 1.0d0 - composition(1,i) - composition(4,i)
            zz=min(composition(3,i)+eq_delta_metal(1),z_max)
            zz2 = zz/composition(3,i)
            composition(3,i) = zz
            do j = 5,11
               composition(j,i) = zz2*composition(j,i)
            end do
            composition(2,i)=1.0d0-composition(1,i)-composition(3,i)-composition(4,i)
         end do
      else
         do i = zone_begin,1,-1
            composition(2,i)=1.0d0-composition(1,i)-composition(3,i)-composition(4,i)
         end do
      endif
! G SOMERS 5/15; ADD LIGHT ELEMENT DIFFUSION
      if(star%ctrl%diffuse_lithium)then
         do kk = 1,num_light
            x_min_floor = 0.0d0
            do i = zone_begin,1,-1
               ii = light_element_id(kk)
               composition(ii,i)=max(composition(ii,i) + eq_delta_light(kk,1),x_min_floor)
            end do
         end do
      endif
!
      jmin=2
      do i=zone_begin+1,zone_end-1
         do j=jmin,num_eq_points
!  FIND EQUALLY SPACED GRID POINTS CLOSEST TO THE MODEL POINT.
            if(eq_radius(j).ge.radius_bl(i))then
!  ENSURE THAT FIRST INTERP. POINT NO LESS THAN FIRST EQUALLY SPACED POINT.
               k0 = max(j-2,1)
!  ENSURE THAT LAST INTERP. POINT NO GREATER THAN LAST EQUALLY SPACED POINT.
               k0 = min(k0,num_eq_points-3)
! JVS fix for NPT = 3?
               if (k0 .eq. 0) k0=1
               jmin=j
               exit
            endif
         end do
         if (j > num_eq_points) then
         k0 = num_eq_points-3
         jmin=num_eq_points
         end if

         do k=1,4
            tabler(k)=eq_radius(k0+k-1)
         end do
         radmod=radius_bl(i)
!  FIND 4 POINT LAGRANGIAN INTERPOLATION FACTORS.
         call intrp2(tabler,facinterp,radmod)
!  PERFORM 4 POINT LAGRANGIAN INTERPOLATION FOR CHANGE IN X.
         dxmod = lagrange4(facinterp, eq_delta_hydrogen, k0)
         x_max = 1.0d0 - composition(3,i) - composition(4,i)
         composition(1,i)=min(composition(1,i) + dxmod,x_max)
! MHP 3/94 ADDED METAL DIFFUSION
         if(star%job%use_diffusion_z)then
            z_max = 1.0d0 - composition(1,i) - composition(4,i)
            dzmod = lagrange4(facinterp, eq_delta_metal, k0)
            zz = min(composition(3,i)+dzmod,z_max)
            zz2 = zz/composition(3,i)
            composition(3,i)=zz
            do j = 5,11
               composition(j,i) = zz2*composition(j,i)
            enddo
            composition(2,i)=1.0d0-composition(1,i)-composition(3,i)-composition(4,i)
         else
            composition(2,i)=1.0d0-composition(1,i)-composition(3,i)-composition(4,i)
         endif
! GES 5/15 ADDED LIGHT ELEMENT DIFFUSION
         if(star%ctrl%diffuse_lithium)then
            do kk = 1,num_light
               ii = light_element_id(kk)
               dxmod = facinterp(1)*eq_delta_light(kk,k0)+ &
                       facinterp(2)*eq_delta_light(kk,k0+1)+ &
                       facinterp(3)*eq_delta_light(kk,k0+2)+ &
                       facinterp(4)*eq_delta_light(kk,k0+3)
               composition(ii,i) = composition(ii,i) + dxmod
            enddo
         endif
      end do
!
      do i = zone_end,num_zones
         x_max = 1.0d0 - composition(3,i) - composition(4,i)
         composition(1,i)=min(composition(1,i) + eq_delta_hydrogen(num_eq_points),x_max)
      end do
! MHP 3/94 ADDED METAL DIFFUSION
      if(star%job%use_diffusion_z)then
         do i = zone_end,num_zones
            zz = max(composition(3,i)+eq_delta_metal(num_eq_points),0.0d0)
            zz2 = zz/composition(3,i)
            composition(3,i) = zz
            do j = 5,11
               composition(j,i) = zz2*composition(j,i)
            end do
            composition(2,i)=1.0d0-composition(1,i)-composition(3,i)-composition(4,i)
         end do
      else
         do i = zone_end,num_zones
            composition(2,i)=1.0d0-composition(1,i)-composition(3,i)-composition(4,i)
         end do
      endif
! GES 5/15 LIGHT ELEMENT DIFFUSION
      if(star%ctrl%diffuse_lithium)then
         do kk = 1,num_light
            do i = zone_end,num_zones
               ii = light_element_id(kk)
               composition(ii,i)=composition(ii,i) + eq_delta_light(kk,num_eq_points)
            end do
         enddo
      endif
!
      do i=1,num_zones
         radius_bl(i)=radius_bl(i)/star%bl_radius_scale
         temperature_bl(i)=temperature_bl(i)/star%bl_temp_scale
         enclosed_mass(i)=enclosed_mass(i)/star%bl_mass_scale
         dlnp_dr(i)=dlnp_dr(i)*star%bl_radius_scale
      end do
      timestep=timestep*star%bl_time_scale
      total_mass=total_mass/star%bl_mass_scale
      return
end subroutine microdiff_etm
