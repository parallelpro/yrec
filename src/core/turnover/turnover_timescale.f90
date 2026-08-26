!----------------------------------------------------------------------
! gettau
!----------------------------------------------------------------------
! G Somers, 3/17; 2026 stitched-model restructure.
!----------------------------------------------------------------------
! Convective turnover timescale, one pressure scale height above the
! base of the surface convection zone -- computed from the STITCHED
! model (core/stitched_model.f90): the converged interior + the
! fixed-step envelope re-integration, walked as one structure by
! turnover_from_interior_new.
!
! History: this routine used to run its own envelope integration and
! hand-stitch env_struct onto the interior ("borrowed from stitch.f")
! on every call -- and its callers (getw mid-step, wrtout at write
! time) each triggered further integrations. The stitched model is
! built once per step in evolve_step, ahead of compute_observables;
! this routine now only walks it. star%pphot is owned by the stitched
! build (set by its atm_get); the mid-step consumers (wind
! saturation, the deuterium limiter) read the stored step-start
! values, which is what the *_old/pphot0 lag bookkeeping always
! assumed.
subroutine compute_turnover_timescale(radius_at_bcz)
      use stitched_model_lib, only: stx_prof, n_ie, &
           ip_mass, ip_logR, ip_logRho, ip_logP, ip_conv, ip_gradr, &
           ip_grada, ip_conv_vel
      use star_info_lib, only: star, json
      use phys_const_lib
      implicit none

      double precision, intent(out) :: radius_at_bcz

      double precision :: combined_radius(json), combined_pressure(json), &
           combined_density(json), combined_mass(json), &
           combined_gravity(json), combined_velocity(json)
      double precision :: combined_grad1(json), combined_grad2(json)
      logical :: combined_convective_flag(json)
      integer :: j, num_points

! TAUCZ = 0.0 (pphot is NOT reset here -- the stitched build owns it)
      star%convective_turnover_timescale = 0.0

! Span selection, as the pre-restructure code: find the top of the
! surface convection zone in the INTERIOR; if it sits more than one
! pressure scale height below the interior's surface (or the star is
! fully convective), the walk stays interior-only -- otherwise the
! stitched envelope is appended (regions 1-2 of the stitched grid;
! the atmosphere is never included). The walker takes json-sized
! arrays; clamp defensively (the historical assembly had the same
! bound).
      do j = star%nz-1, 1, -1
         if (.not. star%convective_flag(j)) exit
      end do
      if (star%logP(j+1) - star%logP(star%nz) .lt. 1.0d0) then
         num_points = min(n_ie, json)
      else
         num_points = star%nz
      end if
      do j = 1, num_points
         combined_mass(j) = stx_prof(ip_mass,j)*star%solar_mass_cgs
         combined_radius(j) = stx_prof(ip_logR,j)
         combined_pressure(j) = stx_prof(ip_logP,j)
         combined_density(j) = stx_prof(ip_logRho,j)
         combined_gravity(j) = combined_mass(j)* &
              exp(ln10*(cgl - 2.0d0*combined_radius(j)))
         combined_velocity(j) = stx_prof(ip_conv_vel,j)
         combined_grad1(j) = stx_prof(ip_gradr,j)
         combined_grad2(j) = stx_prof(ip_grada,j)
         combined_convective_flag(j) = stx_prof(ip_conv,j) .ne. 0.0d0
      end do

!          CALL TAUINTNEW(HCOMPF,HS2,HSF,LCF,HRF,HPF,HDF,HGF,MM,M,HVF,
!      *                  DELF1,DELF2,HSTOT,RBCZ)  ! KC 2025-05-31
      call turnover_from_interior_new(combined_mass, combined_convective_flag, &
           combined_radius, combined_pressure, combined_density, &
           combined_gravity, num_points, star%nz, combined_velocity, &
           combined_grad1, combined_grad2, radius_at_bcz)

      return
end subroutine compute_turnover_timescale
