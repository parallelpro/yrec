!----------------------------------------------------------------------
! tpgrad
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original tpgrad.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! COMPUTES RADIATIVE GRADIENT AND COMPARES WITH ADIABATIC GRADIENT
! COMPUTES CONVECTIVE GRADIENT VIA MIXING LENGTH THEORY IF APPLICABLE
! ASSUMES eos_get AND kap_get HAVE BEEN CALLED (eos_res, kap_res)
! RETURNS is_convective = .TRUE. IF CONVECTIVE
!         radiative_gradient = RADIATIVE GRADIENT
!         actual_gradient = ACTUAL GRADIENT
!          dgrad_dt_component,dgrad_dp_component = NAT-LOG DERIVATIVES OF
!          THE CONVECTIVE GRADIENT
!
! luminosity_lsun (originally B) is the UNLOGGED luminosity in solar
! units: the DELR formula below -- O*B*exp(ln10*(PL-SL-4*TL+CLSUNL-CGL+
! CDELRL))*FTL/FPL, i.e. kappa * (L/Lsun) * Lsun_cgs * P / (G*M*T**4),
! the standard radiative-gradient formula -- only balances dimensionally
! if B multiplies in linearly. The callers pass star%luminosity_lsun.
! 2026 de-tramp (ROADMAP item 3): 33 arguments -> 19. The former
! positional form (14 unpacked eos/kap scalars) merged with its
! result-array wrapper temperature_gradients_r: this single entry now
! takes the eos_get/kap_get result arrays directly and unpacks the
! named-index values into locals at the top -- the same doubles, so
! byte-identical. All per-zone state stays an explicit argument
! (NEVER absorb star% per-zone reads here: this routine is on the
! star-blind envint kernel's integrand path, where the state comes
! from the integration point, not from star%'s zone arrays).
! zone_index (2026 readability W2): the Henyey zone number im of the
! point being evaluated, or -1 for an envelope/atmosphere point; it
! gates the adiabatic-overshoot branch (ladov, iov1 <= im <= iov2).
! Before W2 the callers wrote it to star%iovim just before the call.
! It is optional only because core/read_starting_model.f90 still
! passes it through star%iovim = -1 (see the W2 sidechan report,
! Deferred); when absent, star%iovim is read exactly as before.
module temperature_gradients_lib
      implicit none
contains

subroutine temperature_gradients(log_temperature, log_pressure, &
     eos_res, kap_res, log_radius, log_mass, luminosity_lsun, &
     actual_gradient, radiative_gradient, dgrad_dt_component, &
     dgrad_dp_component, dgrad_dr_component, convective_velocity, &
     want_derivatives, is_convective, pressure_rotation_factor, &
     temperature_rotation_factor, log_teff, ierr, zone_index)

      use eos_lib, only: num_eos_results, i_temperature, i_pressure, &
           i_density, i_dlnrho_dlnt, i_dlnrho_dlnp, i_cp, i_grada, &
           i_dlnrho_dlnt_dt, i_dlnrho_dlnp_dt, i_grada_dt, i_grada_dp, &
           i_cp_dt, i_cp_dp
      use kap_lib, only: num_kap_results, i_kap, i_dlnkap_dlnrho, &
           i_dlnkap_dlnt
      use star_info_lib, only: star
      use luout_lib
      use phys_const_lib
      use math_lib
      implicit none

      double precision, intent(in) :: log_temperature, &
           log_pressure, log_radius, log_mass, luminosity_lsun
      double precision, intent(in) :: eos_res(num_eos_results), &
           kap_res(num_kap_results)
      double precision, intent(out) :: actual_gradient, radiative_gradient
      double precision, intent(out) :: dgrad_dt_component, &
           dgrad_dp_component, dgrad_dr_component
      double precision, intent(out) :: convective_velocity
! named-index values from the result arrays, unpacked below (same
! doubles the 33-argument form received positionally)
      double precision :: temperature, pressure, density, opacity, &
           dlnrho_dlnt, dlnrho_dlnp, dlnkap_dlnt, dlnkap_dlnrho, &
           specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
           dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp
      logical, intent(in) :: want_derivatives
      logical, intent(out) :: is_convective
      double precision, intent(in) :: pressure_rotation_factor, &
           temperature_rotation_factor, log_teff
      double precision, parameter :: vtol=1.0d-10
! max_cubic_iterations: Newton iteration cap for the MLT cubic.
      integer, parameter :: max_cubic_iterations = 25
! superadiabatic_tol: the zone is radiative unless del_rad exceeds del_ad
! by more than this.
      double precision, parameter :: superadiabatic_tol = 1.0d-6
      integer :: iter
      double precision :: deldel, g, presht, phi, phi2, phiphi, test, a1, &
           v, a3, a3p, vp, vd, ddel, qdelat, qdelap, tempot, &
           tempop, qddelt, qddelp, temp1, qa1t, qa1p, qa1r, qa3t, qa3p, &
           qa3r, temp2, temp3, qvt, qvp, qvr, deli, ateffl, deepx

      integer, intent(out) :: ierr
      integer, intent(in), optional :: zone_index
! overshoot_zone_index: zone_index when given, else star%iovim (the
! pre-W2 side channel, still used by core/read_starting_model).
      integer :: overshoot_zone_index

      ierr = 0
      if (present(zone_index)) then
         overshoot_zone_index = zone_index
      else
         overshoot_zone_index = star%iovim
      end if

      temperature = eos_res(i_temperature)
      pressure = eos_res(i_pressure)
      density = eos_res(i_density)
      dlnrho_dlnt = eos_res(i_dlnrho_dlnt)
      dlnrho_dlnp = eos_res(i_dlnrho_dlnp)
      specific_heat_cp = eos_res(i_cp)
      adiabatic_gradient = eos_res(i_grada)
      dlnrho_dlnt_dt = eos_res(i_dlnrho_dlnt_dt)
      dlnrho_dlnp_dt = eos_res(i_dlnrho_dlnp_dt)
      adiabatic_gradient_dt = eos_res(i_grada_dt)
      adiabatic_gradient_dp = eos_res(i_grada_dp)
      specific_heat_cp_dt = eos_res(i_cp_dt)
      specific_heat_cp_dp = eos_res(i_cp_dp)
      opacity = kap_res(i_kap)
      dlnkap_dlnt = kap_res(i_dlnkap_dlnt)
      dlnkap_dlnrho = kap_res(i_dlnkap_dlnrho)

! star%alfmlt/phmlt/cmxmlt are zeroed here and assigned nowhere else
! (they reach the valfmlt/vphmlt/vcmxmlt output columns as zeros); see
! audit/readability-sweep-2026-09-03/SUMMARY.md 1.1 #10 -- author decision.
      star%alfmlt=0.0d0
      star%phmlt=0.0d0
      star%cmxmlt=0.0d0
      radiative_gradient = opacity*luminosity_lsun*exp(ln10*(log_pressure - &
           log_mass - 4d0*log_temperature + star%log10_solar_luminosity - cgl + &
           cdelrl))* &
           temperature_rotation_factor/pressure_rotation_factor
      deldel = radiative_gradient - adiabatic_gradient
      if(deldel.le.superadiabatic_tol) then
! ZONE IS RADIATIVE
       is_convective = .false.
       actual_gradient = radiative_gradient
       convective_velocity=0.0d0
         if (star%ctrl%ladov .and. overshoot_zone_index.ge.star%iov1 .and. &
             overshoot_zone_index.le.star%iov2 .and. overshoot_zone_index.ne.-1) then
            actual_gradient = adiabatic_gradient
         end if
       return
      endif
! ZONE IS CONVECTIVE
      is_convective = .true.
      actual_gradient = adiabatic_gradient
! SKIP MIXING LENGTH THEORY FOR CORES
      if(log_temperature.gt.star%ctrl%tgcut) then
         convective_velocity = 1.0d-11
       if(want_derivatives) then
! DERIVATIVES OF CONVECTIVE GRADIENT
          dgrad_dt_component = adiabatic_gradient_dt
          dgrad_dp_component = adiabatic_gradient_dp
          dgrad_dr_component = 0.0d0
       endif
       return
      endif
! G Somers 9/14, Add the ability to include spots, which alter
! the radiative flux in the convective reigons. This is done by
! setting the flux to Lspotted = L/[f*x^4 + (1-f)], where f is the
! spot filling factor and x is the temperature contrast between
! the spotted surface and the normal surface (0 <= x <= 1). We also
! alter the surface boundary condition used in atm/atm_lib.f90, by looking
! up the pressure at the un-spotted T (ATEFFL) instead of Teff.
!
! We have also included the ability to have the temperature contrast
! SPOTX vary with depth, as is the case in stars. In this case, SPOTX
! is locally determined by X(T) = 1 - (1-Xsurf)(ATeff/T), where ATeff
! is the surface temperature of the unspotted regions. This is the
! analytic expression resulting from a user specified Xsurf, the
! assumption of magnetic equipartition at the surface (Pmag = Psurf/2),
! and the assumption of constant Pmag with depth.
!
! This flux alters DELR, so recalculate DELDEL with the correction. Only
! do this if the spot filling factor is non-zero, and the envelope is
! convective.
      if(star%ctrl%spot_filling_factor .ne. 0.00)then
         if(star%ctrl%spot_depth_varies)then
            ateffl = log_teff - 0.25*log10(star%ctrl%spot_filling_factor * &
                 pow(star%ctrl%spot_temp_contrast, 4.0) + 1.0 - star%ctrl%spot_filling_factor)
            deepx = 1.0 - (1.0 - star%ctrl%spot_temp_contrast)*(exp10(ateffl))/(exp10(log_temperature))
         else
            deepx = star%ctrl%spot_temp_contrast
         endif
         deldel = radiative_gradient/(star%ctrl%spot_filling_factor * pow(deepx, 4.0) + &
              1.0 - star%ctrl%spot_filling_factor) - adiabatic_gradient
      endif
! G Somers END
      g = exp(ln10*(cgl + log_mass - log_radius - log_radius))
      presht = pressure/(density*g)
      phi = star%mixing_length_alpha*density*opacity*presht
      phi2 = phi*phi
      phiphi = 1.0d0/(1.0d0 + cmixl2*phi2)
! SOLVE CUBIC A3*V**3 + V**2 + A1*V - 1.0 = 0.0
! ENSURE THAT THE SQUARE ROOT BEING TAKEN IS THAT OF A
! POSITIVE NUMBER
      test = deldel*g*(-dlnrho_dlnt)/presht
      if(test.gt.0.0d0)then
         a1 = cmixl3*phiphi*temperature**3*opacity/(specific_heat_cp*dsqrt(test))
      else
! ODD CONVERGENCE PROBLEM; WRITE DIAGNOSTIC AND CONTINUE.
         write(*,201)deldel,dlnrho_dlnt,test
 201     format('***WARNING***'/'TPGRAD TRIED TO TAKE SQUARE ROOT', &
                ' OF NEGATIVE NUMBER - ASSUMED ADIABATIC CONVECTION'/ &
                ' DELDEL ',1pe12.3,' QDT ',e12.3,' DISC. ',e12.3)
      write(*,911)opacity,luminosity_lsun,log_pressure,log_mass,log_temperature, &
           temperature_rotation_factor,pressure_rotation_factor,radiative_gradient, &
           adiabatic_gradient,log_radius
 911  format(1p10e12.3)
          ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the driver-side
          ! call sites (core/main, core/crrect, core/starin, setup/hpoint)
          ! preserve the historical stop on a nonzero return.
          ierr = 1
          return
      endif
      v = 1.0d0/a1
      a3 = 0.75d0*phi2*phiphi/a1
      a3p = 3.0d0*a3
      if(a3.gt.1.0d+3) v = pow(a3, (-0.333333333d0))
      do iter = 1,max_cubic_iterations
       v = dmin1(v,1.0d0)
       vp = a1 + v*(2.0d0 + v*a3p)
       vd = (-1.0d0 + v*(a1 + v*(1.0d0 + v*a3)))/vp
       vd = vd*(1.0d0 + vd*(1.0d0 + v*a3p)/vp)
       v = v - vd
       if(dabs(vd).lt.vtol) exit
      end do
      if (iter > max_cubic_iterations) then
      write(run_log_unit,20) log_pressure,log_temperature,opacity, &
           specific_heat_cp,dlnrho_dlnt
   20 format(' -----CUBIC NON-CONVERGENCE(PL,TL,CAPPA,CP,QDT)=' &
             ,2f10.6/1p3e12.6)
! ASSUME ADIABATIC CONVECTION
      convective_velocity = 1.0d-11
      if(want_derivatives) then
! DERIVATIVES OF CONVECTIVE GRADIENT
       dgrad_dt_component = adiabatic_gradient_dt
         dgrad_dp_component = adiabatic_gradient_dp
         dgrad_dr_component = 0.0d0
      endif
      return
      end if
      ddel = deldel*v*(v+a1)
      actual_gradient = adiabatic_gradient + ddel
! CALCULATE CONVECTIVE VELOCITY. THE RADICAND HAS THE SAME SIGN AS THE
! ONE ALREADY CHECKED ABOVE (deldel*g*(-dlnrho_dlnt)/presht > 0), SO NO
! SEPARATE SIGN TEST IS NEEDED HERE.
      convective_velocity = v*dsqrt(g*(-dlnrho_dlnt)*(star%mixing_length_alpha**2)*presht*deldel/8.0d0)
      if(want_derivatives) then
! DERIVATIVES OF CONVECTIVE GRADIENT
       qdelat = adiabatic_gradient_dt*adiabatic_gradient
       qdelap = adiabatic_gradient_dp*adiabatic_gradient
       tempot = dlnkap_dlnrho*dlnrho_dlnt + dlnkap_dlnt
       tempop = dlnkap_dlnrho*dlnrho_dlnp
       qddelt = ((-4.0d0 + tempot)*radiative_gradient - qdelat)/deldel
       qddelp = ((+1.0d0 + tempop)*radiative_gradient - qdelap)/deldel
       temp1 = 2.0d0*cmixl2*phi2*phiphi
       qa1t = tempot*(1d0-temp1)-specific_heat_cp_dt-0.5d0*(dlnrho_dlnt+ &
            qddelt+dlnrho_dlnt_dt/dlnrho_dlnt)+3d0
       qa1p = tempop*(1d0-temp1)-temp1-specific_heat_cp_dp-0.5d0 &
                *(dlnrho_dlnp+qddelp+dlnrho_dlnp_dt/dlnrho_dlnt) + 0.5d0
       qa1r = -2.0d0*temp1 + 2.0d0
       qa3t = 2.0d0*phiphi*tempot - qa1t
       qa3p = 2.0d0*phiphi*(1.0d0+tempop) - qa1p
       qa3r = 4.0d0*phiphi - qa1r
       temp1 = a1*v
       temp2 = a3*v*v*v
       temp3 = 1.0d0/vp
       qvt = -(temp1*qa1t + temp2*qa3t)*temp3
       qvp = -(temp1*qa1p + temp2*qa3p)*temp3
       qvr = -(temp1*qa1r + temp2*qa3r)*temp3
       temp1 = deldel*(a1+v+v)
       temp2 = deldel*a1*v
       deli = 1.0d0/actual_gradient
       dgrad_dt_component = (qdelat+qddelt*ddel+temp1*qvt+temp2*qa1t)*deli
       dgrad_dp_component = (qdelap+qddelp*ddel+temp1*qvp+temp2*qa1p)*deli
       dgrad_dr_component = (temp1*qvr+temp2*qa1r)*deli
      endif

      return
end subroutine temperature_gradients

end module temperature_gradients_lib
