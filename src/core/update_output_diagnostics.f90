!----------------------------------------------------------------------
! update_output_diagnostics
!----------------------------------------------------------------------
! Added 2026 (phase four, step 5 -- ROADMAP.md "Phase four: the star
! layer"). The state-computing blocks that historically lived inside
! io/wrtout.f90, moved verbatim so the writer only reads. In original
! wrtout order:
!   1. the luminosity-breakdown renormalization (MUTATES the model --
!      star%luminosity_breakdown -- skipped during He flash);
!   2. the core-CZ mass and the (dead) core-boundary radius; NOTE the
!      preserved original bug: core_boundary_log_radius reads
!      envelope_boundary_fx, which at that point still holds its
!      SAVEd value from the PREVIOUS call (the FX/FX2 typo in the
!      original wrtout.f); the blanket save below preserves exactly
!      that behavior, and core_boundary_radius remains computed but
!      never consumed, as in the original;
!   3. the central conditions: P/T extrapolated to the center, a real
!      eos_get evaluation for central beta and degeneracy eta, all
!      stored in star%run%central_*;
!   4. the surface-CZ base interpolation filling
!      star%run%envelope_mass/envelope_radius/envelope_cz_*;
!   5. NOT the gettau call: it stays in wrtout at its original spot.
!      gettau's atm-side integration prints progress diagnostics into
!      the .short stream, so hoisting it above wrtout's header writes
!      reorders that stream (values identical, layout not) -- blocked
!      by print interleaving, not by data flow. Documented residual.
! Called from core/main.f90 immediately before wrtout, once per
! output model, exactly the cadence the blocks had inside wrtout.
subroutine update_output_diagnostics(ierr)
      use star_info_lib
      use eos_lib
      use luout_lib
      use const_lib

      implicit none


! --- snu coefficients (as in wrtout) for the MESA-mode block below ---
      double precision :: clsnuf_diag(8), gasnuf_diag(8)
      data gasnuf_diag/1.18D1,2.15D2,7.14D4,7.17D1,2.40D4,6.04D1, &
           1.137D2,1.139D2/
      data clsnuf_diag/0.0D0,1.6D1,4.26D4,2.4D0,1.14D4,1.7D0,6.8D0,6.9D0/
! --- locals carried over from wrtout (names unchanged) ---
      double precision :: total_luminosity_sum, temp_value
      integer :: i
      double precision :: core_boundary_fx2, envelope_boundary_fx
! core_boundary_log_radius/core_boundary_radius (CORERL/CORER) are
! computed but never consumed (original behavior, preserved).
      double precision :: core_boundary_log_radius, core_boundary_radius
      double precision :: pressure_linear, log_pressure_center, &
           log_temperature_center, log_density_center, &
           hydrogen_fraction_center, metal_fraction_center
! temperature_linear_center/density_linear_center are separate
! eqstat/meqos output slots, never read after the call (dead output)
! but kept distinct so no argument slots alias.
      double precision :: temperature_linear_center, density_linear_center
      logical :: is_atmosphere_point, compute_derivatives
      double precision :: beta_center, beta_inverse_center, beta14_center, &
           mean_molecular_weight_center, amu_center, &
           electron_mean_molecular_weight_center, degeneracy_eta_center
      double precision :: qdt_center, qdp_center, qcp_center, dela_center, &
           qdtt_center, qdtp_center, qat_center, qap_center, qcpt_center, &
           qcpp_center
      integer :: ksaha_center
      double precision :: fxion(3)
      double precision :: dd1, dd2, cz_base_mass
      double precision :: envelope_cz_log_temperature, &
           envelope_cz_log_density, envelope_cz_log_pressure
! blanket save, as in wrtout: preserves envelope_boundary_fx's
! previous-call value (the documented bug above) and ksaha_center's
! saha-table state continuity across calls.
      save
   ! INTENTIONAL: incl. the documented FX/FX2 stale carry; reset via output_diag_reset_pending
!  RENORMALIZE LUMINOSITY TERMS TLUMX - SKIPPED FOR HE FLASH
      integer, intent(out) :: ierr

      ierr = 0

! 2026 (phase five, step C): see evolve_step's matching block. This
! includes envelope_boundary_fx, whose previous-call stale value is
! the documented FX/FX2 quirk -- a fresh process starts it at zero,
! so a repeated call must too.
      if (output_diag_reset_pending) then
         total_luminosity_sum = 0.0d0
         temp_value = 0.0d0
         core_boundary_fx2 = 0.0d0
         envelope_boundary_fx = 0.0d0
         core_boundary_log_radius = 0.0d0
         core_boundary_radius = 0.0d0
         pressure_linear = 0.0d0
         log_pressure_center = 0.0d0
         log_temperature_center = 0.0d0
         log_density_center = 0.0d0
         hydrogen_fraction_center = 0.0d0
         metal_fraction_center = 0.0d0
         temperature_linear_center = 0.0d0
         density_linear_center = 0.0d0
         beta_center = 0.0d0
         beta_inverse_center = 0.0d0
         beta14_center = 0.0d0
         mean_molecular_weight_center = 0.0d0
         amu_center = 0.0d0
         electron_mean_molecular_weight_center = 0.0d0
         degeneracy_eta_center = 0.0d0
         qdt_center = 0.0d0
         qdp_center = 0.0d0
         qcp_center = 0.0d0
         dela_center = 0.0d0
         qdtt_center = 0.0d0
         qdtp_center = 0.0d0
         qat_center = 0.0d0
         qap_center = 0.0d0
         qcpt_center = 0.0d0
         qcpp_center = 0.0d0
         dd1 = 0.0d0
         dd2 = 0.0d0
         cz_base_mass = 0.0d0
         envelope_cz_log_temperature = 0.0d0
         envelope_cz_log_density = 0.0d0
         envelope_cz_log_pressure = 0.0d0
         fxion = 0.0d0
         ksaha_center = 0
         is_atmosphere_point = .false.
         compute_derivatives = .false.
         output_diag_reset_pending = .false.
      end if

      if(.not.helium_flash_active) then
       total_luminosity_sum = star%luminosity_breakdown(i_lum_pp1)+star%luminosity_breakdown(i_lum_pp2)+ &
            star%luminosity_breakdown(i_lum_pp3)+star%luminosity_breakdown(i_lum_cno)+star%luminosity_breakdown(i_lum_3alpha)+ &
            star%luminosity_breakdown(i_lum_neu)+star%luminosity_breakdown(i_lum_grav)+star%luminosity_breakdown(i_lum_he_c)
       temp_value = star%luminosity_lsun(star%nz)/total_luminosity_sum
       do i = 1,8
          star%luminosity_breakdown(i) = star%luminosity_breakdown(i)*temp_value
       end do
      endif

!  CALCULATE MASS OF CENTRAL AND SURFACE CONVECTION ZONES
!  THESE MASSES ARE IN SOLAR UNITS
      if(star%core_cz_top_index.gt.1) then
       star%run%core_cz_mass = star%m(star%core_cz_top_index)/solar_mass_cgs
      else
       star%run%core_cz_mass = 0.0D0
      endif

! JVS 10/11 Be more care about the true boundary of the convective core
      if (star%core_cz_top_index.gt.1) then
! JVS 10/11 note: this formula reads envelope_boundary_fx (FX), which
! at this point has not yet been assigned in this call (it is set
! further below, in the JENV block) -- with SAVE it carries over
! whatever value it held at the end of the previous call to this
! subroutine. core_boundary_fx2 (FX2) is computed just above but is
! NOT what is used here -- this looks like a bug (FX2 vs FX typo) in
! the original wrtout.f, preserved exactly, not fixed.
       core_boundary_fx2 = (star%diag%del_grad(i_grad_ad,star%core_cz_top_index+1)-star%diag%del_grad(i_grad_rad,star%core_cz_top_index))/ &
             (star%diag%del_grad(i_grad_ad,star%core_cz_top_index+1)-star%diag%del_grad(i_grad_rad,star%core_cz_top_index))
       core_boundary_log_radius = star%logR(star%core_cz_top_index)+envelope_boundary_fx* &
            (star%logR(star%core_cz_top_index+1)-star%logR(star%core_cz_top_index))-log10_solar_radius
       core_boundary_radius = dexp(ln10*core_boundary_log_radius)
      else
       core_boundary_radius = 0.0D0
      endif
! JVS end

! MHP 02/12 MOVED ABOVE SECTION WHERE THESE ARE USED
!  DETERMINE CENTRAL T,P, AND DENSITY USING THE FIRST SHELL VALUES.
!  CENTRAL ETA AND BETA ARE ALSO CALCULATED.
!  EXTRAPOLATE FROM INNER SHELL P AND T TO CENTRAL P AND T
      temp_value =0.5D0*dexp(ln10*(cc13*(c4pi3l+star%logRho(1)-star%log_mass(1))+star%logRho(1)+cgl+star%log_mass(1)))
      pressure_linear = dexp(ln10*star%logP(1))
      log_pressure_center = dlog10(pressure_linear + temp_value)
!  SDEL(2,1) IS THE ACTUAL T GRADIENT AT POINT 1( = DEL)
      log_temperature_center = star%logT(1) + dlog10(1.0D0+ temp_value*star%diag%del_grad(i_grad_actual,1)/pressure_linear)
      log_density_center = star%logRho(1)
      hydrogen_fraction_center = star%xa(i_h1,1)
      metal_fraction_center = star%xa(i_metals,1)
      is_atmosphere_point = .true.
      compute_derivatives = .false.
!  CALL EQSTAT TO GET TRUE CENTRAL DENSITY, BETA, AND ETA.
      call eos_get(log_temperature_center,temperature_linear_center,log_pressure_center,pressure_linear, &
           log_density_center,density_linear_center,hydrogen_fraction_center,metal_fraction_center, &
           beta_center,beta_inverse_center,beta14_center,fxion,mean_molecular_weight_center, &
           amu_center,electron_mean_molecular_weight_center,degeneracy_eta_center,qdt_center,qdp_center, &
           qcp_center,dela_center,qdtt_center,qdtp_center,qat_center,qap_center,qcpt_center,qcpp_center, &
           compute_derivatives,is_atmosphere_point,ksaha_center, &
           composition_at_zone=star%xa(:,1), ierr=ierr)
      if (ierr /= 0) return
! MHP 02/12 MOVED ABOVE TO WHERE FIRST USED
! STORE CENTRAL RHO,P,T FOR LATER USE
      star%run%central_log10_pressure = log_pressure_center
      star%run%central_log10_temperature = log_temperature_center
      star%run%central_log10_density = log_density_center
      star%run%central_beta = beta_center
      star%run%central_degeneracy_eta = degeneracy_eta_center
! MHP 02/12 FIXED MINOR GLITCH ON BASE OF THE CONVECTION ZONE
! PROPERTIES FOR FULLY CONVECTIVE STARS; TCENTER PCENTER RHOCENTER
! WERE BEING DEFINED AFTER THIS CODE SECTION
      if(star%envelope_cz_bottom_index.lt.star%nz) then
       if(star%envelope_cz_bottom_index.gt.1) then
!  FIND MASS FRACTION OF THE ZONE EDGE AT BASE OF SURFACE C.Z.
! JVS 10/11/13 SDEL(1,JENV) IN DENOMINATOR WAS A TYPO. CHANGED TO SDEL(3,JENV)
!            FX = (SDEL(3,JENV)-SDEL(1,JENV-1))/
!     *           (SDEL(3,JENV)-SDEL(1,JENV-1))
            dd2 = star%diag%del_grad(i_grad_rad,star%envelope_cz_bottom_index-1)-star%diag%del_grad(i_grad_ad,star%envelope_cz_bottom_index-1)
            dd1 = star%diag%del_grad(i_grad_rad,star%envelope_cz_bottom_index)-star%diag%del_grad(i_grad_ad,star%envelope_cz_bottom_index)
            envelope_boundary_fx = dd2/(dd2-dd1)
!            HSB = 0.5D0*(HS1(JENV)+HS1(JENV-1))
            cz_base_mass = star%m(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%m(star%envelope_cz_bottom_index)-star%m(star%envelope_cz_bottom_index-1))
            star%run%envelope_mass = (exp(ln10*star%log_total_mass) - cz_base_mass)/solar_mass_cgs
! MHP 2/98 FIND RADIUS OF CZ BASE
            star%run%envelope_cz_log_radius = star%logR(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%logR(star%envelope_cz_bottom_index)-star%logR(star%envelope_cz_bottom_index-1))-log10_solar_radius
            star%run%envelope_radius = exp(ln10*star%run%envelope_cz_log_radius)
            star%run%envelope_cz_o16 = star%diag%so(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%diag%so(star%envelope_cz_bottom_index)-star%diag%so(star%envelope_cz_bottom_index-1))
            envelope_cz_log_temperature = star%logT(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%logT(star%envelope_cz_bottom_index)-star%logT(star%envelope_cz_bottom_index-1))
            envelope_cz_log_density = star%logRho(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%logRho(star%envelope_cz_bottom_index)-star%logRho(star%envelope_cz_bottom_index-1))
            envelope_cz_log_pressure = star%logP(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%logP(star%envelope_cz_bottom_index)-star%logP(star%envelope_cz_bottom_index-1))
            star%run%envelope_cz_temperature = exp(ln10*envelope_cz_log_temperature)
            star%run%envelope_cz_density = exp(ln10*envelope_cz_log_density)
            star%run%envelope_cz_pressure = exp(ln10*envelope_cz_log_pressure)
       else
          star%run%envelope_mass = star%star_mass
          star%run%envelope_radius = 0.0D0
            star%run%envelope_cz_temperature = 10.0D0**star%run%central_log10_temperature
            star%run%envelope_cz_density = 10.0D0**star%run%central_log10_density
            star%run%envelope_cz_pressure = 10.0D0**star%run%central_log10_pressure
            star%run%envelope_cz_o16 = star%diag%so(1)
       endif
      else
       star%run%envelope_mass = 0.0D0
       star%run%envelope_radius = 0.0D0
         star%run%envelope_cz_temperature = 0.0D0
         star%run%envelope_cz_density = 0.0D0
         star%run%envelope_cz_pressure = 0.0D0
         star%run%envelope_cz_o16 = 0.0D0
      endif



! ---- 2026 MESA-style output: fill the per-model history sources ----
! (star%run% members read by write_history). Formulas are the legacy
! .track v0 branch's, verbatim. Legacy mode skips this: wrtout
! computes the same values internally with pinned print ordering.
      if (.not. use_legacy_output) then
      call gettau(star%xa, star%logR, star%logP, star%logRho, &
           star%m, star%logT, star%fp_rot, star%ft_rot, &
           star%log_Teff, star%log_total_mass, star%log_L, star%nz, &
           star%convective_flag, star%run%envelope_radius)
      star%turnover%convective_turnover_timescale_old = &
           star%turnover%convective_turnover_timescale
      star%turnover%pphot0 = star%turnover%pphot

      star%run%log_R_surface = 0.5d0*(star%log_L + log10_solar_luminosity &
           - c4pil - csigl - 4.0d0*star%log_Teff)
      star%run%log_g_surface = cgl + star%env_comp%stotal &
           - 2.0d0*star%run%log_R_surface
      star%run%log_R_surface = star%run%log_R_surface - log10_solar_radius

      star%run%total_moment_of_inertia = 0.0d0
      if (.not. rotation_active) then
         do i = 1, star%nz
            star%run%total_moment_of_inertia = &
                 star%run%total_moment_of_inertia + &
                 cc23*star%dm(i)*exp(2.0d0*ln10*star%logR(i))
         end do
      else
         do i = 1, star%nz
            star%run%total_moment_of_inertia = &
                 star%run%total_moment_of_inertia + star%i_rot(i)
         end do
      end if

      star%flux%cl37_snu_rate = 0.0d0
      star%flux%ga71_snu_rate = 0.0d0
      if (lsnu) then
         do i = 1, 8
            star%flux%cl37_snu_rate = star%flux%cl37_snu_rate + &
                 clsnuf_diag(i)*star%flux%neutrino_flux_total(i)
            star%flux%ga71_snu_rate = star%flux%ga71_snu_rate + &
                 gasnuf_diag(i)*star%flux%neutrino_flux_total(i)
         end do
      else
         do i = 1, 10
            star%flux%neutrino_flux_total(i) = 0.0d0
         end do
      end if

      star%run%cz_moment_of_inertia = 0.0d0
      if (rotation_active) then
         star%run%rotation_period_days = min(9999.0d0, &
              0.5d0*c4pi/star%omega(star%nz)/8.64d4)
         star%run%surf_velocity_kms = star%omega(star%nz)* &
              exp(ln10*(star%run%log_R_surface+log10_solar_radius))*1.0d-5
         if (star%envelope_cz_bottom_index .lt. star%nz) then
            do i = star%envelope_cz_bottom_index, star%nz
               star%run%cz_moment_of_inertia = &
                    star%run%cz_moment_of_inertia + star%i_rot(i)
            end do
         end if
      else
         star%run%rotation_period_days = 0.0d0
         star%run%surf_velocity_kms = 0.0d0
         if (star%envelope_cz_bottom_index .lt. star%nz) then
            do i = star%envelope_cz_bottom_index, star%nz
               star%run%cz_moment_of_inertia = &
                    star%run%cz_moment_of_inertia + &
                    cc23*star%dm(i)*exp(2.0d0*ln10*star%logR(i))
            end do
         end if
      end if

      if (star%evo%has_h_shell) then
         star%run%h_shell_bot_mass = &
              star%m(star%evo%h_shell_zone_begin)/solar_mass_cgs
         star%run%h_shell_mid_mass = &
              star%m(star%evo%h_shell_midpoint_zone)/solar_mass_cgs
         star%run%h_shell_top_mass = &
              star%m(star%evo%h_shell_end_index)/solar_mass_cgs
         star%run%h_shell_bot_radius = exp(ln10*( &
              star%logR(star%evo%h_shell_zone_begin) &
              - star%run%log_R_surface - log10_solar_radius))
         star%run%h_shell_mid_radius = exp(ln10*( &
              star%logR(star%evo%h_shell_midpoint_zone) &
              - star%run%log_R_surface - log10_solar_radius))
         star%run%h_shell_top_radius = exp(ln10*( &
              star%logR(star%evo%h_shell_end_index) &
              - star%run%log_R_surface - log10_solar_radius))
      else
         star%run%h_shell_bot_mass = 0.0d0
         star%run%h_shell_mid_mass = 0.0d0
         star%run%h_shell_top_mass = 0.0d0
         star%run%h_shell_bot_radius = 0.0d0
         star%run%h_shell_mid_radius = 0.0d0
         star%run%h_shell_top_radius = 0.0d0
      end if
      end if

end subroutine update_output_diagnostics
