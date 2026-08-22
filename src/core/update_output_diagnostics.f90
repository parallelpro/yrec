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
      use star_info_lib, only: star
      use evolve_state_lib, only: output_diag_reset_pending

      use star_info_lib, only: star
      use envelope_comp_lib
      use eos_lib
      use luout_lib
      use const_lib
      implicit none


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
       total_luminosity_sum = star%luminosity_breakdown(1)+star%luminosity_breakdown(2)+ &
            star%luminosity_breakdown(3)+star%luminosity_breakdown(4)+star%luminosity_breakdown(5)+ &
            star%luminosity_breakdown(6)+star%luminosity_breakdown(7)+star%luminosity_breakdown(8)
       temp_value = star%luminosity_lsun(star%num_zones)/total_luminosity_sum
       do i = 1,8
          star%luminosity_breakdown(i) = star%luminosity_breakdown(i)*temp_value
   10    continue
       end do
      endif

!  CALCULATE MASS OF CENTRAL AND SURFACE CONVECTION ZONES
!  THESE MASSES ARE IN SOLAR UNITS
      if(star%core_cz_top_index.gt.1) then
       star%run%core_cz_mass = star%enclosed_mass(star%core_cz_top_index)/solar_mass_cgs
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
       core_boundary_fx2 = (star%diag%del_grad(3,star%core_cz_top_index+1)-star%diag%del_grad(1,star%core_cz_top_index))/ &
             (star%diag%del_grad(3,star%core_cz_top_index+1)-star%diag%del_grad(1,star%core_cz_top_index))
       core_boundary_log_radius = star%log_radius(star%core_cz_top_index)+envelope_boundary_fx* &
            (star%log_radius(star%core_cz_top_index+1)-star%log_radius(star%core_cz_top_index))-log10_solar_radius
       core_boundary_radius = dexp(ln10*core_boundary_log_radius)
      else
       core_boundary_radius = 0.0D0
      endif
! JVS end

! MHP 02/12 MOVED ABOVE SECTION WHERE THESE ARE USED
!  DETERMINE CENTRAL T,P, AND DENSITY USING THE FIRST SHELL VALUES.
!  CENTRAL ETA AND BETA ARE ALSO CALCULATED.
!  EXTRAPOLATE FROM INNER SHELL P AND T TO CENTRAL P AND T
      temp_value =0.5D0*dexp(ln10*(cc13*(c4pi3l+star%log_density(1)-star%log_mass(1))+star%log_density(1)+cgl+star%log_mass(1)))
      pressure_linear = dexp(ln10*star%log_pressure(1))
      log_pressure_center = dlog10(pressure_linear + temp_value)
!  SDEL(2,1) IS THE ACTUAL T GRADIENT AT POINT 1( = DEL)
      log_temperature_center = star%log_temperature(1) + dlog10(1.0D0+ temp_value*star%diag%del_grad(2,1)/pressure_linear)
      log_density_center = star%log_density(1)
      hydrogen_fraction_center = star%composition(1,1)
      metal_fraction_center = star%composition(3,1)
      is_atmosphere_point = .true.
      compute_derivatives = .false.
!  CALL EQSTAT TO GET TRUE CENTRAL DENSITY, BETA, AND ETA.
      call eos_get(log_temperature_center,temperature_linear_center,log_pressure_center,pressure_linear, &
           log_density_center,density_linear_center,hydrogen_fraction_center,metal_fraction_center, &
           beta_center,beta_inverse_center,beta14_center,fxion,mean_molecular_weight_center, &
           amu_center,electron_mean_molecular_weight_center,degeneracy_eta_center,qdt_center,qdp_center, &
           qcp_center,dela_center,qdtt_center,qdtp_center,qat_center,qap_center,qcpt_center,qcpp_center, &
           compute_derivatives,is_atmosphere_point,ksaha_center, &
           composition_at_zone=star%composition(:,1), ierr=ierr)
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
      if(star%envelope_cz_bottom_index.lt.star%num_zones) then
       if(star%envelope_cz_bottom_index.gt.1) then
!  FIND MASS FRACTION OF THE ZONE EDGE AT BASE OF SURFACE C.Z.
! JVS 10/11/13 SDEL(1,JENV) IN DENOMINATOR WAS A TYPO. CHANGED TO SDEL(3,JENV)
!            FX = (SDEL(3,JENV)-SDEL(1,JENV-1))/
!     *           (SDEL(3,JENV)-SDEL(1,JENV-1))
            dd2 = star%diag%del_grad(1,star%envelope_cz_bottom_index-1)-star%diag%del_grad(3,star%envelope_cz_bottom_index-1)
            dd1 = star%diag%del_grad(1,star%envelope_cz_bottom_index)-star%diag%del_grad(3,star%envelope_cz_bottom_index)
            envelope_boundary_fx = dd2/(dd2-dd1)
!            HSB = 0.5D0*(HS1(JENV)+HS1(JENV-1))
            cz_base_mass = star%enclosed_mass(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%enclosed_mass(star%envelope_cz_bottom_index)-star%enclosed_mass(star%envelope_cz_bottom_index-1))
            star%run%envelope_mass = (exp(ln10*star%log_total_mass) - cz_base_mass)/solar_mass_cgs
!           ENVLM = SMASS-HS1(JENV-1)/CMSUN
!          HSR = 0.5D0*(10.0D0**HR(JENV)+10.0D0**HR(JENV-1))
!          ENVX = HSR/(10.0D0**RL)
! MHP 2/98 FIND RADIUS OF CZ BASE
            star%run%envelope_cz_log_radius = star%log_radius(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%log_radius(star%envelope_cz_bottom_index)-star%log_radius(star%envelope_cz_bottom_index-1))-log10_solar_radius
            star%run%envelope_radius = exp(ln10*star%run%envelope_cz_log_radius)
            star%run%envelope_cz_o16 = star%diag%so(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%diag%so(star%envelope_cz_bottom_index)-star%diag%so(star%envelope_cz_bottom_index-1))
            envelope_cz_log_temperature = star%log_temperature(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%log_temperature(star%envelope_cz_bottom_index)-star%log_temperature(star%envelope_cz_bottom_index-1))
            envelope_cz_log_density = star%log_density(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%log_density(star%envelope_cz_bottom_index)-star%log_density(star%envelope_cz_bottom_index-1))
            envelope_cz_log_pressure = star%log_pressure(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%log_pressure(star%envelope_cz_bottom_index)-star%log_pressure(star%envelope_cz_bottom_index-1))
            star%run%envelope_cz_temperature = exp(ln10*envelope_cz_log_temperature)
            star%run%envelope_cz_density = exp(ln10*envelope_cz_log_density)
            star%run%envelope_cz_pressure = exp(ln10*envelope_cz_log_pressure)
       else
          star%run%envelope_mass = star%total_mass_msun
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



      return
end subroutine update_output_diagnostics
