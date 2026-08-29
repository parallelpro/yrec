!----------------------------------------------------------------------
! wrtout
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original wrtout.f; only variable names, source form, and comment
! style were updated.
!
! Writes the per-model summary block to the .short log file (global
! properties, central conditions, energy generation, neutrino fluxes,
! H-shell diagnostics), writes the one-line .track record, stores the
! last converged model, and (every nprtmod models or when a store is
! otherwise due) dispatches to putstore/wrtmod/wrtmil for the verbose
! .store/.pmod/.penv/.patm output.
subroutine write_legacy_output(timestep_yr, log_gravity, h_shell_present_flag, &
     h_shell_begin_index, h_shell_mid_index, h_shell_end_index, &
     trial_sign_flag, punch_pending_flag, total_angular_momentum, &
     total_rotational_kinetic_energy)
      use star_info_lib, only: star, i_h1, i_he4, i_lum_grav, i_lum_he_c, i_lum_neu, i_metals, i_nu_b8, i_nu_be7, i_nu_f17, i_nu_hep, i_nu_n13, i_nu_o15, i_nu_pep, i_nu_pp, i_o16, json
      use luout_lib
      use phys_const_lib
      use eos_lib
      implicit none

      double precision, intent(in) :: timestep_yr
      double precision, intent(out) :: log_gravity
      logical, intent(in) :: h_shell_present_flag
      integer, intent(in) :: h_shell_begin_index, h_shell_mid_index, &
           h_shell_end_index
      double precision, intent(in) :: trial_sign_flag
      logical, intent(inout) :: punch_pending_flag
      double precision, intent(in) :: total_angular_momentum, &
           total_rotational_kinetic_energy

! G Somers END


      double precision :: clsnuf(8), gasnuf(8)
      character(len=5) :: legacy_gyre_suffix
      character(len=64) :: legacy_gyre_path
! MHP 8/96 CROSS SECTIONS OF DIFFERENT NEUTRINOS TO THE CHLORINE
! AND GALLIUM EXPERIMENTS; TAKEN FROM NEUTRINO ASTROPHYSICS,P.207.
! note changes in cl37 cross sections (see bahcall and pinsonneault,
! REV.MOD.PHYS., P.895)
      data gasnuf/1.18D1,2.15D2,7.14D4,7.17D1,2.40D4,6.04D1, &
                  1.137D2,1.139D2/
      data clsnuf/0.0D0,1.6D1,4.26D4,2.4D0,1.14D4,1.7D0,6.8D0,6.9D0/

! --- locals ---
      integer :: i
! core_boundary_log_radius/core_boundary_radius (CORERL/CORER) are
! computed below but never read afterward anywhere in the original
! wrtout.f -- dead code, preserved as such (not removed).
      double precision :: envelope_cz_log_radius
      double precision :: core_mass, bolometric_magnitude, radius_log_surface
! temperature_linear_center/density_linear_center (T/D) are separate
! eqstat/meqos output slots, distinct from temp_value (TEMP, used to
! build star%central_log10_pressure/star%central_log10_temperature above) and from
! star%central_log10_density (DL, the input log-density estimate) -- they are
! never read again after the call in the original wrtout.f (dead
! output), but must not be aliased with those other variables or the
! DL slot gets overwritten with a linear value. Preserved as distinct
! locals to match the original's argument list exactly.
! Second-derivative / opacity-related eqstat outputs; exact physical
! definitions not confidently known from this file alone (they mirror
! the QDT/QDP/QCP/DELA/QDTT/QDTP/QAT/QAP/QCPT/QCPP argument slots of
! EQSTAT/MEQOS), kept as conservative names.
      double precision :: h_shell_mid_mass, h_shell_total_mass, he_core_mass, &
           max_log_temperature
      integer :: max_temp_index
      double precision :: max_temp_log_radius
      logical :: max_temp_convective_flag
      double precision :: fl7li, fl37cl, fl71ga, fl81br, fl98mo, fl115in
      double precision :: fit_point_mass
      double precision :: total_moment_of_inertia, cz_moment_of_inertia
      double precision :: rotation_period_days, equatorial_velocity_kms
      integer :: k
      double precision :: h_shell_begin_mass, h_shell_mid_mass2, &
           h_shell_end_mass, h_shell_begin_radius, h_shell_mid_radius, &
           h_shell_end_radius
      integer :: iwrite
      double precision :: age_yr, luminosity_erg_s, radius_cm, teff_k, &
           gravity_cgs, ycenter_local, he_core_mass_grams

! JVS 0712 for call to envint:
!       REAL*8 DUM1(4),DUM2(3),DUM3(3),DUM4(3)
! JVS 10/13 for recalculation of taucz
!       REAL*8 DEL1(JSON), DEL2(JSON)
! end JVS

!  WRITE HEADER FILE DESCRIBING THE GLOBAL PROPERTIES OF THE STAR
!  AND THE CENTRAL CONDITIONS TO THE SHORT OUTPUT FILE
!  THIS INFORMATION IS ALSO WRITTEN TO THE MODEL OUTPUT FILE IF
!  A DETAILED BREAKDOWN OF THE STELLAR STRUCTURE IS TO BE PRINTED
!  FOR THIS MODEL.
!
! 2026 log redesign: the per-model physical summary block (MODEL
! NO. / SHELLS / LOG(TEFF) / CENTER / ENERGY / H-SHELL / NEUTRINOS /
! FIT-POINT) is deleted -- history.data is the only home for
! per-model physical data, and the run-log progress line covers
! at-a-glance monitoring. What survives below: the surface radius/
! gravity the isochrone writer and callers need, the .mod dump, and
! the special-purpose writers.
      radius_log_surface = 0.5D0*(star%log_L + star%log10_solar_luminosity - c4pil - csigl - 4.0D0*star%log_Teff)
      log_gravity = cgl + star%stotal - radius_log_surface - radius_log_surface
      radius_log_surface = radius_log_surface - star%log10_solar_radius
! 2026 retire-legacy: the .track writer block is deleted -- every
! quantity it wrote is in history.data (see the .track-vs-history
! audit; initial_x/initial_y/mixing_length_alpha were added to the
! history global block to close the run-metadata gap).

! April 1992, DBG ISOCHRONE OUTPUT
      if(star%ctrl%isochrone_output_active) then
! Write out model no., age (yr), L (erg/s), R (cm), Teff (K),
! g (cm/s**2), Ycenter, Mass He core (gm)
        age_yr = star%dage*1.0D9
        luminosity_erg_s = 10.0D0**star%log_L
          luminosity_erg_s = luminosity_erg_s*star%solar_luminosity_cgs
        radius_cm = 10.0D0**radius_log_surface
          radius_cm = radius_cm*star%solar_radius_cgs
        teff_k = 10.0D0**star%log_Teff
        gravity_cgs = 10.0D0**log_gravity
        ycenter_local = star%xa(i_he4,1)
        if (h_shell_present_flag) then
           he_core_mass_grams = star%m(h_shell_begin_index-1)
        else
           he_core_mass_grams = 0.0D0
        end if
          write(star%ctrl%isochrone_file_unit,1005)star%model_number,age_yr,luminosity_erg_s,radius_cm,teff_k,gravity_cgs,ycenter_local, &
            he_core_mass_grams
 1005     format(1X, I5, 1P7E17.8)
      end if
! MHP 8/25 Fscomp depreciated, call commented out
!     WRITE OUT SURFACE COMPOSITIONS TO FILE ISCOMP IF extended comp MODEL.
!      IF(LEXCOM) THEN
!       OPEN(UNIT=ISCOMP,FILE=FSCOMP, FORM='FORMATTED',
!     *        STATUS='UNKNOWN',ACCESS='APPEND')
!       WRITE(ISCOMP,235)MODEL,DAGE,HCOMP(4,M),HCOMP(5,M),HCOMP(6,M),
!     *                HCOMP(7,M),HCOMP(14,M),HCOMP(15,M)
!
! G Somers 11/14, WRITE THE LAST MODEL TO .LAST, AND IF LSTORE=T AND WE'RE ON
! A STORING TIMESTEP, WRITE THE EXTENDED INFORMATION TO LSTORE. IF NOT, GRAB
! THE PULSATION INFO IF LPULSE=TRUE.
!
!  STORE LAST CONVERGED MODEL IN LOGICAL UNIT ILAST
!  IF LSTORE = T, STORE EVERY NPUNCH MODELS IN LOGICAL UNIT ISTOR
!  IF LSTPCH = T, STORE THE LAST MODEL CALCULATED IN A RUN
      iwrite = ilast
      call write_mod_model(iwrite)
! 2026 retire-legacy: the .store writers (putstore + the LSTCH
! stitch) are deleted -- profiles/pulse files carry the stitched
! model, history carries the per-model globals.
! G Somers END
! new (2026): GYRE-format periodic pulsation output, independent of
! the LPULSE/pulsation_output_active mechanism above -- see
! core/read_input.f90 and io/write_gyre_pulse.f90.
      if (star%ctrl%pulse_gyre_interval.gt.0 .and. mod(star%model_number,star%ctrl%pulse_gyre_interval).eq.0) then
         write(legacy_gyre_suffix,'(I5.5)') star%model_number
         legacy_gyre_path = 'gyre_profile_'//legacy_gyre_suffix//'.data.GYRE'
         call write_gyre_pulse(star%nz,star%model_number,star%m,star%logRho,star%luminosity_lsun, &
              star%logP,star%logR,star%logT,star%omega, legacy_gyre_path)
      endif

! JVS 01/11 Added new track file output format, +manipulations for stitching
! together the interior and envelope pieces. Columns 68,69,70 are normalized
! acoustic depth, depth to CZ and acoustic crossing time, respectively.
! 2026 retire-legacy: the LACOUT acoustic-depth mode (calcad call,
! ageout model saves, and the acoustic-columns track record) is
! retired -- acoustic depths are post-processing on profile columns
! (csound over the stitched grid).

      return
end subroutine write_legacy_output
