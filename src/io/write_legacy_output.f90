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
      logical :: time_scaling_disabled
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

      time_scaling_disabled=.false.
!  WRITE HEADER FILE DESCRIBING THE GLOBAL PROPERTIES OF THE STAR
!  AND THE CENTRAL CONDITIONS TO THE SHORT OUTPUT FILE
!  THIS INFORMATION IS ALSO WRITTEN TO THE MODEL OUTPUT FILE IF
!  A DETAILED BREAKDOWN OF THE STELLAR STRUCTURE IS TO BE PRINTED
!  FOR THIS MODEL.
!
      write(short_file_unit,21)
   20 format(1X,127('*'))
   21 format(/,1X,127('*'))
      if(.not.star%ctrl%helium_flash_active) then
       write(short_file_unit,30)star%model_number,star%star_mass,star%xnew,star%znew,star%dage,timestep_yr
   30    format(1X,'MODEL NO.',I5,2X,'MASS',F13.7,2X,'(X,Z)=(',F11.9, &
          ',',F11.9,')',2X,'AGE(GYRS)',F14.8,' STEP(YRS)=',F12.0)
      else
       write(short_file_unit,40)star%model_number,star%star_mass,star%xnew,star%znew,star%dage,timestep_yr
   40    format(1X,'MODEL NO.',I5,2X,'MASS',F13.7,2X,'(X,Z)=(',F11.9, &
          ',',F11.9,')',2X,'AGE(GYRS)',F14.8,' STEP(YRS)=',1PE12.4)
      endif

      bolometric_magnitude = star%solar_bolometric_magnitude-2.5D0*star%log_L
      radius_log_surface = 0.5D0*(star%log_L + star%log10_solar_luminosity - c4pil - csigl - 4.0D0*star%log_Teff)
      log_gravity = cgl + star%stotal - radius_log_surface - radius_log_surface
      write(short_file_unit,50)star%nz,star%job%initial_envelope_x,star%job%initial_envelope_z,star%core_cz_mass,star%envelope_mass, star%envelope_radius
   50 format(1X,'SHELLS=',I5,2X,'(X0,Z0)=(',F9.7,',',F9.7,')',2X, &
       'CONV. ZONE MASSES(MSUN): CORE',F10.7,' ENV.',F10.7, &
       ' RAD. FRAC.',F10.7)
      radius_log_surface = radius_log_surface - star%log10_solar_radius
      write(short_file_unit,60)star%log_Teff,bolometric_magnitude,star%log_L,radius_log_surface,log_gravity
   60 format(1X,'LOG(TEFF)=',F11.8,'  M(BOL)=',F11.7,'  LOG(L/LSUN)=' &
       ,F12.8,'  LOG(R/RSUN)=',F12.8,'  LOG(G) =',F12.8)
! MHP 02/12 MOVED ABOVE SECTION WHERE THESE ARE USED
!  DETERMINE CENTRAL T,P, AND DENSITY USING THE FIRST SHELL VALUES.
!  CENTRAL ETA AND BETA ARE ALSO CALCULATED.
!  EXTRAPOLATE FROM INNER SHELL P AND T TO CENTRAL P AND T
!  SDEL(2,1) IS THE ACTUAL T GRADIENT AT POINT 1( = DEL)
!  CALL EQSTAT TO GET TRUE CENTRAL DENSITY, BETA, AND ETA.
! YC  If LMHD then use MHD equation of state.
!      IF (LMHD) THEN
!         CALL MEQOS(TL,T,PL,P,DL,D,X,Z,BETA,BETAI,BETA14,FXION,RMU,
!     *   AMU,EMU,ETA,QDT,QDP,QCP,DELA,QDTT,QDTP,QAT,QAP,QCPT,QCPP,
!     *   LDERIV,LATMO,KSAHA)
!     *AMU,EMU,ETA,QDT,QDP,QCP,DELA,QDTT,QDTP,QAT,QAP,QCPT,QCPP,LDERIV,
!     *LATMO,KSAHA)
!      END IF
      write(short_file_unit,70)star%central_log10_pressure,star%central_log10_temperature,star%central_log10_density,star%central_beta, &
           star%central_degeneracy_eta,star%xa(i_h1,1),star%xa(i_metals,1),star%xa(i_o16,1)
   70 format(1X,'CENTER: LOG P=',F10.7,' LOG T=',F10.8,' LOG D=', &
       F10.6,' BETA=',F9.7,' ETA=',0PF10.5,'  X=',0PF9.7,' Z=',F9.7, &
       ' O16=',F9.7)
      write(short_file_unit,80)(star%luminosity_breakdown(i),i = 1,5),star%luminosity_breakdown(i_lum_he_c),star%luminosity_breakdown(i_lum_neu),star%luminosity_breakdown(i_lum_grav)
   80 format(1X,'ENERGY: PPI',1PE13.6,'  PPII',E13.6,'  PPIII',E13.6, &
       '  CNO',E13.6,/,9X,'TRIPLE ALPHA',E13.6,'  HE-C',E13.6, &
       '  NEUTRINOS',E13.6,'  GRAV',E13.6)
      h_shell_mid_mass = 0.0D0
      h_shell_total_mass = 0.0D0
      if(h_shell_present_flag) then
! H-SHELL VALUES PRINTED OUT - MASSES IN SOLAR UNITS
! SS1 - MASS INTERIOR TO CENTER OF H SHELL; SS2 = MASS OF H SHELL;
! SS3 = HE CORE MASS; SS4 = MASS INTERIOR TO SHELL WITH MAXIMUM T
       fit_point_mass = star%m(h_shell_mid_index)/star%solar_mass_cgs
       h_shell_total_mass = (star%m(h_shell_end_index) - star%m(h_shell_begin_index-1))/star%solar_mass_cgs
       he_core_mass = star%m(h_shell_begin_index-1)/star%solar_mass_cgs
       max_log_temperature = star%central_log10_temperature
! LOCATE MAXIMUM T - NOTE DIFFERENT METHOD USED FOR HE FLASH
       if(.not.star%ctrl%helium_flash_active) then
          do i = 2,star%nz
             if(star%logT(i).lt.star%logT(i-1))exit
          end do
          if (i > (star%nz)) then
          i = star%nz + 1
          end if
            max_temp_index = i - 1
          if(max_temp_index.gt.1) then
             h_shell_mid_mass = star%m(max_temp_index)/star%solar_mass_cgs
             max_log_temperature = star%logT(max_temp_index)
          else
             h_shell_mid_mass = 0.0D0
             max_log_temperature = star%logT(1)
          endif
          write(short_file_unit,120)fit_point_mass,h_shell_total_mass,he_core_mass,max_log_temperature,h_shell_mid_mass
  120       format(1X,'H-SHELL MID-PT=',F10.7,' MASS TOTAL=', &
                F10.7,2X,'HE-CORE MASS=',F10.7,1X,'MAX-T=',F10.7, &
                ' (MASS=',F9.7,')')
       else
!  HE FLASH
          do i = 2,star%nz
             if(star%logT(i).lt.star%logT(i-1) .and. star%logT(i-1).gt.7.98D0) exit
          end do
          if (i > (star%nz)) then
          i = star%nz + 1
          end if
            max_temp_index = i - 1
          if(max_temp_index.gt.1) then
             h_shell_mid_mass = star%m(max_temp_index)/star%solar_mass_cgs
             max_log_temperature = star%logT(max_temp_index)
!  ADDITIONAL OUTPUT FOR HE FLASH
             max_temp_log_radius = star%logR(max_temp_index)
             max_temp_convective_flag = star%convective_flag(max_temp_index)
             write(short_file_unit,120)fit_point_mass,h_shell_total_mass,he_core_mass,max_log_temperature,h_shell_mid_mass
             write(short_file_unit,150)max_temp_convective_flag,max_temp_log_radius
  150          format(1X,'CONVECTION = ',L1,5X,'LOG(R) MAX-T =',F8.5)
          endif
       endif
!  END H-SHELL SECTION
      endif
!     PRINT OUT NEUTRINO RATES FROM ENGEB CALCULATION
      write(short_file_unit,160) (star%neutrino_flux_total(i),i=1,8)
  160 format(1X,'NEUTRINOS 1E10ERG/CM^2 PP,PEP,HEP,BE7,', &
         'B8,N13,O15,F17:', 1P8E9.2)
! DBG 7/93 from Bahcall's book p 207 table 8.2
      fl7li = 0.0D0*star%neutrino_flux_total(i_nu_pp)+665.0D0*star%neutrino_flux_total(i_nu_pep)+8.4D4*star%neutrino_flux_total(i_nu_hep)+ &
              9.6D0*star%neutrino_flux_total(i_nu_be7)+3.9D4*star%neutrino_flux_total(i_nu_b8)+42.4D0*star%neutrino_flux_total(i_nu_n13)+ &
              246.0D0*star%neutrino_flux_total(i_nu_o15)+249.0D0*star%neutrino_flux_total(i_nu_f17)
      fl37cl = 0.0D0*star%neutrino_flux_total(i_nu_pp)+16.0D0*star%neutrino_flux_total(i_nu_pep)+4.26D4*star%neutrino_flux_total(i_nu_hep)+ &
              2.4D0*star%neutrino_flux_total(i_nu_be7)+1.09D4*star%neutrino_flux_total(i_nu_b8)+1.7D0*star%neutrino_flux_total(i_nu_n13)+ &
              6.8D0*star%neutrino_flux_total(i_nu_o15)+6.9D0*star%neutrino_flux_total(i_nu_f17)
      fl71ga = 11.8D0*star%neutrino_flux_total(i_nu_pp)+215.0D0*star%neutrino_flux_total(i_nu_pep)+7.3D4*star%neutrino_flux_total(i_nu_hep)+ &
              73.2D0*star%neutrino_flux_total(i_nu_be7)+2.43D4*star%neutrino_flux_total(i_nu_b8)+61.8D0*star%neutrino_flux_total(i_nu_n13)+ &
              116.0D0*star%neutrino_flux_total(i_nu_o15)+117.0D0*star%neutrino_flux_total(i_nu_f17)
      fl81br = 0.0D0*star%neutrino_flux_total(i_nu_pp)+75.0D0*star%neutrino_flux_total(i_nu_pep)+9.0D4*star%neutrino_flux_total(i_nu_hep)+ &
              18.3D0*star%neutrino_flux_total(i_nu_be7)+2.7D4*star%neutrino_flux_total(i_nu_b8)+14.5D0*star%neutrino_flux_total(i_nu_n13)+ &
              36.7D0*star%neutrino_flux_total(i_nu_o15)+37.0D0*star%neutrino_flux_total(i_nu_f17)
      fl98mo = 0.0D0*star%neutrino_flux_total(i_nu_pp)+0.0D0*star%neutrino_flux_total(i_nu_pep)+10.0D4*star%neutrino_flux_total(i_nu_hep)+ &
              0.0D0*star%neutrino_flux_total(i_nu_be7)+3.0D4*star%neutrino_flux_total(i_nu_b8)+0.0D0*star%neutrino_flux_total(i_nu_n13)+ &
              0.0D0*star%neutrino_flux_total(i_nu_o15)+0.0D0*star%neutrino_flux_total(i_nu_f17)
      fl115in = 78.0D0*star%neutrino_flux_total(i_nu_pp)+576.0D0*star%neutrino_flux_total(i_nu_pep)+6.1D4*star%neutrino_flux_total(i_nu_hep)+ &
              248.0D0*star%neutrino_flux_total(i_nu_be7)+2.5D4*star%neutrino_flux_total(i_nu_b8)+224.0D0*star%neutrino_flux_total(i_nu_n13)+ &
              355.0D0*star%neutrino_flux_total(i_nu_o15)+356.0D0*star%neutrino_flux_total(i_nu_f17)
      write(short_file_unit,2160) fl7li,fl37cl,fl71ga,fl81br,fl98mo,fl115in
 2160 format(1X,'NEUTRINO ENERGIES (1.E-36ERG): 7Li=', 1PE9.2, &
       ' 37Cl=',1PE9.2,' 71Ga=',1PE9.2,' 81Br=',1PE9.2,' 98Mo=', &
       1PE9.2, ' 115In=', 1PE9.2)
      fit_point_mass = star%m(star%nz)/star%solar_mass_cgs
      write(short_file_unit,170)fit_point_mass,star%logP(star%nz),star%logT(star%nz),star%logR(star%nz)
  170 format(1X,'FIT-POINT    M/MSUN=',F16.12,5X,'(P,T,R) =',3F12.7)
      write(short_file_unit,20)
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
!
!  PRINT OUT MODEL DETAILS IF REQUESTED FOR THIS MODEL. THIS IS ALL DONE
!  IN THE SR PUTSTORE.
!
      if(star%ctrl%lstore.and.mod(star%model_number,star%ctrl%nprtmod).eq.0) then
       call write_store_model(star%xa,star%logRho,star%luminosity_lsun,star%logP,star%logR,star%log_mass,star%logT,star%convective_flag,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure, &
         star%fit_point_temperature,star%fit_point_radius,star%envelope_fit_coeffs,trial_sign_flag,star%luminosity_breakdown,star%core_cz_top_index,star%envelope_cz_bottom_index,star%model_number,star%nz,star%star_mass,star%log_Teff,star%log_L, &
         star%log_total_mass,star%dage,timestep_yr,star%omega,star%m,star%eta_squared,star%mean_radius,star%fp_rot,star%ft_rot,star%j_rot,star%i_rot)
       punch_pending_flag = .false.
      endif
! the call to putstore above creates the necessary pulsation output for LPULSE.
! however, in the event that the above block is not executed and pulsation
! output is desired, call wrtmod.
      if(.not.(star%ctrl%lstore.and.mod(star%model_number,star%ctrl%nprtmod).eq.0) .and. star%job%pulsation_output_active) then
       if(star%ctrl%lmilne) call write_milne(star%xa,star%logRho,star%luminosity_lsun,star%logP,star%logR,star%m,star%nz,star%model_number)
! 2026 retire-legacy: wrtmod (.pmod interior pulse + .FULL sound-
! speed table) deleted; the stitched pulse files carry everything.
      endif
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
