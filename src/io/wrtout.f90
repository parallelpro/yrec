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
subroutine wrtout(timestep_yr, log_gravity, h_shell_present_flag, &
     h_shell_begin_index, h_shell_mid_index, h_shell_end_index, &
     trial_sign_flag, punch_pending_flag, total_angular_momentum, &
     total_rotational_kinetic_energy)
      use star_info_lib, only: star, i_h1, i_he4, i_lum_grav, i_lum_he_c, i_lum_neu, i_metals, i_nu_b8, i_nu_be7, i_nu_f17, i_nu_hep, i_nu_n13, i_nu_o15, i_nu_pep, i_nu_pp, i_o16
      use luout_lib
      use const_lib
      use eos_lib
      implicit none
      integer, parameter :: json = 5000

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

      integer :: icheck
      data icheck/0/
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
! build star%run%central_log10_pressure/star%run%central_log10_temperature above) and from
! star%run%central_log10_density (DL, the input log-density estimate) -- they are
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
      if(.not.helium_flash_active) then
       write(short_file_unit,30)star%model_number,star%star_mass,star%env_comp%xnew,star%env_comp%znew,star%run%dage,timestep_yr
   30    format(1X,'MODEL NO.',I5,2X,'MASS',F13.7,2X,'(X,Z)=(',F11.9, &
          ',',F11.9,')',2X,'AGE(GYRS)',F14.8,' STEP(YRS)=',F12.0)
      else
       write(short_file_unit,40)star%model_number,star%star_mass,star%env_comp%xnew,star%env_comp%znew,star%run%dage,timestep_yr
   40    format(1X,'MODEL NO.',I5,2X,'MASS',F13.7,2X,'(X,Z)=(',F11.9, &
          ',',F11.9,')',2X,'AGE(GYRS)',F14.8,' STEP(YRS)=',1PE12.4)
      endif

      bolometric_magnitude = solar_bolometric_magnitude-2.5D0*star%log_L
      radius_log_surface = 0.5D0*(star%log_L + log10_solar_luminosity - c4pil - csigl - 4.0D0*star%log_Teff)
      log_gravity = cgl + star%env_comp%stotal - radius_log_surface - radius_log_surface
      write(short_file_unit,50)star%nz,initial_envelope_x,initial_envelope_z,star%run%core_cz_mass,star%run%envelope_mass, star%run%envelope_radius
   50 format(1X,'SHELLS=',I5,2X,'(X0,Z0)=(',F9.7,',',F9.7,')',2X, &
       'CONV. ZONE MASSES(MSUN): CORE',F10.7,' ENV.',F10.7, &
       ' RAD. FRAC.',F10.7)
      radius_log_surface = radius_log_surface - log10_solar_radius
      write(short_file_unit,60)star%log_Teff,bolometric_magnitude,star%log_L,radius_log_surface,log_gravity
   60 format(1X,'LOG(TEFF)=',F11.8,'  M(BOL)=',F11.7,'  LOG(L/LSUN)=' &
       ,F12.8,'  LOG(R/RSUN)=',F12.8,'  LOG(G) =',F12.8)
! MHP 02/12 MOVED ABOVE SECTION WHERE THESE ARE USED
!  DETERMINE CENTRAL T,P, AND DENSITY USING THE FIRST SHELL VALUES.
!  CENTRAL ETA AND BETA ARE ALSO CALCULATED.
!  EXTRAPOLATE FROM INNER SHELL P AND T TO CENTRAL P AND T
!      TEMP =0.5D0*DEXP(CLN*(CC13*(C4PI3L+HD(1)-HS(1))+HD(1)+CGL+HS(1)))
!      P = DEXP(CLN*HP(1))
!      PL = DLOG10(P + TEMP)
!  SDEL(2,1) IS THE ACTUAL T GRADIENT AT POINT 1( = DEL)
!      TL = HT(1) + DLOG10(1.0D0+ TEMP*SDEL(2,1)/P)
!      DL = HD(1)
!      X = HCOMP(1,1)
!      Z = HCOMP(3,1)
!      LATMO = .TRUE.
!      LDERIV = .FALSE.
!  CALL EQSTAT TO GET TRUE CENTRAL DENSITY, BETA, AND ETA.
! YC  If LMHD then use MHD equation of state.
!      IF (LMHD) THEN
!         CALL MEQOS(TL,T,PL,P,DL,D,X,Z,BETA,BETAI,BETA14,FXION,RMU,
!     *   AMU,EMU,ETA,QDT,QDP,QCP,DELA,QDTT,QDTP,QAT,QAP,QCPT,QCPP,
!     *   LDERIV,LATMO,KSAHA)
!      ELSE
!      IF (LDH) THEN
!         XXDH = HCOMP(1,1)
!         YYDH = HCOMP(2,1)+HCOMP(4,1)
!         ZZDH = HCOMP(3,1)
!         ZDH(1) = HCOMP(5,1)+HCOMP(6,1)
!         ZDH(2) = HCOMP(7,1)+HCOMP(8,1)
!         ZDH(3) = HCOMP(9,1)+HCOMP(10,1)+HCOMP(11,1)
!      END IF
!      CALL EQSTAT(TL,T,PL,P,DL,D,X,Z,BETA,BETAI,BETA14,FXION,RMU,
!     *AMU,EMU,ETA,QDT,QDP,QCP,DELA,QDTT,QDTP,QAT,QAP,QCPT,QCPP,LDERIV,
!     *LATMO,KSAHA)
!      END IF
      write(short_file_unit,70)star%run%central_log10_pressure,star%run%central_log10_temperature,star%run%central_log10_density,star%run%central_beta, &
           star%run%central_degeneracy_eta,star%xa(i_h1,1),star%xa(i_metals,1),star%xa(i_o16,1)
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
       fit_point_mass = star%m(h_shell_mid_index)/solar_mass_cgs
       h_shell_total_mass = (star%m(h_shell_end_index) - star%m(h_shell_begin_index-1))/solar_mass_cgs
       he_core_mass = star%m(h_shell_begin_index-1)/solar_mass_cgs
       max_log_temperature = star%run%central_log10_temperature
! LOCATE MAXIMUM T - NOTE DIFFERENT METHOD USED FOR HE FLASH
       if(.not.helium_flash_active) then
          do i = 2,star%nz
             if(star%logT(i).lt.star%logT(i-1))exit
          end do
          if (i > (star%nz)) then
          i = star%nz + 1
          end if
            max_temp_index = i - 1
          if(max_temp_index.gt.1) then
             h_shell_mid_mass = star%m(max_temp_index)/solar_mass_cgs
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
             h_shell_mid_mass = star%m(max_temp_index)/solar_mass_cgs
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
      write(short_file_unit,160) (star%flux%neutrino_flux_total(i),i=1,8)
  160 format(1X,'NEUTRINOS 1E10ERG/CM^2 PP,PEP,HEP,BE7,', &
         'B8,N13,O15,F17:', 1P8E9.2)
! DBG 7/93 from Bahcall's book p 207 table 8.2
      fl7li = 0.0D0*star%flux%neutrino_flux_total(i_nu_pp)+665.0D0*star%flux%neutrino_flux_total(i_nu_pep)+8.4D4*star%flux%neutrino_flux_total(i_nu_hep)+ &
              9.6D0*star%flux%neutrino_flux_total(i_nu_be7)+3.9D4*star%flux%neutrino_flux_total(i_nu_b8)+42.4D0*star%flux%neutrino_flux_total(i_nu_n13)+ &
              246.0D0*star%flux%neutrino_flux_total(i_nu_o15)+249.0D0*star%flux%neutrino_flux_total(i_nu_f17)
      fl37cl = 0.0D0*star%flux%neutrino_flux_total(i_nu_pp)+16.0D0*star%flux%neutrino_flux_total(i_nu_pep)+4.26D4*star%flux%neutrino_flux_total(i_nu_hep)+ &
              2.4D0*star%flux%neutrino_flux_total(i_nu_be7)+1.09D4*star%flux%neutrino_flux_total(i_nu_b8)+1.7D0*star%flux%neutrino_flux_total(i_nu_n13)+ &
              6.8D0*star%flux%neutrino_flux_total(i_nu_o15)+6.9D0*star%flux%neutrino_flux_total(i_nu_f17)
      fl71ga = 11.8D0*star%flux%neutrino_flux_total(i_nu_pp)+215.0D0*star%flux%neutrino_flux_total(i_nu_pep)+7.3D4*star%flux%neutrino_flux_total(i_nu_hep)+ &
              73.2D0*star%flux%neutrino_flux_total(i_nu_be7)+2.43D4*star%flux%neutrino_flux_total(i_nu_b8)+61.8D0*star%flux%neutrino_flux_total(i_nu_n13)+ &
              116.0D0*star%flux%neutrino_flux_total(i_nu_o15)+117.0D0*star%flux%neutrino_flux_total(i_nu_f17)
      fl81br = 0.0D0*star%flux%neutrino_flux_total(i_nu_pp)+75.0D0*star%flux%neutrino_flux_total(i_nu_pep)+9.0D4*star%flux%neutrino_flux_total(i_nu_hep)+ &
              18.3D0*star%flux%neutrino_flux_total(i_nu_be7)+2.7D4*star%flux%neutrino_flux_total(i_nu_b8)+14.5D0*star%flux%neutrino_flux_total(i_nu_n13)+ &
              36.7D0*star%flux%neutrino_flux_total(i_nu_o15)+37.0D0*star%flux%neutrino_flux_total(i_nu_f17)
      fl98mo = 0.0D0*star%flux%neutrino_flux_total(i_nu_pp)+0.0D0*star%flux%neutrino_flux_total(i_nu_pep)+10.0D4*star%flux%neutrino_flux_total(i_nu_hep)+ &
              0.0D0*star%flux%neutrino_flux_total(i_nu_be7)+3.0D4*star%flux%neutrino_flux_total(i_nu_b8)+0.0D0*star%flux%neutrino_flux_total(i_nu_n13)+ &
              0.0D0*star%flux%neutrino_flux_total(i_nu_o15)+0.0D0*star%flux%neutrino_flux_total(i_nu_f17)
      fl115in = 78.0D0*star%flux%neutrino_flux_total(i_nu_pp)+576.0D0*star%flux%neutrino_flux_total(i_nu_pep)+6.1D4*star%flux%neutrino_flux_total(i_nu_hep)+ &
              248.0D0*star%flux%neutrino_flux_total(i_nu_be7)+2.5D4*star%flux%neutrino_flux_total(i_nu_b8)+224.0D0*star%flux%neutrino_flux_total(i_nu_n13)+ &
              355.0D0*star%flux%neutrino_flux_total(i_nu_o15)+356.0D0*star%flux%neutrino_flux_total(i_nu_f17)
      write(short_file_unit,2160) fl7li,fl37cl,fl71ga,fl81br,fl98mo,fl115in
 2160 format(1X,'NEUTRINO ENERGIES (1.E-36ERG): 7Li=', 1PE9.2, &
       ' 37Cl=',1PE9.2,' 71Ga=',1PE9.2,' 81Br=',1PE9.2,' 98Mo=', &
       1PE9.2, ' 115In=', 1PE9.2)
      fit_point_mass = star%m(star%nz)/solar_mass_cgs
      write(short_file_unit,170)fit_point_mass,star%logP(star%nz),star%logT(star%nz),star%logR(star%nz)
  170 format(1X,'FIT-POINT    M/MSUN=',F16.12,5X,'(P,T,R) =',3F12.7)
      write(short_file_unit,20)
      if(ltrack) then
! MHP 02/12 MOVED ABOVE TO WHERE FIRST USED
! MHP 8/96
! STORE CENTRAL RHO,P,T FOR LATER USE
!         PCENTER = PL
!         TCENTER = TL
!         DCENTER = DL
!  Total moment of inertia
         total_moment_of_inertia = 0.0D0
         if(.not.rotation_active)then
            do i = 1,star%nz
               total_moment_of_inertia = total_moment_of_inertia + cc23*star%dm(i)*exp(2.0D0*ln10*star%logR(i))
            end do
         else
            do i = 1,star%nz
               total_moment_of_inertia = total_moment_of_inertia + star%i_rot(i)
            end do
         endif
! MHP 12/09 NEW OPTION TO OUTPUT TRACK INFORMATION IN ONE LINE PER MODEL FORMAT.
         if(track_file_version .eq. 0) then
! MHP 8/96 ADD LINE TO COMPUTE SNU's for Cl37 and Ga71.
            star%flux%cl37_snu_rate = 0.0D0
            star%flux%ga71_snu_rate = 0.0D0
            if(lsnu) then
               do i = 1,8
                  star%flux%cl37_snu_rate = star%flux%cl37_snu_rate + clsnuf(i)*star%flux%neutrino_flux_total(i)
                  star%flux%ga71_snu_rate = star%flux%ga71_snu_rate + gasnuf(i)*star%flux%neutrino_flux_total(i)
               end do
            else
               do i = 1,10
                  star%flux%neutrino_flux_total(i) = 0.0D0
               end do
            endif
! MHP 02/12
! ADDED SURFACE CZ MOMENT OF INERTIA CALCULATION TO NON-ROTATING
! MODELS, AND CORRECTLY ZERO OUT TERMS NOT COMPUTED IN SPHERICAL MODELS
! ROTATION I/O
            cz_moment_of_inertia = 0.0D0
            if(rotation_active) then
               rotation_period_days = min(9999.0D0,0.5D0*c4pi/star%omega(star%nz)/8.64D4)
               equatorial_velocity_kms = star%omega(star%nz)*exp(ln10*(radius_log_surface+log10_solar_radius))*1.0D-5
               if(star%envelope_cz_bottom_index.lt.star%nz)then
                  do k = star%envelope_cz_bottom_index,star%nz
                     cz_moment_of_inertia = cz_moment_of_inertia + star%i_rot(k)
                  end do
               endif
            else
               rotation_period_days = 0.0D0
               equatorial_velocity_kms = 0.0D0
               if(star%envelope_cz_bottom_index.lt.star%nz)then
                  do k = star%envelope_cz_bottom_index,star%nz
                     cz_moment_of_inertia = cz_moment_of_inertia + cc23*star%dm(k)*exp(2.0D0*ln10*star%logR(k))
                  end do
               endif
            endif
! END CHANGED SECTION
!        Get location of bottom, middle, top of B burning shell.
            if(h_shell_present_flag) then
             h_shell_begin_mass = star%m(h_shell_begin_index)/solar_mass_cgs
             h_shell_mid_mass2 = star%m(h_shell_mid_index)/solar_mass_cgs
             h_shell_end_mass = star%m(h_shell_end_index)/solar_mass_cgs
             h_shell_begin_radius = exp(ln10*(star%logR(h_shell_begin_index)-radius_log_surface-log10_solar_radius))
             h_shell_mid_radius = exp(ln10*(star%logR(h_shell_mid_index)-radius_log_surface-log10_solar_radius))
             h_shell_end_radius = exp(ln10*(star%logR(h_shell_end_index)-radius_log_surface-log10_solar_radius))
          else
             h_shell_begin_mass = 0.0D0
             h_shell_mid_mass2 = h_shell_begin_mass
             h_shell_end_mass = h_shell_begin_mass
             h_shell_begin_radius = h_shell_begin_mass
             h_shell_mid_radius = h_shell_begin_mass
             h_shell_end_radius = h_shell_begin_mass
             endif
! JVS 0712 Drop sinkline to get pressure at photosphere (updates PPHOT)
!             ABEG0 = ATMBEG
!            AMIN0 = ATMMIN
!            AMAX0 = ATMMAX
!            EBEG0 = ENVBEG
!            EMIN0 = ENVMIN
!            EMAX0 = ENVMAX
!            ATMBEG = ATMSTP
!            ATMMIN = ATMSTP
!            ATMMAX = ATMSTP
!            ENVBEG = ENVSTP
!            ENVMIN = ENVSTP
!            ENVMAX = ENVSTP
!            IDUM = 0
!            B = DEXP(CLN*BL)
!            FPL = FP(M)
!             FTL = FT(M)
!            KATM = 0
!             KENV = 0
!             KSAHA = 0
!CCCC            LPULPT=.FALSE.
!            IXX=0
!            LPRT=.FALSE.
!            LSBC0 = .FALSE.
!            X = HCOMP(1,M)
!             Z = HCOMP(3,M)
!            RLL = 0.5D0*(BL + CLSUNL - 4.0D0*TEFFL - C4PIL - CSIGL)
!            GL = CGL + HSTOT - RLL - RLL
!            PLIM = HP(M)
!! G Somers 10/14, FOR SPOTTED RUNS, FIND THE
!! PRESSURE AT THE AMBIENT TEMPERATURE ATEFFL
!            IF(JENV.EQ.M.AND.SPOTF.NE.0.0.AND.SPOTX.NE.1.0)THEN
!                   ATEFFL = TEFFL - 0.25*LOG10(SPOTF*SPOTX**4.0+1.0-SPOTF)
!            ELSE
!               ATEFFL = TEFFL
!            ENDIF
!            CALL ENVINT(B,FPL,FTL,GL,HSTOT,IXX,LPRT,LSBC0,
!     *            PLIM,RLL,ATEFFL,X,Z,DUM1,IDUM,KATM,KENV,KSAHA,
!     *            DUM2,DUM3,DUM4,LPULPT)
!! G Somers END
!
!! JVS 08/13 TRACKS CZ EVEN WHEN OUTSIDE THE ENVELOPE FITTING POINT
!      IF(JENV.GE.M .AND. .NOT. LC(M)) THEN
!            ENVR = ENVRCZ
!      ELSE IF(JENV.GE.M .AND. LC(M)) THEN
!            DO I=1,M
!                  DEL1(I) = SDEL(1,I)
!                  DEL2(I) = SDEL(3,I)
!            ENDDO
!                  LJVS = .TRUE.
!                  CALL TAUINT(HCOMP,HS2,HS1,LC,HR,HP,HD,HG,M,SVEL,DEL1,DEL2)
!                  LJVS = .FALSE.
!                  DD2 = DEL1(JENV-1)-DEL2(JENV-1)
!                  DD1 = DEL1(JENV)-DEL2(JENV)
!                  FX = DD2/(DD2-DD1)
!                  ENVCZL = HR(JENV-1)+FX*(HR(JENV)-HR(JENV-1))
!                  ENVR = EXP(CLN*ENVCZL)/CRSUN
!      ELSE
!            DO I=1,M
!                  DEL1(I) = SDEL(1,I)
!                  DEL2(I) = SDEL(3,I)
!            ENDDO
!                  LJVS = .TRUE.
!                  CALL TAUINT(HCOMP,HS2,HS1,LC,HR,HP,HD,HG,M,SVEL,DEL1,DEL2)
!                  LJVS = .FALSE.
!      ENDIF
!

! G Somers 3/17, ADDED CALL TO NEW TAUCZ AND PPHOT CALCULATION ROUTINE.

!       CALL GETTAU(HCOMP,HR,HP,HD,HG,HS1,HT,FP,FT,TEFFL,  ! KC 2025-05-31
      call gettau(star%xa,star%logR,star%logP,star%logRho,star%m,star%logT,star%fp_rot,star%ft_rot,star%log_Teff, &
                  star%log_total_mass,star%log_L,star%nz,star%convective_flag,star%run%envelope_radius)
      star%turnover%convective_turnover_timescale_old = star%turnover%convective_turnover_timescale
      star%turnover%pphot0 = star%turnover%pphot

! JVS 02/12 Added PPHOT and SMASS to the output
            write(itrack, 1499) star%model_number,star%nz,star%run%dage,star%log_L,radius_log_surface,log_gravity,star%log_Teff,star%run%core_cz_mass,star%run%envelope_mass, &
            star%run%envelope_radius,star%run%envelope_cz_temperature,star%run%envelope_cz_density,star%run%envelope_cz_pressure,star%run%envelope_cz_o16,star%run%central_log10_temperature,star%run%central_log10_density,star%run%central_log10_pressure,star%run%central_beta,star%run%central_degeneracy_eta,star%xa(i_h1,1),star%xa(i_he4,1), &
            star%xa(i_metals,1),(star%luminosity_breakdown(i),i = 1,5),star%luminosity_breakdown(i_lum_he_c),star%luminosity_breakdown(i_lum_grav),star%luminosity_breakdown(i_lum_neu), &
            star%flux%cl37_snu_rate,star%flux%ga71_snu_rate,(star%flux%neutrino_flux_total(i),i=1,10),(star%xa(i,1),i=4,11), &
            (star%xa(i,star%nz),i=4,15),(star%xa(i,star%nz),i=1,3),star%xa(i_metals,star%nz)/star%xa(i_h1,star%nz), &
            total_angular_momentum,total_rotational_kinetic_energy,total_moment_of_inertia,cz_moment_of_inertia,star%omega(star%nz),star%omega(1),rotation_period_days,equatorial_velocity_kms,star%turnover%convective_turnover_timescale, &
            h_shell_begin_mass,h_shell_mid_mass2,h_shell_end_mass,h_shell_begin_radius,h_shell_mid_radius,h_shell_end_radius,star%turnover%pphot,star%star_mass
! MHP 9/25 added more columns to cz depth to avoid overflow
!     1499       FORMAT(1X,2I8,1P7E16.8,0PF8.4,1P4E12.4,16E16.8,12E10.3,41E16.8)
! MCR 12/25 Preserve precision and 'E' for values w/ 3-digit exponents
 1499       format(1X,2I8,1P7E17.8E3,1P5E12.4,16E17.8E3,12E10.3,41E17.8E3)
         else if(track_file_version .eq.1 .or. track_file_version .eq.2) then
            write(itrack,1500)star%model_number,star%nz,star%run%dage,star%log_L,radius_log_surface,log_gravity,star%log_Teff,star%run%core_cz_mass,star%run%envelope_mass, &
                           star%run%envelope_radius,star%run%envelope_cz_temperature,star%run%envelope_cz_density,star%run%envelope_cz_pressure,star%run%envelope_cz_o16
 1500       format(1X,2I8,1P7E16.8,0PF8.4,1P4E12.4)
            write(itrack,1509)star%run%central_log10_temperature,star%run%central_log10_density,star%run%central_log10_pressure,star%run%central_beta,star%run%central_degeneracy_eta,star%xa(i_h1,1),star%xa(i_he4,1), &
                           star%xa(i_metals,1),total_moment_of_inertia
 1509       format(1X,1P9E16.8)
 1510       format(1X,1P8E16.8)
            write(itrack,1510)(star%luminosity_breakdown(i),i = 1,5), &
                           star%luminosity_breakdown(i_lum_he_c),star%luminosity_breakdown(i_lum_grav),star%luminosity_breakdown(i_lum_neu)

! MHP 8/96 ADD LINE TO COMPUTE SNU's for Cl37 and Ga71.
            star%flux%cl37_snu_rate = 0.0D0
            star%flux%ga71_snu_rate = 0.0D0
            do i = 1,8
               star%flux%cl37_snu_rate = star%flux%cl37_snu_rate + clsnuf(i)*star%flux%neutrino_flux_total(i)
               star%flux%ga71_snu_rate = star%flux%ga71_snu_rate + gasnuf(i)*star%flux%neutrino_flux_total(i)
            end do
          write(itrack, 1515) star%flux%cl37_snu_rate,star%flux%ga71_snu_rate,(star%flux%neutrino_flux_total(i),i=1,10)
 1515       format(1X,2F8.3,1P10E10.3)
            write(itrack,1510)(star%xa(i,1),i=4,11)
! ADD SURFACE X,Y,Z,Z/X.
            write(itrack,1520)(star%xa(i,star%nz),i=4,15), &
              (star%xa(i,star%nz),i=1,3),star%xa(i_metals,star%nz)/star%xa(i_h1,star%nz)
 1520       format(1X,1P8E16.8,/,1X,1P8E16.8)
! ROTATION I/O
            if(rotation_active) then
! MHP 8/25 removed limit on rotation period output
!     ROTP = MIN(9999.0D0,0.5D0*C4PI/OMEGA(M)/8.64D4)
               rotation_period_days = 0.5D0*c4pi/star%omega(star%nz)/8.64D4
               equatorial_velocity_kms = star%omega(star%nz)*exp(ln10*(radius_log_surface+log10_solar_radius))*1.0D-5
               cz_moment_of_inertia = 0.0D0
               if(star%convective_flag(star%nz))then
                  do k = star%envelope_cz_bottom_index,star%nz
                     cz_moment_of_inertia = cz_moment_of_inertia + star%i_rot(k)
                  end do
               else
                  cz_moment_of_inertia = 0.0D0
               endif
            endif
!        Get location of bottom, middle, top of B burning shell.
            if(h_shell_present_flag) then
             h_shell_begin_mass = star%m(h_shell_begin_index)/solar_mass_cgs
             h_shell_mid_mass2 = star%m(h_shell_mid_index)/solar_mass_cgs
             h_shell_end_mass = star%m(h_shell_end_index)/solar_mass_cgs
             h_shell_begin_radius = exp(ln10*(star%logR(h_shell_begin_index)-radius_log_surface-log10_solar_radius))
             h_shell_mid_radius = exp(ln10*(star%logR(h_shell_mid_index)-radius_log_surface-log10_solar_radius))
             h_shell_end_radius = exp(ln10*(star%logR(h_shell_end_index)-radius_log_surface-log10_solar_radius))
             else
             h_shell_begin_mass = 0.0D0
             h_shell_mid_mass2 = h_shell_begin_mass
             h_shell_end_mass = h_shell_begin_mass
             h_shell_begin_radius = h_shell_begin_mass
             h_shell_mid_radius = h_shell_begin_mass
             h_shell_end_radius = h_shell_begin_mass
          endif
            write(itrack,1530)h_shell_begin_mass,h_shell_mid_mass2,h_shell_end_mass,h_shell_begin_radius,h_shell_mid_radius,h_shell_end_radius
 1530       format(1X, 1P6E16.8)
         endif
!        ROTATION STUFF
! 4/09 ADDED TAUCZ TO ROTATION INFORMATION
         if (track_file_version .eq. 2) then

            write(itrack,1540)total_angular_momentum,total_rotational_kinetic_energy,total_moment_of_inertia,cz_moment_of_inertia,star%omega(star%nz), &
                           star%omega(1),rotation_period_days,equatorial_velocity_kms,star%turnover%convective_turnover_timescale
 1540       format(1X, 1P6E13.5,0P,2F11.5,1E13.5)
         end if
         if (track_file_version .eq. 3) then
!        RATIO OF GRAV TO TOTAL ENERGY
!       GROTOT = 100.0*TLUMX(7)/HL(M)
            write(itrack,1501)star%model_number,star%nz,star%run%dage,star%log_L,radius_log_surface,log_gravity,star%log_Teff,star%run%core_cz_mass,star%run%envelope_mass, &
            star%run%envelope_radius, star%env_comp%xnew
 1501       format(1X,2I8,1P5E13.5, 1P2E11.3, 0PF8.4, 1PE13.5)

         end if
      endif

! April 1992, DBG ISOCHRONE OUTPUT
      if(isochrone_output_active) then
! Write out model no., age (yr), L (erg/s), R (cm), Teff (K),
! g (cm/s**2), Ycenter, Mass He core (gm)
        age_yr = star%run%dage*1.0D9
        luminosity_erg_s = 10.0D0**star%log_L
          luminosity_erg_s = luminosity_erg_s*solar_luminosity_cgs
        radius_cm = 10.0D0**radius_log_surface
          radius_cm = radius_cm*solar_radius_cgs
        teff_k = 10.0D0**star%log_Teff
        gravity_cgs = 10.0D0**log_gravity
        ycenter_local = star%xa(i_he4,1)
        if (h_shell_present_flag) then
           he_core_mass_grams = star%m(h_shell_begin_index-1)
        else
           he_core_mass_grams = 0.0D0
        end if
          write(isochrone_file_unit,1005)star%model_number,age_yr,luminosity_erg_s,radius_cm,teff_k,gravity_cgs,ycenter_local, &
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
!  235    FORMAT(I4,F13.9,1P6E10.3)
!       CLOSE(ISCOMP)
!      ENDIF
!
! G Somers 11/14, WRITE THE LAST MODEL TO .LAST, AND IF LSTORE=T AND WE'RE ON
! A STORING TIMESTEP, WRITE THE EXTENDED INFORMATION TO LSTORE. IF NOT, GRAB
! THE PULSATION INFO IF LPULSE=TRUE.
!
!  STORE LAST CONVERGED MODEL IN LOGICAL UNIT ILAST
!  IF LSTORE = T, STORE EVERY NPUNCH MODELS IN LOGICAL UNIT ISTOR
!  IF LSTPCH = T, STORE THE LAST MODEL CALCULATED IN A RUN
      iwrite = ilast
      call wrtlst(iwrite,star%xa,star%logRho,star%luminosity_lsun,star%logP,star%logR,star%log_mass,star%logT,star%convective_flag,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure, &
           star%fit_point_temperature,star%fit_point_radius,star%envelope_fit_coeffs,trial_sign_flag,star%luminosity_breakdown,star%core_cz_top_index,star%envelope_cz_bottom_index,star%model_number,star%nz,star%star_mass,star%log_Teff,star%log_L,star%log_total_mass, &
           star%run%dage,timestep_yr,star%omega)
!
!  PRINT OUT MODEL DETAILS IF REQUESTED FOR THIS MODEL. THIS IS ALL DONE
!  IN THE SR PUTSTORE.
!
      if(lstore.and.mod(star%model_number,nprtmod).eq.0) then
       call putstore(star%xa,star%logRho,star%luminosity_lsun,star%logP,star%logR,star%log_mass,star%logT,star%convective_flag,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure, &
         star%fit_point_temperature,star%fit_point_radius,star%envelope_fit_coeffs,trial_sign_flag,star%luminosity_breakdown,star%core_cz_top_index,star%envelope_cz_bottom_index,star%model_number,star%nz,star%star_mass,star%log_Teff,star%log_L, &
         star%log_total_mass,star%run%dage,timestep_yr,star%omega,star%m,star%eta_squared,star%mean_radius,star%fp_rot,star%ft_rot,star%j_rot,star%i_rot)
       punch_pending_flag = .false.
      endif
! the call to putstore above creates the necessary pulsation output for LPULSE.
! however, in the event that the above block is not executed and pulsation
! output is desired, call wrtmod.
      if(.not.(lstore.and.mod(star%model_number,nprtmod).eq.0) .and. pulsation_output_active) then
       if(lmilne) call wrtmil(star%xa,star%logRho,star%luminosity_lsun,star%logP,star%logR,star%m,star%nz,star%model_number)
!        CALL WRTMOD(M,LSHELL,JXBEG,JXEND,JCORE,JENV,HCOMP,HS1,HD,HL,
!      *   HP,HR,HT,LC,MODEL,BL,TEFFL,OMEGA,FP,FT,ETA2,R0,HJM,HI,HS,
!      *   DAGE)  ! KC 2025-05-31
       call wrtmod(star%nz,star%envelope_cz_bottom_index,star%xa,star%m,star%logRho,star%luminosity_lsun, &
         star%logP,star%logR,star%logT,star%model_number,star%log_L,star%log_Teff,star%fp_rot,star%ft_rot,star%log_mass,star%run%dage)
      endif
! G Somers END
! new (2026): GYRE-format periodic pulsation output, independent of
! the LPULSE/pulsation_output_active mechanism above -- see
! core/parmin.f90 and io/write_gyre_pulse.f90.
      if (pulse_gyre_interval.gt.0 .and. mod(star%model_number,pulse_gyre_interval).eq.0) then
         write(legacy_gyre_suffix,'(I5.5)') star%model_number
         legacy_gyre_path = 'gyre_profile_'//legacy_gyre_suffix//'.data.GYRE'
         call write_gyre_pulse(star%nz,star%model_number,star%m,star%logRho,star%luminosity_lsun, &
              star%logP,star%logR,star%logT,star%omega, legacy_gyre_path)
      endif

! JVS 01/11 Added new track file output format, +manipulations for stitching
! together the interior and envelope pieces. Columns 68,69,70 are normalized
! acoustic depth, depth to CZ and acoustic crossing time, respectively.
        if (acoustic_depth_output) then
            if(star%envelope_cz_bottom_index.gt.1 .and. compute_acoustic_depth) then
                  call calcad(star%logR, star%run%envelope_cz_log_radius, star%nz, star%logRho, star%logP,star%logT,star%log_L, star%fp_rot, star%ft_rot, star%log_total_mass, &
!      *            LPRT, TEFFL, HCOMP, NKK, DAGE, DDAGE, JENV)  ! KC 2025-05-31
                  star%log_Teff, star%xa, star%run%dage, star%envelope_cz_bottom_index)
            else if (star%envelope_cz_bottom_index.eq.1) then
                  taucz_placeholder=0.0D0
            endif
            if (star%convective_flag(star%nz)) then
                  icheck=1
            else if (.not. star%convective_flag(star%nz)) then
                  icheck = 0
            endif

        if (ljlast_placeholder) then
         iwrite = ijlast_placeholder
         call wrtlst(iwrite,star%xa,star%logRho,star%luminosity_lsun,star%logP,star%logR,star%log_mass,star%logT,star%convective_flag,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure, &
         star%fit_point_temperature,star%fit_point_radius,star%envelope_fit_coeffs,trial_sign_flag,star%luminosity_breakdown,star%core_cz_top_index,star%envelope_cz_bottom_index,star%model_number,star%nz,star%star_mass,star%log_Teff,star%log_L,star%log_total_mass, &
         star%run%dage,timestep_yr,star%omega)
        endif

            write(itrack, 1800) star%model_number,star%nz,star%run%dage,star%log_L,radius_log_surface,log_gravity,star%log_Teff,star%run%core_cz_mass,star%run%envelope_mass, &
            star%run%envelope_radius,star%run%envelope_cz_temperature,star%run%envelope_cz_density,star%run%envelope_cz_pressure,star%run%envelope_cz_o16,star%run%central_log10_temperature,star%run%central_log10_density,star%run%central_log10_pressure,star%run%central_beta,star%run%central_degeneracy_eta,star%xa(i_h1,1),star%xa(i_he4,1), &
            star%xa(i_metals,1),(star%luminosity_breakdown(i),i = 1,5),star%luminosity_breakdown(i_lum_he_c),star%luminosity_breakdown(i_lum_grav),star%luminosity_breakdown(i_lum_neu), &
            star%flux%cl37_snu_rate,star%flux%ga71_snu_rate,(star%flux%neutrino_flux_total(i),i=1,10),(star%xa(i,1),i=4,11), &
            (star%xa(i,star%nz),i=4,15),(star%xa(i,star%nz),i=1,3),star%xa(i_metals,star%nz)/star%xa(i_h1,star%nz), &
            total_angular_momentum,taucz_placeholder,tcz_placeholder,tnorm_placeholder,tauhe_placeholder,whe_placeholder,tatmos_placeholder,equatorial_velocity_kms,star%turnover%convective_turnover_timescale, &
            h_shell_begin_mass,h_shell_mid_mass2,h_shell_end_mass,h_shell_begin_radius,h_shell_mid_radius,h_shell_end_radius, icheck
 1800      format(1X,2I8,1P7E16.8,0PF8.4,1P4E12.4,16E16.8,12E10.3, &
                  39E16.8, I8)

! JVS END
       endif

      return
end subroutine wrtout
