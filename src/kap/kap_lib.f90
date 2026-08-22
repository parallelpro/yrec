!----------------------------------------------------------------------
! kap_lib
!----------------------------------------------------------------------
! Added 2026 as part of the YREC readability refactor's phase two
! (disentangling the solver from the physics domains -- see
! GUIDELINES.md's "Physics domains still entangled with the solver").
! Unlike eos_lib.f90's eos_get, this is not a new dispatch consolidation:
! getopac (renamed kap_get) was already the single, clean, explicit-
! interface entry point every external caller used uniformly -- no
! duplicated dispatch logic existed at any of its 8 call sites. This
! rename/module-wrap is purely to give kap/ the same public-facade
! shape as eos_lib (a module named `<domain>_lib`, matching
! GUIDELINES.md's naming rule, which already anticipated `kap_lib` by
! name) and the same numerics_lib/eos_lib precedent for "a module
! hosting a real callable subroutine." kap_get's body, dispatch logic,
! and argument list are otherwise unchanged from getopac.
module kap_lib
      use opacity_table_lib
      implicit none
! the envelope metal fraction the surface-table machinery was
! initialized with (set by kap_init; physics-purity pass 2026 -- the
! kap domain no longer reads star_info)
      double precision, save :: kap_envelope_metal_fraction
contains

! Computes the opacity for a given composition (X, Z), blending
! between molecular/atmosphere tables, interior tables (OPAL/LAOL/
! Kurucz families, optionally interpolated between two Z values or a
! pure-Z table), and a conductive-opacity correction.
subroutine kap_get(log10_density, log10_temperature, hydrogen_fraction, &
     metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, &
     ion_fraction, ierr)
      use const_lib
      use luout_lib
      use opacity_table_lib
      implicit none

      double precision, intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      double precision, intent(inout) :: ion_fraction(3)
! 2026 (ROADMAP.md stage 3): OPTIONAL ierr, the transitional form of
! MESA's ierr-not-stop discipline. When the caller passes ierr, any
! error that historically killed the run deep inside a table lookup
! (its diagnostic message is still printed at the point of failure)
! is returned here instead, ierr /= 0, and no stop occurs. When the
! caller omits ierr, behavior is exactly historical: the same
! message, then a stop -- now located in this facade rather than the
! leaf. Existing call sites need no change.
      integer, intent(out), optional :: ierr







      save

! --- locals ---
      logical :: got_atmosphere_opacity, got_conductive_opacity
      double precision :: atm_opacity, atm_log10_opacity, &
           atm_dlnkap_dlnrho, atm_dlnkap_dlnt
      double precision :: atm_opacity_2, atm_log10_opacity_2, &
           atm_dlnkap_dlnrho_2, atm_dlnkap_dlnt_2
      double precision :: slope
      double precision :: purez_opacity, purez_log10_opacity, &
           purez_dlnkap_dlnrho, purez_dlnkap_dlnt
      double precision :: table_metal_fraction
      double precision :: opacity_2, log10_opacity_2, dlnkap_dlnrho_2, &
           dlnkap_dlnt_2
      double precision :: ramp_weight
      double precision :: conductive_opacity, conductive_log10_opacity, &
           conductive_dlnkap_dlnrho, conductive_dlnkap_dlnt
      double precision :: log10_conductive_opacity
      double precision :: radiative_opacity, radiative_dlnkap_dlnrho, &
           radiative_dlnkap_dlnt
      integer :: jerr

!     THIS SUBROUTINE CALCULATES THE OPACITY FOR A GIVEN X AND Z.
!     IF LDIFZ=T OR LZRAMP=T THEN INTERPOLATE BETWEEN TWO Z TABLES.
!     IN A SMALL T RANGE THE ATMOSPHERE AND INTERIOR OPACITY ARE
!     RAMPED FROM ONE TO THE OTHER.

!     GET ATMOSPHERE OPACITY

      if (present(ierr)) ierr = 0
      jerr = 0
      got_atmosphere_opacity = .false.
      if (log10_temperature.le.molecular_opacity_logt_max) then
         if (use_alex06_tables) then
            call getalex06(log10_density, log10_temperature, &
                 hydrogen_fraction, metal_fraction, atm_opacity, &
                 atm_log10_opacity, atm_dlnkap_dlnrho, atm_dlnkap_dlnt, &
                 jerr)
            if (jerr /= 0) go to 900
            got_atmosphere_opacity = .true.
         else if (use_alex95_tables) then
            call yalo3d(log10_density, log10_temperature, &
                 hydrogen_fraction, metal_fraction, atm_opacity, &
                 atm_log10_opacity, atm_dlnkap_dlnrho, atm_dlnkap_dlnt)
            got_atmosphere_opacity = .true.
         else if (use_kurucz90_tables) then
            call kurucz(log10_density, log10_temperature, atm_opacity, &
                 atm_log10_opacity, atm_dlnkap_dlnrho, atm_dlnkap_dlnt, &
                 jerr, *100)
            if (jerr /= 0) go to 900
            if (use_two_z_tables) then
               call kurucz2(log10_density, log10_temperature, &
                    atm_opacity_2, atm_log10_opacity_2, &
                    atm_dlnkap_dlnrho_2, atm_dlnkap_dlnt_2, jerr, *100)
               if (jerr /= 0) go to 900
               slope = (atm_log10_opacity - atm_log10_opacity_2) / &
                    (kurucz_table_z1 - kurucz_table_z2)
               atm_log10_opacity = atm_log10_opacity_2 + &
                    (metal_fraction - kurucz_table_z2)*slope
               atm_opacity = 10.0d0**atm_log10_opacity
               slope = (atm_dlnkap_dlnrho - atm_dlnkap_dlnrho_2) / &
                    (kurucz_table_z1 - kurucz_table_z2)
               atm_dlnkap_dlnrho = atm_dlnkap_dlnrho_2 + &
                    (metal_fraction - kurucz_table_z2)*slope
               slope = (atm_dlnkap_dlnt - atm_dlnkap_dlnt_2) / &
                    (kurucz_table_z1 - kurucz_table_z2)
               atm_dlnkap_dlnt = atm_dlnkap_dlnt_2 + &
                    (metal_fraction - kurucz_table_z2)*slope
            end if
            got_atmosphere_opacity = .true.
         end if
      end if
  100 continue

!     GET INTERIOR OPACITY IF NEEDED

      if (log10_temperature.lt.molecular_opacity_logt_min .and. &
           got_atmosphere_opacity) goto 1000

!     HELIUM BURNING REGION (HB EVOLUTION) USE PURE Z TABLE
!     mhp 7/12 Altered logic of the opacities in the He burnng
!     regime.  Switched to exclusive usage of OPAL below 50 million K
!     and switched the ramp to above Z = 0.1.
!     JCZ 211125 changing temperature limit to 7.0 to accommodate
!     semiconvection+overshoot HB models, which can reach lower core
!     temperatures
      if ((metal_fraction.gt.0.1d0) .and. (log10_temperature.gt.7.0d0)) then
         if (.not.use_pure_z_table) then
            write(short_file_unit, *) ' ERROR: Z>0.10 T > 5 X 10^7 K', &
                 ' NEED PURE Z TABLE TO CONTINUE. Z,LOG T=', &
                 metal_fraction, log10_temperature
            jerr = 1
            go to 900
         end if
         call gtpurz(log10_density, log10_temperature, purez_opacity, &
              purez_log10_opacity, purez_dlnkap_dlnrho, purez_dlnkap_dlnt, &
              jerr)
         if (jerr /= 0) go to 900
         if (use_opal95_tables) then
!           mhp 7/12 interpolate to maximum z in table
            table_metal_fraction = 0.1d0
            call getopal95(log10_density, log10_temperature, &
                 hydrogen_fraction, table_metal_fraction, opacity, &
                 log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, jerr)
            if (jerr /= 0) go to 900
         else if (use_opal92_tables) then
            call yllo3d(log10_density, log10_temperature, &
                 hydrogen_fraction, opacity, log10_opacity, &
                 dlnkap_dlnrho, dlnkap_dlnt)
            table_metal_fraction = opal_table_z1
         else if (use_laol89_tables) then
            call gtlaol(log10_density, log10_temperature, &
                 hydrogen_fraction, opacity, log10_opacity, &
                 dlnkap_dlnrho, dlnkap_dlnt, jerr)
            if (jerr /= 0) go to 900
            table_metal_fraction = laol_table_z1
         end if
         slope = (log10_opacity - purez_log10_opacity) / &
              (table_metal_fraction - 1.0d0)
         log10_opacity = purez_log10_opacity + (metal_fraction - 1.0d0)*slope
         opacity = 10.0d0**log10_opacity
         slope = (dlnkap_dlnrho - purez_dlnkap_dlnrho) / &
              (table_metal_fraction - 1.0d0)
         dlnkap_dlnrho = purez_dlnkap_dlnrho + &
              (metal_fraction - 1.0d0)*slope
         slope = (dlnkap_dlnt - purez_dlnkap_dlnt) / &
              (table_metal_fraction - 1.0d0)
         dlnkap_dlnt = purez_dlnkap_dlnt + (metal_fraction - 1.0d0)*slope

      else if ((metal_fraction.gt.0.12d0) .or. &
           ((abs(metal_fraction - kap_envelope_metal_fraction).gt. &
           metal_fraction_match_tolerance) .and. .not.use_two_z_tables &
           .and. .not.use_opal95_tables)) then
!        JCZ 211125 changed to 10^7 K in message to reflect above
!        change in logic.
         write(short_file_unit,*) ' Z>0.12 T < 10^7 K', &
              ' OUTSIDE OPAL OPACITY TABLE RANGE OR Z', &
              ' OUTSIDE SINGLE TABLE USED.Z,ZENV,LOG T=', &
              metal_fraction, kap_envelope_metal_fraction, log10_temperature
         jerr = 1
         go to 900

!     NOT HELIUM BURNING REGION (HB EVOLUTION) OR L2Z=T AND
!     Z STILL NOT TOO LARGE IN CORE (<.15) SO CAN USE
!     SECOND Z TABLE RATHER THAN PURE Z TABLE

      else if (use_opal95_tables) then
         call getopal95(log10_density, log10_temperature, &
              hydrogen_fraction, metal_fraction, opacity, log10_opacity, &
              dlnkap_dlnrho, dlnkap_dlnt, jerr)
         if (jerr /= 0) go to 900
      else if (use_opal92_tables) then
         call yllo3d(log10_density, log10_temperature, hydrogen_fraction, &
              opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
         if (use_two_z_tables) then
            call yllo3d2(log10_density, log10_temperature, &
                 hydrogen_fraction, opacity_2, log10_opacity_2, &
                 dlnkap_dlnrho_2, dlnkap_dlnt_2)
            slope = (log10_opacity - log10_opacity_2) / &
                 (opal_table_z1 - opal_table_z2)
            log10_opacity = log10_opacity_2 + &
                 (metal_fraction - opal_table_z2)*slope
            opacity = 10.0d0**log10_opacity
            slope = (dlnkap_dlnrho - dlnkap_dlnrho_2) / &
                 (opal_table_z1 - opal_table_z2)
            dlnkap_dlnrho = dlnkap_dlnrho_2 + &
                 (metal_fraction - opal_table_z2)*slope
            slope = (dlnkap_dlnt - dlnkap_dlnt_2) / &
                 (opal_table_z1 - opal_table_z2)
            dlnkap_dlnt = dlnkap_dlnt_2 + &
                 (metal_fraction - opal_table_z2)*slope
         end if
      else if (use_laol89_tables) then
         call gtlaol(log10_density, log10_temperature, hydrogen_fraction, &
              opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, jerr)
         if (jerr /= 0) go to 900
         if (use_two_z_tables) then
            call gtlaol2(log10_density, log10_temperature, &
                 hydrogen_fraction, opacity_2, log10_opacity_2, &
                 dlnkap_dlnrho_2, dlnkap_dlnt_2, jerr)
            if (jerr /= 0) go to 900
            slope = (log10_opacity - log10_opacity_2) / &
                 (laol_table_z1 - laol_table_z2)
            log10_opacity = log10_opacity_2 + &
                 (metal_fraction - laol_table_z2)*slope
            opacity = 10.0d0**log10_opacity
            slope = (dlnkap_dlnrho - dlnkap_dlnrho_2) / &
                 (laol_table_z1 - laol_table_z2)
            dlnkap_dlnrho = dlnkap_dlnrho_2 + &
                 (metal_fraction - laol_table_z2)*slope
            slope = (dlnkap_dlnt - dlnkap_dlnt_2) / &
                 (laol_table_z1 - laol_table_z2)
            dlnkap_dlnt = dlnkap_dlnt_2 + &
                 (metal_fraction - laol_table_z2)*slope
         end if
!     mhp 7/12 insert final trap - no opacity computed
!     should not be able to get here.
      else
         write(short_file_unit,*) 'NO OPACITY TABLE CHOSEN', &
              ' RUN STOPPED. X Z TL=', hydrogen_fraction, metal_fraction, &
              log10_temperature
         jerr = 1
         go to 900
      end if

 1000 continue
!     DO A RAMP BETWEEN SURFACE AND INTERIOR OPACITY

      if (got_atmosphere_opacity .and. &
           log10_temperature.le.molecular_opacity_logt_max) then
         if (log10_temperature.ge.molecular_opacity_logt_min) then
            ramp_weight = (log10_temperature - molecular_opacity_logt_min) / &
                 (molecular_opacity_logt_max - molecular_opacity_logt_min)
            opacity = ramp_weight*opacity + (1.0d0 - ramp_weight)*atm_opacity
            log10_opacity = dlog10(opacity)
            dlnkap_dlnrho = ramp_weight*dlnkap_dlnrho + &
                 (1.0d0 - ramp_weight)*atm_dlnkap_dlnrho
            dlnkap_dlnt = ramp_weight*dlnkap_dlnt + &
                 (1.0d0 - ramp_weight)*atm_dlnkap_dlnt
         else
            opacity = atm_opacity
            log10_opacity = atm_log10_opacity
            dlnkap_dlnrho = atm_dlnkap_dlnrho
            dlnkap_dlnt = atm_dlnkap_dlnt
         end if
      end if

!     DO CONDUCTIVE OPACITY CORRECTION
      if (use_conductive_opacity) then
!        Get Potekhin conductive opacity
         call condopacpint(log10_density, log10_temperature, &
              hydrogen_fraction, metal_fraction, conductive_opacity, &
              conductive_log10_opacity, conductive_dlnkap_dlnrho, &
              conductive_dlnkap_dlnt, ion_fraction, got_conductive_opacity, &
              jerr)
         if (jerr /= 0) go to 900
      else
         got_conductive_opacity = .false.
      end if

      if (.not.got_conductive_opacity) then
!        If we get here we have no Potekhin opacity, so we
!        try for Hubbard Lampe
         if (log10_temperature.lt.4.2d0) then
            return
         else if (log10_density.lt.(2.0d0*log10_temperature - 13.0d0)) then
            return
         else
!           Do Hubbard Lampe conductive opacity calculation
            log10_conductive_opacity = &
                 dlog10(1.0d0 - 0.6d0*hydrogen_fraction) - 14.6196d0 - &
                 (3.5853d0 + 0.1386d0*log10_density)*log10_density + &
                 (5.1324d0 - 0.3219d0*log10_temperature)*log10_temperature + &
                 0.3901d0*log10_density*log10_temperature
            conductive_opacity = 10.0d0**log10_conductive_opacity
            conductive_dlnkap_dlnrho = 0.3901d0*log10_temperature - &
                 0.2772d0*log10_density - 3.5853d0
            conductive_dlnkap_dlnt = 0.3901d0*log10_density - &
                 0.6438d0*log10_temperature + 5.1324d0
         end if
      end if

      radiative_opacity = opacity   ! Save the radiative opacity stuff
      radiative_dlnkap_dlnrho = dlnkap_dlnrho
      radiative_dlnkap_dlnt = dlnkap_dlnt
!     Add the opacities appropriately
      opacity = radiative_opacity*conductive_opacity / &
           (radiative_opacity + conductive_opacity)  ! e.g. 1/O = 1/OX + 1/OC
      log10_opacity = dlog10(opacity)
      dlnkap_dlnrho = (radiative_dlnkap_dlnrho + conductive_dlnkap_dlnrho - &
           (radiative_opacity*radiative_dlnkap_dlnrho + &
           conductive_opacity*conductive_dlnkap_dlnrho) / &
           (radiative_opacity + conductive_opacity))
      dlnkap_dlnt = (radiative_dlnkap_dlnt + conductive_dlnkap_dlnt - &
           (radiative_opacity*radiative_dlnkap_dlnt + &
           conductive_opacity*conductive_dlnkap_dlnt) / &
           (radiative_opacity + conductive_opacity))

      return

! 2026 (ROADMAP.md stage 3): single failure funnel. The diagnostic
! for whatever failed has already been written at the point of
! failure; here we either hand the error to a caller that asked for
! it, or preserve the historical stop.
  900 continue
      if (present(ierr)) then
         ierr = jerr
         return
      end if
      stop
end subroutine kap_get

!----------------------------------------------------------------------
! kap_init
!----------------------------------------------------------------------
! Added 2026 (phase three, ROADMAP.md stage 1): the kap domain's
! startup-time table-load lifecycle entry, following MESA's
! <mod>_init convention. Wraps setupopac.f90 (which reads whichever
! opacity tables the use_*_tables flags select and builds their
! interpolation splines); setup/setups.f90 previously called
! setupopac directly, making a de-facto-public entry of what is
! really internal loader machinery.
subroutine kap_init(envelope_hydrogen_fraction, &
     envelope_metal_fraction, laol_work_array, &
     alex06_table_path, kurucz_table_path, kurucz_table2_path, &
     laol_table_path, laol_table2_path, opal95_table_path, &
     opal92_table_path, opal92_table2_path, pure_z_table_path, &
     alex95_table_paths, ierr)

      use opacity_table_lib
      implicit none

      double precision, intent(in) :: envelope_hydrogen_fraction, &
           envelope_metal_fraction
      double precision, intent(inout) :: laol_work_array(12)
      character(len=256), intent(in) :: alex06_table_path, &
           kurucz_table_path, kurucz_table2_path, laol_table_path, &
           laol_table2_path, opal95_table_path, opal92_table_path, &
           opal92_table2_path, pure_z_table_path
      character(len=256), intent(in) :: alex95_table_paths(7)
! 2026 (ROADMAP.md stage 3): optional ierr -- same contract as
! kap_get's (see there). Table-load failures print their diagnostic
! at the point of failure, then either surface here or stop here.
      integer, intent(out), optional :: ierr

      integer :: jerr

      if (present(ierr)) ierr = 0
      kap_envelope_metal_fraction = envelope_metal_fraction
      call setupopac(envelope_hydrogen_fraction, laol_work_array, &
           alex06_table_path, kurucz_table_path, kurucz_table2_path, &
           laol_table_path, laol_table2_path, opal95_table_path, &
           opal92_table_path, opal92_table2_path, pure_z_table_path, &
           alex95_table_paths, jerr)
      if (jerr /= 0) then
         if (present(ierr)) then
            ierr = jerr
            return
         end if
         stop
      end if

      return
end subroutine kap_init

!----------------------------------------------------------------------
! kap_update_surface_tables
!----------------------------------------------------------------------
! Added 2026 (phase three, ROADMAP.md stage 1): public lifecycle entry
! for refreshing the cached surface-composition opacity-table slices
! (OPAL95/OPAL92/ALEX95 fixed-X tables) when the envelope hydrogen
! fraction changes. Wraps surfopac.f90; core/starin.f90 and
! setup/hpoint.f90 previously called surfopac directly -- a
! legitimate lifecycle operation that simply had no facade name.
subroutine kap_update_surface_tables(hydrogen_fraction)

      use opacity_table_lib
      implicit none

      double precision, intent(in) :: hydrogen_fraction

      call surfopac(hydrogen_fraction)

      return
end subroutine kap_update_surface_tables

end module kap_lib
