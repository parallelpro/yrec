!----------------------------------------------------------------------
! kap_lib
!----------------------------------------------------------------------
! The public interface of kap/, the opacities: everything code outside
! kap/ may use, here and only here (check_boundaries enforces it).
!     kap_init        install the settings (kap_settings_t, built by the
!                     star layer's opacity_settings) and load the tables
!                     they select
!     kap_get         the opacity of a zone and its log derivatives, in
!                     a result array (i_kap, i_log10_kap,
!                     i_dlnkap_dlnrho, i_dlnkap_dlnt)
!     kap_eval        the same with separate outputs (formerly getopac)
!     kap_update_surface_tables   rebuild the surface-composition tables
!                     when the envelope X changes
!     set_kap_settings   install settings without reloading tables
! kap_eval blends the molecular tables (Alexander 2006 or 1994/95, or
! Kurucz 1990) into the interior ones (OPAL95-format, OPAL92 or LAOL89)
! over molecular_opacity_logt_min..max, blends two Z tables linearly or
! moves toward the pure-Z table in helium-burning cores, and adds the
! conductive opacity harmonically (Potekhin, else Hubbard & Lampe). The
! table families are one module each (opal95_opacity, ...); kap reads
! no star.
module kap_lib
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: kap_settings, kap_settings_t, opacity_table, &
           set_kap_settings
      use opal95_opacity, only: compute_opal95_opacity, read_opal95_tables, &
           build_opal95_surface_table
      use conductive_opacity, only: compute_conductive_opacity
      use kurucz_opacity, only: compute_kurucz_opacity, read_kurucz_tables
      use opal92_opacity, only: compute_opal92_opacity, read_opal92_tables, &
           build_opal92_surface_table
      use laol_opacity, only: compute_laol_opacity, compute_pure_z_opacity, &
           read_laol_tables, build_laol_splines, read_pure_z_table, &
           build_pure_z_splines
      use alex06_opacity, only: compute_alex06_opacity, read_alex06_tables
      use alex94_opacity, only: compute_alex94_opacity, read_alex94_tables, &
           build_alex94_surface_table
      use math_lib
      implicit none
! Everything is private unless exported below (2026: kap_get is the
! single public query -- the named-index result-array form; the
! long-argument engine underneath it, kap_eval, is exported ONLY for
! the domain's own standalone test and stays outside the
! check_boundaries.py cross-domain allowlist). The i_*/num_kap_results
! index parameters carry their own public attributes.
      private
      public :: kap_get, kap_init, kap_update_surface_tables
      public :: kap_eval
! the settings record (opacity_tables) and its setter
      public :: kap_settings_t, set_kap_settings
! the envelope metal fraction the surface-table machinery was
! initialized with (set by kap_init, checked by kap_eval's
! single-Z-table range test)
      real(dp), save :: kap_envelope_metal_fraction
! result-array slots for kap_get (2026, MESA kap-results shape)
      integer, parameter, public :: &
           i_kap = 1, i_log10_kap = 2, i_dlnkap_dlnrho = 3, &
           i_dlnkap_dlnt = 4
      integer, parameter, public :: num_kap_results = 4

contains

! ---------------------------------------------------------------
! Named-index result-array form of kap_eval: the four opacity
! outputs packed into res(num_kap_results). ion_fraction stays an
! explicit inout argument (it is shared eos/kap ionization state,
! typically the caller's eos_res(i_fxion:i_fxion+2) slice).
subroutine kap_get(log10_density, log10_temperature, &
     hydrogen_fraction, metal_fraction, res, ion_fraction, ierr)
      real(dp), intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      real(dp), intent(inout) :: res(num_kap_results)
      real(dp), intent(inout) :: ion_fraction(3)
      integer, intent(out) :: ierr

      call kap_eval(log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction, res(i_kap), &
           res(i_log10_kap), res(i_dlnkap_dlnrho), &
           res(i_dlnkap_dlnt), ion_fraction, ierr)
end subroutine kap_get


! Computes the opacity for a given composition (X, Z), blending
! between molecular/atmosphere tables, interior tables (OPAL/LAOL/
! Kurucz families, optionally interpolated between two Z values or a
! pure-Z table), and a conductive-opacity correction.
subroutine kap_eval(log10_density, log10_temperature, hydrogen_fraction, &
     metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, &
     ion_fraction, ierr)
      use io_units
      real(dp), intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      real(dp), intent(inout) :: ion_fraction(3)
! 2026 (docs/ROADMAP.md stage 3): ierr-not-stop. Any error inside a table
! lookup (its diagnostic is still printed at the point of failure)
! comes back here as ierr /= 0; nothing in the kap domain stops.
      integer, intent(out) :: ierr
! --- locals ---
      logical :: got_atmosphere_opacity, got_conductive_opacity, in_table
      real(dp) :: atm_opacity, atm_log10_opacity, &
           atm_dlnkap_dlnrho, atm_dlnkap_dlnt
      real(dp) :: atm_opacity_2, atm_log10_opacity_2, &
           atm_dlnkap_dlnrho_2, atm_dlnkap_dlnt_2
      real(dp) :: purez_opacity, purez_log10_opacity, &
           purez_dlnkap_dlnrho, purez_dlnkap_dlnt
      real(dp) :: table_metal_fraction
      real(dp) :: opacity_2, log10_opacity_2, dlnkap_dlnrho_2, &
           dlnkap_dlnt_2
      real(dp) :: ramp_weight
      real(dp) :: conductive_opacity, conductive_log10_opacity, &
           conductive_dlnkap_dlnrho, conductive_dlnkap_dlnt
      real(dp) :: log10_conductive_opacity
      real(dp) :: radiative_opacity, radiative_dlnkap_dlnrho, &
           radiative_dlnkap_dlnt
      integer :: jerr

!     This subroutine calculates the opacity for a given X and Z.
!     If ldifz=T or lzramp=T then interpolate between two Z tables.
!     In a small T range the atmosphere and interior opacity are
!     ramped from one to the other.

!     Get atmosphere opacity

      ierr = 0
      jerr = 0
      got_atmosphere_opacity = .false.
      if (log10_temperature <= kap_settings%molecular_opacity_logt_max) then
         if (kap_settings%use_alex06_tables) then
            call compute_alex06_opacity(log10_density, log10_temperature, &
                 hydrogen_fraction, metal_fraction, atm_opacity, &
                 atm_log10_opacity, atm_dlnkap_dlnrho, atm_dlnkap_dlnt, &
                 jerr)
            if (jerr /= 0) then
               ierr = jerr
               return
            end if
            got_atmosphere_opacity = .true.
         else if (kap_settings%use_alex95_tables) then
            call compute_alex94_opacity(log10_density, log10_temperature, &
                 hydrogen_fraction, metal_fraction, atm_opacity, &
                 atm_log10_opacity, atm_dlnkap_dlnrho, atm_dlnkap_dlnt)
            got_atmosphere_opacity = .true.
         else if (kap_settings%use_kurucz90_tables) then
! outside a Kurucz table (in_table false) the interior opacity is used
! alone
            call compute_kurucz_opacity(opacity_table%kurucz(1), log10_density, &
                 log10_temperature, atm_opacity, atm_log10_opacity, &
                 atm_dlnkap_dlnrho, atm_dlnkap_dlnt, in_table, jerr)
            if (jerr /= 0) then
               ierr = jerr
               return
            end if
            if (in_table .and. kap_settings%use_two_z_tables) then
               call compute_kurucz_opacity(opacity_table%kurucz(2), log10_density, &
                    log10_temperature, atm_opacity_2, atm_log10_opacity_2, &
                    atm_dlnkap_dlnrho_2, atm_dlnkap_dlnt_2, in_table, jerr)
               if (jerr /= 0) then
                  ierr = jerr
                  return
               end if
            end if
            if (in_table .and. kap_settings%use_two_z_tables) then
               call blend_in_z(atm_log10_opacity, atm_log10_opacity_2, &
                    metal_fraction, kap_settings%kurucz_table_z1, kap_settings%kurucz_table_z2)
               atm_opacity = exp10(atm_log10_opacity)
               call blend_in_z(atm_dlnkap_dlnrho, atm_dlnkap_dlnrho_2, &
                    metal_fraction, kap_settings%kurucz_table_z1, kap_settings%kurucz_table_z2)
               call blend_in_z(atm_dlnkap_dlnt, atm_dlnkap_dlnt_2, &
                    metal_fraction, kap_settings%kurucz_table_z1, kap_settings%kurucz_table_z2)
            end if
            got_atmosphere_opacity = in_table
         end if
      end if

!     Get interior opacity if needed

      if (.not. (log10_temperature < kap_settings%molecular_opacity_logt_min .and. got_atmosphere_opacity)) then

!     Helium burning region (HB evolution) use pure Z table
!     mhp 7/12 Altered logic of the opacities in the He burnng
!     regime.  Switched to exclusive usage of OPAL below 50 million K
!     and switched the ramp to above Z = 0.1.
!     JCZ 211125 changing temperature limit to 7.0 to accommodate
!     semiconvection+overshoot HB models, which can reach lower core
!     temperatures
      if ((metal_fraction > 0.1d0) .and. (log10_temperature > 7.0d0)) then
         if (.not.kap_settings%use_pure_z_table) then
            write(run_log_unit, *) ' ERROR: Z>0.10 T > 5 X 10^7 K', &
                 ' NEED PURE Z TABLE TO CONTINUE. Z,LOG T=', &
                 metal_fraction, log10_temperature
            jerr = 1
            ierr = jerr
            return
         end if
         call compute_pure_z_opacity(log10_density, log10_temperature, purez_opacity, &
              purez_log10_opacity, purez_dlnkap_dlnrho, purez_dlnkap_dlnt, &
              jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
         if (kap_settings%use_opal95_tables) then
!           mhp 7/12 interpolate to maximum z in table
            table_metal_fraction = 0.1d0
            call compute_opal95_opacity(log10_density, log10_temperature, &
                 hydrogen_fraction, table_metal_fraction, opacity, &
                 log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, jerr)
            if (jerr /= 0) then
               ierr = jerr
               return
            end if
         else if (kap_settings%use_opal92_tables) then
            call compute_opal92_opacity(opacity_table%opal92(1), log10_density, log10_temperature, &
                 hydrogen_fraction, opacity, log10_opacity, &
                 dlnkap_dlnrho, dlnkap_dlnt, jerr)
            if (jerr /= 0) then
               ierr = jerr
               return
            end if
            table_metal_fraction = kap_settings%opal_table_z1
         else if (kap_settings%use_laol89_tables) then
            call compute_laol_opacity(opacity_table%laol(1), log10_density, log10_temperature, &
                 hydrogen_fraction, opacity, log10_opacity, &
                 dlnkap_dlnrho, dlnkap_dlnt, jerr)
            if (jerr /= 0) then
               ierr = jerr
               return
            end if
            table_metal_fraction = kap_settings%laol_table_z1
         end if
!        The pure-Z table plays the role of the second table at Z = 1.
         call blend_in_z(log10_opacity, purez_log10_opacity, &
              metal_fraction, table_metal_fraction, 1.0d0)
         opacity = exp10(log10_opacity)
         call blend_in_z(dlnkap_dlnrho, purez_dlnkap_dlnrho, &
              metal_fraction, table_metal_fraction, 1.0d0)
         call blend_in_z(dlnkap_dlnt, purez_dlnkap_dlnt, &
              metal_fraction, table_metal_fraction, 1.0d0)

      else if ((metal_fraction > 0.12d0) .or. &
           ((abs(metal_fraction - kap_envelope_metal_fraction) > &
           kap_settings%metal_fraction_match_tolerance) .and. .not.kap_settings%use_two_z_tables &
           .and. .not.kap_settings%use_opal95_tables)) then
!        JCZ 211125 changed to 10^7 K in message to reflect above
!        change in logic.
         write(run_log_unit,*) ' Z>0.12 T < 10^7 K', &
              ' OUTSIDE OPAL OPACITY TABLE RANGE OR Z', &
              ' OUTSIDE SINGLE TABLE USED.Z,ZENV,LOG T=', &
              metal_fraction, kap_envelope_metal_fraction, log10_temperature
         jerr = 1
         ierr = jerr
         return
!     Not helium burning region (HB evolution) or l2z=T and
!     Z still not too large in core (<=0.12 above) so can use
!     second Z table rather than pure Z table

      else if (kap_settings%use_opal95_tables) then
         call compute_opal95_opacity(log10_density, log10_temperature, &
              hydrogen_fraction, metal_fraction, opacity, log10_opacity, &
              dlnkap_dlnrho, dlnkap_dlnt, jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
      else if (kap_settings%use_opal92_tables) then
         call compute_opal92_opacity(opacity_table%opal92(1), log10_density, log10_temperature, hydrogen_fraction, &
              opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
         if (kap_settings%use_two_z_tables) then
            call compute_opal92_opacity(opacity_table%opal92(2), log10_density, log10_temperature, &
                 hydrogen_fraction, opacity_2, log10_opacity_2, &
                 dlnkap_dlnrho_2, dlnkap_dlnt_2, jerr)
            if (jerr /= 0) then
               ierr = jerr
               return
            end if
            call blend_in_z(log10_opacity, log10_opacity_2, &
                 metal_fraction, kap_settings%opal_table_z1, kap_settings%opal_table_z2)
            opacity = exp10(log10_opacity)
            call blend_in_z(dlnkap_dlnrho, dlnkap_dlnrho_2, &
                 metal_fraction, kap_settings%opal_table_z1, kap_settings%opal_table_z2)
            call blend_in_z(dlnkap_dlnt, dlnkap_dlnt_2, &
                 metal_fraction, kap_settings%opal_table_z1, kap_settings%opal_table_z2)
         end if
      else if (kap_settings%use_laol89_tables) then
         call compute_laol_opacity(opacity_table%laol(1), log10_density, log10_temperature, hydrogen_fraction, &
              opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
         if (kap_settings%use_two_z_tables) then
            call compute_laol_opacity(opacity_table%laol(2), log10_density, log10_temperature, &
                 hydrogen_fraction, opacity_2, log10_opacity_2, &
                 dlnkap_dlnrho_2, dlnkap_dlnt_2, jerr)
            if (jerr /= 0) then
               ierr = jerr
               return
            end if
            call blend_in_z(log10_opacity, log10_opacity_2, &
                 metal_fraction, kap_settings%laol_table_z1, kap_settings%laol_table_z2)
            opacity = exp10(log10_opacity)
            call blend_in_z(dlnkap_dlnrho, dlnkap_dlnrho_2, &
                 metal_fraction, kap_settings%laol_table_z1, kap_settings%laol_table_z2)
            call blend_in_z(dlnkap_dlnt, dlnkap_dlnt_2, &
                 metal_fraction, kap_settings%laol_table_z1, kap_settings%laol_table_z2)
         end if
!     mhp 7/12 insert final trap - no opacity computed
!     should not be able to get here.
      else
         write(run_log_unit,*) 'NO OPACITY TABLE CHOSEN', &
              ' RUN STOPPED. X Z TL=', hydrogen_fraction, metal_fraction, &
              log10_temperature
         jerr = 1
         ierr = jerr
         return
      end if

      end if
!     Do a ramp between surface and interior opacity

      if (got_atmosphere_opacity .and. &
           log10_temperature <= kap_settings%molecular_opacity_logt_max) then
         if (log10_temperature >= kap_settings%molecular_opacity_logt_min) then
            ramp_weight = (log10_temperature - kap_settings%molecular_opacity_logt_min) / &
                 (kap_settings%molecular_opacity_logt_max - kap_settings%molecular_opacity_logt_min)
            opacity = ramp_weight*opacity + (1.0d0 - ramp_weight)*atm_opacity
            log10_opacity = log10(opacity)
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

!     Do conductive opacity correction
      if (kap_settings%use_conductive_opacity) then
!        Get Potekhin conductive opacity
         call compute_conductive_opacity(log10_density, log10_temperature, &
              hydrogen_fraction, metal_fraction, conductive_opacity, &
              conductive_log10_opacity, conductive_dlnkap_dlnrho, &
              conductive_dlnkap_dlnt, ion_fraction, got_conductive_opacity, &
              jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
      else
         got_conductive_opacity = .false.
      end if

      if (.not.got_conductive_opacity) then
!        If we get here we have no Potekhin opacity, so we
!        try for Hubbard Lampe
         if (log10_temperature < 4.2d0) then
            return
         else if (log10_density < (2.0d0*log10_temperature - 13.0d0)) then
            return
         else
!           Do Hubbard Lampe conductive opacity calculation
            log10_conductive_opacity = &
                 log10(1.0d0 - 0.6d0*hydrogen_fraction) - 14.6196d0 - &
                 (3.5853d0 + 0.1386d0*log10_density)*log10_density + &
                 (5.1324d0 - 0.3219d0*log10_temperature)*log10_temperature + &
                 0.3901d0*log10_density*log10_temperature
            conductive_opacity = exp10(log10_conductive_opacity)
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
      log10_opacity = log10(opacity)
      dlnkap_dlnrho = (radiative_dlnkap_dlnrho + conductive_dlnkap_dlnrho - &
           (radiative_opacity*radiative_dlnkap_dlnrho + &
           conductive_opacity*conductive_dlnkap_dlnrho) / &
           (radiative_opacity + conductive_opacity))
      dlnkap_dlnt = (radiative_dlnkap_dlnt + conductive_dlnkap_dlnt - &
           (radiative_opacity*radiative_dlnkap_dlnt + &
           conductive_opacity*conductive_dlnkap_dlnt) / &
           (radiative_opacity + conductive_opacity))

end subroutine kap_eval

!----------------------------------------------------------------------
! kap_init
!----------------------------------------------------------------------
! Added 2026 (phase three, docs/ROADMAP.md stage 1): the kap domain's
! startup-time table-load lifecycle entry, following MESA's
! <mod>_init convention: installs the settings, reads the conductive
! table when it is on, and loads the tables the settings select
! (load_opacity_tables).
subroutine kap_init(settings, envelope_hydrogen_fraction, &
     envelope_metal_fraction, laol_work_array, &
     alex06_table_path, kurucz_table_path, kurucz_table2_path, &
     laol_table_path, laol_table2_path, opal95_table_path, &
     opal92_table_path, opal92_table2_path, pure_z_table_path, &
     alex95_table_paths, conductive_table_path, ierr)

      use conductive_opacity, only: read_conductive_table

! the run's opacity settings (installed here: kap reads no other
! configuration)
      type(kap_settings_t), intent(in) :: settings
      real(dp), intent(in) :: envelope_hydrogen_fraction, &
           envelope_metal_fraction
      real(dp), intent(inout) :: laol_work_array(12)
      character(len=256), intent(in) :: alex06_table_path, &
           kurucz_table_path, kurucz_table2_path, laol_table_path, &
           laol_table2_path, opal95_table_path, opal92_table_path, &
           opal92_table2_path, pure_z_table_path
      character(len=256), intent(in) :: alex95_table_paths(7)
      character(len=*), intent(in) :: conductive_table_path
! 2026 (docs/ROADMAP.md stage 3): same ierr contract as kap_eval's.
! Table-load failures print their diagnostic at the point of failure
! and surface here as ierr /= 0.
      integer, intent(out) :: ierr

      integer :: jerr

      ierr = 0
      call set_kap_settings(settings)
      kap_envelope_metal_fraction = envelope_metal_fraction
      if (settings%use_conductive_opacity) then
         call read_conductive_table(conductive_table_path, jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
      end if
      call load_opacity_tables(envelope_hydrogen_fraction, laol_work_array, &
           alex06_table_path, kurucz_table_path, kurucz_table2_path, &
           laol_table_path, laol_table2_path, opal95_table_path, &
           opal92_table_path, opal92_table2_path, pure_z_table_path, &
           alex95_table_paths, jerr)
      if (jerr /= 0) then
         ierr = jerr
         return
      end if

end subroutine kap_init

!----------------------------------------------------------------------
! kap_update_surface_tables
!----------------------------------------------------------------------
! Added 2026 (phase three, docs/ROADMAP.md stage 1): public lifecycle entry
! for refreshing the cached surface-composition opacity-table slices
! (OPAL95/OPAL92/ALEX94 fixed-X tables) when the envelope hydrogen
! fraction changes (refresh_surface_tables).
subroutine kap_update_surface_tables(hydrogen_fraction, ierr)
      real(dp), intent(in) :: hydrogen_fraction
! same ierr contract as kap_eval's: table errors come back as
! ierr /= 0.
      integer, intent(out) :: ierr

      call refresh_surface_tables(hydrogen_fraction, ierr)

end subroutine kap_update_surface_tables


! blend_in_z
! Linear interpolation in Z between a quantity evaluated on the first-Z
! table (value_z1, at table_z1) and on the second (value_z2, at
! table_z2); value_z1 is replaced by the blend at metal_fraction. Used
! for log10(opacity) and the two logarithmic derivatives after every
! two-table lookup in kap_eval. 2026 readability: extracted from four
! copies; the arithmetic is token-for-token the original
! (slope, then base + (Z - Z2)*slope).
pure subroutine blend_in_z(value_z1, value_z2, metal_fraction, table_z1, table_z2)
      real(dp), intent(inout) :: value_z1
      real(dp), intent(in) :: value_z2, metal_fraction, table_z1, table_z2
      real(dp) :: slope

      slope = (value_z1 - value_z2) / &
           (table_z1 - table_z2)
      value_z1 = value_z2 + &
           (metal_fraction - table_z2)*slope
end subroutine blend_in_z

!----------------------------------------------------------------------
! load_opacity_tables (formerly setupopac)
!----------------------------------------------------------------------
! Read the tables kap_settings selects and build what their lookups
! need: the interior tables (OPAL95-format and its surface table; OPAL92
! and its surface tables; LAOL89, whose header's Debye-Hueckel metal
! mixture goes to the eos domain; the pure-Z table), then one family of
! molecular tables (Alexander 2006, else Alexander 1994/95 and its
! surface table, else Kurucz 1990). A second Z table is read when
! use_two_z_tables is set.
subroutine load_opacity_tables(envelope_hydrogen_fraction, laol_work_array, &
     alex06_table_path, kurucz_table_path, kurucz_table2_path, &
     laol_table_path, laol_table2_path, opal95_table_path, &
     opal92_table_path, opal92_table2_path, pure_z_table_path, &
     alex95_table_paths, ierr)
      use eos_lib, only: eos_set_debye_huckel_z
      real(dp), intent(in) :: envelope_hydrogen_fraction
      real(dp), intent(inout) :: laol_work_array(12)
      character(len=256), intent(in) :: alex06_table_path, &
           kurucz_table_path, kurucz_table2_path, laol_table_path, &
           laol_table2_path, opal95_table_path, opal92_table_path, &
           opal92_table2_path, pure_z_table_path
      character(len=256), intent(in) :: alex95_table_paths(7)
      integer, intent(out) :: ierr
      real(dp) :: laol_debye_huckel_z(18)

      ierr = 0
! the interior tables
      if (kap_settings%use_opal95_tables) then
         call read_opal95_tables(opal95_table_path, ierr)
         if (ierr /= 0) return
         call build_opal95_surface_table(envelope_hydrogen_fraction)
      end if
      if (kap_settings%use_opal92_tables) then
         call read_opal92_tables(opal92_table_path, opal92_table2_path, ierr)
         if (ierr /= 0) return
         call build_opal92_surface_table(envelope_hydrogen_fraction, ierr)
         if (ierr /= 0) return
      end if
      if (kap_settings%use_laol89_tables) then
         call read_laol_tables(laol_work_array, laol_debye_huckel_z, laol_table_path, &
              laol_table2_path, ierr)
         if (ierr /= 0) return
         call eos_set_debye_huckel_z(laol_debye_huckel_z)
         call build_laol_splines
      end if
      if (kap_settings%use_pure_z_table) then
         call read_pure_z_table(pure_z_table_path, ierr)
         if (ierr /= 0) return
         call build_pure_z_splines
      end if
! the molecular tables
      if (kap_settings%use_alex06_tables) then
         call read_alex06_tables(alex06_table_path, ierr)
         if (ierr /= 0) return
      else if (kap_settings%use_alex95_tables) then
         call read_alex94_tables(alex95_table_paths, ierr)
         if (ierr /= 0) return
         call build_alex94_surface_table(envelope_hydrogen_fraction)
      else if (kap_settings%use_kurucz90_tables) then
         call read_kurucz_tables(kurucz_table_path, kurucz_table2_path, ierr)
         if (ierr /= 0) return
      end if
end subroutine load_opacity_tables

!----------------------------------------------------------------------
! refresh_surface_tables (formerly surfopac)
!----------------------------------------------------------------------
! Rebuild the fixed-composition surface tables at the envelope hydrogen
! fraction (OPAL95-format, OPAL92, Alexander 1994/95), from tables
! already read.
subroutine refresh_surface_tables(hydrogen_fraction, ierr)
      real(dp), intent(in) :: hydrogen_fraction
      integer, intent(out) :: ierr

      ierr = 0
      if (kap_settings%use_opal95_tables) then
         call build_opal95_surface_table(hydrogen_fraction)
      end if
      if (kap_settings%use_opal92_tables) then
         call build_opal92_surface_table(hydrogen_fraction, ierr)
         if (ierr /= 0) return
      end if
      if (kap_settings%use_alex95_tables) then
         call build_alex94_surface_table(hydrogen_fraction)
      end if
end subroutine refresh_surface_tables

end module kap_lib
