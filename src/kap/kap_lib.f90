!----------------------------------------------------------------------
! kap_lib
!----------------------------------------------------------------------
! The public interface of kap/, the opacities: everything code outside
! kap/ may use, here and only here (check_boundaries enforces it).
!     kap_init        check and install the settings (kap_settings_t,
!                     built by the star layer's kap_settings_from_star),
!                     choose the families they select, and load their
!                     tables from the paths record
!     kap_get         the opacity of a zone and its log derivatives, in
!                     a result array (i_kap, i_log10_kap,
!                     i_dlnkap_dlnrho, i_dlnkap_dlnt)
!     kap_update_surface_tables   rebuild the surface-composition tables
!                     when the envelope X changes
!     set_kap_settings   install settings without reloading tables (for
!                     test programs: choosing families the tables of
!                     which were not loaded gives an error at a lookup)
! kap_get evaluates the opacity in named steps (kap_eval): the molecular
! family's opacity below molecular_opacity_logt_max; the interior
! family's above molecular_opacity_logt_min (blended toward the pure-Z
! table in a helium-burning core); the ramp between the two; and the
! conductive opacity added harmonically (Potekhin's table, else Hubbard
! & Lampe's fit). Which family plays each role is opacity_families'
! precedence list; each family is one module (opal95_opacity, ...). kap
! reads no star, and holds its state (opacity_tables' types) here.
module kap_lib
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: kap_settings, kap_settings_t, kap_table_paths_t, &
           opacity_table_state, install_kap_settings
      use opacity_source, only: kap_point, kap_result, kap_source
      use opacity_families, only: select_families
      use opacity_blending, only: blend_results_in_z, ramp_to_molecular, &
           hubbard_lampe_applies, hubbard_lampe_opacity, add_conduction
      use kap_messages, only: kap_log_values, kap_report, message_length, unwritten_line
      implicit none
      private
      public :: kap_get, kap_init, kap_update_surface_tables
! the settings record (opacity_tables), its setter, and the paths record
      public :: kap_settings_t, set_kap_settings, kap_table_paths_t
! result-array slots for kap_get (MESA's kap-results shape)
      integer, parameter, public :: &
           i_kap = 1, i_log10_kap = 2, i_dlnkap_dlnrho = 3, &
           i_dlnkap_dlnt = 4
      integer, parameter, public :: num_kap_results = 4

! A helium-burning core: Z above helium_core_metal_fraction and log T
! above helium_core_logt (10^7 K) take the pure-Z table's blend. Below
! that temperature, Z above max_interior_metal_fraction is outside the
! interior tables.
      real(dp), parameter :: helium_core_metal_fraction = 0.1d0, &
           helium_core_logt = 7.0d0, max_interior_metal_fraction = 0.12d0

! kap's state: every family's tables, their warm-start cursors and
! composition caches (opacity_tables' types), and the envelope metal
! fraction the tables were loaded for. Only this module writes it: the
! family routines get their sets as arguments. initial_tables is never
! written: it holds what a fresh process starts from (the grids' and
! cursors' initial values, zero elsewhere), and kap_init starts each run
! from it, so a second run in one process starts where the first did.
      type(opacity_table_state), save :: tables, initial_tables
! the families the settings select (opacity_families), one per role
      type(kap_source), save :: interior, molecular, helium_core, conduction

contains

!----------------------------------------------------------------------
! kap_get
!----------------------------------------------------------------------
! The opacity at (log10 rho, log10 T, X, Z) and its log derivatives, in
! res (the i_* slots), from kap_eval. ion_fraction is the ionized
! fractions of H, He+ and He++ (the caller's eos_res%ion_fraction),
! which the conductive opacity weights its species by; kap does not
! change them. On an error (ierr /= 0, with its message written) res is
! left as it was.
subroutine kap_get(log10_density, log10_temperature, &
     hydrogen_fraction, metal_fraction, res, ion_fraction, ierr)
      real(dp), intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      real(dp), intent(inout) :: res(num_kap_results)
      real(dp), intent(inout) :: ion_fraction(3)
      integer, intent(out) :: ierr
      type(kap_point) :: point
      type(kap_result) :: result

      point = kap_point(log10_density, log10_temperature, hydrogen_fraction, &
           metal_fraction, ion_fraction)
      call kap_eval(point, result, ierr)
      if (ierr /= 0) return
      res(i_kap) = result%opacity
      res(i_log10_kap) = result%log10_opacity
      res(i_dlnkap_dlnrho) = result%dlnkap_dlnrho
      res(i_dlnkap_dlnt) = result%dlnkap_dlnt
end subroutine kap_get

!----------------------------------------------------------------------
! kap_eval
!----------------------------------------------------------------------
! The opacity at point, in four steps:
!   1. molecular_opacity: the molecular family's, at log T up to its
!      window's top (molecular_opacity_logt_max);
!   2. interior_opacity: the interior family's, unless the point is
!      below the interior window (molecular_opacity_logt_min) and the
!      molecular family has a value there;
!   3. the ramp: between the windows' ends, linear in kappa from the
!      molecular to the interior opacity (below it, the molecular one);
!   4. conduction: the Potekhin conductive opacity, else Hubbard &
!      Lampe's where it applies, added harmonically.
! Each step returns ierr from the family it calls, which has written the
! message.
subroutine kap_eval(point, result, ierr)
      type(kap_point), intent(in) :: point
      type(kap_result), intent(out) :: result
      integer, intent(out) :: ierr
      type(kap_result) :: molecular_result, conductive_result
      logical :: have_molecular, have_conductive

      call molecular_opacity(point, molecular_result, have_molecular, ierr)
      if (ierr /= 0) return
      if (.not. (point%log10_temperature < interior%logt_min .and. have_molecular)) then
         call interior_opacity(point, result, ierr)
         if (ierr /= 0) return
      end if
      if (have_molecular .and. point%log10_temperature <= molecular%logt_max) then
         call ramp_to_molecular(point%log10_temperature, interior%logt_min, &
              molecular%logt_max, molecular_result, result)
      end if
      call conductive_opacity(point, conductive_result, have_conductive, ierr)
      if (ierr /= 0) return
      if (.not. have_conductive) then
         if (.not. hubbard_lampe_applies(point%log10_density, point%log10_temperature)) return
         conductive_result = hubbard_lampe_opacity(point%log10_density, &
              point%log10_temperature, point%hydrogen_fraction)
      end if
      call add_conduction(result, conductive_result)
end subroutine kap_eval

! Step 1: the molecular family's opacity, where its window reaches the
! point; have_molecular is .false. without one, or where its table has
! no value (Kurucz 1990).
subroutine molecular_opacity(point, res, have_molecular, ierr)
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      logical, intent(out) :: have_molecular
      integer, intent(out) :: ierr

      ierr = 0
      have_molecular = .false.
      if (.not. (molecular%selected .and. point%log10_temperature <= molecular%logt_max)) return
      call molecular%lookup(kap_settings, tables, point, res, have_molecular, ierr)
      if (ierr /= 0) have_molecular = .false.
end subroutine molecular_opacity

! Step 2: the interior family's opacity. A helium-burning core blends it
! toward the pure-Z table; above max_interior_metal_fraction, or at a Z
! off the envelope's for a family without a table there, there is none.
subroutine interior_opacity(point, res, ierr)
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      integer, intent(out) :: ierr
      logical :: in_table

      if (point%metal_fraction > helium_core_metal_fraction .and. &
           point%log10_temperature > helium_core_logt) then
         call helium_core_opacity(point, res, ierr)
      else if ((point%metal_fraction > max_interior_metal_fraction) .or. &
           ((abs(point%metal_fraction - tables%envelope_metal_fraction) > &
           kap_settings%metal_fraction_match_tolerance) .and. .not.kap_settings%use_two_z_tables &
           .and. .not.interior%covers_any_z)) then
         call kap_log_values([point%metal_fraction, tables%envelope_metal_fraction, &
              point%log10_temperature], &
              ' Z>0.12 T < 10^7 K' // &
              ' OUTSIDE OPAL OPACITY TABLE RANGE OR Z' // &
              ' OUTSIDE SINGLE TABLE USED.Z,ZENV,LOG T=')
         ierr = 1
      else if (interior%selected) then
         call interior%lookup(kap_settings, tables, point, res, in_table, ierr)
      else
         call no_interior_family(point, ierr)
      end if
end subroutine interior_opacity

! Step 2 in a helium-burning core: the interior family's single-Z table
! (at its Z, table_z) blended linearly in Z toward the pure-Z table,
! placed at Z = 1.
subroutine helium_core_opacity(point, res, ierr)
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      integer, intent(out) :: ierr
      type(kap_result) :: pure_z_result
      real(dp) :: table_z
      logical :: in_table

      if (.not. helium_core%selected) then
         call kap_log_values([point%metal_fraction, point%log10_temperature], &
              ' ERROR: Z>0.10 T > 5 X 10^7 K' // &
              ' NEED PURE Z TABLE TO CONTINUE. Z,LOG T=')
         ierr = 1
         return
      end if
      call helium_core%lookup(kap_settings, tables, point, pure_z_result, in_table, ierr)
      if (ierr /= 0) return
      if (.not. interior%selected) then
! (kap_init refuses these settings; set_kap_settings does not check)
         call no_interior_family(point, ierr)
         return
      end if
      call interior%core_lookup(kap_settings, tables, point, res, table_z, ierr)
      if (ierr /= 0) return
      call blend_results_in_z(res, pure_z_result, point%metal_fraction, table_z, 1.0d0)
end subroutine helium_core_opacity

subroutine no_interior_family(point, ierr)
      type(kap_point), intent(in) :: point
      integer, intent(out) :: ierr

      call kap_log_values([point%hydrogen_fraction, point%metal_fraction, &
           point%log10_temperature], 'NO OPACITY TABLE CHOSEN' // ' RUN STOPPED. X Z TL=')
      ierr = 1
end subroutine no_interior_family

! Step 4's table: the Potekhin conductive opacity, where it has one
! (have_conductive).
subroutine conductive_opacity(point, res, have_conductive, ierr)
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      logical, intent(out) :: have_conductive
      integer, intent(out) :: ierr

      ierr = 0
      have_conductive = .false.
      if (.not. conduction%selected) return
      call conduction%lookup(kap_settings, tables, point, res, have_conductive, ierr)
end subroutine conductive_opacity

!----------------------------------------------------------------------
! kap_init
!----------------------------------------------------------------------
! Once per run (setup's install_run_settings, at the start of run_yrec):
! check and install the settings, start from initial_tables, choose the
! families (opacity_families) and load their tables: the conductive
! table, the interior family's, the pure-Z table and the molecular
! family's. A setting that cannot work, or a table that cannot be read,
! returns ierr = 1 with its message written (kap_messages).
subroutine kap_init(settings, paths, envelope_hydrogen_fraction, &
     envelope_metal_fraction, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(kap_table_paths_t), intent(in) :: paths
      real(dp), intent(in) :: envelope_hydrogen_fraction, &
           envelope_metal_fraction
      integer, intent(out) :: ierr

      call check_kap_settings(settings, ierr)
      if (ierr /= 0) return
      call set_kap_settings(settings)
      tables = initial_tables
      tables%envelope_metal_fraction = envelope_metal_fraction
      call load_family(conduction, ierr)
      if (ierr /= 0) return
      call load_family(interior, ierr)
      if (ierr /= 0) return
      call load_family(helium_core, ierr)
      if (ierr /= 0) return
      call load_family(molecular, ierr)

contains

subroutine load_family(family, ierr)
      type(kap_source), intent(in) :: family
      integer, intent(out) :: ierr

      ierr = 0
      if (family%selected) call family%load(kap_settings, paths, &
           envelope_hydrogen_fraction, tables, ierr)
end subroutine load_family

end subroutine kap_init

!----------------------------------------------------------------------
! set_kap_settings
!----------------------------------------------------------------------
! Install settings and choose their families, without loading tables
! (kap_init's first step; test programs call it to switch a family off).
subroutine set_kap_settings(settings)
      type(kap_settings_t), intent(in) :: settings

      call install_kap_settings(settings)
      call select_families(settings, interior, molecular, helium_core, conduction)
end subroutine set_kap_settings

!----------------------------------------------------------------------
! check_kap_settings
!----------------------------------------------------------------------
! Refuse, before any table is read, settings that cannot give an
! opacity (each would otherwise give a NaN or an undefined opacity with
! ierr = 0 at some lookup): ierr = 1 and a kap_report naming the
! controls.
!   - the pure-Z table blends toward an interior table above Z = 0.1, so
!     it needs one;
!   - with two Z tables (the Z ramp or Z diffusion), a family in use that
!     blends its two tables linearly in Z needs two different table Zs
!     (OPAL 1995 and the Alexander tables interpolate in Z instead);
!   - the molecular ramp needs molecular_opacity_logt_min below
!     molecular_opacity_logt_max when a molecular family is on.
subroutine check_kap_settings(settings, ierr)
      type(kap_settings_t), intent(in) :: settings
      integer, intent(out) :: ierr
      type(kap_source) :: interior_family, molecular_family, helium_core_family, &
           conduction_family
      character(len=message_length) :: lines(1)

      ierr = 0
      lines = unwritten_line
      call select_families(settings, interior_family, molecular_family, &
           helium_core_family, conduction_family)
      if (helium_core_family%selected .and. .not. interior_family%selected) then
         write(lines,'(a)') 'kap_init: use_pure_z_table needs an interior table ' // &
              '(use_opal95_tables, use_opal92_tables or use_laol89_tables)'
      else if (settings%use_two_z_tables .and. same_table_zs(interior_family)) then
         call two_z_message(interior_family)
      else if (settings%use_two_z_tables .and. same_table_zs(molecular_family)) then
         call two_z_message(molecular_family)
      else if (molecular_family%selected .and. .not. (settings%molecular_opacity_logt_min < &
           settings%molecular_opacity_logt_max)) then
         write(lines,'(a,2(1pe12.5,a))') 'kap_init: molecular_opacity_logt_min (', &
              settings%molecular_opacity_logt_min, ') must be below molecular_opacity_logt_max (', &
              settings%molecular_opacity_logt_max, ')'
      end if
      if (lines(1)(1:1) /= unwritten_line) then
         call kap_report(lines)
         ierr = 1
      end if

contains

logical function same_table_zs(family)
      type(kap_source), intent(in) :: family
      same_table_zs = family%selected .and. family%blends_two_z .and. &
           family%table_z(1) == family%table_z(2)
end function same_table_zs

subroutine two_z_message(family)
      type(kap_source), intent(in) :: family
      write(lines,'(a,1pe12.5)') 'kap_init: two ' // trim(family%name) // &
           ' tables (use_z_ramp or use_diffusion_z) need two Zs; ' // &
           trim(family%table_z_controls(1)) // ' = ' // trim(family%table_z_controls(2)) // &
           ' =', family%table_z(1)
end subroutine two_z_message

end subroutine check_kap_settings

!----------------------------------------------------------------------
! kap_update_surface_tables
!----------------------------------------------------------------------
! Rebuild the interior and molecular families' surface-composition
! tables (OPAL 1995-format, OPAL 1992, Alexander 1994) at a new envelope
! X, from the tables kap_init read.
subroutine kap_update_surface_tables(hydrogen_fraction, ierr)
      real(dp), intent(in) :: hydrogen_fraction
      integer, intent(out) :: ierr

      ierr = 0
      if (interior%selected .and. associated(interior%refresh)) then
         call interior%refresh(kap_settings, hydrogen_fraction, tables, ierr)
         if (ierr /= 0) return
      end if
      if (molecular%selected .and. associated(molecular%refresh)) then
         call molecular%refresh(kap_settings, hydrogen_fraction, tables, ierr)
      end if
end subroutine kap_update_surface_tables

end module kap_lib
