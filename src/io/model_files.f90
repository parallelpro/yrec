!----------------------------------------------------------------------
! model_files
!----------------------------------------------------------------------
! The stellar model files. write_mod_model writes the .mod format: the
! final model, which is also the restart file and the solver's recovery
! source, and the Monte Carlo runs' models. read_mod_model reads it back.
! read_yrec7 and read_model2 read the two older starting-model formats
! (read_starting_model, in core, picks the reader from the file's format
! tag, read_model_format_tag):
!   YREC7 (tag NMOD): a header, the physics flags, the luminosities, the
!     three envelope-triangle records, then the shells in four parts: the
!     structure with X and Z, one line per shell; He3 and CNO
!     (yrec7_cno_slots), one value per convection zone and one line per
!     radiative shell; with the extended composition the light elements
!     H2-Be9 (yrec7_light_slots), two shells per line; and in a rotating
!     model -log10(omega), eight shells per line.
!   MODEL2 (tag MOD2, llp 4/16/03): a header, the physics flags and
!     codes, the luminosities, the three envelope records, three heading
!     lines, then one line per shell with every quantity (the
!     composition in model2_slots).
! Both hold the seven luminosities of file_luminosities. Each record's
! layout is a named format below. The three readers fill star% directly
! and return the file's other header values in a model_file_header (the
! module reads no control: a reader takes what the file says).
module model_files
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use star_info_lib, only: i_h1, i_he4, i_metals, i_he3, i_c12, i_c13, &
           i_n14, i_n15, i_o16, i_o17, i_o18, i_h2, i_li6, i_li7, i_be9, &
           i_lum_pp1, i_lum_pp2, i_lum_pp3, i_lum_cno, i_lum_3alpha, &
           i_lum_neu, i_lum_grav
      implicit none
      private
      public :: read_model_format_tag, write_mod_model, read_mod_model, &
           read_yrec7, read_model2

! The header values of a model file that the star does not hold, as
! read_mod_model, read_yrec7 and read_model2 return them: each reader
! returns every value its layout holds, and a layout without one leaves
! its default. Who uses them:
!   read_starting_model     the mixing length and the extended-
!                           composition flag (checked against the
!                           controls), the timestep (the first step's);
!   read_yrec7 itself       the rotation flag and the convection-zone
!                           edges, which lay out the rest of the file;
!   none (the format's      the triangle orientation (reset_state_after_
!   record; test_io checks  read sets the run's to 1 after every read)
!   each is read from its   and MODEL2's codes, flags and disk and wind
!   place)                  constants.
      type, public :: model_file_header
! the model's mixing length (0 when the file has none) and timestep [yr]
         real(dp) :: mixing_length = 0.0d0
         real(dp) :: timestep_yr = 0.0d0
! what the model holds
         logical :: use_extended_composition = .false.
         logical :: rotation_active = .false.
! YREC7 and MODEL2: the convection-zone edges (shell indices), and -1
! when an envelope record number is negative (the triangle was
! reflected), 1 otherwise
         integer :: core_cz_top_index = 0, envelope_cz_bottom_index = 0
         real(dp) :: triangle_orientation = 1.0d0
! MODEL2's physics-option codes, read verbatim (their encoding is not
! interpreted; "?" in the other layouts); its flags, the switches of the
! run that wrote it, under their controls' names; and its disk and wind
! constants (TDISK, PDISK, WMAX): the disk-locking age [Gyr], the disk's
! angular velocity [rad/s] and the wind's saturation angular velocity
! [rad/s], as the controls disk_locking_age_gyr, disk_omega_rad_s and
! wind_saturation_omega hold them (their defaults in the other layouts)
         character(len=6) :: eos_code = '    ? '
         character(len=4) :: atm_code = '  ? ', alok_code = '  ? ', hik_code = '  ? '
         logical :: use_diffusion_y = .false., use_diffusion_z = .false., &
              disk_locking_active = .false., instability_transport_active = .false., &
              use_wind_torque = .false., core_overshoot_active = .false., &
              envelope_overshoot_active = .false., intermediate_overshoot_active = .false., &
              use_pure_z_table = .false., use_semiconvection = .false.
         real(dp) :: disk_locking_age_gyr = 0.0d0, disk_omega_rad_s = 7.2722d-6, &
              wind_saturation_omega = 3.0d-4
      end type model_file_header

! The .mod format's composition columns, by slot, and its column-header
! line (the writer writes it, the reader requires it).
      integer, parameter :: mod_slots(15) = [i_h1, i_he4, i_metals, i_he3, &
           i_c12, i_c13, i_n14, i_n15, i_o16, i_o17, i_o18, i_h2, &
           i_li6, i_li7, i_be9]
      character(len=*), parameter :: mod_columns = &
           'columns: zone log_mass logR luminosity logP logT' &
           // ' logRho conv omega h1 he4 z he3 c12 c13 n14 n15 o16 o17 o18' &
           // ' h2 li6 li7 be9'
! A .mod row: the zone, the six structure columns (log_mass .. logRho),
! the convective flag, then omega and the composition slots, as many as
! the row holds; a global value: its name padded to 20 characters, then
! the value.
      character(len=*), parameter :: mod_row_format = &
           '(i6,6(1x,es26.16e3),1x,l1,*(1x,es26.16e3))', &
           mod_real_format = '(a,1x,es26.16e3)', mod_integer_format = '(a,1x,i12)', &
           mod_logical_format = '(a,1x,l1)'

! The YREC7 layout's composition blocks: part 2's He3 and CNO isotopes,
! part 3's light elements, by slot, in file order.
      integer, parameter :: yrec7_cno_slots(8) = [i_he3, i_c12, i_c13, i_n14, &
           i_n15, i_o16, i_o17, i_o18]
      integer, parameter :: yrec7_light_slots(4) = [i_h2, i_li6, i_li7, i_be9]
! The MODEL2 layout's composition columns, by slot, in file order.
      integer, parameter :: model2_slots(15) = [i_h1, i_he4, i_metals, i_he3, &
           i_c12, i_c13, i_n14, i_n15, i_o16, i_o17, i_o18, i_h2, &
           i_li6, i_li7, i_be9]
! The luminosities both older layouts hold, in file order: pp I-III,
! CNO, 3 alpha, the neutrinos and gravity.
      integer, parameter :: file_luminosities(7) = [i_lum_pp1, i_lum_pp2, &
           i_lum_pp3, i_lum_cno, i_lum_3alpha, i_lum_neu, i_lum_grav]

! The YREC7 records: the header; the flags, convection-zone edges and
! mixing length; the luminosities; an envelope record; a shell of part
! 1; the lines of parts 2 and 3 (eight values, or four for a single
! shell); part 4.
      character(len=*), parameter :: &
           yrec7_header_format = '(A4,1X,2I5,3F8.5,1X,F14.10,F13.9,F13.0)', &
           yrec7_flags_format = '(2(L1,1X),5X,2I4,F7.3)', &
           yrec7_luminosity_format = '(10X,7E10.3)', &
           yrec7_envelope_format = '(3X,I2,F7.5,4F8.5,3E12.5)', &
           yrec7_shell_format = '(0PF13.10,F10.7,1PE14.7,0PF10.7,2F10.7,L1,2I6)', &
           yrec7_composition_format = '(8(1PE9.3,1X))', &
           yrec7_single_shell_format = '(4(1PE9.3,1X))', &
           yrec7_omega_format = '(0P8F10.7)'
! The MODEL2 records: the header, the flags and codes, the
! luminosities ('TLUMX', then the seven values), an envelope record
! ('ENV', the record number, ...), a shell.
      character(len=*), parameter :: &
           model2_header_format = '(A5,2I8,5F16.12,1PE18.10,0P2F16.12)', &
           model2_flags_format = '(2I8,F16.10,1X,A6,1X,3(A4,1X),L1,1X,A4,1X,' // &
                '11(L1,1X),3(1PE18.10))', &
           model2_luminosity_format = '(10X,1P7E17.9)', &
           model2_envelope_format = '(4X,I2,5F16.12,1P3E20.12)', &
           model2_shell_format = '(I6,0P2F18.14,1PE24.16,0P3F18.14,1PE24.16,1x,L1,' // &
                '3(0PF12.9),12(0PE16.8))'

contains

! ---------------------------------------------------------------
! The format tag of the model file on unit iread, its first four
! characters (YMOD, NMOD or MOD2), the unit rewound first. A file it
! cannot read (an empty one) is reported and returns ierr = 1.
subroutine read_model_format_tag(iread, tag, ierr)
      integer, intent(in) :: iread
      character(len=4), intent(out) :: tag
      integer, intent(out) :: ierr
      character(len=256) :: msg
      integer :: ios

      ierr = 0
      rewind iread
      read(iread,'(a4)',iostat=ios,iomsg=msg) tag
      call check_model_read(ios, msg, iread, 'the format tag', ierr)
end subroutine read_model_format_tag

! ---------------------------------------------------------------
! The .mod format holds exactly what a reload uses: the global scalars,
! the per-shell structure, composition and omega, and the two values the
! deck-consistency check needs (mixing length, extended-composition
! flag). Its layout is MESA-style: a version tag, named global values
! (one per line; readers ignore unknown names), a named column-header
! line, then one row per shell at full double precision, so a model
! reloaded between run cards round-trips exactly. Everything is read
! from star%. With replace, the unit is rewound first, so the file holds
! one model, the latest (the final-model file); without it the model is
! appended after the ones already written (the Monte Carlo model file,
! one model per run). The unit stays open. A failed write is reported
! ("cannot write the model file <path>: <reason>") and returns ierr = 1;
! the writes after it are skipped.
!   star reads: age_gyr, convective_flag, log_L_lsun, log_mass,
!     log_Teff, log_total_mass, logP, logR, logRho, logT,
!     luminosity_lsun, mixing_length_alpha, model_number, nz, omega,
!     physics%rotation_active, physics%use_extended_composition,
!     star_mass, timestep_yr, xa
!   star writes: none
subroutine write_mod_model(s, iwrite, replace, ierr)
      use star_info_lib, only: star_info
      use run_messages, only: log_line, check_io, solver_diagnostics
      type(star_info), intent(in) :: s
      integer, intent(in) :: iwrite
      logical, intent(in) :: replace
      integer, intent(inout) :: ierr
      character(len=256) :: msg, path
      character(len=80) :: line
      integer :: k, ios

      ios = 0
      if (replace) rewind(iwrite, iostat=ios, iomsg=msg)
      call put_line('YMOD 1')
      call put_integer('model_number        ', s%model_number)
      call put_integer('num_zones           ', s%nz)
      call put_real(   'star_mass_msun      ', s%star_mass)
      call put_real(   'log_total_mass      ', s%log_total_mass)
      call put_real(   'star_age_gyr        ', s%age_gyr)
      call put_real(   'timestep_yr         ', s%timestep_yr)
      call put_real(   'log_L               ', s%log_L_lsun)
      call put_real(   'log_Teff            ', s%log_Teff)
      call put_real(   'mixing_length_alpha ', s%mixing_length_alpha)
      call put_logical('extended_composition', s%physics%use_extended_composition)
      call put_logical('rotation_active     ', s%physics%rotation_active)
      call put_line(mod_columns)
      do k = 1, s%nz
         if (ios /= 0) exit
         write(iwrite,mod_row_format,iostat=ios,iomsg=msg) k, &
              s%log_mass(k), s%logR(k), s%luminosity_lsun(k), &
              s%logP(k), s%logT(k), s%logRho(k), &
              s%convective_flag(k), s%omega(k), &
              s%xa(mod_slots,k)
      end do
      if (ios == 0) flush(iwrite, iostat=ios, iomsg=msg)
      if (ios /= 0) then
         inquire(unit=iwrite, name=path)
         call check_io(ios, msg, 'write the model file', path, ierr)
         return
      end if

! the solver forensics' record of the write (an appended model, the
! Monte Carlo file's, with its age)
      if (solver_diagnostics()) then
         if (replace) then
            write(line,'(a,i5,a,i3)') ' DUMPED MODEL', s%model_number, '  FILE', iwrite
         else
            write(line,'(a,i5,a,f13.9,a,i3)') ' DUMPED MODEL', s%model_number, &
                 ' AGE', s%age_gyr, '  FILE', iwrite
         end if
         call log_line(trim(line))
      end if

contains

! A line, a named integer, real or logical global: written unless an
! earlier write failed.
subroutine put_line(text)
      character(len=*), intent(in) :: text
      if (ios == 0) write(iwrite,'(a)',iostat=ios,iomsg=msg) text
end subroutine put_line

subroutine put_integer(name, value)
      character(len=*), intent(in) :: name
      integer, intent(in) :: value
      if (ios == 0) write(iwrite,mod_integer_format,iostat=ios,iomsg=msg) name, value
end subroutine put_integer

subroutine put_real(name, value)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: value
      if (ios == 0) write(iwrite,mod_real_format,iostat=ios,iomsg=msg) name, value
end subroutine put_real

subroutine put_logical(name, value)
      character(len=*), intent(in) :: name
      logical, intent(in) :: value
      if (ios == 0) write(iwrite,mod_logical_format,iostat=ios,iomsg=msg) name, value
end subroutine put_logical

end subroutine write_mod_model

! ---------------------------------------------------------------
! Named global values (unknown names are noted in the run log and
! ignored, so the format can grow), the column-header line (it must be
! mod_columns), then one row per shell into star%. The header returns
! the timestep, the mixing length, and the extended-composition and
! rotation flags, which the deck-consistency check and the timestep
! seeding need. The luminosity breakdown is zeroed; the first
! energy-generation call computes it. A malformed file is reported and
! returns ierr = 1.
!   star reads: none
!   star writes: age_gyr, convective_flag, log_L_lsun, log_mass,
!     log_Teff, log_total_mass, logP, logR, logRho, logT,
!     luminosity_breakdown_lsun, luminosity_lsun, model_number, nz,
!     omega, star_mass, xa
subroutine read_mod_model(s, iread, header, ierr)
      use star_info_lib, only: star_info
      use run_messages, only: log_line
      type(star_info), intent(inout) :: s
      integer, intent(in) :: iread
      type(model_file_header), intent(out) :: header
      integer, intent(out) :: ierr

      character(len=1024) :: line
      character(len=64) :: name
      character(len=256) :: msg
      integer :: k, izone, ios

      ierr = 0
      rewind iread
      read(iread,'(a)',iostat=ios,iomsg=msg) line
      call check_model_read(ios, msg, iread, 'the format tag', ierr)
      if (ierr /= 0) return
      if (line(1:6) /= 'YMOD 1') then
         call report_model_error(iread, 'is not a YMOD 1 file', ierr)
         return
      end if

! named globals until the column-header line
      do
         read(iread,'(a)',iostat=ios) line
         if (ios /= 0) then
            call report_model_error(iread, 'ends before its columns line', ierr)
            return
         end if
         if (line(1:8) == 'columns:') exit
         read(line,*,iostat=ios) name
         if (ios /= 0 .or. len_trim(line) == 0) cycle
         select case (trim(name))
         case ('model_number')
            read(line,*,iostat=ios,iomsg=msg) name, s%model_number
         case ('num_zones')
            read(line,*,iostat=ios,iomsg=msg) name, s%nz
         case ('star_mass_msun')
            read(line,*,iostat=ios,iomsg=msg) name, s%star_mass
         case ('log_total_mass')
            read(line,*,iostat=ios,iomsg=msg) name, s%log_total_mass
         case ('star_age_gyr')
            read(line,*,iostat=ios,iomsg=msg) name, s%age_gyr
         case ('timestep_yr')
            read(line,*,iostat=ios,iomsg=msg) name, header%timestep_yr
         case ('log_L')
            read(line,*,iostat=ios,iomsg=msg) name, s%log_L_lsun
         case ('log_Teff')
            read(line,*,iostat=ios,iomsg=msg) name, s%log_Teff
         case ('mixing_length_alpha')
            read(line,*,iostat=ios,iomsg=msg) name, header%mixing_length
         case ('extended_composition')
            read(line,*,iostat=ios,iomsg=msg) name, header%use_extended_composition
         case ('rotation_active')
            read(line,*,iostat=ios,iomsg=msg) name, header%rotation_active
         case default
            call log_line(' read_mod_model: ignoring the unknown global ' // trim(name))
            ios = 0
         end select
         call check_model_read(ios, msg, iread, 'the global '//trim(name), ierr)
         if (ierr /= 0) return
      end do

      if (trim(line) /= mod_columns) then
         call report_model_error(iread, 'has the columns ' // trim(line(9:)), ierr)
         return
      end if
      call check_zone_count(s%nz, iread, ierr)
      if (ierr /= 0) return

      do k = 1, s%nz
         read(iread,*,iostat=ios,iomsg=msg) izone, &
              s%log_mass(k), s%logR(k), s%luminosity_lsun(k), &
              s%logP(k), s%logT(k), s%logRho(k), &
              s%convective_flag(k), s%omega(k), &
              s%xa(mod_slots,k)
         call check_model_read(ios, msg, iread, 'the shells', ierr)
         if (ierr /= 0) return
      end do

! (the energy generation's first call computes the breakdown)
      s%luminosity_breakdown_lsun = 0.0d0
end subroutine read_mod_model

! ---------------------------------------------------------------
! A YREC7-format model (the module header gives its layout) into star%;
! the header returns the file's flags, convection-zone edges, mixing
! length, timestep and triangle orientation, MODEL2's values at their
! defaults. xa and omega are zeroed first, for what the file does not
! hold.
!   star reads: none
!   star writes: age_gyr, convective_flag, envelope_fit_coeffs,
!     fit_point_logP, fit_point_logR, fit_point_logT, log_L_lsun,
!     log_mass, log_Teff, log_total_mass, logP, logR, logRho, logT,
!     luminosity_breakdown_lsun, luminosity_lsun, model_number, nz,
!     omega, star_mass, triangle_log_l_lsun, triangle_log_teff, xa
subroutine read_yrec7(s, iread, header, ierr)
      use run_messages, only: log_line
      use star_info_lib, only: star_info, json
      use math_lib
      type(star_info), intent(inout) :: s
      integer, intent(in) :: iread
      type(model_file_header), intent(out) :: header
      integer, intent(out) :: ierr

      character(len=4) :: header_keyword
      character(len=256) :: msg
      integer :: ios
      real(dp) :: omega_log10(json)
! -log10(omega) at or above this value marks a non-rotating shell
! (omega stored as 0)
      real(dp), parameter :: omega_log10_zero_sentinel = 58.9D0
      integer :: i, j, core_top, envelope_bottom
      integer :: envelope_record_number
      integer :: pair_start_index

      ierr = 0
      s%omega = 0.0d0
      s%xa = 0.0d0
      rewind(iread)

! the header
      read(iread,yrec7_header_format,iostat=ios,iomsg=msg) header_keyword, &
           s%model_number, s%nz, s%star_mass, s%log_Teff, s%log_L_lsun, &
           s%log_total_mass, s%age_gyr, header%timestep_yr
      call check_model_read(ios, msg, iread, 'the header', ierr)
      if (ierr /= 0) return
      call check_zone_count(s%nz, iread, ierr)
      if (ierr /= 0) return

! what the model holds (rotation, the extended composition), the
! convection-zone edges and the mixing length
      read(iread,yrec7_flags_format,iostat=ios,iomsg=msg) header%rotation_active, &
           header%use_extended_composition, header%core_cz_top_index, &
           header%envelope_cz_bottom_index, header%mixing_length
      call check_model_read(ios, msg, iread, 'the flags', ierr)
      if (ierr /= 0) return
! the convection-zone edges index the composition records below
      core_top = header%core_cz_top_index
      envelope_bottom = header%envelope_cz_bottom_index
      call check_zone_index(core_top, s%nz, iread, 'core convection zone top', ierr)
      if (ierr /= 0) return
      call check_zone_index(envelope_bottom, s%nz, iread, &
           'envelope convection zone bottom', ierr)
      if (ierr /= 0) return

      read(iread,yrec7_luminosity_format,iostat=ios,iomsg=msg) &
           s%luminosity_breakdown_lsun(file_luminosities)
      call check_model_read(ios, msg, iread, 'the luminosities', ierr)
      if (ierr /= 0) return
      call luminosities_to_solar_units(s)

! the three envelope-triangle records
      do i = 1, 3
         read(iread,yrec7_envelope_format,iostat=ios,iomsg=msg) envelope_record_number, &
              s%triangle_log_teff(i), s%triangle_log_l_lsun(i), s%fit_point_logP(i), &
              s%fit_point_logT(i), s%fit_point_logR(i), &
              (s%envelope_fit_coeffs(3*i-3+j), j = 1, 3)
         call check_model_read(ios, msg, iread, 'the envelope records', ierr)
         if (ierr /= 0) return
         if (envelope_record_number < 0) header%triangle_orientation = -1.0d0
      end do

! The shells, in four parts. Part 1, one line per shell: log m, log r,
! L, log P, log T, log rho, convective, and X and Z in units of 1e-6.
      do i = 1, s%nz
         call read_yrec7_shell(i)
         if (ierr /= 0) return
      end do

! Part 2, He3 and the CNO isotopes: one value for each convection zone
! (the core's and the envelope's), one line per shell between them.
      if (envelope_bottom <= core_top) then
         read(iread,yrec7_composition_format,iostat=ios,iomsg=msg) s%xa(yrec7_cno_slots,core_top)
         if (ios == 0) read(iread,yrec7_composition_format,iostat=ios,iomsg=msg) &
              s%xa(yrec7_cno_slots,envelope_bottom)
      else
         do j = core_top, envelope_bottom
            read(iread,yrec7_composition_format,iostat=ios,iomsg=msg) s%xa(yrec7_cno_slots,j)
            if (ios /= 0) exit
         end do
      end if
      call check_model_read(ios, msg, iread, 'the CNO abundances', ierr)
      if (ierr /= 0) return
      call fill_convection_zones(yrec7_cno_slots)
! He4 = 1 - X - Z - He3
      do i = 1, s%nz
         s%xa(i_he4,i) = 1.0D0 - s%xa(i_h1,i) - s%xa(i_metals,i) - &
              s%xa(i_he3,i)
      end do

! Part 3, with the extended composition: the light elements, one value
! for a fully convective star, otherwise two shells per line from the
! core convection zone's top to the envelope's bottom, the odd shell
! left over on a line of its own, the convection zones as in part 2.
      if (header%use_extended_composition) then
         if (envelope_bottom == 1) then
            read(iread,yrec7_single_shell_format,iostat=ios,iomsg=msg) &
                 s%xa(yrec7_light_slots,1)
            call check_model_read(ios, msg, iread, 'the light elements', ierr)
            if (ierr /= 0) return
            do j = 1, s%nz
               s%xa(yrec7_light_slots,j) = s%xa(yrec7_light_slots,1)
            end do
         else
            ios = 0
            do pair_start_index = core_top, envelope_bottom-1, 2
               read(iread,yrec7_composition_format,iostat=ios,iomsg=msg) &
                    (s%xa(yrec7_light_slots,j), j = pair_start_index, pair_start_index+1)
               if (ios /= 0) exit
            end do
            call check_model_read(ios, msg, iread, 'the light elements', ierr)
            if (ierr /= 0) return
! an odd number of shells: the last on a line of its own (the loop's
! index has gone one pair past)
            if (mod(envelope_bottom-1 - core_top, 2) /= 0) then
               read(iread,yrec7_single_shell_format,iostat=ios,iomsg=msg) &
                    s%xa(yrec7_light_slots,pair_start_index)
               call check_model_read(ios, msg, iread, 'the light elements', ierr)
               if (ierr /= 0) return
            end if
            call fill_convection_zones(yrec7_light_slots)
         end if
      end if

! Part 4, in a rotating model: -log10(omega), eight shells per line. A
! file that ends before them (a model written before its rotation was
! set up) keeps omega zero, with a note.
      if (header%rotation_active) then
         read(iread,yrec7_omega_format,iostat=ios,iomsg=msg) (omega_log10(i), i = 1, s%nz)
         if (is_iostat_end(ios)) then
            call log_line(' read_yrec7: the model file has no omega records; omega is zero')
            s%omega = 0.0d0
            rewind iread
            return
         end if
         call check_model_read(ios, msg, iread, 'the rotation rates', ierr)
         if (ierr /= 0) return
         do i = 1, s%nz
            if (omega_log10(i) < omega_log10_zero_sentinel) then
               s%omega(i) = pow(10D0, (-omega_log10(i)))
            else
               s%omega(i) = 0.0d0
            end if
         end do
      end if

! (the unit stays open, rewound)
      rewind iread

contains

! Shell i of part 1, with its mass checked against the star's.
subroutine read_yrec7_shell(i)
      integer, intent(in) :: i
      integer :: ix, iz

      read(iread,yrec7_shell_format,iostat=ios,iomsg=msg) s%log_mass(i), s%logR(i), &
           s%luminosity_lsun(i), s%logP(i), s%logT(i), s%logRho(i), &
           s%convective_flag(i), ix, iz
      call check_model_read(ios, msg, iread, 'the shells', ierr)
      if (ierr /= 0) return
      if (s%log_mass(i) < 0D0 .or. s%log_mass(i) > s%log_total_mass) then
         write(msg, '(a,i0,a)') 'shell ', i, &
              ' has a mass below zero or above the star''s'
         call report_model_error(iread, trim(msg), ierr)
         return
      end if
      s%xa(i_h1,i) = 1.0D-6*dble(ix)
      s%xa(i_metals,i) = 1.0D-6*dble(iz)
end subroutine read_yrec7_shell

! The slots' values through the convection zones: the core convection
! zone's shells take its top shell's, the envelope's its bottom shell's.
subroutine fill_convection_zones(slots)
      integer, intent(in) :: slots(:)
      do j = 1, core_top - 1
         s%xa(slots,j) = s%xa(slots,core_top)
      end do
      do j = envelope_bottom + 1, s%nz
         s%xa(slots,j) = s%xa(slots,envelope_bottom)
      end do
end subroutine fill_convection_zones

end subroutine read_yrec7

! ---------------------------------------------------------------
! A MODEL2-format model (the module header gives its layout) into
! star%; the header returns the file's other header values.
!   star reads: none
!   star writes: age_gyr, convective_flag, envelope_fit_coeffs,
!     fit_point_logP, fit_point_logR, fit_point_logT, log_L_lsun,
!     log_mass, log_Teff, log_total_mass, logP, logR, logRho, logT,
!     luminosity_breakdown_lsun, luminosity_lsun, model_number, nz,
!     omega, star_mass, triangle_log_l_lsun, triangle_log_teff, xa
subroutine read_model2(s, iread, header, ierr)
      use star_info_lib, only: star_info
      type(star_info), intent(inout) :: s
      integer, intent(in) :: iread
      type(model_file_header), intent(out) :: header
      integer, intent(out) :: ierr

      character(len=5) :: header_keyword
! the model's initial-composition code, read and not used
      character(len=4) :: composition_code
      character(len=256) :: msg
      integer :: envelope_record_number
      integer :: i, j, k, ios, model_number, nz
      real(dp) :: log_mass_center, log_mass_surface

      ierr = 0
      rewind(iread)

! the header (the zone count checked before it indexes anything)
      read(iread,model2_header_format,iostat=ios,iomsg=msg) header_keyword, model_number, &
           nz, s%star_mass, s%log_Teff, s%log_L_lsun, s%log_total_mass, s%age_gyr, &
           header%timestep_yr, log_mass_center, log_mass_surface
      call check_model_read(ios, msg, iread, 'the header', ierr)
      if (ierr /= 0) return
      call check_zone_count(nz, iread, ierr)
      if (ierr /= 0) return
      s%model_number = model_number
      s%nz = nz
      s%log_mass(1) = log_mass_center
      s%log_mass(nz) = log_mass_surface

! the convection-zone edges, the mixing length, the physics codes and
! flags, and what the model holds (rotation, the extended composition)
      read(iread,model2_flags_format,iostat=ios,iomsg=msg) header%core_cz_top_index, &
           header%envelope_cz_bottom_index, header%mixing_length, header%eos_code, &
           header%atm_code, header%alok_code, header%hik_code, header%use_pure_z_table, &
           composition_code, header%use_extended_composition, header%use_diffusion_y, &
           header%use_diffusion_z, header%use_semiconvection, header%core_overshoot_active, &
           header%envelope_overshoot_active, header%intermediate_overshoot_active, &
           header%rotation_active, header%instability_transport_active, &
           header%use_wind_torque, header%disk_locking_active, &
           header%disk_locking_age_gyr, header%disk_omega_rad_s, &
           header%wind_saturation_omega
      call check_model_read(ios, msg, iread, 'the flags', ierr)
      if (ierr /= 0) return

      read(iread,model2_luminosity_format,iostat=ios,iomsg=msg) &
           s%luminosity_breakdown_lsun(file_luminosities)
      call check_model_read(ios, msg, iread, 'the luminosities', ierr)
      if (ierr /= 0) return
      call luminosities_to_solar_units(s)

! the three envelope-triangle records
      do i = 1, 3
         read(iread,model2_envelope_format,iostat=ios,iomsg=msg) envelope_record_number, &
              s%triangle_log_teff(i), s%triangle_log_l_lsun(i), s%fit_point_logP(i), &
              s%fit_point_logT(i), s%fit_point_logR(i), &
              (s%envelope_fit_coeffs(3*i-3+j), j = 1, 3)
         call check_model_read(ios, msg, iread, 'the envelope records', ierr)
         if (ierr /= 0) return
         if (envelope_record_number < 0) header%triangle_orientation = -1.0d0
      end do

! the three heading lines
      do i = 1, 3
         read(iread,*,iostat=ios,iomsg=msg)
         call check_model_read(ios, msg, iread, 'the column headings', ierr)
         if (ierr /= 0) return
      end do

! one line per shell
      do i = 1, s%nz
         read(iread,model2_shell_format,iostat=ios,iomsg=msg) k, s%log_mass(i), s%logR(i), &
              s%luminosity_lsun(i), s%logP(i), s%logT(i), s%logRho(i), s%omega(i), &
              s%convective_flag(i), s%xa(model2_slots,i)
         call check_model_read(ios, msg, iread, 'the shells', ierr)
         if (ierr /= 0) return
      end do
end subroutine read_model2

! ---------------------------------------------------------------
! The file's seven luminosities in Lsun: a model that stores them in
! erg/s (any above 1e20, |L_nu| for the neutrinos) is converted.
!   star reads: solar_luminosity_cgs
!   star writes: luminosity_breakdown_lsun
subroutine luminosities_to_solar_units(s)
      use star_info_lib, only: star_info
      type(star_info), intent(inout) :: s
! the largest luminosity in Lsun a model can hold [Lsun]
      real(dp), parameter :: largest_solar_units = 1.0D20
      real(dp) :: largest

      largest = max(maxval(s%luminosity_breakdown_lsun(file_luminosities)), &
           abs(s%luminosity_breakdown_lsun(i_lum_neu)))
      if (largest > largest_solar_units) &
           s%luminosity_breakdown_lsun(file_luminosities) = &
           s%luminosity_breakdown_lsun(file_luminosities)/s%solar_luminosity_cgs
end subroutine luminosities_to_solar_units

! ---------------------------------------------------------------
! An error in the model file on unit iread: "the starting model <file>
! <what>", and ierr = 1.
subroutine report_model_error(iread, what, ierr)
      use run_messages, only: report_error
      integer, intent(in) :: iread
      character(len=*), intent(in) :: what
      integer, intent(inout) :: ierr
      character(len=256) :: path

      inquire(unit=iread, name=path)
      call report_error('the starting model ' // trim(path) // ' ' // what)
      ierr = 1
end subroutine report_model_error

! ---------------------------------------------------------------
! After a read of the model file on unit iread with iostat=ios and
! iomsg=msg: on failure, report "cannot read <what> of the starting
! model <file>: <msg>" (check_io) and set ierr = 1.
subroutine check_model_read(ios, msg, iread, what, ierr)
      use run_messages, only: check_io
      integer, intent(in) :: ios, iread
      character(len=*), intent(in) :: msg, what
      integer, intent(inout) :: ierr
      character(len=256) :: path

      if (ios == 0) return
      inquire(unit=iread, name=path)
      call check_io(ios, msg, 'read ' // what // ' of the starting model', path, ierr)
end subroutine check_model_read

! ---------------------------------------------------------------
! A model's zone count must fit the star's arrays (1 to json zones);
! otherwise report it and set ierr = 1.
subroutine check_zone_count(nz, iread, ierr)
      use star_info_lib, only: json
      use run_messages, only: report_error
      integer, intent(in) :: nz, iread
      integer, intent(inout) :: ierr
      character(len=256) :: path
      character(len=64) :: counts

      if (nz >= 1 .and. nz <= json) return
      inquire(unit=iread, name=path)
      write(counts, '(i0,a,i0)') nz, ' zones, not 1 to ', json
      call report_error('the starting model ' // trim(path) // ' has ' // trim(counts))
      ierr = 1
end subroutine check_zone_count

! ---------------------------------------------------------------
! A zone index from a model's header must lie in 1..nz; otherwise
! report it and set ierr = 1.
subroutine check_zone_index(index, nz, iread, what, ierr)
      use run_messages, only: report_error
      integer, intent(in) :: index, nz, iread
      character(len=*), intent(in) :: what
      integer, intent(inout) :: ierr
      character(len=256) :: path
      character(len=64) :: numbers

      if (index >= 1 .and. index <= nz) return
      inquire(unit=iread, name=path)
      write(numbers, '(i0,a,i0)') index, ', not 1 to ', nz
      call report_error('the starting model ' // trim(path) // ': ' // what // &
           ' at zone ' // trim(numbers))
      ierr = 1
end subroutine check_zone_index

end module model_files
