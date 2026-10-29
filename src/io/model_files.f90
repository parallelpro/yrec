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
!     structure with X and Z, one line per shell; He3 and CNO, one value
!     per convection zone and one line per radiative shell; the light
!     elements H2-Be9, two shells per line; and -log10(omega), eight
!     shells per line.
!   MODEL2 (tag MOD2, llp 4/16/03): a header, the physics flags and
!     codes, the luminosities, the three envelope records, three heading
!     lines, then one line per shell with every quantity.
! Both readers take the same arguments, the file's header values (the
! last two lines of arguments are MODEL2's; read_yrec7 gives them their
! defaults). They fill star% directly; several of the header values
! share names with controls_buffer members, so this module never uses
! controls_buffer.
module model_files
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use star_info_lib, only: i_h1, i_he4, i_metals, i_he3, i_c12, i_c13, &
           i_n14, i_n15, i_o16, i_o17, i_o18, i_h2, i_li6, i_li7, i_be9
      implicit none
      private
      public :: read_model_format_tag, write_mod_model, read_mod_model, &
           read_yrec7, read_model2

! The .mod format's composition columns, by slot, and its column-header
! line (the writer writes it, the reader requires it).
      integer, parameter :: mod_slots(15) = [i_h1, i_he4, i_metals, i_he3, &
           i_c12, i_c13, i_n14, i_n15, i_o16, i_o17, i_o18, i_h2, &
           i_li6, i_li7, i_be9]
      character(len=*), parameter :: mod_columns = &
           'columns: zone log_mass logR luminosity logP logT' &
           // ' logRho conv omega h1 he4 z he3 c12 c13 n14 n15 o16 o17 o18' &
           // ' h2 li6 li7 be9'

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
! one model per run).
subroutine write_mod_model(s, iwrite, replace)
      use star_info_lib, only: star_info
      use io_units, only: last_model_unit, run_log_unit
      use run_log, only: solver_diagnostics
      type(star_info), intent(in) :: s
      integer, intent(in) :: iwrite
      logical, intent(in) :: replace
      integer :: k

      if (replace) rewind iwrite
      write(iwrite,'(a)') 'YMOD 1'
      write(iwrite,'(a,1x,i12)')        'model_number        ', s%model_number
      write(iwrite,'(a,1x,i12)')        'num_zones           ', s%nz
      write(iwrite,'(a,1x,es26.16e3)')  'star_mass_msun      ', s%star_mass
      write(iwrite,'(a,1x,es26.16e3)')  'log_total_mass      ', s%log_total_mass
      write(iwrite,'(a,1x,es26.16e3)')  'star_age_gyr        ', s%age_gyr
      write(iwrite,'(a,1x,es26.16e3)')  'timestep_yr         ', s%timestep_yr
      write(iwrite,'(a,1x,es26.16e3)')  'log_L               ', s%log_L_lsun
      write(iwrite,'(a,1x,es26.16e3)')  'log_Teff            ', s%log_Teff
      write(iwrite,'(a,1x,es26.16e3)')  'mixing_length_alpha ', s%mixing_length_alpha
      write(iwrite,'(a,1x,l1)')         'extended_composition', s%ctrl%use_extended_composition
      write(iwrite,'(a,1x,l1)')         'rotation_active     ', s%ctrl%rotation_active
      write(iwrite,'(a)') mod_columns
      do k = 1, s%nz
         write(iwrite,'(i6,6(1x,es26.16e3),1x,l1,16(1x,es26.16e3))') k, &
              s%log_mass(k), s%logR(k), s%luminosity_lsun(k), &
              s%logP(k), s%logT(k), s%logRho(k), &
              s%convective_flag(k), s%omega(k), &
              s%xa(mod_slots,k)
      end do
      flush(iwrite)

! the solver forensics' record of the write
      if (solver_diagnostics()) then
         if (iwrite == last_model_unit) then
          write(run_log_unit,330) s%model_number, iwrite
         else
          write(run_log_unit,340) s%model_number, s%age_gyr, iwrite
         endif
      end if
  330 format(' DUMPED MODEL',I5,'  FILE',I3)
  340 format(' DUMPED MODEL',I5,' AGE',F13.9,'  FILE',I3)
end subroutine write_mod_model

! ---------------------------------------------------------------
! Named global values (unknown names are noted in the run log and
! ignored, so the format can grow), the column-header line (it must be
! mod_columns), then one row per shell into star%. The four header
! values the deck-consistency check and the timestep seeding need are
! returned: the timestep, the mixing length, and the extended-
! composition and rotation flags. The luminosity breakdown is zeroed;
! the first energy-generation call computes it. A malformed file is
! reported and returns ierr = 1.
subroutine read_mod_model(s, iread, timestep_yr, mixing_length0, &
     use_extended_composition0, rotation_active0, ierr)
      use star_info_lib, only: star_info
      use io_units, only: run_log_unit
      type(star_info), intent(inout) :: s
      integer, intent(in) :: iread
      real(dp), intent(out) :: timestep_yr, mixing_length0
      logical, intent(out) :: use_extended_composition0, rotation_active0
      integer, intent(out) :: ierr

      character(len=1024) :: line
      character(len=64) :: name
      character(len=256) :: msg
      integer :: k, izone, ios

      ierr = 0
      timestep_yr = 0.0d0
      mixing_length0 = 0.0d0
      use_extended_composition0 = .false.
      rotation_active0 = .false.

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
            read(line,*,iostat=ios,iomsg=msg) name, timestep_yr
         case ('log_L')
            read(line,*,iostat=ios,iomsg=msg) name, s%log_L_lsun
         case ('log_Teff')
            read(line,*,iostat=ios,iomsg=msg) name, s%log_Teff
         case ('mixing_length_alpha')
            read(line,*,iostat=ios,iomsg=msg) name, mixing_length0
         case ('extended_composition')
            read(line,*,iostat=ios,iomsg=msg) name, use_extended_composition0
         case ('rotation_active')
            read(line,*,iostat=ios,iomsg=msg) name, rotation_active0
         case default
            write(run_log_unit,'(2a)') &
                 ' read_mod_model: ignoring the unknown global ', trim(name)
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
! A YREC7-format model (the module header gives its layout); the
! arguments are the file's header values, MODEL2's at their defaults.
! xa and omega are zeroed first, for what the file does not hold.
subroutine read_yrec7(s, mixing_length, timestep_yr, triangle_orientation, &
     iread, core_cz_top_index, envelope_cz_bottom_index, &
     use_extended_composition, rotation_active, atm_code, eos_code, &
     hik_code, use_diffusion_y, use_diffusion_z, &
     disk_locking_active, instability_transport_active, &
     use_wind_torque, alok_code, core_overshoot_active, &
     envelope_overshoot_active, lovstm, use_pure_z_table, &
     use_semiconvection, disk_pressure, disk_temperature, &
     wind_saturation_omega, ierr)
      use io_units, only: run_log_unit
      use star_info_lib, only: star_info, json
      use math_lib
      type(star_info), intent(inout) :: s

      real(dp), intent(out) :: mixing_length
      real(dp), intent(out) :: timestep_yr
      real(dp), intent(out) :: triangle_orientation
      integer, intent(in) :: iread
      integer, intent(out) :: core_cz_top_index, envelope_cz_bottom_index
      logical, intent(out) :: use_extended_composition, rotation_active
! atm_code, eos_code, hik_code, alok_code: the MODEL2 format's
! physics-option codes, "?" here (their encoding is not interpreted)
      character(len=4), intent(out) :: atm_code
      character(len=6), intent(out) :: eos_code
      character(len=4), intent(out) :: hik_code
      logical, intent(out) :: use_diffusion_y, use_diffusion_z, &
           disk_locking_active, instability_transport_active, use_wind_torque
      character(len=4), intent(out) :: alok_code
      logical, intent(out) :: core_overshoot_active, envelope_overshoot_active, lovstm, &
           use_pure_z_table, use_semiconvection
      real(dp), intent(out) :: disk_pressure, disk_temperature, &
           wind_saturation_omega
      integer, intent(out) :: ierr

      character(len=4) :: header_keyword
      character(len=256) :: msg
      integer :: ios
      real(dp) :: omega_log10(json)
! -log10(omega) at or above this value marks a non-rotating shell
! (omega stored as 0)
      real(dp), parameter :: omega_log10_zero_sentinel = 58.9D0
      integer :: i, j, ix, iz
      integer :: envelope_record_number
      integer :: pair_start_index, parity_test
      integer :: ii

      ierr = 0

      s%omega = 0d0
      s%xa = 0d0

! MODEL2's codes, flags and constants, at their defaults
      atm_code = "  ? "
      eos_code = "    ? "
      hik_code = "  ? "
      alok_code = "  ? "
      use_diffusion_y = .false.
      use_diffusion_z = .false.
      disk_locking_active = .false.
      instability_transport_active = .false.
      use_wind_torque = .false.
      core_overshoot_active = .false.
      envelope_overshoot_active = .false.
      lovstm = .false.
      use_pure_z_table = .false.
      use_semiconvection = .false.
      disk_pressure = 7.2722D-6
      disk_temperature = 0d0
      wind_saturation_omega = 3d-4

      rewind(iread)

! the header
      read(iread,10,iostat=ios,iomsg=msg) header_keyword,s%model_number,s%nz,s%star_mass, &
           s%log_Teff,s%log_L_lsun,s%log_total_mass, &
           s%age_gyr,timestep_yr
 10   format(A4,1X,2I5,3F8.5,1X,F14.10,F13.9,F13.0)
      call check_model_read(ios, msg, iread, 'the header', ierr)
      if (ierr /= 0) return
      call check_zone_count(s%nz, iread, ierr)
      if (ierr /= 0) return

! what the model holds (rotation, the extended composition), the
! convection-zone edges and the mixing length
      read(iread,30,iostat=ios,iomsg=msg) rotation_active,use_extended_composition, &
           core_cz_top_index,envelope_cz_bottom_index,mixing_length
 30   format(2(L1,1X),5X,2I4,F7.3)
      call check_model_read(ios, msg, iread, 'the flags', ierr)
      if (ierr /= 0) return
! the convection-zone edges index the composition records below
      call check_zone_index(core_cz_top_index, s%nz, iread, 'core convection zone top', ierr)
      if (ierr /= 0) return
      call check_zone_index(envelope_cz_bottom_index, s%nz, iread, &
           'envelope convection zone bottom', ierr)
      if (ierr /= 0) return

! the luminosities of pp I-III, CNO, 3 alpha, the neutrinos and gravity
      read(iread,40,iostat=ios,iomsg=msg) (s%luminosity_breakdown_lsun(j),j=1,7)
 40   format(10X,7E10.3)
      call check_model_read(ios, msg, iread, 'the luminosities', ierr)
      if (ierr /= 0) return
      call luminosities_to_solar_units(s)

! the three envelope-triangle records; triangle_orientation is -1 when a
! record number is negative (the triangle was reflected), 1 otherwise
      triangle_orientation = 1D0
      do i = 1,3
       read(iread,70,iostat=ios,iomsg=msg) envelope_record_number,s%triangle_log_teff(i), &
            s%triangle_log_l_lsun(i),s%fit_point_pressure(i), &
            s%fit_point_temperature(i),s%fit_point_radius(i), &
            (s%envelope_fit_coeffs(3*i-3+j),j = 1,3)
 70      format(3X,I2,F7.5,4F8.5,3E12.5)
       call check_model_read(ios, msg, iread, 'the envelope records', ierr)
       if (ierr /= 0) return
       if (envelope_record_number < 0) triangle_orientation = -1D0
      end do

! The shells, in four parts. Part 1, one line per shell: log m, log r,
! L, log P, log T, log rho, convective, and X and Z in units of 1e-6.
      do i = 1,s%nz
       read(iread,100,iostat=ios,iomsg=msg) s%log_mass(i),s%logR(i),s%luminosity_lsun(i), &
            s%logP(i),s%logT(i),s%logRho(i),s%convective_flag(i), &
            ix,iz
 100     format(0PF13.10,F10.7,1PE14.7,0PF10.7,2F10.7,L1,2I6)
       call check_model_read(ios, msg, iread, 'the shells', ierr)
       if (ierr /= 0) return
       if (s%log_mass(i) < 0D0 .or. s%log_mass(i) > s%log_total_mass) then
          write(msg, '(a,i0,a)') 'shell ', i, &
               ' has a mass below zero or above the star''s'
          call report_model_error(iread, trim(msg), ierr)
          return
       endif
       s%xa(i_h1,i) = 1.0D-6*dble(ix)
       s%xa(i_metals,i) = 1.0D-6*dble(iz)
      end do

! Part 2, He3 and the CNO isotopes (slots 4-11): one value for each
! convection zone (the core's and the envelope's), one line per shell
! between them.
      if (envelope_cz_bottom_index <= core_cz_top_index) then
       read(iread,200,iostat=ios,iomsg=msg)(s%xa(i,core_cz_top_index),i=4,11)
       if (ios == 0) read(iread,200,iostat=ios,iomsg=msg)(s%xa(i,envelope_cz_bottom_index),i=4,11)
      else
       do j = core_cz_top_index,envelope_cz_bottom_index
          read(iread,200,iostat=ios,iomsg=msg) (s%xa(i,j),i = 4,11)
 200        format(8(1PE9.3,1X))
          if (ios /= 0) exit
       end do
      endif
      call check_model_read(ios, msg, iread, 'the CNO abundances', ierr)
      if (ierr /= 0) return
! the convective core's value in its shells
      if (core_cz_top_index > 1) then
       do j = 1,core_cz_top_index-1
          do i = 4,11
             s%xa(i,j) = s%xa(i,core_cz_top_index)
          end do
       end do
      endif
! the convective envelope's value in its shells
      if (envelope_cz_bottom_index < s%nz) then
       do j = envelope_cz_bottom_index+1,s%nz
          do i = 4,11
             s%xa(i,j) = s%xa(i,envelope_cz_bottom_index)
          end do
       end do
      endif
! He4 = 1 - X - Z - He3
      do i = 1,s%nz
       s%xa(i_he4,i) = 1.0D0 - s%xa(i_h1,i) - s%xa(i_metals,i) - &
            s%xa(i_he3,i)
      end do

! Part 3, with the extended composition: H2, Li6, Li7 and Be9 (slots
! 12-15), one value for a fully convective star, otherwise two shells
! per line, the convection zones as in part 2.
      if (use_extended_composition) then
       if (envelope_cz_bottom_index == 1) then
          read(iread,300,iostat=ios,iomsg=msg)(s%xa(i,1),i=12,15)
 300        format(4(1PE9.3,1X))
          call check_model_read(ios, msg, iread, 'the light elements', ierr)
          if (ierr /= 0) return
          do j = 1,s%nz
             do i = 12,15
              s%xa(i,j) = s%xa(i,1)
             end do
          end do
       else
          ios = 0
          do pair_start_index = core_cz_top_index,envelope_cz_bottom_index-1,2
             read(iread,200,iostat=ios,iomsg=msg)((s%xa(i,j),i = 12,15), &
                  j = pair_start_index,pair_start_index+1)
             if (ios /= 0) exit
          end do
          call check_model_read(ios, msg, iread, 'the light elements', ierr)
          if (ierr /= 0) return
! an odd number of shells: the last on a line of its own
          parity_test = envelope_cz_bottom_index-1 - core_cz_top_index
          if (mod(parity_test,2) /= 0) then
               read(iread,300,iostat=ios,iomsg=msg)(s%xa(i,pair_start_index),i = 12,15)
               call check_model_read(ios, msg, iread, 'the light elements', ierr)
               if (ierr /= 0) return
          end if
          if (core_cz_top_index > 1) then
             do j = 1,core_cz_top_index-1
              do i = 12,15
                 s%xa(i,j) = s%xa(i,core_cz_top_index)
              end do
             end do
          endif
          if (envelope_cz_bottom_index < s%nz) then
             do j = envelope_cz_bottom_index+1,s%nz
              do i = 12,15
                 s%xa(i,j) = s%xa(i,envelope_cz_bottom_index)
              end do
             end do
          endif
       endif
      endif

! Part 4, in a rotating model: -log10(omega), eight shells per line. A
! file that ends before them (a model written before its rotation was
! set up) keeps omega zero, with a note.
      if (rotation_active) then
      read(iread,500,iostat=ios,iomsg=msg)(omega_log10(ii),ii = 1,s%nz)
 500    format(0P8F10.7)
      if (is_iostat_end(ios)) then
         write(run_log_unit,'(a)') ' read_yrec7: the model file has no omega' // &
              ' records; omega is zero'
         s%omega = 0d0
         rewind iread
         return
      end if
      call check_model_read(ios, msg, iread, 'the rotation rates', ierr)
      if (ierr /= 0) return
      do i = 1,s%nz
           if (omega_log10(i) < omega_log10_zero_sentinel) then
              s%omega(i) = pow(10D0, (-omega_log10(i)))
           else
              s%omega(i) = 0D0
           endif
      end do
      endif

! (the unit stays open, rewound)
      rewind iread
end subroutine read_yrec7

! ---------------------------------------------------------------
! A MODEL2-format model (the module header gives its layout); the
! arguments are the file's header values.
subroutine read_model2(s, mixing_length, timestep_yr, triangle_orientation, &
     iread, core_cz_top_index, envelope_cz_bottom_index, &
     use_extended_composition, rotation_active, atm_code, eos_code, &
     hik_code, use_diffusion_y, use_diffusion_z, &
     disk_locking_active, instability_transport_active, &
     use_wind_torque, alok_code, core_overshoot_active, &
     envelope_overshoot_active, lovstm, use_pure_z_table, &
     use_semiconvection, disk_pressure, disk_temperature, &
     wind_saturation_omega, ierr)
      use star_info_lib, only: star_info
      type(star_info), intent(inout) :: s
      real(dp), intent(out) :: mixing_length
      real(dp), intent(out) :: timestep_yr
      real(dp), intent(out) :: triangle_orientation
      integer, intent(in) :: iread
      integer, intent(out) :: core_cz_top_index, envelope_cz_bottom_index
      logical, intent(out) :: use_extended_composition, rotation_active
! atm_code, eos_code, hik_code, alok_code: short physics-option codes,
! read verbatim from the header (their encoding is not interpreted)
      character(len=4), intent(out) :: atm_code
      character(len=6), intent(out) :: eos_code
      character(len=4), intent(out) :: hik_code
      logical, intent(out) :: use_diffusion_y, use_diffusion_z, &
           disk_locking_active, instability_transport_active, use_wind_torque
      character(len=4), intent(out) :: alok_code
      logical, intent(out) :: core_overshoot_active, envelope_overshoot_active, lovstm, &
           use_pure_z_table, use_semiconvection
      real(dp), intent(out) :: disk_pressure, disk_temperature, &
           wind_saturation_omega
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
      read(iread,10,iostat=ios,iomsg=msg) header_keyword,model_number,nz,s%star_mass, &
           s%log_Teff,s%log_L_lsun,s%log_total_mass,s%age_gyr, &
           timestep_yr,log_mass_center,log_mass_surface
 10   format(A5,2I8,5F16.12,1PE18.10,0P2F16.12)
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
      read(iread,30,iostat=ios,iomsg=msg) core_cz_top_index,envelope_cz_bottom_index, &
           mixing_length,eos_code,atm_code,alok_code,hik_code, &
           use_pure_z_table,composition_code,use_extended_composition, &
           use_diffusion_y,use_diffusion_z,use_semiconvection,core_overshoot_active, &
           envelope_overshoot_active,lovstm,rotation_active, &
           instability_transport_active,use_wind_torque,disk_locking_active, &
           disk_temperature,disk_pressure,wind_saturation_omega
   30 format(2I8,F16.10,1X,A6,1X,3(A4,1X),L1,1X,A4,1X,11(L1,1X), &
           3(1PE18.10))
      call check_model_read(ios, msg, iread, 'the flags', ierr)
      if (ierr /= 0) return

! the luminosities of pp I-III, CNO, 3 alpha, the neutrinos and gravity
! (the line 'TLUMX', then the seven values)
      read(iread,40,iostat=ios,iomsg=msg) (s%luminosity_breakdown_lsun(j),j=1,7)
   40 format(10X,1P7E17.9)
      call check_model_read(ios, msg, iread, 'the luminosities', ierr)
      if (ierr /= 0) return
      call luminosities_to_solar_units(s)

! the three envelope-triangle records ('ENV', the record number, ...);
! triangle_orientation is -1 when a record number is negative (the
! triangle was reflected), 1 otherwise
      triangle_orientation = 1D0
      do i = 1,3
         read(iread,50,iostat=ios,iomsg=msg)envelope_record_number,s%triangle_log_teff(i), &
              s%triangle_log_l_lsun(i),s%fit_point_pressure(i), &
              s%fit_point_temperature(i),s%fit_point_radius(i), &
              (s%envelope_fit_coeffs(3*i-3+j),j=1,3)
         call check_model_read(ios, msg, iread, 'the envelope records', ierr)
         if (ierr /= 0) return
         if (envelope_record_number < 0) triangle_orientation = -1D0
      enddo
   50 format(4X,I2,5F16.12,1P3E20.12)

! the three heading lines
      do i = 1, 3
         read(iread,*,iostat=ios,iomsg=msg)
         call check_model_read(ios, msg, iread, 'the column headings', ierr)
         if (ierr /= 0) return
      end do

! one line per shell
      do i = 1,s%nz
         read(iread,60,iostat=ios,iomsg=msg) k,s%log_mass(i),s%logR(i),s%luminosity_lsun(i), &
              s%logP(i),s%logT(i),s%logRho(i), &
              s%omega(i),s%convective_flag(i),(s%xa(j,i),j=1,15)
         call check_model_read(ios, msg, iread, 'the shells', ierr)
         if (ierr /= 0) return
       enddo
 60   format(I6,0P2F18.14,1PE24.16,0P3F18.14,1PE24.16,1x,L1, &
           3(0PF12.9),12(0PE16.8))
end subroutine read_model2

! ---------------------------------------------------------------
! The seven luminosities in Lsun: a model that stores them in erg/s
! (any above 1e20, |L_nu| for the neutrinos) is converted.
subroutine luminosities_to_solar_units(s)
      use star_info_lib, only: star_info, i_lum_3alpha, i_lum_cno, &
           i_lum_grav, i_lum_neu, i_lum_pp1, i_lum_pp2, i_lum_pp3
      type(star_info), intent(inout) :: s
! the largest luminosity in Lsun a model can hold [Lsun]
      real(dp), parameter :: largest_solar_units = 1.0D20
      real(dp) :: largest
      integer :: j

      largest = max(s%luminosity_breakdown_lsun(i_lum_pp1), &
           s%luminosity_breakdown_lsun(i_lum_pp2),s%luminosity_breakdown_lsun(i_lum_pp3), &
           s%luminosity_breakdown_lsun(i_lum_cno),s%luminosity_breakdown_lsun(i_lum_3alpha), &
           abs(s%luminosity_breakdown_lsun(i_lum_neu)),s%luminosity_breakdown_lsun(i_lum_grav))
      if (largest > largest_solar_units) then
         do j = 1,7
            s%luminosity_breakdown_lsun(j) = s%luminosity_breakdown_lsun(j)/s%solar_luminosity_cgs
         end do
      end if
end subroutine luminosities_to_solar_units

! ---------------------------------------------------------------
! An error in the model file on unit iread: "the starting model <file>
! <what>", and ierr = 1.
subroutine report_model_error(iread, what, ierr)
      use io_units, only: report_error
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
      use io_units, only: check_io
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
      use io_units, only: report_error
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
      use io_units, only: report_error
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
