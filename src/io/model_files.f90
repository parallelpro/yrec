!----------------------------------------------------------------------
! model_files
!----------------------------------------------------------------------
! The stellar model files. write_mod_model writes the .mod format: the
! final model, which is also the restart file and the solver's recovery
! source, and the Monte Carlo runs' models. read_mod_model reads it back.
! read_yrec7 and read_model2 read the two older starting-model formats,
! YREC7 and MODEL2 (llp 4/16/03), which share one layout: a header, the
! luminosities, the envelope records (skipped), then the shells in four
! parts (structure; CNO abundances; the extended composition; log J/M or
! omega). MODEL2 adds the last two lines of arguments. read_starting_model
! (core) picks the reader from the file's format tag. The readers fill
! star% directly and return the file's own header values; several of
! those share names with controls_buffer members, so this module never
! uses controls_buffer.
module model_files
      implicit none
      private
      public :: write_mod_model, read_mod_model, read_yrec7, read_model2

contains

! ---------------------------------------------------------------
! The .mod format holds exactly what a reload uses: the global scalars,
! the per-shell structure, composition and omega, and the two values the
! deck-consistency check needs (mixing length, extended-composition
! flag). Its layout is MESA-style: a version tag, named global values
! (one per line; readers ignore unknown names), a named column-header
! line, then one row per shell at full double precision, so a model
! reloaded between run cards round-trips exactly. Everything is read
! from star%; the unit is rewound first, since the file always holds one
! model, the latest.
subroutine write_mod_model(iwrite)

      use star_info_lib, only: star, i_h1, i_he4, i_metals, i_he3, &
           i_c12, i_c13, i_n14, i_n15, i_o16, i_o17, i_o18, i_h2, &
           i_li6, i_li7, i_be9
      use io_units
      use run_log, only: solver_diagnostics
      implicit none

      integer, intent(in) :: iwrite

      integer :: k
      integer, parameter :: slots(15) = [i_h1, i_he4, i_metals, i_he3, &
           i_c12, i_c13, i_n14, i_n15, i_o16, i_o17, i_o18, i_h2, &
           i_li6, i_li7, i_be9]

      rewind iwrite
      write(iwrite,'(a)') 'YMOD 1'
      write(iwrite,'(a,1x,i12)')        'model_number        ', star%model_number
      write(iwrite,'(a,1x,i12)')        'num_zones           ', star%nz
      write(iwrite,'(a,1x,es26.16e3)')  'star_mass_msun      ', star%star_mass
      write(iwrite,'(a,1x,es26.16e3)')  'log_total_mass      ', star%log_total_mass
      write(iwrite,'(a,1x,es26.16e3)')  'star_age_gyr        ', star%age_gyr
      write(iwrite,'(a,1x,es26.16e3)')  'timestep_yr         ', star%timestep_yr
      write(iwrite,'(a,1x,es26.16e3)')  'log_L               ', star%log_L
      write(iwrite,'(a,1x,es26.16e3)')  'log_Teff            ', star%log_Teff
      write(iwrite,'(a,1x,es26.16e3)')  'mixing_length_alpha ', star%mixing_length_alpha
      write(iwrite,'(a,1x,l1)')         'extended_composition', star%job%use_extended_composition
      write(iwrite,'(a,1x,l1)')         'rotation_active     ', star%job%rotation_active
      write(iwrite,'(a)') 'columns: zone log_mass logR luminosity logP logT' &
           // ' logRho conv omega h1 he4 z he3 c12 c13 n14 n15 o16 o17 o18' &
           // ' h2 li6 li7 be9'
      do k = 1, star%nz
         write(iwrite,'(i6,6(1x,es26.16e3),1x,l1,16(1x,es26.16e3))') k, &
              star%log_mass(k), star%logR(k), star%luminosity_lsun(k), &
              star%logP(k), star%logT(k), star%logRho(k), &
              star%convective_flag(k), star%omega(k), &
              star%xa(slots(1),k),  star%xa(slots(2),k),  star%xa(slots(3),k), &
              star%xa(slots(4),k),  star%xa(slots(5),k),  star%xa(slots(6),k), &
              star%xa(slots(7),k),  star%xa(slots(8),k),  star%xa(slots(9),k), &
              star%xa(slots(10),k), star%xa(slots(11),k), star%xa(slots(12),k), &
              star%xa(slots(13),k), star%xa(slots(14),k), star%xa(slots(15),k)
      end do
      flush(iwrite)

! 2026 log redesign: the per-model terminal echo is replaced by the
! run-log progress line (run_log); the "DUMPED MODEL" bookkeeping
! is solver forensics behind the diagnostics flag.
      if (solver_diagnostics()) then
         if (iwrite .eq. last_model_unit) then
          write(run_log_unit,330) star%model_number, iwrite
         else
          write(run_log_unit,340) star%model_number, star%age_gyr, iwrite
         endif
      end if
  330 format(' DUMPED MODEL',I5,'  FILE',I3)
  340 format(' DUMPED MODEL',I5,' AGE',F13.9,'  FILE',I3)

      return
end subroutine write_mod_model

! ---------------------------------------------------------------
! Named global values (unknown names are reported and ignored, so the
! format can grow), the column-header line (checked against the expected
! layout), then one row per shell into star%. The three values the
! deck-consistency check and the timestep seeding need are returned. The
! luminosity breakdown is zeroed; the first energy-generation call
! computes it.
subroutine read_mod_model(iread, timestep_yr, mixing_length0, &
     use_extended_composition0, rotation_active0, ierr)

      use star_info_lib, only: star, json, i_h1, i_he4, i_metals, &
           i_he3, i_c12, i_c13, i_n14, i_n15, i_o16, i_o17, i_o18, &
           i_h2, i_li6, i_li7, i_be9
      use io_units
      implicit none

      integer, intent(in) :: iread
      double precision, intent(out) :: timestep_yr, mixing_length0
      logical, intent(out) :: use_extended_composition0, rotation_active0
      integer, intent(out) :: ierr

      character(len=1024) :: line
      character(len=64) :: name
      integer :: k, izone, ios
      integer, parameter :: slots(15) = [i_h1, i_he4, i_metals, i_he3, &
           i_c12, i_c13, i_n14, i_n15, i_o16, i_o17, i_o18, i_h2, &
           i_li6, i_li7, i_be9]
      character(len=*), parameter :: expected_columns = &
           'columns: zone log_mass logR luminosity logP logT' &
           // ' logRho conv omega h1 he4 z he3 c12 c13 n14 n15 o16 o17 o18' &
           // ' h2 li6 li7 be9'

      ierr = 0
      timestep_yr = 0.0d0
      mixing_length0 = 0.0d0
      use_extended_composition0 = .false.
      rotation_active0 = .false.

      rewind iread
      read(iread,'(a)') line
      if (line(1:6) .ne. 'YMOD 1') then
         write(run_log_unit,'(a)') 'READ_MOD_MODEL: not a YMOD 1 file'
         ierr = 1
         return
      end if

! named globals until the column-header line
      do
         read(iread,'(a)',iostat=ios) line
         if (ios /= 0) then
            write(run_log_unit,'(a)') &
                 'READ_MOD_MODEL: file ended before columns line'
            ierr = 1
            return
         end if
         if (line(1:8) .eq. 'columns:') exit
         read(line,*,iostat=ios) name
         if (ios /= 0 .or. len_trim(line) == 0) cycle
         select case (trim(name))
         case ('model_number')
            read(line,*) name, star%model_number
         case ('num_zones')
            read(line,*) name, star%nz
         case ('star_mass_msun')
            read(line,*) name, star%star_mass
         case ('log_total_mass')
            read(line,*) name, star%log_total_mass
         case ('star_age_gyr')
            read(line,*) name, star%age_gyr
         case ('timestep_yr')
            read(line,*) name, timestep_yr
         case ('log_L')
            read(line,*) name, star%log_L
         case ('log_Teff')
            read(line,*) name, star%log_Teff
         case ('mixing_length_alpha')
            read(line,*) name, mixing_length0
         case ('extended_composition')
            read(line,*) name, use_extended_composition0
         case ('rotation_active')
            read(line,*) name, rotation_active0
         case default
            write(run_log_unit,'(2a)') &
                 'READ_MOD_MODEL: ignoring unknown global ', trim(name)
         end select
      end do

      if (trim(line) .ne. expected_columns) then
         write(run_log_unit,'(a)') &
              'READ_MOD_MODEL: unexpected column layout:'
         write(run_log_unit,'(a)') trim(line)
         ierr = 1
         return
      end if
      if (star%nz < 1 .or. star%nz > json) then
         write(run_log_unit,'(a,i9)') &
              'READ_MOD_MODEL: bad num_zones ', star%nz
         ierr = 1
         return
      end if

      do k = 1, star%nz
         read(iread,*,iostat=ios) izone, &
              star%log_mass(k), star%logR(k), star%luminosity_lsun(k), &
              star%logP(k), star%logT(k), star%logRho(k), &
              star%convective_flag(k), star%omega(k), &
              star%xa(slots(1),k),  star%xa(slots(2),k),  star%xa(slots(3),k), &
              star%xa(slots(4),k),  star%xa(slots(5),k),  star%xa(slots(6),k), &
              star%xa(slots(7),k),  star%xa(slots(8),k),  star%xa(slots(9),k), &
              star%xa(slots(10),k), star%xa(slots(11),k), star%xa(slots(12),k), &
              star%xa(slots(13),k), star%xa(slots(14),k), star%xa(slots(15),k)
         if (ios /= 0) then
            write(run_log_unit,'(a,i7)') &
                 'READ_MOD_MODEL: read error at shell ', k
            ierr = 1
            return
         end if
      end do

! TLUMX was only a cosmetic seed; compute_zone_energy recomputes it on the first step.
      star%luminosity_breakdown = 0.0d0

      return
end subroutine read_mod_model

! ---------------------------------------------------------------
! A YREC7-format model; the arguments are the file's header values.
subroutine read_yrec7(mixing_length, timestep_yr, triangle_orientation, &
     iread, core_cz_top_index, envelope_cz_bottom_index, &
     use_extended_composition, rotation_active, atm_code, eos_code, &
     hik_code, use_diffusion_y, use_diffusion_z, &
     disk_locking_active, instability_transport_active, &
     use_wind_torque, alok_code, core_overshoot_active, &
     envelope_overshoot_active, lovstm, use_pure_z_table, &
     use_semiconvection, disk_pressure, disk_temperature, &
     wind_saturation_omega, ierr)
      use io_units
      use star_info_lib, only: star, i_lum_3alpha, i_lum_cno, i_lum_grav, i_lum_neu, i_lum_pp1, i_lum_pp2, i_lum_pp3, json
      use math_lib
      implicit none

      double precision, intent(out) :: mixing_length
      double precision, intent(out) :: timestep_yr
      double precision, intent(out) :: triangle_orientation
      integer, intent(in) :: iread
      integer, intent(out) :: core_cz_top_index, envelope_cz_bottom_index
      logical, intent(out) :: use_extended_composition, rotation_active
! atm_code, eos_code, hik_code, alok_code and
! star%initial_composition_code: the format's physics-option codes,
! placeholders here (their encoding is not interpreted)
      character*4, intent(out) :: atm_code
      character*6, intent(out) :: eos_code
      character*4, intent(out) :: hik_code
      logical, intent(out) :: use_diffusion_y, use_diffusion_z, &
           disk_locking_active, instability_transport_active, use_wind_torque
      character*4, intent(out) :: alok_code
      logical, intent(out) :: core_overshoot_active, envelope_overshoot_active, lovstm, &
           use_pure_z_table, use_semiconvection
      double precision, intent(out) :: disk_pressure, disk_temperature, &
           wind_saturation_omega


      character*4 :: header_keyword

      double precision :: omega_log10(json)
! -log10(omega) at or above this value marks a non-rotating shell
! (omega stored as 0) in the legacy model file
      double precision, parameter :: omega_log10_zero_sentinel = 58.9D0
! --- locals ---
      integer :: i, j, ix, iz
      integer :: envelope_record_number
      integer :: pair_start_index, parity_test
      integer :: ii
      double precision :: max_luminosity_component

! Initialization

! Zero the star%xa arrays (star%xa) and the star%omega array.
! This is agreed upon protection in case they are inadvertently used in
! the program, but not present in the old YREC input data.
      integer, intent(out) :: ierr

      ierr = 0

      do i  = 1, json
         star%omega(i) = 0d0
       do j = 1, 15
          star%xa(j,i) = 0d0
       enddo
      enddo

! Set Model2-sspecific strings to "?"
! MHP 4/25 padded character strings to proper length
      atm_code = "  ? "
      eos_code = "    ? "
      hik_code = "  ? "
      alok_code = "  ? "
      star%initial_composition_code = "AG93"

! Set Model2-specific flags false
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

! Set Model2-specific constants to default
      disk_pressure = 7.2722D-6
      disk_temperature = 0d0
      wind_saturation_omega = 3d-4

! End of initialization

! Read in an input model from the LU iread file in YREC7 Format

      rewind(iread)

! get the header record
      read(iread,10) header_keyword,star%model_number,star%nz,star%star_mass, &
           star%log_Teff,star%log_L,star%log_total_mass, &
           star%age_gyr,timestep_yr
 10   format(A4,1X,2I5,3F8.5,1X,F14.10,F13.9,F13.0)

! Get the flags describing the input physics and the information that
! is stored in the model (rotating or not?, extendec compotition, etc.)
      read(iread,30) rotation_active,use_extended_composition, &
           core_cz_top_index,envelope_cz_bottom_index,mixing_length
 30   format(2(L1,1X),5X,2I4,F7.3)

! get the LUMINOSITIES---PP(1-2-3)-CNO-HE-NU-GRAV
      read(iread,40) (star%luminosity_breakdown(j),j=1,7)
 40   format(10X,7E10.3)

! If TLUMX are in ergs, convert to solar units.  Decide by
! comparing to 10**20. If larger, divide by CLSUN.
      max_luminosity_component = dmax1(star%luminosity_breakdown(i_lum_pp1), &
           star%luminosity_breakdown(i_lum_pp2),star%luminosity_breakdown(i_lum_pp3), &
           star%luminosity_breakdown(i_lum_cno),star%luminosity_breakdown(i_lum_3alpha), &
           dabs(star%luminosity_breakdown(i_lum_neu)),star%luminosity_breakdown(i_lum_grav))
      if(max_luminosity_component.gt.1.0D20) then
       do j = 1,7
          star%luminosity_breakdown(j) = star%luminosity_breakdown(j)/star%solar_luminosity_cgs
         enddo
      endif

! Get the ENVELOPE DATA
! FTRI is 1,normally.  It is set to -1  if any of the record numbers
! for the envelope triangle records was set to -1 by WRTLST.
      triangle_orientation = 1D0
      do i = 1,3
       read(iread,70) envelope_record_number,star%triangle_log_teff(i), &
            star%triangle_log_l(i),star%fit_point_pressure(i), &
            star%fit_point_temperature(i),star%fit_point_radius(i), &
            (star%envelope_fit_coeffs(3*i-3+j),j = 1,3)
 70      format(3X,I2,F7.5,4F8.5,3E12.5)
       if(envelope_record_number.lt.0) triangle_orientation = -1D0
      end do

! READ IN HENYEY POINTS IN 4 PARTS
!
! Read FIRST PART:M,R,L,P,T,RHO,CONV(T/F),X,Y,AND Z - one line per shell
      do i = 1,star%nz
       read(iread,100) star%log_mass(i),star%logR(i),star%luminosity_lsun(i), &
            star%logP(i),star%logT(i),star%logRho(i),star%convective_flag(i), &
            ix,iz
 100     format(0PF13.10,F10.7,1PE14.7,0PF10.7,2F10.7,L1,2I6)
       if(star%log_mass(i).lt.0D0.or.star%log_mass(i).gt.star%log_total_mass) then
          write(run_log_unit,1000)
          write(run_log_unit,1050) i
 1000       format(1X,39('>'),40('<')/1X,'RUN STOPPED')
 1050       format(' ERROR IN SUBROUTINE GETY7'/1X,'GLITCH IN SHELL', &
           I3,', SHELL MASS LESS THAN ZERO OR GREATER THAN STAR MASS')
          ! 2026 (phase five, step B): stop converted to ierr; run_yrec
          ! returns the error and the CLI wrapper (main) stops.
          ierr = 1
          return
       endif
       star%xa(1,i) = 1.0D-6*dfloat(ix)
       star%xa(3,i) = 1.0D-6*dfloat(iz)
      end do

! Read SECOND PART:ELEMENT ABUNDANCES: HE3,CNO CYCLE ELEMENTS.
! ABUNDANCES IN SURFACE AND CENTRAL CONVECTION ZONES STORED WITH 1
! VALUE PER ZONE RATHER THAN 1 VALUE PER SHELL
      if(envelope_cz_bottom_index.le.core_cz_top_index) then
       read(iread,200)(star%xa(i,core_cz_top_index),i=4,11)
       read(iread,200)(star%xa(i,envelope_cz_bottom_index),i=4,11)
      else
       do j = core_cz_top_index,envelope_cz_bottom_index
          read(iread,200) (star%xa(i,j),i = 4,11)
 200        format(8(1PE9.3,1X))
       end do
      endif
      if(core_cz_top_index.gt.1) then
! CONVECTIVE CORE- ASSIGN FIRST COMPOSITION VALUE TO SHELLS 1-JCORE
       do j = 1,core_cz_top_index-1
          do i = 4,11
             star%xa(i,j) = star%xa(i,core_cz_top_index)
          end do
       end do
      endif
      if(envelope_cz_bottom_index.lt.star%nz) then
! CONVECTIVE SURFACE- ASSIGN LAST COMPOSITION TO SHELLS JENV-M
       do j = envelope_cz_bottom_index+1,star%nz
          do i = 4,11
             star%xa(i,j) = star%xa(i,envelope_cz_bottom_index)
          end do
       end do
      endif
! DEFINE HE4 = 1 - X - Z - HE3.
      do i = 1,star%nz
       star%xa(2,i) = 1.0D0 - star%xa(1,i) - star%xa(3,i) - &
            star%xa(4,i)
      end do

! READ IN H2,LI6,LI7,BE9 (EXTENDED COMPOSITION)

      if(use_extended_composition) then
       if(envelope_cz_bottom_index.eq.1)then
! FULLY CONVECTIVE MODEL
          read(iread,300)(star%xa(i,1),i=12,15)
 300        format(4(1PE9.3,1X))
          do j = 1,star%nz
             do i = 12,15
              star%xa(i,j) = star%xa(i,1)
             end do
          end do
       else
! GENERAL CASE
! THESE ABUNDANCES ARE READ IN WITH 2 SHELLS PER LINE
          do pair_start_index = core_cz_top_index,envelope_cz_bottom_index-1,2
             read(iread,200)((star%xa(i,j),i = 12,15), &
                  j = pair_start_index,pair_start_index+1)
          end do
! IF AN ODD NUMBER OF ABUNDANCES STORED, READ IN LAST VALUE
          parity_test = envelope_cz_bottom_index-1 - core_cz_top_index
          if(mod(parity_test,2).ne.0) &
               read(iread,300)(star%xa(i,pair_start_index),i = 12,15)
          if(core_cz_top_index.gt.1) then
! CONVECTIVE CORE- ASSIGN FIRST VALUE TO SHELLS 1-JCORE
             do j = 1,core_cz_top_index-1
              do i = 12,15
                 star%xa(i,j) = star%xa(i,core_cz_top_index)
              end do
             end do
          endif
          if(envelope_cz_bottom_index.lt.star%nz) then
! CONVECTIVE SURFACE- ASSIGN LAST VALUE TO SHELLS JENV - M
             do j = envelope_cz_bottom_index+1,star%nz
              do i = 12,15
                 star%xa(i,j) = star%xa(i,envelope_cz_bottom_index)
              end do
             end do
          endif
       endif
      endif

! Read in FOURTH PART:  - LOG J/M STORED , 8 SHELLS PER LINE

      if(rotation_active) then
! READ OMEGA IN. If OMEGA records are missing go to OMEGA BYPASS below
      read(iread,500,end=9999)(omega_log10(ii),ii = 1,star%nz)
 500    format(0P8F10.7)
      do i = 1,star%nz
           if(omega_log10(i) .lt. omega_log10_zero_sentinel) then
              star%omega(i) = pow(10D0, (-omega_log10(i)))
           else
              star%omega(i) = 0D0
           endif
      end do
      endif

! KEEP iread OPEN
      rewind iread

      return

! OMEGA BYPASS
! Come here if OMEGA records are missing -- Can happen if
! a new star%omega file is being generated (LWNEW is true in STARIN)

9999  write(run_log_unit,9998) "GETYREC7: OMEGA records are missing from ", &
        " input model file - OMEGA array zeroed"
9998  format(2A)

      do i  = 1, json
         star%omega(i) = 0d0
      enddo

! KEEP iread OPEN
      rewind iread

      return

end subroutine read_yrec7

! ---------------------------------------------------------------
! A MODEL2-format model: the YREC7 layout plus the MODEL2 additions.
subroutine read_model2(mixing_length, timestep_yr, triangle_orientation, &
     iread, core_cz_top_index, envelope_cz_bottom_index, &
     use_extended_composition, rotation_active, atm_code, eos_code, &
     hik_code, use_diffusion_y, use_diffusion_z, &
     disk_locking_active, instability_transport_active, &
     use_wind_torque, alok_code, core_overshoot_active, &
     envelope_overshoot_active, lovstm, use_pure_z_table, &
     use_semiconvection, disk_pressure, disk_temperature, &
     wind_saturation_omega)
      use star_info_lib, only: star, i_lum_3alpha, i_lum_cno, i_lum_grav, i_lum_neu, i_lum_pp1, i_lum_pp2, i_lum_pp3, json
      implicit none

      double precision, intent(out) :: mixing_length
      double precision, intent(out) :: timestep_yr
      double precision, intent(out) :: triangle_orientation
      integer, intent(in) :: iread
      integer, intent(out) :: core_cz_top_index, envelope_cz_bottom_index
      logical, intent(out) :: use_extended_composition, rotation_active
! atm_code/eos_code/hik_code/alok_code/star%initial_composition_code are short physics-
! option codes read verbatim from the model header; their exact
! encoding is not interpreted in this file, so names are conservative.
      character*4, intent(out) :: atm_code
      character*6, intent(out) :: eos_code
      character*4, intent(out) :: hik_code
      logical, intent(out) :: use_diffusion_y, use_diffusion_z, &
           disk_locking_active, instability_transport_active, use_wind_torque
      character*4, intent(out) :: alok_code
      logical, intent(out) :: core_overshoot_active, envelope_overshoot_active, lovstm, &
           use_pure_z_table, use_semiconvection
      double precision, intent(out) :: disk_pressure, disk_temperature, &
           wind_saturation_omega


      character*5 :: header_keyword
      integer :: envelope_record_number
      integer :: i, j, k
      double precision :: max_luminosity_component

      rewind(iread)

! GET THE HEADER RECORD
      read(iread,10) header_keyword,star%model_number,star%nz,star%star_mass, &
           star%log_Teff,star%log_L,star%log_total_mass,star%age_gyr, &
           timestep_yr,star%log_mass(1),star%log_mass(star%nz)
 10   format(A5,2I8,5F16.12,1PE18.10,0P2F16.12)

! Get the flags describing the input physics and the information that
! is stored in the model (rotating or not?, extendec compotition, etc.)
      read(iread,30) core_cz_top_index,envelope_cz_bottom_index, &
           mixing_length,eos_code,atm_code,alok_code,hik_code, &
           use_pure_z_table,star%initial_composition_code,use_extended_composition, &
           use_diffusion_y,use_diffusion_z,use_semiconvection,core_overshoot_active, &
           envelope_overshoot_active,lovstm,rotation_active, &
           instability_transport_active,use_wind_torque,disk_locking_active, &
           disk_temperature,disk_pressure,wind_saturation_omega
   30 format(2I8,F16.10,1X,A6,1X,3(A4,1X),L1,1X,A4,1X,11(L1,1X), &
           3(1PE18.10))

! Get the  LUMINOSITIES---PP(1-2-3)-CNO-HE-NU-GRAV
! CONVERT TO SOLAR LUMINOSITIES
      read(iread,40) (star%luminosity_breakdown(j),j=1,7)
!   40 FORMAT('TLUMX',5X,1P7E17.9)
   40 format(10X,1P7E17.9)

! If TLUMX are in ergs, convert to solar units.  Decide by
! comparing to 10**20.  IF larger, divide by CLSUN.
      max_luminosity_component = dmax1(star%luminosity_breakdown(i_lum_pp1), &
           star%luminosity_breakdown(i_lum_pp2),star%luminosity_breakdown(i_lum_pp3), &
           star%luminosity_breakdown(i_lum_cno),star%luminosity_breakdown(i_lum_3alpha), &
           dabs(star%luminosity_breakdown(i_lum_neu)),star%luminosity_breakdown(i_lum_grav))
      if (max_luminosity_component.gt.1.0D20) then
       do j = 1,7
          star%luminosity_breakdown(j) = star%luminosity_breakdown(j)/star%solar_luminosity_cgs
         enddo
      endif

! Get the ENVELOPE DATA
! FTRI is 1,normally.  It is set to -1  if any of the record numbers
! for the envelope triangle records was set to -1 by WRTLST.
      triangle_orientation = 1D0
      do i = 1,3
         read(iread,50)envelope_record_number,star%triangle_log_teff(i), &
              star%triangle_log_l(i),star%fit_point_pressure(i), &
              star%fit_point_temperature(i),star%fit_point_radius(i), &
              (star%envelope_fit_coeffs(3*i-3+j),j=1,3)
         if (envelope_record_number .lt. 0) triangle_orientation = -1D0
      enddo
!   50 FORMAT('ENV',I2,5F16.12,1P3E20.12)
   50 format(4X,I2,5F16.12,1P3E20.12)

! Skip over column headings for per-shell data
      read(iread,*)
      read(iread,*)
      read(iread,*)

! Read per-shell data (HENYEY POINGS)
      do i = 1,star%nz
         read(iread,60) k,star%log_mass(i),star%logR(i),star%luminosity_lsun(i), &
              star%logP(i),star%logT(i),star%logRho(i), &
              star%omega(i),star%convective_flag(i),(star%xa(j,i),j=1,15)
       enddo
 60   format(I6,0P2F18.14,1PE24.16,0P3F18.14,1PE24.16,1x,L1, &
           3(0PF12.9),12(0PE16.8))

      return
end subroutine read_model2

end module model_files
