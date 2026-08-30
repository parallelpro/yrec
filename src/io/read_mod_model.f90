!----------------------------------------------------------------------
! read_mod_model
!----------------------------------------------------------------------
! Added 2026 (retire-legacy campaign): reader for the .mod model file
! (see io/write_mod_model.f90 for the format). Named global values --
! unknown names are reported and ignored, so the format can grow --
! then a named column-header line (verified against the expected
! layout), then one row per shell.
!
! Fills star% directly. The three values the loader's deck-
! consistency check and timestep seeding need are returned through
! arguments, matching what read_starting_model's dispatch does for
! the MODEL2/YREC7 readers. The luminosity breakdown is zeroed: the
! old formats' TLUMX record was only ever a cosmetic seed until the
! first energy-generation call recomputes it.
subroutine read_mod_model(iread, timestep_yr, mixing_length0, &
     use_extended_composition0, rotation_active0, ierr)

      use star_info_lib, only: star, json, i_h1, i_he4, i_metals, &
           i_he3, i_c12, i_c13, i_n14, i_n15, i_o16, i_o17, i_o18, &
           i_h2, i_li6, i_li7, i_be9
      use luout_lib
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
            read(line,*) name, star%dage
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

! TLUMX was only a cosmetic seed; engeb recomputes it on the first step.
      star%luminosity_breakdown = 0.0d0

      return
end subroutine read_mod_model
