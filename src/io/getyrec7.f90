!----------------------------------------------------------------------
! getyrec7
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original getyrec7.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! GETY7 - Read a YREC7 format stellar model file into memory
!
! llp  4/16/03
!
! Dummy-argument names for the fields shared with the MODEL2 format
! (essentially the whole signature) are chosen to match getmodel2.f90
! and putyrec7.f90, the sibling MODEL2 reader and YREC7 writer that
! share the same header/envelope/per-shell record layout.
subroutine getyrec7(log_luminosity_lsun, envelope_fit_coeffs, &
     mixing_length, age_gyr, timestep_yr, trial_sign_flag, composition, &
     log_density, log_luminosity, &
     log_pressure, log_radius, log_mass, log_total_mass, log_temperature, &
     iread, short_file_unit, core_cz_top_index, envelope_cz_bottom_index, &
     convective_flag, use_extended_composition, rotation_active, num_shells, &
     model_number, omega, fit_point_pressure, fit_point_radius, &
     total_mass_msun, log_teff, luminosity_breakdown, trial_log_luminosity, &
     trial_log_temperature, fit_point_temperature, &
     atm_code, eos_code, hik_code, use_diffusion_y, use_diffusion_z, &
     disk_locking_active, instability_transport_active, ljdot0, alok_code, &
     lovstc, envelope_overshoot_active, lovstm, use_pure_z_table, lsemic, &
     compmix_code, disk_pressure, disk_temperature, wind_saturation_omega, ierr)
! First three lines above are YREC7 inputs
! Last two lines are MODEL2 add-ons

! Several of this subroutine's own dummy arguments below (values read
! from the old model-file format -- rotation_active,
! envelope_overshoot_active, instability_transport_active, lovstc,
! lovstm, lsemic) happen to share names with unrelated const_lib
! runtime-config module variables, so `use, only:` the one member
! actually needed here rather than a blanket `use const_lib` (which
! would conflict with the dummy-argument declarations below).
      use const_lib, only: solar_luminosity_cgs
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(out) :: log_luminosity_lsun
      double precision, intent(out) :: envelope_fit_coeffs(9), mixing_length
      double precision, intent(out) :: age_gyr, timestep_yr
      double precision, intent(out) :: trial_sign_flag
      double precision, intent(out) :: composition(15,json), log_density(json), &
           log_luminosity(json), log_pressure(json), log_radius(json), &
           log_mass(json)
      double precision, intent(out) :: log_total_mass
      double precision, intent(out) :: log_temperature(json)
      integer, intent(in) :: iread, short_file_unit
      integer, intent(out) :: core_cz_top_index, envelope_cz_bottom_index
      logical, intent(out) :: convective_flag(json)
      logical, intent(out) :: use_extended_composition, rotation_active
      integer, intent(out) :: num_shells, model_number
      double precision, intent(out) :: omega(json)
      double precision, intent(out) :: fit_point_pressure(3), &
           fit_point_radius(3)
      double precision, intent(out) :: total_mass_msun, log_teff
      double precision, intent(out) :: luminosity_breakdown(8)
      double precision, intent(out) :: trial_log_luminosity(3), &
           trial_log_temperature(3), fit_point_temperature(3)
! atm_code/eos_code/hik_code/alok_code/compmix_code are short physics-
! option codes set to placeholder values here; their exact encoding is
! not interpreted in this file, so names are conservative. Naming
! matches getmodel2.f90.
      character*4, intent(out) :: atm_code
      character*6, intent(out) :: eos_code
      character*4, intent(out) :: hik_code
      logical, intent(out) :: use_diffusion_y, use_diffusion_z, &
           disk_locking_active, instability_transport_active, ljdot0
      character*4, intent(out) :: alok_code
      logical, intent(out) :: lovstc, envelope_overshoot_active, lovstm, &
           use_pure_z_table, lsemic
      character*4, intent(out) :: compmix_code
      double precision, intent(out) :: disk_pressure, disk_temperature, &
           wind_saturation_omega


      character*4 :: header_keyword

      double precision :: omega_log10(json)

      save

! --- locals ---
      integer :: i, j, ix, iz
      integer :: envelope_record_number
      integer :: pair_start_index, parity_test
      integer :: ii
      double precision :: max_luminosity_component

! Initialization

! Zero the composition arrays (composition) and the omega array.
! This is agreed upon protection in case they are inadvertently used in
! the program, but not present in the old YREC input data.
      integer, intent(out) :: ierr

      ierr = 0

      do i  = 1, json
         omega(i) = 0d0
       do j = 1, 15
          composition(j,i) = 0d0
       enddo
      enddo

! Set Model2-sspecific strings to "?"
! MHP 4/25 padded character strings to proper length
      atm_code = "  ? "
      eos_code = "    ? "
      hik_code = "  ? "
      alok_code = "  ? "
      compmix_code = "AG93"

! Set Model2-specific flags false
      use_diffusion_y = .false.
      use_diffusion_z = .false.
      disk_locking_active = .false.
      instability_transport_active = .false.
      ljdot0 = .false.
      lovstc = .false.
      envelope_overshoot_active = .false.
      lovstm = .false.
      use_pure_z_table = .false.
      lsemic = .false.

! Set Model2-specific constants to default
      disk_pressure = 7.2722D-6
      disk_temperature = 0d0
      wind_saturation_omega = 3d-4

! End of initialization

! Read in an input model from the LU iread file in YREC7 Format

      rewind(iread)

! get the header record
      read(iread,10) header_keyword,model_number,num_shells,total_mass_msun, &
           log_teff,log_luminosity_lsun,log_total_mass, &
           age_gyr,timestep_yr
 10   format(A4,1X,2I5,3F8.5,1X,F14.10,F13.9,F13.0)

! Get the flags describing the input physics and the information that
! is stored in the model (rotating or not?, extendec compotition, etc.)
      read(iread,30) rotation_active,use_extended_composition, &
           core_cz_top_index,envelope_cz_bottom_index,mixing_length
 30   format(2(L1,1X),5X,2I4,F7.3)

! get the LUMINOSITIES---PP(1-2-3)-CNO-HE-NU-GRAV
      read(iread,40) (luminosity_breakdown(j),j=1,7)
 40   format(10X,7E10.3)

! If TLUMX are in ergs, convert to solar units.  Decide by
! comparing to 10**20. If larger, divide by CLSUN.
      max_luminosity_component = dmax1(luminosity_breakdown(1), &
           luminosity_breakdown(2),luminosity_breakdown(3), &
           luminosity_breakdown(4),luminosity_breakdown(5), &
           dabs(luminosity_breakdown(6)),luminosity_breakdown(7))
      if(max_luminosity_component.gt.1.0D20) then
       do j = 1,7
          luminosity_breakdown(j) = luminosity_breakdown(j)/solar_luminosity_cgs
         enddo
      endif

! Get the ENVELOPE DATA
! FTRI is 1,normally.  It is set to -1  if any of the record numbers
! for the envelope triangle records was set to -1 by WRTLST.
      trial_sign_flag = 1D0
      do i = 1,3
       read(iread,70) envelope_record_number,trial_log_temperature(i), &
            trial_log_luminosity(i),fit_point_pressure(i), &
            fit_point_temperature(i),fit_point_radius(i), &
            (envelope_fit_coeffs(i+i+i-3+j),j = 1,3)
 70      format(3X,I2,F7.5,4F8.5,3E12.5)
       if(envelope_record_number.lt.0) trial_sign_flag = -1D0
 80   continue
      end do

! READ IN HENYEY POINTS IN 4 PARTS
!
! Read FIRST PART:M,R,L,P,T,RHO,CONV(T/F),X,Y,AND Z - one line per shell
      do i = 1,num_shells
       read(iread,100) log_mass(i),log_radius(i),log_luminosity(i), &
            log_pressure(i),log_temperature(i),log_density(i),convective_flag(i), &
            ix,iz
 100     format(0PF13.10,F10.7,1PE14.7,0PF10.7,2F10.7,L1,2I6)
       if(log_mass(i).lt.0D0.or.log_mass(i).gt.log_total_mass) then
          write(short_file_unit,1000)
          write(short_file_unit,1050) i
 1000       format(1X,39('>'),40('<')/1X,'RUN STOPPED')
 1050       format(' ERROR IN SUBROUTINE GETY7'/1X,'GLITCH IN SHELL', &
           I3,', SHELL MASS LESS THAN ZERO OR GREATER THAN STAR MASS')
          ! 2026 (phase five, step B): stop converted to ierr; run_yrec
          ! returns the error and the CLI wrapper (main) stops.
          ierr = 1
          return
       endif
       composition(1,i) = 1.0D-6*dfloat(ix)
       composition(3,i) = 1.0D-6*dfloat(iz)
 110  continue
      end do

! Read SECOND PART:ELEMENT ABUNDANCES: HE3,CNO CYCLE ELEMENTS.
! ABUNDANCES IN SURFACE AND CENTRAL CONVECTION ZONES STORED WITH 1
! VALUE PER ZONE RATHER THAN 1 VALUE PER SHELL
      if(envelope_cz_bottom_index.le.core_cz_top_index) then
       read(iread,200)(composition(i,core_cz_top_index),i=4,11)
       read(iread,200)(composition(i,envelope_cz_bottom_index),i=4,11)
      else
       do j = core_cz_top_index,envelope_cz_bottom_index
          read(iread,200) (composition(i,j),i = 4,11)
 200        format(8(1PE9.3,1X))
 210     continue
       end do
      endif
      if(core_cz_top_index.gt.1) then
! CONVECTIVE CORE- ASSIGN FIRST COMPOSITION VALUE TO SHELLS 1-JCORE
       do j = 1,core_cz_top_index-1
          do i = 4,11
             composition(i,j) = composition(i,core_cz_top_index)
 220        continue
          end do
 230     continue
       end do
      endif
      if(envelope_cz_bottom_index.lt.num_shells) then
! CONVECTIVE SURFACE- ASSIGN LAST COMPOSITION TO SHELLS JENV-M
       do j = envelope_cz_bottom_index+1,num_shells
          do i = 4,11
             composition(i,j) = composition(i,envelope_cz_bottom_index)
 240        continue
          end do
 250     continue
       end do
      endif
! DEFINE HE4 = 1 - X - Z - HE3.
      do i = 1,num_shells
       composition(2,i) = 1.0D0 - composition(1,i) - composition(3,i) - &
            composition(4,i)
 260  continue
      end do

! READ IN H2,LI6,LI7,BE9 (EXTENDED COMPOSITION)

      if(use_extended_composition) then
       if(envelope_cz_bottom_index.eq.1)then
! FULLY CONVECTIVE MODEL
          read(iread,300)(composition(i,1),i=12,15)
 300        format(4(1PE9.3,1X))
          do j = 1,num_shells
             do i = 12,15
              composition(i,j) = composition(i,1)
 305           continue
             end do
 310        continue
          end do
       else
! GENERAL CASE
! THESE ABUNDANCES ARE READ IN WITH 2 SHELLS PER LINE
          do pair_start_index = core_cz_top_index,envelope_cz_bottom_index-1,2
             read(iread,200)((composition(i,j),i = 12,15), &
                  j = pair_start_index,pair_start_index+1)
 315        continue
          end do
! IF AN ODD NUMBER OF ABUNDANCES STORED, READ IN LAST VALUE
          parity_test = envelope_cz_bottom_index-1 - core_cz_top_index
          if(mod(parity_test,2).ne.0) &
               read(iread,300)(composition(i,pair_start_index),i = 12,15)
          if(core_cz_top_index.gt.1) then
! CONVECTIVE CORE- ASSIGN FIRST VALUE TO SHELLS 1-JCORE
             do j = 1,core_cz_top_index-1
              do i = 12,15
                 composition(i,j) = composition(i,core_cz_top_index)
 320              continue
              end do
 330           continue
             end do
          endif
          if(envelope_cz_bottom_index.lt.num_shells) then
! CONVECTIVE SURFACE- ASSIGN LAST VALUE TO SHELLS JENV - M
             do j = envelope_cz_bottom_index+1,num_shells
              do i = 12,15
                 composition(i,j) = composition(i,envelope_cz_bottom_index)
 340              continue
              end do
 350           continue
             end do
          endif
       endif
      endif

! Read in FOURTH PART:  - LOG J/M STORED , 8 SHELLS PER LINE

      if(rotation_active) then
! READ OMEGA IN. If OMEGA records are missing go to OMEGA BYPASS below
      read(iread,500,end=9999)(omega_log10(ii),ii = 1,num_shells)
 500    format(0P8F10.7)
      do i = 1,num_shells
           if(omega_log10(i) .lt. 58.9D0) then
              omega(i) = 10D0**(-omega_log10(i))
           else
              omega(i) = 0D0
           endif
 510    continue
      end do
      endif

! KEEP iread OPEN
      rewind iread

      return

! OMEGA BYPASS
! Come here if OMEGA records are missing -- Can happen if
! a new omega file is being generated (LWNEW is true in STARIN)

9999  write(short_file_unit,9998) "GETYREC7: OMEGA records are missing from ", &
        " input model file - OMEGA array zeroed"
9998  format(2A)

      do i  = 1, json
         omega(i) = 0d0
      enddo

! KEEP iread OPEN
      rewind iread

      return

end subroutine getyrec7
