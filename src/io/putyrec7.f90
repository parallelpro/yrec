!----------------------------------------------------------------------
! putyrec7
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original putyrec7.f; only variable names, source form, and comment
! style were updated.
!
! PUTYREC7 - Write out a stellar model in YREC7 format.
!
! llp 4/16/03
subroutine putyrec7(log_luminosity_lsun, envelope_fit_coeffs, mixing_length, &
     age_gyr, timestep_yr, trial_sign_flag, composition, log_density, &
     log_luminosity, log_pressure, log_radius, log_mass, log_total_mass, &
     log_temperature, iwrite, core_cz_top_index, envelope_cz_bottom_index, &
     convective_flag, use_extended_composition, rotation_active, num_shells, &
     model_number, omega, fit_point_pressure, fit_point_radius, &
     total_mass_msun, log_teff, luminosity_breakdown, trial_log_luminosity, &
     trial_log_temperature, fit_point_temperature)
      use star_info_lib, only: star, i_lum_3alpha, i_lum_cno, i_lum_grav, i_lum_neu, i_lum_pp1, i_lum_pp2, i_lum_pp3, json
!      & ATM,EOS,HIK,LDIFY,LDIFZ,LDISK,LINSTB,LJDOT0,ALOK,
!      & LOVSTC,LOVSTE,LOVSTM,LPUREZ,LSEMIC,COMPMIX,PDISK,TDISK,WMAX)  ! KC 2025-05-31
! First three lines above are YREC7 inputs
! Last two lines are MODEL2 add-ons

!     WRITE MODEL OUT IN ASCII FORMAT
! Several of this subroutine's own dummy arguments below (rotation_active
! etc, values being written out to the model file) happen to share
! names with unrelated const_lib runtime-config module variables, so
! `use, only:` the one member actually needed here rather than a
! blanket `use const_lib`. Same treatment as io/getyrec7.f90.
! (solar_luminosity_cgs now comes from star% -- 2026 phase-A
! eviction; the former `use const_lib, only:` import is gone.)
      implicit none

      double precision, intent(in) :: log_luminosity_lsun
      double precision, intent(in) :: envelope_fit_coeffs(9), mixing_length
      double precision, intent(in) :: age_gyr, timestep_yr, trial_sign_flag
      double precision, intent(in) :: composition(15,json), log_density(json), &
           log_luminosity(json), log_pressure(json), log_radius(json), &
           log_mass(json)
      double precision, intent(in) :: log_total_mass
      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: iwrite, core_cz_top_index, envelope_cz_bottom_index
      logical, intent(in) :: convective_flag(json), use_extended_composition, &
           rotation_active
      integer, intent(in) :: num_shells, model_number
! omega is clamped to a floor value in place below, so it must be
! inout even though the caller does not otherwise rely on that.
      double precision, intent(inout) :: omega(json)
      double precision, intent(in) :: fit_point_pressure(3), fit_point_radius(3)
      double precision, intent(in) :: total_mass_msun, log_teff
      double precision, intent(inout) :: luminosity_breakdown(8)
      double precision, intent(in) :: trial_log_luminosity(3), &
           trial_log_temperature(3), fit_point_temperature(3)


! CHARACTER*4 ATM, LOK, HIK, COMPMIX
! MHP 4/25 changed LOK name to make it unique, used elsewhere
!       CHARACTER*4 ATM,HIK,ALOK,COMPMIX
!       CHARACTER*6 EOS

      double precision :: omega_log10(json)
! --- locals ---
      integer :: i, j, zone_sign_index, ix, iz, num_full_rows, &
           remainder_count, start_index
      integer :: pair_start_index, shell_index, species_index, parity_test
      double precision :: max_luminosity_component

! PUTYREC7 writes the model provided to it out to the LU IWRITE file.
! Any error messages go to LU ISHORT.


!  THE 6 HEADER RECORDS WRITTEN FIRST

      if (timestep_yr .ge. 1D6) then
       write(iwrite,10) model_number,num_shells,total_mass_msun,log_teff, &
            log_luminosity_lsun,log_total_mass,age_gyr,timestep_yr
   10    format('NMOD ',2I5,3F8.5,1X,F14.10,F13.9,F13.0)
      else
       write(iwrite,20) model_number,num_shells,total_mass_msun,log_teff, &
            log_luminosity_lsun,log_total_mass,age_gyr,timestep_yr
   20    format('NMOD ',2I5,3F8.5,1X,F14.10,F13.9,1PE13.6)
      endif

! write PHYSICS FLAGS- ROTATION,EXTENDED COMP,CENTRAL AND SURFACE CZ'S
      write(iwrite,30)rotation_active,use_extended_composition, &
           core_cz_top_index,envelope_cz_bottom_index,mixing_length
   30 format(2(L1,1X),'F F ',2I4,F7.3)

! write LUMINOSITIES
! If TLUMX are in solar units, convert to ergs.  Decide by
! comparing to 10**20.  If smaller, multiply by CLSUN.
      max_luminosity_component = dmax1(luminosity_breakdown(i_lum_pp1), &
           luminosity_breakdown(i_lum_pp2),luminosity_breakdown(i_lum_pp3), &
           luminosity_breakdown(i_lum_cno),luminosity_breakdown(i_lum_3alpha), &
           dabs(luminosity_breakdown(i_lum_neu)),luminosity_breakdown(i_lum_grav))
      if(max_luminosity_component.le.1.0D20) then
       do j = 1,7
          luminosity_breakdown(j) = luminosity_breakdown(j) * star%solar_luminosity_cgs
         enddo
      endif
      write(iwrite,40) (luminosity_breakdown(j),j=1,7)
   40 format('TLUMX',5X,1P7E10.3)

! write ENVELOPE DATA
      do i = 1,3
      zone_sign_index = i
      if(trial_sign_flag.lt.0D0) zone_sign_index = -i
      write(iwrite,50)zone_sign_index,trial_log_temperature(i), &
           trial_log_luminosity(i),fit_point_pressure(i), &
           fit_point_temperature(i),fit_point_radius(i), &
           (envelope_fit_coeffs(3*i-3+j),j=1,3)
   50 format('ENV',I2,F7.5,4F8.5,1P3E12.5)
      end do

! write HENYEY POINTS, one line per shell
      do i = 1,num_shells
       ix = idint(1.0D6*composition(1,i) + 0.50D0)
       iz = idint(1.0D6*composition(3,i) + 0.50D0)
       write(iwrite,100) log_mass(i),log_radius(i),log_luminosity(i), &
            log_pressure(i),log_temperature(i),log_density(i), &
            convective_flag(i),ix,iz
  100    format(0PF13.10,F10.7,1PE14.7,0PF10.7,2F10.7,L1,2I6)
      end do

! WRITE OUT COMPOSITION ARRAY - CENTRAL AND SURFACE
! CONVECTION ZONES HAVE ONLY 1 COMPOSITION STORED PER ZONE
! IF STAR IS FULLY CONVECTIVE, WRITE OUT 1 COMP.
      if(envelope_cz_bottom_index.le.core_cz_top_index) then
       write(iwrite,200)(composition(i,core_cz_top_index),i=4,11)
       write(iwrite,200)(composition(i,envelope_cz_bottom_index),i=4,11)
      else
       do j = core_cz_top_index,envelope_cz_bottom_index
          write(iwrite,200) (composition(i,j),i = 4,11)
  200       format(8(1PE9.3,1X))
       end do
      endif
! EXTENDED COMP - WRITE OUT 1 ABUND IF FULLY CONVECTIVE,
! OTHERWISE WRITE OUT 2 POINTS PER LINE.
      if(use_extended_composition) then
       if(envelope_cz_bottom_index.le.core_cz_top_index) then
          write(iwrite,220)(composition(i,envelope_cz_bottom_index),i=12,15)
  220       format(4(1PE9.3,1X))
       else
          do pair_start_index = core_cz_top_index,envelope_cz_bottom_index-1,2
             write(iwrite,200)((composition(species_index,shell_index), &
                  species_index = 12,15),shell_index = pair_start_index,pair_start_index+1)
          end do
!   IF AN ODD NUMBER OF ABUNDANCES EXISTS, WRITE IN LAST VALUE
          parity_test = envelope_cz_bottom_index-1 - core_cz_top_index
          if(mod(parity_test,2).ne.0) &
               write(iwrite,220)(composition(species_index,pair_start_index), &
               species_index = 12,15)
       endif
      endif
! WRITE OUT RUN OF OMEGA,STORED 8 ELEMENTS PER LINE
! (- THE LOG OF OMEGA IS STORED)
      if(rotation_active) then
       do i = 1,num_shells
          if (omega(i).le.1.D-59) omega(i)=1.D-59
          omega_log10(i) = dabs(dlog10(omega(i)))
       end do
       num_full_rows = int(num_shells/8)
       remainder_count = num_shells - num_full_rows*8
       start_index = 1
       do i = 1,num_full_rows
          write(iwrite,310)(omega_log10(j),j = start_index,start_index + 7)
  310       format(0P8F10.7)
          start_index = start_index + 8
       end do
       if(remainder_count.gt.0) write(iwrite,310) &
            (omega_log10(j),j=start_index,start_index+remainder_count-1)
      endif

      return
end subroutine putyrec7
