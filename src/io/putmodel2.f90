!----------------------------------------------------------------------
! putmodel2
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original putmodel2.f; only variable names, source form, and comment
! style were updated.
!
! PUTMMODEL2 - Write out a stellar model in MODEL2 format.
!
! llp 4/16/03
subroutine putmodel2(log_luminosity_lsun, envelope_fit_coeffs, mixing_length, &
     age_gyr, timestep_yr, trial_sign_flag, composition, log_density, &
     log_luminosity, log_pressure, log_radius, log_mass, log_total_mass, &
     log_temperature, iwrite, ishort, core_cz_top_index, envelope_cz_bottom_index, &
     convective_flag, use_extended_composition, rotation_active, num_shells, &
     model_number, omega, fit_point_pressure, fit_point_radius, &
     total_mass_msun, log_teff, luminosity_breakdown, trial_log_luminosity, &
     trial_log_temperature, fit_point_temperature, &
     atmosphere_flag, eos_flag, high_temp_opacity_flag, diffuse_helium_active, &
     use_diffusion_z, disk_locking_active, instability_transport_active, &
     ljdot0, low_temp_opacity_flag, lovstc, envelope_overshoot_active, lovstm, &
     use_pure_z_table, lsemic, initial_composition_code, disk_pressure, &
     disk_temperature, wind_saturation_omega)
! First three lines above are YREC7 inputs
! Last two lines are MODEL2 add-ons

!  Write output model in MODEL2 format

! Several of this subroutine's own dummy arguments below (rotation_active
! etc, values being written out to the model file) happen to share
! names with unrelated const_lib runtime-config module variables, so
! `use, only:` the one member actually needed here rather than a
! blanket `use const_lib`. Same treatment as io/getyrec7.f90.
      use const_lib, only: solar_luminosity_cgs
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: log_luminosity_lsun
      double precision, intent(in) :: envelope_fit_coeffs(9), mixing_length
      double precision, intent(in) :: age_gyr, timestep_yr, trial_sign_flag
      double precision, intent(in) :: composition(15,json), log_density(json), &
           log_luminosity(json), log_pressure(json), log_radius(json), &
           log_mass(json)
      double precision, intent(in) :: log_total_mass
      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: iwrite, ishort, core_cz_top_index, &
           envelope_cz_bottom_index
      logical, intent(in) :: convective_flag(json), use_extended_composition, &
           rotation_active
      integer, intent(in) :: num_shells, model_number
      double precision, intent(in) :: omega(json)
      double precision, intent(in) :: fit_point_pressure(3), fit_point_radius(3)
      double precision, intent(in) :: total_mass_msun, log_teff
      double precision, intent(inout) :: luminosity_breakdown(8)
      double precision, intent(in) :: trial_log_luminosity(3), &
           trial_log_temperature(3), fit_point_temperature(3)

      character(len=4), intent(in) :: atmosphere_flag
      character(len=6), intent(in) :: eos_flag
      character(len=4), intent(in) :: high_temp_opacity_flag
      logical, intent(in) :: diffuse_helium_active, use_diffusion_z, &
           disk_locking_active, instability_transport_active, ljdot0
      character(len=4), intent(in) :: low_temp_opacity_flag
      logical, intent(in) :: lovstc, envelope_overshoot_active, lovstm, &
           use_pure_z_table, lsemic
      character(len=4), intent(in) :: initial_composition_code
      double precision, intent(in) :: disk_pressure, disk_temperature, &
           wind_saturation_omega
! --- locals ---
      integer :: i, j, zone_sign_index
      double precision :: max_luminosity_component

! PUTYMODEL2 writes the model provided to it out to the LU IWRITE file.
! Any error messages go to LU ISHORT.

      if(atmosphere_flag .eq. ' ? ') then
         write(ishort,7)
  7      format('*** YREC7 input file, flags, etc., have been ', &
                'defaulted.  ***')
      endif

! write header records
      if(age_gyr .lt. 1d3) then
         write(iwrite,10) 'MOD2 ',model_number,num_shells,total_mass_msun, &
              log_teff,log_luminosity_lsun,log_total_mass,age_gyr, &
              timestep_yr,log_mass(1),log_mass(num_shells)
 10      format(A5,2I8,5F16.11,1PE18.10,0P2F16.12)
      else if (age_gyr .lt. 1D4) then
         write(iwrite,11) 'MOD2 ',model_number,num_shells,total_mass_msun, &
              log_teff,log_luminosity_lsun,log_total_mass,age_gyr, &
              timestep_yr,log_mass(1),log_mass(num_shells)
 11      format(A5,2I8,4F16.12,F16.10,1PE18.10,0P2F16.12)
      else if (age_gyr .lt. 1D5) then
         write(iwrite,12) 'MOD2 ',model_number,num_shells,total_mass_msun, &
              log_teff,log_luminosity_lsun,log_total_mass,age_gyr, &
              timestep_yr,log_mass(1),log_mass(num_shells)
 12      format(A5,2I8,4F16.12,F16.9,1PE18.10,0P2F16.12)
      else
         write(iwrite,13) 'MOD2 ',model_number,num_shells,total_mass_msun, &
              log_teff,log_luminosity_lsun,log_total_mass,age_gyr, &
              timestep_yr,log_mass(1),log_mass(num_shells)
 13      format(A5,2I8,4F16.12,F16.8,1PE18.10,0P2F16.12)
      endif

! write physics flags:
      write(iwrite,30) core_cz_top_index,envelope_cz_bottom_index,mixing_length, &
           eos_flag,atmosphere_flag,low_temp_opacity_flag,high_temp_opacity_flag, &
           use_pure_z_table,initial_composition_code,use_extended_composition, &
           diffuse_helium_active,use_diffusion_z,lsemic,lovstc, &
           envelope_overshoot_active,lovstm,rotation_active, &
           instability_transport_active,ljdot0,disk_locking_active, &
           disk_temperature,disk_pressure,wind_saturation_omega
   30 format(2I8,F16.10,1X,A6,1X,3(A4,1X),L1,1X,A4,1X,11(L1,1X), &
           3(1PE18.10))

! write luminosities
! If TLUMX are in solar units, convert to ergs.  Decide by
! comparing to 10**20, if smaller, multiply by CLSUN

      max_luminosity_component = dmax1(luminosity_breakdown(1), &
           luminosity_breakdown(2),luminosity_breakdown(3), &
           luminosity_breakdown(4),luminosity_breakdown(5), &
           dabs(luminosity_breakdown(6)),luminosity_breakdown(7))
      if(max_luminosity_component.le.1.0D20) then
       do j = 1,7
          luminosity_breakdown(j) = luminosity_breakdown(j) * solar_luminosity_cgs
         enddo
      endif
      write(iwrite,40) (luminosity_breakdown(j),j=1,7)
   40 format('TLUMX',5X,1P7E17.9)

! write ENVELOPE DATA
      do i = 1,3
         zone_sign_index = i
         if(trial_sign_flag.lt.0D0) zone_sign_index = -i
         write(iwrite,50)zone_sign_index,trial_log_temperature(i), &
              trial_log_luminosity(i),fit_point_pressure(i), &
              fit_point_temperature(i),fit_point_radius(i), &
              (envelope_fit_coeffs(3*i-3+j),j=1,3)
      enddo
   50 format('ENV',I2,5F16.12,1P3E20.12)

! write column headings for all per shell information
      write(iwrite,*)
      write(iwrite,55) ' SHELL','MASS       ','RADIUS      ', &
        'LUMINOSITY       ','PRESSURE     ','TEMPERATURE   ', &
        'DENSITY     ','OMEGA         ',' C','  H1     ','   He4    ', &
        '   METALS  ','He3      ','C12      ','C13      ','N14      ', &
        'N15      ','O16      ','O17      ','O18      ','H2       ', &
        'Li6      ','Li7      ','Be9      '
 55   format(A6,2A18,A24,3A18,A24,A2,3A12,12A16)
      write(iwrite,*)

! LLP 4/16/03  Output OMEGA rather Log10(OMEGA), a change from YREC7 files

! write per shell information, one line per shell
      do i = 1,num_shells
         write(iwrite,60) i,log_mass(i),log_radius(i),log_luminosity(i), &
              log_pressure(i),log_temperature(i),log_density(i), &
              omega(i),convective_flag(i),(composition(j,i),j=1,15)
       enddo
 60   format(I6,0P2F18.14,1PE24.16,0P3F18.14,1PE24.16,1x,L1, &
           3(0PF12.9),12(0PE16.8))

      return
end subroutine putmodel2
