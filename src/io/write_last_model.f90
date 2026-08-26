!----------------------------------------------------------------------
! wrtlst
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original wrtlst.f; only variable names, source form, and comment
! style were updated.
!
! WRTLST writes the converged model to the "last model" file (stores
! last converged model) and "store" file (every npunch models). It
! decodes the physics flags used elsewhere as short character labels
! (ATM/EOS/ALOK/HIK) and delegates the actual write to putmodel2.
!
! MHP 10/02 NKK removed from declared variable list
subroutine write_last_model(iwrite, composition, log_density, log_luminosity, &
     log_pressure, log_radius, log_mass, log_temperature, convective_flag, &
     trial_log_temperature, trial_log_luminosity, fit_point_pressure, &
     fit_point_temperature, fit_point_radius, envelope_fit_coeffs, &
     trial_sign_flag, luminosity_breakdown, core_cz_top_index, &
     envelope_cz_bottom_index, model_number, num_shells, total_mass_msun, &
     log_teff, log_luminosity_lsun, log_total_mass, age_gyr, timestep_yr, &
     omega)

! WRTLST WRITES THE CONVERGED MODEL TO LAST MODEL A(STORES LAST
! CONVERGED MODEL) AND STORE MODELS D(EVERY NPUNCH MODELS)

!     WRITE MODEL OUT IN ASCII FORMAT
      use luout_lib
      use star_info_lib, only: star, json
      use opacity_table_lib
      use yale_eos_lib
      use scv_eos_lib
      implicit none
      integer, parameter :: nts=63, nps=76

      integer, intent(in) :: iwrite
      double precision, intent(in) :: composition(15,json), log_density(json), &
           log_luminosity(json), log_pressure(json), log_radius(json), &
           log_mass(json), log_temperature(json)
      logical, intent(in) :: convective_flag(json)
      double precision, intent(in) :: trial_log_temperature(3), &
           trial_log_luminosity(3), fit_point_pressure(3), &
           fit_point_temperature(3), fit_point_radius(3), &
           envelope_fit_coeffs(9), trial_sign_flag
      double precision, intent(inout) :: luminosity_breakdown(8)
      integer, intent(in) :: core_cz_top_index, envelope_cz_bottom_index, &
           model_number, num_shells
      double precision, intent(in) :: total_mass_msun, log_teff, &
           log_luminosity_lsun, log_total_mass, age_gyr, timestep_yr
      double precision, intent(in) :: omega(json)

      character(len=6) :: eos_flag
!     CHARACTER*4 ATM, LOK, HIK, COMPMIX
! MHP 4/25 changed LOK name to make it unique, used elsewhere
      character(len=4) :: atmosphere_flag, low_temp_opacity_flag, &
           high_temp_opacity_flag
! MHP 8/25 Removed unused variables
!      CHARACTER*256 FLAOL, FPUREZ
!      CHARACTER*256 FOPALE,FOPALE01,FOPALE06  ! FcondOpacP
! former common/i2o/: initial_composition_code now use-associated from
! run_diag_lib.
! physics flags:
! Determine atmosphere flag, ATM
      if (star%job%atm_choice .eq. 0) then
         atmosphere_flag='EDD '
      elseif (star%job%atm_choice .eq. 1) then
         atmosphere_flag='KS  '
      elseif (star%job%atm_choice .eq. 2) then
         atmosphere_flag='HRA '
      elseif (star%job%atm_choice .eq. 3) then
         atmosphere_flag='KUR '
      elseif (star%job%atm_choice .eq. 4) then
         atmosphere_flag='ALL '
! JNT 06/14
      elseif (star%job%atm_choice .eq. 5) then
         atmosphere_flag='K/C '
      endif
! Determine equation of state flag, EOS
      eos_flag='SAHA  '
      if (use_debye_huckel_correction) eos_flag='SAH+DH'
      if (use_scv_eos) then
         eos_flag='SCV   '
         if (star%ctrl%use_opal95_eos) eos_flag='SCV+OP'
         if (star%ctrl%use_opal2001_eos) eos_flag='SCV+O1'
         if (use_debye_huckel_correction) then
         if (star%ctrl%use_opal2006_eos) eos_flag='SCV+O6'
            eos_flag='SCV+DH'
            if (star%ctrl%use_opal95_eos) eos_flag='SCVDHO'
            if (star%ctrl%use_opal2001_eos) eos_flag='SCDHO1'
            if (star%ctrl%use_opal2006_eos) eos_flag='SCDHO6'
         endif
      else
         if (star%ctrl%use_opal95_eos) then
            eos_flag='OPAL  '
            if (use_debye_huckel_correction) eos_flag='OPA+DH'
         endif
         if (star%ctrl%use_opal2001_eos) then
            eos_flag='OPAL01'
            if (use_debye_huckel_correction) eos_flag='OP1+DH'
         endif
         if (star%ctrl%use_opal2006_eos) then
            eos_flag='OPAL06'
            if (use_debye_huckel_correction) eos_flag='OP6+DH'
         endif
      endif
! Determine low temperature opacities flag, ALOK
      low_temp_opacity_flag='NONE'
      if (star%ctrl%use_alex95_tables) low_temp_opacity_flag='ALEX'
      if (star%ctrl%use_kurucz90_tables) low_temp_opacity_flag='KURZ'
! Determine high temperature opacities flag, HIK
      high_temp_opacity_flag='NONE'
      if (star%ctrl%use_opal95_tables) high_temp_opacity_flag='OP95'
      if (star%ctrl%use_opal92_tables) high_temp_opacity_flag='OP92'
      if (star%ctrl%use_laol89_tables) high_temp_opacity_flag='LL89'

      call write_model2(log_luminosity_lsun,envelope_fit_coeffs,star%mixing_length_alpha, &
           age_gyr,timestep_yr,trial_sign_flag,composition,log_density, &
           log_luminosity,log_pressure,log_radius,log_mass,log_total_mass, &
           log_temperature,iwrite,short_file_unit,core_cz_top_index, &
           envelope_cz_bottom_index,convective_flag,star%job%use_extended_composition, &
           star%job%rotation_active,num_shells,model_number,omega,fit_point_pressure, &
           fit_point_radius,total_mass_msun,log_teff,luminosity_breakdown, &
           trial_log_luminosity,trial_log_temperature,fit_point_temperature, &
           atmosphere_flag,eos_flag,high_temp_opacity_flag, &
           star%job%diffuse_helium_active,star%job%use_diffusion_z,star%job%disk_locking_active, &
           star%job%instability_transport_active,star%job%ljdot0,low_temp_opacity_flag,star%job%lovstc, &
           star%job%envelope_overshoot_active,star%job%lovstm,use_pure_z_table,star%job%lsemic, &
           star%initial_composition_code,star%job%disk_omega_rad_s,star%job%disk_locking_age_gyr, &
           star%job%wind_saturation_omega)
! First three lines above are YREC7 inputs
! Last two lines are MODEL2 add-ons



      write(iowr,360) model_number,num_shells,log_teff,log_luminosity_lsun,age_gyr
  360 format(I6,'  #SHELLS=', I4, '  LogTeff=',F8.5, &
             '  Log(L/Lsun)=',F8.5,'  Age=',F12.5)
      if(iwrite.eq.11) then
       write(short_file_unit,330) model_number,iwrite
      else
       write(short_file_unit,340) model_number,age_gyr,iwrite
      endif
  330 format(' DUMPED MODEL',I5,'  FILE',I3)
  340 format(' DUMPED MODEL',I5,' AGE',F13.9,'  FILE',I3)
      if(iwrite.eq.ilast) rewind ilast
      return
end subroutine write_last_model
