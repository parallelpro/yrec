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
subroutine wrtlst(iwrite, composition, log_density, log_luminosity, &
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
      use const_lib
      use luout_lib
      implicit none
      integer, parameter :: json = 5000
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
           high_temp_opacity_flag, initial_composition_code
! MHP 8/25 Removed unused variables
!      CHARACTER*256 FLAOL, FPUREZ
!      CHARACTER*256 FOPALE,FOPALE01,FOPALE06  ! FcondOpacP




! common/heflsh/: helium_flash_active (originally LKUTHE); not used in
! this file. Naming is local to this batch.
      logical :: helium_flash_active
      common/heflsh/ helium_flash_active


! llp  3/19/03 Add COMMON block /I2O/ for info directly transferred from
!      input to output model - starting with a code for th initial model
!      compostion (COMPMIX)
      common /i2o/ initial_composition_code

!      COMMON/CWIND/WMAX,EXMD,EXW,EXTAU,EXR,EXM,CONSTFACTOR,STRUCTFACTOR,LJDOT0
! MHP 8/17 ADDED EXCEN, C_2 TO COMMON BLOCK FOR MATT ET AL. 2012 CENT. TERM
! common/cwind/: only ljdot0 is used here. Naming matches dadcoeft.f90.
      double precision :: wind_saturation_omega, exmd, &
           wind_law_omega_exponent, extau, exr, exm, exl, expr, &
           constfactor, structfactor, excen, c_2
      logical :: ljdot0
      common/cwind/ wind_saturation_omega, exmd, wind_law_omega_exponent, &
           extau, exr, exm, exl, expr, constfactor, structfactor, excen, &
           c_2, ljdot0
! common/debhu/: only ldh is used here. Naming matches eqstat.f90.
      double precision :: cdh, etadh0, etadh1, zdh(18), xxdy, yydh, zzdh, &
           dhnue(18)
      logical :: ldh
      common/debhu/ cdh, etadh0, etadh1, zdh, xxdy, yydh, zzdh, dhnue, ldh
! OPACITY COMMON BLOCKS - modified 3/09
! common/newopac/: not used in this file; declared only to preserve
! layout. Naming matches getopac.f90.
      double precision :: laol_table_z1, laol_table_z2, opal_table_z1, &
           opal_table_z2, opal95_single_table_z, alex_table_z1, &
           kurucz_table_z1, kurucz_table_z2, molecular_opacity_logt_min, &
           molecular_opacity_logt_max
      logical :: use_alex06_tables, use_laol89_tables, use_opal92_tables, &
           use_opal95_tables, use_kurucz90_tables, use_alex95_tables, &
           use_two_z_tables
      common /newopac/ laol_table_z1, laol_table_z2, opal_table_z1, &
           opal_table_z2, opal95_single_table_z, alex_table_z1, &
           kurucz_table_z1, kurucz_table_z2, molecular_opacity_logt_min, &
           molecular_opacity_logt_max, use_alex06_tables, &
           use_laol89_tables, use_opal92_tables, use_opal95_tables, &
           use_kurucz90_tables, use_alex95_tables, use_two_z_tables
! MHP 8/25 Removed character file names from common block
! common/nwlaol/: only use_pure_z_table is used here. Naming matches
! getopac.f90.
      double precision :: olaol(12,104,52), oxa(12), ot(52), orho(104), &
           tollaol
      integer :: iolaol, numofxyz, numrho, numt, iopurez
      logical :: llaol, use_pure_z_table
      common/nwlaol/ olaol, oxa, ot, orho, tollaol, iolaol, numofxyz, &
           numrho, numt, llaol, use_pure_z_table, iopurez
! KC 2025-05-30 reordered common block elements
!       COMMON/OPALEOS/FOPALE,LOPALE,IOPALE,FOPALE01,LOPALE01,
!      x  FOPALE06,LOPALE06,LNumDeriv
! MHP 8/25 Remove file names from common blocks
! common/opaleos/: use_opal95_eos/use_opal2001_eos/use_opal2006_eos
! (LOPALE/lopale01/lopale06) are used here. Naming matches eqstat.f90.
      logical :: use_opal95_eos, use_opal2001_eos, use_opal2006_eos, &
           use_numerical_derivatives
      integer :: iopale
      common/opaleos/ use_opal95_eos, iopale, use_opal2001_eos, &
           use_opal2006_eos, use_numerical_derivatives
! common/scveos/: only use_scv_eos is used here. Naming matches
! eqscve.f90.
      double precision :: table_log10_temperature(nts), &
           hydrogen_table(nts,nps,12), helium_table(nts,nps,12), &
           entropy_of_mixing_table(nts,nps), metal_table(nts,nps,13), &
           envelope_table(nts,nps,12)
      integer :: num_pressure_points(nts)
      logical :: use_scv_eos
      integer :: scv_temp_index, scv_pressure_index
      common/scveos/ table_log10_temperature, hydrogen_table, &
           helium_table, entropy_of_mixing_table, metal_table, &
           envelope_table, num_pressure_points, use_scv_eos, &
           scv_temp_index, scv_pressure_index

      save

! physics flags:
! Determine atmosphere flag, ATM
      if (atm_choice .eq. 0) then
         atmosphere_flag='EDD '
      elseif (atm_choice .eq. 1) then
         atmosphere_flag='KS  '
      elseif (atm_choice .eq. 2) then
         atmosphere_flag='HRA '
      elseif (atm_choice .eq. 3) then
         atmosphere_flag='KUR '
      elseif (atm_choice .eq. 4) then
         atmosphere_flag='ALL '
! JNT 06/14
      elseif (atm_choice .eq. 5) then
         atmosphere_flag='K/C '
      endif
! Determine equation of state flag, EOS
      eos_flag='SAHA  '
      if (ldh) eos_flag='SAH+DH'
      if (use_scv_eos) then
         eos_flag='SCV   '
         if (use_opal95_eos) eos_flag='SCV+OP'
         if (use_opal2001_eos) eos_flag='SCV+O1'
         if (ldh) then
         if (use_opal2006_eos) eos_flag='SCV+O6'
            eos_flag='SCV+DH'
            if (use_opal95_eos) eos_flag='SCVDHO'
            if (use_opal2001_eos) eos_flag='SCDHO1'
            if (use_opal2006_eos) eos_flag='SCDHO6'
         endif
      else
         if (use_opal95_eos) then
            eos_flag='OPAL  '
            if (ldh) eos_flag='OPA+DH'
         endif
         if (use_opal2001_eos) then
            eos_flag='OPAL01'
            if (ldh) eos_flag='OP1+DH'
         endif
         if (use_opal2006_eos) then
            eos_flag='OPAL06'
            if (ldh) eos_flag='OP6+DH'
         endif
      endif
! Determine low temperature opacities flag, ALOK
      low_temp_opacity_flag='NONE'
      if (use_alex95_tables) low_temp_opacity_flag='ALEX'
      if (use_kurucz90_tables) low_temp_opacity_flag='KURZ'
! Determine high temperature opacities flag, HIK
      high_temp_opacity_flag='NONE'
      if (use_opal95_tables) high_temp_opacity_flag='OP95'
      if (use_opal92_tables) high_temp_opacity_flag='OP92'
      if (use_laol89_tables) high_temp_opacity_flag='LL89'

      call putmodel2(log_luminosity_lsun,envelope_fit_coeffs,cmixl, &
           age_gyr,timestep_yr,trial_sign_flag,composition,log_density, &
           log_luminosity,log_pressure,log_radius,log_mass,log_total_mass, &
           log_temperature,iwrite,short_file_unit,core_cz_top_index, &
           envelope_cz_bottom_index,convective_flag,use_extended_composition, &
           rotation_active,num_shells,model_number,omega,fit_point_pressure, &
           fit_point_radius,total_mass_msun,log_teff,luminosity_breakdown, &
           trial_log_luminosity,trial_log_temperature,fit_point_temperature, &
           atmosphere_flag,eos_flag,high_temp_opacity_flag, &
           diffuse_helium_active,use_diffusion_z,disk_locking_active, &
           instability_transport_active,ljdot0,low_temp_opacity_flag,lovstc, &
           envelope_overshoot_active,lovstm,use_pure_z_table,lsemic, &
           initial_composition_code,disk_pressure,disk_temperature, &
           wind_saturation_omega)
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
end subroutine wrtlst
