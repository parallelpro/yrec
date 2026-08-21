!----------------------------------------------------------------------
! putstore
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original putstore.f; only variable names, source form, and comment
! style were updated.
!
! PUTSTORE - Write out a verbose stellar model to the .store file. This
!            program is mostly a clone of wrtlst.f and putmodel2.f,
!            combined.
!
! G Somers  11/14
!
! 06
subroutine putstore(composition, log_density, log_luminosity, log_pressure, &
     log_radius, log_mass, log_temperature, convective_flag, &
     trial_log_temperature, trial_log_luminosity, fit_point_pressure, &
     fit_point_temperature, fit_point_radius, envelope_fit_coeffs, &
     trial_sign_flag, luminosity_breakdown, core_cz_top_index, &
     envelope_cz_bottom_index, model_number, num_shells, total_mass_msun, &
     log_teff, log_luminosity_lsun, log_total_mass, age_gyr, timestep_yr, &
     omega, mass_coordinate, rotation_eta2, radius_ratio_r0, shape_factor_fp, &
     shape_factor_ft, specific_angular_momentum, shell_moment_of_inertia)

! PUTSTORE PUTS THE MOST RECENT VERBOSE OUTPUT FILE INTO THE .STORE FILE,
! EITHER AT SPECIFIED AGES, EVERY NPRTMOD MODELS, OR AT THE END OF RUNS.

!     WRITE MODEL OUT IN ASCII FORMAT
      use scrtch_lib
      use luout_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000
      integer, parameter :: nts=63, nps=76

      double precision, intent(inout) :: composition(15,json), log_density(json), &
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
      double precision, intent(inout) :: omega(json), mass_coordinate(json), &
           rotation_eta2(json), radius_ratio_r0(json), shape_factor_fp(json), &
           shape_factor_ft(json), specific_angular_momentum(json), &
           shell_moment_of_inertia(json)

      character(len=6) :: eos_flag
!     CHARACTER*4 ATM, LOK, HIK, COMPMIX
! MHP 4/25 changed LOK name to make it unique, used elsewhere
      character(len=4) :: atmosphere_flag, low_temp_opacity_flag, &
           high_temp_opacity_flag, initial_composition_code
! MHP 8/25 Removed unused variables
!      CHARACTER*256 FLAOL, FPUREZ
!      CHARACTER*256 FOPALE,FOPALE01,FOPALE06  ! FcondOpacP




! common/heflsh/: not used in this file. Naming is local to this batch.
      logical :: helium_flash_active
      common/heflsh/ helium_flash_active


! llp  3/19/03 Add COMMON block /I2O/ for info directly transferred from
!      input to output model - starting with a code for th initial model
!      compostion (COMPMIX)
      common /i2o/ initial_composition_code

! MHP 8/17 ADDED EXCEN, C_2 TO COMMON BLOCK FOR MATT ET AL. 2012 CENT. TERM
! common/cwind/: only ljdot0 is used here. Naming matches dadcoeft.f90.
      double precision :: wind_saturation_omega, exmd, &
           wind_law_omega_exponent, extau, exr, exm, exl, expr, &
           constfactor, structfactor, excen, c_2
      logical :: ljdot0
      common/cwind/ wind_saturation_omega, exmd, wind_law_omega_exponent, &
           extau, exr, exm, exl, expr, constfactor, structfactor, excen, &
           c_2, ljdot0
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
! are used here. Naming matches eqstat.f90.
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

! common/sound/: adiabatic_index_gamma1 is used here. Naming is local
! to this batch.
      double precision :: adiabatic_index_gamma1(json)
      logical :: sound_speed_output_active
      common/sound/ adiabatic_index_gamma1, sound_speed_output_active
! common/temp2/: only es_circulation_velocity/gsf_circulation_velocity/
! secular_shear_velocity are used here. Naming matches vcirc.f90.
      double precision :: es_circulation_velocity(json), &
           es_circulation_velocity_prev(json), secular_shear_velocity(json), &
           secular_shear_velocity_prev(json), hle(json), &
           gsf_circulation_velocity(json), gsf_circulation_velocity_prev(json), &
           mu_gradient_velocity(json)
      common/temp2/ es_circulation_velocity, es_circulation_velocity_prev, &
           secular_shear_velocity, secular_shear_velocity_prev, hle, &
           gsf_circulation_velocity, gsf_circulation_velocity_prev, &
           mu_gradient_velocity
! common/roten/: rotational_energy_term, used here. Naming is local to
! this batch.
      double precision :: rotational_energy_term(json)
      common/roten/ rotational_energy_term

      save

! --- locals ---
      integer :: i, j, k, ii, zone_sign_index, idm
      integer :: id(json)
      double precision :: max_luminosity_component
      double precision :: local_g_const, sg, fm, duma, oblateness_a, &
           pole_to_equator_ratio, vtot
! JXBEG, JXEND and LSHELL are not dummy arguments of this subroutine
! and are never assigned anywhere in it -- the CALL PINDEX below reads
! them uninitialized. This looks like a latent bug from an earlier
! refactor (these were presumably meant to come in as arguments, as
! they do in wrtout.f/putyrec7.f), but it is preserved exactly as in
! the original putstore.f; not "fixed" here.
      integer :: jxbeg, jxend
      logical :: lshell
! LMILNE is likewise never declared in any COMMON block of the
! original putstore.f (unlike wrtout.f/wrtlst.f, which get it from
! common/ccout2/) -- it is an implicitly-typed, SAVE'd, always-default
! (.FALSE. in practice) local here, so "IF(LMILNE) CALL WRTMIL(...)"
! below is effectively dead code. Preserved exactly as in the
! original; not "fixed" here. Renamed to lmilne_local (2026) since it
! would otherwise collide with the unrelated const_lib lmilne added
! for former common/ccout2/.
      logical :: lmilne_local

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
      endif
! Determine equation of state flag, EOS
      eos_flag='SAHA  '
      if (use_debye_huckel_correction) eos_flag='SAH+DH'
      if (use_scv_eos) then
         eos_flag='SCV   '
         if (use_opal95_eos) eos_flag='SCV+OP'
         if (use_opal2001_eos) eos_flag='SCV+O1'
         if (use_debye_huckel_correction) then
         if (use_opal2006_eos) eos_flag='SCV+O6'
            eos_flag='SCV+DH'
            if (use_opal95_eos) eos_flag='SCVDHO'
            if (use_opal2001_eos) eos_flag='SCDHO1'
            if (use_opal2006_eos) eos_flag='SCDHO6'
         endif
      else
         if (use_opal95_eos) then
            eos_flag='OPAL  '
            if (use_debye_huckel_correction) eos_flag='OPA+DH'
         endif
         if (use_opal2001_eos) then
            eos_flag='OPAL01'
            if (use_debye_huckel_correction) eos_flag='OP1+DH'
         endif
         if (use_opal2006_eos) then
            eos_flag='OPAL06'
            if (use_debye_huckel_correction) eos_flag='OP6+DH'
         endif
      endif
! Determine low temperature opacities flag, LOK
      low_temp_opacity_flag='NONE'
      if (use_alex95_tables) low_temp_opacity_flag='ALEX'
      if (use_kurucz90_tables) low_temp_opacity_flag='KURZ'
! Determine high temperature opacities flag, HIK
      high_temp_opacity_flag='NONE'
      if (use_opal95_tables) high_temp_opacity_flag='OP95'
      if (use_opal92_tables) high_temp_opacity_flag='OP92'
      if (use_laol89_tables) high_temp_opacity_flag='LL89'

      if(atmosphere_flag .eq. ' ? ') then
         write(short_file_unit,7)
  7      format('*** YREC7 input file, flags, etc., have been ', &
                'defaulted.  ***')
      endif

! 09/25 JvS: Add secondary format option that prints stitched interior and envelope
! points. Output is either in the old format or new format, not both.
      if(lstch)then
          if(lphhd)then
            write(istor,1013) ! header key
          ! write model physics header. Should only happen upon first model output.
          ! write physics flags:
            write(istor,29) core_cz_top_index,envelope_cz_bottom_index,cmixl, &
           eos_flag,atmosphere_flag,low_temp_opacity_flag,high_temp_opacity_flag, &
           use_pure_z_table,initial_composition_code,use_extended_composition, &
           diffuse_helium_active,use_diffusion_z,lsemic,lovstc, &
           envelope_overshoot_active,lovstm,rotation_active, &
           instability_transport_active,ljdot0,disk_locking_active, &
           disk_temperature,disk_pressure,wind_saturation_omega,lstore,lstatm,lstenv, &
           lstmod,lstphys,lstrot
   29      format('#',2I8,F16.10,1X,A6,1X,3(A4,1X),L1,1X,A4,1X,11(L1,1X), &
           3(1PE18.10),1X,6(L1,1X))
           write(istor,1014) ! profile header
           lphhd = .false. ! Turn off the physics header fo the rest of the run.
          endif
        call stitch(composition,log_radius,log_pressure,log_density, &
             log_mass,log_temperature,log_luminosity,mass_coordinate,omega, &
             rotation_eta2,shell_moment_of_inertia,radius_ratio_r0, &
             specific_angular_momentum,shape_factor_fp,shape_factor_ft, &
             log_teff,log_total_mass,log_luminosity_lsun,num_shells, &
             convective_flag,model_number)
 1013     format('# JCORE  JENV  CMIXL  EOS  ATM  ALOK HIK  LPUREZ  COMPMIX', &
      '  LEXCOM  LDIFY  LDIFZ  LSEMIC  LOVSTC  LOVSTE  LOVSTM', &
      '  LROT  LINSTB  LJDOT0  LDISK  TDISK  PDISK  WMAX  LSTORE', &
      '  LSTATM  LSTENV  LSTMOD  LSTPHYS  LSTROT')
 1014     format( &
     'MODEL SHELL MASS RADIUS LUMINOSITY PRESSURE TEMPERATURE DENSITY OMEGA ', &
     'CONVECTIVE INTERIOR_PT ENV_PT ATM_PT H1 He4 METALS He3 C12 C13 N14 N15 O16 ', &
     'O17 O18 H2 Li6 Li7 Be9 OPACITY GRAV DELR DEL DELAD V_CONV GAM1 HII HEII HEIII ', &
     'BETA ETA PPI PPII PPIII CNO TRIPLE_ALPHA E_NUC E_NEU E_GRAV CP DLNRHODLNT A ', &
     'RP/RE FP FT J/M MOMENT DEL_KE V_ES V_GSF V_SS VTOT ')

       else
! write header records
      if(age_gyr .lt. 1d3) then
         write(istor,10) 'MOD2 ',model_number,num_shells,total_mass_msun, &
              log_teff,log_luminosity_lsun,log_total_mass,age_gyr, &
              timestep_yr,log_mass(1),log_mass(num_shells)
 10      format(A5,2I8,5F16.11,1PE18.10,0P2F16.12)
      else if (age_gyr .lt. 1D4) then
         write(istor,11) 'MOD2 ',model_number,num_shells,total_mass_msun, &
              log_teff,log_luminosity_lsun,log_total_mass,age_gyr, &
              timestep_yr,log_mass(1),log_mass(num_shells)
 11      format(A5,2I8,4F16.12,F16.10,1PE18.10,0P2F16.12)
      else if (age_gyr .lt. 1D5) then
         write(istor,12) 'MOD2 ',model_number,num_shells,total_mass_msun, &
              log_teff,log_luminosity_lsun,log_total_mass,age_gyr, &
              timestep_yr,log_mass(1),log_mass(num_shells)
 12      format(A5,2I8,4F16.12,F16.9,1PE18.10,0P2F16.12)
      else
         write(istor,13) 'MOD2 ',model_number,num_shells,total_mass_msun, &
              log_teff,log_luminosity_lsun,log_total_mass,age_gyr, &
              timestep_yr,log_mass(1),log_mass(num_shells)
 13      format(A5,2I8,4F16.12,F16.8,1PE18.10,0P2F16.12)
      endif

! write physics flags:
      write(istor,30) core_cz_top_index,envelope_cz_bottom_index,cmixl, &
           eos_flag,atmosphere_flag,low_temp_opacity_flag,high_temp_opacity_flag, &
           use_pure_z_table,initial_composition_code,use_extended_composition, &
           diffuse_helium_active,use_diffusion_z,lsemic,lovstc, &
           envelope_overshoot_active,lovstm,rotation_active, &
           instability_transport_active,ljdot0,disk_locking_active, &
           disk_temperature,disk_pressure,wind_saturation_omega,lstore,lstatm,lstenv, &
           lstmod,lstphys,lstrot
   30 format(2I8,F16.10,1X,A6,1X,3(A4,1X),L1,1X,A4,1X,11(L1,1X), &
           3(1PE18.10),1X,6(L1,1X))

! write luminosities and output flags
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
      write(istor,40) (luminosity_breakdown(j),j=1,7)
   40 format('TLUMX',5X,1P7E17.9)

! write ENVELOPE DATA
      do i = 1,3
         zone_sign_index = i
         if(trial_sign_flag.lt.0D0) zone_sign_index = -i
         write(istor,50)zone_sign_index,trial_log_temperature(i), &
              trial_log_luminosity(i),fit_point_pressure(i), &
              fit_point_temperature(i),fit_point_radius(i), &
              (envelope_fit_coeffs(3*i-3+j),j=1,3)
      enddo
   50 format('ENV',I2,5F16.12,1P3E20.12)

      call pindex(jxbeg,jxend,lshell,num_shells,id,idm)

      if(lstmod)then
! write column headings for all per shell information
         write(istor,55)
 55      format(/, &
     ' SHELL       MASS             RADIUS             LUMINOSITY            ', &
     'PRESSURE         TEMPERATURE         DENSITY               OMEGA      ', &
     '    C     H1          He4        METALS         He3             C12   ', &
     '          C13             N14             N15             O16         ', &
     '    O17             O18             H2              Li6             Li7', &
     '             Be9           OPAC       GRAV        DELR        DEL      ', &
     '   DELA       V_CONV     GAM1      HII     HEII     HEIII    BETA      ', &
     'ETA       PPI         PPII       PPIII        CNO         3HE         ', &
     'E_NUC        E_NEU       E_GRAV          A           RP/RE           FP', &
     '            FT           J/M          MOMENT        DEL_KE       V_ES ', &
     '      V_GSF      V_SS       VTOT   ',/)

! write out the requested information.
       local_g_const=dexp(ln10*cgl)
         do ii = 1,idm
            i = id(ii)
! write out the basic info
            write(istor,62,advance='no') i,log_mass(i),log_radius(i), &
                 log_luminosity(i),log_pressure(i), &
                 log_temperature(i),log_density(i),omega(i), &
                 convective_flag(i),(composition(j,i),j=1,15)
! write out additional physics if desired
            if(lstphys)then
             sg = dexp(ln10*(cgl - 2.0D0*log_radius(i)))*mass_coordinate(i)
               write(istor,63,advance='no') shell_diag%so(i),sg,shell_diag%del_grad(1,i),shell_diag%del_grad(2,i), &
                 shell_diag%del_grad(3,i),shell_diag%svel(i),adiabatic_index_gamma1(i),0.0,0.0,0.0,shell_diag%sbeta(i),shell_diag%seta(i), &
                 (shell_diag%seg(k,i),k=1,5),shell_diag%sesum(i),shell_diag%seg(6,i),shell_diag%seg(7,i)
!               WRITE(ISTOR,63,ADVANCE='no') SO(I),SG,SDEL(1,I),SDEL(2,I),
!     *           SDEL(3,I),SVEL(I),GAM1(I),SFXION(1,I),SFXION(2,I),SFXION(3,I),
!     *           SBETA(I),SETA(I),(SEG(K,I),K=1,5),SESUM(I),SEG(6,I),SEG(7,I)
            else
               write(istor,63,advance='no') 0.0,0.0,0.0,0.0,0.0,0.0,0.0, &
                 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
            endif
! write out additional rotation info if desired
            if(lstrot.and.rotation_active)then
             fm = dexp(ln10*log_mass(i))
             duma = cc13*omega(i)**2/(local_g_const*fm)*5.D0/(2.D0+rotation_eta2(i))
             oblateness_a = duma * radius_ratio_r0(i)**3
             pole_to_equator_ratio = (1.0D0 - oblateness_a)/(1.0D0 + 0.5D0*oblateness_a)
               vtot = es_circulation_velocity(i)+gsf_circulation_velocity(i)+secular_shear_velocity(i)
               write(istor,64) oblateness_a,pole_to_equator_ratio,shape_factor_fp(i), &
                  shape_factor_ft(i),specific_angular_momentum(i),shell_moment_of_inertia(i), &
                  rotational_energy_term(i),es_circulation_velocity(i), &
                  gsf_circulation_velocity(i),secular_shear_velocity(i),vtot
            else
               write(istor,64) 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
            endif
         enddo
 62   format(I6,0P2F18.14,1PE24.16,0P3F18.14,1PE24.16,1x,L1, &
           3(0PF12.9),12(0PE16.8),2X)

 63   format(1PE10.4,1PE11.3,E12.4,E12.4,E12.4,1PE12.4,0PF9.5,F9.5,F9.5,F9.5, &
!     &     F9.5,F9.5,F9.5,F9.5,F9.5,F9.5,F9.5,E13.5,E13.5,E13.5)
           F9.5,F9.5,E12.4,E12.4,E12.4,E12.4,E12.4,E13.5,E13.5,E13.5)

 64   format(E14.6,E14.6,E14.6,E14.6,E14.6,E14.6,E14.6,E11.3,E11.3,E11.3,E11.3, &
           E11.3)

      endif
! now call wrtmod, with the goal of outputting the envelope and atmosphere, or
! if required by LPULSE.
      if(lstatm.or.lstenv)then
       if(lmilne_local) call wrtmil(composition,log_density,log_luminosity, &
            log_pressure,log_radius,mass_coordinate,num_shells,model_number)
         call wrtmod(num_shells,envelope_cz_bottom_index,composition, &
              mass_coordinate,log_density,log_luminosity,log_pressure, &
              log_radius,log_temperature,model_number,log_luminosity_lsun, &
              log_teff,shape_factor_fp,shape_factor_ft,log_mass,age_gyr)
      endif

      write(istor,65)
 65   format(/,/)

       endif !closes logic block for new vs old store file format

! G Somers END
!CCCCCCCCCCCCCCCCCCCCCCCCC

      return
end subroutine putstore
