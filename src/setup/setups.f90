!----------------------------------------------------------------------
! setups
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original setups.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! MHP 8/25 Passed file names directly instead of through common, as
! character strings in common blocks can cause problems.
!
! One-time model-independent setup performed at program start:
! defines the physical/mathematical constants used throughout the
! code (common/const1//const2//const3//const/debhu), loads the
! opacity tables (setupopac/mhdtbl), the degenerate-electron
! (Fermi-Dirac) equation-of-state table used by eqrelv.f90, the
! Kurucz/Castelli surface-pressure table used for the T-tau surface
! boundary condition, and (if enabled) the SCVH envelope
! equation-of-state tables.
subroutine setups(laol_work_array, alex06_table_path, allard_table_path, &
     atm_table_path, fermi_table_path, kurucz_table_path, &
     kurucz_table2_path, laol_table_path, laol_table2_path, &
     opal95_table_path, opal92_table_path, zams_a_table_path, &
     zams_b_table_path, zams_c_table_path, centre1_table_path, &
     centre2_table_path, centre3_table_path, centre4_table_path, &
     centre5_table_path, opal92_table2_path, pure_z_table_path, &
     scv_h_table_path, scv_he_table_path, scv_z_table_path, &
     alex95_table_paths)
      use luout_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000
! JNT 06/14 ADD NTC FOR KURUCZ/CASTELLI 2004 ATM
      integer, parameter :: nt = 57, ng = 11
      integer, parameter :: ntc = 76, ngc = 11
      integer, parameter :: nts = 63, nps = 76
! MHP 8/97 ADDED NTA AND NGA FOR ALLARD ATMOSPHERE
      integer, parameter :: nta = 54, nga = 5

      double precision, intent(inout) :: laol_work_array(12)
! MHP 8/25 Reduced declared variables to ones actually used here or
! passed to other routine
      character(len=256), intent(in) :: alex06_table_path, &
           allard_table_path, atm_table_path, fermi_table_path, &
           kurucz_table_path, kurucz_table2_path, opal95_table_path, &
           laol_table_path, laol_table2_path, opal92_table_path, &
           zams_a_table_path, zams_b_table_path, zams_c_table_path, &
           centre1_table_path, centre2_table_path, centre3_table_path, &
           centre4_table_path, centre5_table_path, &
           opal92_table2_path, pure_z_table_path, scv_h_table_path, &
           scv_he_table_path, scv_z_table_path
      character(len=256), intent(in) :: alex95_table_paths(7)

      common/lunum/ first_unit, run_unit, standard_unit, fermi_unit, &
           opal_model_unit, opal_envelope_unit, opal_atm_unit, &
           dynamics_unit, laol_table_unit, neutrino_unit, &
           composition_unit, kurucz_table_unit
      integer :: first_unit, run_unit, standard_unit, fermi_unit, &
           opal_model_unit, opal_envelope_unit, opal_atm_unit, &
           dynamics_unit, laol_table_unit, neutrino_unit, &
           composition_unit, kurucz_table_unit
! MHP 8/25 Removed file names from common block
!      COMMON/LUFNM/ FLAST, FFIRST, FRUN, FSTAND, FFERMI,
!     1    FDEBUG, FTRACK, FSHORT, FMILNE,  FMODPT,
!     2    FSTOR, FPMOD, FPENV, FPATM, FDYN,
!     3    FLLDAT, FSNU, FSCOMP, FKUR,
!     4    FMHD1, FMHD2, FMHD3, FMHD4, FMHD5, FMHD6, FMHD7, FMHD8

! DBGLAOL
! MHP 8/25 Removed character file names from common block
! common/nwlaol/: not used in this file; declared only to preserve
! layout. Naming matches getopac.f90.
      double precision :: olaol(12,104,52), oxa(12), ot(52), orho(104), &
           tollaol
      integer :: iolaol, numofxyz, numrho, numt, iopurez
      logical :: llaol, use_pure_z_table
      common/nwlaol/ olaol, oxa, ot, orho, tollaol, iolaol, numofxyz, &
           numrho, numt, llaol, use_pure_z_table, iopurez
! DBG 7/92 COMMON BLOCK ADDED TO COMPUTE DEBYE-HUCKEL CORRECTION.
! common/debhu/: all members assigned here. Naming matches
! eqstat2.f90.
      double precision :: cdh, etadh0, etadh1, zdh(18), xxdy, yydh, zzdh, &
           dhnue(18)
      logical :: ldh
      common/debhu/ cdh, etadh0, etadh1, zdh, xxdy, yydh, zzdh, dhnue, ldh

! common/comp/: only envelope_hydrogen_fraction is used here. Naming
! matches getopac.f90.
      double precision :: envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, amuenv, fxenv(12), xnew, znew, stotal, senv
      common/comp/ envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, amuenv, fxenv, xnew, znew, stotal, senv

! common/const/: all members assigned here. Naming matches
! chkscal.f90/wrthead.f90/vcirc.f90.
      double precision :: solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude
      common/const/ solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude




! common/ccr/: this batch's own block (see eqrelv.f90 for the full
! member-by-member discussion); all read in here from the Fermi table
! file. Naming matches eqrelv.f90.
      double precision :: fermi_table_x_grid(43), fermi_table_eta(43), &
           fermi_table_data(5,43,20)
      integer :: fermi_table_x_lookup(261)
      common/ccr/ fermi_table_x_grid, fermi_table_eta, fermi_table_data, &
           fermi_table_x_lookup

! common/mhd/: only use_mhd_eos is used here. Naming matches
! mhdtbl.f90/qatm.f90.
      logical :: use_mhd_eos
      integer :: unit_zams_a, unit_zams_b, unit_zams_c, unit_centre1, &
           unit_centre2, unit_centre3, unit_centre4, unit_centre5
      common/mhd/use_mhd_eos, unit_zams_a, unit_zams_b, unit_zams_c, &
           unit_centre1, unit_centre2, unit_centre3, unit_centre4, &
           unit_centre5

! common/atmos/: atm_hras/atm_choice are used here; atm_choice_initial/
! use_ttau_relation are unused placeholders. Naming matches
! surfbc.f90/qatm.f90.
      double precision :: atm_hras
      integer :: atm_choice, atm_choice_initial
      logical :: use_ttau_relation
      common/atmos/atm_hras, atm_choice, atm_choice_initial, use_ttau_relation
      double precision :: kurucz_log10_pressure_table(nt,ng), &
           kurucz_teff_table(nt), kurucz_logg_table(ng), kurucz_table_z
      integer :: atm_table_file_unit
! MHP 8/25 Removed file names from common block
! JMH 8/18/91
!      COMMON/ATMOS2/ATMPL(NT,NG),ATMTL(NT),
!     *              ATMGL(NG),ATMZ,IOATM,FATM
! common/atmos2/: the Kurucz surface-pressure table, all used here.
! Naming matches surfp.f90.
      common/atmos2/kurucz_log10_pressure_table, kurucz_teff_table, &
           kurucz_logg_table, kurucz_table_z, atm_table_file_unit
! JNT 06/14
! common/atmos2c/: the Kurucz/Castelli surface-pressure table, all
! used here. Naming matches kcsurfp.f90.
      double precision :: kurucz_castelli_log10_pressure_table(ntc,ngc), &
           kurucz_castelli_teff_table(ntc), kurucz_castelli_logg_table(ngc)
      common/atmos2c/kurucz_castelli_log10_pressure_table, &
           kurucz_castelli_teff_table, kurucz_castelli_logg_table
! MHP 8/92 COMMON BLOCK ADDED FOR LOWER EDGE OF TABLE IN LOG G.
! common/fac/: all members assigned here. Naming matches surfp.f90.
      integer :: kurucz_gmin_index(nt), kurucz_gmax_index(nt), &
           teff_interp_start_index, gravity_interp_indices(4), &
           castelli_gmin_index(ntc), castelli_gmax_index(ntc)
      common/fac/kurucz_gmin_index, kurucz_gmax_index, &
           teff_interp_start_index, gravity_interp_indices, &
           castelli_gmin_index, castelli_gmax_index
! MHP  5/97 ADDED COMMON BLOCK FOR SCV EOS TABLES
! common/scveos/: tlogx/nptsx/tablex/tabley/tablez used here;
! smix/tablenv/idt/idp are unused placeholders. Naming matches
! eqstat2.f90.
      double precision :: tlogx(nts), tablex(nts,nps,12), &
           tabley(nts,nps,12), smix(nts,nps), tablez(nts,nps,13), &
           tablenv(nts,nps,12)
      integer :: nptsx(nts), idt, idp
      logical :: use_scv_eos
      common/scveos/ tlogx, tablex, tabley, smix, tablez, tablenv, nptsx, &
           use_scv_eos, idt, idp
! MHP 8/25 Removed file names from common block
! common/scv2/: not used in this file; declared only to preserve
! layout. Not referenced in any already-converted file.
      integer :: scv_h_unit, scv_he_unit, scv_z_unit
      common/scv2/ scv_h_unit, scv_he_unit, scv_z_unit

      save

! --- locals ---
      double precision :: speed_of_light, electron_mass, boltzmann_constant, &
           planck_constant, hydrogen_atom_mass, electron_charge_esu
      double precision :: hra
      external hra
      integer :: grid_idx, iden_idx, bin_width, bin_start, bin_end, card_idx
      integer :: teff_idx, logg_idx
      logical :: found_valid_pressure
      double precision :: scvhe_dummy_val, scvz_dummy_val
      integer :: scvhe_dummy_npts, scvz_dummy_npts
      integer :: t_idx, p_idx, col_idx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SETUP CONSTANTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      clndp = dlog(10.0d0)
      ln10 = clndp
      clni = 1.0d0/ln10
! Luminosity of Sun
      log10_solar_luminosity = dlog10(solar_luminosity_cgs)
      ln_solar_luminosity = ln10/solar_luminosity_cgs
! Mass of Sun
      solar_mass_cgs = 1.9891d33
      log10_solar_mass = dlog10(solar_mass_cgs)
! Radius of Sun
      log10_solar_radius = dlog10(solar_radius_cgs)
! Bolometric magnitude of the sun
      solar_bolometric_magnitude = 4.79d0
! No. of seconds per year and other mathematical constants
      seconds_per_year = 3.1558d7
      cc13 = 1.0d0/3.0d0
      cc23 = cc13 + cc13
      cpi = 3.1415926535898d0
      c4pi = 4.0d0*cpi
      c4pil = dlog10(c4pi)
      c4pi3l = dlog10(cc13*c4pi)
! Speed of Light
      speed_of_light = 2.99792458d10
! Stefan-Boltzmann constant
      csig = 5.67051d-5
      csigl = dlog10(csig)
! radiation constant/3
      radiation_constant_over_3 = csig*4.0d0/(3.0d0*speed_of_light)
      ca3l = dlog10(radiation_constant_over_3)
! molar gas constant
      gas_constant = 8.314510d7
! log of Gravitational constant G=6.67259D-8
      cgl = -7.17571d0
! mass of the electron
      electron_mass = 9.1093897d-28
! Boltzmann constant
      boltzmann_constant = 1.380658d-16
! Planck's constant
      planck_constant = 6.6260755d-27
      cmkh = 1.50d0*dlog10(2.0d0*cpi*electron_mass)+2.5d0* &
           dlog10(boltzmann_constant)- &
           3.0d0*dlog10(planck_constant)
      cdelrl = -c4pi3l - csigl - dlog10(16.0d0)
      cmixl2 = cc13
      cmixl3 = 16.0d0*dsqrt(2.0d0)*csig
! DBG CALCULATE DEBYE-HUCKEL CONSTAND CDH
! mass hydrogen atom (gm)
      hydrogen_atom_mass = 1.0d0/6.0222137d23
! charge electron (ESU)
      electron_charge_esu=4.802d-10
      cdh = -sqrt(cpi/(boltzmann_constant*hydrogen_atom_mass* &
           hydrogen_atom_mass*hydrogen_atom_mass))*electron_charge_esu* &
           electron_charge_esu*electron_charge_esu/3.0d0
! Electric charge of C,N,O,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,Ca,Ti,Cr,Mn,Fe,Ni
      dhnue(1) = 6.0d0
      dhnue(2) = 7.0d0
      dhnue(3) = 8.0d0
      dhnue(4) = 10.0d0
      dhnue(5) = 11.0d0
      dhnue(6) = 12.0d0
      dhnue(7) = 13.0d0
      dhnue(8) = 14.0d0
      dhnue(9) = 15.0d0
      dhnue(10) = 16.0d0
      dhnue(11) = 17.0d0
      dhnue(12) = 18.0d0
      dhnue(13) = 20.0d0
      dhnue(14) = 22.0d0
      dhnue(15) = 24.0d0
      dhnue(16) = 25.0d0
      dhnue(17) = 26.0d0
      dhnue(18) = 28.0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     EVALUATE TAU = 2/3 TEMPERATURE FOR HRA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      atm_hras = hra(cc23)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     SET UP OPACITY TABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call setupopac(envelope_hydrogen_fraction, laol_work_array, &
           alex06_table_path,kurucz_table_path,kurucz_table2_path, &
           laol_table_path,laol_table2_path, &
           opal95_table_path,opal92_table_path,opal92_table2_path, &
           pure_z_table_path,alex95_table_paths)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     READ IN MHD EOS TABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MHP 8/25 Passed file names directly instead of through common, as
! character strings in common blocks can cause problems
      if(use_mhd_eos) then
           call mhdtbl(zams_a_table_path,zams_b_table_path, &
                zams_c_table_path,centre1_table_path,centre2_table_path, &
                centre3_table_path,centre4_table_path,centre5_table_path)
      endif
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INPUT F-TABLES FOR DEGENERATE EQUATION OF STATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  OPEN DATA FILES
      open(unit=fermi_unit,file=fermi_table_path,form='FORMATTED',status='OLD')
      read(fermi_unit,150)  (fermi_table_x_grid(grid_idx), grid_idx = 1,43)
      read(fermi_unit,150)  (fermi_table_eta(grid_idx), grid_idx = 1,43)
      read(fermi_unit,160) (grid_idx,iden_idx,(fermi_table_data(col_idx,grid_idx,iden_idx), &
           col_idx=1,5),card_idx=1,860)
  150 format(8x,9f8.4)
  160 format(2i5,5f12.7)
      rewind fermi_unit
      bin_start = 1
      do 180 grid_idx=1,41
! MHP 10/02 SHOULD BE INT(DVAL,ETC.)
       if (fermi_table_x_grid(grid_idx+1).le.fermi_table_x_grid(grid_idx)) then
!       IF (INT(DVAL(I+1)).LE.INT(DVAL(I))) THEN
          write(short_file_unit,1000) grid_idx
 1000       format(1x,39('>'),40('<')/1x,'ERROR IN SUBROUTINE SETUPS'/ &
           1x,'GLITCH IN FERMI TABLE ELEMENT',i4/1x,'RUN STOPPED')
          stop
       endif
       bin_width = int((fermi_table_x_grid(grid_idx+1) - fermi_table_x_grid(grid_idx))*20.0d0 + 0.10d0)
       bin_end = bin_start + bin_width - 1
       if (grid_idx.eq.41)  bin_end = 261
       do 170 iden_idx=bin_start,bin_end
          fermi_table_x_lookup(iden_idx) = grid_idx
  170    continue
       bin_start = bin_end + 1
  180 continue

!  CLOSE EQUATION OF STATE FILE.
      close(fermi_unit)
! JMH 8/18/91
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INPUT PRESSURE TABLE FOR SURFACE BOUNDARY CONDITIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ((atm_choice .eq. 3) .or. (atm_choice .eq. 4)) then
! OPEN SURFACE PRESSURE TABLE
        open(atm_table_file_unit,file=atm_table_path, status='OLD')
! GET ABUNDANCE:
        read(atm_table_file_unit,200) kurucz_table_z
  200   format(1x,10x,1pe16.8,/)
! GET VALUES OF LOG Teff:
        read(atm_table_file_unit,201) (kurucz_teff_table(teff_idx),teff_idx=1,nt)
  201   format(1p4e16.8,13(/1p4e16.8),/1pe16.8,/)
! GET VALUES OF LOG G:
        read(atm_table_file_unit,202) (kurucz_logg_table(logg_idx),logg_idx=1,ng)
  202   format(1p4e16.8,/1p4e16.8,/1p3e16.8,/)
! GET GRID OF LOG PRESSURE VALUES:
        do 205 teff_idx=1,nt
          read(atm_table_file_unit,203) (kurucz_log10_pressure_table(teff_idx,logg_idx),logg_idx=1,ng)
  203     format(1p4e16.8,/1p4e16.8,/1p3e16.8)
  205   continue
        rewind atm_table_file_unit
        close(atm_table_file_unit)
!       G Somers 5/15; added a catch that allows the code to work
!                      if the highest gravity P value for a specific
!                      temperature is -999. LCATCH is set to true
!                      once the first valid P is read, and it can
!                      set the minimum P thereafter.
         do teff_idx = 1,nt
            found_valid_pressure = .false.
            do logg_idx = ng,1,-1
               if(kurucz_log10_pressure_table(teff_idx,logg_idx).le.0.0d0)then
!                 check if the first non-999 value has been reached.
!                 if so, set IMIN.
                  if(found_valid_pressure)then
                     kurucz_gmin_index(teff_idx) = logg_idx + 1
                     goto 5
                  endif
               else
!                 if we have reached a non-negative pressure value,
!                 turn off the catch so IMIN can be set. also record
!                 the highest gravity with a pressure for later int.
                  if(.not.found_valid_pressure) kurucz_gmax_index(teff_idx) = logg_idx
                  found_valid_pressure = .true.
               endif
            end do
            kurucz_gmin_index(teff_idx) = 1
    5       continue
!           if all of the P values at a given T are -999, set IMIN
!           to the number of gravity terms. in responce, the code
!           should break when trying to find surface P.
            if(.not.found_valid_pressure) kurucz_gmin_index(teff_idx) = ng
         end do
!        G Somers 5/15 END
! MHP 6/97 ADDED OPTION FOR ALLARD MODEL ATMOSPHERES; USED INSTEAD OF
! KURUCZ FOR TEFF < 10,000 K.
         if(atm_choice .eq. 4)then
!            ATMZA = 0.02D0
          call alfilein(allard_table_path)      ! Get Allard Atmospheres files and
         endif                  ! initialize tables. 9/23/08 LLP


! JNT 6/14 ADD OPTION FOR NEW KURUCZ/CASTELLI ATMOSHPERES
      else if (atm_choice .eq. 5) then
! OPEN SURFACE PRESSURE TABLE
        open(atm_table_file_unit,file=atm_table_path, status='OLD')
! GET ABUNDANCE:
        read(atm_table_file_unit,200) kurucz_table_z
! GET VALUES OF LOG Teff:
        read(atm_table_file_unit,206) (kurucz_castelli_teff_table(teff_idx),teff_idx=1,ntc)
  206   format(1p4e16.8,18(/1p4e16.8),/)
! GET VALUES OF LOG G:
        read(atm_table_file_unit,208) (kurucz_castelli_logg_table(logg_idx),logg_idx=1,ngc)
  208   format(1p4e16.8,/1p4e16.8,/1p3e16.8,/)
! GET GRID OF LOG PRESSURE VALUES:
        do 207 teff_idx=1,ntc
          read(atm_table_file_unit,203) (kurucz_castelli_log10_pressure_table(teff_idx,logg_idx),logg_idx=1,ngc)
  207   continue
        rewind atm_table_file_unit
        close(atm_table_file_unit)
!       G Somers 5/15; added a catch that allows the code to work
!                      if the highest gravity P value for a specific
!                      temperature is -999. LCATCH is set to true
!                      once the first valid P is read, and it can
!                      set the minimum P thereafter.
        do teff_idx = 1,ntc
           found_valid_pressure = .false.
           do logg_idx = ngc,1,-1
              if(kurucz_castelli_log10_pressure_table(teff_idx,logg_idx).le.0.0d0)then
!                check if the first non-999 value has been reached.
!                if so, set IMIN2.
                 if(found_valid_pressure)then
                    castelli_gmin_index(teff_idx) = logg_idx + 1
                    goto 6
                 endif
              else
!                if we have reached a non-negative pressure value,
!                turn off the catch so IMIN2 can be set. also record
!                the highest gravity with a pressure for later int.
                 if(.not.found_valid_pressure) castelli_gmax_index(teff_idx) = logg_idx
                 found_valid_pressure = .true.
              endif
           end do
           castelli_gmin_index(teff_idx) = 1
    6      continue
!          if all of the P values at a given T are -999, set IMIN
!          to the number of gravity terms. in responce, the code
!          should break when trying to find surface P.
           if(.not.found_valid_pressure) castelli_gmin_index(teff_idx) = ng
        end do
! END JNT 6/14
      endif

! MHP 5/97 ADDED OPTION FOR NEW SCV EQUATION OF STATE TABLES.
      if(use_scv_eos)then
         open(unit=scv_h_unit,file=scv_h_table_path,status='OLD')
         open(unit=scv_he_unit,file=scv_he_table_path,status='OLD')
         open(unit=scv_z_unit,file=scv_z_table_path,status='OLD')
!  READ IN EQUATION OF STATE TABLES FOR HYDROGEN AND HELIUM
         do t_idx = 1, nts
            read(scv_h_unit,1) tlogx(t_idx),nptsx(t_idx)
            read(scv_he_unit,1) scvhe_dummy_val,scvhe_dummy_npts
            read(scv_z_unit,1) scvz_dummy_val,scvz_dummy_npts
    1       format(f5.2,i4)
! TABLE GRID POINTS IN T, P(T) ARE THE SAME - NPTSY AND TLOGX
! READ IN TO RETAIN PARALLEL COMMON BLOCK STRUCTURE.
            do p_idx = 1, nptsx(t_idx)
               read(scv_h_unit,2) (tablex(t_idx,p_idx,col_idx),col_idx=1,11)
               read(scv_he_unit,2) (tabley(t_idx,p_idx,col_idx),col_idx=1,11)
    2          format(f6.2,1p2e13.5,0p,8f9.4)
            end do
!  READ IN METAL EQUATION OF STATE TABLE; COMPUTED USING THE PRATHER
! EQUATION OF STATE IN THE OLD YALE CODE.
            do p_idx = nptsx(t_idx),1,-1
               read(scv_z_unit,3)(tablez(t_idx,p_idx,col_idx),col_idx=1,13)
 3             format(f6.2,12f9.4)
            end do
         end do
      endif

      return
end subroutine setups
