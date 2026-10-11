!----------------------------------------------------------------------
! read_controls
!----------------------------------------------------------------------
! Added 2026 (phase five): the controls-read step of the former
! program main, as its own star-layer routine. parmin reads the
! control and physics namelists into the controls_buffer BUFFER and
! fills the job-level file paths (star%job). Errors return via ierr.
!
! 2026 controls->star% campaign, phase B: star%ctrl is the
! authoritative post-read home. The four-step sequence below
! guarantees parmin always reads over pristine defaults (a namelist
! read only overwrites what the input file provides) -- it replaces
! the retired controls_reset_lib snapshot machinery structurally:
!   1. star%ctrl = controls_state()   pristine defaults (component
!                                     initializers, generated from
!                                     the buffer's declarations)
!   2. seed_controls_buffer           buffer <- star%ctrl
!   3. parmin                         namelist reads over the buffer
!   4. store_controls_to_star         star%ctrl <- buffer
! Until phase C migrates them, consumers still read the buffer's
! bare names; star%ctrl is written here but not yet read broadly.
module controls_reader
      use star_info_lib, only: star, control_nml_override, physics_nml_override
      use controls_sync, only: store_controls_to_star
      use physical_constants
      use io_units
      use legacy_controls
      use run_files, only: resolve_output_mode_and_paths, derive_options_and_open_files
      use derived_controls, only: echo_settings
      use run_cards, only: interpret_kind_cards
      implicit none
      private
      public :: read_controls

contains

subroutine read_controls(ierr)

      use star_info_lib, only: star, controls_state
      use controls_sync, only: seed_controls_buffer, &
           store_controls_to_star
      implicit none

      integer, intent(out) :: ierr

      ierr = 0
      star%ctrl = controls_state()
      call seed_controls_buffer
! read in user parameters (read_input fills star%job's table paths
! itself -- see its load/store_table_paths pair below)
      call read_input(ierr)
      if (ierr /= 0) return
      call store_controls_to_star

      return
end subroutine read_controls

!----------------------------------------------------------------------
! parmin
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original parmin.f; only source form was updated -- see the naming
! note below for why variable names are treated differently here.
!
! Reads the CONTROL and PHYSICS namelists (yrec8.nml1/yrec8.nml2, or
! the two files named on the command line) that configure an entire
! YREC run: envelope/atmosphere physics switches, opacity/EOS table
! selection, diffusion and mixing options, rotation and wind-law
! parameters, nuclear cross sections, and the sequence of "kind" run
! cards (evolve / rescale / rescale+evolve) executed by the driver.
! Also opens every logical unit used for the run's input/output files
! and writes a human-readable summary of the active run parameters to
! the .short log.
!
! NAMING NOTE: unlike other files in this modernization project, the
! members of every COMMON block below, and every variable listed in
! NAMELIST /control/ or /physics/, keep their EXACT original spelling
! (only lowercased). Fortran's namelist reader matches input-file
! entries (e.g. "LCALS = .TRUE." in a *.nml1 file) to variables by
! name, against every *.nml1 and *.nml2 file in the repository -- so
! renaming any of these would silently break every existing input
! file that sets it. This also means every dummy argument of this
! particular subroutine happens to be a NAMELIST /control/ member
! too (the namelist read sets these arguments' values directly, which
! are then forwarded to other routines), so none of them could be
! renamed either. Only true local variables (never referenced by
! COMMON, NAMELIST, or the argument list) were renamed for
! readability; see also the two SUBROUTINE EXPAND_VALUE locals below,
! a fully-local scope with no such restriction.
subroutine read_input(ierr)
      integer, intent(out) :: ierr

      call reset_legacy_controls
      ierr = 0
! The table-path variables above keep their legacy NAMELIST /control/
! spellings (fatm, falex06, ...) as purely local read buffers; their
! post-read home is star%job (canonical *_table_path names). The
! load/store pair below replaces the former 24-item argument list
! through which read_controls used to pass the star%job members.
      call load_table_paths
!
! The DATA cards that once set every run parameter's default here are
! gone: defaults now live on the buffer declarations in
! io/controls_buffer.f90 (a use-associated variable cannot be a DATA
! target), and the working mixing length's default (1.4d0) sits on
! its star%mixing_length_alpha component. The historical per-author
! DATA-card commentary is in git history (pre-2026 parmin.f).
!
! MHP 10/24 MIXTURE CONTROL: ISETISO CONTROLS CNO ISOTOPE RATIOS AND
! LIGHT ELEMENT ABUNDANCES D,HE3,LI6,LI7,BE9,B10,B11 (1=USED)
! ISETMIX CONTROLS C+N+O MASS FRACTIONS (1=USED)
! AMIX AND AISO ARE STRINGS IDENTIFYING EITHER A PRESET MIXTURE OR A CUSTOM ONE ('CUS')
! SUPPORTED AMIX AT PRESENT ARE 'GS98','AAG21',M22M','M22P'. SUPPORTED AISO IS 'L21'.
! FOR A CUSTOM MIXTURE YOU CAN ENTER INDIVIDUAL VALUES.
! TO BE ADDED - AMIX FROM ATOMIC OPACITY TABLES (INEWMIX=2) AND OTHER MIXTURES/ISOTOPES
!     L21 DEFAULT ISOTOPE DATA LODDERS ET AL. 2021 SSRV 217,44
!     GS98 DEFAULT CNO FRACTIONS OF METALS GREVESSE&sAUVAL 1998 SSRV 85,161
! MHP 10/24 DATA FOR CNO FRACTIONS AND Z/X OF DIFFERNT SOLAR MIXTURES.
!     AMIXT IS THE LIST OF IDS,EACH OF WHICH HAS A ZX AND CNO FRACS
!     ENTRY 1 =GS98(IN PARMIN),2=ASPLUND ET AL. 2021 A&A 653,141
!     3,4=MAGG ET AL. 2021 (MET,PHOT) A&A 661,140
!
! DBG 1/96 THE ARRAY V, READ IN VIA read_laol_tables, CONTAINED THE MASS FRACTIONS
! OF THE ENVELOPE ELEMENTS. IT WAS USED IN STARIN TO DEFINE FXENV,
! WHICH ARE THE NUMBER FRACTION OF THE ENVELOPE ELEMENTS. FXENV WAS
! THEN UPDATED IN compute_blended_eos AND EQSAHA. VNEW (now a controls_buffer member)
! IS IDENTICAL TO V EXCEPT THAT THE NUMBERS ARE DEFINED EXPLICITLY
! FOR A G&N93 SOLAR MIXTURE. YOU CAN CHANGE THEM VIA THE PHYSICS
! NAMELIST. V IS SET EQUAL TO VNEW IN STARIN EXCEPT WHEN LLAOL=T
! (TO MAINTAIN BACKWARD COMPATIBILITY).
!            Na          Al          Mg          Fe
!            Si          C           H           O
!            N           Ar          Ne          He
!
! G Somers 06/14 ALLOW NEW LI DESTRUCTION CROSS SECTIONS
!           NEW VALUES SHOULD BE IN UNITS OF keV b.
!           DEFAULT LI6 = 5,500 keV b FROM FOWLER ET AL. 1967
!           DEFAULT LI7 = 52 keV b FROM ROLFS & KAVANAGH 1986
!           DEFAULT BE9(P,G)B10 (BE91) = 1.1 keV b INFERRED FROM FOWLER 1988
!               THEY DON'T SHOW THIS REACTION, SO THIS IS APPROXIMATE.
!           DEFAULT BE9(P,D)2HE4 (BE92) = 15,000 keV b FROM FOWLER ET AL. 1967
!           DEFAULT BE9(P,A)LI6 (BE93) = 15,000 keV b FROM FOWLER ET AL. 1967
! MHP 8/14 DEFAULT CROSS-SECTIONS ARE TAKEN FROM THE SOLAR FUSION II PAPER
! REFERENCE ADELBERGER ET AL. 2011. UNITS ARE KeV b
! NOTE: PEP IS THE PROPORTIONALITY CONSTANT RELATIVE TO PP
! NOTE: BE7+E- IS THE PROPORTIONALITY CONSTANT IN THE LINEAR TERM
! THE CODE USES T9, NOT T6, SO ANY EXPRESSION IN TERMS OF T/10^6 K
! NEEDS TO BE DIVIDED BY 1000^0.5 (FOR BOTH PEP AND BE7+E-)
! (first/second derivatives of the cross sections likewise from
! Adelberger et al. 2011, units b)
!
! THIS SUBROUTINE READS ALL USER DEFINED QUANTITIES FROM THE
! FILES yrec8.nml1 and yrec8.nml2
! VALUES FOR LOGICAL UNITS USED IN READS AND WRITES SET IN DATA
! STATEMENT; THEY CAN BE CHANGED IN THE NAMELIST IF NEEDED.
! LOGICAL UNIT 5 = READ FROM SCREEN
! LOGICAL UNIT 6 = WRITE TO SCREEN: FOR BATCH USE STATUS FILE INSTEAD
! SPECIFY ALL LOGICAL UNIT NUMBERS HERE:
! OUTPUT: TERMINAL (STDOUT)
! 2026 re-entrancy fix: these legacy-spelled NAMELIST locals used to be
! initialised by DATA statements and initialised declarations, which
! implicitly SAVE them -- so a second in-process yrec_run inherited
! whatever the previous inlist set for any control the new inlist
! omitted (the canonical star%ctrl buffer is reset by read_controls,
! these were not).  They are now assigned their defaults on every
! entry.  Values are unchanged.
      alphac = 0.0d0
      alphae = 0.0d0
      dtdif = 1.0d-2
      wnew = 0.0d0
      creim = -4.0d-13
      cstmixing = 1.0d0
      cstdiffmix = 1.0d0
      tcut = (/6.5d0,6.5d0,6.82d0,7.7d0,7.5d0/)
      vnew = (/ 0.001999d0, 0.003238d0, 0.037573d0, 0.071794d0, 0.040520d0, 0.173285d0, &
           0.000000d0, 0.482273d0, 0.053152d0, 0.005379d0, 0.098668d0, 0.000000d0 /)
      itdif1 = 1
      itdif2 = 1
      niter1 = 2
      niter2 = 20
      niter3 = 2
      niter4 = 0
      ldifli = .false.
      lsnu = .false.
      lovstc = .false.
      lwnew = .false.
      ljdot0 = .true.
      lthoul = .false.
      lsemic = .false.
      lnews = .false.
      lmonte = .false.
      lreimer = .false.
      imbeg = 1
      imend = 1
      weakscreening = 0.03d0
      sstandard = (/ 0.9828, 1.0485, 0.9815, 0.9241, 1.3818, 1.0542, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &
           1.0, 1.0, 1.0108, 0.7819, 0.2875 /)
      ldifad = .false.
      lnoj = .false.
      tdisk = 0.0d0
      pdisk = 7.2722d-6
      ldisk = .false.
      lsenv0a = .false.
      senv0a = 1.26d-4
      lenvg = .false.
      atmstp = 0.5
      envstp = 0.5
      numrun = 1
      kindrn = 1
      lfirst = .true.
      nmodls = 0
      endage = 0.0
      setdt = 0.0
      rsclm = 0.0
      rsclx = 0.0
      rsclz = 0.0
      rsclcm = 0.0
      end_dcen = 0.0
      end_xcen = 0.0
      end_ycen = 0.0
      element_id = (/ 'HE3', 'C12', 'C13', 'N14', 'N15', 'O16', 'O17', 'O18', 'H2 ', 'LI6', 'LI7', 'BE9' /)
      optol = 1.0d-8
      atmerr = 3.0d-4
      atmd0 = 1.0d-10
      atmbeg = 1.0d-1
      atmmin = 1.0d-1
      atmmax = 5.0d-1
      enverr = 3.0d-4
      envbeg = 1.0d-1
      envmin = 1.0d-1
      envmax = 5.0d-1
      stolr0 = 1.0d-3
      imax = 11
      nuse = 7
      djok = 1.0d-4
      hpttol = (/ 1.0d-8, 8.0d-2, 5.0d-2, 5.0d-2, 1.0d0, 1.0d0, 0.0d0, 5.0d-2, 2.0d-2, 5.0d-2, &
           5.0d-2, 1.0d-1 /)
      lnewcp = .false.
      anewcp = '   '
      xnewcp = 1.3d1
      isetiso = 0
      isetmix = 0
      amix = 'GS98'
      aiso = 'L21'
      r12_13 = 88.27d0
      r14_15 = 411.9d0
      r16_17 = 471.4d0
      r16_18 = 2693.0d0
      xh2_ini = 2.781d-5
      xhe3_ini = 3.461d-5
      xli6_ini = 7.187d-10
      xli7_ini = 1.025d-8
      xbe9_ini = 1.595d-10
      xb10_ini = 1.002d-9
      xb11_ini = 4.405d-9
      zxmix = 0.02292d0
      frac_c = 0.172148d0
      frac_n = 0.050426d0
      frac_o = 0.468195
      mixture_id_table = (/ character(len=8) :: 'GS98', 'AAG21', 'M22P', 'M22M' /)
      zx_mix_table = (/ 0.02292d0, 0.0187d0, 0.0225d0, 0.0226d0 /)
      frac_c_table = (/ 0.172148d0, 0.184156d0, 0.19239d0, 0.191425d0 /)
      frac_n_table = (/ 0.050426d0, 0.050344d0, 0.059012d0, 0.058716d0 /)
      frac_o_table = (/ 0.468195d0, 0.416592d0, 0.415631d0, 0.413545d0 /)
      tridt = 1.0d-2
      tridl = 8.0d-2
      dtwind = 1.0d1
      lptime = .true.
      ltrist = .false.
      lkuthe = .false.
      lovste = .false.
      lexcom = .false.
      lrot = .false.
      linstb = .false.
      alfa = 1.5d0
      fk = 1.0d0
      fw = 1.0d0
      fc = 1.0d0
      fes = 1.0d0
      fgsf = 1.0d0
      fmu = 1.0d0
      fss = 1.0d0
      rcrit = 1.0d3
      wmax = 3.0d-4
      wmax_sun = 1000.0
      kttau = 0
      clsun = 3.8515d33
      crsun = 6.9598d10
      lcore = .false.
      mcore = 0
      fcore = 1.0
      lmhd = .false.
      grtol = 1.0d-8
      ilambda = 1
      niter_gs = 10
      ldify = .false.
      dt_gs = 0.1d0
      xmin = 1.0d-3
      ymin = 1.0d-3
      lthoulfit = .false.
      ldifz = .false.
      lnewdif = .false.
      cmin = 1.0d-20
      abstol = 1.0d-5
      reltol = 1.0d-4
      kemmax = 50
      etadh0 = -1.0d0
      etadh1 = 1.0d0
      ldh = .false.
      fesc = 1.0d0
      fssc = 1.0d0
      fgsfc = 1.0d0
      igsf = 1
      lpurez = .false.
      lrwsh = .false.
      liso = .false.
      lnulos1 = .false.
      toll = 1.0d-5
      tolr = 1.0d-4
      tolz = 1.0d-3
      lcals = .false.
      lcalsolzx = .false.
      calsolage = 4.57d9
      calsolzx = 0.02292d0
      lzramp = .false.
      lcalst = .false.
      lteff = .false.
      lopale = .false.
      lopale01 = .false.
      lcondopacp = .false.
      lopale06 = .false.
      lnumderiv = .false.
      lscv = .false.
      lsolid = .false.
      impjmod = 0
      dmdt0 = -1.0d-14
      compacc = (/ 0.71668d0, 0.265721d0, 0.01757d0, 2.9d-5, 3.013d-3, 3.385d-5, 9.346d-4, 0.0d0, &
           8.462d-3, 0.0d0, 1.696d-5, 0.0d0, 2.0d-9, 2.0d-9, 3.0d-11 /)
      lmdot = .false.
      lsolwind = .false.
      tmolmin = 4.0d0
      tmolmax = 4.1d0
      lalex06 = .false.
      lmwind = .false.
      lrossby = .true.
      lbscale = .false.
      awind = 'K97'
      pmma = 2.0
      pmmb = 1.0
      pmmc = 2.0
      pmmd = 0.0
      pmmm = 0.5
      pmmjd = 1.32e30
      pmmmd = 1.27e12
      pmmsolp = 4.9304d0
      pmmsolw = 2.836e-6
      pmmsoltau = 9.4805d5
      lcodm = .false.
      codm = 2.5d4
      xsli6 = 5.5d3
      xsli7 = 5.2d1
      xsbe91 = 1.1d0
      xsbe92 = 1.5d4
      xsbe93 = 1.5d4
      lxli6 = .false.
      lxli7 = .false.
      lxbe91 = .false.
      lxbe92 = .false.
      lxbe93 = .false.
      sli6 = 1.0d0
      sli7 = 1.0d0
      sbe91 = 1.0d0
      sbe92 = 1.0d0
      sbe93 = 1.0d0
      spotf = 0.00
      spotx = 1.00
      lsdepth = .false.
      s0_1_1 = 4.01d-22
      s0_3_3 = 5.21d3
      s0_3_4 = 5.6d-1
      s0_1_12 = 1.34d0
      s0_1_13 = 7.6d0
      s0_1_14 = 1.66d0
      s0_1_16 = 1.06d1
      s0_1_be7e = 1.7709d-10
      s0_1_be7p = 0.0208d0
      s0_1_15_c12alp = 7.3d4
      s0_1_15_o16 = 3.6d1
      s0p_1_1 = 4.49d-24
      s0p_3_3 = -4.9d0
      s0p_3_4 = -3.6d-4
      s0p_1_12 = 2.6d-3
      s0p_1_13 = -7.83d-3
      s0p_1_14 = -3.3d-3
      s0p_1_16 = -5.4d-2
      s0p_1_be7p = -3.12d-5
      s0pp_1_12 = 8.3d-5
      s0pp_1_13 = 7.29d-4
      s0pp_1_16 = 0.0d0
      s0pp_1_be7p = -2.288d-7
      lnewnuc = .false.
      lcalcenv = .true.

      terminal_unit = 6
! OUTPUT: LAST MODEL (TEXT)
      last_model_unit = 11
! INPUT: FIRST MODEL (TEXT)
      first_unit = 12
! INPUT: PHYSICS NAMELIST
      run_unit = 13
! INPUT: CONTROL NAMELIST
      standard_unit = 14
! INPUT: FERMI TABLES
      fermi_unit = 15
! OUTPUT: THE RUN LOG
      run_log_unit = 20
! OUTPUT: FOR PULSATION CODE, INTERIOR
! OUTPUT: FOR PULSATION CODE, ENVELOPE
! OUTPUT: FOR PULSATION CODE, ATMOSPHERE
! OUTPUT: BINARY OUTPUT OF LAST MODEL
      last_model_binary_lu = 27
! OUTPUT: BINARY OUTPUT OF STORED MODELS
      stored_models_binary_lu = 28
! INPUT: BINARY STARTING MODEL
      first_model_binary_lu = 29
! OUTPUT: INFO RELAVENT TO DYNAMO
      dynamics_unit = 30
! YCK INPUT: OPAL92 OPACITY TABLES
      opal92_table_unit = 32
! YCK INPUT: OPAL95 OPACITY TABLE
      opal95_table_unit = 48
! OUTPUT: EXTENDED COMPOSITION INFO
      composition_unit = 34
! YCK INPUT: KURUCZ LOW T OPACITIES
      kurucz_table_unit = 36
! OUPUT: ISOCHRONE INFORMATION
      isochrone_file_unit = 37
! INPUT: KURUCZ ATMOSPHER TABLE
! YCK INPUT: Alex LOW T OPACITIES
      alex94_table_unit = 39
! INPUT: MHD EQU. OF STATE TABLES (2026 wave 2: these eight stay
! assigned here -- the star%ctrl members have no declaration-time
! default, so they would seed as 0, and mhdst treats a
! unit .le. 0 as "table absent")
      unit_zams_a = 40
      unit_zams_b = 41
      unit_zams_c = 42
      unit_centre1 = 43
      unit_centre2 = 44
      unit_centre3 = 45
      unit_centre4 = 46
      unit_centre5 = 47
! INPUT: OPAL EQUATION OF STATE
      iopale = 49
! (2026: the LAOL and pure-Z LAOL table units, formerly fixed
! iolaol = 61 / iopurez = 62, are newunit locals of kap/laol89/
! read_laol_tables.f90 and read_pure_z_table.f90 now.)
! DBG 4/94
!     INPUT:
! DBG 8/95 SECOND OPOACITY TABLES FOR ZRAMP AND Z DIFFUSION
! (2026: the second LAOL/OPAL92 table units are runtime-allocated
! locals of their readers now, not fixed io_units units.)
      kurucz_table2_unit = 65
! MHP 6/97 ADDED OPTION FOR ALLARD MODEL ATMOSPHERES
      allard_table_unit = 66
! MHP 6/98 MONTE CARLO FOR SNUs
      monte_carlo_unit1 = 70
      monte_carlo_unit2 = 71
! INPUT FILES FOR THE SCV EOS
      scv_h_unit=72
      scv_he_unit=73
      scv_z_unit=74
! Input files for Potekhin conductive opacities. LLP 7/8/06
      icondopacp = 75
!      FcondOpacP = 'condall.d'
! Default FeH and Alpha for new Allard Atmospheres
      alatm_feh = 0d0
      alatm_alpha = 0d0
      laltptau100 = .false.
! 3/09 Input file for 2006 Alexander opacities
      alex06_table_unit = 90

      print *,''
      print *,'Yale Rotating Evolution Code - YREC, v',yrec_version_string(1:len_trim(yrec_version_string)),' (',git_hash_string(1:len_trim(git_hash_string)),')'

! JVS 02/11 Altered the yrec8 input format so that files can be entered
! on the command line, with the *.nml1 as the first argument, and *.nml2 as
! the second. Defaults to yrec8.nml1 and yrec8.nml2 if none are provided.
!      OPEN(UNIT=ISTAND, FILE='yrec8.nml1', STATUS='OLD')
!      OPEN(UNIT=IRUN, FILE='yrec8.nml2', STATUS='OLD')
!      READ(UNIT=ISTAND, NML=CONTROL)
!      READ(UNIT=IRUN, NML=PHYSICS)
!      CLOSE(ISTAND)
!      CLOSE(IRUN)

! Dynamically create format string so version info is nicely spaced
! max(1,...): a zero width (A00) is illegal -- belt and braces on top
! of setversion's empty-hash guard
      write(version_fmt, 315) max(1, len_trim(yrec_version_string)), &
           max(1, len_trim(git_hash_string))
      315 format('(''# YREC v'', A', i2.2, ', '' ('', A', i2.2, &
           &        ', '')'')')

! 2026 (libyrec): when embedded (pyyrec), argv belongs to the host
! process, so the C API injects the paths via star_job_lib overrides;
! blank overrides preserve the CLI getarg behavior exactly.
      if (control_nml_override .ne. ' ') then
         control_nml_file = control_nml_override
      else
      call getarg(1, control_nml_file)
      if (control_nml_file(1:2) .eq. ' ') control_nml_file = 'yrec8.nml1'
      end if
      print *, ' '
      call read_namelist_files
      call adopt_canonical_names
      call resolve_output_mode_and_paths
      call derive_options_and_open_files(ierr)
      call store_table_paths
      if (ierr /= 0) return
      call echo_settings
      call interpret_kind_cards
      return

end subroutine read_input


! ---------------------------------------------------------------
! star%job <-> local NAMELIST-spelled table-path buffers. Load runs
! first so a re-entrant run starts from the previous run's stored
! paths (exactly what the former argument association gave); store
! runs after the reads and expand_value resolution.
subroutine load_table_paths
      falex06 = star%job%alex06_table_path
      fallard = star%job%allard_table_path
      fatm = star%job%atm_table_path
      ffermi = star%job%fermi_table_path
      fkur = star%job%kurucz_table_path
      fkur2 = star%job%kurucz_table2_path
      flaol = star%job%laol_table_path
      flaol2 = star%job%laol_table2_path
      fliv95 = star%job%opal95_table_path
      flldat = star%job%opal92_table_path
      fmhd1 = star%job%zams_a_table_path
      fmhd2 = star%job%zams_b_table_path
      fmhd3 = star%job%zams_c_table_path
      fmhd4 = star%job%centre1_table_path
      fmhd5 = star%job%centre2_table_path
      fmhd6 = star%job%centre3_table_path
      fmhd7 = star%job%centre4_table_path
      fmhd8 = star%job%centre5_table_path
      fopal2 = star%job%opal92_table2_path
      fpurez = star%job%pure_z_table_path
      fscvh = star%job%scv_h_table_path
      fscvhe = star%job%scv_he_table_path
      fscvz = star%job%scv_z_table_path
      opecalex = star%job%alex95_table_paths
end subroutine load_table_paths

subroutine store_table_paths
      star%job%alex06_table_path = falex06
      star%job%allard_table_path = fallard
      star%job%atm_table_path = fatm
      star%job%fermi_table_path = ffermi
      star%job%kurucz_table_path = fkur
      star%job%kurucz_table2_path = fkur2
      star%job%laol_table_path = flaol
      star%job%laol_table2_path = flaol2
      star%job%opal95_table_path = fliv95
      star%job%opal92_table_path = flldat
      star%job%zams_a_table_path = fmhd1
      star%job%zams_b_table_path = fmhd2
      star%job%zams_c_table_path = fmhd3
      star%job%centre1_table_path = fmhd4
      star%job%centre2_table_path = fmhd5
      star%job%centre3_table_path = fmhd6
      star%job%centre4_table_path = fmhd7
      star%job%centre5_table_path = fmhd8
      star%job%opal92_table2_path = fopal2
      star%job%pure_z_table_path = fpurez
      star%job%scv_h_table_path = fscvh
      star%job%scv_he_table_path = fscvhe
      star%job%scv_z_table_path = fscvz
      star%job%alex95_table_paths = opecalex
      star%job%conductive_table_path = fcondopacp
end subroutine store_table_paths

! ---------------------------------------------------------------
! Pick the CONTROL/PHYSICS namelist files (CLI override aware) and
! read them -- through the new-style &star_job/&controls reader when
! the control file carries a star_job group, else the legacy
! NAMELIST /control/ + /physics/ pair.
subroutine read_namelist_files
! 2026: detect the inlist style BEFORE printing the banner, so a
! single-inlist run announces its one file instead of a phantom
! legacy PHYSICS namelist (the getarg(2) default used to print
! 'yrec8.nml2' even though the new-style dispatch never reads it).
      logical :: new_style, single_inlist

      new_style = nml_file_has_group(control_nml_file, 'star_job')
      single_inlist = .false.
      if (new_style) &
           single_inlist = nml_file_has_group(control_nml_file, 'controls')

      if (physics_nml_override .ne. ' ') then
         physics_nml_file = physics_nml_override
      else
      call getarg(2, physics_nml_file)
      if (physics_nml_file(1:2) .eq. ' ') physics_nml_file = 'yrec8.nml2'
      end if
! single-inlist style: both groups live in the one file; point the
! physics slot at it too (copy_inlists_used then dedupes it).
      if (single_inlist) physics_nml_file = control_nml_file

      if (single_inlist) then
         write(*,*) 'inlist           :  ',control_nml_file(1:len_trim(control_nml_file))
      else
      write(*,*) 'CONTROL namelist :  ',control_nml_file(1:len_trim(control_nml_file))
      write(*,*) 'PHYSICS namelist :  ',physics_nml_file(1:len_trim(physics_nml_file))
      end if

! Defaults for the two output paths every run needs: the final
! model file (.mod) and the run log -- whose directory part also
! defines the output directory. Both use {YREC_OUTPUT}, expanded by
! expand_value below ($YREC_OUTPUT, or "output" when unset); a deck
! that sets FLAST/last_model_file or FSHORT/short_output_file
! overrides them. The other f* path controls remain deck-supplied.
      flast = '{YREC_OUTPUT}/final.mod'
      fshort = '{YREC_OUTPUT}/run.log'

! 2026 inlist revamp: dispatch on inlist style. New-style files carry
! &star_job (+ &controls, same file or the second); everything else
! takes the legacy reader path.
      if (new_style) then
      include 'inlist_new_read.inc'
      else
      open(unit=standard_unit, file=control_nml_file, status='OLD')
      open(unit=run_unit, file=physics_nml_file, status='OLD')
      read(unit=standard_unit, nml=control)
      read(unit=run_unit, nml=physics)
      close(standard_unit)
      close(run_unit)
      end if
end subroutine read_namelist_files

! ---------------------------------------------------------------
! Copy the namelist-spelled locals into their canonical const_lib
! homes (see the NAMING NOTE in the header: NAMELIST members must
! keep their original spelling, so renamed configuration lives in
! module variables and is assigned here after the read).
subroutine adopt_canonical_names
! stolr0/imax/nuse must keep their exact NAMELIST /physics/ spelling
! (see this file's naming note at the top), so integration_parameters's
! canonically-named variables are set by copying from them here,
! after the namelist read, rather than by renaming in place.
      tolerance_fraction = stolr0
      max_stage_index = imax
      extrap_order = nuse
! tscut must likewise keep its NAMELIST spelling; copy into const_lib's
! saha_log10t_cutoff here.
      saha_log10t_cutoff = tscut
! use_new_nuclear_rates/weak_screening_threshold: same reasoning,
! copied from their NAMELIST-spelled locals before remap runs, since
! map_user_inputs.f90 reads use_new_nuclear_rates to decide how to compute
! const_lib's cross-section-scale members.
      use_new_nuclear_rates = lnewnuc
      weak_screening_threshold = weakscreening
! 2026 phase-C rename campaign: copy the NAMELIST-spelled locals
! into their canonically-named controls_buffer buffers (see the local
! declaration block above).
      overshoot_alpha_core = alphac
      overshoot_alpha_envelope = alphae
      diffusion_timestep_factor = dtdif
      num_rotation_structure_iters = itdif1
      max_diffusion_iters = itdif2
      diffuse_lithium = ldifli
      calc_neutrinos = lsnu
      core_overshoot_active = lovstc
      set_initial_omega = lwnew
      use_wind_torque = ljdot0
      use_thoul_diffusion = lthoul
      use_semiconvection = lsemic
      improved_first_guess_flag = lnews
      max_iter_level1 = niter1
      max_iter_level2 = niter2
      max_iter_level3 = niter3
      max_iter_level4 = niter4
      nuclear_logT_cutoffs = tcut
      mixture_weights_seed = vnew
      initial_omega = wnew
      monte_carlo_active = lmonte
      mc_run_start = imbeg
      mc_run_end = imend
      reimers_scaling_factor = creim
      use_reimers_wind = lreimer
      constant_mixing_coeff = cstmixing
      constant_settling_reduction = cstdiffmix
! envelope_overshoot_active/rotation_active/instability_transport_
! active/mass_accretion_rate/accreted_composition/use_mass_accretion:
! same reasoning, copied from their NAMELIST-spelled locals. Order
! relative to `call remap` doesn't matter for these -- map_user_inputs.f90
! doesn't read any of them.
      envelope_overshoot_active = lovste
      rotation_active = lrot
      instability_transport_active = linstb
      mass_accretion_rate = dmdt0
      accreted_composition = compacc
      use_mass_accretion = lmdot
! use_itoh_neutrino_loss/calc_envelope_flag: same reasoning, copied
! from their NAMELIST-spelled locals.
      use_itoh_neutrino_loss = lnulos1
      calc_envelope_flag = lcalcenv
! clsun/crsun must likewise keep their NAMELIST spelling; copy into
! const_lib's solar_luminosity_cgs/solar_radius_cgs here (before
! setup/setups.f90 computes the rest of former common/const/ from
! these two).
      star%solar_luminosity_cgs = clsun
      star%solar_radius_cgs = crsun
! lexcom must likewise keep its NAMELIST spelling; copy into
! const_lib's use_extended_composition here.
      use_extended_composition = lexcom
! spotf/spotx/lsdepth, tdisk/pdisk/ldisk, endage/setdt/end_dcen/
! end_xcen/end_ycen, lmhd, and optol must likewise keep their NAMELIST
! spelling; copy into const_lib's canonical names here.
      spot_filling_factor = spotf
      spot_temp_contrast = spotx
      spot_depth_varies = lsdepth
      disk_locking_age_gyr = tdisk
      disk_omega_rad_s = pdisk
      disk_locking_active = ldisk
      target_end_age = endage
      timestep_override = setdt
      central_deuterium_stop = end_dcen
      central_hydrogen_stop = end_xcen
      central_helium_stop = end_ycen
      use_mhd_eos = lmhd
      metal_fraction_match_tolerance = optol
      convergence_tolerance = djok
      settling_tolerance = grtol
      coulomb_log_choice = ilambda
      settling_num_iterations = niter_gs
      diffuse_helium_active = ldify
      settling_timestep_fraction = dt_gs
      hydrogen_diffusion_floor = xmin
      helium_diffusion_min = ymin
      use_thoul_fit = lthoulfit
      use_diffusion_z = ldifz
      use_new_diffusion_routines = lnewdif
      atm_error_tol = atmerr
      atm_step_initial = atmd0
      atm_step_begin = atmbeg
      atm_step_min = atmmin
      atm_step_max = atmmax
      env_error_tol = enverr
      env_step_begin = envbeg
      env_step_min = envmin
      env_step_max = envmax
      atm_choice = kttau
      debye_huckel_eta_min = etadh0
      debye_huckel_eta_max = etadh1
      use_debye_huckel_correction = ldh
      monte_carlo_file1_path = fmonte1
      monte_carlo_file2_path = fmonte2
      num_models = nmodls
      first_call_flag = lfirst
      num_runs = numrun
      max_domega_global = dtwind
      use_structure_dt_limits = lptime
      tri_delta_teffl = tridt
      tri_delta_logl = tridl
      atm_step_size = atmstp
      envelope_step_size = envstp
      envelope_generation_flag = lenvg
      helium_flash_active = lkuthe
      new_species_value = xnewcp
      rescale_species_active = lnewcp
      mixture_change_mode = isetmix
      isotope_change_mode = isetiso
      target_carbon_cno_fraction = frac_c
      target_nitrogen_cno_fraction = frac_n
      target_oxygen_cno_fraction = frac_o
      c12_to_c13_ratio = r12_13
      n14_to_n15_ratio = r14_15
      o16_to_o17_ratio = r16_17
      o16_to_o18_ratio = r16_18
      target_metal_fraction = zxmix
      initial_h2_fraction = xh2_ini
      initial_he3_fraction = xhe3_ini
      initial_li6_fraction = xli6_ini
      initial_li7_fraction = xli7_ini
      initial_be9_fraction = xbe9_ini
      initial_b10_fraction = xb10_ini
      initial_b11_fraction = xb11_ini
      difad_velocity_scale = fw
      mixing_velocity_scale = fc
      es_velocity_scale = fes
      gsf_velocity_scale = fgsf
      mu_gradient_scale = fmu
      secular_shear_velocity_scale = fss
      critical_reynolds = rcrit
      es_mixing_scale = fesc
      secular_shear_mixing_scale = fssc
      gsf_mixing_scale = fgsfc
      gsf_inhibition_mode = igsf
      min_abundance = cmin
      absolute_tolerance = abstol
      relative_tolerance = reltol
      max_burn_iterations = kemmax
      extend_core_inward = lcore
      num_core_shells_added = mcore
      core_mass_reduction_factor = fcore
      use_pure_z_table = lpurez
      rewind_short_file = lrwsh
      isochrone_output_active = liso
      initial_x_array = xenv0a
      initial_z_array = zenv0a
      mixing_length_array = cmixla
      has_senv0_array = lsenv0a
      senv0_array = senv0a
      luminosity_tolerance = toll
      radius_tolerance = tolr
      zx_tolerance = tolz
      calibrate_solar_model = lcals
      calibrate_solar_zx = lcalsolzx
      target_solar_zx = calsolzx
      target_solar_age = calsolage
      use_z_ramp = lzramp
      target_luminosity_lsun = xls
      target_star_luminosity_tolerance = xlstol
      target_teff = steff
      target_radius_rsun = sr
      specify_teff_flag = lteff
      calibrate_star_flag = lcalst
      use_opal95_eos = lopale
      use_opal2001_eos = lopale01
      use_opal2006_eos = lopale06
      use_numerical_derivatives = lnumderiv
      laol_table_z1 = zlaol1
      laol_table_z2 = zlaol2
      opal_table_z1 = zopal1
      opal_table_z2 = zopal2
      opal95_single_table_z = zopal951
      alex_table_z1 = zalex1
      kurucz_table_z1 = zkur1
      kurucz_table_z2 = zkur2
      molecular_opacity_logt_min = tmolmin
      molecular_opacity_logt_max = tmolmax
      use_alex06_tables = lalex06
      use_laol89_tables = llaol89
      use_opal92_tables = lopal92
      use_opal95_tables = lopal95
      use_kurucz90_tables = lkur90
      use_alex95_tables = lalex95
      use_conductive_opacity = lcondopacp
      use_diffusion_advection_transport = ldifad
      no_am_transport_in_core = lnoj
      allard_target_feh = alatm_feh
      allard_target_alpha = alatm_alpha
      allard_use_tau100 = laltptau100
      force_solid_body_rotation = lsolid
      solid_body_mode_flag = impjmod
      use_envelope_triangle_dt = ltrist
      pmm_exponent_a = pmma
      pmm_exponent_b = pmmb
      pmm_exponent_c = pmmc
      pmm_exponent_d = pmmd
      pmm_exponent_m = pmmm
      pmm_norm_jdot = pmmjd
      pmm_norm_mdot = pmmmd
      pmm_solar_pressure = pmmsolp
      pmm_solar_omega = pmmsolw
      pmm_solar_turnover_timescale = pmmsoltau
      use_pmm_wind_law = lmwind
      scale_by_rossby_number = lrossby
      scale_by_b_field = lbscale
      wind_law_name = awind
      wind_saturation_omega = wmax
      constant_background_diffusion_coeff = codm
      use_constant_background_diffusion = lcodm
      s0_pp = s0_1_1
      s0_he3he3 = s0_3_3
      s0_he3he4 = s0_3_4
      s0_p_c12 = s0_1_12
      s0_p_c13 = s0_1_13
      s0_p_n14 = s0_1_14
      s0_p_o16 = s0_1_16
      s0_be7_electron = s0_1_be7e
      s0_be7_p = s0_1_be7p
      s0_n15_p_c12_branch = s0_1_15_c12alp
      s0_n15_p_o16_branch = s0_1_15_o16
      s0p_pp = s0p_1_1
      s0p_he3he3 = s0p_3_3
      s0p_he3he4 = s0p_3_4
      s0p_p_c12 = s0p_1_12
      s0p_p_c13 = s0p_1_13
      s0p_p_n14 = s0p_1_14
      s0p_p_o16 = s0p_1_16
      s0pp_p_c12 = s0pp_1_12
      s0pp_p_c13 = s0pp_1_13
      s0pp_p_o16 = s0pp_1_16
      s0p_be7_p = s0p_1_be7p
      s0pp_be7_p = s0pp_1_be7p
      use_scv_eos = lscv
! LNEWNUC = .false.: the cross sections are the SSTANDARD ratios (to
! Bahcall & Ulrich 1988), as in the F77 code's COMMON /CROSS/;
! map_user_inputs replaces them from the s0_* inputs when LNEWNUC is
! set. (2026-10: the 2026 rename campaign dropped this copy as "never
! referenced", which left every pp/CNO rate zero -- and the Be7
! branching 0/0 -- whenever LNEWNUC was .false.; every shipped deck
! sets it .true., so no pin saw it. SSTANDARD is settable from the
! legacy namelist; &controls inlists get the default ratios.)
      if (.not. lnewnuc) star%cross_section_scale = sstandard
! MHP 8/14 SUBROUTINE TO CONVERT MORE USER-FRIENDLY INPUT VARIABLES
! INTO THE VECTORS USED IN THE CODE (SUPERCEDES OLDER INPUTS)
! (2026 phase C: remap is read-path like this routine -- it reads
! AND writes the controls BUFFER; the stores below/in read_controls
! carry its results into star%ctrl.)
      call map_user_inputs
! 2026 phase C: output_init_mesa (called from
! derive_options_and_open_files below) reads its output controls
! from star%ctrl, but read_controls' final store has not happened
! yet at that point -- store the buffer into star%ctrl here, after
! the namelist reads, the local->canonical copies, and remap.
! read_controls' final store after parmin returns then also captures
! echo_settings' clamps (walpcz etc.); each store is an idempotent
! buffer copy.
      call store_controls_to_star
! MHP 06/13 Added memory of whether the choice of atmospheres has
! been changed during the run, and what the original setting was
      star%atm_choice_initial = kttau
      star%use_ttau_relation = .false.
end subroutine adopt_canonical_names






! Does this namelist file contain the group &<gname>? Used to detect
! new-style inlists (and whether &controls shares the &star_job file).
logical function nml_file_has_group(fname, gname)
      character(len=*), intent(in) :: fname, gname
      character(len=512) :: probe_line
      integer :: probe_unit, ios

      nml_file_has_group = .false.
      open(newunit=probe_unit, file=fname, status='OLD', &
           action='READ', iostat=ios)
      if (ios /= 0) return
      do
         read(probe_unit, '(a)', iostat=ios) probe_line
         if (ios /= 0) exit
         if (index(adjustl(probe_line), '&' // trim(gname)) == 1) then
            nml_file_has_group = .true.
            exit
         end if
      end do
      close(probe_unit)
end function nml_file_has_group

end module controls_reader
