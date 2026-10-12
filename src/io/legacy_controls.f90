!----------------------------------------------------------------------
! legacy_controls
!----------------------------------------------------------------------
! The controls reader's working state: the variables the legacy
! CONTROL and PHYSICS namelists bind to, under their original
! spellings (a namelist matches input entries by variable name, so
! renaming one would break every deck that sets it), the two namelist
! groups, the declarations of the new-style &star_job/&controls
! variables (inlist_new_decl.inc, generated), and the reader's other
! working variables. read_input copies the values into their
! canonical homes (controls_buffer, star%job) after the read;
! reset_legacy_controls clears everything at the start of each read.
module legacy_controls
      use controls_buffer
      use integration_parameters
      implicit none

! The legacy namelists' binding targets, under their original spellings.
! The legacy namelists are frozen: new controls go only into the new
! format (controls_buffer and the registry), so nothing is added here.
! After the read, adopt_canonical_names copies each into its canonical
! home (a controls_buffer variable or star%job); a control spelled the
! same in both reads straight into controls_buffer and needs no copy.
! defaults/controls_registry.tsv maps every legacy name to its
! canonical and new-inlist names. hpttol, ies and imu are read and then
! unused. version_fmt, element_id, the mixture tables, the inlist file
! names are the reader's own working variables.
      character(len=256) :: version_fmt
      character(len=256) :: flaol, fpurez
      character(len=256) :: flaol2, fopal2, fkur2
      character(len=3) :: anewcp, atmp
      character(len=8) :: amix, aiso
      character(len=256) :: fiso
      character(len=256) :: fatm
      character(len=256) :: fstch
      character(len=256) :: fallard, fscvh, fscvhe, fscvz
      character(len=256) :: flast, ffirst, ffermi, &
           fshort, &
           fdyn, &
           flldat, fscomp, fkur, &
           fmhd1, fmhd2, fmhd3, fmhd4, fmhd5, fmhd6, fmhd7, fmhd8
      integer :: kindrn(50)
      double precision :: rsclm(50), rsclx(50), rsclz(50), rsclcm(50)
      character(len=256) :: fliv95
      character(len=256) :: fopale, fopale01, fcondopacp, fopale06
      character(len=256) :: opecalex(7)
      character(len=256) :: falex06

      double precision :: alfa, fk, wmax_sun, pmmwmax, zalex2, zopal952
      logical :: lpmm, lsolwind

      character(len=3) :: element_id(12)
      character(len=8) :: mixture_id_table(4)
      double precision :: zx_mix_table(4), frac_c_table(4), &
           frac_n_table(4), frac_o_table(4)
      character(len=256) :: control_nml_file, physics_nml_file

      character(len=256) :: fmonte1, fmonte2

      character(len=256) :: descrip(2)

      double precision :: tridt, tridl

      integer :: nmodls(50), numrun
      logical :: lfirst(50)

      double precision :: clsun, crsun

      double precision :: dtwind

      logical :: lptime

      double precision :: hpttol(12)

      double precision :: djok

      logical :: lovste

      double precision :: atmstp, envstp
      logical :: lenvg

      logical :: lexcom

      logical :: lkuthe

      double precision :: atmerr, atmd0, atmbeg, atmmin, atmmax

      double precision :: enverr, envbeg, envmin, envmax

      double precision :: stolr0
      integer :: imax, nuse

      double precision :: tscut

      double precision :: xnewcp
      logical :: lnewcp

      integer :: isetmix, isetiso
      double precision :: frac_c, frac_n, frac_o, r12_13, r14_15, r16_17, r16_18, zxmix, xh2_ini, &
           xhe3_ini, xli6_ini, xli7_ini, xbe9_ini, xb10_ini, xb11_ini

      double precision :: optol

      logical :: lrot, linstb

      double precision :: endage(50), setdt(50), end_dcen(50), end_xcen(50), end_ycen(50)

      double precision :: fw, fc, fes, fgsf, fmu, fss, rcrit

      double precision :: etadh0, etadh1
      logical :: ldh

      double precision :: fesc, fssc, fgsfc
      integer :: igsf

      double precision :: grtol
      integer :: ilambda, niter_gs
      logical :: ldify

      double precision :: cmin, abstol, reltol
      integer :: kemmax

      double precision :: dt_gs, xmin, ymin
      logical :: lthoulfit

      logical :: ldifz

      logical :: lnewdif

      integer :: kttau

      logical :: lmhd

      logical :: lcore
      integer :: mcore
      double precision :: fcore

      logical :: lpurez

      logical :: lrwsh, liso

      double precision :: xenv0a(50), zenv0a(50), cmixla(50), senv0a(50)
      logical :: lsenv0a(50)

      logical :: lnulos1

      double precision :: toll, tolr, tolz, calsolzx, calsolage
      logical :: lcals, lcalsolzx

      logical :: lzramp

      double precision :: xls, xlstol, steff, sr
      logical :: lteff, lcalst

      logical :: lopale, lopale01, lopale06, lnumderiv

      double precision :: zlaol1, zlaol2, zopal1, zopal2, zopal951, zalex1, zkur1, zkur2, tmolmin, &
           tmolmax
      logical :: lalex06, llaol89, lopal92, lopal95, lkur90, lalex95

      logical :: lcondopacp

      logical :: ldifad

      logical :: lnoj

      double precision :: sstandard(17)
      logical :: lnewnuc

      double precision :: alphac, alphae
      double precision :: dtdif
      double precision :: wnew, creim
      double precision :: cstmixing, cstdiffmix
      double precision :: tcut(5)
      double precision :: vnew(12)
      integer :: itdif1, itdif2
      integer :: niter1, niter2, niter3, niter4
      logical :: ldifli, lsnu, lovstc
      logical :: lwnew, ljdot0, lthoul
      logical :: lsemic, lnews, lmonte
      logical :: lreimer
      integer :: imbeg, imend

      double precision :: s0_1_1, s0_3_3, s0_3_4, s0_1_12, s0_1_13, s0_1_14, s0_1_16, &
           s0_1_be7e, s0_1_be7p, s0_1_15_c12alp, s0_1_15_o16, s0p_1_1, s0p_3_3, s0p_3_4, &
           s0p_1_12, s0p_1_13, s0p_1_14, s0p_1_16, s0pp_1_12, s0pp_1_13, s0pp_1_16, s0p_1_be7p, &
           s0pp_1_be7p

      logical :: lscv

      double precision :: alatm_feh, alatm_alpha
      logical :: laltptau100

      double precision :: tdisk, pdisk
      logical :: ldisk

      double precision :: weakscreening

      logical :: lsolid
      integer :: impjmod

      double precision :: dmdt0, compacc(15)
      logical :: lmdot

      logical :: ltrist

      double precision :: pmma, pmmb, pmmc, pmmd, pmmm, pmmjd, pmmmd, pmmsolp, pmmsolw, pmmsoltau
      logical :: lmwind, lrossby, lbscale
      character(len=3) :: awind

      double precision :: wmax

      logical :: lcalcenv

      double precision :: codm
      logical :: lcodm

      double precision :: xsli6, xsli7, xsbe91, xsbe92, xsbe93
      logical :: lxli6, lxli7, lxbe91, lxbe92, lxbe93

      double precision :: sli6, sli7, sbe91, sbe92, sbe93

      double precision :: spotf, spotx
      logical :: lsdepth


      namelist /control/ &
           &    cmixla, calsolage, calsolzx, &
           &    descrip, &
           &    endage, &
           &    flaol, fpurez,flaol2, fopal2, &
           &    flast, ffirst, ffermi, fshort, fstch, &
           &    fdyn, flldat, fscomp, fkur, fmhd1, &
           &    fmhd2, fmhd3, fmhd4, fmhd5, fmhd6, fmhd7, fmhd8, fiso, fatm, &
           &    fkur2, fallard, fscvh, fscvhe, fscvz, fopale, fliv95, &
           &    fmonte1,fmonte2, &
           &    kindrn, &
           &    lfirst, set_initial_model_number, &
           &    terminal_interval, report_solver_diagnostics, &
           &    log_L_upper_limit, log_L_lower_limit, &
           &    Teff_upper_limit, Teff_lower_limit, &
           &    log_g_upper_limit, log_g_lower_limit, &
           &    nu_max_upper_limit, nu_max_lower_limit, &
           &    nu_max_sun, delta_nu_sun, Teff_sun, &
           &    inlist_used_file, profile_data_prefix, &
           &    lzramp, lteff, lcalst, lpurez, &
! MHP 9/24 add LCALSOLZX to namelist
           &    liso, lrwsh, lsenv0a,lcals,lcalsolzx, &
           &    llaol89,lopal92,lopal95,lkur90,lalex95, &
           &    npenv, numrun, nmodls, &
           &    opecalex, &
           &    rsclm, rsclx, rsclz, rsclcm, rsclzc, rsclzm1, rsclzm2, &
           &    setdt, senv0a,steff,sr, &
           &    tolr, toll,tolz, &
           &    xenv0a, xls, xlstol, &
           &    zenv0a, &
           &    zlaol1,zlaol2,zopal1,zopal2, zopal951, &
           &    zopal952, zalex1, zalex2, zkur1, zkur2, &
           &      fopale01,fcondopacp,fopale06,falex06,lalex06, &
! MHP 10/24 ADDED END_DCEN,END_XCEN,END_YCEN VECTORS TO NML1, USED IN MAIN
! MHP 10/24 ADDED HEAVY ELEMENT MIXTURE CONTROLS TO NML1,USED IN STARIN
           &  end_dcen,end_xcen,end_ycen,isetmix,isetiso, &
           &  amix,aiso,frac_c,frac_n,frac_o,r12_13,r14_15,r16_17,r16_18,zxmix, &
           &  xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini,xb10_ini,xb11_ini, &
! new (2026): GYRE-format periodic pulsation output interval, additive
! only -- absent from every existing *.nml1 file, so it simply keeps
! its default (off) there.
           &  pulse_gyre_interval
!
      namelist /physics/ &
           &    atmmin, atmbeg, atmerr, atmmax, atmd0, anewcp, atmp, acfpft, &
           &    atime, alphac, alphae, alfa, alpham, atmstp, abstol, betac, &
           &    cmin, clsun, crsun, solar_mass_cgs, G_cgs, &
           &    dpenv, dtdif, dtwind, djok, dt_gs, &
           &    enverr, envmax, envmin, envbeg, envstp,etadh0,etadh1, &
           &    fcorr0, fcorri, fk,  fw, fc, fo, fmu, fes, &
           &    fcore, fgsf, fss, fesc, fssc, fgsfc, fgry, fgrz, &
           &    grtol, &
           &    htoler, hpttol, &
           &    itfp1, itfp2, imax, itdif1, itdif2, ies, igsf, imu, ilambda, &
           &    kttau, kemmax, &
           &    lvfc, ldifad, lnoj, lnewdif, ldify, ldifz, ldifli, lsnu, ldh, &
           &    lnewcp, lkuthe, lovstc, lovste, lovstm, lovmax, &
           &    lexcom, lrot, lnew0, linstb, lwnew, ljdot0, lptime,ladov,ltrist, &
           &    lenvg, lnulos1, lthoul, lthoulfit, &
           &    lopale, lmhd, lcore, lsemic, lnews, &
           &    mcore, &
           &    niter1, niter2, niter3, niter4, nuse, niter_gs, &
           &    optol, &
           &    rcrit, reltol, &
           &    stolr0, &
           &    tcut, tscut, tenv0, tenv1, tgcut, tridt, tridl, &
           &    tollaol, &
           &    vnew, &
           &    walpcz, wnew, weakscreening, &
           &    xnewcp, xmin, &
           &    ymin, tmolmin,tmolmax, &
           &    lmonte,imbeg,imend,sstandard,lscv, &
           &    ldisk,tdisk,pdisk,wmax,lsolid,impjmod,  &  !JNT 09/2025 FOR 05/15
           &    dmdt0,fczdmdt,ftotdmdt,compacc,creim,lreimer,lmdot, &
           &    lopale01,lcondopacp,lopale06,lnumderiv, &
           &    alatm_feh,alatm_alpha,laltptau100, &  ! For new Allard Atmospheres
           &    cstmixing, cstdiffmix,       &  !CFD oct2009 To mimic mixing(reduce settling)
           &    lsolwind,lmwind,lrossby,lpmm,lbscale, &
           &    awind,pmma,pmmb,pmmc,pmmd,pmmm,pmmjd,pmmmd,pmmwmax, &
! MHP 8/17 ADDED WMAX_SUN
           &    pmmsolp,pmmsolw,pmmsoltau,lcodm,codm,wmax_sun, &
! G Somers 6/14
           &    xsli6,xsli7,xsbe91,xsbe92,xsbe93, &
           &    lxli6,lxli7,lxbe91,lxbe92,lxbe93, &
           &    sli6,sli7,sbe91,sbe92,sbe93,lnewnuc, &
           &    spotf, spotx, lsdepth, &
! G Somers END
! MHP 09/14 ADDED CROSS SECTIONS
           &    s0_1_1,s0_3_3,s0_3_4,s0_1_12,s0_1_13,s0_1_14,s0_1_16, &
           &    s0_pep,s0_1_be7e,s0_1_be7p,s0_hep,s0_1_15_c12alp,s0_1_15_o16, &
           &    s0p_1_1,s0p_3_3,s0p_3_4,s0p_1_12,s0p_1_13,s0p_1_14, &
           &    s0p_1_16,s0pp_1_12,s0pp_1_13,s0pp_1_16,s0p_1_be7p,s0pp_1_be7p, &
           &    flag_dx,flag_dw,flag_dz,lstruct_time, &
           &    time_core_min,time_dl,time_dp,time_dr,time_dt,time_dw_global, &
           &    time_dw_mix,time_dx_core_frac,time_dx_core_tot,time_dx_shell, &
           &    time_dx_total,time_dy_core_frac,time_dy_core_tot,time_dy_shell, &
           &    time_dy_total,tol_czbase_fine_width,tol_dl_max,tol_dm_max, &
           &    tol_dm_min,tol_dp_core_max,tol_dp_czbase_max,tol_dp_env_max, &
           &    tol_dx_max,tol_dz_max,time_max_dt_frac,lnewvars, &
           &    lcalcenv

! 2026 inlist revamp: new-style (&star_job/&controls) machinery,
! generated from defaults/controls_registry.tsv.
      include 'inlist_new_decl.inc'


contains

! ---------------------------------------------------------------
! Every variable above, set at the start of each read: a namelist read
! changes only what the deck sets, and an in-process rerun (pyyrec,
! test_reentry) must start from the same state. Each default is stated
! once: a variable copied to a canonical control that has its default
! in controls_buffer takes it from there (alphac = overshoot_alpha_core;
! the buffer holds the pristine star%ctrl at this point), the others
! have theirs here. A variable without a default starts at zero,
! .false. or NUL bytes, as when these were the reader's locals under
! -finit-local-zero.
subroutine set_legacy_defaults
      alphac = overshoot_alpha_core
      alphae = overshoot_alpha_envelope
      dtdif = diffusion_timestep_factor
      wnew = initial_omega
      creim = reimers_scaling_factor
      cstmixing = constant_mixing_coeff
      cstdiffmix = constant_settling_reduction
      tcut = nuclear_logT_cutoffs
      tscut = saha_log10t_cutoff
      vnew = mixture_weights_seed
      itdif1 = num_rotation_structure_iters
      itdif2 = max_diffusion_iters
      niter1 = max_iter_level1
      niter2 = max_iter_level2
      niter3 = max_iter_level3
      niter4 = max_iter_level4
      ldifli = diffuse_lithium
      lsnu = calc_neutrinos
      lovstc = core_overshoot_active
      lwnew = set_initial_omega
      ljdot0 = use_wind_torque
      lthoul = use_thoul_diffusion
      lsemic = use_semiconvection
      lnews = improved_first_guess_flag
      lmonte = monte_carlo_active
      lreimer = use_reimers_wind
      imbeg = mc_run_start
      imend = mc_run_end
      weakscreening = 0.03d0
! sstandard: the cross sections relative to Bahcall & Ulrich 1988, used
! when lnewnuc is .false.
      sstandard = (/ 0.9828, 1.0485, 0.9815, 0.9241, 1.3818, 1.0542, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &
           1.0, 1.0, 1.0108, 0.7819, 0.2875 /)
      ldifad = .false.
      lnoj = .false.
      tdisk = disk_locking_age_gyr
      pdisk = disk_omega_rad_s
      ldisk = disk_locking_active
      lsenv0a = .false.
      senv0a = 1.26d-4
      lenvg = .false.
      atmstp = 0.5
      envstp = 0.5
      numrun = 1
      kindrn = 1
      lfirst = .true.
      nmodls = 0
      endage = target_end_age
      setdt = timestep_override
      rsclm = 0.0
      rsclx = 0.0
      rsclz = 0.0
      rsclcm = 0.0
      end_dcen = central_deuterium_stop
      end_xcen = central_hydrogen_stop
      end_ycen = central_helium_stop
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
! Rescaling one species of the starting model (anewcp names it)
      lnewcp = .false.
      anewcp = '   '
      xnewcp = 1.3d1
! The starting model's mixture: isetmix = 1 sets the C, N and O fractions
! of Z from amix ('GS98', 'AAG21', 'M22P', 'M22M', or 'CUS' for frac_c/n/o
! and zxmix as given); isetiso = 1 sets the CNO isotope ratios and the
! light-element abundances from aiso ('L21', Lodders et al. 2021, or
! 'CUS'). The table rows: GS98 (Grevesse & Sauval 1998), AAG21 (Asplund
! et al. 2021), M22P and M22M (Magg et al. 2022, photospheric and
! meteoritic).
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
      lovste = envelope_overshoot_active
      lexcom = use_extended_composition
      lrot = rotation_active
      linstb = instability_transport_active
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
      lmhd = use_mhd_eos
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
      lpurez = use_pure_z_table
      lrwsh = .false.
      liso = .false.
      lnulos1 = use_itoh_neutrino_loss
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
      dmdt0 = mass_accretion_rate
      compacc = accreted_composition
      lmdot = use_mass_accretion
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
! Light-element cross sections in keV b (used when lxli6 etc are set):
! Li6 5500 (Fowler et al. 1967), Li7 52 (Rolfs & Kavanagh 1986),
! Be9(p,g) 1.1 (approximate, inferred from Fowler 1988), Be9(p,d) and
! Be9(p,a) 15000 (Fowler et al. 1967)
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
      spotf = spot_filling_factor
      spotx = spot_temp_contrast
      lsdepth = spot_depth_varies
! S-factors (keV b) and their first and second derivatives (b), from
! Solar Fusion II (Adelberger et al. 2011). pep is relative to pp, and
! Be7+e- is the constant of the linear term; the code uses T9, so a
! coefficient for T/1e6 K is divided by sqrt(1000).
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
      lcalcenv = calc_envelope_flag
! the final model and the run log, whose directory is the output
! directory ({YREC_OUTPUT} is $YREC_OUTPUT, or output when unset)
      flast = '{YREC_OUTPUT}/final.mod'
      fshort = '{YREC_OUTPUT}/run.log'
      alatm_feh = 0d0
      alatm_alpha = 0d0
      laltptau100 = .false.
! no default
      version_fmt = repeat(achar(0), len(version_fmt))
      flaol = repeat(achar(0), len(flaol))
      fpurez = repeat(achar(0), len(fpurez))
      flaol2 = repeat(achar(0), len(flaol2))
      fopal2 = repeat(achar(0), len(fopal2))
      fkur2 = repeat(achar(0), len(fkur2))
      atmp = repeat(achar(0), len(atmp))
      fiso = repeat(achar(0), len(fiso))
      fatm = repeat(achar(0), len(fatm))
      fstch = repeat(achar(0), len(fstch))
      fallard = repeat(achar(0), len(fallard))
      fscvh = repeat(achar(0), len(fscvh))
      fscvhe = repeat(achar(0), len(fscvhe))
      fscvz = repeat(achar(0), len(fscvz))
      ffirst = repeat(achar(0), len(ffirst))
      ffermi = repeat(achar(0), len(ffermi))
      fdyn = repeat(achar(0), len(fdyn))
      flldat = repeat(achar(0), len(flldat))
      fscomp = repeat(achar(0), len(fscomp))
      fkur = repeat(achar(0), len(fkur))
      fmhd1 = repeat(achar(0), len(fmhd1))
      fmhd2 = repeat(achar(0), len(fmhd2))
      fmhd3 = repeat(achar(0), len(fmhd3))
      fmhd4 = repeat(achar(0), len(fmhd4))
      fmhd5 = repeat(achar(0), len(fmhd5))
      fmhd6 = repeat(achar(0), len(fmhd6))
      fmhd7 = repeat(achar(0), len(fmhd7))
      fmhd8 = repeat(achar(0), len(fmhd8))
      fliv95 = repeat(achar(0), len(fliv95))
      fopale = repeat(achar(0), len(fopale))
      fopale01 = repeat(achar(0), len(fopale01))
      fcondopacp = repeat(achar(0), len(fcondopacp))
      fopale06 = repeat(achar(0), len(fopale06))
      opecalex = repeat(achar(0), len(opecalex))
      falex06 = repeat(achar(0), len(falex06))
      pmmwmax = 0
      zalex2 = 0
      zopal952 = 0
      lpmm = .false.
      control_nml_file = repeat(achar(0), len(control_nml_file))
      physics_nml_file = repeat(achar(0), len(physics_nml_file))
      fmonte1 = repeat(achar(0), len(fmonte1))
      fmonte2 = repeat(achar(0), len(fmonte2))
      descrip = repeat(achar(0), len(descrip))
      xenv0a = 0
      zenv0a = 0
      cmixla = 0
      xls = 0
      xlstol = 0
      steff = 0
      sr = 0
      zlaol1 = 0
      zlaol2 = 0
      zopal1 = 0
      zopal2 = 0
      zopal951 = 0
      zalex1 = 0
      zkur1 = 0
      zkur2 = 0
      llaol89 = .false.
      lopal92 = .false.
      lopal95 = .false.
      lkur90 = .false.
      lalex95 = .false.
end subroutine set_legacy_defaults

end module legacy_controls
