!----------------------------------------------------------------------
! legacy_controls
!----------------------------------------------------------------------
! The controls reader's working state: the variables the legacy
! CONTROL and PHYSICS namelists bind to, under their original
! spellings (a namelist matches input entries by variable name, so
! renaming one would break every deck that sets it), the two legacy
! namelist groups, and the reader's other working variables. (The
! inlist's &control and &physics groups, under the new names, belong
! to read_inlist in read_controls.) read_input copies the values into
! their canonical homes (controls_buffer, star%job) after the read;
! reset_legacy_controls clears everything at the start of each read.
module legacy_controls
      use controls_buffer
      implicit none

! The legacy namelists' binding targets, under their original spellings.
! The legacy namelists are frozen: new controls go only into the new
! format (controls_buffer and the registry), so nothing is added here.
! After the read, adopt_canonical_names copies each into its canonical
! home (a controls_buffer variable or star%job); a control spelled the
! same in both reads straight into controls_buffer and needs no copy.
! defaults/controls_registry.tsv maps every legacy name to its
! canonical and new-inlist names. ies and imu are read and then
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
           &    solberg_hoiland_velocity_scale, solberg_hoiland_mixing_scale, &
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

contains

! ---------------------------------------------------------------
! Every variable above, set at the start of each read: a namelist read
! changes only what the deck sets, and an in-process rerun (pyyrec,
! test_reentry) must start from the same state. Each default is stated
! once: a variable the read copies to a canonical control takes that
! control's default, its controls_buffer initializer (alphac =
! overshoot_alpha_core; the buffer holds the pristine star%ctrl at this
! point); only the variables without a canonical home have theirs here.
! A variable without a default starts at zero, .false. or NUL bytes, as
! when these were the reader's locals under -finit-local-zero.
subroutine set_legacy_defaults
      include 'legacy_mirror.inc'
! The legacy names without a canonical home, and their defaults.
! sstandard: the cross sections relative to Bahcall & Ulrich 1988, used
! when lnewnuc is .false.
      sstandard = (/ 0.9828, 1.0485, 0.9815, 0.9241, 1.3818, 1.0542, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &
           1.0, 1.0, 1.0108, 0.7819, 0.2875 /)
      kindrn = 1
      rsclm = 0.0
      rsclx = 0.0
      rsclz = 0.0
      rsclcm = 0.0
      element_id = (/ 'HE3', 'C12', 'C13', 'N14', 'N15', 'O16', 'O17', 'O18', 'H2 ', 'LI6', 'LI7', 'BE9' /)
! Rescaling one species of the starting model (anewcp names it)
      anewcp = '   '
! The starting model's mixture: isetmix = 1 sets the C, N and O fractions
! of Z from amix ('GS98', 'AAG21', 'M22P', 'M22M', or 'CUS' for frac_c/n/o
! and zxmix as given); isetiso = 1 sets the CNO isotope ratios and the
! light-element abundances from aiso ('L21', Lodders et al. 2021, or
! 'CUS'). The table rows: GS98 (Grevesse & Sauval 1998), AAG21 (Asplund
! et al. 2021), M22P and M22M (Magg et al. 2022, photospheric and
! meteoritic).
      amix = 'GS98'
      aiso = 'L21'
      mixture_id_table = (/ character(len=8) :: 'GS98', 'AAG21', 'M22P', 'M22M' /)
      zx_mix_table = (/ 0.02292d0, 0.0187d0, 0.0225d0, 0.0226d0 /)
      frac_c_table = (/ 0.172148d0, 0.184156d0, 0.19239d0, 0.191425d0 /)
      frac_n_table = (/ 0.050426d0, 0.050344d0, 0.059012d0, 0.058716d0 /)
      frac_o_table = (/ 0.468195d0, 0.416592d0, 0.415631d0, 0.413545d0 /)
      fmu = 1.0d0
      clsun = 3.8515d33
      crsun = 6.9598d10
      lsolwind = .false.
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
! the final model and the run log, whose directory is the output
! directory ({YREC_OUTPUT} is $YREC_OUTPUT, or output when unset)
      flast = '{YREC_OUTPUT}/final.mod'
      fshort = '{YREC_OUTPUT}/run.log'
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
      descrip = repeat(achar(0), len(descrip))
end subroutine set_legacy_defaults

end module legacy_controls
