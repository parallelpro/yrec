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
subroutine parmin(falex06, fallard, fatm, ffermi, fkur, fkur2, flaol, &
     flaol2, fliv95, flldat, fmhd1, fmhd2, fmhd3, fmhd4, fmhd5, fmhd6, &
     fmhd7, fmhd8, fopal2, fpatm, fpenv, fpmod, fpurez, fscvh, fscvhe, &
     fscvz, opecalex)

      use const_lib
      use luout_lib
      use intpar_lib
      implicit none

! PARAMETERS for tabulated surface pressures (n_atm_teff/n_atm_logg),
! Kurucz/Castelli surface pressures (n_katm_teff/n_katm_logg, JNT
! 6/2014), Allard model surface pressures (n_allard_teff/
! n_allard_logg), the SCV EOS tables (n_scv_teff/n_scv_press), and the
! shared array length used by the variable-FC and acoustic-depth
! diagnostics (max_diag_pts).
      integer, parameter :: n_atm_teff = 57, n_atm_logg = 11
      integer, parameter :: n_katm_teff = 76, n_katm_logg = 11
      integer, parameter :: n_allard_teff = 54, n_allard_logg = 5
      integer, parameter :: n_scv_teff = 63, n_scv_press = 76
      integer, parameter :: max_diag_pts = 5000

! --- CONTROL/PHYSICS namelist variables (including this routine's
!     dummy arguments) and their explicit type overrides (CHARACTER
!     length / array DIMENSION); kept at original spelling -- see
!     NAMING NOTE above. Grouped as in the original file; a few
!     members declared here (OLAOL/OXA/OT/ORHO/TOLLAOL, YRECVER,
!     GITHASH, DESCRIP, FMONTE1/FMONTE2, ICLCD/IACAT/IJLAST/IJVS/
!     IJENT/IJDEL, AWIND) are also COMMON members and get their type
!     declared once, below, alongside their COMMON block instead. ---
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
           fdebug, ftrack, fshort, fmilne, fmodpt, &
           fstor, fpmod, fpenv, fpatm, fdyn, &
           flldat, fsnu, fscomp, fkur, &
           fmhd1, fmhd2, fmhd3, fmhd4, fmhd5, fmhd6, fmhd7, fmhd8
      integer :: kindrn(50)
      double precision :: rsclm(50), rsclx(50), rsclz(50), rsclcm(50)
      character(len=256) :: fliv95
      character(len=256) :: fopale, fopale01, fcondopacp, fopale06
      character(len=256) :: opecalex(7)
      character(len=256) :: falex06

! --- NAMELIST-only variables that had no explicit declaration in the
!     original (relying on IMPLICIT typing); kept at original
!     spelling, typed per the implicit rule their first letter implies
!     (A-H,O-Z => double precision; L => logical) ---
      double precision :: alfa, fk, wmax_sun, pmmwmax, zalex2, zopal952
      logical :: lpmm, lsolwind

! --- true local variables (never referenced by COMMON, NAMELIST, or
!     the argument list): freely renamed for readability ---
      character(len=3) :: element_id(12)
      character(len=8) :: mixture_id_table(4)
      double precision :: zx_mix_table(4), frac_c_table(4), &
           frac_n_table(4), frac_o_table(4)
      character(len=256) :: fcalcad
      character(len=256) :: control_nml_file, physics_nml_file
      character(len=256) :: shell_cmd
      integer :: i, j, last_slash_idx
      integer :: short_prefix_len
      double precision :: one_third, two_thirds
! parmin_ln10: this file's own private ln(10) (never read after being
! set -- see the assignment below), distinct from const_lib's ln10
! (which this file now also uses via const3's `use const_lib`, hence
! the rename needed here to avoid a name collision).
      double precision :: parmin_ln10
      double precision :: sum_frac
      integer :: nkind
      integer :: first_model_binary_lu, last_model_binary_lu, &
           stored_models_binary_lu

! common /vnewcb/
      double precision :: vnew(12)
      common /vnewcb/ vnew


! common /lunum/
      integer :: ifirst, irun, istand, ifermi, iopmod, iopenv, iopatm, idyn, illdat, isnu, iscomp, &
           ikur
      common /lunum/ ifirst, irun, istand, ifermi, iopmod, iopenv, iopatm, idyn, illdat, isnu, &
           iscomp, ikur

! common /iomonte/
      character(len=256) :: fmonte1, fmonte2
      integer :: imonte1, imonte2
      common /iomonte/ fmonte1, fmonte2, imonte1, imonte2

! common /desc/
      character(len=256) :: descrip(2)
      common /desc/ descrip

! common /ccout/
      logical :: lstore, lstatm, lstenv, lstmod, lstphys, lstrot, lscrib, lstch, lphhd
      common /ccout/ lstore, lstatm, lstenv, lstmod, lstphys, lstrot, lscrib, lstch, lphhd

! common /ccout1/
      integer :: npenv, nprtmod, nprtpt, npoint
      common /ccout1/ npenv, nprtmod, nprtpt, npoint

! common /pulsegyre/: new (2026) dedicated block for the GYRE-format
! periodic pulsation-structure output feature (io/write_gyre_pulse.f90,
! triggered from io/wrtout.f90 every pulse_gyre_interval converged
! models). Deliberately a separate block from the existing common
! /pulse/ (shared by 9 files already) to keep this purely additive.
      integer :: pulse_gyre_interval
      common /pulsegyre/ pulse_gyre_interval

! common /ccout2/
      logical :: ldebug, lcorr, lmilne, ltrack, lstpch
      common /ccout2/ ldebug, lcorr, lmilne, ltrack, lstpch

! common /cenv/
      double precision :: tridt, tridl, senv0
      logical :: lsenv0, lnew0
      common /cenv/ tridt, tridl, senv0, lsenv0, lnew0

! common /ckind/
      double precision :: rescal(4,50)
      integer :: nmodls(50), iresca(50), numrun
      logical :: lfirst(50)
      common /ckind/ rescal, nmodls, iresca, lfirst, numrun

! common /comp/
      double precision :: xenv, zenv, zenvm, amuenv, fxenv(12), xnew, znew, stotal, senv
      common /comp/ xenv, zenv, zenvm, amuenv, fxenv, xnew, znew, stotal, senv

! clsun/crsun: NAMELIST /physics/ members (must keep this exact
! spelling); copied into const_lib's solar_luminosity_cgs/
! solar_radius_cgs after the namelist read below. The other six former
! common/const/ members (clsunl/clnsun/cmsun/cmsunl/crsunl/cmbol) are
! unused in this file -- setup/setups.f90 computes them from
! solar_luminosity_cgs/solar_radius_cgs at startup -- so they're
! dropped entirely rather than carried forward.
      double precision :: clsun, crsun



! common /ct2/
      double precision :: dtwind
      common /ct2/ dtwind

! common /ct3/
      logical :: lptime
      common /ct3/ lptime

! common /ctol/
      double precision :: htoler(5,2), fcorr0, fcorri, fcorr, hpttol(12)
      integer :: niter1, niter2, niter3
      common /ctol/ htoler, fcorr0, fcorri, fcorr, hpttol, niter1, niter2, niter3

! common /difus/
      double precision :: dtdif, djok
      integer :: itdif1, itdif2
      common /difus/ dtdif, djok, itdif1, itdif2

! lovste: NAMELIST /physics/ member (must keep this exact spelling).
! Former common /dpmix/; every other member already matches
! const_lib's canonical spelling, so only lovste needs to stay local
! (copied into const_lib's envelope_overshoot_active right after the
! namelist read below).
      logical :: lovste

! common /envgen/
      double precision :: atmstp, envstp
      logical :: lenvg
      common /envgen/ atmstp, envstp, lenvg

! lexcom: NAMELIST /physics/ member (must keep this exact spelling);
! copied into const_lib's use_extended_composition after the namelist
! read below.
      logical :: lexcom

! common /heflsh/
      logical :: lkuthe
      common /heflsh/ lkuthe

! common /intatm/
      double precision :: atmerr, atmd0, atmbeg, atmmin, atmmax
      common /intatm/ atmerr, atmd0, atmbeg, atmmin, atmmax

! common /intenv/
      double precision :: enverr, envbeg, envmin, envmax
      common /intenv/ enverr, envbeg, envmin, envmax

! stolr0/imax/nuse: NAMELIST /physics/ members (must keep this exact
! spelling, see this file's naming note at the top). Former common
! /intpar/; copied into intpar_lib's canonically-named variables after
! the namelist read below, since intpar_lib's names differ from these
! namelist-fixed ones.
      double precision :: stolr0
      integer :: imax, nuse

! tscut: NAMELIST /physics/ member (must keep this exact spelling).
! Former common/ctlim/ member alongside atime/tcut/tenv0/tenv1/tgcut;
! those five kept their spelling when ctlim moved to const_lib, but
! tscut did not (const_lib calls it saha_log10t_cutoff), so it stays
! local here and is copied into that canonical name after the
! namelist read below, same treatment as stolr0/imax/nuse above.
      double precision :: tscut

! common /label/
      double precision :: xenv0, zenv0
      common /label/ xenv0, zenv0

! common /newcmp/
      double precision :: xnewcp
      integer :: inewcp
      logical :: lnewcp, lrel
      common /newcmp/ xnewcp, inewcp, lnewcp, lrel

! common /newmx/
      integer :: isetmix, isetiso
      logical :: lmixture, lisotope
      double precision :: frac_c, frac_n, frac_o, r12_13, r14_15, r16_17, r16_18, zxmix, xh2_ini, &
           xhe3_ini, xli6_ini, xli7_ini, xbe9_ini, xb10_ini, xb11_ini
      common /newmx/ isetmix, isetiso, lmixture, lisotope, frac_c, frac_n, frac_o, r12_13, r14_15, &
           r16_17, r16_18, zxmix, xh2_ini, xhe3_ini, xli6_ini, xli7_ini, xbe9_ini, xb10_ini, &
           xb11_ini

! common /optab/
      double precision :: optol, zsi
      integer :: idt, idd(4)
      common /optab/ optol, zsi, idt, idd

! lrot/linstb: NAMELIST /physics/ members (must keep this exact
! spelling). Former common /rot/; every other member already matches
! const_lib's canonical spelling, so only these two need to stay local
! (copied into const_lib's rotation_active/instability_transport_active
! right after the namelist read below).
      logical :: lrot, linstb

! common /sett/
      double precision :: endage(50), setdt(50), end_dcen(50), end_xcen(50), end_ycen(50)
      logical :: lendag(50), lsetdt(50)
      common /sett/ endage, setdt, lendag, lsetdt, end_dcen, end_xcen, end_ycen

! common /vmult/
      double precision :: fw, fc, fo, fes, fgsf, fmu, fss, rcrit
      common /vmult/ fw, fc, fo, fes, fgsf, fmu, fss, rcrit

! common /debhu/
      double precision :: cdh, etadh0, etadh1, zdh(18), xxdh, yydh, zzdh, dhnue(18)
      logical :: ldh
      common /debhu/ cdh, etadh0, etadh1, zdh, xxdh, yydh, zzdh, dhnue, ldh

! common /vmult2/
      double precision :: fesc, fssc, fgsfc
      integer :: ies, igsf, imu
      common /vmult2/ fesc, fssc, fgsfc, ies, igsf, imu

! common /gravst/
      double precision :: grtol
      integer :: ilambda, niter_gs
      logical :: ldify
      common /gravst/ grtol, ilambda, niter_gs, ldify


! common /burtol/
      double precision :: cmin, abstol, reltol
      integer :: kemmax
      common /burtol/ cmin, abstol, reltol, kemmax

! common /lopal95/
      integer :: iliv95
      common /lopal95/ iliv95

! common /gravs2/
      double precision :: dt_gs, xmin, ymin
      logical :: lthoulfit
      common /gravs2/ dt_gs, xmin, ymin, lthoulfit

! common /gravs3/
      double precision :: fgry, fgrz
      logical :: lthoul, ldifz
      common /gravs3/ fgry, fgrz, lthoul, ldifz

! common /gravs4/
      logical :: lnewdif, ldifli
      common /gravs4/ lnewdif, ldifli

! common /pulse/
      double precision :: xmsol
      logical :: lpulse
      integer :: ipver
      common /pulse/ xmsol, lpulse, ipver

! common /po/
      double precision :: poa, pob, poc, pomax
      logical :: lpout
      common /po/ poa, pob, poc, pomax, lpout

! common /track/
      integer :: itrver
      common /track/ itrver

! common /atmos/
      double precision :: hras
      integer :: kttau, kttau0
      logical :: lttau
      common /atmos/ hras, kttau, kttau0, lttau

! common /mhd/
      logical :: lmhd
      integer :: iomhd1, iomhd2, iomhd3, iomhd4, iomhd5, iomhd6, iomhd7, iomhd8
      common /mhd/ lmhd, iomhd1, iomhd2, iomhd3, iomhd4, iomhd5, iomhd6, iomhd7, iomhd8

! common /core/
      logical :: lcore
      integer :: mcore
      double precision :: fcore
      common /core/ lcore, mcore, fcore

! common /nwlaol/
      double precision :: olaol(12,104,52), oxa(12), ot(52), orho(104), tollaol
      integer :: iolaol, numofxyz, numrho, numt, iopurez
      logical :: llaol, lpurez
      common /nwlaol/ olaol, oxa, ot, orho, tollaol, iolaol, numofxyz, numrho, numt, llaol, lpurez, &
           iopurez

! common /chrone/
      logical :: lrwsh, liso
      integer :: iiso
      common /chrone/ lrwsh, liso, iiso

! common /newxym/
      double precision :: xenv0a(50), zenv0a(50), cmixla(50), senv0a(50)
      logical :: lsenv0a(50)
      common /newxym/ xenv0a, zenv0a, cmixla, lsenv0a, senv0a

! common /atmos2/
      double precision :: atmpl(n_atm_teff,n_atm_logg), atmtl(n_atm_teff), atmgl(n_atm_logg), atmz
      integer :: ioatm
      common /atmos2/ atmpl, atmtl, atmgl, atmz, ioatm

! common /atmos2c/
      double precision :: atmplc(n_katm_teff,n_katm_logg), atmtlc(n_katm_teff), atmglc(n_katm_logg)
      common /atmos2c/ atmplc, atmtlc, atmglc

! lnulos1: NAMELIST /physics/ member (must keep this exact spelling);
! copied into const_lib's use_itoh_neutrino_loss after the namelist
! read below. dsnudt/dsnudd (former common/nuloss/'s other two
! members) are unused in this file and confirmed dead everywhere else
! too (nuclear/engeb.f90's neutrino_dlnq_dlnt/neutrino_dlnq_dlnd,
! which shared this block only by position, are genuinely local to
! that file and were made plain locals there), so they're dropped
! entirely rather than carried forward.
      logical :: lnulos1

! common /cals2/
      double precision :: toll, tolr, tolz, calsolzx, calsolage
      logical :: lcals, lcalsolzx
      common /cals2/ toll, tolr, tolz, lcals, lcalsolzx, calsolzx, calsolage

! common /zramp/
      double precision :: rsclzc(50), rsclzm1(50), rsclzm2(50)
      integer :: iolaol2, ioopal2, nk
      logical :: lzramp
      common /zramp/ rsclzc, rsclzm1, rsclzm2, iolaol2, ioopal2, nk, lzramp

! common /calstar/
      double precision :: xls, xlstol, steff, sr, bli, alri, ager, blr, blrp, agei
      logical :: lstar, lteff, lpassr, lcalst
      common /calstar/ xls, xlstol, steff, sr, bli, alri, ager, blr, blrp, agei, lstar, lteff, &
           lpassr, lcalst

! common /opaleos/
      logical :: lopale, lopale01, lopale06, lnumderiv
      integer :: iopale
      common /opaleos/ lopale, iopale, lopale01, lopale06, lnumderiv

! common /newopac/
      double precision :: zlaol1, zlaol2, zopal1, zopal2, zopal951, zalex1, zkur1, zkur2, tmolmin, &
           tmolmax
      logical :: lalex06, llaol89, lopal92, lopal95, lkur90, lalex95, l2z
      common /newopac/ zlaol1, zlaol2, zopal1, zopal2, zopal951, zalex1, zkur1, zkur2, tmolmin, &
           tmolmax, lalex06, llaol89, lopal92, lopal95, lkur90, lalex95, l2z

! common /miscopac/
      integer :: ikur2, icondopacp
      logical :: lcondopacp
      common /miscopac/ ikur2, icondopacp, lcondopacp

! common /alexo/
      integer :: ialxo
      common /alexo/ ialxo

! common /alex06/
      integer :: ialex06
      common /alex06/ ialex06

! common /alexmix/
      double precision :: xalex, zalex
      common /alexmix/ xalex, zalex

! common /varfc/
      double precision :: vfc(max_diag_pts)
      logical :: lvfc, ldifad
      common /varfc/ vfc, lvfc, ldifad

! common /notran/
      logical :: lnoj
      common /notran/ lnoj

! sstandard/lnewnuc: NAMELIST /physics/ members (must keep this exact
! spelling, see this file's naming note at the top). Former common
! /cross/; only these two are actually read anywhere (sstandard's
! namelist value is never referenced outside its own declaration --
! setup/remap.f90 fully recomputes const_lib's cross_section_scale
! from other inputs regardless -- so it's dropped here rather than
! copied). lnewnuc is copied into const_lib's use_new_nuclear_rates
! right after the namelist read below, since remap.f90 needs it.
      double precision :: sstandard(17)
      logical :: lnewnuc

! common /newcross/
      double precision :: s0_1_1, s0_3_3, s0_3_4, s0_1_12, s0_1_13, s0_1_14, s0_1_16, s0_pep, &
           s0_1_be7e, s0_1_be7p, s0_hep, s0_1_15_c12alp, s0_1_15_o16, s0p_1_1, s0p_3_3, s0p_3_4, &
           s0p_1_12, s0p_1_13, s0p_1_14, s0p_1_16, s0pp_1_12, s0pp_1_13, s0pp_1_16, s0p_1_be7p, &
           s0pp_1_be7p
      common /newcross/ s0_1_1, s0_3_3, s0_3_4, s0_1_12, s0_1_13, s0_1_14, s0_1_16, s0_pep, &
           s0_1_be7e, s0_1_be7p, s0_hep, s0_1_15_c12alp, s0_1_15_o16, s0p_1_1, s0p_3_3, s0p_3_4, &
           s0p_1_12, s0p_1_13, s0p_1_14, s0p_1_16, s0pp_1_12, s0pp_1_13, s0pp_1_16, s0p_1_be7p, &
           s0pp_1_be7p

! common /newparam/
      double precision :: flag_dx, flag_dw, flag_dz, time_core_min, time_dl, time_dp, time_dr, &
           time_dt, time_dw_global, time_dw_mix, time_dx_core_frac, time_dx_core_tot, time_dx_shell, &
           time_dx_total, time_dy_core_frac, time_dy_core_tot, time_dy_shell, time_dy_total, &
           tol_czbase_fine_width, tol_dl_max, tol_dm_max, tol_dm_min, tol_dp_core_max, &
           tol_dp_czbase_max, tol_dp_env_max, tol_dx_max, tol_dz_max, time_max_dt_frac
      logical :: lstruct_time, lnewvars
      common /newparam/ flag_dx, flag_dw, flag_dz, time_core_min, time_dl, time_dp, time_dr, time_dt, &
           time_dw_global, time_dw_mix, time_dx_core_frac, time_dx_core_tot, time_dx_shell, &
           time_dx_total, time_dy_core_frac, time_dy_core_tot, time_dy_shell, time_dy_total, &
           tol_czbase_fine_width, tol_dl_max, tol_dm_max, tol_dm_min, tol_dp_core_max, &
           tol_dp_czbase_max, tol_dp_env_max, tol_dx_max, tol_dz_max, time_max_dt_frac, lstruct_time, &
           lnewvars

! common /monte/
      logical :: lmonte
      integer :: imbeg, imend
      common /monte/ lmonte, imbeg, imend

! common /scveos/
      double precision :: tlogx(n_scv_teff), tablex(n_scv_teff,n_scv_press,12), &
           tabley(n_scv_teff,n_scv_press,12), smix(n_scv_teff,n_scv_press), &
           tablez(n_scv_teff,n_scv_press,13), tablenv(n_scv_teff,n_scv_press,12)
      integer :: nptsx(n_scv_teff), idtt, idp
      logical :: lscv
      common /scveos/ tlogx, tablex, tabley, smix, tablez, tablenv, nptsx, lscv, idtt, idp

! common /scv2/
      integer :: iscvh, iscvhe, iscvz
      common /scv2/ iscvh, iscvhe, iscvz

! common /alatm03/
      double precision :: alatm_feh, alatm_alpha
      logical :: laltptau100
      integer :: ioatma
      common /alatm03/ alatm_feh, alatm_alpha, laltptau100, ioatma

! common /alatm04/
      double precision :: dummy1, dummy2, dummy3, dummy4
      common /alatm04/ dummy1, dummy2, dummy3, dummy4

! common /disk/
      double precision :: sage, tdisk, pdisk
      logical :: ldisk
      common /disk/ sage, tdisk, pdisk, ldisk

! weakscreening: NAMELIST /physics/ member (must keep this exact
! spelling). Former common /weak/; copied into const_lib's
! weak_screening_threshold right after the namelist read below.
      double precision :: weakscreening

! common /sbrot/
      logical :: lsolid
      integer :: impjmod
      common /sbrot/ lsolid, impjmod

! dmdt0/compacc/lmdot: NAMELIST /physics/ members (must keep this
! exact spelling). Former common /masschg/; fczdmdt/ftotdmdt/creim/
! lreimer already match const_lib's canonical spelling, so only these
! three need to stay local (copied into const_lib's
! mass_accretion_rate/accreted_composition/use_mass_accretion right
! after the namelist read below).
      double precision :: dmdt0, compacc(15)
      logical :: lmdot

! common /cmixing/
      double precision :: cstmixing, cstdiffmix
      common /cmixing/ cstmixing, cstdiffmix

! common /acdpth/
      double precision :: tauczn, deladj(max_diag_pts), tauhe, tnorm, tcz, whe, acatmr(max_diag_pts), &
           acatmd(max_diag_pts), acatmp(max_diag_pts), acatmt(max_diag_pts), tatmos, ageout(5)
      logical :: lclcd, ljlast, ljwrt, ladon, laoly, lacout
      integer :: iclcd, iacat, ijlast, ijvs, ijent, ijdel
      common /acdpth/ tauczn, deladj, tauhe, tnorm, tcz, whe, acatmr, acatmd, acatmp, acatmt, tatmos, &
           ageout, lclcd, iclcd, iacat, ijlast, ljlast, ljwrt, ladon, laoly, ijvs, ijent, ijdel, &
           lacout

! common /govs/
      logical :: ltrist
      common /govs/ ltrist

! common /pmmwind/
      double precision :: pmma, pmmb, pmmc, pmmd, pmmm, pmmjd, pmmmd, pmmsolp, pmmsolw, pmmsoltau
      logical :: lmwind, lrossby, lbscale
      character(len=3) :: awind
      common /pmmwind/ pmma, pmmb, pmmc, pmmd, pmmm, pmmjd, pmmmd, pmmsolp, pmmsolw, pmmsoltau, &
           lmwind, lrossby, lbscale, awind

! common /cwind/
      double precision :: wmax, exmd, exw, extau, exr, exm, exl, expr, constfactor, structfactor, &
           excen, c_2
      logical :: ljdot0
      common /cwind/ wmax, exmd, exw, extau, exr, exm, exl, expr, constfactor, structfactor, excen, &
           c_2, ljdot0

! lnewtcz/lcalcenv: NAMELIST /physics/ members (must keep this exact
! spelling); copied into const_lib's use_new_turnover_timescale/
! calc_envelope_flag after the namelist read below. taucz/taucz0/
! pphot/pphot0/fracstep (former common/ovrtrn/'s other five members)
! are unused in this file -- they're genuinely evolving per-model
! state read/written by many distant files, now state/turnover_lib.f90
! -- so they're dropped from this file's own declarations entirely.
      logical :: lnewtcz, lcalcenv

! common /mag/
      double precision :: codm
      logical :: lcodm
      common /mag/ codm, lcodm

! common /xsect/
      double precision :: xsli6, xsli7, xsbe91, xsbe92, xsbe93
      logical :: lxli6, lxli7, lxbe91, lxbe92, lxbe93
      common /xsect/ xsli6, xsli7, xsbe91, xsbe92, xsbe93, lxli6, lxli7, lxbe91, lxbe92, lxbe93

! sli6/sli7/sbe91/sbe92/sbe93 are themselves NAMELIST /physics/ members
! (see the "G Somers 6/14" list below) so must keep their exact
! spelling (this file's naming note at the top); the canonical
! const_lib names (li6_rate_scale etc, former common/burnscs/) are set
! via copy-assignment once these are computed below.
      double precision :: sli6, sli7, sbe91, sbe92, sbe93

! common /spots/
      double precision :: spotf, spotx
      logical :: lsdepth
      common /spots/ spotf, spotx, lsdepth

! common /version/
      character(len=10) :: yrecver
      character(len=20) :: githash
      common /version/ yrecver, githash

      save
!
!
! SPLIT NAMELIST INTO TWO: CONTROL and PHYSICS
      namelist /control/ &
           &    cmixla, calsolage, calsolzx, &
           &    descrip, &
           &    endage, &
           &    flaol, fpurez,flaol2, fopal2, &
           &    flast, ffirst, ffermi, fdebug, ftrack, fshort, fstch, &
           &    fmilne, fmodpt, fstor, fpmod, fpatm, fpenv, &
           &    fdyn, flldat, fsnu, fscomp, fkur, fmhd1, &
           &    fmhd2, fmhd3, fmhd4, fmhd5, fmhd6, fmhd7, fmhd8, fiso, fatm, &
           &    fkur2, fallard, fscvh, fscvhe, fscvz, fopale, fliv95, &
           &    fmonte1,fmonte2, &
           &    ipver, itrver, &
           &    kindrn, &
           &    ldebug, lcorr, lmilne, ltrack, lstore, lfirst, &
           &    lstpch, lscrib, lstch, &
! G Somers 11/14
           &    lstatm,lstenv,lstmod,lstphys,lstrot, &
! G Somers END
           &    lpulse, lzramp, lteff, lcalst, lpurez, &
! MHP 9/24 add LCALSOLZX to namelist
           &    liso, lrwsh, lsenv0a, lpout,lcals,lcalsolzx, &
           &    llaol89,lopal92,lopal95,lkur90,lalex95, &
           &    npoint, &
           &    npenv, nprtmod, nprtpt, numrun, nmodls, &
           &    opecalex, &
           &    poa, pob, poc, pomax, &
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
           &    cmin, clsun, crsun, &
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
! G Somers 3/17 USE NEW OVERTURN TIMESCALE CALC?
           &    lnewtcz, lcalcenv
!
! DBG DATA CARDS FOR THE RUN PARAMETERS
! MHP DATA FOR MONTE CARLO OPTION, ETC
      data lmonte,imbeg,imend/.false.,1,1/
! Changed slightly 3He-3He on 9/25/97 to take account of the S'.
!  Previously (6/16/97) used S at Gamow Peak. Agrees with Workshop paper.
!
      data weakscreening/0.03d0/
      data sstandard/0.9828,1.0485,0.9815,0.9241,1.3818,1.0542,1.0, &
           &  1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0108,0.7819,0.2875/
! MHP 7/93 VARIABLE FC OPTION
! MHP 9/94 COMBINED DIFFUSION/ADVECTION OPTION
      data lvfc,ldifad/.false.,.false./
! MHP 9/93
      data lnoj/.false./
      data tdisk,pdisk,ldisk/0.0d0,7.2722d-6,.false./
      data xalex/0.7e0/
      data zalex/0.02e0/
      data lsenv0a, senv0a /50*.false.,50*1.26d-4/
      data xenv0, zenv0/0.7,0.02/
      data ldebug, lcorr, npoint, lmilne, ltrack, lstore, &
           &  lstpch &
           &   /.false., .true., 1, .false., .true., .false.,.false./
      data lscrib/.true./
      data lstch/.false./
      data lenvg, atmstp, envstp/.false.,0.5,0.5/
      data nprtmod, nprtpt/1,1/
      data pulse_gyre_interval/0/
      data numrun, kindrn, lfirst, nmodls &
           &      /1,50*1,50*.true.,50*0/
! MHP 10/24 ADDED NEW DEFAULTS FOR END CONDITIONS ON CENTRAL D,X,Y
      data endage, setdt, rsclm, rsclx, rsclz, rsclcm &
           &      /50*0.0,50*0.0,50*0.0,50*0.0,50*0.0,50*0.0/
      data end_dcen,end_xcen, end_ycen &
           &      /50*0.0,50*0.0,50*0.0/
      data element_id/'HE3','C12','C13','N14','N15','O16','O17','O18','H2 ', &
           & 'LI6','LI7','BE9'/
      data optol,zsi/1.0d-8,0.0d0/
! tcut/saha_log10t_cutoff/tenv0/tenv1/tgcut defaults moved to
! const_lib.f90 (see its own header note): DATA can no longer target
! them here now that they're use-associated from const_lib.
      data atmerr,atmd0,atmbeg,atmmin,atmmax/3.0d-4,1.0d-10,1.0d-1, &
           &      1.0d-1,5.0d-1/
      data enverr,envbeg,envmin,envmax/3.0d-4,1.0d-1,1.0d-1,5.0d-1/
      data stolr0,imax,nuse/1.0d-3,11,7/
      data dtdif,djok,itdif1,itdif2/1.0d-2,1.0d-4,1,1/
      data htoler/6.0d-5,4.5d-5,3.0d-5,9.0d-5,3.0d-5,9.0d-1,5.0d-1, &
           &      5.0d-1,2.0d0,2.5d-6/
      data hpttol/1.0d-8,8.0d-2,5.0d-2,5.0d-2,1.0d0,1.0d0,0.0d0,5.0d-2, &
           &             2.0d-2,5.0d-2,5.0d-2,1.0d-1/
      data lnewcp,anewcp,xnewcp/.false.,'   ',1.3d1/
! MHP 10/24 ADDED NEW MIXTURE CONTROL ISETISO CONTROLS CNO ISOTOPE RATIOS AND
! LIGHT ELEMENT ABUNDANCES D,HE3,LI6,LI7,BE9,B10,B11 (1=USED)
! ISETMIX CONTROLS C+N+O MASS FRACTIONS (1=USED)
! AMIX AND AISO ARE STRINGS IDENTIFYING EITHER A PRESET MIXTURE OR A CUSTOM ONE ('CUS')
! SUPPORTED AMIX AT PRESENT ARE 'GS98','AAG21',M22M','M22P'. SUPPORTED AISO IS 'L21'.
! FOR A CUSTOM MIXTURE YOU CAN ENTER INDIVIDUAL VALUES.
! TO BE ADDED - AMIX FROM ATOMIC OPACITY TABLES (INEWMIX=2) AND OTHER MIXTURES/ISOTOPES
      data isetiso,isetmix,amix,aiso/0,0,'GS98','L21'/
!     L21 DEFAULT ISOTOPE DATA LODDERS ET AL. 2021 SSRV 217,44
      data r12_13,r14_15,r16_17,r16_18,xh2_ini,xhe3_ini,xli6_ini, &
           &      xli7_ini,xbe9_ini,xb10_ini,xb11_ini/88.27d0,411.9d0,471.4d0, &
           &      2693.0d0,2.781d-5,3.461d-5,7.187d-10,1.025d-8,1.595d-10, &
           &      1.002d-9,4.405d-9/
!     GS98 DEFAULT CNO FRACTIONS OF METALS GREVESSE&sAUVAL 1998 SSRV 85,161
      data zxmix,frac_c,frac_n,frac_o/0.02292d0,0.172148d0,0.050426d0, &
           &      0.468195/
! MHP 10/24 DATA FOR CNO FRACTIONS AND Z/X OF DIFFERNT SOLAR MIXTURES.
!     AMIXT IS THE LIST OF IDS,EACH OF WHICH HAS A ZX AND CNO FRACS
!     ENTRY 1 =GS98(IN PARMIN),2=ASPLUND ET AL. 2021 A&A 653,141
!     3,4=MAGG ET AL. 2021 (MET,PHOT) A&A 661,140
      data mixture_id_table,zx_mix_table,frac_c_table,frac_n_table,frac_o_table/'GS98','AAG21', &
           &      'M22P','M22M',0.02292d0,0.0187d0,0.0225d0,0.0226d0, &
           &  0.172148d0,0.184156d0,0.19239d0,0.191425d0, &
           &  0.050426d0,0.050344d0,0.059012d0,0.058716d0, &
           &  0.468195d0,0.416592d0,0.415631d0,0.413545d0/
! acfpft/itfp1/itfp2 defaults moved to const_lib.f90 (former
! common/rot/).
      data tridt,tridl/1.0d-2,8.0d-2/
      data niter1,niter2,niter3,fcorr0,fcorri/2,20,2,0.8d0,0.1d0/
! atime's default moved to const_lib.f90 -- see the tcut/etc. note
! above; ATIME(13) was orginally = 1.5.
      data dtwind /1.0d1/
      data lptime/.true./
! JVS 04/14
      data ltrist/.false./
! cmixl's default (1.4d0) moved to const_lib.f90's own declaration: a
! DATA statement can no longer target it here now that it's
! use-associated from const_lib rather than locally declared/COMMON.
      data lkuthe/.false./
!       DATA DPENV,LNSTDMX,LOVSTC,ALPHAC,LOVSTE,ALPHAE, LOVSTM, ALPHAM
!      */1.0D0,.FALSE.,.FALSE., 0.0D0, .FALSE.,0.0D0, .FALSE., 0.0/
! dpenv/lovstc/alphac/alphae/lovstm/alpham defaults moved to
! const_lib.f90 (former common/dpmix/); lovste stays local (NAMELIST
! spelling).
      data lovste/.false./
! ladov/lovmax/betac defaults moved to const_lib.f90 (former
! common/dpmix/).
! JVS 07/13
! END JVS
      data lnew0,lexcom/.false.,.false./
! walpcz/lwnew/wnew defaults moved to const_lib.f90 (former
! common/rot/); lrot/linstb stay local (NAMELIST spelling).
      data lrot,linstb/.false.,.false./
      data ljdot0,alfa,fk/.true.,1.5d0,1.0d0/
      data fw,fc,fo,fes,fgsf,fmu,fss,rcrit/1.0d0,1.0d0, &
           &      1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d3/
! MHP 8/17 INITIALIZED WMAX_SUN
      data wmax,wmax_sun/3.0d-4,1000.0/
! DBG PULSE DATA CARD FOR PULSATION
      data lpulse/.false./
      data ipver/1/
      data itrver/1/
      data kttau/0/
      data clsun,crsun/3.8515d33,6.9598d10/
! YC  If LMHD is TRUE use MHD equation of state tables.  LU numbers
!     are stored in IOMHDi.
! DBG If LCORE is TRUE then calculate shells interior to start up
!     model's inner most shell.
      data lcore,mcore,fcore/.false.,0,1.0/
      data lmhd/.false./
! MHP 5/90 NEW DATA STATEMENTS FOR NEW PARAMETERS
      data grtol,ilambda,niter_gs,ldify/1.0d-8,1,10, &
           &      .false./
! niter4/lnews/lsnu defaults moved to const_lib.f90 (former
! common/neweng/).
      data dt_gs,xmin,ymin,lthoulfit/0.1d0,1.0d-3,1.0d-3,.false./
      data lthoul,ldifz,fgrz,fgry/.false.,.false.,1.0d0,1.0d0/
      data lnewdif,ldifli/.false.,.false./
      data cmin,abstol,reltol,kemmax/1.0d-20,1.0d-5,1.0d-4,50/
      data etadh0, etadh1, ldh/-1.0d0, 1.0d0, .false./
      data fesc,fssc,fgsfc/1.0d0,1.0d0,1.0d0/
      data ies,igsf,imu/1,1,1/
! lsemic's default moved to const_lib.f90 (former common/dpmix/).
! DBGLAOL
      data tollaol,llaol,lpurez/10.0,.false.,.false./
! DBG 11/11/91
      data lrwsh,liso/.false.,.false./
! 3/92 DBG
      data lnulos1/.false./
! DBG PULSE OUT 7/92
      data pomax,poa,pob,poc,lpout/0.1d0,1.0d0,10.0d0,0.0d0,.false./
! MHP 06/13 ADDED FLAG TO CALIBRATE TO SOLAR Z/X, SOLAR Z/X, SOLAR AGE
      data toll,tolr,tolz,lcals,lcalsolzx,calsolage,calsolzx/1.0d-5, &
           &      1.0d-4,1.0d-3,.false.,.false.,4.57d9,0.02292d0/
!      DATA TOLL,TOLR,LCALS/1.0D-5,1.0D-4,.FALSE./
! DBG 4/94 ZRAMP STUFF
      data rsclzc, rsclzm1, rsclzm2, lzramp/50*-1.0d0, 50*-1.0d0, &
           &        50*-1.0d0, .false./
! DBG 12/94 CALIBRATED STELLAR MODEL STUFF
      data lcalst, lteff/.false., .false./
! YCK >>>  2/95 OPAL eos
! LLP >>> OPAL 2001 EOS, Potekhin Conductive Opacities,
!         OPAL 2006 EOS, Use Numerical Derivitives switches
      data lopale, lopale01,lcondopacp,lopale06,lnumderiv &
           &      /.false.,.false.,.false.,.false.,.false./
! MHP 8/25 Removed hard coded defaults
!     Alex low T opacity
!      DATA OPECALEX/'OPACALEXANDER.X00',
!     +              'OPACALEXANDER.X01',
!     +              'OPACALEXANDER.X02',
!     +              'OPACALEXANDER.X035',
!     +              'OPACALEXANDER.X05',
!     +              'OPACALEXANDER.X07',
!     +              'OPACALEXANDER.X08'/

! DBG 1/96 THE ARRAY V, READ IN VIA RDLAOL, CONTAINED THE MASS FRACTIONS
! OF THE ENVELOPE ELEMENTS. IT WAS USED IN STARIN TO DEFINE FXENV,
! WHICH ARE THE NUMBER FRACTION OF THE ENVELOPE ELEMENTS. FXENV WAS
! THEN UPDATED IN EQSTAT AND EQSAHA. HERE WE DEFINE VNEW PASSED
! IN A COMMON BLOCK VNEWCB. IT IS IDENTICAL TO V EXCEPT THAT THE NUMBERS
! ARE DEFINED HERE EXPLICITY FOR A G&N93 SOLAR MIXTURE. YOU CAN
! CHANGE THEM VIA THE PHYSICS NAMELIST. V IS SET EQUAL TO VNEW IN
!  STARIN EXCEPT WHEN LLAOL=T (TO MAINTAIN BACKWARD COMPATIBILITY.
!            Na          Al          Mg          Fe
!            Si          C           H           O
!            N           Ar          Ne          He
      data vnew/0.001999d0, 0.003238d0, 0.037573d0, 0.071794d0, &
           &           0.040520d0, 0.173285d0, 0.000000d0, 0.482273d0, &
           &           0.053152d0, 0.005379d0, 0.098668d0, 0.000000d0/
! MHP 5/97 OPTION FOR SAUMON, CHABRIER, AND VAN HORN EOS ADDED
      data lscv/.false./
! MHP 3/99 OPTION FOR SB ROTATION ENFORCED IN THE ENTIRE STAR AT
! ALL TIMES
! JNT 09/2025 FOR 05/15 IMPJMOD default set to 0
      data lsolid, impjmod/.false.,0/   !JNT 09/2025
! fczdmdt/ftotdmdt/creim/lreimer defaults moved to const_lib.f90
! (former common/masschg/); dmdt0/compacc/lmdot stay local (NAMELIST
! spelling).
      data dmdt0,compacc,lmdot &
           &      /-1.0d-14, &
           &      0.71668d0,0.265721d0,0.01757d0,2.9d-5, &
           &      3.013d-3,3.385d-5,9.346d-4,0.0d0,8.462d-3,0.0d0, &
           &      1.696d-5,0.0d0,2.0d-9,2.0d-9,3.0d-11, &
           &      .false./
! mhp 8/10 added scaled solar wind mass loss option
!      DATA LSOLWIND,DMSUN,DMWSUN,DMWMAX/.FALSE.,-2.0D-14,2.863E-6,9.054E-5/
      data lsolwind/.false./
! 3/09 Alexander 2006 opacity table options and opacity ramp options
      data tmolmin,tmolmax,lalex06/4.0d0,4.1d0,.false./
!FD 10/09 Mimic mixing options - acting on setling and differential settling
      data cstmixing, cstdiffmix/1.0,1.0/
! JVS 02/11 Initialize acoustic depth common block values appropriately
      data ageout/0.5d0, 1.0d0, 5.0d0, 10.0d0, 20.0d0/
      data lclcd, ljlast, ljwrt, lacout/.false.,.false., .false.,.false./
! JVS end
! MHP 02/12 NEW PARAMETERIZATION OF ANGULAR MOMENTUM AND MASS LOSS
! FROM MAGNETIZED SOLAR-LIKE WINDS
!      DATA LPMM,PMMA,PMMB,PMMC,PMMM,PMMJD,PMMMD,PMMWMAX,
!     *     PMMSOLP,PMMSOLW,PMMSOLTAU
!     *   /.FALSE.,2.0D0,1.0D0,0.0D0,0.22D0,1.32E30,1.27E12,
!     *     2.836E-5,4.9304D0,2.836E-6,1.065E6/
! G Somers 6/16 NEW PARAMETERIZATION OF ANGULAR MOMENTUM AND MASS LOSS FROM
! SOLAR WINDS. FOLLOW MATT ET AL (2012) FORMULATION, BUT DEFAULT TO KAWALER
! TYPE LAW.
      data lmwind,lrossby,lbscale/.false.,.true.,.false./
! MHP 8/17 CHANGED DEFAULT FOR PMMA TO 2 FROM 0
      data awind,pmma,pmmb,pmmc,pmmd,pmmm/'K97',2.0,1.0,2.0,0.0,0.5/
      data pmmjd,pmmmd,pmmsolp,pmmsolw,pmmsoltau/1.32e30,1.27e12, &
           &      4.9304d0,2.836e-6,9.4805d5/
! MHP 02/12 PERMIT CONSTANT DIFFUSION COEFFICIENT
      data lcodm,codm/.false.,2.5d4/
!
! G Somers 06/14 ALLOW NEW LI DESTRUCTION CROSS SECTIONS
!           NEW VALUES SHOULD BE IN UNITS OF keV b.
!           DEFAULT LI6 = 5,500 keV b FROM FOWLER ET AL. 1967
!           DEFAULT LI7 = 52 keV b FROM ROLFS & KAVANAGH 1986
!           DEFAULT BE9(P,G)B10 (BE91) = 1.1 keV b INFERRED FROM FOWLER 1988
!               THEY DON'T SHOW THIS REACTION, SO THIS IS APPROXIMATE.
!           DEFAULT BE9(P,D)2HE4 (BE92) = 15,000 keV b FROM FOWLER ET AL. 1967
!           DEFAULT BE9(P,A)LI6 (BE93) = 15,000 keV b FROM FOWLER ET AL. 1967
      data xsli6, xsli7, xsbe91, xsbe92, xsbe93 &
           &      /5.5d3, 5.2d1, 1.1d0, 1.5d4, 1.5d4/
      data lxli6, lxli7, lxbe91, lxbe92, lxbe93, &
           &      sli6, sli7, sbe91, sbe92, sbe93 &
           &      /.false.,.false.,.false.,.false.,.false., &
           &      1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/
      data spotf, spotx, lsdepth/0.00, 1.00, .false./
! G Somers END
! MHP 8/14 DEFAULT CROSS-SECTIONS ARE TAKEN FROM THE SOLAR FUSION II PAPER
! REFERENCE ADELBERGER ET AL. 2011. UNITS ARE KeV b
      data s0_1_1,s0_3_3,s0_3_4,s0_1_12,s0_1_13/ &
           &      4.01d-22,5.21d3,5.6d-1,1.34d0,7.6d0/
      data s0_1_14,s0_1_16,s0_pep,s0_1_be7e/ &
           &      1.66d0,1.06d1,3.5734d-6,1.7709d-10/
! NOTE: PEP IS THE PROPORTIONALITY CONSTANT RELATIVE TO PP
! NOTE: BE7+E- IS THE PROPORTIONALITY CONSTANT IN THE LINEAR TERM
! THE CODE USES T9, NOT T6, SO ANY EXPRESSION IN TERMS OF T/10^6 K
! NEEDS TO BE DIVIDED BY 1000^0.5 (FOR BOTH PEP AND BE7+E-)
      data s0_1_be7p,s0_hep,s0_1_15_c12alp,s0_1_15_o16/ &
           &      0.0208d0,8.6d-20,7.3d4,3.6d1/
! REFERENCE FIRST DERRIVATIVES OF CROSS-SECTIONS (ADELBERGER ET AL. 2011)
! UNITS ARE b
      data s0p_1_1,s0p_3_3,s0p_3_4,s0p_1_12/ &
           &      4.49d-24,-4.9d0,-3.6d-4,2.6d-3/
      data s0p_1_13,s0p_1_14,s0p_1_16,s0p_1_be7p/ &
           &      -7.83d-3,-3.3d-3,-5.4d-2,-3.12d-5/
! REFERENCE SECOND DERIVATIVES OF CROSS SECTIONS (ADELBERGER ET AL. 2011)
      data s0pp_1_12,s0pp_1_13,s0pp_1_16,s0pp_1_be7p/ &
           &      8.3d-5,7.29d-4,0.0d0,-2.288d-7/
      data lnewnuc /.false./
      data time_core_min,time_dl,time_dp,time_dr,time_dt,time_dw_global, &
           &      time_dw_mix,time_dx_core_frac,time_dx_core_tot,time_dx_shell, &
           &      time_dx_total,time_dy_core_frac,time_dy_core_tot, &
           &      time_dy_shell,time_dy_total,time_max_dt_frac,lnewvars/ &
           &      1.0d-3,2.0d-2,4.0d-2,2.0d-2,2.0d-2,8.0d-2, &
           &      8.0d-2,0.5d0,2.0d-2,0.1d0, &
           &      1.5d-3,0.5d0,2.0d-2,0.1d0,1.5d-3,1.5d0,.false./
      data flag_dx,flag_dw,flag_dz,lstruct_time,tol_czbase_fine_width, &
           &      tol_dl_max,tol_dm_max,tol_dm_min,tol_dp_core_max, &
           &      tol_dp_czbase_max,tol_dp_env_max,tol_dx_max,tol_dz_max/ &
           &      0.05d0,0.10d0,0.05d0,.false.,0.0d0, &
           &      0.02d0,0.08d0,1.0d-8,0.05d0,0.05d0,0.05d0,1.0d0,1.0d0/
! DEFAULT TO YES FOR NEW TAUCZ CALCULATION
      data lnewtcz,lcalcenv /.true.,.true./
!
! THIS SUBROUTINE READS ALL USER DEFINED QUANTITIES FROM THE
! FILES yrec8.nml1 and yrec8.nml2
! VALUES FOR LOGICAL UNITS USED IN READS AND WRITES SET IN DATA
! STATEMENT; THEY CAN BE CHANGED IN THE NAMELIST IF NEEDED.
! LOGICAL UNIT 5 = READ FROM SCREEN
! LOGICAL UNIT 6 = WRITE TO SCREEN: FOR BATCH USE STATUS FILE INSTEAD
! SPECIFY ALL LOGICAL UNIT NUMBERS HERE:
! OUTPUT: STATUS FILE
      iowr = 6
! OUTPUT: LAST MODEL (TEXT)
      ilast = 11
! INPUT: FIRST MODEL (TEXT)
      ifirst = 12
! INPUT: PHYSICS NAMELIST
      irun = 13
! INPUT: CONTROL NAMELIST
      istand = 14
! INPUT: FERMI TABLES
      ifermi = 15
! OUTPUT: RESERVED for DEBUGGING
      idebug = 18
! OUTPUT: TRACK
      itrack = 19
! OUTPUT: ALL DIAGNOSTIC INFO
      short_file_unit = 20
! OUTPUT: MILNE INVARIANT VARIABLES
      imilne = 21
! OUTPUT: SHELL BY SHELL INFO ON MODELS
      imodpt = 22
! OUTPUT: SAVED MODELS, CAN BE USED AS STARTING MODEL
      istor = 23
! OUTPUT: FOR PULSATION CODE, INTERIOR
      iopmod = 24
! OUTPUT: FOR PULSATION CODE, ENVELOPE
      iopenv = 25
! OUTPUT: FOR PULSATION CODE, ATMOSPHERE
      iopatm = 26
! OUTPUT: BINARY OUTPUT OF LAST MODEL
      last_model_binary_lu = 27
! OUTPUT: BINARY OUTPUT OF STORED MODELS
      stored_models_binary_lu = 28
! INPUT: BINARY STARTING MODEL
      first_model_binary_lu = 29
! OUTPUT: INFO RELAVENT TO DYNAMO
      idyn = 30
! YCK INPUT: OPAL92 OPACITY TABLES
      illdat = 32
! YCK INPUT: OPAL95 OPACITY TABLE
      iliv95 = 48
! OUTPUT: SNU FLUXES
      isnu = 33
! OUTPUT: EXTENDED COMPOSITION INFO
      iscomp = 34
! YCK INPUT: KURUCZ LOW T OPACITIES
      ikur = 36
! OUPUT: ISOCHRONE INFORMATION
      iiso = 37
! INPUT: KURUCZ ATMOSPHER TABLE
      ioatm = 38
! YCK INPUT: Alex LOW T OPACITIES
      ialxo = 39
! INPUT: MHD EQU. OF STATE TABLES
      iomhd1 = 40
      iomhd2 = 41
      iomhd3 = 42
      iomhd4 = 43
      iomhd5 = 44
      iomhd6 = 45
      iomhd7 = 46
      iomhd8 = 47
! INPUT: OPAL EQUATION OF STATE
      iopale = 49
! INPUT LAOL OPACITIES IN DENSE GRID FORMAT
      iolaol = 61
! INPUT: LAOL OPACITIES FOR PURE CN IN DENSE GRID FORMAT
      iopurez = 62
! DBG 4/94
!     INPUT:
! DBG 8/95 SECOND OPOACITY TABLES FOR ZRAMP AND Z DIFFUSION
      iolaol2 = 63
      ioopal2 = 64
      ikur2 = 65
! MHP 6/97 ADDED OPTION FOR ALLARD MODEL ATMOSPHERES
      ioatma = 66
! MHP 6/98 MONTE CARLO FOR SNUs
      imonte1 = 70
      imonte2 = 71
! INPUT FILES FOR THE SCV EOS
      iscvh=72
      iscvhe=73
      iscvz=74
! Input files for Potekhin conductive opacities. LLP 7/8/06
      icondopacp = 75
!      FcondOpacP = 'condall.d'
! Default FeH and Alpha for new Allard Atmospheres
      alatm_feh = 0d0
      alatm_alpha = 0d0
      laltptau100 = .false.
! 3/09 Input file for 2006 Alexander opacities
      ialex06 = 90

      print *,''
      print *,'Yale Rotating Evolution Code - YREC, v',yrecver(1:len_trim(yrecver)),' (',githash(1:len_trim(githash)),')'

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
      write(version_fmt, 315) len_trim(yrecver), len_trim(githash)
      315 format('(''# YREC v'', A', i2.2, ', '' ('', A', i2.2, &
           &        ', '')'')')

      call getarg(1, control_nml_file)
      if (control_nml_file(1:2) .eq. ' ') control_nml_file = 'yrec8.nml1'
      print *, ' '
      write(*,*) 'CONTROL namelist :  ',control_nml_file(1:len_trim(control_nml_file))

      call getarg(2, physics_nml_file)
      if (physics_nml_file(1:2) .eq. ' ') physics_nml_file = 'yrec8.nml2'
      write(*,*) 'PHYSICS namelist :  ',physics_nml_file(1:len_trim(physics_nml_file))

      open(unit=istand, file=control_nml_file, status='OLD')
      open(unit=irun, file=physics_nml_file, status='OLD')
      read(unit=istand, nml=control)
      read(unit=irun, nml=physics)
      close(istand)
      close(irun)
! stolr0/imax/nuse must keep their exact NAMELIST /physics/ spelling
! (see this file's naming note at the top), so intpar_lib's
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
! remap.f90 reads use_new_nuclear_rates to decide how to compute
! const_lib's cross-section-scale members.
      use_new_nuclear_rates = lnewnuc
      weak_screening_threshold = weakscreening
! envelope_overshoot_active/rotation_active/instability_transport_
! active/mass_accretion_rate/accreted_composition/use_mass_accretion:
! same reasoning, copied from their NAMELIST-spelled locals. Order
! relative to `call remap` doesn't matter for these -- remap.f90
! doesn't read any of them.
      envelope_overshoot_active = lovste
      rotation_active = lrot
      instability_transport_active = linstb
      mass_accretion_rate = dmdt0
      accreted_composition = compacc
      use_mass_accretion = lmdot
! use_itoh_neutrino_loss/use_new_turnover_timescale/calc_envelope_flag:
! same reasoning, copied from their NAMELIST-spelled locals.
      use_itoh_neutrino_loss = lnulos1
      use_new_turnover_timescale = lnewtcz
      calc_envelope_flag = lcalcenv
! clsun/crsun must likewise keep their NAMELIST spelling; copy into
! const_lib's solar_luminosity_cgs/solar_radius_cgs here (before
! setup/setups.f90 computes the rest of former common/const/ from
! these two).
      solar_luminosity_cgs = clsun
      solar_radius_cgs = crsun
! lexcom must likewise keep its NAMELIST spelling; copy into
! const_lib's use_extended_composition here.
      use_extended_composition = lexcom
! MHP 8/14 SUBROUTINE TO CONVERT MORE USER-FRIENDLY INPUT VARIABLES
! INTO THE VECTORS USED IN THE CODE (SUPERCEDES OLDER INPUTS)
      call remap
! MHP 06/13 Added memory of whether the choice of atmospheres has
! been changed during the run, and what the original setting was
      kttau0 = kttau
      lttau = .false.
! DBG WRITE OUT ENTIRE NAMELIST TO ISHORT
      write(short_file_unit,nml=physics)
      write(short_file_unit,nml=control)

! Post-process all CONTROL namelist vars that hold path values.
! Expand any placeholders found in the string with the value taken from a
! corresponding environment variable, if one is defined.
      call expand_value(falex06)
      call expand_value(fallard)
      call expand_value(fatm)
      call expand_value(fcondopacp)
      call expand_value(ffermi)
      call expand_value(ffirst)
      call expand_value(flast)
      call expand_value(fliv95)
      call expand_value(fmodpt)
      call expand_value(fopale06)
      call expand_value(fpatm)
      call expand_value(fpenv)
      call expand_value(fpmod)
      call expand_value(fpurez)
      call expand_value(fscomp)
      call expand_value(fscvh)
      call expand_value(fscvhe)
      call expand_value(fscvz)
      call expand_value(fshort)
      call expand_value(fsnu)
      call expand_value(fstor)
      call expand_value(ftrack)

! Create output directory as specified in the FTRACK value of CONTROL
! namelist if it doesn't already exist.
      shell_cmd = 'mkdir -p '
        ! find index of last '/' char. Use that to snip out the directory name.
      do i = len_trim(ftrack), 1, -1
          if (ftrack(i:i) .eq. '/') then
              last_slash_idx = i
              goto 1250
          endif
      end do
      1250 continue
      shell_cmd(len_trim(shell_cmd)+2:len_trim(ftrack(1:last_slash_idx))+len_trim(shell_cmd)) = ftrack(1:last_slash_idx)
      print *,"OUTPUT placed in :  ",ftrack(1:last_slash_idx)
      print *, ''
      call system(shell_cmd)


! JVS 02/11 Acoustic depth/ Asteroseismic glitch output. Puts output
! in the same directory as all other output, and names it with the
! same conventions
      if (lacout) then
            iclcd = 91
            short_prefix_len=index(fshort,'short')
            fcalcad=fshort(1:short_prefix_len-1)//'calcad'
            open(unit=iclcd, file=fcalcad, status='UNKNOWN')
            write(iclcd,*) 'Acoustic depth calculation output file'
            write(iclcd,*) 'age (Gyr),radius(cm),1/sound speed(s/cm),radius (CZ), 1/cs (CZ)                                         &
&       delad,gamma1,P, T, X'


!      IACAT = 92
!      FACAT=FSHORT(1:MRK-1)//'acatm'
!      OPEN(UNIT=IACAT, FILE=FACAT, STATUS='UNKNOWN')
!      WRITE(IACAT,*) 'Acoustic depth calculation output file: atmosphere integration'
!      WRITE(IACAT,*) 'age (Gyr),radius(cm),1/sound speed(s/cm),delad,gamma1,
!     * P, T, X'

!            IJLAST = 93
!            FJLAST=FSHORT(1:MRK-1)//'jlast'
!            OPEN(UNIT=IJLAST, FILE=FJLAST, STATUS='UNKNOWN')

!      IJVS = 94
!      FJVS=FSHORT(1:MRK-1)//'jvs'
!      OPEN(UNIT=IJVS, FILE=FJVS, STATUS='UNKNOWN')
!      WRITE(IJVS,*) 'The dels'

!      IJENT = 95
!      FJENT=FSHORT(1:MRK-1)//'ent'
!      OPEN(UNIT=IJENT, FILE=FJENT, STATUS='UNKNOWN')
!      WRITE(IJENT,*) 'Profiles: Age(Gyr), log(R), log(L), Log(LHe3),
!     * Conv flag, log(T), log(P), log(D)'

!      IJDEL = 96
!      FJDEL=FSHORT(1:MRK-1)//'del'
!      OPEN(UNIT=IJDEL, FILE=FJDEL, STATUS='UNKNOWN')
!      WRITE(IJDEL,*) 'Profiles: Age(Gyr), depth, del, del, del'
      endif


! JVS END
!
! G Somers 6/14, DEFINE SCALING COEFFICIENT FOR LI/BE CROSS SECTIONS
!          DEFAULT LI6 = 5.5 MeV b FROM FOWLER ET AL. 1967
!          DEFAULT LI7 = 52 keV b FROM ROLFS & KAVANAGH 1986
!           DEFAULT BE9(P,G)B10 (BE91) = 1.1 keV b INFERRED FROM FOWLER 1988
!           DEFAULT BE9(P,D)2HE4 (BE92) = 15,000 keV b FROM FOWLER ET AL. 1967
!           DEFAULT BE9(P,A)LI6 (BE93) = 15,000 keV b FROM FOWLER ET AL. 1967
      if (lxli6) then
          sli6 = xsli6/5.5d3
      endif
      if (lxli7) then
          sli7 = xsli7/5.2d1
      endif
      if (lxbe91) then
          sbe91 = xsbe91/1.1d0
      endif
      if (lxbe92) then
          sbe92 = xsbe92/1.5d4
      endif
      if (lxbe93) then
          sbe93 = xsbe93/1.5d4
      endif
! sli6/sli7/sbe91/sbe92/sbe93 must keep their NAMELIST spelling (see
! declaration above), so copy into const_lib's canonical names here.
      li6_rate_scale = sli6
      li7_rate_scale = sli7
      be9_pg_rate_scale = sbe91
      be9_pd_rate_scale = sbe92
      be9_palpha_rate_scale = sbe93
! G Somers END
! MHP 8/25 open relevant table as well as ensuring that only one is selected
!  Disable Older OPAL EOS's if a newer one is specified
      if (lopale06) then
          lopale01 = .false.
          lopale = .false.
         open(iopale, file=fopale06,status='OLD')
      endif
      if (lopale01) then
         open(iopale, file=fopale01,status='OLD')
         lopale = .false.
      else if(lopale) then
         open(iopale, file=fopale,status='OLD')
      endif
! 3/09 Disable older Alexander opacities if a newer one is specified
      if(lalex06)then
         lalex95 = .false.
         lkur90 = .false.
      endif
      if(lalex95)then
         lkur90 = .false.
      endif

      open(istor,file=fstor,form='FORMATTED',status='UNKNOWN')
      rewind(istor)
! G Somers 11/14 write the new header for the .store file, if LSTORE = TRUE.
      if(lstore)then
! JvS 08/25 Added stitched interior and envelope option      
         if(lstch)then
            lphhd = .true.
         else
            write(istor,version_fmt) yrecver, githash
            write(istor,1012)
         endif   
      endif
      1012 format('# Header Key',/,'# ModType    ModNum    #Shells    ', &
           &  'M/Msun    log(Teff)    log(L/Lsun)    log(M/gram)    Age/Gyr', &
           &  '    Timestep/yr    log(M_inner/gram)    log(M_outer/gram)',/, &
           &  '# JCORE  JENV  CMIXL  EOS  ATM  ALOK HIK  LPUREZ  COMPMIX', &
           &  '  LEXCOM  LDIFY  LDIFZ  LSEMIC  LOVSTC  LOVSTE  LOVSTM', &
           &  '  LROT  LINSTB  LJDOT0  LDISK  TDISK  PDISK  WMAX  LSTORE', &
           &  '  LSTATM  LSTENV  LSTMOD  LSTPHYS  LSTROT',/,'# PPI_lum', &
           &  '    PPII_lum    PPIII_lum    CNO_lum    3He_lum    ', &
           &  'OldNeu_lum    Grav_lum',/,'# ENV1-3   log(Teff)', &
           &  '    log(L/Lsun)    P_base    T_base    R_base    MatrixElements' &
           &  ,/)
! 1013 FORMAT('# JCORE  JENV  CMIXL  EOS  ATM  ALOK HIK  LPUREZ  COMPMIX',
!     1 '  LEXCOM  LDIFY  LDIFZ  LSEMIC  LOVSTC  LOVSTE  LOVSTM',
!     1 '  LROT  LINSTB  LJDOT0  LDISK  TDISK  PDISK  WMAX  LSTORE',
!     1 '  LSTATM  LSTENV  LSTMOD  LSTPHYS  LSTROT'
!     1 ,/)
! 1014     FORMAT(
!     1'MODEL SHELL MASS RADIUS LUMINOSITY PRESSURE TEMPERATURE DENSITY OMEGA ',
!     1'CONVECTIVE INTERIOR_PT H1 He4 METALS He3 C12 C13 N14 N15 O16 O17 O18 H2 Li6 Li7 ',
!     1'Be9 OPACITY GRAV DELR DEL DELAD V_CONV GAM1 HII HEII HEIII BETA ',
!     1'ETA PPI PPII PPIII CNO TRIPLE_ALPHA E_NUC E_NEU E_GRAV CP DLNRHODLNT A RP/RE FP ',
!     1'FT J/M MOMENT DEL_KE V_ES V_GSF V_SS VTOT ')
! G Somers END

      open(short_file_unit,file=fshort,form='FORMATTED',status='UNKNOWN')
      rewind(short_file_unit)
      if (ltrack) then
          open(unit=itrack,file=ftrack, form='FORMATTED', &
           &         status='UNKNOWN')
          rewind(itrack)
      endif

! SNU OUTPUT
      if(lsnu) then
          open(unit=isnu,file=fsnu, form='FORMATTED', &
           &         status='UNKNOWN')
          rewind(isnu)
      endif
      if(ldebug) then
            open(idebug,file=fdebug,form='FORMATTED', &
           &          status='UNKNOWN')
      end if
!     MHP Moved opening of FMILNE here, no need to do elsewhere
      if(lmilne)then
         open(unit=imilne,file=fmilne,form='FORMATTED', &
           &      status='UNKNOWN',access='APPEND')
      endif
!     MHP 10/02 LBNIN never set, ignore loop
!      IF (.NOT.LBNIN) THEN
         open(unit=ifirst,file=ffirst,form='FORMATTED',status='OLD')
!      END IF
      open(unit=ilast,file=flast,form='FORMATTED',status='UNKNOWN')
      open(unit=imodpt,file=fmodpt,form='FORMATTED',status='UNKNOWN')
!     OPEN ALL PULSE FILES
      if(lpulse) then
      open(iopmod, file=fpmod,status='UNKNOWN',form='FORMATTED')
      open(iopenv, file=fpenv,status='UNKNOWN',form='FORMATTED')
      open(iopatm, file=fpatm,status='UNKNOWN',form='FORMATTED')
      end if
! MHP 6/98
! MHP 8/25 Moved call from main to here for opening idyn
      if(lmonte)then
         open(unit=idyn,file=fdyn,form='FORMATTED',status='OLD')
         open(imonte1, file=fmonte1,status='UNKNOWN',form='FORMATTED')
         open(imonte2, file=fmonte2,status='UNKNOWN',form='FORMATTED')
      endif
!     MHP 8/25 Moved opening of conductive opacity and EoS tables here, to avoid complicated passages of declared variables.
      if(lcondopacp)then
         open(icondopacp,file=fcondopacp,status='OLD')
      endif
      if(liso) then
         open(iiso, file=fiso,status='UNKNOWN', form='FORMATTED')
      endif
      if(lsemic)then
         if(lovstc.or.lovste.or.lovstm)then
            write(short_file_unit,2)lsemic,lovste,lovstc,lovstm
      2       format(1x,'ERROR IN SUBROUTINE PARMIN'/'SEMI-CONVECTION', &
           &  ' AND OVERSHOOT FLAGS BOTH TURNED ON'/'FLAGS LSEMIC',l2, &
           &  ' OVERSHOOT - CORE,ENVELOPE,INTERMEDIATE-',3l2/'RUN STOPPED')
            stop
         endif
      endif
      write(short_file_unit,1)(hpttol(i),i=1,12),alphae,alphac,linstb,ljdot0, &
           &               alfa,fk,fw,fc,fo,fmu,rcrit
      1 format(1x,'PT TOL',12f6.3/1x,'O.S.ENV',f6.3,' O.S.CORE',f6.3, &
           &        ' LINSTB ',l1,' LJDOT ',l1,' WIND IND.',f6.3,' FK', &
           &        1pe8.2/1x,' FV',0pf5.2,' FC',f5.2,' COUPLING' &
           &        ,f6.3, ' F MU',f5.2,' RCRIT',f9.1)
      tenv = 0.5d0*(tenv0 + tenv1)
      if(lrot) then
         lnew0 = .true.
         if(walpcz.lt.-2.0d0) walpcz = -2.0d0
         if(walpcz.gt.0.0d0) walpcz = 0.0d0
! JNT 09/2025 LSOLID OVERWRITES IMPJMOD
         if(lsolid) impjmod = 1
!CCCCC OLD OR NEW WINDLAW.
         if(.not.lmwind)then
!CCCCC INSTRUCTIONS FOR THE OLD WINDLAW
!CCCCC
!CCCCC SET UP COEFFICIENTS FOR ANGULAR MOMENTUM LOSS VIA WINDS.
!CCCCC GIVEN THE INDEX ALFA, THE FORMULA FOR JDOT IS
!CCCCC JDOT = FK*2.036D33*1.452D9**ALFA*(MDOT/1.0D-14)**(1-2ALFA/3)
!CCCCC *OMEGA**(1+4ALFA/3)*(R/RSUN)**(2-ALFA)*(M/MSUN)**-ALFA/3
!CCCCC EXMD = EXPONENT OF MDOT TERM; EXW SAME FOR OMEGA;EXR FOR R;EXM FOR M.
!CCCCC FK IS A FUDGE FACTOR,SET TO 1 TO REPRODUCE THE OBSERVED SOLAR ANGULAR
!CCCCC MOMENTUM LOSS.
            one_third = 1.0d0/3.0d0
            two_thirds = 2.0d0/3.0d0
            constfactor = fk*2.036d33*1.452d9**alfa
            exmd = 1.0d0 - two_thirds*alfa
            exw = 1.0d0 + 2.0d0*two_thirds*alfa
            exr = 2.0d0 - alfa
            exm = -one_third*alfa
!
! INSTRUCTIONS FOR THE NEW WINDLAW
!
! SET UP INDICES FOR LOSS LAW IN TERMS OF PMM A, B, C, M.
! INCLUDE A PMMD AS WELL, FOR TURNING TAUCZ and CENTRIFUG DEPENDENCE ON/OFF.
! EXCEN = EXPONENT FOR CENTRIFUGAL TERM, {K2^2/(K2^2+0.5*W^2 R^3/GM)}^PMMM
!
! AWIND = 'K97' ENFORCES KRISHNAMURTHI (1997) WIND LAW.
! AWIND = 'V13' ENFORCES VAN SADERS + PINSONNEAULT (2013) WIND LAW, AS ADAPTED
!         FROM MATT ET AL. (2008,2012).
! AWIND = 'CUS' ADOPTS THE M-D VALUES GIVEN IN THE NAMELIST.
!
!
! THEN CONVERT THE PMM VALUES INTO EXPONENTS FOR THE TORQUE CALCULATION.
! BASIC EQUATIONS, WHERE EVERY TERM IS SCALED TO SOLAR, ARE:
!
! Jdot = w * Bmag^4m * Mdot^1-2m * R^5m+2 * M^-m * Fcen^md
! Bmag = Pphot^0.5 * (w * R^c)^b * Tcz^d
! GM/R * Mdot = (Lx/Lbol * Lbol)^a -> Mdot = R * M^-1 * Lbol * w^a
! ...where Lx/Lbol is assumed to scale with w to the power of a.
!
! THIS GIVES:
! JDOT ~ w ^ 1+a-2ma+4mb
!      ~ Tcz ^ 4md
!      ~ R ^ 3+3m+4mbc
!      ~ M ^ m-1
!      ~ Lbol ^ 1-2m
!      ~ Pphot ^ 2m
!      ~ Fcen ^ md
!
! For moment, ignore GM/R term in Mdot. This makes:
! JDOT ~ R ^ 2+5m+4mbc
! JDOT ~ M ^ -m
!
         else
            if(awind.eq.'K97')then
               pmmm = alfa/3.0
! MHP 8/17 CORRECT DEFAULT FOR A = 2, NOT 0
               pmma = 2.0
               pmmb = 1.0
               pmmc = 2.0
               pmmd = 0.0
               lbscale = .false.
            elseif(awind.eq.'V13')then
               pmmm = 0.22
               pmma = 2.0
               pmmb = 1.0
               pmmc = 0.0
               pmmd = 1.0
               lbscale = .true.
            endif
            exw   = 1.0d0 + pmma - 2.0d0*pmma*pmmm + 4.0d0*pmmm*pmmb
! G Somers 8/17 ZERO'D OUT EXTAU. TAUCZ TERM NOW COMPUTED IN
! MWIND/MCOWIND, NOT IN AMCALC.
!            EXTAU = 4.0D0*PMMM
            extau = 0.0
! JvS 09/25 REMOVED TYPO IN EXR = 2.0D0+5.0D0*PMMM-4.0D0*PMMM*PMMB*PMMC
            exr   = 2.0d0+5.0d0*pmmm-4.0d0*pmmm*pmmc
            exm   = -pmmm
            exl   = 1.0d0 - 2.0d0*pmmm
! G Somers 11/17, ADDED PMMD TO SWTICH OFF IN K97 FORM.
            expr  = 2.0d0*pmmm*pmmd
            excen = pmmm*pmmd
! INITIALIZE CONSTANT FACTOR FOR CENTRIFUGAL TERM
            c_2 = 0.0506
! SET THE CONSTANT FACTOR
            constfactor = fk*pmmjd/pmmsolw**exw
! IF RELEVANT RESET THE SATURATION THRESHOLD IN
! TERMS OF THE SOLAR ROTATION RATE.  WMAX_SUN<1000
! INDICATES SATURATION (AT THE SUN), SO
! WSAT = WMAX_SUN*PMMSOLW
            if(wmax_sun.lt.1.0e3)then
               wmax = wmax_sun*pmmsolw
            endif
         endif
      endif
!     WINDLAW END
!
      parmin_ln10 = dlog(10.0d0)
      if(lnewcp) then
       lrel = .true.
       if(atmp.eq.'ABS') lrel = .false.
! DECIDE WHICH ELEMENT IN ARRAY HCOMP TO BE RESCALED
! USING CHARACTER ARRAY AID AND INPUT CHARACTER VARIABLE ANEWCP
       do 10 i = 1,12
          if(anewcp.eq.element_id(i)) then
! INEWCP IS THE INDEX OF THE ELEMENT BEING ALTERED
            inewcp = i + 3
            goto 30
          endif
      10    continue
! ANEWCP NOT A RECOGNIZED ELEMENT
       lnewcp = .false.
       write(short_file_unit,20) anewcp
      20    format(1x,'VARIABLE',a4,1x,'NOT A RECOGNIZED ELEMENT'/1x, &
           &    'RESCALING NOT PERFORMED')
      30    continue
      endif
      lmixture = .false.
      lisotope = .false.
      if(isetmix.eq.1)then
! IF DEFAULT MIX (GS98) CNOFRACS ARE ALREADY SET.
!         IF(AMIX.EQ.'GS98')THEN
!            LMIXTURE = .TRUE.
!            GOTO 606
!         ENDIF
! FOR A CUSTOM MIX,DISABLE IF THE SUM OF CNO MASS FRACTIONS EXCEEDS ONE
! OR IF ANY MASS FRACTION IS NEGATIVE
         if(amix.eq.'CUS')then
            if(frac_c.lt.0.0d0.or.frac_n.lt.0.0d0.or.frac_o.lt.0.0d0)then
               write(*,591)frac_c,frac_n,frac_o
               write(short_file_unit,591)frac_c,frac_n,frac_o
      591          format('NEGATIVE INPUT CNO FRACTION ',3e12.4, &
           &              ' MIX NOT MODIFIED')
               goto 602
            endif
            sum_frac=frac_c+frac_n+frac_o
            if(sum_frac.ge.1.0d0)then
               write(*,598)frac_c,frac_n,frac_o
               write(short_file_unit,598)frac_c,frac_n,frac_o
      598          format('INPUT CNO FRACTION ',3e12.4, &
           &              ' EXCEEDS 1. MIX NOT MODIFIED')
               goto 602
            endif
! VALID MIXTURE, USE CUSTOM ENTRIES FROM .NML1
            lmixture = .true.
            goto 606
         endif
! SEARCH THROUGH OTHER VALID MIXTURE ENTRIES;IF FOUND,ASSIGN
!         DO I = 2,4
         do i = 1,4
            if(amix.eq.mixture_id_table(i))then
               zxmix = zx_mix_table(i)
               frac_c = frac_c_table(i)
               frac_n = frac_n_table(i)
               frac_o = frac_o_table(i)
               lmixture = .true.
               goto 606
            endif
         end do
!     NO VALID MIX SPECIFIED
         write(*,589)amix
         write(short_file_unit,589)amix
      589    format('DESIRED CNO MIXTURE ',a8,' NOT FOUND. MIX NOT ALTERED.')
      endif
      606 if(lmixture)then
         write(*,604)amix,frac_c,frac_n,frac_o
         write(short_file_unit,604)amix,frac_c,frac_n,frac_o
      604    format('CNO MIXTURE ',a8,' C ',e12.4,' N ',e12.4,' O ', &
           &         e12.4,' APPLIED TO STARTING MODEL.')
      endif
!     CHECK IF ISOTOPE RATIOS NEED TO BE ALTERED
      602 if(isetiso.eq.1)then
! FOR A CUSTOM MIX,DISABLE IF THE SUM OF CNO MASS FRACTIONS EXCEEDS ONE
! OR IF ANY MASS FRACTION IS NEGATIVE
         if(aiso.eq.'CUS')then
            if(r12_13.lt.0.0d0 .or. r14_15.lt.0.0d0 .or. r16_17.lt.0.0d0 &
           &  .or. r16_18.lt.0.0d0 .or. xh2_ini.lt.0.0d0 .or. xhe3_ini.lt.0.0d0 &
           &  .or. xli6_ini.lt.0.0d0 .or. xli7_ini.lt.0.0d0 .or. &
           &  xbe9_ini.lt.0.0d0.or.xb10_ini.lt.0.0d0 .or.xb11_ini.lt.0.0d0)then
               write(*,596)r12_13,r14_15,r16_17,r16_18, &
           &  xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini,xb10_ini,xb11_ini
               write(short_file_unit,596)r12_13,r14_15,r16_17,r16_18, &
           &  xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini,xb10_ini,xb11_ini
      596          format('NEGATIVE INPUT ISOTOPE RATIO OR LIGHT ELEMENT' &
           &   ' MASS FRACTION ',11e12.4, &
           &              ' MIX NOT MODIFIED')
               goto 603
            endif
            sum_frac= xh2_ini+xhe3_ini+xli6_ini+xli7_ini+xbe9_ini+ &
           &                 xb10_ini+xb11_ini
            if(sum_frac.ge.1.0d0)then
               write(*,595)xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini, &
           &  xb10_ini,xb11_ini
               write(short_file_unit,595)xh2_ini,xhe3_ini,xli6_ini,xli7_ini, &
           &  xbe9_ini,xb10_ini,xb11_ini
      595          format('SUM OF LIGHT ELEMENT MASS FRACTIONS EXCEEDS 1', &
           &  11e12.4,' MIX NOT MODIFIED')
               goto 603
            endif
!     CURRENTLY THERE ARE ONLY 2 VALID OPTIONS - THE DEFAULT (L21) OR
! A CUSTOM MIXTURE (CUS) - IF NEITHER IS TRUE, EXIT
         else if(aiso.ne.'L21')then
            goto 603
         endif
!     PASSED ALL CHECKS - EITHER THE DEFAULT OR THE CUSTOM SETTINGS WILL BE APPLIED
         lisotope = .true.
      endif
      if(lisotope)then
         write(*,605)aiso,r12_13,r16_18, &
           &  xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini
         write(short_file_unit,605)aiso,r12_13,r16_18, &
           &  xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini
      605    format('ISOTOPE AND LIGHT ELEMENT MIXTURE ',a8,' C12/C13 ', &
           &    e12.4,' O16/O18 ',e12.4,' H2 ',e12.4,' HE3 ',e12.4,' LI6 ', &
           &    e12.4,' LI7 ',e12.4,' BE9 ',e12.4,' APPLIED TO STARTING MODEL.')
      endif
! DBG 12/95 ENSURE CORRECT PARAMETERS FOR Z DIFFUSION
      603 if (ldifz) then
           ldify=.true.
         lthoul=.true.
      end if

! G SOMERS 04/15; ENSURE CORRECT PARAMETERS FOR LIGHT ELEMENT DIFFUSION.
      if(ldifli)then
          lnewdif=.true.
          ldify=.true.
          ldifz=.true.
          lthoul=.true.
          lthoulfit=.false.
          ilambda=4
      endif

!     WRITE OUT RUN PARAMETERS.

      write(short_file_unit,50)xenv0,zenv0
      50 format(30x,'RUN DATA VALUES'/3x, &
           &        'LINE  1     XENV0     ZENV0       ZSI'/2x, &
           &        'STANDARD       N/A       N/A  0.00E+00'/3x, &
           &        'CURRENT',1p2e10.0,2x,e8.2)
      if(npoint.le.0) npoint = 9999
      write(short_file_unit,70) ldebug,lcorr,npoint,lmilne,ltrack,lstore,lstpch
      70 format(3x,'LINE  2    LDEBUG     LCORR    NPOINT    LMILNE    LTRA                                                            &
&CK    LSTORE    LSTPCH'/2x,'STANDARD',2(9x,'T'),6x, &
           & '9999',9x,'F',2(9x,'T'),9x,'T'/3x,'CURRENT',2(9x,l1),6x, &
           & i4,3(9x,l1),9x,l1)
      if(npenv.le.0) npenv = 9999
      write(short_file_unit,90) lscrib,lstatm,lstenv,lstmod,lstphys,lstrot
      90 format(3x,'LINE  3    LSCRIB    LSTATM    LSTENV    LSTMOD',5x, &
           & 'LSTPHYS    LSTROT'/2x,'STANDARD',9x,'T',2(9x,'F'), &
           & 9x,'T',2(9x,'F'),9x,'2'/3x,'CURRENT',5(9x,l1),9x,l1)
      write(short_file_unit,110)lenvg,atmstp,envstp
      110 format(3x,'LINE 4     LENVG    ATMSTP    ENVSTP'/2x,'STANDARD',7x, &
           &        'N/A',6x,'0.50',6x,'0.50'/3x,'CURRENT',9x,l1,2(4x,f6.3))
      if(nprtmod.le.0) nprtmod = 9999
      if(nprtpt.le.0) nprtpt = 9999
      if(pulse_gyre_interval.lt.0) pulse_gyre_interval = 0
      write(short_file_unit,130) nprtmod,nprtpt
      130 format(3x,'LINE  4  NPRTMOD    NPRTPT'/2x,'STANDARD', &
           & 1(7x,'N/A'),9x,'5'/3x,'CURRENT',2(6x,i4))

!     SPIT OUT NAMELIST VARIABLES TO ISHORT

      write(short_file_unit,25) (tcut(j),j=1,5),tscut,tenv0,tenv1,tgcut
      25 format(3x,'LINE  2    TCUT-  E  TCUT- PP  TCUT-CNO  TCUT-                                                                     &
&3A  TCUT- NU TCUT-SAHA     TENV0     TENV1     TGCUT'/2x, 'STANDAR                                                            &
&D',9x,'6.50',6x,'6.50',6x,'6.82',6x,'7.70',6x,'7.50',6x, &
           & '6.00',6x,'3.00',6x,'9.00',6x,'6.90'/3x,'CURRENT', &
           & 9(5x,f5.2))
      write(short_file_unit,35) atmerr,atmmax,atmd0,enverr,envmax,envmin
! MHP 10/02 obsolete variables removed
!      WRITE(ISHORT,35) NIATM,ATMERR,ATMMAX,ATMD0,NIENV,ENVERR,ENVMAX,
!     *ENVMIN
      35 format(3x,'LINE  3     NIATM    ATMERR    ATMMAX     ATMD0',5x, &
           & 'NIENV    ENVERR    ENVMAX    ENVMIN'/2x,'STANDARD',9x,'8',2x, &
           & '3.00E-04  5.00E-01  1.00E-10',8x,'10  3.00E-04  5.00E-01  2.50E-0                                                            &
&1'/3x,'CURRENT',10x,3(1pe10.2),10x,3(1pe10.2))
!     *1'/3X,'CURRENT',7X,I3,3(1PE10.2),7X,I3,3(1PE10.2))
      write(short_file_unit,45) tridt,tridl/tridt/tridl
      45 format(3x,'LINE  4     TRIDT     TRIDL    LSENV0',5x, &
           & 'SENV0'/2x,'STANDARD',6x,'0.01',6x,'0.08',9x,'F',2x, &
           & '1.00E-07'/3x,'CURRENT',2f10.4)
      write(short_file_unit,55)(htoler(5,j),j=1,2),((htoler(i,j),i=1,4),j=1,2)
      55 format(3x,'LINE  5 TOL.RHS-P TOL.RHS-T MIN.COR-P MIN.COR-T MIN.COR                                                            &
&-R MIN.COR-L MAX.COR-P MAX.COR-T MAX.COR-R MAX.COR-L'/2x,'STANDARD                                                            &
&  3.00E-05  2.50E-06  6.00E-05  4.50E-05  3.00E-05  9.00E-05  9.00                                                            &
&E-01  5.00E-01  5.00E-01  2.00E+00'/3x,'CURRENT',10(2x,1pe8.2))
      write(short_file_unit,65)(hpttol(j),j=1,8)
      65 format(3x,'LINE  6  D(S)-MIN  D(S)-MAX   FLAG-DX   FLAG-DZ    MAX                                                             &
&DP MAX DL/LT    MAX DX    MAX DZ'/2x,'STANDARD  1.00E-08  8.00E-02                                                            &
&  5.00E-02  1.00E+00  5.00E-02  2.00E-02  1.00E+00  1.00E+00'/3x, &
           & 'CURRENT',8(2x,1pe8.2))
      write(short_file_unit,75)lnewcp,anewcp,lrel,xnewcp
      75 format(3x,'LINE  7    LNEWCP    ANEWCP      LREL    XNEWCP'/2x, &
           & 'STANDARD',9x,'F',7x,'N/A',9x,'T',7x,'N/A'/3x,'CURRENT',9x,l1,7x, &
           & a3,9x,l1,1pe10.2)
      write(short_file_unit,85)acfpft,itfp1,itfp2
      85 format(3x,'LINE  8    ACFPFT     ITFP1     ITFP2'/2x,'STANDARD', &
           & 1x,'1.000E-20',9x,'2',8x,'20'/3x,'CURRENT',1pe10.3,6x,i4,6x,i4)
      write(short_file_unit,105) niter1,niter2,fcorr0,fcorri
      105 format(32x,'RUN DATA VALUES'/3x, &
           & 'LINE  1    NITER1    NITER2    FCORR0    FCORRI'/2x, &
           & 'STANDARD',9x,'2',8x,'20',6x,'0.80',6x,'0.10'/3x, &
           & 'CURRENT',2(6x,i4),2(5x,f5.2))
      write(short_file_unit,145)(atime(i),i=1,3),atime(7)
      145 format(3x,'LINE  5 XCORE MIN  DEL.XCORE FRAC.XCORE DEL.XSHELL'/2x, &
           & 'STANDARD     0.001     0.020     0.500     0.100'/3x,'CURRENT', &
           & 4(4x,f6.3))
      write(short_file_unit,155)(atime(i),i=4,6)
      155 format(3x,'LINE  7 DEL.YCORE  FRAC.YCORE DEL.YSHELL'/2x,'STANDARD                                                             &
&    0.020     0.300    0.0015'/3x,'CURRENT',2(5x,f5.3),4x,f6.4)
      write(short_file_unit,165) lkuthe
      165 format(3x,'LINE  8    LKUTHE'/,2x, &
           & 'STANDARD',9x,'F'/3x,'CURRENT',9x,l1)
      write(short_file_unit,175) cmixl,dpenv,lovstc,alphac,lovste,alphae
      175 format(3x,'LINE  9   CMIXL     DPENV    LOVSTC    ALPHAC    LOVSTE                                                            &
&    ALPHAE'/2x,'STANDARD',7x,'N/A',6x,'1.00',2(9x,'F',6x,'0.00')/ &
           & 3x,'CURRENT',2(5x,f5.2),2(9x,l1,5x,f5.2))
      write(short_file_unit,185) lnew0,lexcom
      185 format(3x,'LINE 10   LNEW0    LEXCOM'/2x, &
           & 'STANDARD',2(9x,'F')/3x,'CURRENT',2(9x,l1))
      write(short_file_unit,195) lrot,walpcz,linstb
      195 format(3x,'LINE 11      LROT    WALPCZ    LINSTB'/2x,'STANDARD', &
           & 7x,'N/A',6x,'0.00',7x,'N/A'/3x,'CURRENT',9x,l1,5x,f5.2,9x,l1)
      if(kttau .eq. 0) then
           write(short_file_unit, 197)
      else if (kttau .eq. 1) then
           write(short_file_unit, 198)
      else if (kttau .eq. 2) then
           write(short_file_unit, 1999)
      else if (kttau .eq. 3) then
           write(short_file_unit, 1888)
      else if (kttau .eq. 4) then
           write(short_file_unit, 1889)
! JNT 6/14 ADD FOR NEW KURUCZ/CASTELLI ATMOSPHERE TABLES
      else if (kttau .eq. 5) then
           write(short_file_unit, 1887)
      end if
      197 format(' USING EDDINGTON T-TAU RELATION.')
      198 format(' USING KRISHNA-SWAMY T-TAU RELATION.')
      1999 format(' USING HARVARD-SMITHSONIAN REFERENCE ATMOSPHERE')
      1888 format(' USING KURUCZ ATMOSPHERE TABLE')
      1889 format(' USING ALLARD ATMOSPHERE TABLE')
      1887 format(' USING KURUCZ/CASTELLI ATMOSPHERE TABLE')

! DBG PULSE
      if (lpulse) then
          write(short_file_unit,196)
      196     format(/,' CALCULATE PULSATION OUTPUT ON LAST MODEL')
      end if
      if (lpurez) then
          write(short_file_unit,*) ' USING PURE C AND N OPACITY TABLES'
      end if

      write(imodpt,310) descrip(1),  descrip(2)

      write(short_file_unit,314)
      write(short_file_unit,version_fmt) yrecver, githash
      write(short_file_unit,310) descrip(1),  descrip(2)
      310    format('# DESCRIPTION OF RUN:',a80,/, '#',9x,'  ',8x,': ', &
           &           a80,/,'#', 100('='))

      if(ltrack) then

         write(itrack,314)
      314    format('#',/,'#',100('='))
         write(itrack,version_fmt) yrecver, githash
         write(itrack,320) descrip(1), descrip(2)
      320    format('# DESCRIPTION OF RUN:',a80,/, '#',9x,'  ',8x,': ', &
           &           a80,/,'#', 100('='))

      endif
      write(short_file_unit,323)
      323 format(' USING OSCILATORY SPLINE INTERPOLATION IN HPOINT')

!     INTERPRET RUN FROM SEQUENCE OF "KIND" CARDS

      write(short_file_unit,200)
      200 format(/35x,'RUN CARDS'/)

      lfirst(1) = .true.

!     RUN LOOP
      do 1000 nkind=1, numrun
! READ IN NMODLS AND MODEL SOURCE(MEMORY OR FIRST MODEL)-SAME FOR ALL
       iresca(nkind) = kindrn(nkind)
       if(kindrn(nkind).eq.1) then
! EVOLVE CARD
! MHP 10/24 GENERALIZE STOP CONDITIONS
!          LENDAG(NKIND) = ENDAGE(NKIND).GT.0D0
            if(endage(nkind).gt.0.0d0 .or. end_dcen(nkind).gt.0.0d0 &
           &  .or. end_xcen(nkind).gt.0.0d0 .or. end_ycen(nkind).gt.0.0d0)then
               lendag(nkind)=.true.
            else
               lendag(nkind)=.false.
            endif
          lsetdt(nkind) = setdt(nkind).gt.0d0
            if (nmodls(nkind).gt.0) then
          if (lfirst(nkind)) then
             write(iowr,350) nkind,nmodls(nkind)
             write(short_file_unit,350) nkind,nmodls(nkind)
      350          format(/1x,'RUN #',i3,'   EVOLVE ',i5, &
           &          ' MODELS, STARTING', &
           &          ' WITH THE INPUT "FIRST MODEL".')
          else
             write(iowr,351) nkind, nmodls(nkind)
             write(short_file_unit,351) nkind, nmodls(nkind)
      351          format(/1x,'RUN #',i3,'   EVOLVE ',i5, &
           &          ' MODELS, STARTING', &
           &          ' WITH THE PREVIOUS RUN''S LAST MODEL.')
          end if
! GENERALIZE STOP CONDITIONS
          if(lendag(nkind).or.lsetdt(nkind)) then
             write(iowr,370)lendag(nkind),lsetdt(nkind), &
           &               endage(nkind), setdt(nkind),end_dcen(nkind), &
           &          end_xcen(nkind),end_ycen(nkind)
             write(short_file_unit,370)lendag(nkind),lsetdt(nkind), &
           &          endage(nkind), setdt(nkind),end_dcen(nkind), &
           &          end_xcen(nkind),end_ycen(nkind)
      370          format(1x,'EVOLVE TO AGE ',l1,' SET DELT ', &
           &               l1,' FINAL AGE ', e9.2,' FIXED TSTEP ',e9.2, &
           &   ' CENTRAL D ',e10.4,' CENTRAL X ',e12.4,' CENTRAL Y ',e12.4)
          endif
            end if
       else if(kindrn(nkind).eq.2) then
! RESCALE CARD:  RESCALE STARTING MODEL
! QUANTITIES TO BE RESCALED STORED IN ARRAY RESCALE(4,50)
! WHERE THE ELEMENTS MASS,X,Z,CORE MASS ARE STORED IN ORDER
          rescal(1,nkind) = rsclm(nkind)
          rescal(2,nkind) = rsclx(nkind)
          rescal(3,nkind) = rsclz(nkind)
          rescal(4,nkind) = rsclcm(nkind)
            if (nmodls(nkind) .gt. 0) then
          if (lfirst(nkind)) then
             write(iowr,450) nkind
             write(short_file_unit,450) nkind
      450          format(/1x,'RUN #',i3, &
           &          '   RESCALE THE INPUT MODEL: "FIRST MODEL".')
          else
             write(iowr,451) nkind
             write(short_file_unit,451) nkind
      451          format(/1x,'RUN #',i3, &
           &          '   RESCALE THE PREVIOUS RUN''S LAST MODEL.')
          end if
          write(iowr,452) nmodls(nkind),(rescal(i,nkind),i = 1,4)
          write(short_file_unit,452) nmodls(nkind), &
           &       (rescal(i,nkind),i = 1,4)
      452       format(1x,'RELAX RESCALED MODEL',i3, &
           &       ' TIMES. RESCALE THE FOLLOW', &
           &       'ING(0=USE CURRENT VALUE):'/1x,'MASS ', &
           &       f9.6,3x,'X',f9.6,3x,'Z', &
           &       f9.6,3x,'CORE MASS',f9.6)
            end if
         else if(kindrn(nkind).eq.3) then
! RESCALE AND EVOLVE CARD:  RESCALE STARTING MODEL
! QUANTITIES TO BE RESCALED STORED IN ARRAY RESCALE(4,50)
! WHERE THE ELEMENTS MASS,X,Z,CORE MASS ARE STORED IN ORDER
            rescal(1,nkind) = rsclm(nkind)
            rescal(2,nkind) = rsclx(nkind)
            rescal(3,nkind) = rsclz(nkind)
            rescal(4,nkind) = rsclcm(nkind)
            if (lfirst(nkind)) then
               write(iowr,550) nkind
               write(short_file_unit,550) nkind
      550          format(/1x,'RUN #',i3, &
           &          '   RESCALE & EVOLVE THE INPUT MODEL: "FIRST MODEL".')
            else
               write(iowr,451) nkind
               write(short_file_unit,451) nkind
!   551          FORMAT(/1X,'RUN #',I3,
!      1         '   RESCALE & EVOLVE THE PREVIOUS RUN''S LAST MODEL.')
            end if
            write(iowr,452) nmodls(nkind),(rescal(i,nkind),i = 1,4)
            write(short_file_unit,452) nmodls(nkind), &
           &       (rescal(i,nkind),i = 1,4)
       end if
         if(rescal(3,nkind).ge.0.0d0)  zenv=rescal(3,nkind)
      1000 continue
      return

end subroutine parmin


!----------------------------------------------------------------------
! expand_value
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original parmin.f (this subroutine lived in the same source file);
! only variable names, source form, and comment style were updated.
!
! Replaces any defined "{YREC_XXX}" placeholder string in the passed
! namelist path variable with the value of the corresponding
! environment variable (falling back to a built-in default path if
! that environment variable is not set).
subroutine expand_value(path_value)

      implicit none

      integer, parameter :: n_env_vars = 3  ! Number of possible env. vars.
      character(len=256) :: path_value
      character(len=256) :: temp_value
      character(len=256) :: placeholder_names(n_env_vars)
      character(len=256) :: default_values(n_env_vars)
      character(len=256) :: placeholder
      integer :: env_val_len, orig_val_len, placeholder_len, default_len
      integer :: i
! Each placeholder (env var name enclosed in curly braces) to be
! supported for expansion in namelists, along with default value in the
! case where the var is not defined in the execution environment.
! The number of assignment pairs here must match the value of the
! parameter n_env_vars, above.
! Each placeholder (env var name enclosed in curly braces) to be
! supported for expansion in namelists, along with default value in the
! case where the var is not defined in the execution environment.
! The number of assignment pairs here must match the value of the
! parameter NUM_ENVVARS, above.
      placeholder_names(1) = "{YREC_INPUT}"
      default_values(1) = "../../input"

      placeholder_names(2) = "{YREC_START}"
      default_values(2) = "../../startmodels"

      placeholder_names(3) = "{YREC_OUTPUT}"
      default_values(3) = "output"

      do 5000 i=1, n_env_vars
        placeholder_len = len_trim(placeholder_names(i))
        default_len = len_trim(default_values(i))
        placeholder = placeholder_names(i)(1:placeholder_len)

        if (path_value(1:placeholder_len) .eq. placeholder) then    ! If placeholder
            call getenv(placeholder(2:placeholder_len-1), temp_value)
            env_val_len = len_trim(temp_value)
            orig_val_len = len_trim(path_value)
            if (env_val_len .ne. 0) then       ! If env var set
                temp_value(env_val_len+1:env_val_len+orig_val_len-placeholder_len) = path_value(placeholder_len+1:orig_val_len)
                  !print *,"Override: ",TEMP(1:LEN_TRIM(TEMP))
            else                          ! No env var. Use default path.
                temp_value(1:default_len) = default_values(i)(1:len_trim(default_values(i)))
                temp_value(default_len+1:default_len+orig_val_len) = path_value(placeholder_len+1:orig_val_len)
                  !print *,"default: ",TEMP(1:LEN_TRIM(TEMP))
            end if
            path_value = temp_value
        end if
      5000 continue

end subroutine expand_value
