!----------------------------------------------------------------------
! remap
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original remap.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! MHP 6/2014
! Subroutine to remap new, more intuitively named, namelist parameters
! onto existing code parameters.
subroutine remap
      implicit none

! PARAMETERS NT AND NG FOR TABULATED SURFACE PRESSURES.
      integer, parameter :: nt = 57, ng = 11
! PARAMETERS NTA AND NGA FOR TABULATED ALLARD MODEL SURFACE PRESSURES.
      integer, parameter :: nta = 54, nga = 5
      integer, parameter :: nts = 63, nps = 76
      integer, parameter :: json = 5000

! common/ctol/: only chi_grid_scale (originally HPTTOL) is used here,
! being assigned from the TOL_* namelist parameters below. Naming
! matches mixgrid.f90.
      double precision :: htoler(5,2), fcorr0, fcorri, fcorr, &
           chi_grid_scale(12)
      integer :: niter1, niter2, niter3
      common/ctol/ htoler, fcorr0, fcorri, fcorr, chi_grid_scale, niter1, &
           niter2, niter3

! common/ctlim/: only atime is used (set) here. Naming matches
! eqburn.f90.
      double precision :: atime(14), tcut(5), saha_log10t_cutoff, tenv0, &
           tenv1, tenv, tgcut
      common/ctlim/ atime, tcut, saha_log10t_cutoff, tenv0, tenv1, tenv, tgcut

! common/difus/: only dtdif is used (set) here. Naming matches
! dadcoeft.f90.
      double precision :: dtdif, convergence_tolerance
      integer :: itdif1, max_iterations
      common/difus/ dtdif, convergence_tolerance, itdif1, max_iterations

! common/ct2/: max permitted change in omega per global timestep. Not
! referenced in any already-converted file.
      double precision :: max_domega_global
      common/ct2/ max_domega_global

! common/ct3/: master flag for whether the ATIME(8)-(11) structure-change
! timestep limits are applied (only used in the pre-MS; disabled
! automatically for MS stars). Not referenced in any already-converted
! file.
      logical :: use_structure_dt_limits
      common/ct3/ use_structure_dt_limits

! MHP 8/96 CROSS SECTIONS PUT IN COMMON BLOCK.
! MHP 6/14 DERIVATIVES ADDED
! common/cross/: cross-section scale factors, all set here. Naming
! matches rates.f90.
      double precision :: cross_section_scale(17), qs0e_scale(8), &
           qqs0ee_scale(8), o16_gamma_scale, c12_alpha_scale
      logical :: use_new_nuclear_rates
      common/cross/ cross_section_scale, qs0e_scale, qqs0ee_scale, &
           o16_gamma_scale, c12_alpha_scale, use_new_nuclear_rates

! common/newcross/: user-supplied nuclear reaction S-factors (and
! first/second derivative ratios relative to the Adelberger et al. 1998
! Solar Fusion I values) that feed the cross_section_scale/qs0e_scale/
! qqs0ee_scale vectors above. Not referenced in any already-converted
! file.
      double precision :: s0_pp, s0_he3he3, s0_he3he4, s0_p_c12, s0_p_c13, &
           s0_p_n14, s0_p_o16, s0_pep, s0_be7_electron, s0_be7_p, s0_hep, &
           s0_n15_p_c12_branch, s0_n15_p_o16_branch, &
           s0p_pp, s0p_he3he3, s0p_he3he4, s0p_p_c12, s0p_p_c13, s0p_p_n14, &
           s0p_p_o16, s0pp_p_c12, s0pp_p_c13, s0pp_p_o16, s0p_be7_p, &
           s0pp_be7_p
      common/newcross/ s0_pp, s0_he3he3, s0_he3he4, s0_p_c12, s0_p_c13, &
           s0_p_n14, s0_p_o16, s0_pep, s0_be7_electron, s0_be7_p, s0_hep, &
           s0_n15_p_c12_branch, s0_n15_p_o16_branch, &
           s0p_pp, s0p_he3he3, s0p_he3he4, s0p_p_c12, s0p_p_c13, s0p_p_n14, &
           s0p_p_o16, s0pp_p_c12, s0pp_p_c13, s0pp_p_o16, s0p_be7_p, &
           s0pp_be7_p

! 10/14 MHP NEW PARAMETERS - REPLACING DTDIF,DTWIND, HPTTOL & ATIME VECTORS
! KC 2025-05-30 reordered common block elements
! common/newparam/: the "intuitively named" namelist parameters this
! routine remaps onto the legacy HPTTOL/ATIME/DTDIF/DTWIND/LPTIME
! storage above. Not referenced in any already-converted file.
      double precision :: flag_dx, flag_dw, flag_dz, time_core_min, &
           time_dl, time_dp, time_dr, time_dt, time_dw_global, time_dw_mix, &
           time_dx_core_frac, time_dx_core_tot, time_dx_shell, &
           time_dx_total, time_dy_core_frac, time_dy_core_tot, &
           time_dy_shell, time_dy_total, tol_czbase_fine_width, &
           tol_dl_max, tol_dm_max, tol_dm_min, tol_dp_core_max, &
           tol_dp_czbase_max, tol_dp_env_max, tol_dx_max, tol_dz_max, &
           time_max_dt_frac
      logical :: lstruct_time, lnewvars
      common/newparam/ flag_dx, flag_dw, flag_dz, &
           time_core_min, time_dl, time_dp, time_dr, time_dt, &
           time_dw_global, time_dw_mix, time_dx_core_frac, &
           time_dx_core_tot, time_dx_shell, time_dx_total, &
           time_dy_core_frac, time_dy_core_tot, time_dy_shell, &
           time_dy_total, tol_czbase_fine_width, tol_dl_max, tol_dm_max, &
           tol_dm_min, tol_dp_core_max, tol_dp_czbase_max, tol_dp_env_max, &
           tol_dx_max, tol_dz_max, time_max_dt_frac, lstruct_time, lnewvars

      save

      integer :: i

      double precision :: s0_pp_ref, s0_he3he3_ref, s0_he3he4_ref, &
           s0_p_c12_ref, s0_p_c13_ref
      data s0_pp_ref,s0_he3he3_ref,s0_he3he4_ref,s0_p_c12_ref,s0_p_c13_ref/ &
           4.07D-22,5.15D3,5.4D-1,1.45D0,5.5D0/
! NOTE: PEP IS THE PROPORTIONALITY CONSTANT RELATIVE TO PP
! NOTE: BE7+E- IS THE PROPORTIONALITY CONSTANT IN THE LINEAR TERM
! THE CODE USES T9, NOT T6, SO ANY EXPRESSION IN TERMS OF T/10^6 K
! NEEDS TO BE DIVIDED BY 1000^0.5 (FOR BOTH PEP AND BE7+E-)
      double precision :: s0_p_n14_ref, s0_p_o16_ref, s0_pep_ref, &
           s0_be7_electron_ref
      data s0_p_n14_ref,s0_p_o16_ref,s0_pep_ref,s0_be7_electron_ref/ &
           3.32D0,9.4E0,3.4848D-6,1.752D-10/
      double precision :: s0_be7_p_ref, s0_hep_ref, s0_n15_p_c12_branch_ref, &
           s0_n15_p_o16_branch_ref
      data s0_be7_p_ref,s0_hep_ref,s0_n15_p_c12_branch_ref, &
           s0_n15_p_o16_branch_ref/0.0243D0,8.0D-20,6.75D4,6.4D1/
! REFERENCE ADELBERGER ET AL. 1998 CROSS SECTIONS - USED IN HARD-CODED
! FIRST AND SECOND DERIVATIVE TERMS (EXPRESSED AS S'/S AND S''/S)
! UNITS ARE KeV b
!       DATA S0_1_1_A98,S0_3_3_A98,S0_3_4_A98,S0_1_12_A98,S0_1_13_A98/
!      *     4.00D-22,5.4D3,5.3D-1,1.34D0,7.6D0/  ! KC 2025-05-31
      double precision :: s0_pp_a98
      data s0_pp_a98/4.00D-22/
!       DATA S0_1_14_A98,S0_1_16_A98,S0_PEP_A98,S0_1_BE7E_A98/
!      *     3.5D0,9.4D0,3.5734D-6,1.77088D-10/
!       DATA S0_1_BE7P_A98,S0_HEP_A98,S0_1_15_C12ALP_A98,S0_1_15_O16_A98/
!      *     0.019D0,2.3D-20,6.75D4,6.4D1/
      double precision :: s0_be7_p_a98
      data s0_be7_p_a98/0.019D0/
! REFERENCE FIRST DERRIVATIVES OF CROSS-SECTIONS (ADELBERGER ET AL. 1998)
! UNITS ARE b
      double precision :: qs0e_pp_a98, qs0e_he3he3_a98, qs0e_he3he4_a98, &
           qs0e_p_c12_a98
      data qs0e_pp_a98,qs0e_he3he3_a98,qs0e_he3he4_a98,qs0e_p_c12_a98/ &
           4.48D-24,-4.1D0,-3.0D-4,2.6D-3/
      double precision :: qs0e_p_c13_a98, qs0e_p_n14_a98, qs0e_p_o16_a98, &
           qs0e_be7_p_a98
      data qs0e_p_c13_a98,qs0e_p_n14_a98,qs0e_p_o16_a98,qs0e_be7_p_a98/ &
           -7.8D-3,-1.28D-2,-2.4D-2,-1.33D-5/
! REFERENCE SECOND DERIVATIVES OF CROSS SECTIONS (ADELBERGER ET AL. 1998)
      double precision :: qqs0ee_p_c12_a98, qqs0ee_p_c13_a98, &
           qqs0ee_p_o16_a98, qqs0ee_be7_p_a98
      data qqs0ee_p_c12_a98,qqs0ee_p_c13_a98,qqs0ee_p_o16_a98, &
           qqs0ee_be7_p_a98/8.3D-5,7.3D-4,5.7D-5,7.22D-8/
! The code explicitly treats non-equilibrium C12,C13,N14,O16.
! N15 burning is assumed to immediately follow N14+p, with the
! branching ratio between C12+Alpha and O16+Gamma calculated using
! the relative S0 for the two.  To change the branching ratio, change
! S0_1_15_C12ALP and S0_1_15_O16.  O17 proton capture is assumed to
! immediately follow O16 proton capture, producing N14.
! He-burning reactions triple-alpha,C12+He4,C13+He4,N14+He4 are tracked
! in the code, using the Caughlin&Fowler(1988) cross-sections.  As resonances
! are important for these reactions, a simple scaling by S0 is not sufficient,
! and updates should instead use different functional forms (e.g. NACRE II).
! IF THE LNEWNUC FLAG IS SET TRUE THE SSTANDARD VECTOR IS NORMALIZED TO UNITY
! AND THE INDIVIDUAL CROSS-SECTION VALUES AND THEIR DERIVATIVES ARE USED INSTEAD.
      do i = 1,8
         qs0e_scale(i) = 1.0d0
         qqs0ee_scale(i) = 1.0d0
      end do
      c12_alpha_scale = 1.0d0
      o16_gamma_scale = 1.0d0
!      WRITE(*,911)(QS0E(I),I=1,8),(QQS0EE(I),I=1,8),
!     * (SSTANDARD(I),I=1,17),LNEWNUC
! 911  FORMAT(8E12.3/8E12.3/8E12.3/9E12.3/L5)
      if (use_new_nuclear_rates) then
         do i = 1,17
            cross_section_scale(i) = 1.0d0
         end do
! The energy generation subroutine engeb has a vector, SSTANDARD,
! which is the ratio of the desired S0 to that of Bahcall & Ulrich (1988).
! If this option is true, the namelist SSTANDARD vector is replaced
! by explicit cross-sections for individual reactions, and the implied
! Sstandard vector is derived here.  The first and second derivatives
! of the reactions are permitted to be separately modified, unlike in the
! original version of engeb.  In this case the inputs are the ratios
! S'/S0 and S''/S0, relative to the Solar Fusion I (Adelberger et al. 1998) ratios.
         cross_section_scale(1) = s0_pp/s0_pp_ref
         cross_section_scale(2) = s0_he3he3/s0_he3he3_ref
         cross_section_scale(3) = s0_he3he4/s0_he3he4_ref
         cross_section_scale(4) = s0_p_c12/s0_p_c12_ref
         cross_section_scale(5) = s0_p_c13/s0_p_c13_ref
         cross_section_scale(6) = s0_p_n14/s0_p_n14_ref
         cross_section_scale(7) = s0_p_o16/s0_p_o16_ref
!      SSTANDARD(8) = S0_4_13/S0_4_13_REF resonant; not fit by simple S0
!      SSTANDARD(9) = S0_4_16/S0_4_16_REF not used
!      SSTANDARD(10) = S0_4_12/S0_4_12_REF resonant; not fit by simple S0
!      SSTANDARD(11) = S0_4_14/S0_4_14_REF resonant; not fit by simple S0
!      SSTANDARD(12) = S0_3A/S0_3A_REF resonant; not fit by simple S0
!      SSTANDARD(13) = S0_12_12/S0_12_12_REF not used
         cross_section_scale(14) = s0_pep/s0_pep_ref
         cross_section_scale(15) = s0_be7_electron/s0_be7_electron_ref
         cross_section_scale(16) = s0_be7_p/s0_be7_p_ref
         cross_section_scale(17) = s0_hep/s0_hep_ref
! ABILITY TO CHANGE THE HARD-CODED BRANCHING RATIOS FOR THE OUTCOME OF N15+P
         c12_alpha_scale = s0_n15_p_c12_branch/s0_n15_p_c12_branch_ref
         o16_gamma_scale = s0_n15_p_o16_branch/s0_n15_p_o16_branch_ref
! NUCLEAR REACTION DERIVATIVES ARE USED IN THE VECTORS Q2-Q5 IN ENGEB AND NRATES.
! THE Q2 AND Q3 FACTORS ARE DEFINED PROPORTIONAL TO S'/S (REFERENCE: SOLAR FUSION I IN 1998)
! THE Q4 AND Q5 FACTORS ARE DEFINTED PROPORTIONAL TO S''/S (SAME REFERENCE)
! FIRST DERIVATIVE TERMS - TO BE USED TO MULTIPLY Q2 AND Q3 TERMS IN ENGEB
! MHP 4/25 FIXED TYPOS AND TRAPPED OUT POTENTIAL DIVIDE BY ZERO ERRORS
         if (s0_pp.gt.0.0d0) then
            qs0e_scale(1) = (s0p_pp/s0_pp)/(qs0e_pp_a98/s0_pp_a98)
         else
            qs0e_scale(1) = 0.0d0
         end if
         if (s0_he3he3.gt.0.0d0) then
            qs0e_scale(2) = (s0p_he3he3/s0_he3he3)/(qs0e_he3he3_a98/s0_pp_a98)
         else
            qs0e_scale(2) = 0.0d0
         end if
         if (s0_he3he4.gt.0.0d0) then
            qs0e_scale(3) = (s0p_he3he4/s0_he3he4)/(qs0e_he3he4_a98/s0_pp_a98)
         else
            qs0e_scale(3) = 0.0d0
         end if
         if (s0_p_c12.gt.0.0d0) then
            qs0e_scale(4) = (s0p_p_c12/s0_p_c12)/(qs0e_p_c12_a98/s0_pp_a98)
         else
            qs0e_scale(4) = 0.0d0
         end if
         if (s0_p_c13.gt.0.0d0) then
            qs0e_scale(5) = (s0p_p_c13/s0_p_c13)/(qs0e_p_c13_a98/s0_pp_a98)
         else
            qs0e_scale(5) = 0.0d0
         end if
         if (s0_p_n14.gt.0.0d0) then
            qs0e_scale(6) = (s0p_p_n14/s0_p_n14)/(qs0e_p_n14_a98/s0_pp_a98)
         else
            qs0e_scale(6) = 0.0d0
         end if
         if (s0_p_o16.gt.0.0d0) then
            qs0e_scale(7) = (s0p_p_o16/s0_p_o16)/(qs0e_p_o16_a98/s0_pp_a98)
         else
            qs0e_scale(7) = 0.0d0
         end if
!         QS0E(8)=(S0P_1_BE7P/S0_1_BE7P)/(QS0E_1_BE7P_A98/S0_1_BE7P__A98)
         if (s0_be7_p.gt.0.0d0) then
            qs0e_scale(8) = (s0p_be7_p/s0_be7_p)/ &
                 (qs0e_be7_p_a98/s0_be7_p_a98)
         else
            qs0e_scale(8) = 0.0d0
         end if
! SECOND DERIVATIVE TERMS - TO BE USED TO MULTIPLY Q4 AND Q5 TERMS IN ENGEB
!         QQS0EE(1) = (S0PP_1_1/S0_1_1)/(QQS0EE_1_1_A98/S0_1_1_A98) ZEROED OUT IN 2003 VERSION
!         QQS0EE(2) = (S0PP_3_3/S0_3_3)/(QQS0EE_3_3_A98/S0_3_3_A98) ZEROED OUT IN 2003 VERSION
!         QQS0EE(3) = (S0PP_3_4/S0_3_4)/(QQS0EE_3_4_A98/S0_3_4_A98) ZEROED OUT IN 2003 VERSION
         if (s0_p_c12.gt.0.0d0) then
            qqs0ee_scale(4) = (s0pp_p_c12/s0_p_c12)/ &
                 (qqs0ee_p_c12_a98/s0_pp_a98)
         else
            qqs0ee_scale(4) = 0.0d0
         end if
         if (s0_p_c13.gt.0.0d0) then
            qqs0ee_scale(5) = (s0pp_p_c13/s0_p_c13)/ &
                 (qqs0ee_p_c13_a98/s0_pp_a98)
         else
            qqs0ee_scale(5) = 0.0d0
         end if
!         QQS0EE(6) = (S0PP_1_14/S0_1_14)/(QQS0EE_1_14_A98/S0_1_14_A98) ZEROED OUT IN 2003 VERSION
         if (s0_p_o16.gt.0.0d0) then
            qqs0ee_scale(7) = (s0pp_p_o16/s0_p_o16)/ &
                 (qqs0ee_p_o16_a98/s0_pp_a98)
         else
            qqs0ee_scale(7) = 0.0d0
         end if
         if (s0_be7_p.gt.0.0d0) then
            qqs0ee_scale(8) = (s0pp_be7_p/s0_be7_p)/ &
                 (qqs0ee_be7_p_a98/s0_be7_p_a98)
         else
            qqs0ee_scale(8) = 0.0d0
         end if
      end if
! OPTION TO OVERWRITE SPATIAL AND TEMPORAL TOLERANCES WITH MORE INTUITIVE VARIABLE NAMES.
! NOTE: THIS OVERWRITES THE HPTTOL AND ATIME VECTORS.
      if (lnewvars) then
! MINIMUM LOG MASS SPACING BETWEEN SHELLS
         chi_grid_scale(1) = tol_dm_min
! MAXIMUM LOG MASS SPACING BETWEEN SHELLS
         chi_grid_scale(2) = tol_dm_max
! MAXIMUM L/LSURF SPACING BETWEEN SHELLS
         chi_grid_scale(9) = tol_dl_max
! MAXIMUM LOG P SPACING BETWEEN SHELLS, RADIATIVE INTERIOR
         chi_grid_scale(11) = tol_dp_core_max
! MAXIMUM LOG P SPACING BETWEEN SHELLS, CONVECTIVE ENVELOPE
         chi_grid_scale(8) = tol_dp_env_max
! MAXIMUM LOG P SPACING BETWEEN SHELLS, NEAR CZ BASE
         chi_grid_scale(10) = tol_dp_czbase_max
! WIDTH (IN LOG P) OF FINELY ZONED REGION NEAR THE CZ BASE
         chi_grid_scale(7) = tol_czbase_fine_width
! MAXIMUM X SPACING BETWEEN SHELLS (USUALLY DISABLED)
         chi_grid_scale(5) = tol_dx_max
! MAXIMUM Z SPACING BETWEEN SHELLS (USUALLY DISABLED)
         chi_grid_scale(6) = tol_dz_max
! X DIFFERENCE BETWEEN SHELLS REQUIRED TO FLAG (AVOID NUMERICAL DIFFUSION)
         chi_grid_scale(3) = flag_dx
! Z DIFFERENCE BETWEEN SHELLS REQUIRED TO FLAG (AVOID NUMERICAL DIFFUSION)
         chi_grid_scale(4) = flag_dz
! MAXIMUM LOG OMEGA SPACING BETWEEN SHELLS (ROTATING MODELS ONLY)
         chi_grid_scale(12) = flag_dw
! MINIMUM CENTRAL ABUNDANCE TO USE CORE BURNING DT CRITERION
         atime(1) = time_core_min
! MAX PERMITTED ABSOLUTE DELTA X BURNED IN CORE IN TIMESTEP
         atime(2) = time_dx_core_tot
! MAX PERMITTED FRACTION OF X BURNED IN CORE IN TIMESTEP
         atime(3) = time_dx_core_frac
! MAX PERMITTED SOLAR MASSES OF X BURNED IN TIMESTEP
         atime(6) = time_dx_total
! MAX PERMITTED FRACTION OF X BURNED IN SHELL IN TIMESTEP
         atime(7) = time_dx_shell
! MAX PERMITTED ABSOLUTE DELTA Y BURNED IN CORE IN TIMESTEP
         atime(4) = time_dy_core_tot
! MAX PERMITTED FRACTION OF Y BURNED IN CORE IN TIMESTEP
         atime(5) = time_dy_core_frac
! MAX PERMITTED FRACTION OF Y BURNED IN SHELL IN TIMESTEP
         atime(12) = time_dy_shell
! MAX PERMITTED SOLAR MASSES OF Y BURNED IN SHELL IN TIMESTEP
         atime(14) = time_dy_total
! MAX PERMITTED CHANGES IN T, P, R, L AT ANY GIVEN SHELL IN TIMESTEP
         atime(8) = time_dt
         atime(9) = time_dp
         atime(10) = time_dr
         atime(11) = time_dl
! MASTER FLAG THAT USES/DOES NOT USE ATIME(8)-(11)
! NOTE: THIS IS USED ONLY IN THE PRE-MS AND DISABLED FOR MS STARS AUTOMATICALLY.
         use_structure_dt_limits = lstruct_time
! MAX PERMITTED CHANGE IN OMEGA PER GLOBAL TIMESTEP
         max_domega_global = time_dw_global
! THE CODE CAN MIX IN A SERIES OF SMALLER TIMESTEPS; THIS IS
! THE MAX PERMITTED CHANGE IN OMEGA IN A SMALL, ROTATIONAL MIXING STEP
         dtdif = time_dw_mix
! THE TIMESTEP IS NOT PERMITTED TO CHANGE FROM ONE MODEL TO THE NEXT BY
! MORE THAN THIS SCALE FACTOR
!          ATIME(13) = TIME_MAX_DT_FAC
         atime(13) = time_max_dt_frac
      end if
      return
end subroutine remap
