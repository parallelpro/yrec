!----------------------------------------------------------------------
! engeb
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original engeb.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! 3/92 DBG Added new neutrino loss calculation routines
! NOVEMBER 6, 1990 (JNB)
! THE FOLLOWING SUBROUTINE CALCULATES THE NUCLEAR REACTION RATES WITH
!   SPECIAL ATTENTION TO DETAIL REQUIRED FOR CALCULATING SOLAR NEUTRINO
!   FLUXES. THE NEUTRINO FLUXES ARE EVALUATED ALSO IN THIS SUBROUTINE.
! ALL PHYSICAL QUANTITIES ARE CGS. THE REACTION RATES ARE IN
!  GM^{-1}S^{-1}.
!  MHP 10/97
!
!   On 10/13/97, JNB converted the nuclear masses from neutral nuclear
!   masses to bare nuclear masses by subtracting Z(I)*(m_e)*c^2 from
!   the neutral nuclear masses. This caused changes in a number of
!   places: in ANUC(I), in Q1(I)-Q7(I), in the calculation of the Be7electron
!   Be7proton rates, and in the calculation of the N15p branching ratio.
!
!   JNB made some purely cosmetic changes on 1/20/96. Revised all of the
!   input Q's on 9/23-25/97 to agree with submitted version of Solar
!   Fusion Workshop paper. The SStandard are fixed to agree with
!   the Workshop paper. JNB recalculated all of the EG(I) to determine
!   the best values for the energy generation for all the reactions,
!   taking account of my improved calculations of neutrino energy loss.
!   The calculations are documented in Vol. 19, 132-141, 1997 of my notes.
!   The neutral atom mass differences are taken from Table of Isotopes,
!   8th Ed, 1996 and the neutrino energy losses from Bahcall, Gallium
!   solar neutrino experiments, Phys Rev C, in press, 1997.
!
!  All numbers in this subroutine have been calculated by John Bahcall
!  so that they agree with the modern numbers in Neutrino
!  Astrophysics (1989) or much more recent results, as indicated in
!  comments lines.
!
!  PREVIOUS VERSIONS OF THE YALE CODE
!  USED THE NUCLEAR REACTION RATES PRIMARILY FROM FOWLER, CAUGHLAN, AND
!  ZIMMERMAN IN ANN. REV. ASTR. AND ASTROPHYS. (1975) OR, FOR THE BE7 AND
!  B8 RATES, THE VALUES FROM BAHCALL AND SEARS (1972), ALSO. ANN. REV. ASTR.
!  AND ASTROPHYS. SOME VALUES WERE FROM HARRIS ET AL (1983, ANN. REV. ASTR.
!  IN A FEW CASES, NUMERICAL ERRORS (NOT DUE TO REVISIONS IN THE NUCLEAR DATA)
!  HAVE BEEN CORRECTED.  IN A SMALL NUMBER OF OTHER PLACES, THE PREVIOUS
!  PROGRAM CONTAINED PHYSICALLY INCORRECT STATEMENTS, WHICH HAVE BEEN CORRECTED
!  IN THE PRESENT VERSION. IN THE COURSE OF THE CHECKING AND REVISIONS, I
!  HAVE ADDED MANY COMMENT STATEMENTS IN ORDER TO MAKE THE PROGRAM MORE
!  TRANSPARENT AND SIMPLER TO REVISE.  SOME STATEMENTS HAVE BEEN ADDED
!  TO SIMPLIFY THE PROCESS OF REVISING THE SUBROUTINE AS IMPROVED DATA
!  BECOME AVAILABLE.  FOR EXAMPLE, I HAVE INSERTED EXPLICITLY
!  THE STANDARD CROSS SECTION FACTORS IN A WAY THAT IS EASY TO REVISE;
!  SSTANDARD(I) = 1.0 WHEN THE ITH CROSS SECTION HAS THE STANDARD VALUE
!  MEASUREMENTS USED IN NEUTRINO ASTROPHYSICS.
!  MHP 10/97
!
!  The PREVIOUS standard values of SStandard(I) used were those that
!   were given in Table 1 of Bahcall and Pinsonneault (1992).  This work
!   was  published in the Reviews of Modern Physics, 64,885, 1992. The
!   current version of SStandard(I) refers to Bahcall and
!   Pinsonneault(1995) , Rev. Mod. Phys. 67, 781 (1995), Table 1, if
!   changes have occured.  Otherwise, they are the same as in the
!   Bahcall-Pinsonnealt (1992) Rev. Mod. Phys. article.
!
!  Weakscreening is a parameter set in this subroutine. To obtain
!  the Graboske et al. and Salpeter standard results, use:
!  weakscreening = 0.03.  To investigate the effect of always using
!  weak screening, use a value for weakscreening greater than unity,
!  e. g., 30. The transition region between weak and strong screening
!  is defined by Graboske et al. as $\Lambda_{12} = 0.1$, see page
!  ApJ 181, 465 (1973) [5/15/97].
!
!  The value of SStandard(17) for hep is taken from Carlson et al (1991)
!  Phys. Rev. C 44, 619. It corresponds to an S sub 0 =
!  1.3 E-20 keV-barns, a factor of 0.1625  smaller than indicated by the
!  older measurements and calculations used in Neutrino Astrophysics.
!
! THE NUCLEAR ENERGY RELEASE TO THE STAR FROM EACH REACTION
!  IS TAKEN FROM BAHCALL AND ULRICH, RMP 60, 297 (1988) AND TAKES ACCOUNT
!  ACCURATELY OF NEUTRINO ENERGY LOSS.
!
! I HAVE ADDED A NEW SECTION AT THE END OF THE SUBROUTINE THAT
!  CALCULATES THE SOLAR NEUTRINO FLUXES AT THE EARTH.  THESE FLUXES
!  ARE IN THE UNITS OF CM^-2 SEC^-2 PER GM. TO GET THE FLUX
!  FROM A SHELL, MULTIPLY BY THE MASS OF THE SHELL IN UNITS OF
!  GRAMS.  THE FLUXES ARE IN A COMMON BLOCK, FLUXES.  I ALSO CALCULATE
!  THE FICTIONAL NEUTRINO FLUXES ASSOCIATED WITH THE HE3 + HE3 AND WITH
!  THE HE3 + HE4 RECTIONS; THESE FICTIONAL FLUXES ARE USEFUL DIAGNOSTICS
!  OF THE SOLAR MODEL.
!
! *****************************************
! CHANGING NUCLEAR REACTION CROSS SECTIONS.
! *****************************************
! NUCLEAR CROSS SECTIONS CAN BE CHANGED SIMPLY BY INSERTING NEW NUMBERS
!  FOR THE DATA VALUES OF SSTANDARD(I), WHICH ARE THE RATIOS OF THE
!  DESIRED CROSS SECTION FACTORS TO THE VALUES GIVEN IN TABLES 3.2 AND
!  3.4 OF NEUTRINO ASTROPHYSICS (1989). IF THE VALUE GIVEN IN NEUTRINO
!  ASTROPHYSICS IS USED, THEN SSTANDARD(I) = 1.0 .  TO INCREASE THE CROSS
!  SECTION FOR REACTION K BY A FACTOR OF TWO COMPARED TO THE STANDARD
!  VALUE, SET SSTANDARD(K) = 2.0 .  TO DETERMINE WHICH VALUE OF I GOES
!  WITH WHICH REACTION, SEE THE SECTION JUST BELOW.
!  THE ENERGY DERIVATIVES ENTER IN A FORM IN WHICH THEY ARE DIVIDED BY
!  THE ABSOLUTE VALUES OF THE CROSS SECTIONS AT ZERO ENERGY.  THUS IF
!  THE SHAPE OF THE CROSS SECTION EXTRAPOLATION IS UNCHANGED AND ONLY
!  THE INTERCEPT OF S(E) AT ZERO ENERGY IS CHANGED, THEN NO CORRECTION
!  NEED BE MADE FOR THE DERIVATIVES.  THEY ARE AUTOMATICALLY SCALED
!  CORRECTLY.  THE EXACT WAY THAT THE DERIVATIVES ENTER THE
!  CALCULATIONS IS DESCRIBED IN THE SECTION LABELED ``DEFINING THE
!  Q(I)'' THAT IS PRESENTED BELOW.
! ************************************
! IDENTIFYING THE REACTIONS.
! ************************************
!  THE VALUE OF J DENOTES WHICH OF THE REACTIONS THE COEFFICIENTS
!  REFER TO:
!  J = 1, PP; J = 2, HE3+HE3; J = 3, HE3+ HE4; J =4, P + C12;  J = 5, P+C13;
!  J = 6. P + N14; J = 7, P + O16.
!  REACTIONS J = 8, 13 ARE NOT RELEVANT FOR THE SOLAR INTERIOR; THEY ARE
!   HOLDOVERS FROM THE EARLIER YALE CODE.
!  REACTION 14 IS PEP; REACTION 15 IS BE7 ELECTRON CAPTURE; REACTION 16 IS
!   BE7 PROTON CAPTURE; REACTION 17 IS THE HEP REACTION.
!  do not change sstandard(14) unless you want to change the ratio of
!  pep to pp.
!   REACTIONS 14-17 WERE NOT EXPLICITLY INCLUDED IN THE YALE
!    PREVIOUS VERSION OF THE CODE, BUT THEY ARE IMPORTANT FOR
!    SOLAR NEUTRINO CALCULATIONS.
!   THE BRANCHING OF THE N15 + P REACTIONS IS TREATED IN A SERIES OF
!    SEPARATE STATEMENTS FOLLOWING THE CALCULATION OF THE BE7 + P
!    REACTION. SEE THE DEFINITIONS OF F3 AND F4.  IF THE CROSS-SECTION
!    FACTORS OF THE N15 + P REACTIONS ARE REVISED, THEN THE NUMERICAL
!    COEFFICIENTS MUST BE CHANGED IN THE DEFINITION OF C12ALPHA AND
!    O16GAMMA.
! FOR Q1(I), ...,Q(5(I), I = 8 CORRESPONDS TO THE BE7 +  P REACTION.
!  THIS ASSIGNMENT FOR I = 8 IS ONLY VALID FOR THE LISTED Q'S AND NOT
!  FOR OTHER ARRAYS IN THE PROGRAM.
!  IU IS THE SHELL NUMBER.
!
! Dummy-argument renaming used throughout this modernized version:
!   EPP1 -> pp_chain_energy_gen        EPP2 -> he3he4_be7_electron_energy_gen
!   EPP3 -> he3he4_be7_proton_energy_gen   ECN -> cno_cycle_energy_gen
!   E3AL -> triple_alpha_energy_gen    PEP -> dlnepsilon_dlnrho
!   PET -> dlnepsilon_dlnt             SUM1 -> total_energy_gen_rate
!   DL -> log_density                  TL -> log_temperature
!   X -> hydrogen_fraction             Y -> helium_fraction
!   XHE3 -> he3_fraction               XC12/XC13 -> c12_fraction/c13_fraction
!   XN14 -> n14_fraction               XO16/XO18 -> o16_fraction/o18_fraction
!   XH2 -> deuterium_fraction          IU -> shell_index
!   HR1..HR13 -> reaction_rate_1..13 (yr^-1, amu^-1, output for kemcom)
!   HF1 -> n15_alpha_branch_fraction   HF2 -> be7_electron_capture_fraction
subroutine engeb(pp_chain_energy_gen, he3he4_be7_electron_energy_gen, &
     he3he4_be7_proton_energy_gen, cno_cycle_energy_gen, &
     triple_alpha_energy_gen, dlnepsilon_dlnrho, dlnepsilon_dlnt, &
     total_energy_gen_rate, log_density, &
     log_temperature, hydrogen_fraction, helium_fraction, he3_fraction, &
     c12_fraction, c13_fraction, n14_fraction, o16_fraction, &
     o18_fraction, deuterium_fraction, shell_index, reaction_rate_1, &
     reaction_rate_2, reaction_rate_3, reaction_rate_4, reaction_rate_5, &
     reaction_rate_6, reaction_rate_7, reaction_rate_8, reaction_rate_9, &
     reaction_rate_10, reaction_rate_11, reaction_rate_12, &
     reaction_rate_13, n15_alpha_branch_fraction, &
     be7_electron_capture_fraction)

      use luout_lib
      use const_lib
      use nuclear_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(out) :: pp_chain_energy_gen, &
           he3he4_be7_electron_energy_gen, he3he4_be7_proton_energy_gen, &
           cno_cycle_energy_gen, triple_alpha_energy_gen, &
           dlnepsilon_dlnrho, dlnepsilon_dlnt
! total_energy_gen_rate (originally SUM1) is intent(inout), not
! intent(out): in the log_temperature.le.tcut(1) early-return branch
! below (preserved verbatim from the original), it is never assigned,
! so the caller's incoming value is left untouched on that path -- a
! real property of the original F77 code, not a bug introduced here.
      double precision, intent(inout) :: total_energy_gen_rate
      double precision, intent(in) :: log_density, log_temperature, &
           hydrogen_fraction, helium_fraction, he3_fraction, c12_fraction, &
           c13_fraction, n14_fraction, o16_fraction, o18_fraction, &
           deuterium_fraction
      integer, intent(in) :: shell_index
      double precision, intent(out) :: reaction_rate_1(json), &
           reaction_rate_2(json), reaction_rate_3(json), &
           reaction_rate_4(json), reaction_rate_5(json), &
           reaction_rate_6(json), reaction_rate_7(json), &
           reaction_rate_8(json), reaction_rate_9(json), &
           reaction_rate_10(json), reaction_rate_11(json), &
           reaction_rate_12(json), reaction_rate_13(json), &
           n15_alpha_branch_fraction(json), &
           be7_electron_capture_fraction(json)

      double precision :: deuterium_burning_rate(json), &
           deuterium_burning_rate_start(json), accreted_mass_fraction
      integer :: jcz
      common/deuter/ deuterium_burning_rate, deuterium_burning_rate_start, &
           accreted_mass_fraction, jcz





! 7/91 COMMON BLOCK ADDED TO SKIP FLUX CALCULATIONS IF LSNU=F
! common/neweng/: only lsnu is used here. Naming matches mix.f90.
      integer :: niter4
      logical :: lnews, lsnu
      common/neweng/ niter4, lnews, lsnu

! common/neweps/: alpha_capture_energy/neutrino_loss_rate, both set
! here. Naming matches ytime.f90.
      double precision :: alpha_capture_energy, neutrino_loss_rate
      common/neweps/alpha_capture_energy, neutrino_loss_rate

! common/fluxes/: only neutrino_flux is used here. Naming matches
! wrtmonte.f90/wrtout.f90.
      double precision :: neutrino_flux(10), neutrino_flux_total(10), &
           cl37_snu_rate, ga71_snu_rate
      common/fluxes/ neutrino_flux, neutrino_flux_total, cl37_snu_rate, &
           ga71_snu_rate

! common/be7/: single-member block holding the Be7 mass fraction (not
! needed for the neutrino fluxes since Be7 is always in equilibrium,
! but of interest diagnostically). First appearance of this common
! block in the converted sources.
      double precision :: be7_mass_fraction
      common/be7/be7_mass_fraction

! 9/06 GN --- New neutrino loss common block
! KC 2025-05-30 reordered common block elements
!       COMMON/NULOSS/LNULOS1,DSNUDT,DSNUDD
! common/nuloss/: switch (use_itoh_neutrino_loss) selecting the Itoh
! 1996 neutrino-loss routines used below, and the log-derivatives of
! the resulting loss rate w.r.t. T/rho (neutrino_dlnq_dlnt/
! neutrino_dlnq_dlnd), all set/used here. First appearance of this
! common block in the converted sources.
      double precision :: neutrino_dlnq_dlnt, neutrino_dlnq_dlnd
      logical :: use_itoh_neutrino_loss
      common/nuloss/neutrino_dlnq_dlnt, neutrino_dlnq_dlnd, &
           use_itoh_neutrino_loss



! G Somers 3/17, ADDING NEW TAUCZ COMMON BLOCK
! common/ovrtrn/: only convective_turnover_timescale is used here.
! Naming matches mixcz.f90.
      logical :: use_new_turnover_timescale, calc_envelope_flag
      double precision :: convective_turnover_timescale, &
           convective_turnover_timescale_old, pphot, pphot0, fracstep
      common/ovrtrn/use_new_turnover_timescale, calc_envelope_flag, &
           convective_turnover_timescale, convective_turnover_timescale_old, &
           pphot, pphot0, fracstep

! JVS 10/11 Common block for He3+He3 luminosity
! common/grab/: he3_luminosity_placeholder/he3_total_placeholder are
! actively SET here (despite the "placeholder" name, inherited from
! wrtout.f90 where this block first appeared unused); the rate arrays
! are not used in this file. Naming matches wrtout.f90.
      double precision :: he3_luminosity_placeholder, he3_total_placeholder, &
           he3_he3_rate_placeholder(json), he3_he4_rate_placeholder(json)
      common/grab/ he3_luminosity_placeholder, he3_total_placeholder, &
           he3_he3_rate_placeholder, he3_he4_rate_placeholder

      double precision :: mass_fraction(13), reaction_rate(13), &
           dlnrate_dlnrho(13), dlnrate_dlnt(13), screening_factor(13), &
           dscreen_dlnrho(13), dscreen_dlnt(13), reaction_energy_gen(13), &
           charge_product(13), z53(13), z43(13), z23(13), &
           z86(13), q1(8), q2(8), q3(8), q4(8), q5(8), q6(7), q7(7), &
           q8(7), eg(50)
      double precision :: atomic_mass_amu(13), atomic_number(13)
      double precision :: v1(7), v2(7), v3(7)
      double precision :: years_per_sec_over_amu

      integer :: num_isotopes, nrxns

! ***MHP 5/91 STATEMENTS ADDED FOR EVOLVED STAR NEUTRINO LOSSES.
! NEUTRINO COEFFICIENTS

      data v1,v2,v3/ &
           6.002d19,2.084d20,1.872d21,9.383d-1,-4.141d-1,5.829d-2,5.5924, &
           4.886d10,7.580d10,6.023d10,6.290d-3,7.483d-3,3.061d-4,1.5654, &
           2.320d-7,8.449d-8,1.787d-8,2.581d-2,1.734d-2,6.990d-4,0.56457/
! ***************************
! ANUC ARE ATOMIC MASS UNITS.
! ***************************
!  On 10/13/97, JNB converted ANUC(I) from neutral atomic masses to
!  bare nuclear masses.
!
!  The scale is the mass of C12 divided by 12 or 931.49432 MeV,
!  which is 1.6605402 times 10^{-24} gm. The values are obtained
!  by dividing the mass excess (expressed in MeV) by 931.49432 MeV and
!  adding to this the atomic mass number, A.  The value for Be7, which
!  is used implicitly in this subroutine was, 7.016930, until 13/13/97.
!  We now use the bare nuclear mass of 7.014735 .
!
      data atomic_mass_amu/1.008665,1.007276,2.013553,3.015501,3.014933, &
           4.001506,11.996709,13.000064,13.999233,15.990526,17.994772, &
           19.986954,23.978458/, &
           atomic_number/0.,1.,1.,1.,2.,2.,6.,6.,7.,8.,8.,10.,12./, &
           num_isotopes/13/

! THE ISOTOPES ARE NEUTRON, H1, D, H3, HE3, HE4, C12, C13, N14, O16, O18,
!  NE20, MG24, RESPECTIVELY. ALL OF THESE NUMBERS WERE CHECKED.
! NELEM IS THE NUMBER OF ISOTOPES INCLUDED.
! **************************************************************************
! THE QUANTITIES Q1(J), Q2(J), ...,Q5(J) ARE THE TERMS IN EQUATION 3.14 OF
!  NEUTRINO ASTROPHYSICS AND IN EQUATION 53 OF FOWLER, CAUGHLAN, AND
!  ZIMMERMAN (1967) EQ. 53, IN BOTH CASES MULTIPLIED BY T SUB 9 ^(-2/3).
!  THE REACTIONS CORRESPONDING TO EACH J ARE LISTED ABOVE, UNDER:
!   IDENTIFYING THE RECTIONS.
!  FOR THIS SET OF PARAMETERS, AND ONLY FOR THIS SET OF PARAMETERS,
!  J = 8 CORRESPONDS TO THE BE7 + P REACTION.
! **************************************************************************
! GENERAL EXPRESSION FOR Q'S:
! T_9^(-2/3)[S_EFF/S(0)] =
!   [T_9^(-2/3) + Q1(I)T_9^(-1/3) + Q2(I) + Q3(I)T_9^(1/3) +
!    Q4(I)T_9^(2/3) + Q(5)T_9 ]
!
! BY COMPARISON WITH EQUATION 3.14 OF NEUTRINO ASTROPHYSICS, WE SEE THAT
!
!  Q1 = (5/(12*TAU))*T_9^(-1/3) .
!  Q2 = (S'/S)(E_0))*T_9^(-2/3) .
!  Q3 = (S'/S)(35/36)(K*10^9 K)
!  Q4 = (S''/2S)(E_0^2)(T_9^(-4/3)
!  Q5 = (89/72)(S''/S)(E_0)(KT)(T_9^(-5/3)
!
! EACH OF THE Q'S IS INDEPENDENT OF TEMPERATURE (T), AS CAN BE SEEN FROM
!  EQUATIONS 3.10 AND 3.11 .
!
! MHP 10/97
! All of the values of the Q1, ...., Q5 have been recalculated, using
!  where needed nuclear cross sections given in Tables 3.2 and 3.4 of
!  Neutrino Astrophysics.  They have been updated on 9/16/97 to reflect
!  the derivatives and cross section factors in Adelberger et al. 1998,
!  the paper on the Solar Fusion Workshop.
!
! ******************************************************************
! Q6 IS THE COEFFICIENT OF THE TEMPERATURE TERM IN THE DEFINITION OF
!  TAU, EQUATION 3.10 OF NEUTRINO ASTROPHYSICS.
!  TAU = Q6*(T SUB 9 TO THE (-1/3) POWER ).
! ******************************************************************
! SLIGHT CHANGES HAVE BEEN MADE IN THE PREVIOUS VALUES OF Q6 TO MAKE
!  THE DATA MORE ACCURATE.
! NOTE THAT Q6 IS NEGATIVE.
! ********************************************************************
! Q7 IS THE CONSTANT IN FRONT OF THE REACTION RATE. THE LOGARITHMS ARE
!  RELATED TO S FACTORS BY A LOGARITHM AND VARIOUS NUMERICAL FACTORS.
! ********************************************************************
! THE GENERAL RELATION IS:
!  Q7 = 70.62860 + ((LN(Z0*Z1/A))/3) -LN(A0*A1) - LN(S SUB 0)
!       -LN(1 + DELTA_01)
!  HERE S SUB 0 IS THE CROSS SECTION FACTOR IN UNITS OF KEV-BARNS.
!  THE NUMERICAL VALUES USED HERE ARE TAKEN FROM TABLES 3.4 AND 3.2
!   OF NUCLEAR ASTROPHYSICS.
!  THE QUANTITY DELTA_01 IS NON-ZERO (EQUAL TO UNITY) ONLY WHEN THE
!   TWO REACTING NUCLEI, 0 AND 1, ARE IDENTICAL.
! Q8 REFLECTS A TERM IN THE EXPONENTIAL THAT OCCURS IN THE RATES, THE
!  TERM BEING PROPORTIONAL TO E^( CONSTANT*T_9^2). SEE HARRIS, ET AL.
!  1983, ANNUAL REV. ASTRON. ASTROPHYS.21, 165 (1983) FOR THE MEANING
!  OF THIS OBSCURE TERM.
! THE VALUES OF Q1(I),...Q7(I) GIVEN BELOW WERE OBTAINED USING THE DATA
!  IN NEUTRINO ASTROPHYSICS WITH THE AID OF AUXILARY COMPUTER CODES
!  THAT GENERATED ACCURATE EVALUATIONS. THE NUMBERS GIVEN HERE ARE
!  IN MANY CASES COMPLETELY DIFFERENT FROM THE VALUES IN THE ORIGINAL
!  YALE CODE.
! NRXNS IS THE NUMBER OF REACTIONS BEING TRACKED.
! DBG 8/94 APPLIED MHP UPDATE TO NUCLEAR REACTIONS
! 10/13/97. Changed Q(I) so that are now calculated for the bare nuclear
! masses. The calculates were made cues.f . Last date on the calculations
! in my notes is 10/10/97.  Note Q6 is -tau(T_9 = 1); in cues.f, I calculate
! tau(T_6 = 1). The connection is:  Q6 = 0.1*tau(cues.f).
!
!
! BP00 values
!      DATA Q1/0.12317,.03392,.0325,.0304,.03035,.0273,.02494,.040572/,
!     1Q2/1.08749,-.273,-.2085,.7630,-0.4044,-1.60,-1.224, -0.2095/,
!     2Q3/.93833,-.0648,-.0474,.1626,-.08598,-.3064,-.2139,-0.0595/,
!     3Q4/0.,0.,0.,4.79,7.456,0.0,.69703,.16762/,
!     4Q5/0.,0.,0.,2.595,4.032,0.0,0.3097,.12114/,
!     5Q6/-3.3804,-12.2757,-12.826,-13.6899,-13.7173,-15.2281,-16.6925/,
!     6Q7/20.8964,76.6003,67.8036,69.130,70.3809,69.8517,70.8012/,
!     7Q8/0.,0.,0.,0.0,0.0,0.0,0./,
!     8NRXNS/13/
! Revised values 7/21/03
!      DATA Q1/0.12317,.03392,.0325,.0304,.03035,.0273,.02494,.040572/,
!     1Q2/1.08749,-.273,-.2085,.7630,-0.4044,-1.60,-1.224, 0.0/,
!     2Q3/.93833,-.0648,-.0474,.1626,-.08598,-.3064,-.2139,0.0/,
!     3Q4/0.,0.,0.,4.79,7.456,0.0,.69703,0.0/,
!     4Q5/0.,0.,0.,2.595,4.032,0.0,0.3097,0.0/,
!     5Q6/-3.3804,-12.2757,-12.826,-13.6899,-13.7173,-15.2281,-16.6925/,
!     6Q7/20.8964,76.6003,67.8036,69.130,70.3809,69.8517,70.8012/,
!     7Q8/0.,0.,0.,0.0,0.0,0.0,0./,
!     8NRXNS/13/
! MHP 9/14 RESTORED BE7+P DERIVATIVES THAT WERE ZEROED OUT
      data q1/0.12317,.03392,.0325,.0304,.03035,.0273,.02494,.040572/, &
           q2/1.08749,-.273,-.2085,.7630,-0.4044,-1.60,-1.224,-0.2095/, &
           q3/.93833,-.0648,-.0474,.1626,-.08598,-.3064,-.2139,-0.0595/, &
           q4/0.,0.,0.,4.79,7.456,0.0,.69703,0.16762/, &
           q5/0.,0.,0.,2.595,4.032,0.0,0.3097,0.12114/, &
           q6/-3.3804,-12.2757,-12.826,-13.6899,-13.7173,-15.2281,-16.6925/, &
           q7/20.8964,76.6003,67.8036,69.130,70.3809,69.8517,70.8012/, &
           q8/0.,0.,0.,0.0,0.0,0.0,0./, &
           nrxns/13/
! ***NOTE THAT SSTANDARD IS AN INPUT PARAMETER SET IN THE NAMELIST;
! PREVIOUS PUBLISHED SETS OF SSTANDARD ARE INDICATED BELOW.
! Changed slightly 3He-3He on 9/25/97 to take account of the S'.
!
!****************FD Feb 09 ************************************
! Original data from Neutrino Astrophysics 1989 Bahcall
! Table 3.2 and 3.4
! Table 3.2 pp
!         Q(MeV)  S_0(KeV barns)
! H-p     1.44    4.07E-022
! He3-He3 12.86   5.15E+003
! He3-He4 1.59    0.54
! 7Be-p   0.14    0.02
! He-p    19.8    8.00E-020
!
!
! Table 3.4 cno   Q(MeV)  S_0(MeV barns)
!
! C12-p   1.94    1.45E-003
! C13-p   7.55    5.50E-003
! N14-p   7.3     3.32E-003
! N15-p O16       12.13   6.40E-002
! N15-p C12       4.97    78
! O16-p   0.6     9.40E-003
!
!****************************************************************************
! FD-MP Fev 2009 values changed to included the results presented at Gran Sasso
! by A. Formicola
! Cross sections are:
!  s11      pp      = 3.92 +- 0.08 D-25 MeVb      Physun talk Formicola
!  s13      Hep     = 8.6  +- 1.3  D-20 KeVb
!  s33      He3-He3 = 5.32 +- 0.08 MeVb
!  s34      He3_He4 = 0.568 +- 0.014 KeVb         mean LUNA 2007 - Brown 2007 - Singh 2004
!  s17      Be7+p   = 22.1 +-0.6 +-0.6 eVb        A.Junghans PRC70(2004)045501
!  s114     N14-p   = 1.57 +- 0.13 KeVb           Marta et al. Luna coll. PRC 78 (2008) 022802)
!
! This cross section factor gives the following SSTANDARD
!
!  SSTANDARD=0.9631,1.0330,1.0519,0.9241,1.3818,0.47290,1.0,
! 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0108,0.84770,1.0750
!**********************************************************************
! BP03 values, changed 7/21/03
!     SStandard/0.9681,1.0485,0.9815,0.9241,1.3818,1.0542,1.0,
!    $          1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0108,0.8807,1.075/
!  Previously (6/16/97) used S at Gamow Peak. Agrees with Workshop paper.
!
!     SStandard/0.9828,1.0485,0.9815,0.9241,1.3818,1.0542,1.0,
!    $          1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0108,0.7819,0.2875/
!  This is the revised version created from the Seattle Workshop paper
!  on Solar Fusion Reactions. Corrections inserted 9/16/1997. See pg.
!  132, Volume 19, of notes.
!     9  SStandard/0.9828,1.0291,0.9815,0.9241,1.3818,1.0542,1.0,
!     $ 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.018,0.7819,0.2875/
! VALUES BASED ON PARKER REVIEW ARTICLE.
!    9  SSTANDARD/1.0049,0.9709,0.9870,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
!    $  1.0,1.0,1.0,1.0,1.0,1.913/
! BP 92 VALUES
!     9  SSTANDARD/0.9828,0.9709,0.9870,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
!     $  1.0,1.0,1.0,1.0,0.9218,0.1625/
! BP 95 VALUES
!     9  SSTANDARD/0.9558,0.9690,0.9712,1.0,1.0,0.992,1.0,1.0,1.0,1.0,
!     $  1.0,1.0,1.0,1.0,1.0,0.92088,0.1625/
! LATEST VALUES, INCLUDING NEW PP RATE
!    9  SSTANDARD/0.9558,0.9709,0.9870,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
!    $  1.0,1.0,1.0,1.0,0.9218,0.1625/
! Dar & Shaviv Values
!    9  SSTANDARD/0.9828,0.9709,0.8333,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
!    $  1.0,1.0,1.0,1.0,0.6996,0.1625/
! Schramm & Shi Values
!    9  SSTANDARD/0.9828,1.0874,0.9870,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
!    $  1.0,1.0,1.0,1.0,0.8313,0.1625/
! *********************************************************************
! THE VALUES OF SSTANDARD(I) ARE TO BE CHANGED FROM UNITY IF THE CROSS
!  SECTION FACTORS ARE NOT THE ONES GIVEN IN NEUTRINO ASTROPHYSICS,
!  TABLE 3.2 AND TABLE 3.4 .
! THE CURRENT VALUES ARE CHANGED TO BE THE PREFERRED VALUES LISTED IN
!  THE LAST COLUMN OF TABLE 1 OF BAHCALL AND PINSONNEAULT (1991),
!   6/14/91.
! *********************************************************************
! ZPRD IS USED IN SCREENING CALCULATIONS. IT IS THE PRODUCT OF THE
!  CHARGES OF THE INTERACTING IONS. ZPRD WAS CHECKED. Z86 IS USED
!  IN CALCULATING INTERMEDIATE SCREENING AND IS DEFINED BY GRABOSKE ET
!  AL, AP. J. 181, PAGE 465 (1973), IN TABLE 4.  Z86 WAS CHECKED AND
!  SOME NUMERICAL VALUES WERE MADE SLIGHTLY MORE ACCURATE.  Z53, Z43,
!  AND Z23 ARE ALSO DEFINED IN TABLE 4 (SEE ABOVE).  SINCE THEY ARE
!  ONLY USED IN STRONG SCREENING, THE VALUES OF Z53, Z43, AND Z23
!  WERE NOT CHECKED.
      data charge_product/1.,4.,4.,6.,6.,7.,8.,12.,16.,12.,14.,12.,36./, &
           z53/1.175, &
           3.73,3.73,4.804,4.804,5.385,5.941,9.014,11.24,9.014,10.15, &
           9.104,23.28/,z43/0.52,1.31,1.31,1.488,1.488,1.61,1.721,2.577, &
           3.025,2.577,2.81,2.577,5.668/,z23/-0.413,-0.655,-0.655,-0.643, &
           -0.643,-0.659,-0.673,-0.889,-0.946,-0.889,-0.92,-0.889,-1.36/, &
           z86/1.630,5.917,5.917,8.302,8.302,9.520,10.716,16.192,20.978, &
           16.192,18.606,16.192,45.6635/, &
           years_per_sec_over_amu/5.240358d-8/
      save

      double precision :: term, ion_mean_weight_inverse, &
           electron_mean_weight_inverse, xtr, zeta_sum, &
           electron_number_density_na, dd, density, t9, t9_p13, t9_p23, &
           t9_m13, t9_m23, t9_m1, t9_m2, t9_m12, t9_m32
      double precision :: dgdeut, qrtdeut, en, rdeut, qdeut, zz, &
           tfacdeut, tfacdeut2, rdeutmax, rdeut2
      double precision :: pfmc2, efmkt, fprf, degd
      double precision :: xxl, xxl6, xxl8, zcurl, zbar, z58, z28, z33, tm1
      double precision :: uwk, uint, ustr
      double precision :: r1, r2, a1, a2, a3, a4, a5, dr1, da1
      double precision :: be7electron, be7proton, temp3, qrbe7, &
           zprdbe7p, z86be7p, utotbe7p, camube7
      double precision :: f1, f2, f3, f4, o16gamma, c12alpha
      double precision :: convert, egdeut
      double precision :: sum2, sum3
      double precision :: fourpiau2, q6hep, zprdhe3p, z86he3p, utothe3p
      double precision :: flux_value
      double precision :: carbon_fraction_total, oxygen_fraction_total, &
           neutrino_temp, neutrino_density, neutrino_loss_snu
      double precision :: el, eli, ez, ez3, emue, ex1, ex2, ex3, el2
      double precision :: polx10, polx11, polx12, polx21, polx22, polx31, &
           polx32, qedn, qetn, qetnx, qednx
      integer :: i, k, nz

! ZERO OUT THE ENERGY YIELDS FROM NEUTRINOS(ENU) AND ALPHA CAPTURE
! REACTIONS (EALPCA).
      neutrino_loss_rate = 0.0d0
      alpha_capture_energy = 0.0d0
! DEFINE NEXT THE FRACTIONAL ABUNDANCES BY MASS OF THE IMPORTANT
!  ISOTOPES.
! X, Y, Z, XHE3,..., XBE9 ARE THE MASS FRACTIONS OF THE ISOTOPES.
!  THE ABUNDANCES OF NEUTRONS, H2, H3, NE20,AND MG24, WHICH ARE,
!  RESPECTIVELY, XFRAC(I) FOR I = 1,3,4,12,13, ARE NO LONGER USED.
      mass_fraction(1) = 0.0
      mass_fraction(2) = hydrogen_fraction
! MHP 5/02 ADDED DEUTERIUM
!      XFRAC(3) = 0.0
      mass_fraction(3) = deuterium_fraction
      mass_fraction(4) = 0.0
      mass_fraction(5) = he3_fraction
      mass_fraction(6) = helium_fraction
      mass_fraction(7) = c12_fraction
      mass_fraction(8) = c13_fraction
      mass_fraction(9) = n14_fraction
      mass_fraction(10) = o16_fraction
      mass_fraction(11) = o18_fraction
      mass_fraction(12) = 0.0
      mass_fraction(13) = 0.0
! *******************************************************************
! BEGIN CALCULATION OF SCREENING CORRECTION.
! *******************************************************************
!  THE BASIC REFERENCES ARE SALPETER, AUSTRALIAN JOURNAL OF PHYSICS,
!  VOL. 7, 373 (1954). THE FORMULA FOR WEAK SCREENING THAT IS BEING
!  PROGRAMMED IS EQUATION (25) OF THIS PAPER. THE OTHER IMPORTANT
!  REFERENCES ARE: DEWITT, GRABOSKE, AND COOPER, AP. J. 181, 439 (1973)
!  AND GRABOSKE ET AL., AP. J. 181, 457 (1973). THE VALUES OF EMU AND
!  ZET ARE ESSENTIAL FOR COMPUTING WEAK SCREENING; THE VALUE OF AMU IS
!  USED IN AN NON-ESSENTIAL WAY IN THIS COMPUTATION. XTR IS USED IN
!  COMPUTING INTERMEDIATE SCREENING.
! AMU IS ONE OVER THE MEAN MOLECULAR WEIGHT OF THE IONS, MU SUB I .
! EMU IS ONE OVER THE ELECTRON MEAN MOLECULAR WEIGHT, MU SUB E.
!  EMU IS USED HERE AS THE NAME FOR THE SECOND PART OF THE ZETA FUNCTION
!  IN THE SALPETER EXPRESSION FOR WEAK SCREENING.
! XTR IS USED LATER IN THE INTERMEDIATE SCREENING CALCULATION.  THE AVERAGE
!  OF THE QUANTITY Z**(3B -1) IS EQUAL TO XTR/AMU.
! ZET IS THE FIRST PART OF THE SALPETER SCREENING ZETA VARIABLE.
! MU = SUM OVER I OF [X(I)/A(I)].
! MU SUB E = SUM OVER I [ Z(I)*X(I)/A(I)].
      ion_mean_weight_inverse = 0.
      electron_mean_weight_inverse = 0.
      xtr = 0.
      zeta_sum = 0.
      do 10 i = 1,num_isotopes
         term = mass_fraction(i)/atomic_mass_amu(i)
         ion_mean_weight_inverse = ion_mean_weight_inverse+term
         electron_mean_weight_inverse = electron_mean_weight_inverse+ &
              term*atomic_number(i)
         xtr = xtr+term*atomic_number(i)**1.58
         zeta_sum = zeta_sum+term*atomic_number(i)**2
   10 continue
! DL AND DT ARE THE THE LOG10 OF THE DENSITY AND TEMPERATURE.
!  THE UNIT OF TEMPERATURE IS 10^9 K AND THE UNIT OF DENSITY IS
!  GM PER CM^3 .
! PDT AND PDP ARE THE DERIVATIVES OF THE DENSITY WITH RESPECT TO
!  TEMPERATURE AND DENSITY.
! DD = LOG RHO TO THE BASE 10.
! CLN = LN10.  CLN IS CONVERSION BETWEEN LOG10 AND LN.
! CONVERT DENSITY TO UNLOGGED FORM.
! RWE = RHO/(MU SUB E), I. E., THE NUMBER OF ELECTRONS DIVIDED BY
!  AVOGADRO'S NUMBER.
      electron_number_density_na = ( exp(ln10*log_density) )* &
           electron_mean_weight_inverse
! THE EXPRESSION FOR RWE WAS INCORRECT IN THE ORIGINAL YALE SUBROUTINE.
!  THE ORIGINAL VERSION HAD ( EXP(CLN*DL) ) DIVIDED BY EMU INSTEAD OF
!  MULTIPLIED BY EMU.  RWE IS USED LATER IN COMPUTING THE SCREENING
!  CORRECTION.
      dd = log_density
! SET RATES EQUAL TO ZERO FOR THE LOG_10(T) < 6.0.
! REPLACED FIXED 1 MILLION K THRESHOLD WITH TCUT(1).
!      IF(TL.LE.6.0) THEN
      if (log_temperature.le.tcut(1)) then
! MHP 5/02 DEUTERIUM BURNING
         dgdeut = 0.0d0
         qrtdeut = 0.0d0
         en = -20.
         dlnepsilon_dlnrho = 0.
         dlnepsilon_dlnt = 0.
         do 20 i = 1,nrxns
            eg(i) = 0.
            reaction_rate(i) = 0.
            reaction_energy_gen(i) = 0.
   20    continue
         go to 200
      end if
! T9P13 IS THE TEMPERATURE IN UNITS OF 10^9 DEGREES K TO THE PLUS 1/3
!  POWER.  MINUS IS DENOTED BY M.  HERE T9 IS THE TEMPERATURE IN UNITS
!  OF 10^9 K, CONVERTED FROM THE LOG_10 (T) AND RHO IS THE DENSITY IN
!  CGS UNITS.
      density=exp(ln10*dd)
      t9 = exp(ln10*(log_temperature - 9.0d0))
      t9_p13 = t9**cc13
      t9_p23 = t9_p13**2
      t9_m13=1./t9_p13
      t9_m23=t9_m13**2
      t9_m1=1./t9
      t9_m2=t9_m1**2
      t9_m12=1./dsqrt(t9)
      t9_m32=t9_m1*t9_m12
! MHP 5/02 ADD DEUTERIUM BURNING TERM TO THE CODE
! IF DEUTERIUM IS ABOVE A MINIMUM THRESHOLD VALUE.
! RDEUT IS THE RATE (EXCLUDING FACTORS OF THE
! ABUNDANCES) AND QRTDEUT IS THE DERIVATIVE W/R/T T.
! NOTE THAT SCREENING IS EXCLUDED - REASONABLE GIVEN
! THE LOW TEMPERATURES INVOLVED.
      if (deuterium_fraction.le.1.0d-11) then
         rdeut = 0.0d0
         qrtdeut = 0.0d0
      else
! ENERGY YIELD FOR DEUTERIUM BURNING
        qdeut = 5.494d0
        zz = -3.72d0*t9_m13
        tfacdeut = 1.0d0+0.112d0*t9_p13+3.38d0*t9_p23+2.65d0*t9
! FACTOR OF 6.023D23/ REFLECTS AVAGADROS NUMBER DIVIDED BY THE
! MASS OF THE DEUTERON IN AMU
        rdeut = density*2.240d3*t9_p23*exp(zz)*tfacdeut*6.023d23/ &
             atomic_mass_amu(3)
        tfacdeut2 = 0.112d0*t9_p13+6.76d0*t9_p23+7.95d0*t9
        qrtdeut = cc13*((tfacdeut2/tfacdeut) -2.0d0 - zz)
! NOW LIMIT DEUTERIUM BURNING IN A SURFACE CZ TO BE ON A TIME SCALE
! NO SHORTER THAN THE CONVECTIVE OVERTURN TIMESCALE.
        if (shell_index.ge.jcz .and. convective_turnover_timescale.gt.1.0d0) then
           rdeutmax = 6.023d23/atomic_mass_amu(3)/convective_turnover_timescale
           rdeut2 = rdeut*hydrogen_fraction
           if (rdeut2.gt.rdeutmax) then
! JVS 0712 Commented out write command
!              WRITE(*,*)RDEUT2,RDEUTMAX
              if (hydrogen_fraction.gt.1.0d-6) then
                 rdeut = rdeutmax/hydrogen_fraction
              end if
           end if
         end if
      end if
! ***********************************
! F PRIME/F
! ***********************************
! THE NEXT PIECE OF CODE COMPUTES FIRST THE FERMI ENERGY DIVIDED BY KT, WHERE
!  PFMC2 IS THE FERMI MOMENTUM DIVIDED BY MASS OF ELECTRON TIMES C, ALL
!  SQUARED AND EFMKT IS THE FERMI ENERGY DIVIDED BY KT.  THE QUANTITY
!  FPRF IS THE RATIO OF F PRIME TO F IN SALPETER'S SCREENING CORRECTION.
!  THE VALUE OF FPRF IS DETERMINED IN THE INTERMEDIATE CASE BY AN
!  INTERPOLATION FORMULA DEPENDING UPON THE DEGREE OF DEGENERACY, AS
!  MEASURED BY DEGD = LOG10(E_F/KT). THE NUMERICAL VALUES FOR THE FIT
!  WERE TAKEN FROM SALPETER'S ORIGINAL PAPER, FIGURE 1. THE ONLY CHANGES
!  IN THIS PART OF THE SUBROUTINE WERE THE CORRECTION OF THE ERROR IN
!  THE DEFINITION OF RWE (SEE ABOVE) AND REFINEMENTS OF THE COEFFICIENTS
!  IN THE EXPRESSIONS FOR PFMC2 AND EFMKT.
      pfmc2=1.017677e-4*electron_number_density_na**0.6666667
      efmkt=5.92986*t9_m1*(dsqrt(1.+pfmc2)-1.)
      if (efmkt.le.1.e-2) then
         fprf=1.0
      else
         degd=dlog10(efmkt)
         if (degd.ge.1.5) then
            fprf=0.0
         else
            fprf=0.75793-0.54621*degd-0.30964*degd**2+0.12535*degd**3+ &
            0.1203*degd**4-0.012857*degd**5-0.014768*degd**6
         end if
      end if
! END OF CALCULATION OF FERMI ENERGY DIVIDED BY KT AND THE INTERPOLATION
!  FORMULA FOR (F PRIME/F), WHICH APPEARS IN SALPETER'S SCREENING
!  FORMULA.
! NOW WE GET TO THE COMPUTATION OF WEAK SCREENING (SEE ALSO PAGE 61 OF
!  NEUTRINO ASTROPHYSICS, WHICH GIVES ONLY A SIMPLIFIED FORMULA) AND
!  ALSO THE MORE COMPLICATED INTERMEDIATE AND STRONG SCREENING CASES
!  (SEE REFERENCES TO AP. J. 181 ABOVE FOR INTERMEDIATE
!  AND STRONG SCREENING), ESPECIALLY PAGE 465.
! XXL IS USED IN ALL THE SCREENING FORMULAE.  THIS QUANTITY IS THE
!  IS THE FUNCTION CALLED LAMBDA SUB ZERO BY GRABOSKE ET AL. THE
!  QUANTITY XXL**0.86 = XXL8 IS USED IN CALCULATING INTERMEDIATE
!  SCREENING.
! ZCURL IS THE QUANTITY DEFINED BY EQUATION (4) OF DEWITT ET AL; IT
!  IS THEIR Z WITH A CURLY SYMBOL ON ITS TOP. ZCURL IS USED IN WEAK
!  AND IN INTERMEDIATE SCREENING AND WAS FIRST DEFINED BY SALPETER.
!  ZCURL IS THE SAME AS SALPETER'S ZETA EXCEPT FOR THE FACTOR OF 1/AMU.
! (Z SUB 1 TIMES Z SUB 2)*XXL*ZCURL GIVES THE WEAK SCREENING FACTOR,
!  THE SAME AS SALPETER OR AS EQUATION (19) OF GRABOSKE ET AL.
! Z BAR IS THE SAME AS Z BAR OF DE WITT ET AL AND IS THE AVERAGE CHARGE
!  OF THE IONS.  IT IS EQUAL TO EMU/AMU.  THE QUANTITY (Z BAR)**0.28 =
!  Z28 OCCURS IN THE COMPUTATION OF INTERMEDIATE SCREENING.
! XXL6 IS USED FOR COMPUTING STRONG STRONG SCREENING.
! THE NOTATION USED HERE IS EXPLAINED IN LARGE PART BY
!  EQUATION (4) OF DEWITT ET AL., AP. J. 181, PAGE 439.
! THE FINAL EXPRESSION FOR WEAK SCREENING IS EXACTLY EQUAL TO SALPETER'S
!  FORMULA, WHICH INCLUDES A DEGENERACY CORRECTION. THE MORE GENERAL
!  EXPRESSIONS ARE GIVEN IN TABLE 4 AND EQUATION (19) OF GRABOSKE ET AL.
      xxl=5.9426e-6*t9_m32*dsqrt(density*ion_mean_weight_inverse)
      xxl6=xxl**0.666667
      xxl8=xxl**0.86
      zcurl=dsqrt((zeta_sum+fprf*electron_mean_weight_inverse)/ &
           ion_mean_weight_inverse)
      zbar=electron_mean_weight_inverse/ion_mean_weight_inverse
      z58=zcurl**0.58
      z28=zbar**0.28
      z33=zbar**cc13
      tm1=xxl*zcurl
! COMPUTE SCREENING FOR EACH OF THE REACTIONS.
      do 30 i=1,nrxns
         uwk=tm1*charge_product(i)
         if (uwk.le.weak_screening_threshold) then
! WEAKSCREENING IS A NUMERICAL PARAMETER PASSED IN THE FLUX COMMON
!  BLOCK. TO OBTAIN THE GRABOSKE ET AL. AND SALPETER STANDARD RESULTS,
!  USE: WEAKSCREENING = 0.03.  FOR THE STANDARD SOLAR MODEL, THIS IS THE
!  VALUE THAT SHOULD BE ADOPTED. TO INVESTIGATE THE EFFECT OF ALWAYS USING
!  WEAK SCREENING, USE A LARGE VALUE FOR WEAKSCREENING, E. G., 30.  AS
!  LONG AS WEAKSCREENING IS ASSUMED TO BE BIGGER THAN ONE, THE PROGRAM
!  WILL ALWAYS CALCULATE FOR THE SUN WITH THE WEAK SCREENING
!  APPROXIMATION.
! UTOT IS THE FINAL SCREENING CORRECTION WHICH APPEARS IN THE
!  RATE EXPRESSION AS: EXP(UTOT) .
! DSCR IS THE LOGARITHMIC DERIVATIVE WITH RESPECT TO DENSITY OF THE
!  SCREENING CORRECTION, D LOG (E^(U_TOT)) /D LOG RHO .
! DSCT IS THE LOGARITHMIC DERIVATIVE OF THE SCREENING WITH RESPECT TO T,
! D LOG (E^(U_TOT)) /D LOG T .  THE LIMITING FORMULAE GIVEN IN THE FIRST
! OPTION ARE OBVIOUS SINCE IN THE WEAK LIMIT SALPETER'S FORMULA SHOWS
! THAT U IS PROPORTIONAL TO THE SQUARE ROOT OF RHO TIMES T^(-3/2).
            screening_factor(i)=uwk
            dscreen_dlnrho(i)=0.5*uwk
            dscreen_dlnt(i)=-1.5*uwk
         else
            uint=0.38*xxl8*xtr*z86(i)/(ion_mean_weight_inverse*z58*z28)
            if (uwk.le.2.) then
               screening_factor(i)=uint
               dscreen_dlnrho(i)=0.43*uint
               dscreen_dlnt(i)=-1.29*uint
            else
               ustr=0.624*z33*xxl6*(z53(i)+0.316*z33*z43(i)+0.737* &
               z23(i)/(zbar*xxl6))
               if (ustr.lt.uint.or.uwk.ge.5.) then
                  screening_factor(i)=ustr
                  dscreen_dlnrho(i)=0.208*z33*(z53(i)+0.316*z33*z43(i))*xxl6
                  dscreen_dlnt(i)=-3.*dscreen_dlnrho(i)
               else
                  screening_factor(i)=uint
                  dscreen_dlnrho(i)=0.43*uint
                  dscreen_dlnt(i)=-1.29*uint
               end if
            end if
         end if
   30 continue
! ****************************************************************
! END OF SCREENING CALCULATION. WEAK AND INTERMEDIATE SCREENING FORMS
!  ARE GIVEN CORRECTLY.  STRONG SCREENING WAS NOT CHECKED BECAUSE IT IS
!  NOT RELEVANT FOR THE SUN.
! ****************************************************************
      nz=1
      if (hydrogen_fraction.eq.0.0) then
         f1=0.
         f2=0.
         f3=0.
         f4=0.
         goto 50
      end if
      nz=8
! **************************************************************
!  CALCULATE REACTION RATES FOR THE THREE PRINCIPAL RECTIONS OF
!   THE PP CHAIN: PP, HE3+HE3, HE3 +HE4, AND THE FOUR PROTON
!   BURNING REACTIONS ON C12: C13, N14, AND O16.
! **************************************************************
! R1 IS (T SUB 9)^(-3/2) TIMES (S SUB EFF)/(S SUB 0). THE CORRECT
!  EXPRESSION FOR (S SUB EFF)/(S SUB 0) IS GIVEN IN EQUATION 3.14 OF
!  NEUTRINO ASTROPHYSICS.  THE NUMERICAL FORM THAT IS USED IS EQUATION
!  52 OF FOWLER, CAUHLAN, AND ZIMMERMAN, VOL. 5, 1967.
! RATE(I) IS THE RATE OF THE DIFFERENT REACTIONS PER SECOND PER GRAM,
!  EXCEPT THAT THE MASS FRACTIONS ARE OMITTED AT THIS POINT AND PUT IN
!  LATER.
! DRATT IS LOGARITHMIC DERIVATIVE OF RATE WITH RESPECT TO TEMPERATURE,
!  D LOG RATE DIVIDED BY D LOG T, LOG TO BASE 10.
! DRATRO IS THE LOGARITHMIC DERIVATIVE OF THE RATE WITH RESPECT TO
!  DENSITY, D LOG RATE/D LOG RHO, LOG TO BASE 10.
! THE PREVIOUS YALE VERSION HAD THE RATE FOR THE O16 + P REACTION
!  MULTIPLIED BY T9**(-1/7). THIS FACTOR IS INCORRECT AND HAS BEEN
!  REMOVED; IT APPEARED BEFORE AS AN IF STATEMENT REFERRING ONLY TO
!  RATE(7).
      do 40 i=1,7
!         R1=T9M23+Q1(I)*T9M13+Q2(I)+Q3(I)*T9P13+Q4(I)*T9P23+Q5(I)*T9
! MHP 8/14 RATES CORRECTED TO PERMIT USER MODIFICATION OF REACTION
! RATE DERIVATIVES
         r1=t9_m23+q1(i)*t9_m13+qs0e_scale(i)*(q2(i)+q3(i)*t9_p13)+ &
              qqs0ee_scale(i)*(q4(i)*t9_p23+q5(i)*t9)
         reaction_rate(i)=density*r1*exp(q6(i)*t9_m13+q7(i)+(q8(i)*t9)**2+ &
              screening_factor(i))
         reaction_rate(i) = reaction_rate(i)*cross_section_scale(i)
         if (reaction_rate(i).lt.1.e-30) then
            reaction_rate(i)=0.
            dlnrate_dlnt(i)=0.
         else
            dlnrate_dlnrho(i)=1.+dscreen_dlnrho(i)
            dlnrate_dlnt(i)=dscreen_dlnt(i)-(q6(i)*t9_m13+(2.*t9_m23+ &
            q1(i)*t9_m13-q3(i)*t9_p13-2.*q4(i)*t9_p23-3.*q5(i)*t9)/r1)/3.+ &
            2.*(q8(i)*t9)**2
         end if
   40 continue
! ***************************************************************
! END OF CALCULATION OF REACTION RATES FOR FIRST 7 REACTIONS.
! ***************************************************************
! ***********************************************
! BE7 BURNING
! ***********************************************
! CALCULATE THE BURNING OF BE7 BY PROTONS (WHICH PRODUCES THE MOST
!  EXPERIMENTALLY ACCESSIBLE SOLAR NEUTRINOS) AND THE BURNING OF BE7
!  BY ELECTRON CAPTURE (THE DOMINANT PROCESS).
!  THE ELECTRON CAPTURE RATE IN SEC^(-1) IS GIVEN BY EQUATION (3.18) OF
!  NEUTRINO ASTROPHYSICS.  CAN OMIT THE FACTOR OF RHO (IN CGS UNITS)
!  WHICH ALSO APPEARS IN THE BE7PROTON RATE. THE BE7 PROTON CAPTURE
!  RATE IS GIVEN BY TABLE 3.2 AND EQUATION (3.12).  NOTE THAT 1 OVER MU
!  SUB E IS EQUAL TO EMU.
! HERE WE USE THE NOTATION BE7 + PROTON = BE7PROTON AND
!  BE7 + E = BE7ELECTRON.
! F1 IS THE FRACTION OF THE BE7 THAT IS BURNED BY ELECTRON CAPTURE.
! F2 IS THE FRACTION OF THE BE7 THAT IS BURNED BY PROTON CAPTURE.
! F3 IS THE FRACTION OF THE N14 THAT IS BURNED BY P, ALPHA REACTION.
! F4 IS THE FRACTION OF THE N15 THAT IS BURNED BY P, GAMMA REACTION.
! SEE TABLE 21 OF BAHCALL AND ULRICH (1988), REV. MOD. PHYS. 60.
! 10/13/97. I changed Temp3 (i.e., tau) by 5/10^5 as a result of
! using bare nuclear masses. Previously coefficient was -10.26202.
! I did not change Be7electron for the slightly different S0 for
! electron capture, since that is done in SStandard. Previous coefficient
! in Be7electron expression was (3.126571E+5). 10/14/97.
!
         be7electron = (1.752e-10)*t9_m12*(1.0 + 0.004*(1000.*t9 - 16.))
         be7electron = be7electron*electron_mean_weight_inverse*cross_section_scale(15)
         temp3 = (-10.2625*t9_m13)
         be7proton = (3.128813e+5)*hydrogen_fraction*cross_section_scale(16)*exp(temp3)
! INCLUDE FOR BE7PROTON THE T9M23 FACTOR AND ALL CORRECTIONS PROPORTIONAL TO
!  Q1,...,Q5 FROM EQUATION 3.14 OF NEUTRINO ASTROPHYSICS. THESE
!  CORRECTIONS ARE DEFINED EARLIER IN THIS SUBROUTINE.
!         QRBE7 = T9M23 + Q1(8)*T9M13 + Q2(8)+ Q3(8)*T9P13
!     $          + Q4(8)*T9P23 + Q5(8)*T9
! MHP 9/14 ADDED THE ABILITY TO ALTER DERIVATIVES INDEPENDENTLY
         qrbe7 = t9_m23 + q1(8)*t9_m13 + qs0e_scale(8)*(q2(8)+ q3(8)*t9_p13) &
              + qqs0ee_scale(8)*(q4(8)*t9_p23 + q5(8)*t9)
         be7proton = be7proton*qrbe7
! CALCULATE THE SCREENING CORRECTION FOR BE7 + P REACTION.  USE WEAK AND
!  INTERMEDIATE SCREENING FORMULAE.
         zprdbe7p = 4.0
         z86be7p = 5.7790
         uwk = tm1*zprdbe7p
         if (uwk.le.weak_screening_threshold) then
            utotbe7p = uwk
         else
            uint = 0.38*xxl8*xtr*z86be7p/(ion_mean_weight_inverse*z58*z28)
            utotbe7p = uint
          end if
         be7proton = be7proton*exp(utotbe7p)
! END OF CALCULATION OF SCREENING CORRECTION FOR BE7 + P REACTION.
! MULTIPLY RATES BY FACTOR OF RHO/[(ATOMIC MASS UNIT)*A(BE7)] TO GET
!  IN UNITS OF GM^{-1}.  CALL FACTOR CAMUBE7. Corrected for Be7 bare
! mass number. Previously used 8.582295E+22 in the expression below.
! These corrections were made on 10/14/97.
!
         camube7 = density*8.584981e+22
         be7proton = camube7*be7proton
         be7electron = camube7*be7electron
! END OF MULTIPLICATION INSERTED NOVEMBER 6, 1990.
         f1 = be7electron/(be7electron + be7proton)
         f2 = be7proton/(be7electron + be7proton)
! *********************************************************************
! END OF CALCULATION OF CRUCIAL BE7 ELECTRON CAPTURE AND PROTON CAPTURE
!  RATES AND THEIR RATIO.
! *********************************************************************
! ***************************
! N15 + P BRANCHING.
! ***************************
! THE FOLLOWING STATEMENTS COMPUTE THE BRANCHING OF N15(P,ALPHA)C12
!  AND N15(P,GAMMA)O16. THESE STATEMENTS REPLACE OUTDATED STATEMENTS
!  IN THE YALE CODE. THE CNO CROSS-SECTION FACTORS ARE FROM TABLE 3.4
!  OF NEUTRINO ASTROPHYSICS.  THE RATIO OF THE REACTIONS DEPENDS ONLY
!  UPON THE EFFECTIVE ZERO ENERGY S-FACTOR, WHICH IS S0(ZERO ENERGY)
!  TIMES THE COMBINATION OF TEMPERATURE AND S-FACTOR DERIVATIVES/S0
!  THAT WAS USED PREVIOUSLY AS R1 IN THE RATE CALCULATIONS. THE
!  NUMERICAL COEFFICIENTS THAT APPEAR IN THE RATE WERE REPRESENTED
!  BY THE Q1(J),...Q5(J) FOR THE OTHER REACTIONS.
!  THE QVALUES FOR THE N15 REACTIONS HAVE BEEN COMPUTED SEPARATELY.
! DBG 8/94 APPLIED MHP UPDATE TO NUCLEAR REACTIONS
! To agree with Solar Fusion Workshop paper, the value of S0 in keV-b
! has been changed from 78000. to 675000. on 9/25/97. On 10/14/97,
! JNB changed cues.f so as to compute the Q-coefficients for the N15 + p
! reactions.  Also, checked that the coefficients are the same as the
! ones given earlier in energy.f when we use the older CNO data.
!
      o16gamma = t9_m23 + 0.0273016*t9_m13 + 0.14374 + 0.027490*t9_p13 &
                + 6.14685*t9_p23 + 2.98940*t9
! MULTIPLY BY THE VALUE OF S0 IN KEV-B.
!      O16GAMMA = O16GAMMA*64.
! MHP 9/14 ADDED THE OPTION TO MODIFY THE RELATIVE CROSS SECTIONS
! FOR N15+P -> C12+ALPHA AND O16+GAMMA
      o16gamma = o16gamma*64*o16_gamma_scale
!
      c12alpha = t9_m23 + 0.0273016*t9_m13 + 2.01186 + 0.384763*t9_p13 &
                + 17.0579*t9_p23 + 8.29580*t9
!      C12ALPHA = C12ALPHA*67500
      c12alpha = c12alpha*67500*c12_alpha_scale
      f3 = c12alpha/(c12alpha + o16gamma)
      f4 = 1.0d0 - f3
! END OF NEW ROUTINE FOR THE BRANCHING OF N15 + P .
   50 do 60 i=nz,nrxns
         reaction_rate(i)=0.
         dlnrate_dlnrho(i)=0.
         dlnrate_dlnt(i)=0.
   60 continue
! ***MHP 3/91 ALPHA CAPTURE REACTIONS UPDATED TO CAUGHLAN AND FOWLER(1988)
!    RATES.  THE RATES ARE EXPRESSED IN THE SAME TERMS USED BY CZ, WITH
!    THE CONVERSION FACTOR IN THE FRONT OBTAINED FROM VANDENBERG'S
!    NOTES ON THE REACTION RATES.
!  RATE(8)  HE4+C13
!  RATE(10) HE4+C12=>O16
!  RATE(11) HE4+N14=>O18
!  RATE(12) TRIPLE ALPHA
      if (log_temperature.lt.tcut(4)) go to 100
! C13(ALPHA,N) O16
      r1=t9_m23+0.0129d0*t9_m13+2.04d0+0.184d0*t9_p13
      a1 = 6.77d15*exp(-32.329d0*t9_m13-(t9/1.284d0)**2)
      a2 = 3.82d5*exp(-9.373*t9_m1)
      a3 = 1.41d6*exp(-11.873*t9_m1)
      a4 = 2.00d9*exp(-20.409*t9_m1)
      a5 = 2.92d9*exp(-29.283*t9_m1)
      reaction_rate(8) = 1.157126d22*density*exp(screening_factor(8))* &
           (a1*r1+t9_m32*(a2+a3+a4+a5))
      dlnrate_dlnrho(8)=1.0d0+dscreen_dlnrho(8)
      dr1 = cc13*(-2.0d0*t9_m23-0.0129d0*t9_m13+0.184d0*t9_p13)
      da1 = a1*(cc13*32.329d0*t9_m13 - 2.0d0*(t9/1.284d0)**2)
      dlnrate_dlnt(8) = dscreen_dlnt(8)+density/reaction_rate(8)*(dr1*a1 + &
           r1*da1 + &
           a2*(9.373*t9_m1-1.5d0)+a3*(11.873*t9_m1-1.5d0)+ &
           a4*(20.409*t9_m1-1.5d0)+a5*(29.283*t9_m1-1.5d0))
! C12(ALPHA,GAMMA)O16
      r1 = 1.0d0/(1.0d0+0.0489d0*t9_m23)
      r2 = 1.0d0/(1.0d0+0.2654d0*t9_m23)
      a1 = t9_m2*exp(-32.120*t9_m13)
      a2 = 1.04d8*r1**2*exp(-(t9/3.496)**2)
      a3 = 1.76d8*r2**2
      a4 = 1.25d3*t9_m32*exp(-27.499*t9_m1)
      a5 = 1.43d-2*t9**5*exp(-15.541*t9_m1)
      reaction_rate(10) = 1.25388d22*density*exp(screening_factor(10))* &
           (a1*(a2+a3)+a4+a5)
      dlnrate_dlnrho(10) = 1.0d0+dscreen_dlnrho(10)
      dlnrate_dlnt(10) = dscreen_dlnt(10)+density/reaction_rate(10)* &
           (a1*((cc13*32.120*t9_m13-2.0d0)* &
           (a2+a3)+a2*(r1*cc13*0.1956-2.0d0*(t9/3.496)**2)+a3* &
           (cc13*1.0616d0*r2))+a4*(27.499*t9_m1-1.5d0)+a5* &
           (5.0d0+15.541*t9_m1))
! N14(ALPHA,GAMMA)F18 + F18=>O18+EPLUS+NU
      r1 = t9_m23+0.012d0*t9_m13+1.45d0+0.177d0*t9_p13+1.97d0*t9_p23 &
           +0.406d0*t9
      a1 = 7.78d9*exp(-36.031d0*t9_m13-(t9/0.881d0)**2)
      a2 = t9_m32*2.36d-10*exp(-2.798d0*t9_m1)
      a3 = t9_m32*2.03d0*exp(-5.054d0*t9_m1)
      a4 = t9_m23*1.15d4*exp(-12.310*t9_m1)
      reaction_rate(11)= 1.07452d22*density*exp(screening_factor(11))* &
           (a1*r1+a2+a3+a4)
      dlnrate_dlnrho(11)=1.+dscreen_dlnrho(11)
      dr1 = cc13*(-2.0d0*t9_m23-0.012d0*t9_m13+0.177d0*t9_p13+ &
            3.94d0*t9_p23)+0.406d0*t9
      da1 = a1*(cc13*36.031d0*t9_m13-2.0d0*(t9/0.881d0)**2)
      dlnrate_dlnt(11) = dscreen_dlnt(11)+density/reaction_rate(11)*(dr1*a1+ &
              r1*da1+a2* &
              (2.798d0*t9_m1-1.5d0)+a3*(5.054d0*t9_m1-1.5d0)+ &
              a4*(12.310d0*t9_m1-cc23))
! TRIPLE ALPHA
      reaction_rate(12) = 1.565315d21*density**2*t9_m1*t9_m2*2.79e-8* &
                 exp(-4.4027*t9_m1+screening_factor(12))
      dlnrate_dlnrho(12) = 2.0d0+dscreen_dlnrho(12)
      dlnrate_dlnt(12) = -3.0d0+dscreen_dlnt(12)+4.4027d0*t9_m1
! *******************
! EG(I)
! *******************
! MULTIPLY THE RATES PER GRAM, RATE(I), BY THE ABUNDANCES OF THE
!  REACTING SPECIES BY MASS, TO GET THE TOTAL RATES PER GRAM, EG.
  100 eg(1)=reaction_rate(1)*hydrogen_fraction*hydrogen_fraction
! MHP 5/02 ADD DEUTERIUM BURNING IF RELEVANT
      if (deuterium_fraction.gt.1.0d-11) then
         egdeut = rdeut*hydrogen_fraction*deuterium_fraction
      else
         egdeut = 0.0d0
      end if
      eg(2)=reaction_rate(2)*he3_fraction*he3_fraction
      eg(3)=reaction_rate(3)*he3_fraction*helium_fraction
      eg(4)=reaction_rate(4)*hydrogen_fraction*c12_fraction
      eg(5)=reaction_rate(5)*hydrogen_fraction*c13_fraction
      eg(6)=reaction_rate(6)*hydrogen_fraction*n14_fraction
      eg(7)=reaction_rate(7)*hydrogen_fraction*o16_fraction
      eg(8)=reaction_rate(8)*helium_fraction*c13_fraction
!     EG(9)=RATE(9)*Y*X016
      eg(10)=reaction_rate(10)*helium_fraction*c12_fraction
      eg(11)=reaction_rate(11)*helium_fraction*n14_fraction
      eg(12)=reaction_rate(12)*helium_fraction**3
!     EG(13)=RATE(13)*XC12*XC12
! ******************************************************************
! ****************************************
! ENERGY GENERATION.
! ****************************************
! CALCULATE ENERGY GENERATION BY MULTIPLYING RATES PER GRAM PER SEC BY
!  THE ENERGY RELEASE.  THE ENERGIES ARE TAKEN FROM TABLE 21 OF BAHCALL AND
!  ULRICH (1988), REV. MOD. PHYS. 60, 297. THIS TABLE IS BASED UPON A CAREFUL
!  CALCULATION OF THE AVERAGE AMOUNT OF ENERGY LOSS BY NEUTRINOS FOR
!  EACH REACTION. THE NUMBERS FOR THE C12 + P REACTION SEQUENCE AND THE
!  C13 + P REACTION ARE BROKEN DOWN SEPARATELY FOR THIS
!  SUBROUTINE.
! THE FINAL NUMBERS ARE IN ERG PER GM PER SECOND.
! DEFINE THE CONSTANT TO CONVERT MEV'S TO ERGS. THE NUMBERS THAT APPEAR
!  ARE IN MEV SO THEY CAN BE EASILY IDENTIFIED.
      convert = 1.602177e-6
! THE MULTIPLYING CONSTANTS BELOW ARE IN MEV.
! JNB changed the pp energy release by 0.002 MeV because of a better
! estimate of the neutrino energy loss on 9/25/97. See pg. 139 of
! Vol. 19 of my notes.  On pgs. 139-141, I document other small changes
! to this energy generation. No large changes; all of order keV changes
! except for the rare 8B reaction.  9/28/97.
!
      reaction_energy_gen(1)=eg(1)*6.664*convert
! MHP 5/02 ADD DEUTERIUM BURNING
      if (deuterium_fraction.gt.1.0d-11) then
         dgdeut = egdeut*qdeut*convert
      else
         dgdeut = 0.0d0
      end if
      reaction_energy_gen(2)=eg(2)*12.860*convert
      reaction_energy_gen(3)=eg(3)*(1.586+f1*17.394+f2*11.499)*convert
      reaction_energy_gen(4)=eg(4)*3.457372*convert
      reaction_energy_gen(5)=eg(5)*7.550628*convert
      reaction_energy_gen(6)=eg(6)*(9.054+f3*4.966+f4*12.128)*convert
      reaction_energy_gen(7)=eg(7)*3.553*convert
      reaction_energy_gen(8)=eg(8)*2.216*convert
      reaction_energy_gen(10)=eg(10)*7.162*convert
      reaction_energy_gen(11)=eg(11)*5.815*convert
      reaction_energy_gen(12)=eg(12)*7.275*convert
! JVS 10/11 Need to grab He3 energy generation
      he3_luminosity_placeholder = reaction_energy_gen(2)
      he3_total_placeholder = reaction_energy_gen(2)+reaction_energy_gen(3)
! JVS end

! *******************************************************************
! END OF CALCULATION OF ENERGY RELEASE.
! *******************************************************************
! SET TO ZERO O16+ALPHA AND C12+C12 RATES.
      dlnrate_dlnrho(9) = 0.0
      dlnrate_dlnt(9) = 0.0
      reaction_energy_gen(9) = 0.0
      dlnrate_dlnrho(13) = 0.0
      dlnrate_dlnt(13) = 0.0
      reaction_energy_gen(13) = 0.0
! END OF XEROING OUT OF REACTIONS 9 AND 13.
      total_energy_gen_rate=0.0
      sum2=0.0
      sum3=0.0
      do 110 i=1,nrxns
! *******************************************************************
! SUM OF THE TOTAL ENERGY GENERATION IN ERGS PER GRM PER SECOND WITH
! DERIVATIVES WITH RESPECT TO DENSITY AND TO TEMPERATURE.
! *******************************************************************
! SUM1 = SUM OF ALL ENERGY GENERATION RATES. NOTE THAT THE BURNING OF
!  BE7 IS INCLUDED IN DG(3) ABOVE.
! SUM2 = SUM OVER I OF DG(I)* [D LOG RATE(I) / D LOG RHO ].
! SUM3 = SUM OVER I OF DG(I)* [D LOG RATE(I) / D LOG T ].
         total_energy_gen_rate=total_energy_gen_rate+reaction_energy_gen(i)
         sum2=sum2+reaction_energy_gen(i)*dlnrate_dlnrho(i)
         sum3=sum3+reaction_energy_gen(i)*dlnrate_dlnt(i)
  110 continue
! MHP 5/02 ADD DEUTERIUM BURNING
!      IF(IU.LE.10)THEN
!         WRITE(*,911)IU,SUM1,DGDEUT,SUM2,SUM3,QRTDEUT
!      ENDIF
! 911  FORMAT(I5,1P5E13.3)
      total_energy_gen_rate = total_energy_gen_rate + dgdeut
      sum2 = sum2 + dgdeut
      sum3 = sum3 + dgdeut*qrtdeut
      if (total_energy_gen_rate.le.1.e-12) then
         en=-20.
         dlnepsilon_dlnrho=0.
         dlnepsilon_dlnt=0.
         do 120 i=1,nrxns
            eg(i)=0.
  120    continue
      else
! ******************************************************
! GLOBAL QUANTITIES THAT ARE RETURNED BY THE SUBROUTINE.
! ******************************************************
! PEP AND PET ARE THE DERIVATIVES OF THE TOTAL ENERGY GENERATION RATE
!  WITH RESPECT TO DENSITY AND TEMPERATURE.
! MHP 5/90 CHANGE DERIVATIVES TO BE D LN EPS/D LN RHO AND D LN EPS/D LN T
! TO PUT THEM IN THE SAME FORM AS PRATHER DERIVATIVES.
!        EN=DLOG10(SUM1)
!        PEPD=SUM2/SUM1
!        PEP=PDP*PEPD
!        PET=SUM3/SUM1+PDT*PEPD
         dlnepsilon_dlnrho = sum2
         dlnepsilon_dlnt = sum3
      end if
! PDP = D LOG RHO/ D LOG P; PDT = D LOG RHO/ D LOG T.
! *****************************************************
! END OF COMPUTATION OF THE GLOBAL QUANTITIES.
! *****************************************************
      do 130 i=1,nrxns
         if (reaction_rate(i).le.1.e-5) reaction_rate(i) = 0.0
  130 continue
! ******************************************************
! RATES PER 10^9 YEARS PER ATOMIC MASS UNIT: HRK(IU)
! ******************************************************
! HR1, ..., HR13 ARE THE RATES OF THE INDIVIDUAL REACTIONS.
!  THE INTERPRETATION OF WHICH REACTION GOES WITH WHICH SYMBOL CAN BE
!  MADE EASILY BY LOOKING AT THE DEFINITIONS OF THE EG(I)'S.
!  THE ABUNDANCES ARE UPDATED IN SUBROUTINE KEMCOM USING THESE MATRICES.
! C21 IS THE PRODUCT OF (10^9 YEARS/1 SECOND)*(1 ATOMIC MASS UNIT/1
!  GRAM). I HAVE USED HERE SIDEREAL YEAR IN CONVERTING TO SECONDS.
  200 reaction_rate_1(shell_index)=reaction_rate(1)*years_per_sec_over_amu
      reaction_rate_2(shell_index)=reaction_rate(2)*years_per_sec_over_amu
      reaction_rate_3(shell_index)=reaction_rate(3)*years_per_sec_over_amu
      reaction_rate_4(shell_index)=reaction_rate(4)*years_per_sec_over_amu
      reaction_rate_5(shell_index)=reaction_rate(5)*years_per_sec_over_amu
      reaction_rate_6(shell_index)=reaction_rate(6)*years_per_sec_over_amu
      reaction_rate_7(shell_index)=reaction_rate(7)*years_per_sec_over_amu
      reaction_rate_8(shell_index)=reaction_rate(8)*years_per_sec_over_amu
      reaction_rate_9(shell_index)=reaction_rate(9)*years_per_sec_over_amu
      reaction_rate_10(shell_index)=reaction_rate(10)*years_per_sec_over_amu
      reaction_rate_11(shell_index)=reaction_rate(11)*years_per_sec_over_amu
      reaction_rate_12(shell_index)=reaction_rate(12)*years_per_sec_over_amu
      reaction_rate_13(shell_index)=reaction_rate(13)*years_per_sec_over_amu
      n15_alpha_branch_fraction(shell_index)=f3
      be7_electron_capture_fraction(shell_index)=f1
! ****************************************
! END OF COMPUTATION OF HRK(IU).
! ****************************************
! ****************************************
! CALCULATING THE TOTAL ENERGY GENERATION.
! ****************************************
! THE SUMMATION OF THE ENERGIES IS GIVEN IN TABLE 21 OF NEUTRINO
!  ASTROPHYSICS.
! THE ORIGINAL YALE SUBROUTINE CONTAINED SERIOUS ERRORS.  THE
!  CALCULATION OF THE RATE OF ENERGY THROUGH EPP2 AND EPP3 (SEE
!  BELOW) CONTAINED TWO FACTORS OF BRANCHING RATIOS, RATHER THAN
!  THE SINGLE FACTOR THAT SHOULD BE PRESENT. THIS HAD THE EFFECT
!  OF REDUCING ARTIFICIALLY THE ENERGY CALCULATED FROM THESE
!  REACTIONS.
! EPP1 INCLUDES THE ENERGY GENERATED BY THE PP REACTION, BY THE H2 + P
!  REACTION, AND BY THE HE3 + HE3 REACTION.  SEE TABLE 21 OF NEUTRINO
!  ASTROPHYSICS.
      pp_chain_energy_gen = reaction_energy_gen(1)+reaction_energy_gen(2)+dgdeut
!      EPP1 = DG(1)+DG(2)
! EPP3 INCLUDES THE ENERGY GENERATED BY THE HE3 + HE4 REACTION AND BY
!  THE BURNING OF BE7 THROUGH PROTON CAPTURE.
!      EPP3 = EG(3)*(1.586 + F2*11.499)*CONVERT
      he3he4_be7_proton_energy_gen = eg(3)*f2*(1.586 + 11.499)*convert
! EPP2 INCLUDES THE ENERGY GENERATED BY THE HE3 + HE4 REACTION AND BY
!  THE BURNING OF BE7 THROUGH ELECTRON CAPTURE.
!      EPP2 = EG(3)*(1.586 + F1*17.394)*CONVERT
      he3he4_be7_electron_energy_gen = reaction_energy_gen(3) - &
           he3he4_be7_proton_energy_gen
! ECN IS THE ENERGY GENERATED THROUGH THE CNO CYCLE.
      cno_cycle_energy_gen=reaction_energy_gen(4)+reaction_energy_gen(5)+ &
           reaction_energy_gen(6)+reaction_energy_gen(7)
! E3AL IS THE ENERGY GENERATED THROUGH THE TRIPLE-ALPHA REACTION AND
!  IS NEGLIGIBLE FOR THE SUN.
      triple_alpha_energy_gen = reaction_energy_gen(12)


! ENERGY FROM ALPHA CAPTURE REACTIONS.
      alpha_capture_energy=reaction_energy_gen(8)+reaction_energy_gen(10)+ &
           reaction_energy_gen(11)
      if (lsnu) then
! MHP 9/91 CHANGE TO TURN OFF NEUTRINO CALC FOR HYDROGEN-EXHAUSTED CORE.
         if (hydrogen_fraction.le.1.0d-6) then
            do i=1,10
               neutrino_flux(i)=0.0d0
            end do
         else
! ****************************************************************
! CALCULATION OF NEUTRINO FLUXES
! ****************************************************************
! THIS PART OF THE SUBROUTINE CALCULATES THE NEUTRINOS FLUXES IN
!  NUMBER PER GRAM PER SQUARE CENTIMETER PER SECOND AT THE EARTH'S SURFACE
!  (ASSUMING NOTHING HAPPENS TO THE NEUTRINOS AFTER THEY ARE CREATED).
! SEE TABLES 3.1 AND 3.3 OF NEUTRINOS ASTROPHYSICS OR EQUATIONS 6.1-6.8
!  FOR THE REACTIONS. THE ORDER OF THE REACTIONS IS THE SAME AS IN
!  EQUATIONS 6.1-6.8 .
! DEFINE 4*PI*(AU)**2 .
         fourpiau2 = 2.812295e+27
! FLUX OF PP NEUTRINOS.
         neutrino_flux(1) = eg(1)/fourpiau2
! FLUX OF PEP NEUTRINOS. USE EQUATION 3.17 OF NEUTRINO ASTROPHYSICS.
! Note that should not change SStandard(14) unless the ratio of pep to pp
!  is changed.  Pep rate is explicitly scaled here with respect to the pp
!  rate.
         neutrino_flux(2) = (3.4848e-6)*electron_number_density_na*t9_m12* &
              (1.0 + 20.*t9)*eg(1)
         neutrino_flux(2) = neutrino_flux(2)*cross_section_scale(14)/fourpiau2
! FLUX OF HEP NEUTRINOS.  USE EQUATION 3.12 DIRECTLY.
         q6hep = -6.1399
! Q6 IS THE NEGATIVE OF THE COEFFICIENT OF T9M13 IN TAU, EQUATION 3.10.
         neutrino_flux(3) = (1.71724e+11)*density*t9_m23*exp(q6hep*t9_m13)
! THE DERIVATIVES OF THE CROSS SECTION FACTOR ARE NOT KNOWN AND ARE
!  TAKEN TO BE ZERO.  THE ONLY TERM FROM EQUATION 3.14 THAT SURVIVES
!  IS 5/(12*TAU).
         neutrino_flux(3) = (1.0 + 0.067862*t9_p13)*cross_section_scale(17)* &
              neutrino_flux(3)
! CALCULATE WEAK OR INTERMEDIATE SCREENING FOR HEP NEUTRINOS.
         zprdhe3p = 2.0
         z86he3p = 3.08687
         uwk = tm1*zprdhe3p
         if (uwk.le.weak_screening_threshold) then
            utothe3p = uwk
         else
            uint = 0.38*xxl8*xtr*z86he3p/(ion_mean_weight_inverse*z58*z28)
            utothe3p = uint
         end if
! END OF CALCULATION OF SCREENING CORRECTION FOR HE3 + P REACTION.
         neutrino_flux(3) = neutrino_flux(3)*exp(utothe3p)
         neutrino_flux(3) = neutrino_flux(3)*hydrogen_fraction*he3_fraction/ &
              fourpiau2
! COMPUTE BE7MASSFRACTION. THIS IS NOT REQUIRED FOR THE NEUTRINO
!  FLUXES SINCE BE7 IS ALWAYS IN EQUILIBRIUM WITH THE SLOWER PRODUCTION
!  RATE OF HE3 + HE4.  HOWEVER, IT IS OF INTEREST IN SOME APPLICATIONS
!  TO KNOW THE BE7 MASS FRACTION, SO I COMPUTE IT HERE AND IT CAN BE
!  EXTRACTED WITH A COMMON STATEMENT IF DESIRED.
         be7_mass_fraction = eg(3)/(be7proton + be7electron)
! END OF NOVEMBER 6, 1990  ADDITION.
! FLUX OF BE7 NEUTRINOS.
         neutrino_flux(4) = eg(3)*f1/fourpiau2
! FLUX OF B8 NEUTRINOS.
         neutrino_flux(5) = eg(3)*f2/fourpiau2
! FLUX OF N13 NEUTRINOS.
         neutrino_flux(6) = eg(4)/fourpiau2
! FLUX OF O15 NEUTRINOS.
         neutrino_flux(7) = eg(6)/fourpiau2
! FLUX OF F17 NEUTRINOS.
         neutrino_flux(8) = eg(7)/fourpiau2
! FLUX OF FICTIONAL HE3 + HE3 NEUTRINOS.
         neutrino_flux(9) = eg(2)/fourpiau2
! FLUX OF FICTIONAL HE3 + HE4 NEUTRINOS.
         neutrino_flux(10) = eg(3)/fourpiau2
! SET UNITS OF NEUTRINO FLUXES TO BE 10**10 PER CM^2 PER SEC PER GM AT THE
!  EARTH. MULTIPLY BY 10**-10.
!  IF THE VALUE FOR THIS SHELL IS NEGLIGIBLY SMALL, SET EQUAL TO ZERO.
         do k = 1,10
            neutrino_flux(k) = (1.0e-10)*neutrino_flux(k)
            flux_value = neutrino_flux(k)
! KC 2025-05-30 CHANGED 1.E-50 TO 0.0 TO AVOID UNDERFLOW
            if (flux_value.le.0.0) then
              neutrino_flux(k) = 0.0
            end if
         end do
! MHP 9/91 ENDIF INSERTED HERE.
         end if
      end if
! END OF NEUTRINO FLUX ROUTINE.
!
! ***MHP 5/91
!C CALCULATE NEUTRINO LOSSES FOR NEUTRINO-COOLED CORES OF EVOLVED STARS.
! 3/92 DBG Added option to use new (more sophisticated) neutrino loss
! routines.  See subroutine NEUTR for complete description.


      if (log_temperature.le.tcut(5)) return


          carbon_fraction_total = c12_fraction+c13_fraction
          oxygen_fraction_total = o16_fraction+o18_fraction
          neutrino_temp=10.0**log_temperature
          neutrino_density=10.0**log_density


!**** Itoh 1996 Neutrino loss routines - Grant Newsham 9/06 *****


      if (use_itoh_neutrino_loss) then



          call neutrino(neutrino_temp,neutrino_density,hydrogen_fraction, &
               helium_fraction,carbon_fraction_total,oxygen_fraction_total, &
               neutrino_loss_snu,neutrino_dlnq_dlnt,neutrino_dlnq_dlnd)


          neutrino_loss_rate = -neutrino_loss_snu
          neutrino_dlnq_dlnt = -neutrino_dlnq_dlnt*neutrino_temp/neutrino_loss_rate
          neutrino_dlnq_dlnd = -neutrino_dlnq_dlnd*neutrino_density/neutrino_loss_rate


          total_energy_gen_rate = total_energy_gen_rate + neutrino_loss_rate


          dlnepsilon_dlnrho = dlnepsilon_dlnrho + neutrino_dlnq_dlnd
          dlnepsilon_dlnt = dlnepsilon_dlnrho + neutrino_dlnq_dlnt


!****************************************************************


       else



!     THESE ARE OLD NEUTRINO LOSS ROUTINES


         el = t9/5.9302
         eli = 1.0/el
         ez = exp(cc13*ln10*(dd-9.0))*eli*(0.7937+0.2063*hydrogen_fraction)
         ez3 = ez**3
         emue = 0.5*(1.0+hydrogen_fraction)
         ex1 = 0.0
         if (t9.ge.0.2) then
!C PAIR NEUTRINOS
            el2 = el*el
            polx10=(1.+el2*(-13.04+el2*(133.5+el2*(1534.+el2*918.6))))
            polx11 = v1(1) + ez*(v1(2) + ez*v1(3))
            polx12 = ez3 + eli*(v1(4) + eli*(v1(5) + eli*v1(6)))
            ex1 = dexp(-ez*v1(7)-eli-eli-ln10*dd)*polx10*polx11/polx12
         end if
!C PHOTO NEUTRINOS
         polx21 = v2(1) + ez*(v2(2) + ez*v2(3))
         polx22 = ez3 + eli*(v2(4) + eli*(v2(5) + eli*v2(6)))
         ex2 = emue*el**5*dexp(-ez*v2(7))*polx21/polx22
!C PLASMA NEUTRINOS
         polx31 = v3(1) + ez*(v3(2) + ez*v3(3))
         polx32 = ez3 + eli*(v3(4) + eli*(v3(5) + eli*v3(6)))
         ex3 = emue**3*dexp(-ez*v3(7)+ln10*(dd+dd))*polx31/polx32
         neutrino_loss_rate = -(ex1 + ex2 + ex3)
         total_energy_gen_rate = total_energy_gen_rate + neutrino_loss_rate
         qetnx = 0.0
         qednx = 0.0
         if (t9.ge.0.2) then
! MHP 10/02 fixed column 72 problem
            qedn= ez*(v1(2)+2.*ez*v1(3))/polx11 - ez*v1(7) &
                  - 3.0d0*ez3/polx12
            qetn= -qedn+ eli*(v1(4)+eli*(2.*v1(5)+3.*eli*v1(6)))/polx12
            qetnx= el2*(-26.08+el2*(534.+el2*(9204.+el2*7348.8)))/polx10
            qetnx = (qetnx + qetn + eli+eli)*ex1
            qednx = (-1.0 +cc13*qedn)*ex1
         end if
         qedn= ez*(v2(2)+2.*ez*v2(3))/polx21 - ez*v2(7) - 3.*ez3/polx22
         qetn= -qedn+ eli*(v2(4)+eli*(2.*v2(5)+3.*eli*v2(6)))/polx22
         qetnx = qetnx + (5.0 + qetn)*ex2
         qednx = qednx +cc13*qedn*ex2
         qedn= ez*(v3(2)+2.*ez*v3(3))/polx31 - ez*v3(7) - 3.*ez3/polx32
         qetn=-qedn + eli*(v3(4)+eli*(2.*v3(5)+3.*eli*v3(6)))/polx32
         qetnx = qetnx + qetn*ex3
         qednx = qednx + (2.0 +cc13*qedn)*ex3



         dlnepsilon_dlnt = dlnepsilon_dlnt - qetnx
         dlnepsilon_dlnrho = dlnepsilon_dlnrho - qednx



      end if


      return
end subroutine engeb
