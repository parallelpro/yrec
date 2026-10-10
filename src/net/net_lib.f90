!----------------------------------------------------------------------
! net_lib
!----------------------------------------------------------------------
! Aggregates YREC's COMMON-free nuclear-physics utility routines into a
! single module, following the numerics_lib/const_lib convention (see
! GUIDELINES.md, and MESA's own one small public *_lib module per
! domain). Each routine below was originally its own file
! (nuclear/<name>.f90, with per-file conversion notes preserved below);
! bodies are unchanged, only concatenated under one module. Callers use
! this via `use net_lib`. neutrino_losses_for_mixture.f90 contained two subroutines
! (neutrino_losses_for_mixture, azbar); both are included here.
!
! lir and ratext, originally also housed in nuclear/ (and formerly
! part of this module), moved out to numerics_lib (2026): both are
! generic numerical utilities with no nuclear-physics content -- they
! were only filed here because burn_li_be (now in burn_lib) happens to be
! one of their callers. See numerics_lib.f90 and GUIDELINES.md's rule
! that folder/module placement should track function, not caller.
! ---------------------------------------------------------------------
! MAP OF THE MODULE (2026; updated at the physics-purity split).
! net_lib now holds ONLY the pure kernels -- functions of
! (logRho, logT, composition, controls), no model state, the surface
! test_nuclear pins:
!     rates      13 reaction rates + branching fractions at one point
!                (consumes cross_section_scale from setup/map_user_inputs.f90).
!     sneut      Itoh et al. (1996) neutrino losses (pure fits).
!     neutrino_losses_for_mixture   azbar + sneut wrapper for a full mixture.
!     neutrino   legacy interface to the loss rates.
!     azbar      abar/zbar/molar-abundance bookkeeping.
!     safe_divide_by_exp division-with-exponent guard.
!     ifermi12, zfermim12   Fermi-Dirac integral inverses (rates).
! The BURN DRIVERS (compute_zone_energy/equilibrium_burn_rates, dburn/dburnm/deutrate,
! burn_li_be/burn_li_be_midpoint/store_li_be_rates) -- everything that reads or writes
! star_info -- moved to core/burn_lib.f90: they advance the model and
! are star-layer code (MESA: struct_burn_mix, not net).
! ---------------------------------------------------------------------
module net_lib
! the thermal neutrino losses (Itoh et al. 1996 via Timmes) live in
! net/thermal_neutrinos.f90; re-exported here so net_lib stays the
! domain's one public entry
      use thermal_neutrino_lib
      implicit none

! ---- named indices for the positional array slots (2026) ----
! The 13-reaction arrays of rates (below) and compute_zone_energy (core/burn_lib.f90):
! rate/reaction_rate, screening_factor, charge_product, z53/z43/z23/z86,
! eg, reaction_energy_gen, ... and star%cross_section_scale(1:13).
! Slots 9 (O16 + alpha) and 13 (C12 + C12) are carried but never
! computed. num_reactions was NRXNS.
      integer, parameter, public :: &
           r_pp = 1, r_he3he3 = 2, r_he3he4 = 3, r_pc12 = 4, r_pc13 = 5, &
           r_pn14 = 6, r_po16 = 7, r_ac13 = 8, r_ao16_unused = 9, &
           r_ac12 = 10, r_an14 = 11, r_3alpha = 12, r_c12c12_unused = 13, &
           num_reactions = 13
! num_rate_rows: length of the one-zone rate vector written by rates
! (rows 1..13 are the reactions above in r_* order, times
! gyr_amu_per_sec_gram; rows 14 and 15 are the C12+alpha and Be7
! electron-capture branching fractions). The row indices are the rr_*
! parameters of rotation_scratch_lib, and num_rate_rows must equal the
! leading extent of rot_scr%reaction_rate_by_zone(15,json).
      integer, parameter, public :: num_rate_rows = 15
! rows of the composition rate vector filled by rates (rate_vec,
! rate_by_zone(:,k), rot_scr%reaction_rate_by_zone): the rates of the
! table reactions per Gyr per amu (originally HR1-HR13) and the two
! branching fractions (HF1, HF2). Rows rr_zero9 and rr_zero13 are
! zero placeholders (reactions 9 and 13 are never computed). (2026-10:
! moved here from rotation_scratch_lib -- they are the network's.)
      integer, parameter, public :: rr_pp = 1, rr_he3_he3 = 2, &
           rr_he3_he4 = 3, rr_c12_p = 4, rr_c13_p = 5, rr_n14_p = 6, &
           rr_o16_p = 7, rr_c13_alpha = 8, rr_zero9 = 9, &
           rr_c12_alpha = 10, rr_n14_alpha = 11, rr_triple_alpha = 12, &
           rr_zero13 = 13, rr_frac_c12_alpha = 14, &
           rr_frac_be7_electron = 15
! star%cross_section_scale(14:17) (SSTANDARD) carries four more
! cross-section factors that are not rate-array slots.
      integer, parameter, public :: &
           s_pep = 14, s_be7e = 15, s_be7p = 16, s_hep = 17
! q1..q5(8), star%qs0e_scale(8), star%qqs0ee_scale(8): slot 8 holds the
! Be7 + p coefficients (only for these arrays; rate slot 8 is He4+C13).
      integer, parameter, public :: iq_be7p = 8
! The 13-isotope mass tables (atomic_mass/atomic_charge) and the local
! mass-fraction arrays of rates and compute_zone_energy. num_isotopes was NELEM.
      integer, parameter, public :: &
           iso_n = 1, iso_h1 = 2, iso_h2 = 3, iso_h3 = 4, iso_he3 = 5, &
           iso_he4 = 6, iso_c12 = 7, iso_c13 = 8, iso_n14 = 9, &
           iso_o16 = 10, iso_o18 = 11, iso_ne20 = 12, iso_mg24 = 13, &
           num_isotopes = 13
! Species count of the H/He/C/O mixture handed to neutrino_losses_for_mixture/sneut.
! gyr_amu_per_sec_gram (the C21 of compute_zone_energy/rates/deutrate) lives in
! phys_const_lib.
! S0 in keV-b of the two N15 + p branches, multiplying the S_eff/S0
! polynomials o16_gamma_rate / c12_alpha_n15p_rate (and o16gamma /
! c12alpha in compute_zone_energy). Kept as integers, the kind of the original
! literals 64 and 67500 (both exact when promoted to double).
      integer, parameter, public :: s0_n15pg_kevb = 64, s0_n15pa_kevb = 67500

! d(p,g)3He (2026-10): Avogadro's number over the deuteron mass in amu
! turns CF88's N_A<sigma v> into the per-gram rate that deutrate (the
! deuterium depletion rate) and compute_zone_energy (its energy release) both use,
! and sets both convective-overturn caps. The two copies used to
! differ by 0.68%: deutrate had 3.0115d23 = 6.023d23/2 (A = 2), compute_zone_energy
! 6.023d23/2.013553 (the deuteron mass, as in compute_zone_energy's isotope table).
      double precision, parameter, public :: avogadro_dpg = 6.023d23
      double precision, parameter, public :: deuteron_mass_amu = 2.013553d0
! ---- the rate data shared by every consumer (2026-10: rates and
! burn_lib's compute_zone_energy each carried an identical copy; see
! nuclear_reaction_rates). Coefficient names and comments: network_rate_vector's header.
! S-factor fits of reactions 1-7 and of Be7(p,g) (q1..q5 slot iq_be7p),
! FCZ67 eq. 53 / Neutrino Astrophysics eq. 3.14 times T9^(-2/3); q6 the
! tau coefficient, q7 the log prefactor, q8 the Harris et al. T9^2 term.
      double precision, parameter, private :: &
           fit_q1(iq_be7p) = [0.12317d0, .03392d0, .0325d0, .0304d0, .03035d0, .0273d0, .02494d0, .040572d0], &
           fit_q2(iq_be7p) = [1.08749d0, -.273d0, -.2085d0, .7630d0, -0.4044d0, -1.60d0, -1.224d0, -0.2095d0], &
           fit_q3(iq_be7p) = [.93833d0, -.0648d0, -.0474d0, .1626d0, -.08598d0, -.3064d0, -.2139d0, -0.0595d0], &
           fit_q4(iq_be7p) = [0.d0, 0.d0, 0.d0, 4.79d0, 7.456d0, 0.0d0, .69703d0, 0.16762d0], &
           fit_q5(iq_be7p) = [0.d0, 0.d0, 0.d0, 2.595d0, 4.032d0, 0.0d0, 0.3097d0, 0.12114d0], &
           fit_q6(r_po16) = [-3.3804d0, -12.2757d0, -12.826d0, -13.6899d0, -13.7173d0, -15.2281d0, -16.6925d0], &
           fit_q7(r_po16) = [20.8964d0, 76.6003d0, 67.8036d0, 69.130d0, 70.3809d0, 69.8517d0, 70.8012d0], &
           fit_q8(r_po16) = [0.d0, 0.d0, 0.d0, 0.0d0, 0.0d0, 0.0d0, 0.d0]
! screening: charge products and the Graboske et al. (1973) table-4
! factors of each reaction slot
      double precision, parameter, private :: &
           screen_zprod(num_reactions) = [1.d0, 4.d0, 4.d0, 6.d0, 6.d0, 7.d0, 8.d0, 12.d0, 16.d0, 12.d0, 14.d0, 12.d0, 36.d0], &
           screen_z53(num_reactions) = [1.175d0, 3.73d0, 3.73d0, 4.804d0, 4.804d0, 5.385d0, 5.941d0, 9.014d0, 11.24d0, &
                9.014d0, 10.15d0, 9.104d0, 23.28d0], &
           screen_z43(num_reactions) = [0.52d0, 1.31d0, 1.31d0, 1.488d0, 1.488d0, 1.61d0, 1.721d0, 2.577d0, 3.025d0, &
                2.577d0, 2.81d0, 2.577d0, 5.668d0], &
           screen_z23(num_reactions) = [-0.413d0, -0.655d0, -0.655d0, -0.643d0, -0.643d0, -0.659d0, -0.673d0, -0.889d0, &
                -0.946d0, -0.889d0, -0.92d0, -0.889d0, -1.36d0], &
           screen_z86(num_reactions) = [1.630d0, 5.917d0, 5.917d0, 8.302d0, 8.302d0, 9.520d0, 10.716d0, 16.192d0, &
                20.978d0, 16.192d0, 18.606d0, 16.192d0, 45.6635d0]
! bare nuclear masses (amu) and charges of the iso_* isotopes
      double precision, parameter, public :: &
           isotope_mass(num_isotopes) = [1.008665d0, 1.007276d0, 2.013553d0, 3.015501d0, 3.014933d0, 4.001506d0, &
                11.996709d0, 13.000064d0, 13.999233d0, 15.990526d0, 17.994772d0, 19.986954d0, 23.978458d0], &
           isotope_charge(num_isotopes) = [0.d0, 1.d0, 1.d0, 1.d0, 2.d0, 2.d0, 6.d0, 6.d0, 7.d0, 8.d0, 8.d0, 10.d0, 12.d0]

! The rates of one zone, from nuclear_reaction_rates: rate(r_*) per gram
! per second with the abundance factors omitted (screened, scaled), and
! its logarithmic derivatives; the Be7 and N15 branchings; the two Be7
! channel rates and the screened hep factor (flux diagnostics); and the
! zone values the energy and flux consumers reuse. All zero below
! ctrl%nuclear_logT_cutoffs(1).
      type, public :: nuclear_rates
         double precision :: rate(num_reactions) = 0.0d0
         double precision :: dlnrate_dlnt(num_reactions) = 0.0d0
         double precision :: dlnrate_dlnrho(num_reactions) = 0.0d0
         double precision :: frac_be7_ecap = 0.0d0, frac_be7_pcap = 0.0d0
         double precision :: frac_n15_pa = 0.0d0, frac_n15_pg = 0.0d0
         double precision :: be7_electron = 0.0d0, be7_proton = 0.0d0
         double precision :: hep = 0.0d0
         double precision :: density = 0.0d0, t9 = 0.0d0, t9_p13 = 0.0d0, &
              t9_p23 = 0.0d0, t9_m13 = 0.0d0, t9_m23 = 0.0d0, t9_m12 = 0.0d0
         double precision :: electron_density_na = 0.0d0
      end type nuclear_rates

contains


!----------------------------------------------------------------------
! deuterium_proton_capture_rate
!----------------------------------------------------------------------
! The d(p,g)3He rate excluding the abundance factors, rho * N_A<sigma v>
! / m_d (Caughlan & Fowler 1988, the T9**(-2/3) form), for burn_lib's
! deutrate and nuclear_energy alike -- previously coded once in each. The arithmetic is
! compute_zone_energy's, token for token, so compute_zone_energy's value is unchanged. Screening is
! excluded, as before: the burning temperatures are low.
double precision function deuterium_proton_capture_rate(density, t9, t9_p13, t9_p23, t9_m13, t9_m23)
      use math_lib
      implicit none
      double precision, intent(in) :: density, t9, t9_p13, t9_p23, &
           t9_m13, t9_m23
      double precision :: zz, tfacdeut

      zz = -3.72d0*t9_m13
      tfacdeut = 1.0d0+0.112d0*t9_p13+3.38d0*t9_p23+2.65d0*t9
      deuterium_proton_capture_rate = density*2.240d3*t9_m23*exp(zz)*tfacdeut*avogadro_dpg/ &
           deuteron_mass_amu
end function deuterium_proton_capture_rate

!----------------------------------------------------------------------
! nuclear_reaction_rates
!----------------------------------------------------------------------
! The H/He reaction rates of one zone and their logarithmic derivatives
! with respect to T and rho (2026-10): the one evaluation behind both
! the composition rate vector (network_rate_vector, below) and the energy generation
! (burn_lib's compute_zone_energy), which used to carry identical copies of this
! block. mass_frac(iso_*) is the composition the screening sees. The
! arithmetic is compute_zone_energy's, token for token; the formulas, references and
! coefficient history are documented in network_rate_vector's header below.
!
! Each consumer applies its own floor afterwards (rates zeroes rates
! <= 1e-5, compute_zone_energy zeroes reactions 1-7 below 1e-30 before summing eps).
! No hydrogen: reactions 1-7 and the branchings are zero. Below
! nuclear_logT_cutoffs(1) everything is zero; the alpha captures need
! nuclear_logT_cutoffs(4).
!
! Provenance of the formulas and coefficients (moved verbatim from the
! body of rates, 2026-10, when this routine became the one copy):
! ANUC ARE ATOMIC MASS UNITS.
! On 10/13/97, JNB converted ANUC(I) from neutral atomic masses to
! bare nuclear masses.
! The scale is the mass of C12 divided by 12 or 931.49432 MeV,
! which is 1.6605402 times 10^{-24} gm. The values are obtained
! by dividing the mass excess (expressed in MeV) by 931.49432 MeV and
! adding to this the atomic mass number, A.  The value for Be7, which
! is used implicitly in this subroutine was, 7.016930, until 13/13/97.
! We now use the bare nuclear mass of 7.014735 .
! THE ISOTOPES ARE NEUTRON, H1, D, H3, HE3, HE4, C12, C13, N14, O16, O18,
! NE20, MG24, RESPECTIVELY. ALL OF THESE NUMBERS WERE CHECKED.
! THE QUANTITIES Q1(J), Q2(J), ...,Q5(J) ARE THE TERMS IN EQUATION 3.14 OF
! NEUTRINO ASTROPHYSICS AND IN EQUATION 53 OF FOWLER, CAUGHLAN, AND
! ZIMMERMAN (1967) EQ. 53, IN BOTH CASES MULTIPLIED BY T SUB 9 ^(-2/3).
! THE REACTIONS CORRESPONDING TO EACH J ARE LISTED ABOVE, UNDER:
! IDENTIFYING THE RECTIONS.
! FOR THIS SET OF PARAMETERS, AND ONLY FOR THIS SET OF PARAMETERS,
! J = 8 CORRESPONDS TO THE BE7 + P REACTION.
! GENERAL EXPRESSION FOR Q'S:
! T_9^(-2/3)[S_EFF/S(0)] =
! [T_9^(-2/3) + Q1(I)T_9^(-1/3) + Q2(I) + Q3(I)T_9^(1/3) +
! Q4(I)T_9^(2/3) + Q(5)T_9 ]
! BY COMPARISON WITH EQUATION 3.14 OF NEUTRINO ASTROPHYSICS, WE SEE THAT
! EACH OF THE Q'S IS INDEPENDENT OF TEMPERATURE (T), AS CAN BE SEEN FROM
! EQUATIONS 3.10 AND 3.11 .
! ALL OF THE VALUES OF THE Q1, ...., Q5 HAVE BEEN RECALCULATED, USING
! WHERE NEEDED NUCLEAR CROSS SECTIONS GIVEN IN TABLES 3.2 AND 3.4 OF NEUTRINO
! ASTROPHYSICS.
! Q6 IS THE COEFFICIENT OF THE TEMPERATURE TERM IN THE DEFINITION OF
! TAU, EQUATION 3.10 OF NEUTRINO ASTROPHYSICS.
! TAU = Q6*(T SUB 9 TO THE (-1/3) POWER ).
! SLIGHT CHANGES HAVE BEEN MADE IN THE PREVIOUS VALUES OF Q6 TO MAKE
! THE DATA MORE ACCURATE.
! NOTE Q6 IS NEGATIVE.
! Q7 IS THE CONSTANT IN FRONT OF THE REACTION RATE. THE LOGARITHMS ARE
! RELATED TO S FACTORS BY A LOGARITHM AND VARIOUS NUMERICAL FACTORS.
! THE GENERAL RELATION IS:
! Q7 = 70.62860 + ((LN(Z0*Z1/A))/3) -LN(A0*A1) - LN(S SUB 0)
! -LN(1 + DELTA_01)
! HERE S SUB 0 IS THE CROSS SECTION FACTOR IN UNITS OF KEV-BARNS.
! THE NUMERICAL VALUES USED HERE ARE TAKEN FROM TABLES 3.4 AND 3.2
! OF NUCLEAR ASTROPHYSICS.
! THE QUANTITY DELTA_01 IS NON-ZERO (EQUAL TO UNITY) ONLY WHEN THE
! TWO REACTING NUCLEI, 0 AND 1, ARE IDENTICAL.
! Q8 REFLECTS A TERM IN THE EXPONENTIAL THAT OCCURS IN THE RATES, THE
! TERM BEING PROPORTIONAL TO E^( CONSTANT*T_9^2). SEE HARRIS, ET AL.
! 1983, ANNUAL REV. ASTRON. ASTROPHYS.21, 165 (1983) FOR THE MEANING
! OF THIS OBSCURE TERM.
! THE VALUES OF Q1(I),...Q7(I) GIVEN BELOW WERE OBTAINED USING THE DATA
! IN NEUTRINO ASTROPHYSICS WITH THE AID OF AUXILARY COMPUTER CODES
! THAT GENERATED ACCURATE EVALUATIONS. THE NUMBERS GIVEN HERE ARE
! IN MANY CASES COMPLETELY DIFFERENT FROM THE VALUES IN THE ORIGINAL
! YALE CODE.
! DBG 8/94 APPLIED MHP UPDATE TO NUCLEAR REACTIONS
! 10/13/97. Changed Q(I) so that are now calculated for the bare nuclear
! masses. The calculates were made cues.f . Last date on the calculations
! in my notes is 10/10/97.  Note Q6 is -tau(T_9 = 1); in cues.f, I calculate
! tau(T_6 = 1). The connection is:  Q6 = 0.1*tau(cues.f).
! BP00 values
! DATA Q1/0.12317,.03392,.0325,.0304,.03035,.0273,.02494,.040572/,
! 1 Q2/1.08749,-.273,-.2085,.7630,-0.4044,-1.60,-1.224, -0.2095/,
! 2 Q3/.93833,-.0648,-.0474,.1626,-.08598,-.3064,-.2139,-0.0595/,
! 3 Q4/0.,0.,0.,4.79,7.456,0.0,.69703,.16762/,
! 4 Q5/0.,0.,0.,2.595,4.032,0.0,0.3097,.12114/,
! 5 Q6/-3.3804,-12.2757,-12.826,-13.6899,-13.7173,-15.2281,
! 5    -16.6925/,
! 6 Q7/20.8964,76.6003,67.8036,69.130,70.3809,69.8517,70.8012/,
! 7 Q8/0.,0.,0.,0.0,0.0,0.0,0.0/,
! 8NRXNS/13/
! Revised values 7/21/03
! DATA Q1/0.12317,.03392,.0325,.0304,.03035,.0273,.02494,.040572/,
! 1Q2/1.08749,-.273,-.2085,.7630,-0.4044,-1.60,-1.224, 0.0/,
! 2Q3/.93833,-.0648,-.0474,.1626,-.08598,-.3064,-.2139,0.0/,
! 3Q4/0.,0.,0.,4.79,7.456,0.0,.69703,0.0/,
! 4Q5/0.,0.,0.,2.595,4.032,0.0,0.3097,0.0/,
! 5Q6/-3.3804,-12.2757,-12.826,-13.6899,-13.7173,-15.2281,-16.6925/,
! 6Q7/20.8964,76.6003,67.8036,69.130,70.3809,69.8517,70.8012/,
! 7Q8/0.,0.,0.,0.0,0.0,0.0,0./,
! 8NRXNS/13/
! MHP 9/14 RESTORED BE7+P DERIVATIVES THAT WERE ZEROED OUT
! For different values of SSTANDARD, check the comments in ENGEB
! THE VALUES OF SSTANDARD(I) ARE TO BE CHANGED FROM UNITY IF THE CROSS
! SECTION FACTORS ARE NOT THE ONES GIVEN IN NEUTRINO ASTROPHYSICS,
! TABLE 3.2 AND TABLE 3.4 .
! THE CURRENT VALUES ARE CHANGED TO BE THE PREFERRED VALUES LISTED IN
! THE LAST COLUMN OF TABLE 1 OF BAHCALL AND PINSONNEAULT (1991),
! 6/14/91.
! THE VALUE OF SSTANDARD(17) FOR HEP IS TAKEN FROM WOLFS ET AL. , PHYS.
! REV. LETTERS, 63, 2721 (1989).  IT CORRESPONDS TO AN S SUB 0 =
! 15.3E-20 KEV-BARNS, A FACTOR OF 1.913 LARGER THAN FOR THE OLDER
! MEASUREMENTS USED IN NEUTRINO ASTROPHYSICS.
! ZPRD IS USED IN SCREENING CALCULATIONS. IT IS THE PRODUCT OF THE
! CHARGES OF THE INTERACTING IONS. ZPRD WAS CHECKED. Z86 IS USED
! IN CALCULATING INTERMEDIATE SCREENING AND IS DEFINED BY GRABOSKE ET
! AL, AP. J. 181, PAGE 465 (1973), IN TABLE 4.  Z86 WAS CHECKED AND
! SOME NUMERICAL VALUES WERE MADE SLIGHTLY MORE ACCURATE.  Z53, Z43,
! AND Z23 ARE ALSO DEFINED IN TABLE 4 (SEE ABOVE).  SINCE THEY ARE
! ONLY USED IN STRONG SCREENING, THE VALUES OF Z53, Z43, AND Z23
! WERE NOT CHECKED.
! DEFINE NEXT THE FRACTIONAL ABUNDANCES BY MASS OF THE IMPORTANT
! ISOTOPES.
! X, Y, Z, XHE3,..., XBE9 ARE THE MASS FRACTIONS OF THE ISOTOPES.
! THE ABUNDANCES OF NEUTRONS, H2, H3, NE20,AND MG24, WHICH ARE,
! RESPECTIVELY, XFRAC(I) FOR I = 1,3,4,12,13, ARE NO LONGER USED.
! BEGIN CALCULATION OF SCREENING CORRECTION.
! THE BASIC REFERENCES ARE SALPETER, AUSTRALIAN JOURNAL OF PHYSICS,
! VOL. 7, 373 (1954). THE FORMULA FOR WEAK SCREENING THAT IS BEING
! PROGRAMMED IS EQUATION (25) OF THIS PAPER. THE OTHER IMPORTANT
! REFERENCES ARE: DEWITT, GRABOSKE, AND COOPER, AP. J. 181, 439 (1973)
! AND GRABOSKE ET AL., AP. J. 181, 457 (1973). THE VALUES OF EMU AND
! ZET ARE ESSENTIAL FOR COMPUTING WEAK SCREENING; THE VALUE OF AMU IS
! USED IN AN NON-ESSENTIAL WAY IN THIS COMPUTATION. XTR IS USED IN
! COMPUTING INTERMEDIATE SCREENING.
! AMU IS ONE OVER THE MEAN MOLECULAR WEIGHT OF THE IONS, MU SUB I .
! EMU IS ONE OVER THE ELECTRON MEAN MOLECULAR WEIGHT, MU SUB E.
! EMU IS USED HERE AS THE NAME FOR THE SECOND PART OF THE ZETA FUNCTION
! IN THE SALPETER EXPRESSION FOR WEAK SCREENING.
! XTR IS USED LATER IN THE INTERMEDIATE SCREENING CALCULATION.  THE AVERAGE
! OF THE QUANTITY Z**(3B -1) IS EQUAL TO XTR/AMU.
! ZET IS THE FIRST PART OF THE SALPETER SCREENING ZETA VARIABLE.
! MU = SUM OVER I OF [X(I)/A(I)].
! MU SUB E = SUM OVER I [ Z(I)*X(I)/A(I)].
! DL AND DT ARE THE THE LOG10 OF THE DENSITY AND TEMPERATURE.
! THE UNIT OF TEMPERATURE IS 10^9 K AND THE UNIT OF DENSITY IS
! GM PER CM^3 .
! DD = LOG RHO TO THE BASE 10.
! CLN = LN10.  CLN IS CONVERSION BETWEEN LOG10 AND LN.
! CONVERT DENSITY TO UNLOGGED FORM.
! RWE = RHO/(MU SUB E), I. E., THE NUMBER OF ELECTRONS DIVIDED BY
! AVOGADRO'S NUMBER.
! THE EXPRESSION FOR RWE WAS INCORRECT IN THE ORIGINAL YALE SUBROUTINE.
! THE ORIGINAL VERSION HAD ( EXP(CLN*DL) ) DIVIDED BY EMU INSTEAD OF
! MULTIPLIED BY EMU.  RWE IS USED LATER IN COMPUTING THE SCREENING
! CORRECTION.
! SET RATES EQUAL TO ZERO FOR THE LOG_10(T) < TCUT(1)
! T9P13 IS THE TEMPERATURE IN UNITS OF 10^9 DEGREES K TO THE PLUS 1/3
! POWER.  MINUS IS DENOTED BY M.  HERE T9 IS THE TEMPERATURE IN UNITS
! OF 10^9 K, CONVERTED FROM THE LOG_10 (T) AND RHO IS THE DENSITY IN
! CGS UNITS.
! F PRIME/F
! THE NEXT PIECE OF CODE COMPUTES FIRST THE FERMI ENERGY DIVIDED BY KT, WHERE
! PFMC2 IS THE FERMI MOMENTUM DIVIDED BY MASS OF ELECTRON TIMES C, ALL
! SQUARED AND EFMKT IS THE FERMI ENERGY DIVIDED BY KT.  THE QUANTITY
! FPRF IS THE RATIO OF F PRIME TO F IN SALPETER'S SCREENING CORRECTION.
! THE VALUE OF FPRF IS DETERMINED IN THE INTERMEDIATE CASE BY AN
! INTERPOLATION FORMULA DEPENDING UPON THE DEGREE OF DEGENERACY, AS
! MEASURED BY DEGD = LOG10(E_F/KT). THE NUMERICAL VALUES FOR THE FIT
! WERE TAKEN FROM SALPETER'S ORIGINAL PAPER, FIGURE 1. THE ONLY CHANGES
! IN THIS PART OF THE SUBROUTINE WERE THE CORRECTION OF THE ERROR IN
! THE DEFINITION OF RWE (SEE ABOVE) AND REFINEMENTS OF THE COEFFICIENTS
! IN THE EXPRESSIONS FOR PFMC2 AND EFMKT.
! END OF CALCULATION OF FERMI ENERGY DIVIDED BY KT AND THE INTERPOLATION
! FORMULA FOR (F PRIME/F), WHICH APPEARS IN SALPETER'S SCREENING
! FORMULA.
! NOW WE GET TO THE COMPUTATION OF WEAK SCREENING (SEE ALSO PAGE 61 OF
! NEUTRINO ASTROPHYSICS, WHICH GIVES ONLY A SIMPLIFIED FORMULA) AND
! ALSO THE MORE COMPLICATED INTERMEDIATE AND STRONG SCREENING CASES
! (SEE REFERENCES TO AP. J. 181 ABOVE FOR INTERMEDIATE
! AND STRONG SCREENING), ESPECIALLY PAGE 465.
! XXL IS USED IN ALL THE SCREENING FORMULAE.  THIS QUANTITY IS THE
! IS THE FUNCTION CALLED LAMBDA SUB ZERO BY GRABOSKE ET AL. THE
! QUANTITY XXL**0.86 = XXL8 IS USED IN CALCULATING INTERMEDIATE
! SCREENING.
! ZCURL IS THE QUANTITY DEFINED BY EQUATION (4) OF DEWITT ET AL; IT
! IS THEIR Z WITH A CURLY SYMBOL ON ITS TOP. ZCURL IS USED IN WEAK
! AND IN INTERMEDIATE SCREENING AND WAS FIRST DEFINED BY SALPETER.
! ZCURL IS THE SAME AS SALPETER'S ZETA EXCEPT FOR THE FACTOR OF 1/AMU.
! (Z SUB 1 TIMES Z SUB 2)*XXL*ZCURL GIVES THE WEAK SCREENING FACTOR,
! THE SAME AS SALPETER OR AS EQUATION (19) OF GRABOSKE ET AL.
! Z BAR IS THE SAME AS Z BAR OF DE WITT ET AL AND IS THE AVERAGE CHARGE
! OF THE IONS.  IT IS EQUAL TO EMU/AMU.  THE QUANTITY (Z BAR)**0.28 =
! Z28 OCCURS IN THE COMPUTATION OF INTERMEDIATE SCREENING.
! XXL6 IS USED FOR COMPUTING STRONG STRONG SCREENING.
! THE NOTATION USED HERE IS EXPLAINED IN LARGE PART BY
! EQUATION (4) OF DEWITT ET AL., AP. J. 181, PAGE 439.
! THE FINAL EXPRESSION FOR WEAK SCREENING IS EXACTLY EQUAL TO SALPETER'S
! FORMULA, WHICH INCLUDES A DEGENERACY CORRECTION. THE MORE GENERAL
! EXPRESSIONS ARE GIVEN IN TABLE 4 AND EQUATION (19) OF GRABOSKE ET AL.
! COMPUTE SCREENING FOR EACH OF THE REACTIONS.
! star%ctrl%weak_screening_threshold (WEAKSCREENING) IS A NUMERICAL
! PARAMETER READ FROM THE INPUT. TO OBTAIN THE GRABOSKE ET AL. AND SALPETER STANDARD RESULTS,
! USE: WEAKSCREENING = 0.03.  FOR THE STANDARD SOLAR MODEL, THIS IS THE
! VALUE THAT SHOULD BE ADOPTED. TO INVESTIGATE THE EFFECT OF ALWAYS USING
! WEAK SCREENING, USE A LARGE VALUE FOR WEAKSCREENING, E. G., 30.  AS
! LONG AS WEAKSCREENING IS ASSUMED TO BE BIGGER THAN ONE, THE PROGRAM
! WILL ALWAYS CALCULATE FOR THE SUN WITH THE WEAK SCREENING
! APPROXIMATION.
! UTOT IS THE FINAL SCREENING CORRECTION WHICH APPEARS IN THE
! RATE EXPRESSION AS: EXP(UTOT) .
! END OF SCREENING CALCULATION. WEAK AND INTERMEDIATE SCREENING FORMS
! ARE GIVEN CORRECTLY.  STRONG SCREENING WAS NOT CHECKED BECAUSE IT IS
! NOT RELEVANT FOR THE SUN.
! Reactions first_zeroed_rxn..num_reactions are zeroed after this block:
! only the alpha captures when hydrogen burns, all 13 when it is gone.
! CALCULATE REACTION RATES FOR THE THREE PRINCIPAL REACTIONS OF
! THE PP CHAIN: PP, HE3+HE3, HE3 +HE4, AND THE FOUR PROTON
! BURNING REACTIONS ON C12: C13, N14, AND O16.
! R1 IS (T SUB 9)^(-3/2) TIMES (S SUB EFF)/(S SUB 0). THE CORRECT
! EXPRESSION FOR (S SUB EFF)/(S SUB 0) IS GIVEN IN EQUATION 3.14 OF
! NEUTRINO ASTROPHYSICS.  THE NUMERICAL FORM THAT IS USED IS EQUATION
! 52 OF FOWLER, CAUHLAN, AND ZIMMERMAN, VOL. 5, 1967.
! RATE(I) IS THE RATE OF THE DIFFERENT REACTIONS PER SECOND PER GRAM,
! EXCEPT THAT THE MASS FRACTIONS ARE OMITTED.
! THE PREVIOUS YALE VERSION HAD THE RATE FOR THE O16 + P REACTION
! MULTIPLIED BY T9**(-1/7). THIS FACTOR IS INCORRECT AND HAS BEEN
! REMOVED; IT APPEARED BEFORE AS AN IF STATEMENT REFERRING ONLY TO
! RATE(7).
! R1=T9M23+Q1(I)*T9M13+Q2(I)+Q3(I)*T9P13+Q4(I)*T9P23+Q5(I)*T9
! MHP 8/14 RATES CORRECTED TO PERMIT USER MODIFICATION OF REACTION
! RATE DERIVATIVES
! END OF CALCULATION OF REACTION RATES FOR FIRST 7 REACTIONS.
! BE7 BURNING
! CALCULATE THE BURNING OF BE7 BY PROTONS (WHICH PRODUCES THE MOST
! EXPERIMENTALLY ACCESSIBLE SOLAR NEUTRINOS) AND THE BURNING OF BE7
! BY ELECTRON CAPTURE (THE DOMINANT PROCESS). THE PREVIOUSLY USED
! EXPRESSION CONTAINED A TERM THAT WAS PHYSICALLY INCORRECT
! AND ALSO USED CROSS SECTION FACTORS THAT WERE OUT OF DATE.
! THE ELECTRON CAPTURE RATE IN SEC^(-1) IS GIVEN BY EQUATION (3.18) OF
! NEUTRINO ASTROPHYSICS.  CAN OMIT THE FACTOR OF RHO (IN CGS UNITS)
! WHICH ALSO APPEARS IN THE BE7PROTON RATE. THE BE7 PROTON CAPTURE
! RATE IS GIVEN BY TABLE 3.2 AND EQUATION (3.12).  NOTE THAT 1 OVER MU
! SUB E IS EQUAL TO EMU.
! HERE WE USE THE NOTATION BE7 + PROTON = BE7PROTON AND
! BE7 + E = BE7ELECTRON.
! F1 IS THE FRACTION OF THE BE7 THAT IS BURNED BY ELECTRON CAPTURE.
! 1-F1 IS THE FRACTION OF THE BE7 THAT IS BURNED BY PROTON CAPTURE.
! F3 IS THE FRACTION OF THE N14 THAT IS BURNED BY P, ALPHA REACTION.
! 1-F3 IS THE FRACTION OF THE N15 THAT IS BURNED BY P, GAMMA REACTION.
! SEE TABLE 21 OF BAHCALL AND ULRICH (1988), REV. MOD. PHYS. 60.
! 10/13/97. I changed Temp3 (i.e., tau) by 5/10^5 as a result of
! using bare nuclear masses. Previously coefficient was -10.26202.
! I did not change Be7electron for the slightly different S0 for
! electron capture, since that is done in SStandard. Previous coefficient
! in Be7electron expression was (3.126571E+5). 10/14/97.
! INCLUDE FOR BE7PROTON THE T9M23 FACTOR AND ALL CORRECTIONS PROPORTIONAL TO
! Q1,...,Q5 FROM EQUATION 3.14 OF NEUTRINO ASTROPHYSICS. THESE
! CORRECTIONS ARE DEFINED EARLIER IN THIS SUBROUTINE.
! QRBE7 = T9M23 + Q1(8)*T9M13 + Q2(8)+ Q3(8)*T9P13
! $          + Q4(8)*T9P23 + Q5(8)*T9
! MHP 9/14 ADDED THE ABILITY TO ALTER DERIVATIVES INDEPENDENTLY
! CALCULATE THE SCREENING CORRECTION FOR BE7 + P REACTION.  USE WEAK AND
! INTERMEDIATE SCREENING FORMULAE.
! END OF CALCULATION OF SCREENING CORRECTION FOR BE7 + P REACTION.
! MULTIPLY RATES BY FACTOR OF RHO/[(ATOMIC MASS UNIT)*A(BE7)] TO GET
! IN UNITS OF GM^{-1}.  CALL FACTOR CAMUBE7. Corrected for Be7 bare
! mass number. Previously used 8.582295E+22 in the expression below.
! These corrections were made on 10/14/97.
! END OF MULTIPLICATION INSERTED NOVEMBER 6, 1990.
! END OF CALCULATION OF CRUCIAL BE7 ELECTRON CAPTURE AND PROTON CAPTURE
! RATES AND THEIR RATIO.
! N15 + P BRANCHING.
! THE FOLLOWING STATEMENTS COMPUTE THE BRANCHING OF N15(P,ALPHA)C12
! AND N15(P,GAMMA)O16. THESE STATEMENTS REPLACE OUTDATED STATEMENTS
! IN THE YALE CODE. THE CNO CROSS-SECTION FACTORS ARE FROM TABLE 3.4
! OF NEUTRINO ASTROPHYSICS.  THE RATIO OF THE REACTIONS DEPENDS ONLY
! UPON THE EFFECTIVE ZERO ENERGY S-FACTOR, WHICH IS S0(ZERO ENERGY)
! TIMES THE COMBINATION OF TEMPERATURE AND S-FACTOR DERIVATIVES/S0
! THAT WAS USED PREVIOUSLY AS R1 IN THE RATE CALCULATIONS. THE
! NUMERICAL COEFFICIENTS THAT APPEAR IN THE RATE WERE REPRESENTED
! BY THE Q1(J),...Q5(J) FOR THE OTHER REACTIONS.
! THE QVALUES FOR THE N15 REACTIONS HAVE BEEN COMPUTED SEPARATELY.
! DBG 8/94 APPLIED MHP UPDATE TO NUCLEAR REACTIONS
! To agree with Solar Fusion Workshop paper, the value of S0 in keV-b
! has been changed from 78000. to 675000. on 9/25/97. On 10/14/97,
! JNB changed cues.f so as to compute the Q-coefficients for the N15 + p
! reactions.  Also, checked that the coefficients are the same as the
! ones given earlier in energy.f when we use the older CNO data.
! MULTIPLY BY THE VALUE OF S0 IN KEV-B.
! MHP 9/14 ADDED THE OPTION TO MODIFY THE RELATIVE CROSS SECTIONS
! FOR N15+P -> C12+ALPHA AND O16+GAMMA
! END OF NEW ROUTINE FOR THE BRANCHING OF N15 + P .
! ***MHP 3/91 ALPHA CAPTURE REACTIONS UPDATED TO CAUGHLAN AND FOWLER(1988)
! RATES.  THE RATES ARE EXPRESSED IN THE SAME TERMS USED BY CZ, WITH
! THE CONVERSION FACTOR IN THE FRONT OBTAINED FROM VANDENBERG'S
! NOTES ON THE REACTION RATES.
! RATE(8)  HE4+C13
! RATE(10) HE4+C12=>O16
! RATE(11) HE4+N14=>O18
! RATE(12) TRIPLE ALPHA
! C13(ALPHA,N) O16
! C12(ALPHA,GAMMA)O16
! N14(ALPHA,GAMMA)F18 + F18=>O18+EPLUS+NU
! 2026 (bugsweep Batch 3): the linear coefficient of the CF88
! N14(a,g)F18 polynomial is 0.117 (1 + 0.012 T9^1/3 + 1.45 T9^2/3
! + 0.117 T9 + 1.97 T9^4/3 + 0.406 T9^5/3); both copies of the rate
! carried 0.177 -- a digit transposition inherited from the F77
! source (~0.5% of the rate at T9 = 0.15).
! TRIPLE ALPHA
! END OF XEROING OUT OF REACTIONS 9 AND 13.
! RATES PER 10^9 YEARS PER ATOMIC MASS UNIT: HRK(IU)
! HR1, ..., HR13 ARE THE RATES OF THE INDIVIDUAL REACTIONS.
! THE INTERPRETATION OF WHICH REACTION GOES WITH WHICH SYMBOL CAN BE
! MADE EASILY BY LOOKING AT THE DEFINITIONS OF THE EG(I)'S.
! THE ABUNDANCES ARE UPDATED IN SUBROUTINE KEMCOM USING THESE MATRICES.
! gyr_amu_per_sec_gram (C21) IS THE PRODUCT OF (10^9 YEARS/1 SECOND)*
! (1 ATOMIC MASS UNIT/1 GRAM). I HAVE USED HERE SIDEREAL YEAR IN
! CONVERTING TO SECONDS.
! rows rr_zero9/rr_zero13: zero-filled placeholders (rate(9), rate(13)
! are always 0 -- see the XEROING block above).
! END OF COMPUTATION OF HRK(IU).
subroutine nuclear_reaction_rates(log_density, log_temperature, mass_frac, nr)
      use star_info_lib, only: star
      use phys_const_lib
      use math_lib
      implicit none
      double precision, intent(in) :: log_density, log_temperature, &
           mass_frac(num_isotopes)
      type(nuclear_rates), intent(out) :: nr

      double precision :: screening_factor(num_reactions), &
           dscreen_dlnrho(num_reactions), dscreen_dlnt(num_reactions)
      double precision :: term, ion_mean_weight_inverse, &
           electron_mean_weight_inverse, xtr, zeta_sum, hydrogen_fraction
      double precision :: density, t9, t9_p13, t9_p23, t9_m13, t9_m23, &
           t9_m1, t9_m2, t9_m12, t9_m32
      double precision :: pfmc2, efmkt, fprf, degd
      double precision :: xxl, xxl6, xxl8, zcurl, zbar, z58, z28, z33, tm1
      double precision :: uwk, uint, ustr
      double precision :: r1, r2, a1, a2, a3, a4, a5, dr1, da1, s_sum
      double precision :: temp3, qrbe7, camube7, o16gamma, c12alpha
      integer :: i, first_zeroed_rxn

! composition moments for the screening
      hydrogen_fraction = mass_frac(iso_h1)
      ion_mean_weight_inverse = 0.d0
      electron_mean_weight_inverse = 0.d0
      xtr = 0.d0
      zeta_sum = 0.d0
      do i = 1,num_isotopes
         term = mass_frac(i)/isotope_mass(i)
         ion_mean_weight_inverse = ion_mean_weight_inverse+term
         electron_mean_weight_inverse = electron_mean_weight_inverse+ &
              term*isotope_charge(i)
         xtr = xtr+term*pow(isotope_charge(i), 1.58d0)
         zeta_sum = zeta_sum+term*isotope_charge(i)**2
      end do
      nr%electron_density_na = ( exp(ln10*log_density) )* &
           electron_mean_weight_inverse
      if (log_temperature.le.star%ctrl%nuclear_logT_cutoffs(1)) return

      density=exp(ln10*log_density)
      t9 = exp(ln10*(log_temperature - 9.0d0))
      t9_p13 = pow(t9, cc13)
      t9_p23 = t9_p13**2
      t9_m13=1.d0/t9_p13
      t9_m23=t9_m13**2
      t9_m1=1.d0/t9
      t9_m2=t9_m1**2
      t9_m12=1.d0/dsqrt(t9)
      t9_m32=t9_m1*t9_m12
      nr%density = density
      nr%t9 = t9
      nr%t9_p13 = t9_p13
      nr%t9_p23 = t9_p23
      nr%t9_m13 = t9_m13
      nr%t9_m23 = t9_m23
      nr%t9_m12 = t9_m12

! screening (Salpeter weak; Graboske et al. intermediate and strong),
! with its logarithmic derivatives
      pfmc2=1.017677d-4*pow(nr%electron_density_na, (2.0d0/3.0d0))
      efmkt=5.92986d0*t9_m1*(dsqrt(1.d0+pfmc2)-1.d0)
      if (efmkt.le.1.d-2) then
         fprf=1.0d0
      else
         degd=log10(efmkt)
         if (degd.ge.1.5d0) then
            fprf=0.0d0
         else
            fprf=0.75793d0-0.54621d0*degd-0.30964d0*degd**2+0.12535d0*degd**3+ &
            0.1203d0*degd**4-0.012857d0*degd**5-0.014768d0*degd**6
         end if
      end if
      xxl=5.9426d-6*t9_m32*dsqrt(density*ion_mean_weight_inverse)
      xxl6=pow(xxl, (2.0d0/3.0d0))
      xxl8=pow(xxl, 0.86d0)
      zcurl=dsqrt((zeta_sum+fprf*electron_mean_weight_inverse)/ &
           ion_mean_weight_inverse)
      zbar=electron_mean_weight_inverse/ion_mean_weight_inverse
      z58=pow(zcurl, 0.58d0)
      z28=pow(zbar, 0.28d0)
      z33=pow(zbar, cc13)
      tm1=xxl*zcurl
      do i=1,num_reactions
         uwk=tm1*screen_zprod(i)
         if (uwk.le.star%ctrl%weak_screening_threshold) then
            screening_factor(i)=uwk
            dscreen_dlnrho(i)=0.5d0*uwk
            dscreen_dlnt(i)=-1.5d0*uwk
         else
            uint=0.38d0*xxl8*xtr*screen_z86(i)/(ion_mean_weight_inverse*z58*z28)
            if (uwk.le.2.d0) then
               screening_factor(i)=uint
               dscreen_dlnrho(i)=0.43d0*uint
               dscreen_dlnt(i)=-1.29d0*uint
            else
               ustr=0.624d0*z33*xxl6*(screen_z53(i)+0.316d0*z33*screen_z43(i)+0.737d0* &
               screen_z23(i)/(zbar*xxl6))
               if (ustr.lt.uint.or.uwk.ge.5.d0) then
                  screening_factor(i)=ustr
                  dscreen_dlnrho(i)=0.208d0*z33*(screen_z53(i)+0.316d0*z33*screen_z43(i))*xxl6
                  dscreen_dlnt(i)=-3.d0*dscreen_dlnrho(i)
               else
                  screening_factor(i)=uint
                  dscreen_dlnrho(i)=0.43d0*uint
                  dscreen_dlnt(i)=-1.29d0*uint
               end if
            end if
         end if
      end do

! reactions 1-7 (S-factor form) and the two branchings, with hydrogen
      first_zeroed_rxn = r_pp
      if (hydrogen_fraction.ne.0.0d0) then
         first_zeroed_rxn = r_ac13
         do i=r_pp,r_po16
            r1=t9_m23+fit_q1(i)*t9_m13+star%qs0e_scale(i)*(fit_q2(i)+fit_q3(i)*t9_p13)+ &
                 star%qqs0ee_scale(i)*(fit_q4(i)*t9_p23+fit_q5(i)*t9)
            nr%rate(i)=density*r1*exp(fit_q6(i)*t9_m13+fit_q7(i)+(fit_q8(i)*t9)**2+ &
                 screening_factor(i))
            nr%rate(i) = nr%rate(i)*star%cross_section_scale(i)
            nr%dlnrate_dlnrho(i)=1.d0+dscreen_dlnrho(i)
! the exact derivative of the rate above, scales included (2026-10:
! the q3 term lacked qs0e_scale and the q4, q5 terms qqs0ee_scale, so
! once those scales were ~1 instead of ~1e-22 it no longer matched)
            nr%dlnrate_dlnt(i)=dscreen_dlnt(i)-(fit_q6(i)*t9_m13+(2.d0*t9_m23+ &
            fit_q1(i)*t9_m13-star%qs0e_scale(i)*fit_q3(i)*t9_p13- &
            star%qqs0ee_scale(i)*(2.d0*fit_q4(i)*t9_p23+3.d0*fit_q5(i)*t9))/r1)/3.d0+ &
            2.d0*(fit_q8(i)*t9)**2
         end do
! Be7: electron capture vs proton capture
         nr%be7_electron = (1.752d-10)*t9_m12*(1.0d0 + 0.004d0*(1000.d0*t9 - 16.d0))
         nr%be7_electron = nr%be7_electron*electron_mean_weight_inverse*star%cross_section_scale(s_be7e)
         temp3 = (-10.2625d0*t9_m13)
         nr%be7_proton = (3.128813d+5)*hydrogen_fraction*star%cross_section_scale(s_be7p)*exp(temp3)
         qrbe7 = t9_m23 + fit_q1(iq_be7p)*t9_m13 + star%qs0e_scale(iq_be7p)*(fit_q2(iq_be7p)+ fit_q3(iq_be7p)*t9_p13) &
              + star%qqs0ee_scale(iq_be7p)*(fit_q4(iq_be7p)*t9_p23 + fit_q5(iq_be7p)*t9)
         nr%be7_proton = nr%be7_proton*qrbe7
         nr%be7_proton = nr%be7_proton*exp(weak_intermediate_screening(4.0d0, 5.7790d0))
         camube7 = density*8.584981d+22
         nr%be7_proton = camube7*nr%be7_proton
         nr%be7_electron = camube7*nr%be7_electron
         nr%frac_be7_ecap = nr%be7_electron/(nr%be7_electron + nr%be7_proton)
         nr%frac_be7_pcap = nr%be7_proton/(nr%be7_electron + nr%be7_proton)
! N15: (p,alpha)C12 vs (p,gamma)O16
         o16gamma = t9_m23 + 0.0273016d0*t9_m13 + 0.14374d0 + 0.027490d0*t9_p13 &
                   + 6.14685d0*t9_p23 + 2.98940d0*t9
         o16gamma = o16gamma*s0_n15pg_kevb*star%o16_gamma_scale
         c12alpha = t9_m23 + 0.0273016d0*t9_m13 + 2.01186d0 + 0.384763d0*t9_p13 &
                   + 17.0579d0*t9_p23 + 8.29580d0*t9
         c12alpha = c12alpha*s0_n15pa_kevb*star%c12_alpha_scale
         nr%frac_n15_pa = c12alpha/(c12alpha + o16gamma)
         nr%frac_n15_pg = 1.0d0 - nr%frac_n15_pa
      end if
      do i=first_zeroed_rxn,num_reactions
         nr%rate(i)=0.d0
         nr%dlnrate_dlnrho(i)=0.d0
         nr%dlnrate_dlnt(i)=0.d0
      end do

! alpha captures (Caughlan & Fowler 1988) and triple alpha
      if (log_temperature.ge.star%ctrl%nuclear_logT_cutoffs(4)) then
! C13(alpha,n)O16
      r1=t9_m23+0.0129d0*t9_m13+2.04d0+0.184d0*t9_p13
      a1 = 6.77d15*exp(-32.329d0*t9_m13-(t9/1.284d0)**2)
      a2 = 3.82d5*exp(-9.373d0*t9_m1)
      a3 = 1.41d6*exp(-11.873d0*t9_m1)
      a4 = 2.00d9*exp(-20.409d0*t9_m1)
      a5 = 2.92d9*exp(-29.283d0*t9_m1)
      s_sum = a1*r1+t9_m32*(a2+a3+a4+a5)
      nr%rate(r_ac13) = 1.157126d22*density*exp(screening_factor(r_ac13))*s_sum
      nr%dlnrate_dlnrho(r_ac13)=1.0d0+dscreen_dlnrho(r_ac13)
      dr1 = cc13*(-2.0d0*t9_m23-0.0129d0*t9_m13+0.184d0*t9_p13)
      da1 = a1*(cc13*32.329d0*t9_m13 - 2.0d0*(t9/1.284d0)**2)
      nr%dlnrate_dlnt(r_ac13) = dscreen_dlnt(r_ac13)+(dr1*a1 + r1*da1 + &
           t9_m32*(a2*(9.373d0*t9_m1-1.5d0)+a3*(11.873d0*t9_m1-1.5d0)+ &
           a4*(20.409d0*t9_m1-1.5d0)+a5*(29.283d0*t9_m1-1.5d0)))/s_sum
! C12(alpha,gamma)O16
      r1 = 1.0d0/(1.0d0+0.0489d0*t9_m23)
      r2 = 1.0d0/(1.0d0+0.2654d0*t9_m23)
      a1 = t9_m2*exp(-32.120d0*t9_m13)
      a2 = 1.04d8*r1**2*exp(-(t9/3.496d0)**2)
      a3 = 1.76d8*r2**2
      a4 = 1.25d3*t9_m32*exp(-27.499d0*t9_m1)
      a5 = 1.43d-2*t9**5*exp(-15.541d0*t9_m1)
      s_sum = a1*(a2+a3)+a4+a5
      nr%rate(r_ac12) = 1.25388d22*density*exp(screening_factor(r_ac12))*s_sum
      nr%dlnrate_dlnrho(r_ac12) = 1.0d0+dscreen_dlnrho(r_ac12)
      nr%dlnrate_dlnt(r_ac12) = dscreen_dlnt(r_ac12)+ &
           (a1*((cc13*32.120d0*t9_m13-2.0d0)*(a2+a3)+ &
           a2*(cc13*0.1956d0*t9_m23*r1-2.0d0*(t9/3.496d0)**2)+ &
           a3*(cc13*1.0616d0*t9_m23*r2))+ &
           a4*(27.499d0*t9_m1-1.5d0)+a5*(5.0d0+15.541d0*t9_m1))/s_sum
! N14(alpha,gamma)F18(e+ nu)O18; linear coefficient 0.117 (bugsweep
! Batch 3: both former copies carried 0.177)
      r1 = t9_m23+0.012d0*t9_m13+1.45d0+0.117d0*t9_p13+1.97d0*t9_p23 &
           +0.406d0*t9
      a1 = 7.78d9*exp(-36.031d0*t9_m13-(t9/0.881d0)**2)
      a2 = t9_m32*2.36d-10*exp(-2.798d0*t9_m1)
      a3 = t9_m32*2.03d0*exp(-5.054d0*t9_m1)
      a4 = t9_m23*1.15d4*exp(-12.310d0*t9_m1)
      s_sum = a1*r1+a2+a3+a4
      nr%rate(r_an14)= 1.07452d22*density*exp(screening_factor(r_an14))*s_sum
      nr%dlnrate_dlnrho(r_an14)=1.d0+dscreen_dlnrho(r_an14)
      dr1 = cc13*(-2.0d0*t9_m23-0.012d0*t9_m13+0.117d0*t9_p13+ &
            3.94d0*t9_p23)+0.406d0*t9
      da1 = a1*(cc13*36.031d0*t9_m13-2.0d0*(t9/0.881d0)**2)
      nr%dlnrate_dlnt(r_an14) = dscreen_dlnt(r_an14)+(dr1*a1+r1*da1+ &
              a2*(2.798d0*t9_m1-1.5d0)+a3*(5.054d0*t9_m1-1.5d0)+ &
              a4*(12.310d0*t9_m1-cc23))/s_sum
! triple alpha
      nr%rate(r_3alpha) = 1.565315d21*density**2*t9_m1*t9_m2*2.79d-8* &
                 exp(-4.4027d0*t9_m1+screening_factor(r_3alpha))
      nr%dlnrate_dlnrho(r_3alpha) = 2.0d0+dscreen_dlnrho(r_3alpha)
      nr%dlnrate_dlnt(r_3alpha) = -3.0d0+dscreen_dlnt(r_3alpha)+4.4027d0*t9_m1
      end if

! hep, He3(p,e+ nu)He4: flux only (screened, scaled; times X*X(He3) by
! the consumer)
      nr%hep = (1.71724d+11)*density*t9_m23*exp(-6.1399d0*t9_m13)
      nr%hep = (1.0d0 + 0.067862d0*t9_p13)*star%cross_section_scale(s_hep)*nr%hep
      nr%hep = nr%hep*exp(weak_intermediate_screening(2.0d0, 3.08687d0))

contains

! weak screening, or intermediate above the weak threshold (the Be7+p
! and hep treatment): charge product zprod, Graboske z86
double precision function weak_intermediate_screening(zprod, z86)
      double precision, intent(in) :: zprod, z86
      uwk = tm1*zprod
      if (uwk.le.star%ctrl%weak_screening_threshold) then
         weak_intermediate_screening = uwk
      else
         weak_intermediate_screening = 0.38d0*xxl8*xtr*z86/(ion_mean_weight_inverse*z58*z28)
      end if
end function weak_intermediate_screening

end subroutine nuclear_reaction_rates

!----------------------------------------------------------------------
! network_rate_vector (formerly rates)
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original rates.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! JULY 3, 1991 (MHP)
! THIS SUBROUTINE COMPUTES THE NUCLEAR BURNING RATES FOR USE IN KEMCOM.
! IT IS A STRIPPED-DOWN VERSION OF ENGEB; SEE ENGEB FOR DETAILED NOTES
! ON THE REACTIONS.
! ALL PHYSICAL QUANTITIES ARE CGS. THE REACTION RATES ARE IN
!  GM^{-1}S^{-1}.
!  MHP 10/97
!
! On 10/13/97, JNB converted the nuclear masses from neutral nuclear
! masses to bare nuclear masses by subtracting Z(I)*(m_e)*c^2 from
! the neutral nuclear masses. This caused changes in a number of
! places: in ANUC(I), in Q1(I)-Q7(I), in the calculation of the Be7electron
! Be7proton rates, and in the calculation of the N15p branching ratio.
!
! JNB made some purely cosmetic changes on 1/20/96. Revised all of the
! input Q's on 9/23-25/97 to agree with submitted version of Solar
! Fusion Workshop paper. The SStandard are fixed to agree with
! the Workshop paper. JNB recalculated all of the EG(I) to determine
! the best values for the energy generation for all the reactions,
! taking account of my improved calculations of neutrino energy loss.
! The calculations are documented in Vol. 19, 132-141, 1997 of my notes.
! The neutral atom mass differences are taken from Table of Isotopes,
! 8th Ed, 1996 and the neutrino energy losses from Bahcall, Gallium
! solar neutrino experiments, Phys Rev C, in press, 1997.
!
! ALL NUMBERS IN THIS SUBROUTINE HAVE BEEN CHECKED AND REVISED, WHERE
!  NECESSARY, BY JOHN BAHCALL SO THAT THEY AGREE WITH THE MODERN NUMBERS
!  IN NEUTRINO ASTROPHYSICS (1989).
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
! THE ENERGY DERIVATIVES ENTER IN A FORM IN WHICH THEY ARE DIVIDED BY
!  THE ABSOLUTE VALUES OF THE CROSS SECTIONS AT ZERO ENERGY.  THUS IF
!  THE SHAPE OF THE CROSS SECTION EXTRAPOLATION IS UNCHANGED AND ONLY
!  THE INTERCEPT OF S(E) AT ZERO ENERGY IS CHANGED, THEN NO CORRECTION
!  NEED BE MADE FOR THE DERIVATIVES.  THEY ARE AUTOMATICALLY SCALED
!  CORRECTLY.  THE EXACT WAY THAT THE DERIVATIVES ENETER THE
!  CALCUALTIONS IS DESCRIBED IN THE SECTION LABELED ``DEFINING THE
!  Q(I)' THAT IS PRESENTED BELOW.
! ************************************
! IDENTIFYING THE REACTIONS.
! ************************************
!  THE VALUE OF J DENOTES WHICH OF THE REACTIONS THE COEFFICIENTS
!  REFER TO:
!  J = 1, PP; J = 2, HE3+HE3; J = 3, HE3+ HE4; J =4, P + C12;  J = 5, P+C13;
!  J = 6. P + N14; J = 7, P + O16 ; J = 8, HE4+C13; J = 10, HE4+C12;
!  J = 11, HE4 + N14; J = 12, TRIPLE ALPHA.
!  REACTION 14 IS PEP; REACTION 15 IS BE7 ELECTRON CAPTURE; REACTION 16 IS
!   BE7 PROTON CAPTURE; REACTION 17 IS THE HEP REACTION.
!  do not change sstandard(14) unless you want to change the ratio of
!  pep to pp.
!   REACTIONS 14-17 WERE NOT EXPLICITLY INCLUDED IN THE YALE
!    PREVIOUS VERSION OF THE CODE, BUT THEY ARE MOST OF THE STORY FOR
!    THE SOLAR NEUTRINO PROBLEM.
!   THE BRANCHING OF THE N15 + P REACTIONS IS TREATED IN A SERIES OF
!    SEPARATE STATEMENTS FOLLOWING THE CALCULATION OF THE BE7 + P
!    REACTION. SEE THE DEFINITIONS OF F3 AND F4.  IF THE CROSS-SECTION
!    FACTORS OF THE N15 + P REACTIONS ARE REVISED, THEN THE NUMERICAL
!    COEFFICIENTS MUST BE CHANGED IN THE DEFINITION OF C12ALPHA AND
!    O16GAMMA.
! FOR Q1(I), ...,Q(5(I), I = 8 CORRESPONDS TO THE BE7 +  P REACTION.
!  THIS ASSIGNMENT FOR I = 8 IS ONLY VALID FOR THE LISTED Q'S AND NOT
!  FOR OTHER ARRAYS IN THE PROGRAM.
! rate_vec (2026 W3): the one-zone output vector, rows rr_* of
! rotation_scratch_lib -- formerly the 15 separate per-zone arrays
! HR1..HR13, HF1, HF2 written at element IU (the shell number). Callers
! pass the column rate_by_zone(:, k) of their (num_rate_rows, json)
! array. Rows rr_zero9 and rr_zero13 are zero-filled placeholders
! (reactions 9 and 13 are never computed).
subroutine network_rate_vector(log_density, log_temperature, xa, rate_vec)
      use star_info_lib, only: i_h1, i_h2, i_he4, i_he3, i_c12, i_c13, &
           i_n14, i_o16, i_o18
      use phys_const_lib
      implicit none

! xa: the zone's composition in the star%xa layout (i_* slots)
      double precision, intent(in) :: log_density, log_temperature, xa(:)
      double precision, intent(out) :: rate_vec(num_rate_rows)

      double precision :: mass_frac(num_isotopes)
      type(nuclear_rates) :: nr
      integer :: i

! the same screening composition as compute_zone_energy's, deuterium included
! (2026-10; this path used to leave it out, a ~1e-5 effect while
! deuterium survives)
      mass_frac = 0.0d0
      mass_frac(iso_h1) = xa(i_h1)
      mass_frac(iso_h2) = xa(i_h2)
      mass_frac(iso_he3) = xa(i_he3)
      mass_frac(iso_he4) = xa(i_he4)
      mass_frac(iso_c12) = xa(i_c12)
      mass_frac(iso_c13) = xa(i_c13)
      mass_frac(iso_n14) = xa(i_n14)
      mass_frac(iso_o16) = xa(i_o16)
      mass_frac(iso_o18) = xa(i_o18)
      call nuclear_reaction_rates(log_density, log_temperature, mass_frac, nr)
      do i=1,num_reactions
         if(nr%rate(i).le.1.d-5) nr%rate(i) = 0.0d0
      end do
! per 10^9 yr per amu (gyr_amu_per_sec_gram, sidereal year)
      rate_vec(rr_pp)=nr%rate(r_pp)*gyr_amu_per_sec_gram
      rate_vec(rr_he3_he3)=nr%rate(r_he3he3)*gyr_amu_per_sec_gram
      rate_vec(rr_he3_he4)=nr%rate(r_he3he4)*gyr_amu_per_sec_gram
      rate_vec(rr_c12_p)=nr%rate(r_pc12)*gyr_amu_per_sec_gram
      rate_vec(rr_c13_p)=nr%rate(r_pc13)*gyr_amu_per_sec_gram
      rate_vec(rr_n14_p)=nr%rate(r_pn14)*gyr_amu_per_sec_gram
      rate_vec(rr_o16_p)=nr%rate(r_po16)*gyr_amu_per_sec_gram
      rate_vec(rr_c13_alpha)=nr%rate(r_ac13)*gyr_amu_per_sec_gram
      rate_vec(rr_zero9)=nr%rate(r_ao16_unused)*gyr_amu_per_sec_gram
      rate_vec(rr_c12_alpha)=nr%rate(r_ac12)*gyr_amu_per_sec_gram
      rate_vec(rr_n14_alpha)=nr%rate(r_an14)*gyr_amu_per_sec_gram
      rate_vec(rr_triple_alpha)=nr%rate(r_3alpha)*gyr_amu_per_sec_gram
      rate_vec(rr_zero13)=nr%rate(r_c12c12_unused)*gyr_amu_per_sec_gram
      rate_vec(rr_frac_c12_alpha)=nr%frac_n15_pa
      rate_vec(rr_frac_be7_electron)=nr%frac_be7_ecap
end subroutine network_rate_vector



end module net_lib
