!----------------------------------------------------------------------
! rates
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original rates.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model). Numeric literals (including
! ones missing the customary D0/D-3/etc suffix, e.g. many of the
! DATA-statement values below) are copied character-for-character from
! the original -- do not "correct" them, see project notes.
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
! IU IS THE SHELL NUMBER.
subroutine rates(log_density,log_temperature,hydrogen_fraction, &
     helium_fraction,he3_fraction,c12_fraction,c13_fraction,n14_fraction, &
     o16_fraction,o18_fraction,zone_idx,rate_pp,rate_he3_he3,rate_he3_he4, &
     rate_c12_p,rate_c13_p,rate_n14_p,rate_o16_p,rate_c13_alpha,rate_zero9, &
     rate_c12_alpha,rate_n14_alpha,rate_triple_alpha,rate_zero13, &
     frac_c12_alpha,frac_be7_electron)

      use const_lib
      implicit none
      integer, parameter :: json=5000

      double precision, intent(in) :: log_density, log_temperature, &
           hydrogen_fraction, helium_fraction, he3_fraction, c12_fraction, &
           c13_fraction, n14_fraction, o16_fraction, o18_fraction
      integer, intent(in) :: zone_idx
      double precision, intent(out) :: rate_pp(json), rate_he3_he3(json), &
           rate_he3_he4(json), rate_c12_p(json), rate_c13_p(json), &
           rate_n14_p(json), rate_o16_p(json), rate_c13_alpha(json), &
           rate_zero9(json), rate_c12_alpha(json), rate_n14_alpha(json), &
           rate_triple_alpha(json), rate_zero13(json)
      double precision, intent(out) :: frac_c12_alpha(json), &
           frac_be7_electron(json)





      double precision :: mass_frac(13), rate(13), screening_factor(13), &
           charge_product(13), z53(13), z43(13), z23(13), z86(13)
      double precision :: q1(8), q2(8), q3(8), q4(8), q5(8), q6(7), q7(7), &
           q8(7)
      double precision :: atomic_mass(13), atomic_charge(13)
      integer :: num_isotopes, num_reactions

! ***************************
! ANUC ARE ATOMIC MASS UNITS.
! ***************************
! On 10/13/97, JNB converted ANUC(I) from neutral atomic masses to
! bare nuclear masses.
!
! The scale is the mass of C12 divided by 12 or 931.49432 MeV,
! which is 1.6605402 times 10^{-24} gm. The values are obtained
! by dividing the mass excess (expressed in MeV) by 931.49432 MeV and
! adding to this the atomic mass number, A.  The value for Be7, which
! is used implicitly in this subroutine was, 7.016930, until 13/13/97.
! We now use the bare nuclear mass of 7.014735 .
!
      data atomic_mass/1.008665,1.007276,2.013553,3.015501,3.014933,4.001506, &
           11.996709,13.000064,13.999233,15.990526,17.994772,19.986954, &
           23.978458/,atomic_charge/0.,1.,1.,1.,2.,2.,6.,6.,7.,8.,8.,10.,12./, &
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
! ALL OF THE VALUES OF THE Q1, ...., Q5 HAVE BEEN RECALCULATED, USING
!  WHERE NEEDED NUCLEAR CROSS SECTIONS GIVEN IN TABLES 3.2 AND 3.4 OF NEUTRINO
!  ASTROPHYSICS.
! ******************************************************************
! Q6 IS THE COEFFICIENT OF THE TEMPERATURE TERM IN THE DEFINITION OF
!  TAU, EQUATION 3.10 OF NEUTRINO ASTROPHYSICS.
!  TAU = Q6*(T SUB 9 TO THE (-1/3) POWER ).
! ******************************************************************
! SLIGHT CHANGES HAVE BEEN MADE IN THE PREVIOUS VALUES OF Q6 TO MAKE
!  THE DATA MORE ACCURATE.
! NOTE Q6 IS NEGATIVE.
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
!     1 Q2/1.08749,-.273,-.2085,.7630,-0.4044,-1.60,-1.224, -0.2095/,
!     2 Q3/.93833,-.0648,-.0474,.1626,-.08598,-.3064,-.2139,-0.0595/,
!     3 Q4/0.,0.,0.,4.79,7.456,0.0,.69703,.16762/,
!     4 Q5/0.,0.,0.,2.595,4.032,0.0,0.3097,.12114/,
!     5 Q6/-3.3804,-12.2757,-12.826,-13.6899,-13.7173,-15.2281,
!     5    -16.6925/,
!     6 Q7/20.8964,76.6003,67.8036,69.130,70.3809,69.8517,70.8012/,
!     7 Q8/0.,0.,0.,0.0,0.0,0.0,0.0/,
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
      num_reactions/13/
! For different values of SSTANDARD, check the comments in ENGEB
! *********************************************************************
! THE VALUES OF SSTANDARD(I) ARE TO BE CHANGED FROM UNITY IF THE CROSS
!  SECTION FACTORS ARE NOT THE ONES GIVEN IN NEUTRINO ASTROPHYSICS,
!  TABLE 3.2 AND TABLE 3.4 .
! THE CURRENT VALUES ARE CHANGED TO BE THE PREFERRED VALUES LISTED IN
!  THE LAST COLUMN OF TABLE 1 OF BAHCALL AND PINSONNEAULT (1991),
!   6/14/91.
!  THE VALUE OF SSTANDARD(17) FOR HEP IS TAKEN FROM WOLFS ET AL. , PHYS.
!   REV. LETTERS, 63, 2721 (1989).  IT CORRESPONDS TO AN S SUB 0 =
!   15.3E-20 KEV-BARNS, A FACTOR OF 1.913 LARGER THAN FOR THE OLDER
!   MEASUREMENTS USED IN NEUTRINO ASTROPHYSICS.
! *********************************************************************
! ZPRD IS USED IN SCREENING CALCULATIONS. IT IS THE PRODUCT OF THE
!  CHARGES OF THE INTERACTING IONS. ZPRD WAS CHECKED. Z86 IS USED
!  IN CALCULATING INTERMEDIATE SCREENING AND IS DEFINED BY GRABOSKE ET
!  AL, AP. J. 181, PAGE 465 (1973), IN TABLE 4.  Z86 WAS CHECKED AND
!  SOME NUMERICAL VALUES WERE MADE SLIGHTLY MORE ACCURATE.  Z53, Z43,
!  AND Z23 ARE ALSO DEFINED IN TABLE 4 (SEE ABOVE).  SINCE THEY ARE
!  ONLY USED IN STRONG SCREENING, THE VALUES OF Z53, Z43, AND Z23
!  WERE NOT CHECKED.
      double precision :: c21
      data charge_product/1.,4.,4.,6.,6.,7.,8.,12.,16.,12.,14.,12.,36./,z53/ &
         1.175,3.73,3.73,4.804,4.804,5.385,5.941,9.014,11.24,9.014, &
         10.15,9.104,23.28/,z43/0.52,1.31,1.31,1.488,1.488,1.61,1.721, &
         2.577,3.025,2.577,2.81,2.577,5.668/,z23/-0.413,-0.655,-0.655, &
         -0.643,-0.643,-0.659,-0.673,-0.889,-0.946,-0.889,-0.92,-0.889, &
         -1.36/,z86/1.630,5.917,5.917,8.302,8.302,9.520,10.716,16.192, &
         20.978,16.192,18.606,16.192,45.6635/, &
         c21/5.240358E-8/
      save
! DEFINE NEXT THE FRACTIONAL ABUNDANCES BY MASS OF THE IMPORTANT
!  ISOTOPES.
! X, Y, Z, XHE3,..., XBE9 ARE THE MASS FRACTIONS OF THE ISOTOPES.
!  THE ABUNDANCES OF NEUTRONS, H2, H3, NE20,AND MG24, WHICH ARE,
!  RESPECTIVELY, XFRAC(I) FOR I = 1,3,4,12,13, ARE NO LONGER USED.

      integer :: i, nz
      double precision :: mu_ion_inv, mu_e_inv, xtr, zeta0, trm
      double precision :: rho_over_mu_e, log_rho_local, density, t9, t9p13, &
           t9p23, t9m13, t9m23, t9m1, t9m2, t9m12, t9m32
      double precision :: fermi_mom_sq, fermi_energy_over_kt, f_prime_over_f, &
           log_degeneracy
      double precision :: lambda0, lambda0_23, lambda0_86, z_curl, z_bar, &
           z_curl_58, z_bar_28, z_bar_13, lambda0_zcurl
      double precision :: weak_screening_u, intermediate_screening_u, &
           strong_screening_u
      double precision :: r1, r2, a1, a2, a3, a4, a5
      double precision :: be7_electron_rate, be7_proton_rate, be7_temp_factor, &
           be7_q_factor, be7p_charge_product, be7p_z86, be7p_screening_u, &
           be7_mass_factor
      double precision :: be7_electron_frac, c12_alpha_frac, o16_gamma_frac
      double precision :: o16_gamma_rate, c12_alpha_n15p_rate

      mass_frac(1) = 0.0
      mass_frac(2) = hydrogen_fraction
      mass_frac(3) = 0.0
      mass_frac(4) = 0.0
      mass_frac(5) = he3_fraction
      mass_frac(6) = helium_fraction
      mass_frac(7) = c12_fraction
      mass_frac(8) = c13_fraction
      mass_frac(9) = n14_fraction
      mass_frac(10) = o16_fraction
      mass_frac(11) = o18_fraction
      mass_frac(12) = 0.0
      mass_frac(13) = 0.0
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
      mu_ion_inv = 0.
      mu_e_inv = 0.
      xtr = 0.
      zeta0 = 0.
      do 10 i = 1,num_isotopes
         trm = mass_frac(i)/atomic_mass(i)
         mu_ion_inv = mu_ion_inv+trm
         mu_e_inv = mu_e_inv+trm*atomic_charge(i)
         xtr = xtr+trm*atomic_charge(i)**1.58
         zeta0 = zeta0+trm*atomic_charge(i)**2
   10 continue
! DL AND DT ARE THE THE LOG10 OF THE DENSITY AND TEMPERATURE.
!  THE UNIT OF TEMPERATURE IS 10^9 K AND THE UNIT OF DENSITY IS
!  GM PER CM^3 .
! DD = LOG RHO TO THE BASE 10.
! CLN = LN10.  CLN IS CONVERSION BETWEEN LOG10 AND LN.
! CONVERT DENSITY TO UNLOGGED FORM.
! RWE = RHO/(MU SUB E), I. E., THE NUMBER OF ELECTRONS DIVIDED BY
!  AVOGADRO'S NUMBER.
      rho_over_mu_e = ( exp(ln10*log_density) )*mu_e_inv
! THE EXPRESSION FOR RWE WAS INCORRECT IN THE ORIGINAL YALE SUBROUTINE.
!  THE ORIGINAL VERSION HAD ( EXP(CLN*DL) ) DIVIDED BY EMU INSTEAD OF
!  MULTIPLIED BY EMU.  RWE IS USED LATER IN COMPUTING THE SCREENING
!  CORRECTION.
      log_rho_local = log_density
! SET RATES EQUAL TO ZERO FOR THE LOG_10(T) < TCUT(1)
      if(log_temperature.le.tcut(1)) then
         do 20 i = 1,num_reactions
            rate(i) = 0.
   20    continue
         go to 200
      endif
! T9P13 IS THE TEMPERATURE IN UNITS OF 10^9 DEGREES K TO THE PLUS 1/3
!  POWER.  MINUS IS DENOTED BY M.  HERE T9 IS THE TEMPERATURE IN UNITS
!  OF 10^9 K, CONVERTED FROM THE LOG_10 (T) AND RHO IS THE DENSITY IN
!  CGS UNITS.
      density=exp(ln10*log_rho_local)
      t9 = exp(ln10*(log_temperature - 9.0d0))
      t9p13 = t9**cc13
      t9p23 = t9p13**2
      t9m13=1./t9p13
      t9m23=t9m13**2
      t9m1=1./t9
      t9m2=t9m1**2
      t9m12=1./dsqrt(t9)
      t9m32=t9m1*t9m12
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
      fermi_mom_sq=1.017677E-4*rho_over_mu_e**0.6666667
      fermi_energy_over_kt=5.92986*t9m1*(dsqrt(1.+fermi_mom_sq)-1.)
      if(fermi_energy_over_kt.le.1.E-2) then
         f_prime_over_f=1.0
      else
         log_degeneracy=dlog10(fermi_energy_over_kt)
         if(log_degeneracy.ge.1.5) then
            f_prime_over_f=0.0
         else
            f_prime_over_f=0.75793-0.54621*log_degeneracy-0.30964*log_degeneracy**2+0.12535*log_degeneracy**3+ &
            0.1203*log_degeneracy**4-0.012857*log_degeneracy**5-0.014768*log_degeneracy**6
         endif
      endif
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
      lambda0=5.9426E-6*t9m32*dsqrt(density*mu_ion_inv)
      lambda0_23=lambda0**0.666667
      lambda0_86=lambda0**0.86
      z_curl=dsqrt((zeta0+f_prime_over_f*mu_e_inv)/mu_ion_inv)
      z_bar=mu_e_inv/mu_ion_inv
      z_curl_58=z_curl**0.58
      z_bar_28=z_bar**0.28
      z_bar_13=z_bar**cc13
      lambda0_zcurl=lambda0*z_curl
! COMPUTE SCREENING FOR EACH OF THE REACTIONS.
      do 30 i=1,num_reactions
         weak_screening_u=lambda0_zcurl*charge_product(i)
         if(weak_screening_u.le.weak_screening_threshold) then
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
            screening_factor(i)=weak_screening_u
         else
            intermediate_screening_u=0.38*lambda0_86*xtr*z86(i)/(mu_ion_inv*z_curl_58*z_bar_28)
            if(weak_screening_u.le.2.) then
               screening_factor(i)=intermediate_screening_u
            else
               strong_screening_u=0.624*z_bar_13*lambda0_23*(z53(i)+0.316*z_bar_13*z43(i)+0.737* &
               z23(i)/(z_bar*lambda0_23))
               if(strong_screening_u.lt.intermediate_screening_u.or.weak_screening_u.ge.5.) then
                  screening_factor(i)=strong_screening_u
               else
                  screening_factor(i)=intermediate_screening_u
               endif
            endif
         endif
   30 continue
! ****************************************************************
! END OF SCREENING CALCULATION. WEAK AND INTERMEDIATE SCREENING FORMS
!  ARE GIVEN CORRECTLY.  STRONG SCREENING WAS NOT CHECKED BECAUSE IT IS
!  NOT RELEVANT FOR THE SUN.
! ****************************************************************
      nz=1
      if(hydrogen_fraction.eq.0.0) then
         be7_electron_frac=0.
         c12_alpha_frac=0.
         goto 50
      endif
      nz=8
! **************************************************************
!  CALCULATE REACTION RATES FOR THE THREE PRINCIPAL REACTIONS OF
!   THE PP CHAIN: PP, HE3+HE3, HE3 +HE4, AND THE FOUR PROTON
!   BURNING REACTIONS ON C12: C13, N14, AND O16.
! **************************************************************
! R1 IS (T SUB 9)^(-3/2) TIMES (S SUB EFF)/(S SUB 0). THE CORRECT
!  EXPRESSION FOR (S SUB EFF)/(S SUB 0) IS GIVEN IN EQUATION 3.14 OF
!  NEUTRINO ASTROPHYSICS.  THE NUMERICAL FORM THAT IS USED IS EQUATION
!  52 OF FOWLER, CAUHLAN, AND ZIMMERMAN, VOL. 5, 1967.
! RATE(I) IS THE RATE OF THE DIFFERENT REACTIONS PER SECOND PER GRAM,
!  EXCEPT THAT THE MASS FRACTIONS ARE OMITTED.
! THE PREVIOUS YALE VERSION HAD THE RATE FOR THE O16 + P REACTION
!  MULTIPLIED BY T9**(-1/7). THIS FACTOR IS INCORRECT AND HAS BEEN
!  REMOVED; IT APPEARED BEFORE AS AN IF STATEMENT REFERRING ONLY TO
!  RATE(7).
      do 40 i=1,7
!         R1=T9M23+Q1(I)*T9M13+Q2(I)+Q3(I)*T9P13+Q4(I)*T9P23+Q5(I)*T9
! MHP 8/14 RATES CORRECTED TO PERMIT USER MODIFICATION OF REACTION
! RATE DERIVATIVES
         r1=t9m23+q1(i)*t9m13+qs0e_scale(i)*(q2(i)+q3(i)*t9p13)+ &
            qqs0ee_scale(i)*(q4(i)*t9p23+q5(i)*t9)
         rate(i)=density*r1*exp(q6(i)*t9m13+q7(i)+(q8(i)*t9)**2+screening_factor(i))
         rate(i) = rate(i)*cross_section_scale(i)
         if(rate(i).lt.1.E-30) rate(i)=0.0d0
   40 continue
! ***************************************************************
! END OF CALCULATION OF REACTION RATES FOR FIRST 7 REACTIONS.
! ***************************************************************
! ***********************************************
! BE7 BURNING
! ***********************************************
! CALCULATE THE BURNING OF BE7 BY PROTONS (WHICH PRODUCES THE MOST
!  EXPERIMENTALLY ACCESSIBLE SOLAR NEUTRINOS) AND THE BURNING OF BE7
!  BY ELECTRON CAPTURE (THE DOMINANT PROCESS). THE PREVIOUSLY USED
!  EXPRESSION CONTAINED A TERM THAT WAS PHYSICALLY INCORRECT
!  AND ALSO USED CROSS SECTION FACTORS THAT WERE OUT OF DATE.
! THE ELECTRON CAPTURE RATE IN SEC^(-1) IS GIVEN BY EQUATION (3.18) OF
!  NEUTRINO ASTROPHYSICS.  CAN OMIT THE FACTOR OF RHO (IN CGS UNITS)
!  WHICH ALSO APPEARS IN THE BE7PROTON RATE. THE BE7 PROTON CAPTURE
!  RATE IS GIVEN BY TABLE 3.2 AND EQUATION (3.12).  NOTE THAT 1 OVER MU
!  SUB E IS EQUAL TO EMU.
! HERE WE USE THE NOTATION BE7 + PROTON = BE7PROTON AND
!  BE7 + E = BE7ELECTRON.
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
!
         be7_electron_rate = (1.752E-10)*t9m12*(1.0 + 0.004*(1000.*t9-16.))
         be7_electron_rate = be7_electron_rate*mu_e_inv*cross_section_scale(15)
         be7_temp_factor = (-10.2625*t9m13)
         be7_proton_rate = (3.128813E+5)*hydrogen_fraction*cross_section_scale(16)*exp(be7_temp_factor)

! INCLUDE FOR BE7PROTON THE T9M23 FACTOR AND ALL CORRECTIONS PROPORTIONAL TO
!  Q1,...,Q5 FROM EQUATION 3.14 OF NEUTRINO ASTROPHYSICS. THESE
!  CORRECTIONS ARE DEFINED EARLIER IN THIS SUBROUTINE.
!         QRBE7 = T9M23 + Q1(8)*T9M13 + Q2(8)+ Q3(8)*T9P13
!     $          + Q4(8)*T9P23 + Q5(8)*T9
! MHP 9/14 ADDED THE ABILITY TO ALTER DERIVATIVES INDEPENDENTLY
         be7_q_factor = t9m23 + q1(8)*t9m13 + qs0e_scale(8)*(q2(8)+ q3(8)*t9p13) &
                + qqs0ee_scale(8)*(q4(8)*t9p23 + q5(8)*t9)
         be7_proton_rate = be7_proton_rate*be7_q_factor
! CALCULATE THE SCREENING CORRECTION FOR BE7 + P REACTION.  USE WEAK AND
!  INTERMEDIATE SCREENING FORMULAE.
         be7p_charge_product = 4.0
         be7p_z86 = 5.7790
         weak_screening_u = lambda0_zcurl*be7p_charge_product
         if(weak_screening_u.le.weak_screening_threshold) then
            be7p_screening_u = weak_screening_u
         else
            intermediate_screening_u = 0.38*lambda0_86*xtr*be7p_z86/(mu_ion_inv*z_curl_58*z_bar_28)
            be7p_screening_u = intermediate_screening_u
          endif
         be7_proton_rate = be7_proton_rate*exp(be7p_screening_u)
! END OF CALCULATION OF SCREENING CORRECTION FOR BE7 + P REACTION.
! MULTIPLY RATES BY FACTOR OF RHO/[(ATOMIC MASS UNIT)*A(BE7)] TO GET
!  IN UNITS OF GM^{-1}.  CALL FACTOR CAMUBE7. Corrected for Be7 bare
! mass number. Previously used 8.582295E+22 in the expression below.
! These corrections were made on 10/14/97.
!
         be7_mass_factor = density*8.584981E+22
         be7_proton_rate = be7_mass_factor*be7_proton_rate
         be7_electron_rate = be7_mass_factor*be7_electron_rate
! END OF MULTIPLICATION INSERTED NOVEMBER 6, 1990.
         be7_electron_frac = be7_electron_rate/(be7_electron_rate + be7_proton_rate)
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
      o16_gamma_rate = t9m23 + 0.0273016*t9m13 + 0.14374 + 0.027490*t9p13 &
                + 6.14685*t9p23 + 2.98940*t9
! MULTIPLY BY THE VALUE OF S0 IN KEV-B.
!      O16GAMMA = O16GAMMA*64.
! MHP 9/14 ADDED THE OPTION TO MODIFY THE RELATIVE CROSS SECTIONS
! FOR N15+P -> C12+ALPHA AND O16+GAMMA
      o16_gamma_rate = o16_gamma_rate*64*o16_gamma_scale
!
      c12_alpha_n15p_rate = t9m23 + 0.0273016*t9m13 + 2.01186 + 0.384763*t9p13 &
                + 17.0579*t9p23 + 8.29580*t9
!      C12ALPHA = C12ALPHA*67500.
      c12_alpha_n15p_rate = c12_alpha_n15p_rate*67500*c12_alpha_scale
      c12_alpha_frac = c12_alpha_n15p_rate/(c12_alpha_n15p_rate + o16_gamma_rate)
      o16_gamma_frac = 1.0d0 - c12_alpha_frac
! END OF NEW ROUTINE FOR THE BRANCHING OF N15 + P .
   50 do 60 i=nz,num_reactions
         rate(i)=0.
   60 continue
!***MHP 3/91 ALPHA CAPTURE REACTIONS UPDATED TO CAUGHLAN AND FOWLER(1988)
!   RATES.  THE RATES ARE EXPRESSED IN THE SAME TERMS USED BY CZ, WITH
!   THE CONVERSION FACTOR IN THE FRONT OBTAINED FROM VANDENBERG'S
!   NOTES ON THE REACTION RATES.
!  RATE(8)  HE4+C13
!  RATE(10) HE4+C12=>O16
!  RATE(11) HE4+N14=>O18
!  RATE(12) TRIPLE ALPHA
      if(log_temperature.lt.tcut(4)) go to 100
! C13(ALPHA,N) O16
      r1=t9m23+0.0129d0*t9m13+2.04d0+0.184d0*t9p13
      a1 = 6.77d15*exp(-32.329d0*t9m13-(t9/1.284d0)**2)
      a2 = 3.82d5*exp(-9.373*t9m1)
      a3 = 1.41d6*exp(-11.873*t9m1)
      a4 = 2.00d9*exp(-20.409*t9m1)
      a5 = 2.92d9*exp(-29.283*t9m1)
      rate(8) = 1.157126d22*density*exp(screening_factor(8))*(a1*r1+t9m32* &
           (a2+a3+a4+a5))
! C12(ALPHA,GAMMA)O16
      r1 = 1.0d0/(1.0d0+0.0489d0*t9m23)
      r2 = 1.0d0/(1.0d0+0.2654d0*t9m23)
      a1 = t9m2*exp(-32.120*t9m13)
      a2 = 1.04d8*r1**2*exp(-(t9/3.496)**2)
      a3 = 1.76d8*r2**2
      a4 = 1.25d3*t9m32*exp(-27.499*t9m1)
      a5 = 1.43d-2*t9**5*exp(-15.541*t9m1)
      rate(10) = 1.25388d22*density*exp(screening_factor(10))*(a1*(a2+a3)+a4+a5)
! N14(ALPHA,GAMMA)F18 + F18=>O18+EPLUS+NU
      r1 = t9m23+0.012d0*t9m13+1.45d0+0.177d0*t9p13+1.97d0*t9p23 &
           +0.406d0*t9
      a1 = 7.78d9*exp(-36.031d0*t9m13-(t9/0.881d0)**2)
      a2 = t9m32*2.36d-10*exp(-2.798d0*t9m1)
      a3 = t9m32*2.03d0*exp(-5.054d0*t9m1)
      a4 = t9m23*1.15d4*exp(-12.310*t9m1)
      rate(11)= 1.07452d22*density*exp(screening_factor(11))*(a1*r1+a2+a3+a4)
! TRIPLE ALPHA
      rate(12) = 1.565315d21*density**2*t9m1*t9m2*2.79E-8* &
                 exp(-4.4027*t9m1+screening_factor(12))
  100 continue
      rate(9) = 0.0d0
      rate(13) = 0.0d0
! END OF XEROING OUT OF REACTIONS 9 AND 13.
      do 130 i=1,num_reactions
         if(rate(i).le.1.E-5) rate(i) = 0.0
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
  200 rate_pp(zone_idx)=rate(1)*c21
      rate_he3_he3(zone_idx)=rate(2)*c21
      rate_he3_he4(zone_idx)=rate(3)*c21
      rate_c12_p(zone_idx)=rate(4)*c21
      rate_c13_p(zone_idx)=rate(5)*c21
      rate_n14_p(zone_idx)=rate(6)*c21
      rate_o16_p(zone_idx)=rate(7)*c21
      rate_c13_alpha(zone_idx)=rate(8)*c21
      rate_zero9(zone_idx)=rate(9)*c21
      rate_c12_alpha(zone_idx)=rate(10)*c21
      rate_n14_alpha(zone_idx)=rate(11)*c21
      rate_triple_alpha(zone_idx)=rate(12)*c21
      rate_zero13(zone_idx)=rate(13)*c21
      frac_c12_alpha(zone_idx)=c12_alpha_frac
      frac_be7_electron(zone_idx)=be7_electron_frac
! ****************************************
! END OF COMPUTATION OF HRK(IU).
! ****************************************
      return
end subroutine rates
