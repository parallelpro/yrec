!----------------------------------------------------------------------
! nuclear_energy_lib
!----------------------------------------------------------------------
! The nuclear energy generation of one zone (2026-10: the physics of
! the former engeb, moved here and made free of model-state side
! effects; core/burn_step's compute_zone_energy is now the star-layer
! adapter that stores the results).
! nuclear_energy returns, in a nuclear_energy_result: eps and its
! absolute derivatives d eps/dln rho, d eps/dln T for the Henyey
! Jacobian; the chain shares (pp, Be7 electron/proton branches, CNO,
! triple alpha); the alpha-capture energy, thermal neutrino losses and
! He3 energy rates; and, when ctrl%calc_neutrinos, the solar neutrino
! fluxes and the Be7 abundance. Everything it computes comes back in
! the record. It reads only configuration from star (ctrl, the
! cross-section scales).
module nuclear_energy_lib
      use net_lib
      implicit none
      private

      type, public :: nuclear_energy_result
! .false. at or below nuclear_logT_cutoffs(1): nothing else is set
         logical :: burning = .false.
! erg/g/s; derivatives absolute (d eps/dln rho, d eps/dln T)
         double precision :: total_energy_gen_rate = 0.0d0
         double precision :: deps_dlnrho = 0.0d0, deps_dlnt = 0.0d0
         double precision :: pp_chain_energy_gen = 0.0d0
         double precision :: he3he4_be7_electron_energy_gen = 0.0d0
         double precision :: he3he4_be7_proton_energy_gen = 0.0d0
         double precision :: cno_cycle_energy_gen = 0.0d0
         double precision :: triple_alpha_energy_gen = 0.0d0
         double precision :: eps_alpha_capture = 0.0d0
         double precision :: eps_neutrino_loss = 0.0d0
         double precision :: eps_he3_he3 = 0.0d0, eps_he3_burning = 0.0d0
! set only when ctrl%calc_neutrinos (fluxes_set) and X > 1e-6 (be7_set)
         logical :: fluxes_set = .false., be7_set = .false.
         double precision :: neutrino_flux(10) = 0.0d0
         double precision :: be7_mass_fraction = 0.0d0
      end type nuclear_energy_result

      public :: nuclear_energy

contains

!----------------------------------------------------------------------
! nuclear_energy -- history and sources (the header of engeb.f, kept
! with the physics it describes)
!----------------------------------------------------------------------
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
!  GRAMS.  THE FLUXES ARE STORED IN star%neutrino_flux.  I ALSO CALCULATE
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
! Legacy names (engeb.f) of the result components and arguments:
!   EPP1 -> pp_chain_energy_gen        EPP2 -> he3he4_be7_electron_energy_gen
!   EPP3 -> he3he4_be7_proton_energy_gen   ECN -> cno_cycle_energy_gen
!   E3AL -> triple_alpha_energy_gen    PEP -> deps_dlnrho
!   PET -> deps_dlnt                   SUM1 -> total_energy_gen_rate
!   DL -> log_density                  TL -> log_temperature
!   X -> hydrogen_fraction             Y -> helium_fraction
!   XHE3 -> he3_fraction               XC12/XC13 -> c12_fraction/c13_fraction
!   XN14 -> n14_fraction               XO16/XO18 -> o16_fraction/o18_fraction
!   XH2 -> deuterium_fraction
subroutine nuclear_energy(log_density, &
     log_temperature, hydrogen_fraction, helium_fraction, he3_fraction, &
     c12_fraction, c13_fraction, n14_fraction, o16_fraction, &
     o18_fraction, deuterium_fraction, in_surface_cz, turnover_timescale, e)

      use star_info_lib, only: star, i_nu_b8, i_nu_be7, i_nu_f17, i_nu_hep, i_nu_n13, i_nu_o15, i_nu_pep, i_nu_pp, &
           i_nu_he3he3, i_nu_he3he4
      use phys_const_lib
      use math_lib
      implicit none

      double precision, intent(in) :: log_density, log_temperature, &
           hydrogen_fraction, helium_fraction, he3_fraction, c12_fraction, &
           c13_fraction, n14_fraction, o16_fraction, o18_fraction, &
           deuterium_fraction
! the zone lies in the surface convection zone (deuterium burning is
! capped at the overturn rate there when turnover_timescale > 1 s)
      logical, intent(in) :: in_surface_cz
      double precision, intent(in) :: turnover_timescale
      type(nuclear_energy_result), intent(out) :: e

! neutrino_dsnu_dt/neutrino_dsnu_drho are the absolute derivatives
! d snu/dT, d snu/d rho returned by thermal_neutrino_losses (the Itoh
! 1996 path, selected by star%ctrl%use_itoh_neutrino_loss); they are set
! and consumed only inside thermal_neutrino_loss below.
      double precision :: neutrino_dsnu_dt, neutrino_dsnu_drho

! nr: this zone's reaction rates, derivatives and branchings, from
! net_lib's nuclear_reaction_rates (2026-10: the one copy of the rate
! block that compute_zone_energy and rates each used to carry)
      type(nuclear_rates) :: nr
      double precision :: mass_fraction(num_isotopes), &
           reaction_rate(num_reactions), dlnrate_dlnrho(num_reactions), &
           dlnrate_dlnt(num_reactions), reaction_energy_gen(num_reactions), &
           eg(num_reactions)
      double precision :: v1(7), v2(7), v3(7)
      data v1,v2,v3/ &
           6.002d19,2.084d20,1.872d21,9.383d-1,-4.141d-1,5.829d-2,5.5924d0, &
           4.886d10,7.580d10,6.023d10,6.290d-3,7.483d-3,3.061d-4,1.5654d0, &
           2.320d-7,8.449d-8,1.787d-8,2.581d-2,1.734d-2,6.990d-4,0.56457d0/
      double precision :: dd, density, t9, t9_p13, t9_p23, t9_m13, t9_m23, t9_m12
      double precision :: dgdeut, qrtdeut, rdeut, qdeut, zz, &
           tfacdeut, tfacdeut2, rdeutmax, rdeut2
      double precision :: frac_be7_ecap, frac_be7_pcap, frac_n15_pa, frac_n15_pg
      double precision :: convert, egdeut
      double precision :: sum2, sum3
      double precision :: fourpiau2
      double precision :: flux_value
      double precision :: carbon_fraction_total, oxygen_fraction_total, &
           neutrino_temp, neutrino_density, neutrino_loss_snu
      double precision :: el, eli, ez, ez3, emue, ex1, ex2, ex3, el2
      double precision :: polx10, polx11, polx12, polx21, polx22, polx31, &
           polx32, qedn, qetn, qetnx, qednx
      integer :: i, k

      call fill_mass_fractions
! SET RATES EQUAL TO ZERO FOR THE LOG_10(T) < TCUT(1)
      if (log_temperature.le.star%ctrl%nuclear_logT_cutoffs(1)) then
! MHP 5/02 DEUTERIUM BURNING
         dgdeut = 0.0d0
         qrtdeut = 0.0d0
         e%deps_dlnrho = 0.d0
         e%deps_dlnt = 0.d0
         do i = 1,num_reactions
            eg(i) = 0.d0
            reaction_rate(i) = 0.d0
            reaction_energy_gen(i) = 0.d0
         end do
         return
      end if
      e%burning = .true.
! the rates, their logarithmic derivatives and the branchings
      call nuclear_reaction_rates(log_density, log_temperature, mass_fraction, nr)
      density = nr%density
      t9 = nr%t9
      t9_p13 = nr%t9_p13
      t9_p23 = nr%t9_p23
      t9_m13 = nr%t9_m13
      t9_m23 = nr%t9_m23
      t9_m12 = nr%t9_m12
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
! 2026 (bugsweep Batch 2): T9**(-2/3), as in CF88's d(p,g)3He fit and
! in qrtdeut below (which already differentiates the -2/3 form); was
! T9**(+2/3) (engeb.f). 2026-10: the rate itself is deuterium_proton_capture_rate, shared
! with deutrate; zz and tfacdeut stay for qrtdeut.
        rdeut = deuterium_proton_capture_rate(density, t9, t9_p13, t9_p23, t9_m13, t9_m23)
        tfacdeut2 = 0.112d0*t9_p13+6.76d0*t9_p23+7.95d0*t9
        qrtdeut = cc13*((tfacdeut2/tfacdeut) -2.0d0 - zz)
! NOW LIMIT DEUTERIUM BURNING IN A SURFACE CZ TO BE ON A TIME SCALE
! NO SHORTER THAN THE CONVECTIVE OVERTURN TIMESCALE.
        if (in_surface_cz .and. turnover_timescale.gt.1.0d0) then
           rdeutmax = avogadro_dpg/deuteron_mass_amu/turnover_timescale
           rdeut2 = rdeut*hydrogen_fraction
           if (rdeut2.gt.rdeutmax) then
              if (hydrogen_fraction.gt.1.0d-6) then
                 rdeut = rdeutmax/hydrogen_fraction
              end if
           end if
         end if
      end if
! the energy path drops reactions 1-7 below 1e-30 (with their
! derivatives) before summing eps; rates floors at 1e-5 instead
      reaction_rate = nr%rate
      dlnrate_dlnrho = nr%dlnrate_dlnrho
      dlnrate_dlnt = nr%dlnrate_dlnt
      do i=r_pp,r_po16
         if (reaction_rate(i).lt.1.d-30) then
            reaction_rate(i)=0.d0
            dlnrate_dlnrho(i)=0.d0
            dlnrate_dlnt(i)=0.d0
         end if
      end do
      frac_be7_ecap = nr%frac_be7_ecap
      frac_be7_pcap = nr%frac_be7_pcap
      frac_n15_pa = nr%frac_n15_pa
      frac_n15_pg = nr%frac_n15_pg
! *****************************************************************
! EG(I) = rate times the reactant mass fractions; then energies, the
! derivatives for the Henyey Jacobian, and the neutrino fluxes.
! *****************************************************************
      eg(r_pp)=reaction_rate(r_pp)*hydrogen_fraction*hydrogen_fraction
! MHP 5/02 ADD DEUTERIUM BURNING IF RELEVANT
      if (deuterium_fraction.gt.1.0d-11) then
         egdeut = rdeut*hydrogen_fraction*deuterium_fraction
      else
         egdeut = 0.0d0
      end if
      eg(r_he3he3)=reaction_rate(r_he3he3)*he3_fraction*he3_fraction
      eg(r_he3he4)=reaction_rate(r_he3he4)*he3_fraction*helium_fraction
      eg(r_pc12)=reaction_rate(r_pc12)*hydrogen_fraction*c12_fraction
      eg(r_pc13)=reaction_rate(r_pc13)*hydrogen_fraction*c13_fraction
      eg(r_pn14)=reaction_rate(r_pn14)*hydrogen_fraction*n14_fraction
      eg(r_po16)=reaction_rate(r_po16)*hydrogen_fraction*o16_fraction
      eg(r_ac13)=reaction_rate(r_ac13)*helium_fraction*c13_fraction
!     EG(9)=RATE(9)*Y*X016
      eg(r_ac12)=reaction_rate(r_ac12)*helium_fraction*c12_fraction
      eg(r_an14)=reaction_rate(r_an14)*helium_fraction*n14_fraction
      eg(r_3alpha)=reaction_rate(r_3alpha)*helium_fraction**3
!     EG(13)=RATE(13)*XC12*XC12
! ******************************************************************
! ****************************************
      call sum_energy_generation
      call chain_energy_shares
      call solar_neutrino_fluxes
      call thermal_neutrino_loss
      return

contains


! ---------------------------------------------------------------
! Zero the neutrino/alpha-capture yields, form the fractional
! abundances of the burning species, the ion and electron mean
! molecular weights, and the screening precursors (xtr, zet).
subroutine fill_mass_fractions
      use math_lib
! ZERO OUT THE ENERGY YIELDS FROM NEUTRINOS(ENU) AND ALPHA CAPTURE
! REACTIONS (EALPCA).
! DEFINE NEXT THE FRACTIONAL ABUNDANCES BY MASS OF THE IMPORTANT
!  ISOTOPES.
! X, Y, Z, XHE3,..., XBE9 ARE THE MASS FRACTIONS OF THE ISOTOPES.
!  THE ABUNDANCES OF NEUTRONS, H2, H3, NE20,AND MG24, WHICH ARE,
!  RESPECTIVELY, XFRAC(I) FOR I = 1,3,4,12,13, ARE NO LONGER USED.
      mass_fraction(iso_n) = 0.0d0
      mass_fraction(iso_h1) = hydrogen_fraction
! MHP 5/02 ADDED DEUTERIUM
!      XFRAC(3) = 0.0
      mass_fraction(iso_h2) = deuterium_fraction
      mass_fraction(iso_h3) = 0.0d0
      mass_fraction(iso_he3) = he3_fraction
      mass_fraction(iso_he4) = helium_fraction
      mass_fraction(iso_c12) = c12_fraction
      mass_fraction(iso_c13) = c13_fraction
      mass_fraction(iso_n14) = n14_fraction
      mass_fraction(iso_o16) = o16_fraction
      mass_fraction(iso_o18) = o18_fraction
      mass_fraction(iso_ne20) = 0.0d0
      mass_fraction(iso_mg24) = 0.0d0
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
      dd = log_density
end subroutine fill_mass_fractions

! ---------------------------------------------------------------
! Energy release per reaction (MeV -> erg), reactions 9 and 13
! zeroed as in the original Yale code, the total energy generation
! with d ln eps / d ln rho and d ln eps / d ln T, and the global
! outputs returned to the caller.
subroutine sum_energy_generation
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
      convert = 1.602177d-6
! THE MULTIPLYING CONSTANTS BELOW ARE IN MEV.
! JNB changed the pp energy release by 0.002 MeV because of a better
! estimate of the neutrino energy loss on 9/25/97. See pg. 139 of
! Vol. 19 of my notes.  On pgs. 139-141, I document other small changes
! to this energy generation. No large changes; all of order keV changes
! except for the rare 8B reaction.  9/28/97.
!
      reaction_energy_gen(r_pp)=eg(r_pp)*6.664d0*convert
! MHP 5/02 ADD DEUTERIUM BURNING
      if (deuterium_fraction.gt.1.0d-11) then
         dgdeut = egdeut*qdeut*convert
      else
         dgdeut = 0.0d0
      end if
      reaction_energy_gen(r_he3he3)=eg(r_he3he3)*12.860d0*convert
      reaction_energy_gen(r_he3he4)=eg(r_he3he4)*(1.586d0+frac_be7_ecap*17.394d0+frac_be7_pcap*11.499d0)*convert
      reaction_energy_gen(r_pc12)=eg(r_pc12)*3.457372d0*convert
      reaction_energy_gen(r_pc13)=eg(r_pc13)*7.550628d0*convert
      reaction_energy_gen(r_pn14)=eg(r_pn14)*(9.054d0+frac_n15_pa*4.966d0+frac_n15_pg*12.128d0)*convert
      reaction_energy_gen(r_po16)=eg(r_po16)*3.553d0*convert
      reaction_energy_gen(r_ac13)=eg(r_ac13)*2.216d0*convert
      reaction_energy_gen(r_ac12)=eg(r_ac12)*7.162d0*convert
      reaction_energy_gen(r_an14)=eg(r_an14)*5.815d0*convert
      reaction_energy_gen(r_3alpha)=eg(r_3alpha)*7.275d0*convert
! JVS 10/11 Need to grab He3 energy generation
      e%eps_he3_he3 = reaction_energy_gen(r_he3he3)
      e%eps_he3_burning = reaction_energy_gen(r_he3he3)+reaction_energy_gen(r_he3he4)
! JVS end

! *******************************************************************
! END OF CALCULATION OF ENERGY RELEASE.
! *******************************************************************
! SET TO ZERO O16+ALPHA AND C12+C12 RATES.
      dlnrate_dlnrho(r_ao16_unused) = 0.0d0
      dlnrate_dlnt(r_ao16_unused) = 0.0d0
      reaction_energy_gen(r_ao16_unused) = 0.0d0
      dlnrate_dlnrho(r_c12c12_unused) = 0.0d0
      dlnrate_dlnt(r_c12c12_unused) = 0.0d0
      reaction_energy_gen(r_c12c12_unused) = 0.0d0
! END OF XEROING OUT OF REACTIONS 9 AND 13.
      e%total_energy_gen_rate=0.0d0
      sum2=0.0d0
      sum3=0.0d0
      do i=1,num_reactions
! *******************************************************************
! SUM OF THE TOTAL ENERGY GENERATION IN ERGS PER GRM PER SECOND WITH
! DERIVATIVES WITH RESPECT TO DENSITY AND TO TEMPERATURE.
! *******************************************************************
! SUM1 = SUM OF ALL ENERGY GENERATION RATES. NOTE THAT THE BURNING OF
!  BE7 IS INCLUDED IN DG(3) ABOVE.
! SUM2 = SUM OVER I OF DG(I)* [D LOG RATE(I) / D LOG RHO ].
! SUM3 = SUM OVER I OF DG(I)* [D LOG RATE(I) / D LOG T ].
         e%total_energy_gen_rate=e%total_energy_gen_rate+reaction_energy_gen(i)
         sum2=sum2+reaction_energy_gen(i)*dlnrate_dlnrho(i)
         sum3=sum3+reaction_energy_gen(i)*dlnrate_dlnt(i)
      end do
! MHP 5/02 ADD DEUTERIUM BURNING
      e%total_energy_gen_rate = e%total_energy_gen_rate + dgdeut
      sum2 = sum2 + dgdeut
      sum3 = sum3 + dgdeut*qrtdeut
      if (e%total_energy_gen_rate.le.1.d-12) then
         e%deps_dlnrho=0.d0
         e%deps_dlnt=0.d0
         do i=1,num_reactions
            eg(i)=0.d0
         end do
      else
! ******************************************************
! GLOBAL QUANTITIES THAT ARE RETURNED BY THE SUBROUTINE.
! ******************************************************
! PEP AND PET ARE THE DERIVATIVES OF THE TOTAL ENERGY GENERATION RATE
!  WITH RESPECT TO DENSITY AND TEMPERATURE: THE ABSOLUTE
!  D EPS / D LN RHO AND D EPS / D LN T (ERG/G/S, NATURAL LOG), I.E.
!  SUM OVER I OF EPS(I) * D LN RATE(I) / D LN (RHO OR T) -- NOT THE
!  LOGARITHMIC D LN EPS. NO CALLER DIVIDES BY EPS: THE HENYEY JACOBIAN
!  AND THE PULSATION OUTPUT (GYRE'S EPS*EPS_T) BOTH WANT THIS FORM.
         e%deps_dlnrho = sum2
         e%deps_dlnt = sum3
      end if
! *****************************************************
! END OF COMPUTATION OF THE GLOBAL QUANTITIES.
! *****************************************************
end subroutine sum_energy_generation

! ---------------------------------------------------------------
! The energy of each chain (pp, the two Be7 branches, CNO, triple
! alpha) and of the alpha captures, from the per-reaction energies
! summed in sum_energy_generation. (2026-10: the per-Gyr-per-amu rate
! copies and branchings this routine also produced, once relayed to
! kemcom and later stored per shell, had no reader and are gone.)
subroutine chain_energy_shares
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
      e%pp_chain_energy_gen = reaction_energy_gen(r_pp)+reaction_energy_gen(r_he3he3)+dgdeut
! EPP3 INCLUDES THE ENERGY GENERATED BY THE HE3 + HE4 REACTION AND BY
!  THE BURNING OF BE7 THROUGH PROTON CAPTURE (THE frac_be7_pcap BRANCH
!  OF EG(3), INCLUDING ITS SHARE OF THE 1.586 MEV HE3 + HE4 RELEASE).
      e%he3he4_be7_proton_energy_gen = eg(r_he3he4)*frac_be7_pcap*(1.586d0 + 11.499d0)*convert
! EPP2 INCLUDES THE ENERGY GENERATED BY THE HE3 + HE4 REACTION AND BY
!  THE BURNING OF BE7 THROUGH ELECTRON CAPTURE: THE REST OF DG(3).
      e%he3he4_be7_electron_energy_gen = reaction_energy_gen(r_he3he4) - &
           e%he3he4_be7_proton_energy_gen
! ECN IS THE ENERGY GENERATED THROUGH THE CNO CYCLE.
      e%cno_cycle_energy_gen=reaction_energy_gen(r_pc12)+reaction_energy_gen(r_pc13)+ &
           reaction_energy_gen(r_pn14)+reaction_energy_gen(r_po16)
! E3AL IS THE ENERGY GENERATED THROUGH THE TRIPLE-ALPHA REACTION AND
!  IS NEGLIGIBLE FOR THE SUN.
      e%triple_alpha_energy_gen = reaction_energy_gen(r_3alpha)


! ENERGY FROM ALPHA CAPTURE REACTIONS.
      e%eps_alpha_capture=reaction_energy_gen(r_ac13)+reaction_energy_gen(r_ac12)+ &
           reaction_energy_gen(r_an14)
end subroutine chain_energy_shares

subroutine solar_neutrino_fluxes
      use math_lib
! the solar neutrino fluxes and the Be7 abundance (ctrl%calc_neutrinos)
      if (star%ctrl%calc_neutrinos) then
         e%fluxes_set = .true.
! MHP 9/91 CHANGE TO TURN OFF NEUTRINO CALC FOR HYDROGEN-EXHAUSTED CORE.
         if (hydrogen_fraction.le.1.0d-6) then
            do i=1,10
               e%neutrino_flux(i)=0.0d0
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
         fourpiau2 = 2.812295d+27
! FLUX OF PP NEUTRINOS.
         e%neutrino_flux(i_nu_pp) = eg(r_pp)/fourpiau2
! FLUX OF PEP NEUTRINOS. USE EQUATION 3.17 OF NEUTRINO ASTROPHYSICS.
! Note that should not change SStandard(14) unless the ratio of pep to pp
!  is changed.  Pep rate is explicitly scaled here with respect to the pp
!  rate.
         e%neutrino_flux(i_nu_pep) = (3.4848d-6)*nr%electron_density_na*t9_m12* &
              (1.0d0 + 20.d0*t9)*eg(r_pp)
         e%neutrino_flux(i_nu_pep) = e%neutrino_flux(i_nu_pep)*star%cross_section_scale(s_pep)/fourpiau2
! FLUX OF HEP NEUTRINOS.  USE EQUATION 3.12 DIRECTLY.
! hep: the screened rate factor from nuclear_reaction_rates
         e%neutrino_flux(i_nu_hep) = nr%hep
         e%neutrino_flux(i_nu_hep) = e%neutrino_flux(i_nu_hep)*hydrogen_fraction*he3_fraction/ &
              fourpiau2
! COMPUTE BE7MASSFRACTION. THIS IS NOT REQUIRED FOR THE NEUTRINO
!  FLUXES SINCE BE7 IS ALWAYS IN EQUILIBRIUM WITH THE SLOWER PRODUCTION
!  RATE OF HE3 + HE4.  HOWEVER, IT IS OF INTEREST IN SOME APPLICATIONS
!  TO KNOW THE BE7 MASS FRACTION, SO I COMPUTE IT HERE AND STORE IT
!  IN e%be7_mass_fraction.
         e%be7_set = .true.
         e%be7_mass_fraction = eg(r_he3he4)/(nr%be7_proton + nr%be7_electron)
! END OF NOVEMBER 6, 1990  ADDITION.
! FLUX OF BE7 NEUTRINOS.
         e%neutrino_flux(i_nu_be7) = eg(r_he3he4)*frac_be7_ecap/fourpiau2
! FLUX OF B8 NEUTRINOS.
         e%neutrino_flux(i_nu_b8) = eg(r_he3he4)*frac_be7_pcap/fourpiau2
! FLUX OF N13 NEUTRINOS.
         e%neutrino_flux(i_nu_n13) = eg(r_pc12)/fourpiau2
! FLUX OF O15 NEUTRINOS.
         e%neutrino_flux(i_nu_o15) = eg(r_pn14)/fourpiau2
! FLUX OF F17 NEUTRINOS.
         e%neutrino_flux(i_nu_f17) = eg(r_po16)/fourpiau2
! FLUX OF FICTIONAL HE3 + HE3 NEUTRINOS.
         e%neutrino_flux(i_nu_he3he3) = eg(r_he3he3)/fourpiau2
! FLUX OF FICTIONAL HE3 + HE4 NEUTRINOS.
         e%neutrino_flux(i_nu_he3he4) = eg(r_he3he4)/fourpiau2
! SET UNITS OF NEUTRINO FLUXES TO BE 10**10 PER CM^2 PER SEC PER GM AT THE
!  EARTH. MULTIPLY BY 10**-10.
!  IF THE VALUE FOR THIS SHELL IS NEGLIGIBLY SMALL, SET EQUAL TO ZERO.
         do k = 1,10
            e%neutrino_flux(k) = (1.0d-10)*e%neutrino_flux(k)
            flux_value = e%neutrino_flux(k)
! KC 2025-05-30 CHANGED 1.E-50 TO 0.0 TO AVOID UNDERFLOW
            if (flux_value.le.0.0d0) then
              e%neutrino_flux(k) = 0.0d0
            end if
         end do
! MHP 9/91 ENDIF INSERTED HERE.
         end if
      end if
end subroutine solar_neutrino_fluxes

subroutine thermal_neutrino_loss
      use math_lib
! the thermal neutrino losses above nuclear_logT_cutoffs(5), with their
! absolute derivatives, added to eps
! END OF NEUTRINO FLUX ROUTINE.
!
! ***MHP 5/91
!C CALCULATE NEUTRINO LOSSES FOR NEUTRINO-COOLED CORES OF EVOLVED STARS.
! 3/92 DBG Added option to use new (more sophisticated) neutrino loss
! routines.  See subroutine NEUTR for complete description.


      if (log_temperature.le.star%ctrl%nuclear_logT_cutoffs(5)) return


          carbon_fraction_total = c12_fraction+c13_fraction
          oxygen_fraction_total = o16_fraction+o18_fraction
          neutrino_temp=exp10(log_temperature)
          neutrino_density=exp10(log_density)


!**** Itoh 1996 Neutrino loss routines - Grant Newsham 9/06 *****


      if (star%ctrl%use_itoh_neutrino_loss) then
! ITOH ET AL. (1996) LOSSES (neutrino/neutrino_losses_for_mixture/sneut IN net_lib).
          call thermal_neutrino_losses(neutrino_temp,neutrino_density,hydrogen_fraction, &
               helium_fraction,carbon_fraction_total,oxygen_fraction_total, &
               neutrino_loss_snu,neutrino_dsnu_dt,neutrino_dsnu_drho)
          e%eps_neutrino_loss = -neutrino_loss_snu
          e%total_energy_gen_rate = e%total_energy_gen_rate + e%eps_neutrino_loss

! 2026 (bugsweep Batch 2): neutrino() returns snu >= 0 and the
! ABSOLUTE derivatives d snu/dT, d snu/d rho. deps_dlnt/dlnrho
! are the absolute d eps/d ln T, d eps/d ln rho accumulated in
! sum2/sum3 (see the qetnx/qednx terms of the Beaudet branch below),
! so the loss contributes d(-snu)/d ln T = -T d snu/dT. The original
! (engeb.f) first converted to dimensionless log-derivatives -- the
! wrong units for these sums -- and then built the T line from the
! already-updated rho line (PET = PEP + DSNUDT), discarding the
! nuclear sum3 altogether.
          e%deps_dlnrho = e%deps_dlnrho - neutrino_density*neutrino_dsnu_drho
          e%deps_dlnt = e%deps_dlnt - neutrino_temp*neutrino_dsnu_dt
      else
! LEGACY BEAUDET, PETROSIAN & SALPETER (1967) FITS FOR THE PAIR, PHOTO
!  AND PLASMA NEUTRINO LOSSES (COEFFICIENTS v1, v2, v3 ABOVE), WITH
!  THEIR D/D LN T AND D/D LN RHO ACCUMULATED IN qetnx/qednx.
         el = t9/5.9302d0
         eli = 1.0d0/el
         ez = exp(cc13*ln10*(dd-9.0d0))*eli*(0.7937d0+0.2063d0*hydrogen_fraction)
         ez3 = ez**3
         emue = 0.5d0*(1.0d0+hydrogen_fraction)
         ex1 = 0.0d0
         if (t9.ge.0.2d0) then
!C PAIR NEUTRINOS
            el2 = el*el
            polx10=(1.d0+el2*(-13.04d0+el2*(133.5d0+el2*(1534.d0+el2*918.6d0))))
            polx11 = v1(1) + ez*(v1(2) + ez*v1(3))
            polx12 = ez3 + eli*(v1(4) + eli*(v1(5) + eli*v1(6)))
            ex1 = exp(-ez*v1(7)-eli-eli-ln10*dd)*polx10*polx11/polx12
         end if
!C PHOTO NEUTRINOS
         polx21 = v2(1) + ez*(v2(2) + ez*v2(3))
         polx22 = ez3 + eli*(v2(4) + eli*(v2(5) + eli*v2(6)))
         ex2 = emue*el**5*exp(-ez*v2(7))*polx21/polx22
!C PLASMA NEUTRINOS
         polx31 = v3(1) + ez*(v3(2) + ez*v3(3))
         polx32 = ez3 + eli*(v3(4) + eli*(v3(5) + eli*v3(6)))
         ex3 = emue**3*exp(-ez*v3(7)+ln10*(dd+dd))*polx31/polx32
         e%eps_neutrino_loss = -(ex1 + ex2 + ex3)
         e%total_energy_gen_rate = e%total_energy_gen_rate + e%eps_neutrino_loss
         qetnx = 0.0d0
         qednx = 0.0d0
         if (t9.ge.0.2d0) then
! MHP 10/02 fixed column 72 problem
            qedn= ez*(v1(2)+2.d0*ez*v1(3))/polx11 - ez*v1(7) &
                  - 3.0d0*ez3/polx12
            qetn= -qedn+ eli*(v1(4)+eli*(2.d0*v1(5)+3.d0*eli*v1(6)))/polx12
            qetnx= el2*(-26.08d0+el2*(534.d0+el2*(9204.d0+el2*7348.8d0)))/polx10
            qetnx = (qetnx + qetn + eli+eli)*ex1
            qednx = (-1.0d0 +cc13*qedn)*ex1
         end if
         qedn= ez*(v2(2)+2.d0*ez*v2(3))/polx21 - ez*v2(7) - 3.d0*ez3/polx22
         qetn= -qedn+ eli*(v2(4)+eli*(2.d0*v2(5)+3.d0*eli*v2(6)))/polx22
         qetnx = qetnx + (5.0d0 + qetn)*ex2
         qednx = qednx +cc13*qedn*ex2
         qedn= ez*(v3(2)+2.d0*ez*v3(3))/polx31 - ez*v3(7) - 3.d0*ez3/polx32
         qetn=-qedn + eli*(v3(4)+eli*(2.d0*v3(5)+3.d0*eli*v3(6)))/polx32
         qetnx = qetnx + qetn*ex3
         qednx = qednx + (2.0d0 +cc13*qedn)*ex3
         e%deps_dlnt = e%deps_dlnt - qetnx
         e%deps_dlnrho = e%deps_dlnrho - qednx
      end if

end subroutine thermal_neutrino_loss

end subroutine nuclear_energy

end module nuclear_energy_lib
