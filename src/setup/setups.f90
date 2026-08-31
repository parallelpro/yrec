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
! (Fermi-Dirac) equation-of-state table used by fully_ionized_eos.f90, the
! Kurucz/Castelli surface-pressure table used for the T-tau surface
! boundary condition, and (if enabled) the SCVH envelope
! equation-of-state tables.
subroutine setups(ierr)
      use eos_lib
      use kap_lib
      use atm_lib
      use atm_table_lib
      use star_info_lib, only: star, json
      use phys_const_lib
      use yale_eos_lib
      use math_lib
      use ttau_lib, only: hsra_t_tau_offset
      implicit none
      integer, intent(out) :: ierr
! JNT 06/14 ADD NTC FOR KURUCZ/CASTELLI 2004 ATM
      integer, parameter :: nt = 57, ng = 11
      integer, parameter :: ntc = 76, ngc = 11
      integer, parameter :: nts = 63, nps = 76
! MHP 8/97 ADDED NTA AND NGA FOR ALLARD ATMOSPHERE
      integer, parameter :: nta = 54, nga = 5


! MHP 8/25 Reduced declared variables to ones actually used here or
! passed to other routine
! 2026 de-tramp: the 24 table paths and the mixture work array come
! from star%job directly (they were only ever star%job members passed
! positionally through star_setup) -- callers set star%job first.

! former common/lunum/: all 12 members now use-associated from
! const_lib (see that file's header note) rather than locally
! declared/common'd here.
! MHP 8/25 Removed file names from common block
!      COMMON/LUFNM/ FLAST, FFIRST, FRUN, FSTAND, FFERMI,
!     1    FDEBUG, FTRACK, FSHORT, FMILNE,  FMODPT,
!     2    FSTOR, FPMOD, FPENV, FPATM, FDYN,
!     3    FLLDAT, FSNU, FSCOMP, FKUR,
!     4    FMHD1, FMHD2, FMHD3, FMHD4, FMHD5, FMHD6, FMHD7, FMHD8
! --- locals ---
      double precision :: speed_of_light, electron_mass, boltzmann_constant, &
           planck_constant, hydrogen_atom_mass, electron_charge_esu
      integer :: teff_idx, logg_idx
      logical :: found_valid_pressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SETUP CONSTANTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      clndp = log(10.0d0)
      ln10 = clndp
      clni = 1.0d0/ln10
! Luminosity of Sun
      star%log10_solar_luminosity = log10(star%solar_luminosity_cgs)
      star%ln_solar_luminosity = ln10/star%solar_luminosity_cgs
! Mass of Sun (namelist control since 2026; default 1.9891d33)
      star%solar_mass_cgs = star%ctrl%solar_mass_cgs
      star%log10_solar_mass = log10(star%solar_mass_cgs)
! Radius of Sun
      star%log10_solar_radius = log10(star%solar_radius_cgs)
! Bolometric magnitude of the sun
      star%solar_bolometric_magnitude = 4.79d0
! No. of seconds per year and other mathematical constants
      seconds_per_year = 3.1558d7
      cc13 = 1.0d0/3.0d0
      cc23 = cc13 + cc13
      cpi = 3.1415926535898d0
      c4pi = 4.0d0*cpi
      c4pil = log10(c4pi)
      c4pi3l = log10(cc13*c4pi)
! Speed of Light
      speed_of_light = 2.99792458d10
! Stefan-Boltzmann constant
      csig = 5.67051d-5
      csigl = log10(csig)
! radiation constant/3
      radiation_constant_over_3 = csig*4.0d0/(3.0d0*speed_of_light)
      ca3l = log10(radiation_constant_over_3)
! molar gas constant
      gas_constant = 8.314510d7
! log of the gravitational constant. G_cgs < 0 (default) keeps the
! historical hard-coded log10 G = -7.17571 (G=6.6726D-8) bit-for-bit;
! a positive G_cgs (namelist control since 2026) overrides it.
      if (star%ctrl%G_cgs > 0.0d0) then
         cgl = log10(star%ctrl%G_cgs)
      else
         cgl = -7.17571d0
      end if
! mass of the electron
      electron_mass = 9.1093897d-28
! Boltzmann constant
      boltzmann_constant = 1.380658d-16
! Planck's constant
      planck_constant = 6.6260755d-27
      cmkh = 1.50d0*log10(2.0d0*cpi*electron_mass)+2.5d0* &
           log10(boltzmann_constant)- &
           3.0d0*log10(planck_constant)
      cdelrl = -c4pi3l - csigl - log10(16.0d0)
      cmixl2 = cc13
      cmixl3 = 16.0d0*dsqrt(2.0d0)*csig
! DBG CALCULATE DEBYE-HUCKEL CONSTAND CDH
! mass hydrogen atom (gm)
      hydrogen_atom_mass = 1.0d0/6.0222137d23
! charge electron (ESU)
      electron_charge_esu=4.802d-10
      debye_huckel_coefficient = -sqrt(cpi/(boltzmann_constant*hydrogen_atom_mass* &
           hydrogen_atom_mass*hydrogen_atom_mass))*electron_charge_esu* &
           electron_charge_esu*electron_charge_esu/3.0d0
! Electric charge of C,N,O,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,Ca,Ti,Cr,Mn,Fe,Ni
      debye_huckel_nu(1) = 6.0d0
      debye_huckel_nu(2) = 7.0d0
      debye_huckel_nu(3) = 8.0d0
      debye_huckel_nu(4) = 10.0d0
      debye_huckel_nu(5) = 11.0d0
      debye_huckel_nu(6) = 12.0d0
      debye_huckel_nu(7) = 13.0d0
      debye_huckel_nu(8) = 14.0d0
      debye_huckel_nu(9) = 15.0d0
      debye_huckel_nu(10) = 16.0d0
      debye_huckel_nu(11) = 17.0d0
      debye_huckel_nu(12) = 18.0d0
      debye_huckel_nu(13) = 20.0d0
      debye_huckel_nu(14) = 22.0d0
      debye_huckel_nu(15) = 24.0d0
      debye_huckel_nu(16) = 25.0d0
      debye_huckel_nu(17) = 26.0d0
      debye_huckel_nu(18) = 28.0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     EVALUATE TAU = 2/3 TEMPERATURE FOR HRA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      star%atm_hras = hsra_t_tau_offset()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     SET UP OPACITY TABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ierr = 0
      call kap_init(star%envelope_hydrogen_fraction, &
           star%envelope_metal_fraction, star%job%mixture_weights, &
           star%job%alex06_table_path,star%job%kurucz_table_path,star%job%kurucz_table2_path, &
           star%job%laol_table_path,star%job%laol_table2_path, &
           star%job%opal95_table_path,star%job%opal92_table_path,star%job%opal92_table2_path, &
           star%job%pure_z_table_path,star%job%alex95_table_paths, ierr)
      if (ierr /= 0) return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     READ IN MHD EOS TABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MHP 8/25 Passed file names directly instead of through common, as
! character strings in common blocks can cause problems
! 2026 (ROADMAP.md stage 1): the use_mhd_eos gate moved inside
! eos_lib's eos_init; a no-op when MHD is off, exactly as before.
! eos_init also reads the Fermi-Dirac degenerate-electron table
! (the F-tables block that used to sit inline just below).
      call eos_init(star%job%fermi_table_path,star%job%scv_h_table_path,star%job%scv_he_table_path, &
           star%job%scv_z_table_path,star%job%zams_a_table_path,star%job%zams_b_table_path, &
           star%job%zams_c_table_path,star%job%centre1_table_path,star%job%centre2_table_path, &
           star%job%centre3_table_path,star%job%centre4_table_path,star%job%centre5_table_path, ierr)
      if (ierr /= 0) return
!
!
! (The F-tables read for the degenerate equation of state that lived
! here moved into eos_lib's eos_init above -- 2026, ROADMAP.md
! stage 1 -- along with its state's relocation from atm_table_lib to
! yale_eos_lib.)
! JMH 8/18/91
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INPUT PRESSURE TABLE FOR SURFACE BOUNDARY CONDITIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2026 (ROADMAP.md stage 1): the Kurucz (atm_choice 3/4), Allard
! (atm_choice 4), and Kurucz/Castelli (atm_choice 5) surface table
! loads that lived inline here moved into atm_lib's atm_init.
      call atm_init(star%job%atm_table_path,star%job%allard_table_path, ierr)
      if (ierr /= 0) return

! (The SCV equation-of-state table reads that lived here moved into
! eos_lib's eos_init -- 2026, ROADMAP.md stage 1 -- alongside the
! Fermi-Dirac table read.)

      return
end subroutine setups
