!----------------------------------------------------------------------
! wrthead
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original wrthead.f; only variable names, source form, and comment
! style were updated.
!
! Write the headers for all the appropriate output files.
subroutine wrthead(total_mass_msun)

      use const_lib
      use luout_lib
      use star_info_lib, only: star
      implicit none

      double precision, intent(in) :: total_mass_msun
! --- locals ---
      double precision :: total_mass_grams

      if (rescale_kind(nk) .eq. 1) then
         write(iowr, 47) nk, initial_envelope_x, initial_envelope_z, &
              star%mixing_length_alpha, num_models(nk)
      else if (rescale_kind(nk) .eq. 2) then
         write(iowr, 48) nk, initial_envelope_x, initial_envelope_z, &
              star%mixing_length_alpha, num_models(nk)
      else if (rescale_kind(nk) .eq. 3) then
         write(iowr, 49) nk, initial_envelope_x, initial_envelope_z, &
              star%mixing_length_alpha, num_models(nk)
      end if
  47  format(/, ' RUN=',I2,' EVOLVE  ', ' X=',F8.6, &
             ' Z=',F8.6,' CMIXL=', F8.6, ' NO.MODS=', I5)
  48  format(/, ' RUN=',I2,' RESCALE ', ' X=',F8.6, &
             ' Z=',F8.6,' CMIXL=', F8.6, ' NO.MODS=', I5)
  49  format(/, ' RUN=',I2,' RESCALE&EVOLVE ', ' X=',F8.6, &
             ' Z=',F8.6,' CMIXL=', F8.6, ' NO.MODS=', I5)

      if (isochrone_output_active) then
! header stuff for isochrone output
         total_mass_grams = total_mass_msun*star%solar_mass_cgs
         write(isochrone_file_unit, 1495) total_mass_grams, &
              initial_envelope_x,initial_envelope_z,star%mixing_length_alpha,star%solar_bolometric_magnitude
 1495    format(7X, 1P5E16.8)
      end if

      if (ltrack .and. first_call_flag(nk)) then
! ITRVER identifies version of track out file.  If you change
! the track out file then change this version number.
         write(itrack, 1500)track_file_version,total_mass_msun,initial_envelope_x, &
              initial_envelope_z,star%mixing_length_alpha
 1500    format('#Version=',i3,'  Mtot/Msun =',1PE16.8, &
              '  Initial: X =',1PE16.8,' Z =',1PE16.8, &
              '  Mix. length =', 1PE16.8)
         if(track_file_version .eq. 0) then
!            WRITE(ITRACK, 1503)
! 1503       FORMAT(
!     1'# Model #, shells, AGE(Gyr), log(L/Lsun), log(R/Rsun), log(g),',
!     1' log(Teff), Mconv. core (Msun), Mconv. env, R,T,Rho,P,cappa env',/,
!     2'# Central: log(T), log(RHO), log(P), BETA, ETA, X, Y, Z',/,
!     3'# Luminosity: ppI, ppII, ppIII, CNO, triple-alpha,',
!     3' He-C, gravity, neutrinos (old)',/,
!     3'# Cl SNU, Ga SNU, Neutrinos (1E10 erg/CM^2 at earth): pp, pep, hep, Be7,',
!     3' B8, N13, O15, F17 2xdiagnostic',/,
!     4'# Central Abundances: He3, C12, C13, N14, N15, O16,',
!     4' O17, O18',/,
!     5'# Surface Abundances: He3, C12, C13, N14, N15, O16,',
!     5' O17,O18 H2, Li6, Li7, Be9 X Y Z Z/X',/,
!     6'# Jtot, KE rot tot, total I, CZ I, Omega center, surface, Prot (day), Vrot (km/s), TauCZ (sec) ',/,
!     7'# H shell loc: mass frac-base, midpoint, top; radius frac-',
!     7'base, midpoint, top, Pphot, mass (msun)')
! G Somers 11/14; Added option to create .track file header. Uncomment the following
! block if desired.
            write(itrack, 1504)
! 1504       FORMAT(
!     1'# ',/,
!     2'#    Step    Shls     Age (Gyr)      log(L/Lsun)    log(R/Rsun)        log(g)        log(Teff)',
!     2'       Mconv.core      Mconv.env     Rcore     Tcore      Rho_core     P_core     kappa_env   ',
!     2' log(T)_cen     log(Rho)_cen     log(P)_cen        BETA             ETA            X_cen      ',
!     2'     Y_cen           Z_cen          ppI_lum         ppII_lum       ppIII_lum        CNO_lum   ',
!     2'    3-alpha_lum       He-C_lum        gravity       OLD NUTRINOS   Cl SNU    Ga SNU    **pp** ',
!     2'   **pep**   **hep**   **Be7**   **B8**    **N13**   **O15**   **F17**  **diag1** **diag2**   ',
!     2'   He3_cen        C12_cen         C13_cen         N14_cen         B10_cen         O16_cen     ',
!     2'    B11_cen         O18_cen         He3_surf        C12_surf        C13_surf        N14_surf  ',
!     2'      B10_surf        O16_surf        B11_surf        O18_surf        H2_surf         Li6_surf',
!     2'        Li7_surf        Be9_surf         X_surf          Y_surf          Z_surf         Z/X_su',
!     2'rf          Jtot         KE_rot_tot       total I           CZ I         Omega_surf      Omega',
!     2'_cen       Prot (days)     Vrot (km/s)     TauCZ (s)       Mfrac_base      Mfrac_midp      Mfr',
!     2'ac_top       Rfrac_base      Rfrac_midp      Rfrac_top         P_phot           Mass     ',/,
!     2'# ')
 1504       format( &
     '# ',/, &
     '     Step    Shls         Age_gyr       LogL_lsun       LogR_rsun           Log_g', &
     '        log_Teff        Mco_core         Mco_env     Rco_env     Tco_env     Dco_env', &
     '     Pco_env     Oco_env        LogT_cen        LogD_cen        logP_cen        Beta_cen', &
     '         Eta_cen           X_cen           Y_cen           Z_cen        ppI_lsun', &
     '       ppII_lsun      ppIII_lsun        CNO_lsun         3a_lsun        HeC_lsun', &
     '      Egrav_lsun       Neut_lsun    Cl_snu    Ga_snu   pp_neut  pep_neut  hep_neut', &
     '  Be7_neut   B8_neut  N13_neut  O15_neut  F17_neut     diag1     diag2         He3_cen', &
     '         C12_cen         C13_cen         N14_cen         N15_cen         O16_cen', &
     '         O17_cen         O18_cen         He3_sur         C12_sur         C13_sur', &
     '         N14_sur         N15_sur         O16_sur         O17_sur         O18_sur', &
     '          H2_sur         Li6_sur         Li7_sur         Be9_sur           X_sur', &
     '           Y_sur           Z_sur         Z_X_sur            Jtot      KE_rot_tot', &
     '           I_tot            I_cz       Omega_sur       Omega_cen      Prot_sur_d', &
     '        Vrot_kms         TauCZ_s    MHshell_base     MHshell_mid     MHshell_top', &
     '    RHshell_base     RHShell_mid     RHshell_top       logP_phot       Mass_msun',/, &
     '# ')
! G Somers END.
         else if(track_file_version .eq. 1) then
            write(itrack, 1505)
 1505       format( &
     '# Model #, shells, AGE, log(L/Lsun), log(R/Rsun), log(g),', &
     ' log(Teff), Mconv. core, Mconv. env.)' ,/, &
     '# Central: log(T), log(RHO), log(P), BETA, ETA, X, Y, Z',/, &
     '# Luminosity: ppI, ppII, ppIII, CNO, triple-alpha,', &
     ' He-C, gravity, neutrinos (old)',/, &
     ' Neutrinos (1E10 erg/CM^2 at earth): pp, pep, hep, Be7,', &
     ' B8, N13, O15, F17',/, &
     '# Central Abundances: He3, C12, C13, N14, N15, O16,', &
     ' O17, O18',/, &
     '# Surface Abundances: He3, C12, C13, N14, N15, O16,', &
     ' O17,O18',/, &
     '    "        " cont: H2, Li6, Li7, Be9',/, &
     '# H shell loc: mass frac-base, midpoint, top; radius frac-', &
     'base, midpoint, top')
         else if(track_file_version .eq. 2) then
            write(itrack, 1505)
            write(itrack, 1510)
 1510       format( &
     ' Jtot, K.E. Rotation, OMEGAsurf, OMEGAcenter')
         else if(track_file_version .eq. 3) then
          write(itrack, 1515)
 1515       format( &
     '#Model #, shells, AGE, log(L/Lsun), log(R/Rsun), log(g),', &
     ' log(Teff), Mconv. core, Mconv. env., % Grav Energy, X env')
         end if
      end if

      return
end subroutine wrthead
