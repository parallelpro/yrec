!----------------------------------------------------------------------
! wrtmod
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original wrtmod.f; only variable names, source form, and comment
! style were updated.
!
! Writes the sound-speed diagnostic profile (if requested), integrates
! a fresh envelope/atmosphere from the surface down to the converged
! model for printout, and writes the pulsation model (.pmod-style)
! records for each requested point.
!
!       SUBROUTINE WRTMOD(M,LSHELL,JXBEG,JXEND,JCORE,JENV,HCOMP,HS1,HD,
!      *HL,HP,HR,HT,LC,MODEL,BL,TEFFL,OMEGA,FP,FT,ETA2,R0,HJM,HI,HS,
!      * DAGE)  ! KC 2025-05-31
subroutine wrtmod(num_shells, envelope_cz_bottom_index, composition, &
     mass_coordinate, log_density, log_luminosity, log_pressure, log_radius, &
     log_temperature, model_number, log_luminosity_lsun, log_teff, &
     shape_factor_fp, shape_factor_ft, log_mass, age_gyr)

      use luout_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: num_shells, envelope_cz_bottom_index
      double precision, intent(in) :: composition(15,json)
      double precision, intent(in) :: mass_coordinate(json), log_density(json), &
           log_luminosity(json), log_pressure(json), log_radius(json), &
           log_temperature(json)
      integer, intent(in) :: model_number
      double precision, intent(in) :: log_luminosity_lsun, log_teff
      double precision, intent(in) :: shape_factor_fp(json), shape_factor_ft(json)
      double precision, intent(in) :: log_mass(json)
      double precision, intent(in) :: age_gyr


! common/lunum/: only opal_model_unit (IOPMOD) is used here. Naming
! matches setkrz.f90.
      integer :: first_unit, run_unit, standard_unit, fermi_unit, &
           opal_model_unit, opal_envelope_unit, opal_atm_unit, &
           dynamics_unit, laol_table_unit, neutrino_unit, &
           composition_unit, kurucz_table_unit
      common/lunum/ first_unit, run_unit, standard_unit, fermi_unit, &
           opal_model_unit, opal_envelope_unit, opal_atm_unit, &
           dynamics_unit, laol_table_unit, neutrino_unit, &
           composition_unit, kurucz_table_unit

! common/label/: initial_envelope_x/initial_envelope_z, used here.
! Naming matches wrthead.f90.
      double precision :: initial_envelope_x, initial_envelope_z
      common/label/ initial_envelope_x, initial_envelope_z

! common/ccout/: only lstore/lstatm/lstenv are used here. Naming
! matches ccoeft.f90.
      logical :: lstore, lstatm, lstenv, lstmod, lstphys, lstrot, lscrib, &
           lstch, lphhd
      common/ccout/ lstore, lstatm, lstenv, lstmod, lstphys, lstrot, &
           lscrib, lstch, lphhd

! common/ccout1/: not used in this file. Naming matches wrtmil.f90.
      integer :: npenv, nprtmod, print_point_interval, npoint
      common/ccout1/ npenv, nprtmod, print_point_interval, npoint

! common/ccout2/: not used in this file. Naming matches meqos.f90.
      logical :: ldebug, lcorr, lmilne, ltrack, lstpch
      common/ccout2/ ldebug, lcorr, lmilne, ltrack, lstpch

! common/comp/: only stotal is used here. Naming matches getopac.f90.
      double precision :: envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, amuenv, fxenv(12), xnew, znew, stotal, senv
      common/comp/ envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, amuenv, fxenv, xnew, znew, stotal, senv

! common/const/: only clsunl/crsunl are used here. Naming matches
! vcirc.f90.
      double precision :: solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude
      common/const/ solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude



! common/const3/: only mixing_length is used here. Naming matches
! mix.f90.
      double precision :: cdelrl, mixing_length, cmixl2, cmixl3, clndp, &
           seconds_per_year
      common/const3/ cdelrl, mixing_length, cmixl2, cmixl3, clndp, &
           seconds_per_year

! common/envgen/: atmosphere/envelope step sizes for printout, used
! here. Naming is local to this batch.
      double precision :: atm_step_size, envelope_step_size
      logical :: envelope_generation_flag
      common/envgen/ atm_step_size, envelope_step_size, envelope_generation_flag

! common/flag/: not used in this file. Naming matches mixcz.f90.
      logical :: use_extended_composition
      common/flag/ use_extended_composition

! common/intatm/: atmosphere integration control, all used here.
! Naming is local to this batch.
      double precision :: atm_error_tol, atm_step_initial, atm_step_begin, &
           atm_step_min, atm_step_max
      common/intatm/ atm_error_tol, atm_step_initial, atm_step_begin, &
           atm_step_min, atm_step_max

! common/intenv/: envelope integration control, all used here. Naming
! is local to this batch.
      double precision :: env_error_tol, env_step_begin, env_step_min, &
           env_step_max
      common/intenv/ env_error_tol, env_step_begin, env_step_min, env_step_max

! common/rot/: not used in this file. Naming matches momi.f90.
      double precision :: wnew, walpcz, acfpft
      integer :: itfp1, itfp2
      logical :: rotation_active, instability_transport_active, lwnew
      common/rot/ wnew, walpcz, acfpft, itfp1, itfp2, rotation_active, &
           instability_transport_active, lwnew

! common/scrtch/: only sesum/seg/so/del_grad are used here. Naming
! matches microdiff_setup.f90.
      double precision :: sesum(json), seg(7,json), sbeta(json), seta(json)
      logical :: locons(json)
      double precision :: so(json), del_grad(3,json), sfxion(3,json), &
           svel(json), scp(json)
      common/scrtch/ sesum, seg, sbeta, seta, locons, so, del_grad, &
           sfxion, svel, scp

! DBG PULSE
! common/pulse/: pulsation_mass_msun/pulsation_output_active/
! pulsation_file_version, all used here. Naming is local to this
! batch.
      double precision :: pulsation_mass_msun
      logical :: pulsation_output_active
      integer :: pulsation_file_version
      common/pulse/ pulsation_mass_msun, pulsation_output_active, &
           pulsation_file_version

! common/pulse1/: per-shell partial derivatives used to build the
! pulsation model output. The exact physical definition of each is not
! confidently known from this file alone (they parallel the eqstat
! QDP/QDT/QCP/RMU/EMU family, saved per shell for pulsation output);
! names below are conservative guesses, flagged accordingly. lpumod is
! an unused placeholder.
      double precision :: pulse_dlnrho_dlnp(json), pulse_dlneps_dlnrho(json), &
           pulse_dlneps_dlnt(json), pulse_dlnkap_dlnrho(json), &
           pulse_dlnkap_dlnt(json), pulse_specific_heat(json), &
           pulse_mean_molecular_weight(json), pulse_dlnrho_dlnt(json), &
           pulse_electron_mean_molecular_weight(json)
      logical :: lpumod
      common/pulse1/ pulse_dlnrho_dlnp, pulse_dlneps_dlnrho, &
           pulse_dlneps_dlnt, pulse_dlnkap_dlnrho, pulse_dlnkap_dlnt, &
           pulse_specific_heat, pulse_mean_molecular_weight, &
           pulse_dlnrho_dlnt, pulse_electron_mean_molecular_weight, lpumod

! common/pulse2/: not used in this file; declared only to preserve
! layout. Naming is local to this batch (kept close to the original
! cryptic names since the members are unused here).
      double precision :: qqdp, qqed, qqet, qqod, qqot, qdel, qdela, qqcp, &
           qrmu, qtl, qpl, qdl, qo, qol, qt, qp, qqdt, qemu, qd, qfs
      common/pulse2/ qqdp, qqed, qqet, qqod, qqot, qdel, qdela, qqcp, qrmu, &
           qtl, qpl, qdl, qo, qol, qt, qp, qqdt, qemu, qd, qfs

! MHP 7/96 COMMON BLOCK INSERTED FOR SOUND SPEED
! common/sound/: adiabatic_index_gamma1/sound_speed_output_active,
! both used here. Naming matches putstore.f90.
      double precision :: adiabatic_index_gamma1(json)
      logical :: sound_speed_output_active
      common/sound/ adiabatic_index_gamma1, sound_speed_output_active
! DBG
! DBG 7/92 COMMON BLOCK ADDED TO COMPUTE DEBYE-HUCKEL CORRECTION.
! common/debhu/: only ldh/xxdy(=XXDH)/yydh/zzdh/zdh are used here.
! Naming matches eqstat.f90/wrtlst.f90 (xxdy is the converted name for
! the original XXDH slot).
      double precision :: cdh, etadh0, etadh1, zdh(18), xxdy, yydh, zzdh, &
           dhnue(18)
      logical :: ldh
      common/debhu/ cdh, etadh0, etadh1, zdh, xxdy, yydh, zzdh, dhnue, ldh
! DBG 7/95 To store variables for pulse output
! common/pualpha/: valfmlt/vphmlt/vcmxmlt are used here. Naming
! matches tpgrad.f90.
      double precision :: alfmlt, phmlt, cmxmlt
      double precision :: valfmlt(json), vphmlt(json), vcmxmlt(json)
      common/pualpha/ alfmlt, phmlt, cmxmlt, valfmlt, vphmlt, vcmxmlt
! common/roten/: not used in this file. Naming matches putstore.f90.
      double precision :: rotational_energy_term(json)
      common/roten/ rotational_energy_term

! G Somers 10/14, Add spot common block, and store common block.
! common/spots/: spot_filling_factor/spot_temp_contrast, used here.
! Naming matches tpgrad.f90.
      double precision :: spot_filling_factor, spot_temp_contrast
      logical :: spot_depth_varies
      common/spots/ spot_filling_factor, spot_temp_contrast, spot_depth_varies
! common/temp2/: not used in this file. Naming matches vcirc.f90.
      double precision :: es_circulation_velocity(json), &
           es_circulation_velocity_prev(json), secular_shear_velocity(json), &
           secular_shear_velocity_prev(json), hle(json), &
           gsf_circulation_velocity(json), gsf_circulation_velocity_prev(json), &
           mu_gradient_velocity(json)
      common/temp2/ es_circulation_velocity, es_circulation_velocity_prev, &
           secular_shear_velocity, secular_shear_velocity_prev, hle, &
           gsf_circulation_velocity, gsf_circulation_velocity_prev, &
           mu_gradient_velocity
! common/quadd/: not used in this file. Naming matches vcirc.f90.
      double precision :: phisp(json), phirot(json), phidis(json), &
           circulation_correction_ratio(json)
      common/quadd/ phisp, phirot, phidis, circulation_correction_ratio
! G Somers END

      double precision :: dum1(4), dum2(3), dum3(3), dum4(3)
! KC 2025-05-30 addressed warning messages from Makefile.legacy
!      DATA IHEADR/4H****/

! pulsation_point_id/num_pulsation_points (ID/IDM) are never assigned
! anywhere in this subroutine (no CALL PINDEX or equivalent here,
! unlike putstore.f90/wrtlst.f90). With the implicit SAVE below they
! keep their static-storage default (in practice 0 on essentially all
! platforms), so "DO 220 J = 1,IDM" never executes and the WRITE(
! IOPMOD,5001) statement always reports IDM=0. This looks like a
! latent bug (perhaps a lost CALL PINDEX from an earlier refactor),
! but it is preserved exactly as in the original wrtmod.f; not "fixed"
! here.
      integer :: pulsation_point_id(json)
      integer :: num_pulsation_points

      save

! --- locals ---
      integer :: i
      double precision :: fr, fm, xx1, xxx, xx2, xx3, xx4, xx5, sound_velocity
      double precision :: rmid, dr, divp, divr, qpr1, qdr1, qpr2, qdr2, &
           qqpr, qqdr
      logical :: lsbc0, lprt
      integer :: katm, kenv, ksaha
      double precision :: abeg0, amin0, amax0, ebeg0, emin0, emax0
      double precision :: b, rl, gl, x, z, fpl, ftl
      integer :: ixx
      double precision :: hstot, plim
      logical :: lpulpt
      integer :: idum
      double precision :: ateffl
      double precision :: fsi
      double precision :: rsurfl, tempr, qsmass
      integer :: j
      double precision :: fs, pelpf

         print*,'wrtmod LSOUND 1: ',sound_speed_output_active

      if(sound_speed_output_active)then

         print*,'wrtmod LSOUND 2: ',sound_speed_output_active

!CFD 10/09 Add an extra output to plot the sound speed profile easyly
         open(unit=500,file='Csound.dat',status='unknown')
         write(500,9124)
 9123    format(1X,I6,1X,2F12.8,1P,3E16.8)
 9124    format(1X,'Shell#',5x,'R/Rsun',7x,'M/Msun',7x,'Cs',15x,'Rho')
!CFD end
         do i = 1,num_shells
            fr = exp(ln10*log_radius(i))/solar_radius_cgs
            fm = exp(ln10*log_mass(i))/solar_mass_cgs
            xx1 = fm/fr**3
            xxx = -exp(ln10*(cgl+log_mass(i)+log_density(i)-log_pressure(i)-log_radius(i)))
            xx2 = -xxx/adiabatic_index_gamma1(i)
            xx3 = adiabatic_index_gamma1(i)
            xx4 = -xx2-xxx*(pulse_dlnrho_dlnp(i)+del_grad(2,i)*pulse_dlnrho_dlnt(i))
            xx5 = exp(ln10*(c4pil+log_density(i)+3.0D0*log_radius(i)-log_mass(i)))
            sound_velocity = 1.0D-5*sqrt(adiabatic_index_gamma1(i)*exp(ln10*(log_pressure(i)-log_density(i))))
            write(imodpt,123)fr,fm,xx1,xx2,xx3,xx4,xx5,sound_velocity
 123        format(1X,2F12.8,1P6E16.8)
!CFD 10/09 Add an extra output to plot the sound speed profile easyly
            write(500,9123)i,fr,fm,sound_velocity,log_density(i),log_temperature(i)
!CFD end

!  DERIVATIVES OF DP/DR, DRHO/DR
            if(i.lt.num_shells) then
              rmid = 0.5D0*(exp(ln10*log_radius(i))+exp(ln10*log_radius(i+1)))
              dr = exp(ln10*log_radius(i+1)) - exp(ln10*log_radius(i))
              divp = 0.5D0*(adiabatic_index_gamma1(i)*exp(ln10*log_pressure(i))+ &
                    adiabatic_index_gamma1(i+1)*exp(ln10*log_pressure(i+1)))*dr
              divr = 0.5D0*(exp(ln10*log_density(i))+exp(ln10*log_density(i+1)))*dr
              qpr1 = exp(ln10*(cgl+log_mass(i)+log_density(i)-2.0D0*log_radius(i)))
              qdr1 = exp(ln10*(log_density(i)-log_pressure(i)))*pulse_dlnrho_dlnp(i)*qpr1
              qpr2 = exp(ln10*(cgl+log_mass(i+1)+log_density(i+1)-2.0D0*log_radius(i+1)))
              qdr2 = exp(ln10*(log_density(i+1)-log_pressure(i+1)))*pulse_dlnrho_dlnp(i+1)*qpr1
              qqpr = (qpr1-qpr2)/divp
              qqdr = (qdr1-qdr2)/divr
              write(imodpt,124)rmid/solar_radius_cgs,qpr1,qdr1,qqpr,qqdr
 124          format(1X,F11.7,1P4E16.8)
            endif
         end do
      endif
!
! G Somers 11/14 LCHEMO BLOCK REMOVED, AS THIS INFO IS ALREADY IN .STORE.
!
! DBG PULSE: PRINT OUT PULSATION ENV AND ATM IN ENVINT
! G Somers 11/14 CHANGE TO NEW I/O FLAGS.
      if(lstatm .or. lstenv .or. pulsation_output_active) then
!  INTEGRATE AN ENVELOPE FROM THE SURFACE TO THE CONVERGED MODEL,
!  PRINTING OUT THE RESULTS.
!  SET UP FLAGS AND COUNTERS.
       lsbc0 = .false.
       if(lstore)lprt = .true.
       katm = 0
       kenv = 0
       ksaha = 0
!  SAVE THE INTEGRATION STEP PARAMETERS AND ENFORCE THE SPACING
!  REQUESTED FOR PRINTOUT PURPOSES.
       abeg0 = atm_step_begin
       amin0 = atm_step_min
       amax0 = atm_step_max
       ebeg0 = env_step_begin
       emin0 = env_step_min
       emax0 = env_step_max
       atm_step_begin = atm_step_size
       atm_step_min = atm_step_size
       atm_step_max = atm_step_size
       env_step_begin = envelope_step_size
       env_step_min = envelope_step_size
       env_step_max = envelope_step_size
       b = dexp(ln10*log_luminosity_lsun)
       rl = 0.5D0*(log_luminosity_lsun + log10_solar_luminosity - 4.0D0*log_teff - c4pil - csigl)
       gl = cgl + stotal - rl - rl
       x = composition(1,num_shells)
       z = composition(3,num_shells)
       fpl = shape_factor_fp(num_shells)
       ftl = shape_factor_ft(num_shells)
       ixx=0
       hstot = stotal
       plim = log_pressure(num_shells)
! DBG PULSE: ADDED ARGUEMENT TO ENVINT TO TURN ON/OFF PULSE OUTPUT
         lpulpt = pulsation_output_active
            if (ldh) then
               xxdy = composition(1,num_shells)
               yydh = composition(2,num_shells)+composition(4,num_shells)
               zzdh = composition(3,num_shells)
               zdh(1) = composition(5,num_shells)+composition(6,num_shells)
               zdh(2) = composition(7,num_shells)+composition(8,num_shells)
               zdh(3) = composition(9,num_shells)+composition(10,num_shells)+composition(11,num_shells)
            end if
! MHP 10/02  define ISTORE - used in ENVINT
         idum = 0
! G Somers 10/14, FOR SPOTTED RUNS, FIND THE
! PRESSURE AT THE AMBIENT TEMPERATURE ATEFFL
        if(envelope_cz_bottom_index.eq.num_shells.and.spot_filling_factor.ne.0.0.and.spot_temp_contrast.ne.1.0)then
            ateffl = log_teff - 0.25*log10(spot_filling_factor * spot_temp_contrast**4.0 + 1.0 - spot_filling_factor)
       else
          ateffl = log_teff
       endif
       call envint(b,fpl,ftl,gl,hstot,ixx,lprt,lsbc0,plim,rl, &
                     ateffl,x,z,dum1,idum,katm,kenv,ksaha,dum2, &
                     dum3,dum4,lpulpt)
! G Somers END
       atm_step_begin = abeg0
       atm_step_min = amin0
       atm_step_max = amax0
       env_step_begin = ebeg0
       env_step_min = emin0
       env_step_max = emax0
      endif
!
! G Somers 11/14 LCONZO (convection zone info) block deleted.
!
! G Somers 11/14 LJOUT (rotation info) block deleted.
!
      fsi = dexp(-ln10*stotal)
! DBG PULSE: WRITE HEADER INFORMATION FOR PULSE MODEL
      if(pulsation_output_active) then
         rsurfl = 0.5D0*(log_luminosity_lsun - c4pil - csigl - 4.0D0*log_teff + log10_solar_luminosity)
         tempr = rsurfl - log10_solar_radius
         qsmass = pulsation_mass_msun
         write (opal_model_unit, 5001) model_number,num_pulsation_points,pulsation_file_version,qsmass, &
               log_teff,log_luminosity_lsun,tempr, age_gyr, mixing_length, initial_envelope_x, initial_envelope_z
 5001    format(' MODEL#=', I5, '  NUMBER OF SHELLS IN MODEL=',I5, &
                ' VER=',I2,/, &
                ' MASS=',F8.5, '  LOG(TEFF)=',F8.5,/, ' LOG(L/LSUN)=', &
                F16.10, '  LOG(R/RSUN)=',F16.10, /, &
                ' AGE=', 1PE12.5,' GYR',/, &
                ' MIXING LENGTH PARAMETER=', 0PF16.10,/, &
                ' ZAMS (X,Z)=', 2F18.10)
      end if
      do 220 j = 1,num_pulsation_points
       i = pulsation_point_id(j)
       fs = fsi*mass_coordinate(i)
!
! G Somers 11/14 FINAL LSCRIB BLOCK DELETED
!
! DBG WRITE PULSE MODEL
!       PRINT*, 'LPULSE=',LPULSE
         if (pulsation_output_active.and.lstore) then
! MHP 10/02 uncommented pelpf statement, used later in i/o
!         PELPF = CGAS * DEXP(CLN*(HT(I) + HD(I)))* PEMU(I)
!          ADDED X AND Z TO OUTPUT
         if ((j.eq.2).and.(i.eq.1)) goto 5003
         if(pulsation_file_version.eq.1) then
         pelpf = gas_constant * dexp(ln10*(log_temperature(i) + log_density(i)))* pulse_electron_mean_molecular_weight(i)
         write(opal_model_unit, 5052)log_radius(i),fs,log_luminosity(i),log_temperature(i),log_density(i), &
                log_pressure(i), sesum(i),so(i), pulse_dlnrho_dlnp(i), pulse_dlneps_dlnrho(i), &
                pulse_dlneps_dlnt(i), pulse_dlnkap_dlnrho(i), pulse_dlnkap_dlnt(i), del_grad(2,i),del_grad(3,i), &
                pulse_specific_heat(i), pulse_mean_molecular_weight(i), pulse_dlnrho_dlnt(i), pelpf
         else if (pulsation_file_version.eq.2) then
         write(opal_model_unit, 6052)log_radius(i),fs,log_luminosity(i),log_temperature(i),log_density(i), &
            log_pressure(i), sesum(i),so(i), pulse_dlnrho_dlnp(i), pulse_dlneps_dlnrho(i), &
            pulse_dlneps_dlnt(i), pulse_dlnkap_dlnrho(i), pulse_dlnkap_dlnt(i), del_grad(2,i),del_grad(3,i), &
            pulse_specific_heat(i), pulse_mean_molecular_weight(i), pulse_dlnrho_dlnt(i), composition(1,i),composition(3,i)
         else if (pulsation_file_version.eq.3) then
! DBG 7/95 Modifed to include mixing length variables
         write(opal_model_unit, 6053)log_radius(i),fs,log_luminosity(i),log_temperature(i),log_density(i),valfmlt(i), &
            log_pressure(i), sesum(i),so(i), pulse_dlnrho_dlnp(i), pulse_dlneps_dlnrho(i),vphmlt(i), &
            pulse_dlneps_dlnt(i), pulse_dlnkap_dlnrho(i), pulse_dlnkap_dlnt(i), del_grad(2,i),del_grad(3,i),vcmxmlt(i), &
            pulse_specific_heat(i), pulse_mean_molecular_weight(i), pulse_dlnrho_dlnt(i), composition(1,i),composition(3,i)
         end if
 5003      continue
 5052      format(5E16.9,/,5E16.9,/,5E16.9,/,5E16.9)
 6052      format(5E23.16,/,5E23.16,/,5E23.16,/,5E23.16)
 6053      format(6E23.16,/,6E23.16,/,6E23.16,/,5E23.16)
      end if
! DBG END
  220 continue
!
! G Somers 11/14 REMOVED LONG BLOCK
!
      return
end subroutine wrtmod
