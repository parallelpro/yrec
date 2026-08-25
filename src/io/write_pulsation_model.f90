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
subroutine write_pulsation_model(num_shells, envelope_cz_bottom_index, composition, &
     mass_coordinate, log_density, log_luminosity, log_pressure, log_radius, &
     log_temperature, model_number, log_luminosity_lsun, log_teff, &
     shape_factor_fp, shape_factor_ft, log_mass, age_gyr)

      use atm_lib
      use envint_lib, only: atm_get
      use star_info_lib, only: star, i_grad_actual, i_grad_ad, json
      use luout_lib
      use phys_const_lib
      use yale_eos_lib

      implicit none

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






















! G Somers END

      double precision :: dum1(4), dum2(3), dum3(3), dum4(3)
! KC 2025-05-30 addressed warning messages from Makefile.legacy
!      DATA IHEADR/4H****/

! pulsation_point_id/num_pulsation_points (ID/IDM) are never assigned
! anywhere in this subroutine (no CALL PINDEX or equivalent here,
! unlike write_store_model.f90/write_last_model.f90). With the implicit SAVE below they
! keep their static-storage default (in practice 0 on essentially all
! platforms), so "DO 220 J = 1,IDM" never executes and the WRITE(
! IOPMOD,5001) statement always reports IDM=0. This looks like a
! latent bug (perhaps a lost CALL PINDEX from an earlier refactor),
! but it is preserved exactly as in the original wrtmod.f; not "fixed"
! here.
      integer :: pulsation_point_id(json)
      integer :: num_pulsation_points
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

         print*,'wrtmod LSOUND 1: ',star%sound_speed_output_active

      if(star%sound_speed_output_active)then

         print*,'wrtmod LSOUND 2: ',star%sound_speed_output_active

!CFD 10/09 Add an extra output to plot the sound speed profile easyly
         open(unit=500,file='Csound.dat',status='unknown')
         write(500,9124)
 9123    format(1X,I6,1X,2F12.8,1P,3E16.8)
 9124    format(1X,'Shell#',5x,'R/Rsun',7x,'M/Msun',7x,'Cs',15x,'Rho')
!CFD end
         do i = 1,num_shells
            fr = exp(ln10*log_radius(i))/star%solar_radius_cgs
            fm = exp(ln10*log_mass(i))/star%solar_mass_cgs
            xx1 = fm/fr**3
            xxx = -exp(ln10*(cgl+log_mass(i)+log_density(i)-log_pressure(i)-log_radius(i)))
            xx2 = -xxx/star%adiabatic_index_gamma1(i)
            xx3 = star%adiabatic_index_gamma1(i)
            xx4 = -xx2-xxx*(star%pulse_dlnrho_dlnp(i)+star%del_grad(i_grad_actual,i)*star%pulse_dlnrho_dlnt(i))
            xx5 = exp(ln10*(c4pil+log_density(i)+3.0D0*log_radius(i)-log_mass(i)))
            sound_velocity = 1.0D-5*sqrt(star%adiabatic_index_gamma1(i)*exp(ln10*(log_pressure(i)-log_density(i))))
            write(imodpt,123)fr,fm,xx1,xx2,xx3,xx4,xx5,sound_velocity
 123        format(1X,2F12.8,1P6E16.8)
!CFD 10/09 Add an extra output to plot the sound speed profile easyly
            write(500,9123)i,fr,fm,sound_velocity,log_density(i),log_temperature(i)
!CFD end

!  DERIVATIVES OF DP/DR, DRHO/DR
            if(i.lt.num_shells) then
              rmid = 0.5D0*(exp(ln10*log_radius(i))+exp(ln10*log_radius(i+1)))
              dr = exp(ln10*log_radius(i+1)) - exp(ln10*log_radius(i))
              divp = 0.5D0*(star%adiabatic_index_gamma1(i)*exp(ln10*log_pressure(i))+ &
                    star%adiabatic_index_gamma1(i+1)*exp(ln10*log_pressure(i+1)))*dr
              divr = 0.5D0*(exp(ln10*log_density(i))+exp(ln10*log_density(i+1)))*dr
              qpr1 = exp(ln10*(cgl+log_mass(i)+log_density(i)-2.0D0*log_radius(i)))
              qdr1 = exp(ln10*(log_density(i)-log_pressure(i)))*star%pulse_dlnrho_dlnp(i)*qpr1
              qpr2 = exp(ln10*(cgl+log_mass(i+1)+log_density(i+1)-2.0D0*log_radius(i+1)))
              qdr2 = exp(ln10*(log_density(i+1)-log_pressure(i+1)))*star%pulse_dlnrho_dlnp(i+1)*qpr1
              qqpr = (qpr1-qpr2)/divp
              qqdr = (qdr1-qdr2)/divr
              write(imodpt,124)rmid/star%solar_radius_cgs,qpr1,qdr1,qqpr,qqdr
 124          format(1X,F11.7,1P4E16.8)
            endif
         end do
      endif
!
! G Somers 11/14 LCHEMO BLOCK REMOVED, AS THIS INFO IS ALREADY IN .STORE.
!
! DBG PULSE: PRINT OUT PULSATION ENV AND ATM IN ENVINT
! G Somers 11/14 CHANGE TO NEW I/O FLAGS.
      if(star%job%lstatm .or. star%ctrl%lstenv .or. star%job%pulsation_output_active) then
!  INTEGRATE AN ENVELOPE FROM THE SURFACE TO THE CONVERGED MODEL,
!  PRINTING OUT THE RESULTS.
!  SET UP FLAGS AND COUNTERS.
       lsbc0 = .false.
       if(star%ctrl%lstore)lprt = .true.
       katm = 0
       kenv = 0
       ksaha = 0
!  SAVE THE INTEGRATION STEP PARAMETERS AND ENFORCE THE SPACING
!  REQUESTED FOR PRINTOUT PURPOSES.
       abeg0 = star%job%atm_step_begin
       amin0 = star%job%atm_step_min
       amax0 = star%job%atm_step_max
       ebeg0 = star%job%env_step_begin
       emin0 = star%job%env_step_min
       emax0 = star%job%env_step_max
       star%job%atm_step_begin = star%ctrl%atm_step_size
       star%job%atm_step_min = star%ctrl%atm_step_size
       star%job%atm_step_max = star%ctrl%atm_step_size
       star%job%env_step_begin = star%ctrl%envelope_step_size
       star%job%env_step_min = star%ctrl%envelope_step_size
       star%job%env_step_max = star%ctrl%envelope_step_size
       b = dexp(ln10*log_luminosity_lsun)
       rl = 0.5D0*(log_luminosity_lsun + star%log10_solar_luminosity - 4.0D0*log_teff - c4pil - csigl)
       gl = cgl + star%stotal - rl - rl
       x = composition(1,num_shells)
       z = composition(3,num_shells)
       fpl = shape_factor_fp(num_shells)
       ftl = shape_factor_ft(num_shells)
       ixx=0
       hstot = star%stotal
       plim = log_pressure(num_shells)
! DBG PULSE: ADDED ARGUEMENT TO ENVINT TO TURN ON/OFF PULSE OUTPUT
         lpulpt = star%job%pulsation_output_active
            if (use_debye_huckel_correction) then
               debye_huckel_x = composition(1,num_shells)
               debye_huckel_y = composition(2,num_shells)+composition(4,num_shells)
               debye_huckel_z_total = composition(3,num_shells)
               debye_huckel_z(1) = composition(5,num_shells)+composition(6,num_shells)
               debye_huckel_z(2) = composition(7,num_shells)+composition(8,num_shells)
               debye_huckel_z(3) = composition(9,num_shells)+composition(10,num_shells)+composition(11,num_shells)
            end if
! MHP 10/02  define ISTORE - used in ENVINT
         idum = 0
! G Somers 10/14, FOR SPOTTED RUNS, FIND THE
! PRESSURE AT THE AMBIENT TEMPERATURE ATEFFL
        if(envelope_cz_bottom_index.eq.num_shells.and.star%ctrl%spot_filling_factor.ne.0.0.and.star%ctrl%spot_temp_contrast.ne.1.0)then
            ateffl = log_teff - 0.25*log10(star%ctrl%spot_filling_factor * star%ctrl%spot_temp_contrast**4.0 + 1.0 - star%ctrl%spot_filling_factor)
       else
          ateffl = log_teff
       endif
       call atm_get(b,fpl,ftl,gl,hstot,ixx,lprt,lsbc0,plim,rl, &
                     ateffl,x,z,dum1,idum,katm,kenv,ksaha,dum2, &
                     dum3,dum4,lpulpt)
! G Somers END
       star%job%atm_step_begin = abeg0
       star%job%atm_step_min = amin0
       star%job%atm_step_max = amax0
       star%job%env_step_begin = ebeg0
       star%job%env_step_min = emin0
       star%job%env_step_max = emax0
      endif
!
! G Somers 11/14 LCONZO (convection zone info) block deleted.
!
! G Somers 11/14 LJOUT (rotation info) block deleted.
!
      fsi = dexp(-ln10*star%stotal)
! DBG PULSE: WRITE HEADER INFORMATION FOR PULSE MODEL
      if(star%job%pulsation_output_active) then
         rsurfl = 0.5D0*(log_luminosity_lsun - c4pil - csigl - 4.0D0*log_teff + star%log10_solar_luminosity)
         tempr = rsurfl - star%log10_solar_radius
         qsmass = star%pulsation_mass_msun
         write (star%ctrl%opal_model_unit, 5001) model_number,num_pulsation_points,star%ctrl%pulsation_file_version,qsmass, &
               log_teff,log_luminosity_lsun,tempr, age_gyr, star%mixing_length_alpha, star%job%initial_envelope_x, star%job%initial_envelope_z
 5001    format(' MODEL#=', I5, '  NUMBER OF SHELLS IN MODEL=',I5, &
                ' VER=',I2,/, &
                ' MASS=',F8.5, '  LOG(TEFF)=',F8.5,/, ' LOG(L/LSUN)=', &
                F16.10, '  LOG(R/RSUN)=',F16.10, /, &
                ' AGE=', 1PE12.5,' GYR',/, &
                ' MIXING LENGTH PARAMETER=', 0PF16.10,/, &
                ' ZAMS (X,Z)=', 2F18.10)
      end if
      do j = 1,num_pulsation_points
       i = pulsation_point_id(j)
       fs = fsi*mass_coordinate(i)
!
! G Somers 11/14 FINAL LSCRIB BLOCK DELETED
!
! DBG WRITE PULSE MODEL
!       PRINT*, 'LPULSE=',LPULSE
         if (star%job%pulsation_output_active.and.star%ctrl%lstore) then
! MHP 10/02 uncommented pelpf statement, used later in i/o
!         PELPF = CGAS * DEXP(CLN*(HT(I) + HD(I)))* PEMU(I)
!          ADDED X AND Z TO OUTPUT
         if (j.ne.2 .or. i.ne.1) then
         if(star%ctrl%pulsation_file_version.eq.1) then
         pelpf = gas_constant * dexp(ln10*(log_temperature(i) + log_density(i)))* star%pulse_electron_mean_molecular_weight(i)
         write(star%ctrl%opal_model_unit, 5052)log_radius(i),fs,log_luminosity(i),log_temperature(i),log_density(i), &
                log_pressure(i), star%sesum(i),star%so(i), star%pulse_dlnrho_dlnp(i), star%pulse_dlneps_dlnrho(i), &
                star%pulse_dlneps_dlnt(i), star%pulse_dlnkap_dlnrho(i), star%pulse_dlnkap_dlnt(i), star%del_grad(i_grad_actual,i),star%del_grad(i_grad_ad,i), &
                star%pulse_specific_heat(i), star%pulse_mean_molecular_weight(i), star%pulse_dlnrho_dlnt(i), pelpf
         else if (star%ctrl%pulsation_file_version.eq.2) then
         write(star%ctrl%opal_model_unit, 6052)log_radius(i),fs,log_luminosity(i),log_temperature(i),log_density(i), &
            log_pressure(i), star%sesum(i),star%so(i), star%pulse_dlnrho_dlnp(i), star%pulse_dlneps_dlnrho(i), &
            star%pulse_dlneps_dlnt(i), star%pulse_dlnkap_dlnrho(i), star%pulse_dlnkap_dlnt(i), star%del_grad(i_grad_actual,i),star%del_grad(i_grad_ad,i), &
            star%pulse_specific_heat(i), star%pulse_mean_molecular_weight(i), star%pulse_dlnrho_dlnt(i), composition(1,i),composition(3,i)
         else if (star%ctrl%pulsation_file_version.eq.3) then
! DBG 7/95 Modifed to include mixing length variables
         write(star%ctrl%opal_model_unit, 6053)log_radius(i),fs,log_luminosity(i),log_temperature(i),log_density(i),star%rot%valfmlt(i), &
            log_pressure(i), star%sesum(i),star%so(i), star%pulse_dlnrho_dlnp(i), star%pulse_dlneps_dlnrho(i),star%rot%vphmlt(i), &
            star%pulse_dlneps_dlnt(i), star%pulse_dlnkap_dlnrho(i), star%pulse_dlnkap_dlnt(i), star%del_grad(i_grad_actual,i),star%del_grad(i_grad_ad,i),star%rot%vcmxmlt(i), &
            star%pulse_specific_heat(i), star%pulse_mean_molecular_weight(i), star%pulse_dlnrho_dlnt(i), composition(1,i),composition(3,i)
         end if
         end if
 5052      format(5E16.9,/,5E16.9,/,5E16.9,/,5E16.9)
 6052      format(5E23.16,/,5E23.16,/,5E23.16,/,5E23.16)
 6053      format(6E23.16,/,6E23.16,/,6E23.16,/,5E23.16)
      end if
! DBG END
      end do
!
! G Somers 11/14 REMOVED LONG BLOCK
!
      return
end subroutine write_pulsation_model
