!----------------------------------------------------------------------
! rscale
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original rscale.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! DBG 5/94 Added rescaling of Z in core ZRAMP stuff.
!
! Rescales a model's composition and/or mass at the start of a new
! kind-card run (run_index), per the targets in common/ckind/'s
! rescale_params(:,run_index): surface X and interior Z are rescaled
! multiplicatively/additively across the whole model for main-sequence
! stars (composition(1,1) > 1e-10), while horizontal-branch models
! (helium-core/H-shell stars) instead rescale the core mass and/or
! total mass, holding the H-burning shell mass fixed and redistributing
! the difference into the envelope. Also optionally rescales one extra
! species (common/newcmp/) and ramps interior Z toward a target core
! value (common/zramp/).
subroutine rscale(luminosity_array, composition, shell_mass_log, &
     total_mass_log, num_zones, run_index, star_mass, convective_flag, ierr)
! DBG 5/94 Added rescaling of Z in core ZRAMP stuff.
      use const_lib
      use star_info_lib, only: star
      use luout_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: luminosity_array(json)
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(inout) :: shell_mass_log(json)
      double precision, intent(inout) :: total_mass_log
      integer, intent(in) :: num_zones, run_index
      double precision, intent(inout) :: star_mass
      logical, intent(in) :: convective_flag(json)








      save

! --- locals ---
      double precision :: x_rescale_factor
      integer :: zone_idx, species_idx, error_species_index, icomp
      double precision :: delta_z, z_rescale_factor
      double precision :: log_mass_shift
      double precision :: core_mass_old, env_mass_old, env_mass_old2, &
           shell_mass_old, delta_env_mass, delta_core_mass, env_mass_new, &
           mass_scale_factor, env_mass_check, shell_mass_prev, &
           mass_shift_grams, unlogged_mass_temp, core_mass_new, &
           shell_mass_new, env_mass_total_check
      double precision :: z_ramp_slope, mass_fraction_local, z_ramp_value
      integer :: core_edge, envelope_edge, shell_begin, shell_end, shell_mid
      logical :: has_h_shell

! ************
!      write(*,*)
!      write(*,*)'Entering rscale ',HSTOT-HS(M)
! ************

      integer, intent(out) :: ierr

      ierr = 0

      if(rescale_params(2,run_index).gt.0.0d0) then
!  RESCALE X BY MULTIPLYING ALL SHELL X VALUES BY THE RATIO (XNEW/XOLD)
!  WHERE XOLD = OLD SURFACE X VALUE
!  THIS METHOD OK FOR BOTH HORIZONTAL BRANCH AND MAIN SEQUENCE STARS
         if(rescale_params(2,run_index).le.1.0d0) then
            x_rescale_factor = rescale_params(2,run_index)/dmax1(star%env_comp%xnew,1.0d-20)
            do 10 zone_idx = 1,num_zones
               composition(1,zone_idx) = dmin1(rescale_params(2,run_index), &
                    composition(1,zone_idx)*x_rescale_factor)
   10       continue
            star%env_comp%xnew = rescale_params(2,run_index)
            initial_envelope_x = rescale_params(2,run_index)
! DBG 4/95 BUG FIX XENV IS USED IN SOME ROUTINES AND NOT XENV0 SO CHANGE
!     XENV WHENEVER X IS CHANGED.
! MHP 7/99 THIS IS NOT A BUG, IT IS NECESSARY.
! XENV IS USED TO SET THE INITIAL MU FOR THE YALE EQUATION OF STATE,
! AND IT SHOULD NOT BE CHANGED AFTER THE FIRST RUN.
!            XENV = XENV0
         else
!  DESIRED X >100%; X NOT CHANGED
!            ACOMP = ' X '
            error_species_index = 2
            write(short_file_unit,1000)run_index,rescale_params(error_species_index,run_index)
 1000       format(1x,'ERROR IN SUBROUTINE RSCALE'/1x,'RESCALING OF X', &
           ' IN KIND CARD #',i3,1x,'FAILED - DESIRED COMP',f9.6/1x, &
           'GREATER THAN UNITY.  X NOT RESCALED')
         endif
      endif
      if(rescale_params(3,run_index).ge.0.0d0) then
!  RESCALE Z BY ADDING (RESCAL(3,NK)-OLD Z) TO EACH Z VALUE IN THE STAR
!  THE CNO CYCLE ELEMENTS AND HE3 ARE MULTIPLIED BY THE RATIO OF THE
!  DESIRED NEW Z TO THE OLD Z - LIGHT ELEMENT ABUNDANCES ARE LEFT ALONE
         if(rescale_params(3,run_index).le.1.0d0) then
            delta_z = rescale_params(3,run_index) - star%env_comp%znew
            do 30 species_idx = 1,num_zones
               z_rescale_factor=dmax1(0.d0,composition(3,species_idx)+delta_z)/ &
                    (composition(3,species_idx)+1.d-30)
               composition(3,species_idx) = dmax1(0.0d0,composition(3,species_idx)+delta_z)
               do 20 zone_idx = 5,11
                  composition(zone_idx,species_idx) = composition(zone_idx,species_idx)*z_rescale_factor
   20          continue
   30       continue
            star%env_comp%znew = rescale_params(3,run_index)
            initial_envelope_z = rescale_params(3,run_index)
! DBG 4/95 BUG FIX ZENV IS USED IN MANY ROUTINES AND NOT ZENV0 SO CHANGE
!     ZENV WHENEVER Z IS CHANGED.
            star%env_comp%envelope_metal_fraction = initial_envelope_z
         else
!  DESIRED Z >100%; Z NOT CHANGED
!            ACOMP = ' Z '
            error_species_index = 3
            write(short_file_unit,1002)run_index,rescale_params(error_species_index,run_index)
 1002       format(1x,'ERROR IN SUBROUTINE RSCALE'/1x,'RESCALING OF Z', &
           ' IN KIND CARD #',i3,1x,'FAILED - DESIRED COMP',f9.6/1x, &
           'GREATER THAN UNITY.  Z NOT RESCALED')
         endif
      endif
      if(composition(1,1).gt.1.0d-10) then
! ************
!      write(*,*)'Entering MS change ',HSTOT-HS(M)
! ************

!  MAIN SEQUENCE RESCALING - MASS AND SINGLE ELEMENT.
         if(rescale_species_active.and.new_species_value.ge.0.0d0) then
!  RESCALE THE ABUNDANCE OF ONE ELEMENT OTHER THAN X,Y,Z
!  JNEWCP = INDEX OF ELEMENT TO BE CHANGED IN MATRIX HCOMP
!  IF LREL = T, ABUNDANCE IS RELATIVE TO SURFACE HYDROGEN ABUNDANCE
!  ON A LOGARITHMIC SCALE WHERE X ABUNDANCE = 12.0
!  E.G. AN ABUNDANCE OF 3.0 MEANS 1.0D-9* SURFACE H ABUNDANCE
            if(value_relative_to_h) new_species_value = &
                 dexp(ln10*(new_species_value-12.0d0))*composition(1,num_zones)
            if(new_species_value.lt.1.0d0) then
               do 35 zone_idx = 1,num_zones
                  composition(new_species_index,zone_idx) = new_species_value
   35          continue
            else
!  ERROR - RESCALED ABUNDANCE >100% - ABUNDANCE NOT CHANGED
!               ACOMP = ANEWCP
               write(short_file_unit,1004)icomp,run_index,new_species_value
 1004       format(1x,'ERROR IN SUBROUTINE RSCALE'/1x,'RESCALING OF ', &
           'SPECIES ',i3,' IN KIND CARD # ',i3,' FAILED - DESIRED COMP',f9.6/1x, &
           'GREATER THAN UNITY.  Z NOT RESCALED')
            endif
         endif
!  RESCALE STAR MASS BY MULTIPLYING ALL MASS POINTS IN THE STAR BY
!  A SCALE FACTOR(MNEW/MOLD) WHILE LEAVING ALL OTHER LOCAL VARIABLES
!  UNCHANGED.
         if(rescale_params(1,run_index).gt.0.0d0) then
            log_mass_shift = dlog10(rescale_params(1,run_index)/star_mass)
            total_mass_log = total_mass_log + log_mass_shift
            star_mass = rescale_params(1,run_index)
            star%env_comp%stotal = total_mass_log
            do 40 zone_idx = 1,num_zones
               shell_mass_log(zone_idx) = shell_mass_log(zone_idx) + log_mass_shift
   40       continue
         endif
      else
!  HORIZONTAL BRANCH RESCALING;CORE MASS AND TOTAL MASS.
!  FOR TOTAL MASS,TREAT AS MAIN SEQUENCE BUT ONLY CHANGE
!  THE MASS POINTS OUTSIDE THE H-BURNING SHELL;
!  FOR CORE MASS CHANGE THE POINTS INSIDE THE H-BURNING
!  SHELL AND CORRESPONDINGLY REDUCE/INCREASE THE MASS OUTSIDE
!  THE SHELL.
! ************
!      write(*,*)'Entering HB scale ',HSTOT-HS(M)
! ************


         if(rescale_params(1,run_index).gt.0.0d0.or.rescale_params(4,run_index).gt.0.0d0) then
!  RESCALE MASS AND/OR CORE MASS
!  FIND H-BURNING SHELL.
!

! ************* Added call to FINDSHELL (gn - 10/05)**************

            call findsh(composition,luminosity_array,convective_flag,num_zones, &
                 core_edge,envelope_edge,shell_begin,shell_end,shell_mid, &
                 has_h_shell)
! ****************************************************************

! FIND CORE, ENVELOPE, AND H-BURNING SHELL MASSES.
            core_mass_old = exp(ln10*shell_mass_log(shell_begin-1))/solar_mass_cgs
            env_mass_old = (exp(ln10*total_mass_log)-exp(ln10*shell_mass_log(shell_end)))/solar_mass_cgs
            env_mass_old2 = (exp(ln10*shell_mass_log(num_zones))-exp(ln10*shell_mass_log(shell_end)))/solar_mass_cgs
            shell_mass_old = (exp(ln10*shell_mass_log(shell_end))-exp(ln10*shell_mass_log(shell_begin-1)))/solar_mass_cgs
            delta_env_mass = 0.0d0
            delta_core_mass = 0.0d0
            write(short_file_unit,63)core_mass_old,shell_mass_old,env_mass_old,star_mass
   63       format(1x,'HB-OLD MASSES: CORE ', &
           f9.6,' SHELL ',f9.6,' ENV ',f9.6,' TOTAL',f9.6)
         endif

         if(rescale_params(1,run_index).gt.0.0d0) then

! MASS RESCALING : ADD OR SUBTRACT MASS FROM THE ENVELOPE OUTSIDE
! THE H-BURNING SHELL (ONLY).
!
! ENSURE THAT THE NEW ENVELOPE MASS IS POSITIVE.


            delta_env_mass = rescale_params(1,run_index) - star_mass
            env_mass_new = env_mass_old + delta_env_mass

! ** Reduce either total mass outside core OR the standard envelope only ***
            if((((10d0**total_mass_log)/solar_mass_cgs)-core_mass_old).ge.0.01d0) then
! ************
      write(*,*)'Entering old method ',(10**total_mass_log-10**shell_mass_log(num_zones))/ &
                                        (10**shell_mass_log(num_zones)-10**shell_mass_log(shell_begin-1))
      write(*,*)'SENV ',10**total_mass_log - 10**shell_mass_log(num_zones)
      write(*,*)'Envelope ',10**shell_mass_log(num_zones) - 10**shell_mass_log(shell_begin-1)
! ************

! **************************************************************************
!               write(*,*)JXBEG-1,JXMID,JXEND,M
            if(env_mass_new.le.0.0d0)then
               write(short_file_unit,69)env_mass_old,env_mass_old+delta_env_mass,rescale_params(1,run_index),star_mass
   69          format(1x,'ERROR IN SUBROUTINE RSCALE'/1x, &
           'DESIRED NEW ENVELOPE MASS LESS THAN ZERO'/1x, &
           'OLD ENVELOPE MASS ',1pe9.2,' NEW ENVELOPE ',e9.2, &
           ' NEW AND OLD TOTAL MASS ',2e10.2)
               ! 2026 (phase five, step B): stop converted to ierr; run_yrec
               ! returns the error and the CLI wrapper (main) stops.
               ierr = 1
               return
            endif

! ***** Calculate scale factor for mass rescaling *****

!            HSTOT1 = DLOG10(RESCAL(1,NK)/SMASS)
        mass_scale_factor =(rescale_params(1,run_index)-exp(ln10*shell_mass_log(shell_end))/solar_mass_cgs)/ &
                          (star_mass-exp(ln10*shell_mass_log(shell_end))/solar_mass_cgs)

! *****************************************************

!            HSTOT = HSTOT + HSTOT1
        total_mass_log=dlog10(10**shell_mass_log(shell_end)+mass_scale_factor*(10**total_mass_log-10**shell_mass_log(shell_end)))

            star_mass = rescale_params(1,run_index)
            star%env_comp%stotal = total_mass_log
            do 70 zone_idx = shell_end+1,num_zones
        shell_mass_log(zone_idx)=dlog10(10**shell_mass_log(shell_end)+ &
             mass_scale_factor*(10**shell_mass_log(zone_idx)-10**shell_mass_log(shell_end)))
   70       continue
            env_mass_old = env_mass_new
! *****
!            write(*,*)'Mass difference ',DENV
!            write(*,*)'HSTOT1 ',HSTOT1
!            write(*,*)'SFACTOR ',DLOG10(SFACTOR)
! ************
!      write(*,*)'leaving old method ',HSTOT-HS(M)
! ************

         else

! ************
      write(*,*)'Entering new method ',(10**total_mass_log-10**shell_mass_log(num_zones))/ &
                                        (10**shell_mass_log(num_zones)-10**shell_mass_log(shell_begin-1))
      write(*,*)'SENV ',10**total_mass_log - 10**shell_mass_log(num_zones)
      write(*,*)'Envelope ',10**shell_mass_log(num_zones) - 10**shell_mass_log(shell_begin-1)
! ************
!           *** print debug info ***
!            write(*,*)JXBEG-1,JXMID,JXEND,M
!            write(*,*)(10**HS(JXBEG-1))/CMSUN,' core'
!            write(*,*)(10**HS(JXMID))/CMSUN,' mid'
!            write(*,*)(10**HS(JXEND))/CMSUN,' end'
!            write(*,*)(10**HS(M))/CMSUN,' M'
!            write(*,*)(10**HSTOT)/CMSUN,' total'
            env_mass_check = ((10**total_mass_log)/solar_mass_cgs)-core_mass_old
            if(env_mass_check.le.0.0d0)then
               write(short_file_unit,69)env_mass_old,env_mass_old+delta_env_mass,rescale_params(1,run_index),star_mass
               ! 2026 (phase five, step B): stop converted to ierr; run_yrec
               ! returns the error and the CLI wrapper (main) stops.
               ierr = 1
               return
            endif

! **** Calculate scale factor *****

!           HSTOT1 = DLOG10(RESCAL(1,NK)/SMASS)
        mass_scale_factor =(rescale_params(1,run_index)-exp(ln10*shell_mass_log(shell_begin-1))/solar_mass_cgs)/ &
                          (star_mass-exp(ln10*shell_mass_log(shell_begin-1))/solar_mass_cgs)
!        write(*,*)'hstot1',(10**HSTOT1)
!        write(*,*)'smass',SMASS

!            HSTOT = HSTOT + HSTOT1
       total_mass_log=dlog10(10**shell_mass_log(shell_begin-1)+ &
            mass_scale_factor*(10**total_mass_log-10**shell_mass_log(shell_begin-1)))

            star_mass = rescale_params(1,run_index)
            star%env_comp%stotal = total_mass_log
            do 78 zone_idx = shell_begin,num_zones
      shell_mass_log(zone_idx)=dlog10(10**shell_mass_log(shell_begin-1)+ &
           mass_scale_factor*(10**shell_mass_log(zone_idx)-10**shell_mass_log(shell_begin-1)))
 78     continue
            env_mass_old = (exp(ln10*total_mass_log)-exp(ln10*shell_mass_log(shell_end)))/solar_mass_cgs
            env_mass_total_check=(exp(ln10*shell_mass_log(num_zones))-exp(ln10*shell_mass_log(shell_begin-1)))/solar_mass_cgs
!            write(*,*)'total envelope ',ENVTOTAL
!            write(*,*)'HSTOT1 ',HSTOT1
!            write(*,*)'SFACTOR ',DLOG10(SFACTOR)

! ************
!      write(*,*)'leaving new method ',HSTOT-HS(M)
! ************

         endif
         endif
! ******************************************************************

         if(rescale_params(4,run_index).gt.0.0d0) then

! RESCALE CORE MASS.
! THE MASS OF THE H-BURNING SHELL IS HELD FIXED, AND MASS IS TRANSFERRED
! FROM THE CORE TO THE ENVELOPE (OR VICE VERSA).
!

            shell_mass_prev= (exp(ln10*shell_mass_log(shell_end)))
            delta_core_mass = rescale_params(4,run_index) - core_mass_old
            delta_env_mass = - delta_core_mass
            env_mass_new = env_mass_old2 + delta_env_mass
            if(env_mass_new.le.0.0d0)then
               write(short_file_unit,71)env_mass_old,env_mass_new,rescale_params(4,run_index),core_mass_old
   71          format(1x,'ERROR IN SUBROUTINE RSCALE'/1x, &
           'NEW ENVELOPE MASS LESS THAN ZERO BECAUSE OF CORE', &
           ' RESCALING'/1x, &
           'OLD ENVELOPE MASS ',1pe10.2,' NEW ENVELOPE ',e10.2, &
           ' NEW AND OLD CORE MASS ',2e10.2)
               ! 2026 (phase five, step B): stop converted to ierr; run_yrec
               ! returns the error and the CLI wrapper (main) stops.
               ierr = 1
               return
            endif
            log_mass_shift = dlog10(rescale_params(4,run_index)/core_mass_old)
            do 80 zone_idx = 1,shell_begin-1
               shell_mass_log(zone_idx) = shell_mass_log(zone_idx) + log_mass_shift
   80       continue

!  HOLD H-SHELL MASS CONSTANT;SHIFT IS THE CHANGE IN THE
!  UNLOGGED MASS OF EACH POINT IN THE H SHELL.
            mass_shift_grams = delta_core_mass*solar_mass_cgs
            do 90 zone_idx = shell_begin,shell_end
               unlogged_mass_temp = exp(ln10*shell_mass_log(zone_idx))
               shell_mass_log(zone_idx) = log10(unlogged_mass_temp + mass_shift_grams)
   90       continue
!  NOW SHRINK OR EXPAND THE ENVELOPE MASS TO RETAIN TOTAL CONSTANT MASS.
!

        mass_scale_factor =(exp(ln10*shell_mass_log(num_zones))-exp(ln10*shell_mass_log(shell_end)))/ &
                          (exp(ln10*shell_mass_log(num_zones))-shell_mass_prev)

            do 100 zone_idx = shell_end+1,num_zones
        shell_mass_log(zone_idx)=dlog10(10**shell_mass_log(num_zones)- &
             mass_scale_factor*(10**shell_mass_log(num_zones)-10**shell_mass_log(zone_idx)))
  100       continue

         endif

         if(rescale_params(1,run_index).gt.0.0d0.or.rescale_params(4,run_index).gt.0.0d0)then
            core_mass_new = exp(ln10*shell_mass_log(shell_begin-1))/solar_mass_cgs
            env_mass_new = (exp(ln10*total_mass_log)-exp(ln10*shell_mass_log(shell_end)))/solar_mass_cgs
            shell_mass_new = (exp(ln10*shell_mass_log(shell_end))-exp(ln10*shell_mass_log(shell_begin-1)))/solar_mass_cgs
            write(short_file_unit,101)core_mass_new,shell_mass_new,env_mass_new,star_mass
  101       format(1x,'HB-RESCALED MASSES: CORE ', &
           f9.6,' SHELL ',f9.6,' ENV ',f9.6,' TOTAL ',f9.6)

         endif
      endif
!  CHANGE Y TO REFLECT NEW X,Z, AND HE3 VALUES.
      do 110 zone_idx = 1,num_zones
         composition(2,zone_idx) = 1.0d0-composition(1,zone_idx)-composition(3,zone_idx)-composition(4,zone_idx)
  110 continue
!
! DBG 5/94 rescale interior Z if LZRAMP flag is T.
! Z is linearly adjusted from RSCLZC at the center to surface Z at
! mass fraction RSCLZM.  Compensate changing Z with X.
      if(use_z_ramp.and.(rsclzc(run_index).gt.0d0).and.(rsclzm1(run_index).gt.0d0) &
           .and.(rsclzm2(run_index).gt.0d0))then
            z_ramp_slope = (star%env_comp%znew - rsclzc(run_index))/(rsclzm2(run_index)-rsclzm1(run_index))
            do zone_idx = 1,num_zones
               mass_fraction_local = 10.0d0**shell_mass_log(zone_idx)/(star_mass*solar_mass_cgs)
               if (mass_fraction_local .lt. rsclzm1(run_index)) then
                  z_rescale_factor = (composition(3,zone_idx)-rsclzc(run_index))/composition(3,zone_idx)
                  composition(3,zone_idx) = rsclzc(run_index)
                  composition(1, zone_idx) = 1.0d0-composition(3,zone_idx)-composition(4,zone_idx)-composition(2,zone_idx)
                  do species_idx = 5,11
                      composition(species_idx,zone_idx) = composition(species_idx,zone_idx)*z_rescale_factor
                  enddo
               else if (mass_fraction_local .lt. rsclzm2(run_index)) then
                  z_ramp_value = (mass_fraction_local-rsclzm1(run_index))*z_ramp_slope+rsclzc(run_index)
                  z_rescale_factor = (composition(3,zone_idx)-z_ramp_value)/composition(3,zone_idx)
                  composition(3, zone_idx) = z_ramp_value
                  composition(1, zone_idx) = 1.0d0-composition(3,zone_idx)-composition(4,zone_idx)-composition(2,zone_idx)
                  do species_idx = 5,11
                      composition(species_idx,zone_idx) = composition(species_idx,zone_idx)*z_rescale_factor
                  enddo
               end if
            enddo
      end if
! ************
!      write(*,*)'Leaving rscale ',HSTOT-HS(M)
! ************
!      IF(HSTOT.LT.HS(M)) STOP
      return
end subroutine rscale
