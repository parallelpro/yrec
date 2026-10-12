!----------------------------------------------------------------------
! model_rescaling
!----------------------------------------------------------------------
! Rescaling a model at the start of a rescale run card (rescale_model):
! its X and Z (every model), its total mass and one species
! (main-sequence models) or its core and total mass (horizontal-branch
! models), and optionally a Z ramp toward the center.
module model_rescaling
      use shell_boundaries, only: locate_shell_boundaries
      implicit none
      private
      public :: rescale_model

contains

!----------------------------------------------------------------------
! rescale_model (legacy rscale)
!----------------------------------------------------------------------
! Rescales the model at the start of run card run_index to the targets
! in star%job%rescale_params(:,run_index): (1) the total mass in Msun,
! (2) the surface X, (3) the surface Z, (4) the core mass in Msun. A
! target is applied when positive; Z also when zero. In order:
! - X: every shell's X scaled by X_new/X_surf (every model);
! - Z: every shell's Z shifted by Z_new - Z_surf, the CNO isotopes
!   scaled with it (every model);
! - main-sequence models (central X > 1e-10): one species set to a
!   given abundance, and the mass scaled by M_new/M_old throughout;
! - horizontal-branch models (no central hydrogen, an H-burning shell
!   around a helium core): the total mass changed by stretching the
!   mass outside the H shell, and the core mass changed at fixed
!   H-shell and total mass, the envelope taking up the difference;
! - Y = 1 - X - Z - X(He3) in every shell;
! - with use_z_ramp, Z set to rsclzc in the core and ramped to the
!   surface Z (DBG 5/94).
! luminosity_array is L/Lsun (it locates the H shell); shell_mass_log
! and total_mass_log are log10 m and log10 M (g); star_mass is M in
! Msun. xnew and znew hold the old surface X and Z on entry, and the
! targets on return when X or Z is rescaled.
subroutine rescale_model(luminosity_array, composition, shell_mass_log, &
     total_mass_log, num_zones, run_index, star_mass, convective_flag, ierr)
      use physical_constants
      use star_info_lib, only: star, json
      use io_units
      use eos_lib, only: eos_set_mixture
      use math_lib
      implicit none

      double precision, intent(in) :: luminosity_array(json)
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(inout) :: shell_mass_log(json)
      double precision, intent(inout) :: total_mass_log
      integer, intent(in) :: num_zones, run_index
      double precision, intent(inout) :: star_mass
      logical, intent(in) :: convective_flag(json)
! --- locals ---
      double precision :: x_rescale_factor
      integer :: zone_idx, species_idx, error_species_index
      double precision :: delta_z, z_rescale_factor
      double precision :: log_mass_shift
      double precision :: core_mass_old, env_mass_old, env_mass_old2, &
           shell_mass_old, delta_env_mass, delta_core_mass, env_mass_new, &
           mass_scale_factor, env_mass_check, shell_mass_prev, &
           mass_shift_grams, unlogged_mass_temp, core_mass_new, &
           shell_mass_new
      double precision :: z_ramp_slope, mass_fraction_local, z_ramp_value
      integer :: core_edge, envelope_edge, shell_begin, shell_end, shell_mid
      logical :: has_h_shell

      integer, intent(out) :: ierr

      ierr = 0

      call rescale_surface_x
      call rescale_interior_z
      if(composition(1,1).gt.1.0d-10) then
! Main-sequence model: one species, and the mass.
         call rescale_one_species
         call rescale_main_sequence_mass
      else
! Horizontal-branch model: the total mass changes only the points
! outside the H-burning shell; the core mass changes the points inside
! it, and the mass outside the shell shrinks or grows to match.
         call measure_hb_masses
         call rescale_hb_total_mass
         if (ierr /= 0) return
         call rescale_hb_core_mass
         if (ierr /= 0) return
         call report_hb_masses
      endif
      call set_helium_from_x_and_z
      call ramp_interior_z
      return

contains

! ---------------------------------------------------------------
! Scale every shell's X by X_new/X_surf, X_surf = xnew being the old
! surface X, capping each at X_new:
!   X_i -> min(X_new, X_i X_new/max(X_surf, 1e-20)).
! Valid for main-sequence and horizontal-branch models alike. A
! requested X above 1 is refused with a message.
subroutine rescale_surface_x
      if(star%job%rescale_params(2,run_index).gt.0.0d0) then
         if(star%job%rescale_params(2,run_index).le.1.0d0) then
            x_rescale_factor = star%job%rescale_params(2,run_index)/dmax1(star%xnew,1.0d-20)
            do zone_idx = 1,num_zones
               composition(1,zone_idx) = dmin1(star%job%rescale_params(2,run_index), &
                    composition(1,zone_idx)*x_rescale_factor)
            end do
            star%xnew = star%job%rescale_params(2,run_index)
            star%job%initial_envelope_x = star%job%rescale_params(2,run_index)
! The envelope X, star%envelope_hydrogen_fraction (legacy XENV), is
! deliberately left alone: it sets the initial mean molecular weight
! for the Yale EOS and must not change after the first run (MHP 7/99).
         else
! Requested X above 1: X is not changed.
            error_species_index = 2
            write(run_log_unit,1000)run_index,star%job%rescale_params(error_species_index,run_index)
 1000       format(1x,'ERROR IN SUBROUTINE RSCALE'/1x,'RESCALING OF X', &
           ' IN KIND CARD #',i3,1x,'FAILED - DESIRED COMP',f9.6/1x, &
           'GREATER THAN UNITY.  X NOT RESCALED')
         endif
      endif
end subroutine rescale_surface_x

! ---------------------------------------------------------------
! Shift every shell's Z by dZ = Z_new - Z_surf, Z_surf = znew being the
! old surface Z, and scale the CNO isotopes (slots C12..O18) by each
! shell's new over old Z:
!   Z_i -> max(0, Z_i + dZ),
!   X_k,i -> X_k,i max(0, Z_i + dZ)/(Z_i + 1e-30),  k = C12..O18.
! He3 and the light elements are left alone. The target is applied when
! >= 0, so a target Z of 0 sets Z = 0. A requested Z above 1 is refused
! with a message.
subroutine rescale_interior_z
      if(star%job%rescale_params(3,run_index).ge.0.0d0) then
         if(star%job%rescale_params(3,run_index).le.1.0d0) then
            delta_z = star%job%rescale_params(3,run_index) - star%znew
            do zone_idx = 1,num_zones
               z_rescale_factor=dmax1(0.d0,composition(3,zone_idx)+delta_z)/ &
                    (composition(3,zone_idx)+1.d-30)
               composition(3,zone_idx) = dmax1(0.0d0,composition(3,zone_idx)+delta_z)
               do species_idx = 5,11
                  composition(species_idx,zone_idx) = composition(species_idx,zone_idx)*z_rescale_factor
               end do
            end do
            star%znew = star%job%rescale_params(3,run_index)
            star%job%initial_envelope_z = star%job%rescale_params(3,run_index)
! Unlike X, the envelope Z (star%envelope_metal_fraction, legacy ZENV)
! follows the new Z, since many routines read it (DBG 4/95), and the
! EOS mixture with it.
            star%envelope_metal_fraction = star%job%initial_envelope_z
            call eos_set_mixture(star%envelope_hydrogen_fraction, &
                 star%envelope_metal_fraction, star%amuenv, &
                 star%fxenv)
         else
! Requested Z above 1: Z is not changed.
            error_species_index = 3
            write(run_log_unit,1002)run_index,star%job%rescale_params(error_species_index,run_index)
 1002       format(1x,'ERROR IN SUBROUTINE RSCALE'/1x,'RESCALING OF Z', &
           ' IN KIND CARD #',i3,1x,'FAILED - DESIRED COMP',f9.6/1x, &
           'GREATER THAN UNITY.  Z NOT RESCALED')
         endif
      endif
end subroutine rescale_interior_z

! ---------------------------------------------------------------
! With rescale_species_active (legacy LNEWCP), set one species other
! than X, Y and Z, slot new_species_index, to new_species_value in
! every shell (main-sequence models).
subroutine rescale_one_species
         if(star%ctrl%rescale_species_active.and.star%job%new_species_value.ge.0.0d0) then
! new_species_value is a mass fraction or, with value_relative_to_h
! (legacy LREL), a log abundance relative to the surface X on a scale
! where X is 12, by mass:
!   X_k = 10^(value - 12) X_surf,
! so 3.0 means 1e-9 times the surface X. The converted value
! overwrites star%job%new_species_value.
            if(star%ctrl%value_relative_to_h) star%job%new_species_value = &
                 exp(ln10*(star%job%new_species_value-12.0d0))*composition(1,num_zones)
            if(star%job%new_species_value.lt.1.0d0) then
               do zone_idx = 1,num_zones
                  composition(star%ctrl%new_species_index,zone_idx) = star%job%new_species_value
               end do
            else
! An abundance of 1 or more is refused; the species is not changed.
               write(run_log_unit,1004)star%ctrl%new_species_index,run_index,star%job%new_species_value
 1004       format(1x,'ERROR IN SUBROUTINE RSCALE'/1x,'RESCALING OF ', &
           'SPECIES ',i3,' IN KIND CARD # ',i3,' FAILED - DESIRED COMP',f9.6/1x, &
           'GREATER THAN UNITY.  Z NOT RESCALED')
            endif
         endif
end subroutine rescale_one_species

! ---------------------------------------------------------------
! Scale the total mass and every mass point by M_new/M_old
! (main-sequence models): log M and every log m shift by
! log10(M_new/M_old); nothing else in the model changes.
subroutine rescale_main_sequence_mass
         if(star%job%rescale_params(1,run_index).gt.0.0d0) then
            log_mass_shift = log10(star%job%rescale_params(1,run_index)/star_mass)
            total_mass_log = total_mass_log + log_mass_shift
            star_mass = star%job%rescale_params(1,run_index)
            star%stotal = total_mass_log
            do zone_idx = 1,num_zones
               shell_mass_log(zone_idx) = shell_mass_log(zone_idx) + log_mass_shift
            end do
         endif
end subroutine rescale_main_sequence_mass

! ---------------------------------------------------------------
! When a horizontal-branch model's mass or core mass is to be rescaled,
! find its H-burning shell, points shell_begin..shell_end
! (locate_shell_boundaries, gn 10/05), and its old masses in Msun:
!   core           m(shell_begin-1),
!   shell          m(shell_end) - m(shell_begin-1),
!   envelope       M - m(shell_end),
!   env_mass_old2  m(nz) - m(shell_end), the envelope inside the
!                  fitting point.
subroutine measure_hb_masses
         if(star%job%rescale_params(1,run_index).gt.0.0d0.or.star%job%rescale_params(4,run_index).gt.0.0d0) then
            call locate_shell_boundaries(composition,luminosity_array,convective_flag,num_zones, &
                 core_edge,envelope_edge,shell_begin,shell_end,shell_mid, &
                 has_h_shell)
            core_mass_old = exp(ln10*shell_mass_log(shell_begin-1))/star%solar_mass_cgs
            env_mass_old = (exp(ln10*total_mass_log)-exp(ln10*shell_mass_log(shell_end)))/star%solar_mass_cgs
            env_mass_old2 = (exp(ln10*shell_mass_log(num_zones))-exp(ln10*shell_mass_log(shell_end)))/star%solar_mass_cgs
            shell_mass_old = (exp(ln10*shell_mass_log(shell_end))-exp(ln10*shell_mass_log(shell_begin-1)))/star%solar_mass_cgs
            delta_env_mass = 0.0d0
            delta_core_mass = 0.0d0
            write(run_log_unit,63)core_mass_old,shell_mass_old,env_mass_old,star_mass
   63       format(1x,'HB-OLD MASSES: CORE ', &
           f9.6,' SHELL ',f9.6,' ENV ',f9.6,' TOTAL',f9.6)
         endif
end subroutine measure_hb_masses

! ---------------------------------------------------------------
! Rescale a horizontal-branch model's total mass to M_new by stretching
! the mass above a fixed base m_b:
!   m_i -> m_b + s (m_i - m_b),  s = (M_new - m_b)/(M_old - m_b),
! and likewise M. The base is the top of the H shell, m(shell_end),
! when at least 0.01 Msun lies outside the core ("old method"), else
! the top of the core, m(shell_begin-1) ("new method"). Each method
! prints diagnostics to the terminal. Sets ierr if the envelope mass
! would not be positive.
subroutine rescale_hb_total_mass
         if(star%job%rescale_params(1,run_index).gt.0.0d0) then

            delta_env_mass = star%job%rescale_params(1,run_index) - star_mass
            env_mass_new = env_mass_old + delta_env_mass

            if((((pow(10d0, total_mass_log))/star%solar_mass_cgs)-core_mass_old).ge.0.01d0) then
      write(*,*)'Entering old method ',(exp10(total_mass_log)-exp10(shell_mass_log(num_zones)))/ &
                                        (exp10(shell_mass_log(num_zones))-exp10(shell_mass_log(shell_begin-1)))
      write(*,*)'SENV ',exp10(total_mass_log) - exp10(shell_mass_log(num_zones))
      write(*,*)'Envelope ',exp10(shell_mass_log(num_zones)) - exp10(shell_mass_log(shell_begin-1))
! The new envelope mass outside the H shell must be positive.
            if(env_mass_new.le.0.0d0)then
               write(run_log_unit,69)env_mass_old,env_mass_old+delta_env_mass,star%job%rescale_params(1,run_index),star_mass
   69          format(1x,'ERROR IN SUBROUTINE RSCALE'/1x, &
           'DESIRED NEW ENVELOPE MASS LESS THAN ZERO'/1x, &
           'OLD ENVELOPE MASS ',1pe9.2,' NEW ENVELOPE ',e9.2, &
           ' NEW AND OLD TOTAL MASS ',2e10.2)
               ! run_yrec returns the error; the command-line wrapper
               ! (main) stops.
               ierr = 1
               return
            endif

! Base m_b = m(shell_end); s from the masses in Msun.
        mass_scale_factor =(star%job%rescale_params(1,run_index)-exp(ln10*shell_mass_log(shell_end))/star%solar_mass_cgs)/ &
                          (star_mass-exp(ln10*shell_mass_log(shell_end))/star%solar_mass_cgs)

        total_mass_log=log10(exp10(shell_mass_log(shell_end))+mass_scale_factor*(exp10(total_mass_log)-exp10(shell_mass_log(shell_end))))

            star_mass = star%job%rescale_params(1,run_index)
            star%stotal = total_mass_log
            do zone_idx = shell_end+1,num_zones
        shell_mass_log(zone_idx)=log10(exp10(shell_mass_log(shell_end))+ &
             mass_scale_factor*(exp10(shell_mass_log(zone_idx))-exp10(shell_mass_log(shell_end))))
            end do
            env_mass_old = env_mass_new

         else

      write(*,*)'Entering new method ',(exp10(total_mass_log)-exp10(shell_mass_log(num_zones)))/ &
                                        (exp10(shell_mass_log(num_zones))-exp10(shell_mass_log(shell_begin-1)))
      write(*,*)'SENV ',exp10(total_mass_log) - exp10(shell_mass_log(num_zones))
      write(*,*)'Envelope ',exp10(shell_mass_log(num_zones)) - exp10(shell_mass_log(shell_begin-1))
! The mass outside the core must be positive. (This tests the old
! mass, M_old - M_core, not the new one.)
            env_mass_check = ((exp10(total_mass_log))/star%solar_mass_cgs)-core_mass_old
            if(env_mass_check.le.0.0d0)then
               write(run_log_unit,69)env_mass_old,env_mass_old+delta_env_mass,star%job%rescale_params(1,run_index),star_mass
               ! run_yrec returns the error; the command-line wrapper
               ! (main) stops.
               ierr = 1
               return
            endif

! Base m_b = m(shell_begin-1); s from the masses in Msun.
        mass_scale_factor =(star%job%rescale_params(1,run_index)-exp(ln10*shell_mass_log(shell_begin-1))/star%solar_mass_cgs)/ &
                          (star_mass-exp(ln10*shell_mass_log(shell_begin-1))/star%solar_mass_cgs)
       total_mass_log=log10(exp10(shell_mass_log(shell_begin-1))+ &
            mass_scale_factor*(exp10(total_mass_log)-exp10(shell_mass_log(shell_begin-1))))

            star_mass = star%job%rescale_params(1,run_index)
            star%stotal = total_mass_log
            do zone_idx = shell_begin,num_zones
      shell_mass_log(zone_idx)=log10(exp10(shell_mass_log(shell_begin-1))+ &
           mass_scale_factor*(exp10(shell_mass_log(zone_idx))-exp10(shell_mass_log(shell_begin-1))))
            end do
            env_mass_old = (exp(ln10*total_mass_log)-exp(ln10*shell_mass_log(shell_end)))/star%solar_mass_cgs
         endif
         endif
end subroutine rescale_hb_total_mass

! ---------------------------------------------------------------
! Rescale a horizontal-branch model's core mass to M_core,new at fixed
! H-shell mass and total mass, moving mass between the core and the
! envelope:
! - the core's log m shift by log10(M_core,new/M_core,old);
! - the H-shell points shift by dM_core = M_core,new - M_core,old in m;
! - the envelope points above the shell move toward m(nz), which stays
!   fixed (as does M):
!     m_i -> m_nz - s (m_nz - m_i),
!     s = (m_nz - m_new(shell_end))/(m_nz - m_old(shell_end)).
! Sets ierr if the envelope inside the fitting point, env_mass_old2,
! would not keep a positive mass.
subroutine rescale_hb_core_mass
         if(star%job%rescale_params(4,run_index).gt.0.0d0) then

            shell_mass_prev= (exp(ln10*shell_mass_log(shell_end)))
            delta_core_mass = star%job%rescale_params(4,run_index) - core_mass_old
            delta_env_mass = - delta_core_mass
            env_mass_new = env_mass_old2 + delta_env_mass
            if(env_mass_new.le.0.0d0)then
               write(run_log_unit,71)env_mass_old,env_mass_new,star%job%rescale_params(4,run_index),core_mass_old
   71          format(1x,'ERROR IN SUBROUTINE RSCALE'/1x, &
           'NEW ENVELOPE MASS LESS THAN ZERO BECAUSE OF CORE', &
           ' RESCALING'/1x, &
           'OLD ENVELOPE MASS ',1pe10.2,' NEW ENVELOPE ',e10.2, &
           ' NEW AND OLD CORE MASS ',2e10.2)
               ! run_yrec returns the error; the command-line wrapper
               ! (main) stops.
               ierr = 1
               return
            endif
            log_mass_shift = log10(star%job%rescale_params(4,run_index)/core_mass_old)
            do zone_idx = 1,shell_begin-1
               shell_mass_log(zone_idx) = shell_mass_log(zone_idx) + log_mass_shift
            end do

! Every H-shell point moves by dM_core (g), which keeps the shell's
! mass.
            mass_shift_grams = delta_core_mass*star%solar_mass_cgs
            do zone_idx = shell_begin,shell_end
               unlogged_mass_temp = exp(ln10*shell_mass_log(zone_idx))
               shell_mass_log(zone_idx) = log10(unlogged_mass_temp + mass_shift_grams)
            end do
! Shrink or stretch the envelope toward m(nz) to keep the total mass.
        mass_scale_factor =(exp(ln10*shell_mass_log(num_zones))-exp(ln10*shell_mass_log(shell_end)))/ &
                          (exp(ln10*shell_mass_log(num_zones))-shell_mass_prev)

            do zone_idx = shell_end+1,num_zones
        shell_mass_log(zone_idx)=log10(exp10(shell_mass_log(num_zones))- &
             mass_scale_factor*(exp10(shell_mass_log(num_zones))-exp10(shell_mass_log(zone_idx))))
            end do

         endif
end subroutine rescale_hb_core_mass

! ---------------------------------------------------------------
! Echo the rescaled horizontal-branch core, shell, envelope and total
! masses.
subroutine report_hb_masses
         if(star%job%rescale_params(1,run_index).gt.0.0d0.or.star%job%rescale_params(4,run_index).gt.0.0d0)then
            core_mass_new = exp(ln10*shell_mass_log(shell_begin-1))/star%solar_mass_cgs
            env_mass_new = (exp(ln10*total_mass_log)-exp(ln10*shell_mass_log(shell_end)))/star%solar_mass_cgs
            shell_mass_new = (exp(ln10*shell_mass_log(shell_end))-exp(ln10*shell_mass_log(shell_begin-1)))/star%solar_mass_cgs
            write(run_log_unit,101)core_mass_new,shell_mass_new,env_mass_new,star_mass
  101       format(1x,'HB-RESCALED MASSES: CORE ', &
           f9.6,' SHELL ',f9.6,' ENV ',f9.6,' TOTAL ',f9.6)

         endif
end subroutine report_hb_masses

! ---------------------------------------------------------------
! Set every shell's Y from its new X, Z and He3: Y = 1 - X - Z - X(He3).
subroutine set_helium_from_x_and_z
      do zone_idx = 1,num_zones
         composition(2,zone_idx) = 1.0d0-composition(1,zone_idx)-composition(3,zone_idx)-composition(4,zone_idx)
      end do
end subroutine set_helium_from_x_and_z

! ---------------------------------------------------------------
! With use_z_ramp (legacy LZRAMP; DBG 5/94) and rsclzc, rsclzm1 and
! rsclzm2 all positive, rescale the interior Z as a function of the
! mass fraction q = m/M: Z = rsclzc for q < rsclzm1, then linear up to
! the surface Z (znew) at q = rsclzm2,
!   Z(q) = rsclzc + (q - rsclzm1)(znew - rsclzc)/(rsclzm2 - rsclzm1);
! Z beyond rsclzm2 is left alone. X takes up the change, Y being kept:
! X = 1 - Z - X(He3) - Y. The CNO isotopes (C12..O18) are multiplied by
! (Z_old - Z_new)/Z_old.
subroutine ramp_interior_z
      if(star%ctrl%use_z_ramp.and.(star%job%rsclzc(run_index).gt.0d0).and.(star%job%rsclzm1(run_index).gt.0d0) &
           .and.(star%job%rsclzm2(run_index).gt.0d0))then
            z_ramp_slope = (star%znew - star%job%rsclzc(run_index))/(star%job%rsclzm2(run_index)-star%job%rsclzm1(run_index))
            do zone_idx = 1,num_zones
               mass_fraction_local = exp10(shell_mass_log(zone_idx))/(star_mass*star%solar_mass_cgs)
               if (mass_fraction_local .lt. star%job%rsclzm1(run_index)) then
                  z_rescale_factor = (composition(3,zone_idx)-star%job%rsclzc(run_index))/composition(3,zone_idx)
                  composition(3,zone_idx) = star%job%rsclzc(run_index)
                  composition(1, zone_idx) = 1.0d0-composition(3,zone_idx)-composition(4,zone_idx)-composition(2,zone_idx)
                  do species_idx = 5,11
                      composition(species_idx,zone_idx) = composition(species_idx,zone_idx)*z_rescale_factor
                  enddo
               else if (mass_fraction_local .lt. star%job%rsclzm2(run_index)) then
                  z_ramp_value = (mass_fraction_local-star%job%rsclzm1(run_index))*z_ramp_slope+star%job%rsclzc(run_index)
                  z_rescale_factor = (composition(3,zone_idx)-z_ramp_value)/composition(3,zone_idx)
                  composition(3, zone_idx) = z_ramp_value
                  composition(1, zone_idx) = 1.0d0-composition(3,zone_idx)-composition(4,zone_idx)-composition(2,zone_idx)
                  do species_idx = 5,11
                      composition(species_idx,zone_idx) = composition(species_idx,zone_idx)*z_rescale_factor
                  enddo
               end if
            enddo
      end if
end subroutine ramp_interior_z

end subroutine rescale_model

end module model_rescaling
