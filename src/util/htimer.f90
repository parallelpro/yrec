!----------------------------------------------------------------------
! htimer
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original htimer.f; only variable names, source form, and comment
! style were updated.
!
! htimer determines the timestep for the next model.
! There are several possible means of determining the timestep; it
! chooses the most stringent one.
! First of all, the timestep may not be more than atime(13) times the
!      previous one.
! In sr ptime, the timestep based on structure changes from one
!      model to the next is computed if desired (structure_dt).
! In sr wtime, the timestep based on changes in the angular velocity
!      from one model to the next is computed (rotating models only)
!      (rotation_dt).
! In sr xtime, the timestep based on hydrogen burning is computed
!      (hydrogen_dt).
! In sr ytime, the timestep based on helium burning is computed
!      (helium_dt).
! Details of the computation of these timesteps can be found in the
!      appropriate subroutines.
! The final timestep (previous_timestep) is the minimum of the above
!      timesteps.
! **note: if a criterion is inapplicable or turned off by user
!      parameters then the 'timestep' for that criterion is
!      arbitrarily set to a large number.
! If a fixed timestep is being enforced (timestep_override_active=T),
! the above criteria are short-circuited and a timestep
! (timestep_override) is enforced. Caveat user.
! If the model is being evolved to a fixed age (end_age_stop_active=T),
! the timestep will be reduced if the one used would put the model
! past the desired age (target_end_age).
!
! A negative timestep indicates that a model is being relaxed
! and not evolved. Check for negative dt and exit if found.
subroutine htimer(previous_timestep, hydrogen_dt, num_points, log_density, &
     luminosity, enclosed_mass, shell_mass, log_temperature, composition, &
     convective_core_edge_zone, h_shell_midpoint_zone, &
     luminosity_components, age_gyr, timestep_years, kind_card_index, &
     log_pressure, log_radius, omega, max_domega_frac, h_shell_zone_begin, &
     log_teff)
      use star_info_lib, only: star
      use star_info_lib, only: star, json
      use phys_const_lib
      implicit none

      double precision, intent(inout) :: previous_timestep
      double precision, intent(out) :: hydrogen_dt
      integer, intent(in) :: num_points
      double precision, intent(in) :: log_density(json), luminosity(json), &
           enclosed_mass(json), shell_mass(json), log_temperature(json)
      double precision, intent(in) :: composition(15,json)
      integer, intent(in) :: convective_core_edge_zone, &
           h_shell_midpoint_zone
      double precision, intent(in) :: luminosity_components(8)
      double precision, intent(in) :: age_gyr
      double precision, intent(out) :: timestep_years
      integer, intent(in) :: kind_card_index
      double precision, intent(in) :: log_pressure(json), log_radius(json), &
           omega(json)
      double precision, intent(out) :: max_domega_frac
      integer, intent(in) :: h_shell_zone_begin
      double precision, intent(in) :: log_teff
      double precision :: structure_dt, rotation_dt, helium_dt, &
           hydrogen_luminosity, envelope_dt, time_left_years
      double precision :: energy_gen_terms(6)
      double precision :: rate_pp(json), rate_he3_he3(json), &
           rate_he3_he4(json), rate_c12_p(json), rate_c13_p(json), &
           rate_n14_p(json), rate_o16_p(json), rate_c13_alpha(json), &
           rate_zero9(json), rate_c12_alpha(json), rate_n14_alpha(json), &
           rate_triple_alpha(json), rate_zero13(json)
      double precision :: frac_c12_alpha(json), frac_be7_electron(json)

      if (previous_timestep.ge.0.0d0) then
! if user is fixing tstep, set dt to given value and exit
      if(star%job%timestep_override_active(kind_card_index)) then
       hydrogen_dt = star%job%timestep_override(kind_card_index)*seconds_per_year
       previous_timestep = hydrogen_dt
      else
! mhp 9/01  turn off structure-based timestep setting above a critical
!           temperature; this is done when
      if(star%job%use_structure_dt_limits) then
         if(log_temperature(1).gt.7.1d0 .and. luminosity_components(7).lt. 0.0d0) then
            write(*, 100) log_temperature(1), luminosity_components(7)
 100        format('LPTIME SET FALSE - TC ',F7.4,' EGRAV ',E10.2)
            star%job%use_structure_dt_limits = .false.
         endif
      endif
!  find timestep based on changes in structure variables from
!       one model to the next.
!  note that this returns the timestep stored in the model on the
!       first call to htimer for each kind card.
      if(star%job%use_structure_dt_limits) then
       call ptime(previous_timestep,luminosity,log_pressure,log_radius,log_temperature,num_points,structure_dt)
      else
       structure_dt = 1.0d20
      endif
!  mhp 6/90 a zero timestep in the model will fritz out the code, so
!  avoid this by setting the structure timestep arbitrarily large.
      if(structure_dt.lt.1.0d0) structure_dt=1.0d20
!  find timestep for rotating models based on changes in angular velocity
!       from one model to the next.
!  this sr also returns the maximum change in omega from the previous
!       model to the current one, which is used in the diffusion routine
!       to determine the number of diffusion timesteps required.
!  note that this returns the timestep stored in the model on the
!       first call to htimer for each kind card.
      if(star%job%rotation_active) then
       call wtime(previous_timestep,num_points,omega,rotation_dt,max_domega_frac)
      else
       rotation_dt = 1.0d20
      endif
!  mhp 6/90 a zero timestep in the model will fritz out the code, so
!  avoid this by setting the rotation timestep arbitrarily large.
      if(rotation_dt.lt.1.0d0) rotation_dt=1.0d20
!  find the timestep based on hydrogen burning.
!  determine the total luminosity due to hydrogen burning (in solar units).
      hydrogen_luminosity = luminosity_components(1) + luminosity_components(2) + &
           luminosity_components(3) + luminosity_components(4)
!  skip h-burning timestep criteria for stars without h burning.
      if(hydrogen_luminosity.gt.1.0d-34) then
       call xtime(log_density,composition,luminosity,enclosed_mass,shell_mass, &
            log_temperature,hydrogen_luminosity,convective_core_edge_zone, &
            h_shell_midpoint_zone,num_points,hydrogen_dt, &
            rate_pp,rate_he3_he3,rate_he3_he4,rate_c12_p,rate_c13_p, &
            rate_n14_p,rate_o16_p,rate_c13_alpha,rate_zero9,rate_c12_alpha, &
            rate_n14_alpha,rate_triple_alpha,rate_zero13, &
            frac_c12_alpha,frac_be7_electron)
      else
       hydrogen_dt = 1.0d20
      endif
!  find the timestep based on helium burning.
!  skip for stars without he burning.
      if(luminosity_components(5).gt.1.0d-34) then
       call ytime(energy_gen_terms,composition,log_density,luminosity, &
            enclosed_mass,log_temperature,convective_core_edge_zone, &
            num_points,helium_dt, &
            rate_pp,rate_he3_he3,rate_he3_he4,rate_c12_p,rate_c13_p, &
            rate_n14_p,rate_o16_p,rate_c13_alpha,rate_zero9,rate_c12_alpha, &
            rate_n14_alpha,rate_triple_alpha,rate_zero13, &
            frac_c12_alpha,frac_be7_electron,h_shell_zone_begin)
      else
       helium_dt = 1.0d20
      endif

!  04/14 jvs added timestep governor based on the size of the envelope
!  triangle. timesteps are restricted such that the model does not move
!  more than tri_delta_logl or tri_delta_teffl
      if(star%ctrl%use_envelope_triangle_dt .and. .not. star%job%use_structure_dt_limits) then
       call entime(previous_timestep,luminosity,log_teff,num_points,envelope_dt)
      else
       envelope_dt = 1.0d20
      endif

!     limit increase in time step from one model to the next for models
!     with a non-zero timestep.
!      write(*,299)structure_dt,rotation_dt,hydrogen_dt,helium_dt,previous_timestep
! 299  format('P ',E10.5,' W ',E10.5,' H ',E10.5,' He ',E10.5,' prev. ',
!     * E10.5)
!         WRITE(*,*)'H time ',hydrogen_dt/seconds_per_year,'   He Time ',helium_dt/seconds_per_year

      if(previous_timestep.gt.1.0d0) then
         previous_timestep = star%ctrl%atime(13)*previous_timestep
!  now set the timestep to be the minimum of the entropy based timestep,
!  the nuclear burning timesteps, and the previous timestep*atime(13).
!  04/14 jvs added envelope_dt
         previous_timestep = min(previous_timestep,hydrogen_dt,helium_dt,rotation_dt,structure_dt,envelope_dt)
      else
!  set the timestep to be the minimum of the the above timesteps,
!  but neglect the previous timestep.
!  04/14 jvs added envelope_dt
         previous_timestep = min(hydrogen_dt,helium_dt,rotation_dt,structure_dt,envelope_dt)
      endif

! c ******** Grant added ********
! c      IF(TLUMX(5).GT.1.0D-2) THEN
! c         ATIME(6)=0.000025
! c         ATIME(5)=0.01
! c       IF(TLUMX(5).GT.5.0D-1) THEN
! c         ATIME(6)=0.0000075
! c         ATIME(5)=0.003
! c       IF(TLUMX(5).GT.5.0D0) THEN
! c         ATIME(6)=0.00000075
! c         ATIME(5)=0.0003
! c       IF(TLUMX(5).GT.5.0D01) THEN
! c         ATIME(6)=0.000000025
! c         ATIME(5)=0.00001
! c       IF(TLUMX(5).GT.5.0D02) THEN
! c         ATIME(6)=0.000000005
! c         ATIME(5)=0.000005
! c       IF(TLUMX(5).GT.2.0D03) THEN
! c         ATIME(6)=0.0000000005
! c         ATIME(5)=0.0000005
! c         ENDIF
! c         ENDIF
! c         ENDIF
! c         ENDIF
! c         ENDIF
! c      ENDIF
! c *****************************

      hydrogen_dt = previous_timestep
      end if
      end if
      previous_timestep = abs(previous_timestep)
      timestep_years = previous_timestep/seconds_per_year
!     mhp 10/24 flag includes other stop conditions,
!  so only do this if there is a true stop age.
!     if evolving to a given age, ensure age not exceeded
!      IF(LENDAG(NK)) THEN
      if(star%job%end_age_stop_active(kind_card_index) .and. star%job%target_end_age(kind_card_index).gt.0.0d0) then
       time_left_years = star%job%target_end_age(kind_card_index) - age_gyr*1.0d9
       if(time_left_years.lt. timestep_years) then
         timestep_years = time_left_years
         hydrogen_dt = timestep_years*seconds_per_year
         previous_timestep = hydrogen_dt
       endif
      endif


      return
end subroutine htimer
