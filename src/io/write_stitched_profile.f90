!----------------------------------------------------------------------
! stitch
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original stitch.f; only variable names, source form, and comment
! style were updated.
!
! STITCH: an alternate file format for .store that provides profiles
! for each desired model. Splices the envelope and atmosphere
! solutions onto the interior when LSTENV and LSTATM are true. Will
! not provide atmosphere information when atmosphere tables are used
! (atmo_struct%num_atm_points is 0 then).
!
! 2026 (.store convergence): pure reader. The envelope/atmosphere
! come from env_struct/atmo_struct as build_stitched_model filled
! them for the converged model (evolve_step runs the stitch before
! any output); this writer no longer integrates anything itself.
!
! The output columns in the new .store format are:
! 1 MODEL, 2 SHELL, 3 log(mass[g]), 4 log(r[cm]), 5 L/Lsun, 6 log(P[cgs]), 7 log(T[K])',
! 8 log(DENSITY[cgs]),9 OMEGA(rad/s),10 CONVECTIVE?, 11 INTERIOR_POINT?, 12 ENVLELOPE_PT?
! 13 ATMOSPHERE_POINT?, 14 H1(mass frac), 15 He4(mass frac),16 METALS(mass frac),
! 17 He3(mass frac), 18 C12(mass frac), 19 C13(mass frac), 20 N14(mass frac),
! 21 N15(mass frac), 22 O16(mass frac), 23 O17(mass frac), 24 O18(mass frac),
! 25 H2(mass frac), 26 Li6(mass frac),27 Li7(mass frac),28 Be9(mass frac),29 OPACITY[cgs]
! 30 GRAVITY(cgs), 31 DELR(Rad. temp. grad), 32 DEL(actual temp grad),
! 33 DELA(adiabatic temp grad), 34 CONVECTIVE _VELOCITY[cm/s],35 GAM1(adiabatic exponent),
! 36 HII, 37 HEII, 38 HEIII, 39 BETA, 40 ETA, 41 PPI, 42 PPII, 43 PPIII, 44, CNO, 45 3HE
! 46 E_NUC,47 E_NEU,48 E_GRAV,49 Cp,50 dlnrho/dlnT,51 A, 52 RP/RE, 53 FP, 54 FT, 55 J/M,
! 56 MOMENT, 57 DEL_KE, 58 V_ES, 59 V_GSF, 60 V_SS, 61 VTOT   '
subroutine write_stitched_profile(composition, log_radius, log_pressure, log_density, &
     log_mass, log_temperature, log_luminosity, mass_coordinate, omega, &
     rotation_eta2, shell_moment_of_inertia, radius_ratio_r0, &
     specific_angular_momentum, shape_factor_fp, shape_factor_ft, &
     log_teff, log_total_mass, log_luminosity_lsun, m, convective_flag, &
     model)
      use rotation_scratch_lib

      use star_info_lib, only: star, i_eps_grav, i_eps_neu, i_grad_actual, i_grad_ad, i_grad_rad, json
      use atmstruct_lib
      use envstruct_lib
      use luout_lib
      use phys_const_lib

      implicit none

      double precision, intent(in) :: composition(15,json)
      double precision, intent(in) :: log_radius(json), log_pressure(json), &
           log_density(json), log_mass(json), log_temperature(json), &
           log_luminosity(json), mass_coordinate(json)
      double precision, intent(in) :: omega(json), rotation_eta2(json), &
           shell_moment_of_inertia(json), radius_ratio_r0(json), &
           specific_angular_momentum(json)
      double precision, intent(in) :: shape_factor_fp(json), &
           shape_factor_ft(json)
      double precision, intent(in) :: log_teff, log_total_mass, &
           log_luminosity_lsun
      integer, intent(in) :: m
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: model

      double precision :: envs1(json)
! --- locals ---
      integer :: i, j, k
      double precision :: cg, sg, fm, duma, a_val, rpoleq, vtot
      double precision :: b
      double precision :: rad

!
!
! STITCH: and alternate file format for .store that provides profiles for each
! desired model. Stitches the envelope and atmosphere solutions onto the interior
! when LSTENV and LSTATM are true. Will not provide atmosphere information when
! atmosphere tables are used.
!
! ****************************  WRITE OUT INTERIOR INFORMATION   **********************

      cg=dexp(ln10*cgl)
      do i = 1,m
            sg = dexp(ln10*(cgl - 2.0d0*log_radius(i)))*mass_coordinate(i)
! write out the basic info
            write(istor,62,advance='no') model,i,log_mass(i),log_radius(i), &
                 log_luminosity(i),log_pressure(i),log_temperature(i), &
                 log_density(i),omega(i),convective_flag(i),.true.,.false., &
                 .false.,(composition(j,i),j=1,15)
! write out additional physics if desired
            write(istor,63,advance='no') star%opacity_zone(i),sg,star%gradr(i),star%gradT(i), &
                 star%grada(i),star%conv_vel(i),star%adiabatic_index_gamma1(i), &
                 star%fxion_zone(1,i),star%fxion_zone(2,i),star%fxion_zone(3,i), &
                 star%beta(i),star%eta(i),(star%eps_channels(k,i),k=1,5),star%eps_total(i),star%eps_channels(i_eps_neu,i),star%eps_channels(i_eps_grav,i), &
                 star%scp(i),star%pulse_dlnrho_dlnt(i)
! write out additional rotation info if rotation is on
            if(star%job%rotation_active)then
              fm = dexp(ln10*log_mass(i))
              duma = cc13*omega(i)**2/(cg*fm)*5.d0/(2.d0+rotation_eta2(i))
              a_val = duma * radius_ratio_r0(i)**3
              rpoleq = (1.0d0 - a_val)/(1.0d0 + 0.5d0*a_val)
              vtot = star%es_circulation_velocity(i)+star%gsf_circulation_velocity(i)+ &
                   star%secular_shear_velocity(i)
              write(istor,64) a_val,rpoleq,shape_factor_fp(i), &
                   shape_factor_ft(i),specific_angular_momentum(i), &
                   shell_moment_of_inertia(i),rot_scr%rotational_energy_term(i), &
                   star%es_circulation_velocity(i),star%gsf_circulation_velocity(i), &
                   star%secular_shear_velocity(i),vtot
            else
               write(istor,64) 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
            end if
         end do


! **************************   WRITE OUT ENVELOPE INFORMATION   **********************

! 2026 (.store convergence): this writer no longer "drops a
! sinkline" of its own -- build_stitched_model already re-integrated
! the envelope/atmosphere at the converged model with the fixed
! output step sizes, in evolve_step before any output, and
! env_struct/atmo_struct hold that integration. Reading them here
! makes the writer a pure reader and also fixes a long-standing bug:
! the old sinkline overrode the adaptive envelope/atmosphere step
! sizes with the fixed output ones and never restored them, so from
! the first .store write onward the SOLVER's envelope integrations
! silently ran at output resolution.
      b = dexp(ln10*log_luminosity_lsun)
      if(star%ctrl%lstenv)then ! only provide an envelope if asked to do so
! DEFINE SOME ARRAYS WE NEED
      do i=1,env_struct%num_env_points
          envs1(i) = dexp(ln10*(env_struct%env_log10_mass(i)+log_total_mass))
      end do
         do i=m+1,m+env_struct%num_env_points
            sg = dexp(ln10*(cgl - 2.0d0*env_struct%env_log10_radius(i-m)))*envs1(i-m)
! write out the basic info. Omega and abundances take value of last interior point.
            write(istor,62,advance='no') model,i,env_struct%env_log10_mass(i-m)+log_total_mass, &
            env_struct%env_log10_radius(i-m),env_struct%env_luminosity(i-m), &
            env_struct%env_log10_pressure(i-m),env_struct%env_log10_temperature(i-m), &
            env_struct%env_log10_density(i-m),omega(m),env_struct%env_convective_flag(i-m), &
            .false.,.true.,.false., &
            (composition(j,m),j=1,15)
! write out additional physics
               write(istor,63,advance='no') env_struct%env_opacity(i-m),sg,env_struct%env_gradients(1,i-m), &
                 env_struct%env_gradients(2,i-m), &
                 env_struct%env_gradients(3,i-m),env_struct%env_convective_velocity(i-m), &
                 env_struct%env_gamma1(i-m),env_struct%env_ion_fraction(1,i-m),env_struct%env_ion_fraction(2,i-m), &
                 env_struct%env_ion_fraction(3,i-m),env_struct%env_beta(i-m),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, &
                 env_struct%env_specific_heat_cp(i-m),env_struct%env_dlnrho_dlnt(i-m)
! zero out rotation columns for envelope
               write(istor,64) 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
         end do
       end if



! *************************** WRITE OUT ATMOSPHERE INFORMATION  ************************
! Finish with the atmosphere, if the atmosphere was computed.
! num_atm_points is 0 when the surface boundary came from a
! tabulated atmosphere (no gray integration to stitch) -- that count
! replaces the old side effect in envint that forced LSTATM off.
       if(star%job%lstatm .and. atmo_struct%num_atm_points.gt.0)then
            do i=atmo_struct%num_atm_points,1,-1
! write out the basic info. Omega and abundances take value of last interior point.
            rad = dlog10(dexp(ln10*env_struct%env_log10_radius(env_struct%num_env_points)) + &
                 atmo_struct%atmo_delta_depth(i))
            write(istor,62,advance='no') model,atmo_struct%num_atm_points-i+m+env_struct%num_env_points, &
                 log_total_mass,rad,b, &
            atmo_struct%atmo_log10_pressure(i),atmo_struct%atmo_log10_temperature(i), &
            atmo_struct%atmo_log10_density(i),omega(m), &
            .false.,.false.,.false.,.true.,(composition(j,m),j=1,15)
! write out additional physics
            write(istor,63,advance='no') atmo_struct%atmo_opacity(i),sg,atmo_struct%atmo_gradients(1,i), &
                 atmo_struct%atmo_gradients(2,i),atmo_struct%atmo_gradients(3,i),0.0,atmo_struct%atmo_gamma1(i), &
                 atmo_struct%atmo_ion_fraction(1,i),atmo_struct%atmo_ion_fraction(2,i), &
                 atmo_struct%atmo_ion_fraction(3,i),atmo_struct%atmo_beta(i),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, &
                 atmo_struct%atmo_specific_heat_cp(i),atmo_struct%atmo_dlnrho_dlnt(i)
!  zero placeholders for rotation output
            write(istor,64) 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
         end do
       end if



! **************************    Output format codes   ******************************

 62   format(I6,1X,I6,0P2F18.14,1PE24.16,0P3F18.14,1PE24.16,1X,L1,1X,L1,1X,L1,1X,L1, &
     &     3(0PF12.9),12(0PE16.8),2X)

 63   format(1PE10.4,1PE11.3,E12.4,E12.4,E12.4,1PE12.4,0PF9.5,F9.5,F9.5,F9.5, &
!!     &     F9.5,F9.5,F9.5,F9.5,F9.5,F9.5,F9.5,E13.5,E13.5,E13.5)
!     &     F9.5,F9.5,E12.4,E12.4,E12.4,E12.4,E12.4,E13.5,E13.5,E13.5,E13.5)
     &  F9.5,F9.5,E12.4,E12.4,E12.4,E12.4,E12.4,E13.5,E13.5,E13.5,E13.5,E13.5,E13.5,E13.5, &
     &  E13.5)
 64   format(E14.6,E14.6,E14.6,E14.6,E14.6,E14.6,E14.6,E11.3,E11.3,E11.3,E11.3)



      return
end subroutine write_stitched_profile
