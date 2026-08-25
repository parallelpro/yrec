!----------------------------------------------------------------------
! stitch
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original stitch.f; only variable names, source form, and comment
! style were updated.
!
! STITCH: an alternate file format for .store that provides profiles
! for each desired model. Stitches the envelope and atmosphere
! solutions onto the interior when LSTENV and LSTATM are true. Will
! not provide atmosphere information when atmosphere tables are used.
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

      use atm_lib
      use envint_lib, only: atm_get
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

      double precision :: dum1(4), dum2(3), dum3(3), dum4(3)
      double precision :: envs1(json)
! --- locals ---
      integer :: i, j, k
      double precision :: cg, sg, fm, duma, a_val, rpoleq, vtot
      double precision :: abeg0, amin0, amax0, ebeg0, emin0, emax0
      integer :: idum
      double precision :: b, fpl, ftl
      integer :: katm, kenv, ksaha, ixx
      logical :: lprt, lsbc0
      double precision :: x, z, rl, gl, plim, ateffl
      logical :: lpulpt
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
            write(istor,63,advance='no') star%so(i),sg,star%del_grad(i_grad_rad,i),star%del_grad(i_grad_actual,i), &
                 star%del_grad(i_grad_ad,i),star%svel(i),star%adiabatic_index_gamma1(i), &
                 star%sfxion(1,i),star%sfxion(2,i),star%sfxion(3,i), &
                 star%sbeta(i),star%seta(i),(star%seg(k,i),k=1,5),star%sesum(i),star%seg(i_eps_neu,i),star%seg(i_eps_grav,i), &
                 star%scp(i),star%pulse_dlnrho_dlnt(i)
! write out additional rotation info if rotation is on
            if(star%job%rotation_active)then
              fm = dexp(ln10*log_mass(i))
              duma = cc13*omega(i)**2/(cg*fm)*5.d0/(2.d0+rotation_eta2(i))
              a_val = duma * radius_ratio_r0(i)**3
              rpoleq = (1.0d0 - a_val)/(1.0d0 + 0.5d0*a_val)
              vtot = star%circ%es_circulation_velocity(i)+star%circ%gsf_circulation_velocity(i)+ &
                   star%circ%secular_shear_velocity(i)
              write(istor,64) a_val,rpoleq,shape_factor_fp(i), &
                   shape_factor_ft(i),specific_angular_momentum(i), &
                   shell_moment_of_inertia(i),star%rot%rotational_energy_term(i), &
                   star%circ%es_circulation_velocity(i),star%circ%gsf_circulation_velocity(i), &
                   star%circ%secular_shear_velocity(i),vtot
            else
               write(istor,64) 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
            end if
         end do


! **************************   WRITE OUT ENVELOPE INFORMATION   **********************

      if(star%ctrl%lstenv)then ! only provide an envelope if asked to do so
! Begin by "dropping a sinkline" with the envelope integrator
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
      idum = 0
      b = dexp(ln10*log_luminosity_lsun)
      fpl = shape_factor_fp(m)
      ftl = shape_factor_ft(m)
      katm = 0
      kenv = 0
      ksaha = 0
      ixx=0
      lprt=.true.
      lsbc0 = .false.
      x = composition(1,m)
      z = composition(3,m)
      rl = 0.5d0*(log_luminosity_lsun + star%log10_solar_luminosity - &
           4.0d0*log_teff - c4pil - csigl)
      gl = cgl + log_total_mass - rl - rl
      plim = log_pressure(m)
! G Somers 10/14, FOR SPOTTED RUNS, FIND THE
! PRESSURE AT THE AMBIENT TEMPERATURE ATEFFL
      if(convective_flag(m).and.star%ctrl%spot_filling_factor.ne.0.0.and. &
           star%ctrl%spot_temp_contrast.ne.1.0)then
          ateffl = log_teff - 0.25*log10(star%ctrl%spot_filling_factor* &
               star%ctrl%spot_temp_contrast**4.0+1.0-star%ctrl%spot_filling_factor)
      else
          ateffl = log_teff
      end if
      call atm_get(b,fpl,ftl,gl,log_total_mass,ixx,lprt,lsbc0, &
         plim,rl,ateffl,x,z,dum1,idum,katm,kenv,ksaha, &
         dum2,dum3,dum4,lpulpt)

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
! Finish with the atmosphere, if the atmosphere was computed
       if(star%job%lstatm)then
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
