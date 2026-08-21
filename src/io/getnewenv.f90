!----------------------------------------------------------------------
! getnewenv
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original getnewenv.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! MHP 07/02 new subroutine to move the outer fitting point when
! the total mass is being decreased; add new points by dropping
! a sinkline from the surface down and including new points up to
! the desired envelope mass
subroutine getnewenv(target_envelope_mass, composition, log_density, &
     log_luminosity, log_pressure, log_radius, log_mass, enclosed_mass, &
     shell_mass, &
     log_total_mass, log_temperature, convective_flag, eta_squared, &
     moment_of_inertia, specific_angular_momentum, qiw, mean_radius, &
     rotational_kinetic_energy, log_luminosity_lsun, total_angular_momentum, &
     total_rotational_ke, log_teff, num_zones, new_points_added_flag)
      use run_diag_lib
      use envstruct_lib
      use envelope_comp_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(inout) :: target_envelope_mass
      double precision, intent(inout) :: composition(15,json), &
           log_density(json), log_luminosity(json), log_pressure(json), &
           log_radius(json), log_mass(json), enclosed_mass(json)
      double precision, intent(inout) :: shell_mass(json)
      double precision, intent(in) :: log_total_mass
      double precision, intent(inout) :: log_temperature(json)
      logical, intent(inout) :: convective_flag(json)
      double precision, intent(inout) :: eta_squared(json), &
           moment_of_inertia(json), specific_angular_momentum(json), &
           qiw(json), mean_radius(json)
      double precision, intent(inout) :: rotational_kinetic_energy(json)
      double precision, intent(in) :: log_luminosity_lsun
      double precision, intent(inout) :: total_angular_momentum, &
           total_rotational_ke
      double precision, intent(in) :: log_teff
      integer, intent(inout) :: num_zones
      logical, intent(out) :: new_points_added_flag













      double precision :: envint_dummy1(4), envint_dummy2(3), &
           envint_dummy3(3), envint_dummy4(3)

! G Somers END
      save

      integer :: species_end_index
      integer :: old_num_zones
      double precision :: envelope_mass_before
      double precision :: env_max_saved, env_min_saved, env_begin_saved
      logical :: surface_bc_flag, print_flag
      integer :: katm, kenv, ksaha
      double precision :: luminosity_linear
      double precision :: log_radius_surface, log_gravity_surface
      double precision :: hydrogen_fraction, metal_fraction
      double precision :: fp_surface, ft_surface
      double precision :: pressure_limit
      logical :: pulsation_output_flag
      integer :: ixx_flag
      integer :: envint_unused_flag
      double precision :: spot_adjusted_log_teff
      double precision :: pressure_offset, density_offset, temperature_offset, &
           radius_offset
      double precision :: mass_test, mass_test2
      integer :: num_new_env_points
      integer :: zone_index, species_index, k
      double precision :: mass_interp_x0, mass_interp_x1, mass_interp_x2
      double precision :: interp_fraction
      double precision :: omega_ref
      double precision :: sum_angular_momentum, sum_rotational_ke
      double precision :: angular_momentum_shell
      double precision :: mass_at_base
! omega (originally OMEGA) is declared only via the original file's
! DIMENSION statement, never appearing in the SUBROUTINE argument
! list -- so, unlike every other file in this project, it is genuine
! local scratch here, not a dummy argument aliasing the caller's
! angular-velocity array. Because of the blanket SAVE below, its
! value nonetheless persists across calls to this subroutine within
! one run (matching the original's behavior, however accidental).
! Preserved exactly as such: omega(old_num_zones) is read in the
! WALPCZ.GE.0 and general-case branches below before this call ever
! writes it, using whatever this array held from a previous call (or
! its initial zero-fill) -- not the caller's actual rotation state.
      double precision :: omega(json)

      if(use_extended_composition)then
         species_end_index = 15
      else
         species_end_index = 11
      endif
! RESTRICT THE ENVELOPE MASS TO A MINIMUM OF 10**-12.
      target_envelope_mass = min(target_envelope_mass,-1.0D-12)
! SAVE CURRENT VALUES OF THE TOTAL NUMBER OF POINTS AND ENVELOPE MASS.
      old_num_zones = num_zones
      envelope_mass_before = env_comp%senv
! SET NUMERICAL PARAMETERS OF THE ENVELOPE INTEGRATION
      env_max_saved = env_step_max
      env_min_saved = env_step_min
      env_begin_saved = env_step_begin
      env_step_max = chi_grid_scale(8)
      env_step_min = chi_grid_scale(8)
      env_step_begin = chi_grid_scale(8)
      surface_bc_flag = .false.
      print_flag = .true.
      katm = 0
      kenv = 0
      ksaha = 0
! SET THE PHYSICAL PROPERTIES OF THE ENVELOPE SOLUTION
! LUMINOSITY
      luminosity_linear = dexp(ln10*log_luminosity_lsun)
! RADIUS
      log_radius_surface = 0.5D0*(log_luminosity_lsun + log10_solar_luminosity &
           - 4.0D0*log_teff - c4pil - csigl)
! SURFACE GRAVITY
      log_gravity_surface = cgl + env_comp%stotal - log_radius_surface - log_radius_surface
! COMPOSITION
      hydrogen_fraction = composition(1,num_zones)
      metal_fraction = composition(3,num_zones)
! NEGLECT ROTATIONAL DISTORTION
      fp_surface = 1.0D0
      ft_surface = 1.0D0
! PREVENT THE INTEGRATION FROM SHOOTING PAST THE
! PRESSURE OF THE OUTER MODEL POINT, AT LEAST THE
! FIRST TIME IT TRIES TO DO SO.
      pressure_limit = log_pressure(num_zones)
! DO NOT DO SOLAR PULSATION OUTPUT
      pulsation_output_flag = .false.
! SET UP VALUES FOR THE EQUATION OF STATE CALCULATION
      ixx_flag = 0
      if (use_debye_huckel_correction) then
         debye_huckel_x = composition(1,num_zones)
         debye_huckel_y = composition(2,num_zones)+composition(4,num_zones)
         debye_huckel_z_total = composition(3,num_zones)
         debye_huckel_z(1) = composition(5,num_zones)+composition(6,num_zones)
         debye_huckel_z(2) = composition(7,num_zones)+composition(8,num_zones)
         debye_huckel_z(3) = composition(9,num_zones)+composition(10,num_zones)+ &
              composition(11,num_zones)
      end if
! INTEGRATE DOWN TO THE CURRENT FITTING POINT USING THE SURFACE L AND TEFF.
! MHP 10/02  define ISTORE - used in ENVINT
      envint_unused_flag = 0
! G Somers 10/14, FOR SPOTTED RUNS, FIND THE
! PRESSURE AT THE AMBIENT TEMPERATURE ATEFFL
      if(convective_flag(num_zones).and.spot_filling_factor.ne.0.0.and. &
           spot_temp_contrast.ne.1.0)then
         spot_adjusted_log_teff = log_teff - 0.25*log10(spot_filling_factor * &
              spot_temp_contrast**4.0 + 1.0 - spot_filling_factor)
      else
         spot_adjusted_log_teff = log_teff
      endif
      call envint(luminosity_linear,fp_surface,ft_surface,log_gravity_surface, &
           log_total_mass,ixx_flag,print_flag,surface_bc_flag,pressure_limit, &
           log_radius_surface, &
           spot_adjusted_log_teff,hydrogen_fraction,metal_fraction, &
           envint_dummy1,envint_unused_flag,katm,kenv,ksaha,envint_dummy2, &
           envint_dummy3,envint_dummy4,pulsation_output_flag)
! G Somers END
! RESET THE NUMERICAL PARAMETERS FOR THE ENVELOPE INTEGRATION
      env_step_max = env_max_saved
      env_step_min = env_min_saved
      env_step_begin = env_begin_saved
      env_comp%senv = target_envelope_mass
! STOP IF THE DESIRED NUMBER OF POINTS EXCEEDS THE ARRAY DIMENSIONS
      if(num_zones+env_struct%num_env_points.ge.json) stop 9999
! THE FIRST POINT IN THE ENVELOPE SOLUTION IS THE SET OF PROPERTIES
! OBTAINED AT A MASS EQUAL TO THE LAST INTERIOR POINT.  USE THIS TO
! ENFORCE CONSISTENCY WITH THE INTERIOR SOLUTION;
! ADJUST THE (P, RHO, T, R) POINTS TO BE CONSISTENT
! WITH THE LAST MODEL POINT AND SHIFT ALL OF THE POINTS DOWN BY ONE.
      pressure_offset = log_pressure(num_zones) - env_struct%env_log10_pressure(1)
      density_offset = log_density(num_zones) - env_struct%env_log10_density(1)
      temperature_offset = log_temperature(num_zones) - env_struct%env_log10_temperature(1)
      radius_offset = log_radius(num_zones) - env_struct%env_log10_radius(1)
      mass_test = log_mass(num_zones)
      num_new_env_points = 0
      do zone_index = 1,env_struct%num_env_points - 1
         mass_test2 = env_comp%stotal+env_struct%env_log10_mass(zone_index+1)
         if(mass_test2-mass_test.gt.1.0D-10)then
            env_struct%env_log10_density(zone_index) = env_struct%env_log10_density(zone_index+1)+density_offset
            env_struct%env_log10_pressure(zone_index) = env_struct%env_log10_pressure(zone_index+1)+pressure_offset
            env_struct%env_log10_radius(zone_index) = env_struct%env_log10_radius(zone_index+1)+radius_offset
            env_struct%env_log10_mass(zone_index) = env_struct%env_log10_mass(zone_index+1)
            env_struct%env_log10_temperature(zone_index) = env_struct%env_log10_temperature(zone_index+1)+temperature_offset
            env_struct%env_hydrogen_fraction(zone_index) = env_struct%env_hydrogen_fraction(zone_index+1)
            env_struct%env_metal_fraction(zone_index) = env_struct%env_metal_fraction(zone_index+1)
            mass_test = mass_test2
            num_new_env_points = num_new_env_points + 1
         endif
      end do
      env_struct%num_env_points = num_new_env_points
! ASSIGN NEW POINTS
      do zone_index = num_zones+1,num_zones+env_struct%num_env_points
         species_index = zone_index-num_zones
! LUMINOSITY ASSUMED CONSTANT
         log_luminosity(zone_index) = log_luminosity(num_zones)
! INCLUDE NEW POINTS UP TO THE DIFFERENT DESIRED FITTING POINT
         if(env_struct%env_log10_mass(species_index).lt.env_comp%senv)then
            log_density(zone_index) = env_struct%env_log10_density(species_index)
            log_pressure(zone_index) = env_struct%env_log10_pressure(species_index)
            log_radius(zone_index) = env_struct%env_log10_radius(species_index)
            log_mass(zone_index) = env_struct%env_log10_mass(species_index) + env_comp%stotal
            log_temperature(zone_index) = env_struct%env_log10_temperature(species_index)
            composition(1,zone_index) = env_struct%env_hydrogen_fraction(species_index)
            composition(3,zone_index) = env_struct%env_metal_fraction(species_index)
            do k = 4,species_end_index
               composition(k,zone_index) = composition(k,num_zones)
            end do
            composition(2,zone_index)=1.0D0-composition(1,zone_index)- &
                 composition(3,zone_index)-composition(4,zone_index)
            convective_flag(zone_index) = env_struct%env_convective_flag(species_index)
         else
! POINTS BEYOND THIS ARE ABOVE THE NEW DESIRED FITTING POINT;
! INTERPOLATE LINEARLY, SET NEW NUMBER OF TOTAL POINTS, AND EXIT
            if(species_index.eq.1)then
! INTERPOLATE BETWEEN THE LAST INTERIOR POINT AND THE FIRST ENVELOPE POINT
               mass_interp_x0 = log_mass(num_zones)
               mass_interp_x1 = env_comp%stotal + env_comp%senv
               mass_interp_x2 = env_struct%env_log10_mass(species_index) + env_comp%stotal
               if(mass_interp_x2-mass_interp_x0.lt.1.0D-14) stop 9998
               interp_fraction = (mass_interp_x1-mass_interp_x0)/ &
                    (mass_interp_x2-mass_interp_x0)
               log_density(zone_index) = log_density(num_zones)+interp_fraction* &
                    (env_struct%env_log10_density(species_index)-log_density(num_zones))
               log_pressure(zone_index) = log_pressure(num_zones)+interp_fraction* &
                    (env_struct%env_log10_pressure(species_index)-log_pressure(num_zones))
               log_radius(zone_index) = log_radius(num_zones)+interp_fraction* &
                    (env_struct%env_log10_radius(species_index)-log_radius(num_zones))
               log_mass(zone_index) = mass_interp_x1
               log_temperature(zone_index) = log_temperature(num_zones)+interp_fraction* &
                    (env_struct%env_log10_temperature(species_index)-log_temperature(num_zones))
               composition(1,zone_index) = composition(1,num_zones)+interp_fraction* &
                    (composition(1,num_zones)-env_struct%env_hydrogen_fraction(species_index))
               composition(3,zone_index) = composition(3,num_zones)+interp_fraction* &
                    (composition(3,num_zones)-env_struct%env_metal_fraction(species_index))
               do k = 4,species_end_index
                  composition(k,zone_index) = composition(k,num_zones)
               end do
               composition(2,zone_index)=1.0D0-composition(1,zone_index)- &
                    composition(3,zone_index)-composition(4,zone_index)
               if(env_struct%env_convective_flag(species_index).or.convective_flag(num_zones))then
                  convective_flag(zone_index) = .true.
               else
                  convective_flag(zone_index) = .false.
               endif
            else
! INTERPOLATE BETWEEN THE LAST 2 ENVELOPE POINTS
               mass_interp_x0 = env_struct%env_log10_mass(species_index-1) + env_comp%stotal
               mass_interp_x1 = env_comp%stotal + env_comp%senv
               mass_interp_x2 = env_struct%env_log10_mass(species_index) + env_comp%stotal
               if(mass_interp_x2-mass_interp_x0.lt.1.0D-14) stop 9998
               interp_fraction = (mass_interp_x1-mass_interp_x0)/ &
                    (mass_interp_x2-mass_interp_x0)
               log_density(zone_index) = env_struct%env_log10_density(species_index-1)+interp_fraction* &
                    (env_struct%env_log10_density(species_index)-env_struct%env_log10_density(species_index-1))
               log_pressure(zone_index) = env_struct%env_log10_pressure(species_index-1)+interp_fraction* &
                    (env_struct%env_log10_pressure(species_index)-env_struct%env_log10_pressure(species_index-1))
               log_radius(zone_index) = env_struct%env_log10_radius(species_index-1)+interp_fraction* &
                    (env_struct%env_log10_radius(species_index)-env_struct%env_log10_radius(species_index-1))
               log_mass(zone_index) = mass_interp_x1
               log_temperature(zone_index) = env_struct%env_log10_temperature(species_index-1)+interp_fraction* &
                    (env_struct%env_log10_temperature(species_index)-env_struct%env_log10_temperature(species_index-1))
               composition(1,zone_index) = env_struct%env_hydrogen_fraction(species_index-1)+interp_fraction* &
                    (env_struct%env_hydrogen_fraction(species_index)-env_struct%env_hydrogen_fraction(species_index-1))
               composition(3,zone_index) = env_struct%env_metal_fraction(species_index-1)+interp_fraction* &
                    (env_struct%env_metal_fraction(species_index)-env_struct%env_metal_fraction(species_index-1))
               do k = 4,species_end_index
                  composition(k,zone_index) = composition(k,num_zones)
               end do
               composition(2,zone_index)=1.0D0-composition(1,zone_index)- &
                    composition(3,zone_index)-composition(4,zone_index)
               if(env_struct%env_convective_flag(species_index).or.env_struct%env_convective_flag(species_index-1))then
                  convective_flag(zone_index) = .true.
               else
                  convective_flag(zone_index) = .false.
               endif
            endif
            num_zones = zone_index
            goto 587
         endif
      end do
! ASSIGN THE BOUNDARY AT THE PHOTOSPHERE FOR ENVELOPE MASS BELOW 1.0D-12.
      num_zones = num_zones + env_struct%num_env_points
 587  continue
! ADD THE UNLOGGED MASSES OF THE NEW SHELLS (HS1) AND COMPUTE THE
! MASS CONTENTS OF THE NEW SHELLS (HS2).
      do zone_index = old_num_zones,num_zones
         enclosed_mass(zone_index) = 10.0D0**log_mass(zone_index)
      end do
      do zone_index = old_num_zones,num_zones-1
         shell_mass(zone_index) = 0.5D0*(enclosed_mass(zone_index+1)-enclosed_mass(zone_index-1))
      end do
      mass_at_base = 0.5D0*(enclosed_mass(num_zones)+enclosed_mass(num_zones-1))
      shell_mass(num_zones) = 10.0D0**log_total_mass - mass_at_base
! RECOMPUTE TERMS RELATED TO ROTATION.
      if(rotation_active)then
! FIRST GUESS AT THE ROTATION RATES; ASSIGN A
! VECTOR OF OMEGA SUCH THAT
! OMEGA*R**WALPCZ = CONSTANT.
         if(walpcz.le.-2.0D0)then
! CONSTANT J/M
            do zone_index = old_num_zones+1,num_zones
               specific_angular_momentum(zone_index) = specific_angular_momentum(old_num_zones)
               moment_of_inertia(zone_index) = cc23*10.0D0**(2.0D0*log_radius(zone_index))
               omega(zone_index) = specific_angular_momentum(zone_index)/moment_of_inertia(zone_index)
            end do
         else if(walpcz.ge.0.0D0)then
! SOLID BODY ROTATION
            do zone_index = old_num_zones+1,num_zones
               omega(zone_index) = omega(old_num_zones)
               moment_of_inertia(zone_index) = cc23*10.0D0**(2.0D0*log_radius(zone_index))
               specific_angular_momentum(zone_index) = omega(old_num_zones)*moment_of_inertia(zone_index)
            end do
         else
! GENERAL CASE
            omega_ref = omega(old_num_zones)*10.0D0**(log_radius(old_num_zones)*walpcz)
            do zone_index = old_num_zones+1,num_zones
               omega(zone_index) = omega_ref/10.0D0**(log_radius(zone_index)*walpcz)
               moment_of_inertia(zone_index) = cc23*10.0D0**(2.0D0*log_radius(zone_index))
               specific_angular_momentum(zone_index) = omega(zone_index)*moment_of_inertia(zone_index)
            end do
         endif
         call getrot(log_density,specific_angular_momentum,log_radius,log_mass, &
              shell_mass,convective_flag,num_zones,eta_squared, &
              moment_of_inertia,omega,qiw,mean_radius)
! GIVEN OMEGA AND I, FIND ANGULAR MOMENTUM AND ROTATIONAL K.E.
       sum_angular_momentum = 0.0D0
       sum_rotational_ke = 0.0D0
       do 550 zone_index = 1,num_zones
! MHP 10/02 logic reversed!
!          HJM(I) = HJ/HS2(I)
          angular_momentum_shell = specific_angular_momentum(zone_index)*shell_mass(zone_index)
          rotational_kinetic_energy(zone_index) = 0.5D0*omega(zone_index)*angular_momentum_shell
          sum_angular_momentum = sum_angular_momentum + angular_momentum_shell
          sum_rotational_ke = sum_rotational_ke + rotational_kinetic_energy(zone_index)
 550     continue
       total_angular_momentum = sum_angular_momentum
       total_rotational_ke = sum_rotational_ke
      endif
      write(*,910)
 910  format(1X,'NEW INTERIOR POINTS FROM CHANGE IN ENVELOPE MASS'/ &
     ' J,LOG RHO, LOG L, LOG P, LOG R, LOG M, LOG T, CONV T/F')
      write(*,911)(zone_index,log_density(zone_index),log_luminosity(zone_index), &
           log_pressure(zone_index),log_radius(zone_index), &
           log_mass(zone_index)-env_comp%stotal, &
           log_temperature(zone_index),convective_flag(zone_index), &
           zone_index = old_num_zones,num_zones)
 911  format(I5,1P6E16.8,L2)
      new_points_added_flag = .true.
       write(*,597)envelope_mass_before,env_comp%senv
 597  format(5X,'***** NEW ENVELOPE MASS CALCULATED *****'/8X, &
     'OLD SENV ',1PE22.13,'  NEW SENV',E22.13)
      return
end subroutine getnewenv
