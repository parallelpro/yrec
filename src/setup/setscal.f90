!----------------------------------------------------------------------
! setscal
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original setscal.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Sets up the run-list to evolve a model to a target point on the
! HR diagram (target luminosity plus either target Teff or target
! R/Rsun, per common/calstar/'s specify_teff_flag): expands the two
! seed runs into 50 runs by copying their rescaling/mixing-length
! parameters forward (odd runs from 3 rerun the seed rescaling; even
! runs from 4 reuse run 2's age/timestep stop criteria), mirroring
! setcal.f90's solar-calibration run-list expansion.
subroutine setscal

      use star_info_lib, only: star
      use luout_lib
      use const_lib
      implicit none
!     LSTAR     T - have got a star at Teff and L
!     LPASSR    T - on run have just passed Teff
!     XLS       Luminosity (L/Lsun) wanted by adjusting Y
!     XLSTOL    tolerance wanted for luminosity
!     LTEFF     T - specify Teff for star
!               F - specify R/Rsun for star
!     STEFF     Effective temperature of star (K) or...
!     SR        Radius of star (R/Rsun)
!     TEFF      Teff of current model
!     ALR       log(R/Rsun) of current model
!     ALRI      log(R/Rsun) of previous model
!     DAGE      age of current model (Gyr)
!     AGEI      age of previous model (Gyr)
!     AGER      age of model at R*
!     BL        luminosity of current model
!     BLI       luminosity of previous model
!     BLR       luminosity of model at R
!     BLRP      luminosity of model at R* of previous run
!     XP        X of previous run = RESCAL(2, NK-1)
! mhp 10/02 linct not used, commented out
!      LINCT = .TRUE.

! --- locals ---
      integer :: i, j

      star_found_flag = .false.
      just_passed_target_radius_flag = .false.
      if (specify_teff_flag) then
         target_radius_rsun = sqrt(target_luminosity_lsun* &
              star%solar_luminosity_cgs/(c4pi*csig))/(target_teff*target_teff* &
              star%solar_radius_cgs)
      else
         target_teff = ((target_luminosity_lsun*star%solar_luminosity_cgs)/ &
              (c4pi*csig*target_radius_rsun*target_radius_rsun* &
              star%solar_radius_cgs*star%solar_radius_cgs))**0.25d0
      end if
      log_r_prev_model = 0
!     SET UP RUN TO EVOLVE TO L, Teff IN HR-DIAGRM.
!     THIS CONSISTS OF SETTING THE NUMBER OF RUNS TO THE MAXIMUM (50),
!     AND COPYING THE RELEVANT PARAMETERS FROM THE FIRST TWO RUNS TO
!     THE NEXT SERIES OF 24 CALIBRATING RUNS.
      num_runs = 50
      do i = 2,50
         initial_x_array(i) = initial_x_array(1)
         initial_z_array(i) = initial_z_array(1)
         mixing_length_array(i) = mixing_length_array(1)
         has_senv0_array(i) = has_senv0_array(1)
         senv0_array(i) = senv0_array(1)
      end do
      do i = 3,49,2
         rescale_kind(i) = rescale_kind(1)
         first_call_flag(i) = .true.
         num_models(i) = num_models(1)
         rsclzc(i) = rsclzc(1)
         rsclzm1(i) = rsclzm1(1)
         rsclzm2(i) = rsclzm2(1)
         do j = 1,4
            rescale_params(j,i) = rescale_params(j,1)
         end do
      end do
      do i = 4,50,2
         rescale_kind(i) = 1
         first_call_flag(i) = .false.
         num_models(i) = num_models(2)
         target_end_age(i) = target_end_age(2)
         end_age_stop_active(i) = end_age_stop_active(2)
         timestep_override(i) = timestep_override(2)
         timestep_override_active(i) = timestep_override_active(2)
      end do
      write(*,*) ' Evolve to R*, L* = ', target_radius_rsun, &
           target_luminosity_lsun
      write(itrack,*) '#Evolve to R*, L* = ', target_radius_rsun, &
           target_luminosity_lsun
      return
end subroutine setscal
