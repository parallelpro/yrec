!----------------------------------------------------------------------
! chkcal
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original chkcal.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Checks whether the current model is a converged solar calibration
! (log L and log R, and optionally log Z/X, within tolerance of the
! solar values) and, if not, applies a fixed 2x2 Newton correction to
! the envelope hydrogen fraction (X) and mixing-length alpha using
! pre-determined empirical partial derivatives of log L and log R with
! respect to X and alpha (dlum_dx/drad_dx/dlum_dalpha/drad_dalpha).
! ONLY CALLED FOR EVEN NK, ASSUMES RESCALING ON ODD NK AND EVOLVING
! ON EVEN NK
subroutine chkcal(log_l_lsun, log_r_rsun, run_index, current_zx)

      use const_lib
      use luout_lib
      implicit none

      double precision, intent(in) :: log_l_lsun, log_r_rsun
      integer, intent(in) :: run_index
      double precision, intent(in) :: current_zx


!      COMMON/SETT/ENDAGE(50),SETDT(50),LENDAG(50),LSETDT(50)

! common/newxym/: initial_x_array/initial_z_array/mixing_length_array
! are used here; has_senv0_array/senv0_array are unused placeholders.
! Naming matches wrtmonte.f90.
      double precision :: initial_x_array(50), initial_z_array(50), &
           mixing_length_array(50)
      logical :: has_senv0_array(50)
      double precision :: senv0_array(50)
      common /newxym/ initial_x_array, initial_z_array, &
           mixing_length_array, has_senv0_array, senv0_array

! common/calsun/: dlum_dx/drad_dx/dlum_dalpha/drad_dalpha/log_l_prev/
! log_r_prev/delta_x/delta_alpha/solar_calibration_active, all used
! here. Naming matches wrtmonte.f90.
      double precision :: dlum_dx, drad_dx, dlum_dalpha, drad_dalpha, &
           log_l_prev, log_r_prev, delta_x, delta_alpha
      logical :: solar_calibration_active
      common/calsun/ dlum_dx, drad_dx, dlum_dalpha, drad_dalpha, &
           log_l_prev, log_r_prev, delta_x, delta_alpha, &
           solar_calibration_active

! MHP 6/13 ADD OPTION TO CALIBRATE SOLAR Z/X, SOLAR Z/X, SOLAR AGE
! common/cals2/: luminosity_tolerance/radius_tolerance/zx_tolerance/
! calibrate_solar_zx/target_solar_zx are used here; calibrate_solar_model/
! target_solar_age are unused placeholders. Not referenced in any
! already-converted file.
      double precision :: luminosity_tolerance, radius_tolerance, &
           zx_tolerance
      logical :: calibrate_solar_model, calibrate_solar_zx
      double precision :: target_solar_zx, target_solar_age
      common/cals2/ luminosity_tolerance, radius_tolerance, zx_tolerance, &
           calibrate_solar_model, calibrate_solar_zx, target_solar_zx, &
           target_solar_age
!      COMMON/CALS2/TOLL,TOLR,LCALS


      save
!     DATA TOLL,TOLR/1.0D-5,1.0D-4/

! locals -- log_zx_mismatch/delta_z/log_zx_mismatch_prev (originally
! ZL/DZ/ZLP); ZL is set only in the early-return Z/X check below but
! (via the blanket SAVE, matching the original) is still read again
! from its previous-call value at the "New Z" report near label 1234.
      double precision :: log_zx_mismatch, delta_z, log_zx_mismatch_prev

! CHECK TO SEE IF THE MODEL IS A CALIBRATED SOLAR MODEL.
! IF NOT, ESTIMATE CORRECTIONS TO X AND ALPHA.
      if(abs(log_l_lsun).lt.luminosity_tolerance .and. &
           abs(log_r_rsun).lt.radius_tolerance)then
! ADD CHECK FOR Z/X
         if(calibrate_solar_zx .and.target_solar_zx.gt.0.0d0)then
            log_zx_mismatch = log10(current_zx)-log10(target_solar_zx)
            if(abs(log_zx_mismatch).lt.zx_tolerance)then
               solar_calibration_active = .true.
               goto 9999
            endif
         else
! CALIBRATED SOLAR MODEL.  SET UP OUTPUT FLAGS AND EXIT
            solar_calibration_active = .true.
            goto 9999
         endif
      endif

! Added code to use pre-determined partial derivatives    LLP  6/28/09
      dlum_dx = -3.78d0            ! empirical result:  -3.783    RMS error .070
      drad_dx = -0.89d0             ! empirical result:  -0.890    RMS error .048
      dlum_dalpha = 0.0139d0            ! empirical result:  +0.139    RMS error .0022
      drad_dalpha = -0.050d0            ! empirical result:  -0.0504   RMS error .0059
      goto 1234            ! Bypass partial derivative code
! mhp 5/96 added change to compute solar calibration for 3 kind cards
      if(run_index.eq.3)then
!     SET UP RUN TO DETERMINE DERIVATIVE OF L AND R WITH RESPECT TO X.
!
!     3.7 is empirical average dL/dX, so trial DX is BL / (dL/dX)   ! llp 6/18/09
!         DX = BL/3.7D0
!         DX = +.01
         rescale_params(2,run_index+1) = rescale_params(2,1)+delta_x
! STORE PREVIOUS L AND R.
         initial_x_array(run_index+1) = rescale_params(2,run_index+1)
         initial_x_array(run_index+2)=initial_x_array(run_index+1)
         initial_x_array(run_index+3)=initial_x_array(run_index+1)
         log_l_prev = log_l_lsun
         log_r_prev = log_r_rsun
         goto 9999
      else if(run_index.eq.6)then
!     EVALUATE DERIVATIVE OF L AND R WITH RESPECT TO X.
!         DLDX = (BL - BLP)/DX
!         DRDX = (RL - RLP)/DX
!         WRITE(IOWR,*) "DX,DLDX,BL,BLP,DRDX,RL,RLP: ",
!    *      DX,DLDX,BL,BLP,DRDX,RL,RLP
!     SET UP RUN TO DETERMINE DERIVATIVE OF L AND R WITH RESPECT TO ALPHA.

!     .042 is typical average dR/dA, so trial DA is RL / (dR/dA)   ! llp 6/18/09
!         DA = RL/0.042D0
         mixing_length_array(run_index+1) = mixing_length_array(1)+delta_alpha
         mixing_length_array(run_index+2) = mixing_length_array(run_index+1)
         mixing_length_array(run_index+3) = mixing_length_array(run_index+1)
         rescale_params(2,run_index+1) = rescale_params(2,1)
         initial_x_array(run_index+1)=initial_x_array(1)
         initial_x_array(run_index+2)=initial_x_array(run_index+1)
         initial_x_array(run_index+3)=initial_x_array(run_index+1)
         log_l_prev = log_l_lsun
         log_r_prev = log_r_rsun
         goto 9999
      else if(run_index.eq.9)then
!     EVALUATE DERIVATIVE OF L AND R WITH RESPECT TO ALPHA.
         dlum_dalpha = (log_l_lsun - log_l_prev)/delta_alpha
         drad_dalpha = (log_r_rsun - log_r_prev)/delta_alpha
         write(iowr,*) "DA,DLDA,BL,BLP,DRDA,RL,RLP: ", &
             delta_alpha,dlum_dalpha,log_l_lsun,log_l_prev,drad_dalpha, &
             log_r_rsun,log_r_prev
!     USE DERIVATIVES OF L AND R WITH RESPECT TO X AND ALPHA TO
!     GET IMPROVED GUESSES FOR ALPHA AND X.
         delta_alpha = ((log_l_lsun*drad_dx/dlum_dx-log_r_rsun)/ &
              (drad_dalpha-dlum_dalpha*drad_dx/dlum_dx))
         delta_x = -(log_l_lsun + dlum_dalpha*delta_alpha)/dlum_dx
         mixing_length_array(run_index+1) = mixing_length_array(1)+delta_alpha
         mixing_length_array(run_index+2) = mixing_length_array(run_index+1)
         mixing_length_array(run_index+3) = mixing_length_array(run_index+1)
         rescale_params(2,run_index+1) = rescale_params(2,1)+delta_x
         initial_x_array(run_index+1) = rescale_params(2,run_index+1)
         initial_x_array(run_index+2) = rescale_params(2,run_index+1)
         initial_x_array(run_index+3) = rescale_params(2,run_index+1)
         log_l_prev = log_l_lsun
         log_r_prev = log_r_rsun
         write(iowr,*) "New X, Old X, Calc DX: ", &
             rescale_params(2,run_index+1), rescale_params(2,1), delta_x
         write(iowr,*) "New A, Old A, Calc DA: ", &
             mixing_length_array(run_index+1), mixing_length_array(1), &
             delta_alpha
         goto 9999
!      ELSE
      endif   ! terrminate old partial derivative code

 1234 continue
!     USE DERIVATIVES OF L AND R WITH RESPECT TO X AND ALPHA TO
!     GET IMPROVED GUESSES FOR ALPHA AND X.
         delta_alpha = ((log_l_lsun*drad_dx/dlum_dx-log_r_rsun)/ &
              (drad_dalpha-dlum_dalpha*drad_dx/dlum_dx))
         delta_x = -(log_l_lsun + dlum_dalpha*delta_alpha)/dlum_dx
         mixing_length_array(run_index+1) = &
              mixing_length_array(run_index-2)+delta_alpha
         mixing_length_array(run_index+2) = mixing_length_array(run_index+1)
         mixing_length_array(run_index+3) = mixing_length_array(run_index+1)
         if(calibrate_solar_zx)then
            rescale_params(3,run_index+1) = &
                 rescale_params(3,run_index-2)*target_solar_zx/current_zx
            delta_z = rescale_params(3,run_index+1) - &
                 rescale_params(3,run_index-2)
            initial_z_array(run_index+1) = rescale_params(3,run_index+1)
            initial_z_array(run_index+2) = rescale_params(3,run_index+1)
            initial_z_array(run_index+3) = rescale_params(3,run_index+1)
         endif
         rescale_params(2,run_index+1) = rescale_params(2,run_index-2)+delta_x
         initial_x_array(run_index+1) = rescale_params(2,run_index+1)
         initial_x_array(run_index+2) = rescale_params(2,run_index+1)
         initial_x_array(run_index+3) = rescale_params(2,run_index+1)
         write(iowr,*) "New BL, Old BL, Delta BL: ", &
             log_l_lsun, log_l_prev, log_l_lsun-log_l_prev
         write(iowr,*) "New RL, Old RL, Delta RL: ", &
             log_r_rsun, log_r_prev, log_r_rsun-log_r_prev
         write(iowr,*) "New X, Old X, DX: ", &
             rescale_params(2,run_index+1), rescale_params(2,run_index-2), &
             delta_x
         write(iowr,*) "New A, Old A, DA: ", &
             mixing_length_array(run_index+1), &
             mixing_length_array(run_index-2), delta_alpha
         if(calibrate_solar_zx)then
            write(iowr,*) "New Z, Old Z, DZ: ", &
             rescale_params(3,run_index+1), rescale_params(3,run_index-2), &
             delta_z
         endif
         log_l_prev = log_l_lsun
         log_r_prev = log_r_rsun
         log_zx_mismatch_prev = log_zx_mismatch
         goto 9999
!      ENDIF
 9999 continue
      return
end subroutine chkcal
