!----------------------------------------------------------------------
! eq2mod
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original eq2mod.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! The diffusion equation for angular momentum transport is solved using
! an equally spaced grid for each unstable region. Once the equation
! has been solved, we need the new run of angular momentum at the
! original model points.
!
! eq2mod transforms back to the original grid and calculates the
! change in specific angular momentum in each shell.
!
! Input variables:
!
!   Variables pertaining to the equally spaced grid points used to
!   solve the diffusion equation:
!
!   delta_angular_momentum - change in J at each point due to
!       rotationally induced instabilities.
!   angular_momentum - contains the run of angular momentum prior to
!       angular momentum redistribution.
!   rot_diff%ntot (common/egrid/) - number of equally spaced grid points in the
!       unstable region.
!   total_delta_angular_momentum - total change in angular momentum
!       across grid.
!
!   Variables pertaining to the original model points:
!
!   shell_mass - the unlogged mass content of the model points.
!   zone_begin, zone_end - the first and last model points contained
!       in the unstable region.
!   convective_flag - a logical array which is T if the model point is
!       convective for angular momentum purposes and F otherwise.
!   num_points - the total number of model points.
!
! Variables in common blocks:
!
! common/splin/ contains dummy arrays used to perform an osculatory
! spline interpolation.
! common/egrid/ contains the vector of equally spaced grid points.
!
! Output variable:
!   specific_angular_momentum - the new run of specific angular
!       momentum at the model points.
!-------------------------------------------------------------------------
subroutine eq2mod(delta_angular_momentum, angular_momentum, shell_mass, &
     zone_begin, zone_end, convective_flag, num_points, &
     total_delta_angular_momentum, specific_angular_momentum)
      use rotdiff_lib
      use numerics_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: delta_angular_momentum(json), &
           angular_momentum(json), shell_mass(json)
      integer, intent(in) :: zone_begin, zone_end
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_points
      double precision, intent(in) :: total_delta_angular_momentum
      double precision, intent(inout) :: specific_angular_momentum(json)


      save

      integer :: i0, i1, i, nmod, ii
      double precision :: sumjmod, sumdjmod, test, val, ratio

! CHECK FOR LOWER AND UPPER CONVECTION ZONES.
      if (.not.convective_flag(zone_begin) .or. zone_begin.eq.1) then
         i0 = zone_begin
      else
         do i = zone_begin - 1, 1, -1
            if (.not.convective_flag(i)) then
               i0 = i + 1
               goto 10
            end if
         end do
         i0 = 1
   10    continue
      end if
      if (.not.convective_flag(zone_end) .or. zone_end.eq.num_points) then
         i1 = zone_end
      else
         do i = zone_end+1, num_points
            if (.not.convective_flag(i)) then
               i1 = i - 1
               goto 20
            end if
         end do
         i1 = num_points
   20    continue
      end if
! INTERPOLATE IN DJ/J AS A FUNCTION OF ECHI
      nmod = zone_end - zone_begin + 1
      do i = 1, nmod
         rot_diff%xval(i) = rot_diff%chi(i)
      end do
      do i = 1, rot_diff%ntot
         rot_diff%xtab(i) = rot_diff%echi(i)
         rot_diff%ytab(i) = delta_angular_momentum(i)/angular_momentum(i)
      end do
      call osplin(rot_diff%xval, rot_diff%yval, rot_diff%xtab, rot_diff%ytab, rot_diff%ntot, nmod)
! APPLY THE FRACTIONAL CHANGE IN J/M ACROSS THE MODEL.
      sumjmod = 0.0d0
      sumdjmod = 0.0d0
      if (i0.lt.zone_begin) then
         do i = i0, zone_begin-1
            sumjmod = sumjmod + specific_angular_momentum(i)*shell_mass(i)
            sumdjmod = sumdjmod + rot_diff%yval(1)*specific_angular_momentum(i)* &
                 shell_mass(i)
            specific_angular_momentum(i) = specific_angular_momentum(i)* &
                 (1.0d0+rot_diff%yval(1))
         end do
      end if
      do i = zone_begin, zone_end
         ii = i - zone_begin + 1
         sumjmod = sumjmod + specific_angular_momentum(i)*shell_mass(i)
         sumdjmod = sumdjmod + rot_diff%yval(ii)*specific_angular_momentum(i)* &
              shell_mass(i)
         specific_angular_momentum(i) = specific_angular_momentum(i)* &
              (1.0d0+rot_diff%yval(ii))
      end do
      if (i1.gt.zone_end) then
         do i = zone_end+1, i1
            sumjmod = sumjmod + specific_angular_momentum(i)*shell_mass(i)
            sumdjmod = sumdjmod + rot_diff%yval(nmod)*specific_angular_momentum(i)* &
                 shell_mass(i)
            specific_angular_momentum(i) = specific_angular_momentum(i)* &
                 (1.0d0+rot_diff%yval(nmod))
         end do
      end if
! CHECK THAT ANGULAR MOMENTUM HAS BEEN CONSERVED.  THE TOTAL ANGULAR
! MOMENTUM CAN BE REDUCED FROM SURFACE ANGULAR MOMENTUM LOSS; CARE IS
! USED HERE BECAUSE DUMDJ SHOULD BE ZERO IN MODELS WITHOUT LOSS.
      test = abs(sumdjmod - total_delta_angular_momentum)
      val = 1.0d-5*abs(total_delta_angular_momentum)
!      WRITE(*,*)SUMDJ,SUMDJMOD,SUMJMOD
!      WRITE(*,911)(DJ(I)/EJ(I),I=1,NTOT)
!  911  FORMAT(1P8E10.2)
! IF THE TEST SUM IS NOT THE SAME WITHIN THE RELATIVE TOLERANCE VAL,
! ADJUST THE TOTAL ANGULAR MOMENTUM OF THE ENTIRE REGION BY A CONSTANT
! FACTOR TO ENFORCE ANGULAR MOMENTUM CONSERVATION.
      if (test.gt.val) then
         ratio = (sumjmod + total_delta_angular_momentum)/ &
              (sumjmod + sumdjmod)
         do i = i0, i1
            specific_angular_momentum(i) = specific_angular_momentum(i)*ratio
         end do
      end if
      return
end subroutine eq2mod
