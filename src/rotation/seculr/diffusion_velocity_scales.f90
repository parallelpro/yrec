!----------------------------------------------------------------------
! codiff
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original codiff.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! 2/91 SUBROUTINE ALTERED TO ALLOW DIFFERENT FC'S FOR DIFFERENT
!      MECHANISMS (COMMON BLOCK VMULT2)
!
!  THIS SR DETERMINES THE RUN OF CHARACTERISTIC VELOCITY LENGTH
!  SCALES FOR THE DIFFUSION EQUATIONS.
!
!  INPUT VARIABLES
!
!       DECLARED:
!  radius_mid_prev (HRU) : RUN OF RADIUS (UNLOGGED) OF THE MIDPOINTS OF
!     THE SHELLS.
!  num_zones (M) : NUMBER OF SHELLS.
!  FC,FO,FW : USER PARAMETERS.  THE DIFFUSION COEFFICIENTS
!  CAN BE MODIFIED BY THESE THREE PARAMETERS.
!  VELM : RUN OF CONVECTIVE VELOCITIES.
!
!  OUTPUT VARIABLES
!
!       DECLARED
!  am_diffusion_coeff (COD) : RUN OF DIFFUSION COEFFICENTS FOR ANGULAR
!     MOMENTUM TRANSPORT.
!  mixing_diffusion_coeff (COD2) : SAME, WITH RESPECT TO COMPOSITION
!     TRANSPORT.
!  radius_mid (RMID) : radius_mid(I) IS DEFINED AS THE AVERAGE OF
!     radius_mid_prev(I) AND radius_mid_prev(I-1).
!       FROM COMMON BLOCKS:
!       IN COMMON BLOCKS
!  HLE : RUN OF CHARACTERSITIC VELOCITY LENGTH SCALES.
!
!  NOTE ON STORAGE: COD,V, AND HLE ARE CALCULATED AT THE MIDPOINT
!  BETWEEN MASS POINTS. ELEMENT I CONTAINS THE INFORMATION FOR THE
!  (I,I-1) INTERFACE.
!
! Note: the original dummy argument name HRU is reused (positionally)
! both for the radius run and, historically, other quantities in
! commented-out callers -- here it is unambiguously the run of
! (unlogged) shell-midpoint radii, matching circulation_velocities.f90's radius(json).
subroutine diffusion_velocity_scales(radius_mid_prev, num_zones, radius_mid, &
     am_diffusion_coeff, mixing_diffusion_coeff)
      use rotation_scratch_lib
      use star_info_lib, only: star

      use star_info_lib
      use phys_const_lib
      implicit none
      double precision, intent(in) :: radius_mid_prev(json)
      integer, intent(in) :: num_zones
      double precision, intent(out) :: radius_mid(json)
      double precision, intent(out) :: am_diffusion_coeff(json), &
           mixing_diffusion_coeff(json)
      integer :: i
      double precision :: con1, con2

!  THIS SR DETERMINES THE RUN OF CHARACTERISTIC VELOCITY LENGTH
!  SCALES FOR THE DIFFUSION EQUATIONS.
      circ_scr%hle(1) = 0.0d0
      radius_mid(1) = 0.0d0
      do i = 2,num_zones
         circ_scr%hle(i) = 0.0d0
         radius_mid(i) = 0.5d0*(radius_mid_prev(i)+radius_mid_prev(i-1))
      end do
! MHP 9/14 ADDED LOOP TO ALLOW A CONSTANT BACKGROUND DIFFUSION COEFFICIENT
      con1 = c4pi*star%ctrl%difad_velocity_scale
      if (.not.star%ctrl%use_constant_background_diffusion) then
         if (.not.star%ctrl%lvfc) then
            con2 = con1*star%ctrl%mixing_velocity_scale
            do i = 2,num_zones
               am_diffusion_coeff(i)=con1*(star%es_circulation_velocity(i)+ &
                    star%gsf_circulation_velocity(i)+star%secular_shear_velocity(i))* &
                    radius_mid(i)
               mixing_diffusion_coeff(i)=con2*(star%ctrl%es_mixing_scale* &
                    star%es_circulation_velocity(i)+star%ctrl%gsf_mixing_scale* &
                    star%gsf_circulation_velocity(i)+star%ctrl%secular_shear_mixing_scale* &
                    star%secular_shear_velocity(i))*radius_mid(i)
!               HLE(I)=RMID(I)
            end do
         else
            do i = 2,num_zones
               am_diffusion_coeff(i) = con1*(star%es_circulation_velocity(i)+ &
                    star%gsf_circulation_velocity(i)+star%secular_shear_velocity(i))* &
                    radius_mid(i)
               mixing_diffusion_coeff(i) = am_diffusion_coeff(i)*star%vfc(i)
!               HLE(I) = RMID(I)
            end do
         end if
      else
         if (.not.star%ctrl%lvfc) then
            con2 = con1*star%ctrl%mixing_velocity_scale
            do i = 2,num_zones
               am_diffusion_coeff(i)=con1*(star%ctrl%constant_background_diffusion_coeff+ &
                    (star%es_circulation_velocity(i)+star%gsf_circulation_velocity(i)+ &
                    star%secular_shear_velocity(i))*radius_mid(i))
               mixing_diffusion_coeff(i)=con2*(star%ctrl%es_mixing_scale* &
                    star%es_circulation_velocity(i)+star%ctrl%gsf_mixing_scale* &
                    star%gsf_circulation_velocity(i)+star%ctrl%secular_shear_mixing_scale* &
                    star%secular_shear_velocity(i))*radius_mid(i)
            end do
         else
            do i = 2,num_zones
               am_diffusion_coeff(i) = con1*(star%es_circulation_velocity(i)+ &
                    star%gsf_circulation_velocity(i)+star%secular_shear_velocity(i))* &
                    radius_mid(i)
               mixing_diffusion_coeff(i) = am_diffusion_coeff(i)*star%vfc(i)
! MHP 8/13 ADD D.C. AFTER SCALE FACTOR FOR MIXING APPLIED
               am_diffusion_coeff(i) = am_diffusion_coeff(i)+ &
                    con1*star%ctrl%constant_background_diffusion_coeff
            end do
         end if
      end if
! OPTION TO SUPPRESS ANGULAR MOMENTUM TRANSPORT (BUT PERMIT MIXING).
      if (star%ctrl%no_am_transport_in_core) then
         do i = 2,num_zones
            am_diffusion_coeff(i) = 0.0d0
         end do
      end if
      return
end subroutine diffusion_velocity_scales
