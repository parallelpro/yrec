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
! (unlogged) shell-midpoint radii, matching vcirc.f90's radius(json).
subroutine codiff(radius_mid_prev, num_zones, radius_mid, &
     am_diffusion_coeff, mixing_diffusion_coeff)

      use temp2_lib
      use mdphy_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: radius_mid_prev(json)
      integer, intent(in) :: num_zones
      double precision, intent(out) :: radius_mid(json)
      double precision, intent(out) :: am_diffusion_coeff(json), &
           mixing_diffusion_coeff(json)


! common/intvar/: not used in this file. Naming matches vcirc.f90.
      double precision :: interface_luminosity(json), delami(json), &
           delmi(json), dm(json), epsilm(json), interface_gravity_factor(json), &
           hs3(json), pm(json), qdtmi(json), interface_radius(json), tm(json)
      common/intvar/ interface_luminosity, delami, delmi, dm, epsilm, &
           interface_gravity_factor, hs3, pm, qdtmi, interface_radius, tm





! common/advec/: not used in this file (no references anywhere in the
! executable code of the original). Not referenced in any already-
! converted file; kept as lowercased placeholders pending a confirmed
! source.
      double precision :: fadv(json), fadv0(json)
      common/advec/ fadv, fadv0



! MHP 02/12 PERMIT CONSTANT DIFFUSION COEFFICIENT
! KC 2025-05-30 reordered common block elements
!       COMMON/MAG/LCODM,CODM
! common/mag/: constant_background_diffusion_coeff/use_constant_
! background_diffusion (originally CODM/LCODM), both used here. Naming
! is local to this batch; not referenced in any already-converted file.
      double precision :: constant_background_diffusion_coeff
      logical :: use_constant_background_diffusion
      common/mag/ constant_background_diffusion_coeff, &
           use_constant_background_diffusion

      save

      integer :: i
      double precision :: con1, con2

!  THIS SR DETERMINES THE RUN OF CHARACTERISTIC VELOCITY LENGTH
!  SCALES FOR THE DIFFUSION EQUATIONS.
      circ_vel%hle(1) = 0.0d0
      radius_mid(1) = 0.0d0
      do i = 2,num_zones
         circ_vel%hle(i) = 0.0d0
         radius_mid(i) = 0.5d0*(radius_mid_prev(i)+radius_mid_prev(i-1))
      end do
! MHP 9/14 ADDED LOOP TO ALLOW A CONSTANT BACKGROUND DIFFUSION COEFFICIENT
      con1 = c4pi*difad_velocity_scale
      if (.not.use_constant_background_diffusion) then
         if (.not.lvfc) then
            con2 = con1*mixing_velocity_scale
            do i = 2,num_zones
               am_diffusion_coeff(i)=con1*(circ_vel%es_circulation_velocity(i)+ &
                    circ_vel%gsf_circulation_velocity(i)+circ_vel%secular_shear_velocity(i))* &
                    radius_mid(i)
               mixing_diffusion_coeff(i)=con2*(es_mixing_scale* &
                    circ_vel%es_circulation_velocity(i)+gsf_mixing_scale* &
                    circ_vel%gsf_circulation_velocity(i)+secular_shear_mixing_scale* &
                    circ_vel%secular_shear_velocity(i))*radius_mid(i)
!               HLE(I)=RMID(I)
            end do
         else
            do i = 2,num_zones
               am_diffusion_coeff(i) = con1*(circ_vel%es_circulation_velocity(i)+ &
                    circ_vel%gsf_circulation_velocity(i)+circ_vel%secular_shear_velocity(i))* &
                    radius_mid(i)
               mixing_diffusion_coeff(i) = am_diffusion_coeff(i)*vfc(i)
!               HLE(I) = RMID(I)
            end do
         end if
      else
         if (.not.lvfc) then
            con2 = con1*mixing_velocity_scale
            do i = 2,num_zones
               am_diffusion_coeff(i)=con1*(constant_background_diffusion_coeff+ &
                    (circ_vel%es_circulation_velocity(i)+circ_vel%gsf_circulation_velocity(i)+ &
                    circ_vel%secular_shear_velocity(i))*radius_mid(i))
               mixing_diffusion_coeff(i)=con2*(es_mixing_scale* &
                    circ_vel%es_circulation_velocity(i)+gsf_mixing_scale* &
                    circ_vel%gsf_circulation_velocity(i)+secular_shear_mixing_scale* &
                    circ_vel%secular_shear_velocity(i))*radius_mid(i)
            end do
         else
            do i = 2,num_zones
               am_diffusion_coeff(i) = con1*(circ_vel%es_circulation_velocity(i)+ &
                    circ_vel%gsf_circulation_velocity(i)+circ_vel%secular_shear_velocity(i))* &
                    radius_mid(i)
               mixing_diffusion_coeff(i) = am_diffusion_coeff(i)*vfc(i)
! MHP 8/13 ADD D.C. AFTER SCALE FACTOR FOR MIXING APPLIED
               am_diffusion_coeff(i) = am_diffusion_coeff(i)+ &
                    con1*constant_background_diffusion_coeff
            end do
         end if
      end if
! OPTION TO SUPPRESS ANGULAR MOMENTUM TRANSPORT (BUT PERMIT MIXING).
      if (no_am_transport_in_core) then
         do i = 2,num_zones
            am_diffusion_coeff(i) = 0.0d0
         end do
      end if
      return
end subroutine codiff
