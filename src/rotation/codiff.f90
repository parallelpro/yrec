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

!       old version: COMMON/INTVR2/AMUM(JSON),THDIFM(JSON),VISCM(JSON),WM(JSON)
! common/mdphy/: not used in this file. Naming matches vcirc.f90/
! rotgrid.f90.
      double precision :: amum(json), cpm(json), delm(json), &
           del_adiabatic_mix(json), del_radiative_mix(json), esumm(json), &
           om(json), qdtm(json), thdifm(json), velm(json), viscm(json), &
           epsm(json)
      common/mdphy/ amum, cpm, delm, del_adiabatic_mix, del_radiative_mix, &
           esumm, om, qdtm, thdifm, velm, viscm, epsm

! common/temp2/: es_circulation_velocity/gsf_circulation_velocity/
! secular_shear_velocity/hle are used here; the "_prev" pair and
! mu_gradient_velocity are unused placeholders. Naming matches
! vcirc.f90.
      double precision :: es_circulation_velocity(json), &
           es_circulation_velocity_prev(json), secular_shear_velocity(json), &
           secular_shear_velocity_prev(json), hle(json), &
           gsf_circulation_velocity(json), gsf_circulation_velocity_prev(json), &
           mu_gradient_velocity(json)
      common/temp2/ es_circulation_velocity, es_circulation_velocity_prev, &
           secular_shear_velocity, secular_shear_velocity_prev, hle, &
           gsf_circulation_velocity, gsf_circulation_velocity_prev, &
           mu_gradient_velocity

! common/vmult/: only difad_velocity_scale/mixing_velocity_scale (FW/FC)
! are used here. Naming matches rotgrid.f90/vcirc.f90/dadcoeft.f90.
      double precision :: difad_velocity_scale, mixing_velocity_scale, fo, &
           es_velocity_scale, gsf_velocity_scale, mu_gradient_scale, &
           secular_shear_velocity_scale, critical_reynolds
      common/vmult/ difad_velocity_scale, mixing_velocity_scale, fo, &
           es_velocity_scale, gsf_velocity_scale, mu_gradient_scale, &
           secular_shear_velocity_scale, critical_reynolds

! common/vmult2/: only es_mixing_scale/secular_shear_mixing_scale/
! gsf_mixing_scale (FESC/FSSC/FGSFC) are used here; ies/gsf_inhibition_
! mode/imu are unused placeholders. Naming matches vcirc.f90/
! dadcoeft.f90.
      double precision :: es_mixing_scale, secular_shear_mixing_scale, &
           gsf_mixing_scale
      integer :: ies, gsf_inhibition_mode, imu
      common/vmult2/ es_mixing_scale, secular_shear_mixing_scale, &
           gsf_mixing_scale, ies, gsf_inhibition_mode, imu

! common/advec/: not used in this file (no references anywhere in the
! executable code of the original). Not referenced in any already-
! converted file; kept as lowercased placeholders pending a confirmed
! source.
      double precision :: fadv(json), fadv0(json)
      common/advec/ fadv, fadv0

! MHP 7/93
! common/varfc/: only vfc/lvfc are used here (under the LVFC branch);
! use_diffusion_advection_transport (LDIFAD) is an unused placeholder
! in this file. Naming matches rotgrid.f90/vcirc.f90.
      double precision :: vfc(json)
      logical :: lvfc, use_diffusion_advection_transport
      common/varfc/ vfc, lvfc, use_diffusion_advection_transport

! MHP 9/93
! common/notran/: no_am_transport_in_core (originally LNOJ) is used
! here to optionally suppress angular-momentum transport. Naming
! matches vcirc.f90.
      logical :: no_am_transport_in_core
      common/notran/ no_am_transport_in_core

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
      hle(1) = 0.0d0
      radius_mid(1) = 0.0d0
      do i = 2,num_zones
         hle(i) = 0.0d0
         radius_mid(i) = 0.5d0*(radius_mid_prev(i)+radius_mid_prev(i-1))
      end do
! MHP 9/14 ADDED LOOP TO ALLOW A CONSTANT BACKGROUND DIFFUSION COEFFICIENT
      con1 = c4pi*difad_velocity_scale
      if (.not.use_constant_background_diffusion) then
         if (.not.lvfc) then
            con2 = con1*mixing_velocity_scale
            do i = 2,num_zones
               am_diffusion_coeff(i)=con1*(es_circulation_velocity(i)+ &
                    gsf_circulation_velocity(i)+secular_shear_velocity(i))* &
                    radius_mid(i)
               mixing_diffusion_coeff(i)=con2*(es_mixing_scale* &
                    es_circulation_velocity(i)+gsf_mixing_scale* &
                    gsf_circulation_velocity(i)+secular_shear_mixing_scale* &
                    secular_shear_velocity(i))*radius_mid(i)
!               HLE(I)=RMID(I)
            end do
         else
            do i = 2,num_zones
               am_diffusion_coeff(i) = con1*(es_circulation_velocity(i)+ &
                    gsf_circulation_velocity(i)+secular_shear_velocity(i))* &
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
                    (es_circulation_velocity(i)+gsf_circulation_velocity(i)+ &
                    secular_shear_velocity(i))*radius_mid(i))
               mixing_diffusion_coeff(i)=con2*(es_mixing_scale* &
                    es_circulation_velocity(i)+gsf_mixing_scale* &
                    gsf_circulation_velocity(i)+secular_shear_mixing_scale* &
                    secular_shear_velocity(i))*radius_mid(i)
            end do
         else
            do i = 2,num_zones
               am_diffusion_coeff(i) = con1*(es_circulation_velocity(i)+ &
                    gsf_circulation_velocity(i)+secular_shear_velocity(i))* &
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
