!----------------------------------------------------------------------
! vcirc
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original vcirc.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
!   SKIPPING CONVECTIVE REGIONS, CHECK FOR INSTABILITY AGAINST EACH
!   SECULAR MECHANISM.IF UNSTABLE,COMPUTE A CIRCULATION VELOCITY.
!   THE PROCEDURE FOR CALCULATING STABILITY IS AS FOLLOWS:
!      STABILITY IS CALCULATED FOR THE I/I-1 INTERFACE.
!      EACH STABILITY CRITERION IS EXPRESSED AS FOLLOWS:
!       A COMBINATION OF VARIABLES(RHO,P,ETC.)STORED AT EACH SHELL
!       MULTIPLIED BY A GRADIENT BETWEEN THE SHELLS MUST BE LESS THAN
!       SOME CRITICAL NUMBER.
! Computes, at each unstable interface, the Eddington-Sweet meridional
! circulation velocity (circ_vel%es_circulation_velocity, Zahn 1991/1992, with
! a quadrupole correction), the GSF-instability circulation velocity
! (circ_vel%gsf_circulation_velocity, Kippenhahn 1980 estimate), and the
! diffusive/dynamical shear velocity (circ_vel%secular_shear_velocity), then
! combines them into total_circulation_velocity (HV).
!
!       SUBROUTINE VCIRC(HJM,HR,HRU,IMIN,IMAX,IT,LCZ,M,OMEGA,LDO,  ! KC 2025-05-31
!                        QWLNR,QWRMAX,HV,CLTOT,DT,HP)
subroutine vcirc(log_radius, radius, zone_min, zone_max, iteration, &
     am_transport_convective_flag, num_zones, omega, any_transport_active, &
     dlnomega_dlnr, dlnomega_dlnr_max, total_circulation_velocity, &
     total_luminosity, timestep, log_pressure)

      use temp2_lib
      use mdphy_lib
      use luout_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: log_radius(json), radius(json)
      integer, intent(in) :: zone_min, zone_max, iteration
      logical, intent(in) :: am_transport_convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(in) :: omega(json)
      logical, intent(out) :: any_transport_active
      double precision, intent(in) :: dlnomega_dlnr(json), &
           dlnomega_dlnr_max(json)
      double precision, intent(out) :: total_circulation_velocity(json)
      double precision, intent(in) :: total_luminosity, timestep, &
           log_pressure(json)





! common/intfac/: lagrange_interp_weights (originally FACI), the 4-point
! Lagrangian interpolation weights used to interpolate viscosity/
! thermal diffusivity onto zone interfaces. Naming is local to this
! batch.
      double precision :: lagrange_interp_weights(4,json)
      common/intfac/ lagrange_interp_weights

! MHP 3/92 COMMON BLOCKS ADDED FOR SHEAR VELOCITY AND MU INHIBITION.
! common/taukh/: es_velocity_coeff1/es_velocity_coeff2/es_shear_coeff
! (FES1/FES2/FES3, leading/next-order/shear-driven-term coefficients of
! the Eddington-Sweet velocity expansion) and gsf_kippenhahn_coeff
! (FGSFK, the Kippenhahn 1980 GSF velocity coefficient) are used here;
! fact6/fgsfj are unused placeholders (their only references are in
! commented-out code -- the KM1974 and James & Kahn 1971 velocity
! estimates, both superseded per the comments below). Naming is local
! to this batch.
      double precision :: fact6(json), es_velocity_coeff1(json), &
           es_velocity_coeff2(json), fgsfj(json), gsf_kippenhahn_coeff(json), &
           es_shear_coeff(json)
      common/taukh/ fact6, es_velocity_coeff1, es_velocity_coeff2, fgsfj, &
           gsf_kippenhahn_coeff, es_shear_coeff

! common/intvar/: interface (interpolated) counterparts of various
! model quantities. Only interface_luminosity (ALM, used to test
! whether a zone is within the "no transport" core in the LNOJ branch)
! and interface_gravity_factor (HGM, an interface gravity-like factor
! used in the GSF mu-gradient stability check) and interface_radius
! (RM) are used here; the rest are unused placeholders (their only
! references, if any, are in commented-out code). Naming is local to
! this batch.
      double precision :: interface_luminosity(json), delami(json), &
           delmi(json), dm(json), epsilm(json), interface_gravity_factor(json), &
           hs3(json), pm(json), qdtmi(json), interface_radius(json), tm(json)
      common/intvar/ interface_luminosity, delami, delmi, dm, epsilm, &
           interface_gravity_factor, hs3, pm, qdtmi, interface_radius, tm

!      JvS 09/25 RENAMED ARRAYS IN INTVR2 TO PLAY NICELY WITH MDPHY
!      COMMON/INTVR2/AMUM(JSON),THDIFM(JSON),VISCM(JSON),WM(JSON)
! common/intvr2/: mean_molecular_weight_interface/
! thermal_diffusivity_interface/kinematic_viscosity_interface/
! omega_interface (originally AMUMI/THDIFMI/VISCMI/WM), all used here
! -- interface-averaged counterparts of mix_phys%amum/mix_phys%thdifm/
! mix_phys%viscm (former common/mdphy/) and of the omega dummy
! argument. Naming is local to this batch.
      double precision :: mean_molecular_weight_interface(json), &
           thermal_diffusivity_interface(json), &
           kinematic_viscosity_interface(json), omega_interface(json)
      common/intvr2/ mean_molecular_weight_interface, &
           thermal_diffusivity_interface, kinematic_viscosity_interface, &
           omega_interface




! MHP 06/02 ADDED FACT7 AND FACT8 FOR DIF/AD TREATMENT
! common/vfact/: mu_gradient_richardson_coeff (FACT5, used in the GSF
! mu-gradient Richardson-number stability check) and
! difad_shear_coeff1/difad_shear_coeff2 (FACT7/FACT8, MHP 06/02, added
! for the diffusion+advection treatment) are used here; fact1-fact4 are
! unused placeholders (their only references are in commented-out
! code). Naming is local to this batch.
      double precision :: fact1(json), fact2(json), fact3(json), fact4(json), &
           mu_gradient_richardson_coeff(json), difad_shear_coeff1(json), &
           difad_shear_coeff2(json)
      common/vfact/ fact1, fact2, fact3, fact4, mu_gradient_richardson_coeff, &
           difad_shear_coeff1, difad_shear_coeff2

! common/vmult/: only es_velocity_scale/gsf_velocity_scale/
! secular_shear_velocity_scale/critical_reynolds are used here. Naming
! matches rotgrid.f90.
      double precision :: difad_velocity_scale, mixing_velocity_scale, fo, &
           es_velocity_scale, gsf_velocity_scale, mu_gradient_scale, &
           secular_shear_velocity_scale, critical_reynolds
      common/vmult/ difad_velocity_scale, mixing_velocity_scale, fo, &
           es_velocity_scale, gsf_velocity_scale, mu_gradient_scale, &
           secular_shear_velocity_scale, critical_reynolds

! common/vmult2/: only gsf_inhibition_mode (IGSF) is used here; the
! rest are unused placeholders (their only references are in
! commented-out code). Naming matches dadcoeft.f90.
      double precision :: es_mixing_scale, secular_shear_mixing_scale, &
           gsf_mixing_scale
      integer :: ies, gsf_inhibition_mode, imu
      common/vmult2/ es_mixing_scale, secular_shear_mixing_scale, &
           gsf_mixing_scale, ies, gsf_inhibition_mode, imu

! common/prevmu/: mu_gradient_velocity_prev (originally VMUP), the
! previous-iteration circ_vel%mu_gradient_velocity, saved (but not subsequently
! read back) here. Naming is local to this batch.
      double precision :: mu_gradient_velocity_prev(json)
      common/prevmu/ mu_gradient_velocity_prev

! 3/93 QUADRUPOLE DRIVING TERM ADDED.
! common/quadd/: only circulation_correction_ratio (RAT) is used here;
! phisp/phirot/phidis are unused placeholders. Naming is local to this
! batch.
      double precision :: phisp(json), phirot(json), phidis(json), &
           circulation_correction_ratio(json)
      common/quadd/ phisp, phirot, phidis, circulation_correction_ratio

      double precision :: quadrupole_moment(json), local_gravity(json)
      common/quadru/ quadrupole_moment, local_gravity

! MHP 9/93
! common/notran/: no_am_transport_in_core (originally LNOJ), "9/93
! MIXING WITHOUT TRANSPORT ADDED" -- zeroes out the AM-transport
! velocity/diffusion estimates in the core to avoid numerical problems
! in the H-burning shell.
      logical :: no_am_transport_in_core
      common/notran/ no_am_transport_in_core

! common/varfc/: only use_diffusion_advection_transport (LDIFAD) is
! used here; vfc/lvfc are unused placeholders. Naming matches
! rotgrid.f90.
      double precision :: vfc(json)
      logical :: lvfc, use_diffusion_advection_transport
      common/varfc/ vfc, lvfc, use_diffusion_advection_transport

! MHP 11/94 ADVECTIVE VELOCITY AND DIFFUSIVE COMPONENT.
! common/difad2/: es_advective_velocity/es_advective_velocity_prev/
! es_diffusive_velocity/es_diffusive_velocity_prev (VESA/VESA0/VESD/
! VESD0), all used here. Naming matches rotgrid.f90.
      double precision :: es_advective_velocity(json), &
           es_advective_velocity_prev(json), es_diffusive_velocity(json), &
           es_diffusive_velocity_prev(json)
      common/difad2/ es_advective_velocity, es_advective_velocity_prev, &
           es_diffusive_velocity, es_diffusive_velocity_prev

! MHP 06/02
! HIGHER ORDER TERMS IN M.C. VELOCITY
! common/difad3/: see dadcoeft.f90 for the full member-by-member
! discussion; here facd2/facd3 (used to build vesd2/vesd3, the second/
! third-order velocity terms), fgeom (the geometric/Jacobian factor,
! used to scale advective_term1/diffusive_term2 into am_advective_coeff/
! am_diffusive_coeff -- note this file does not itself declare
! common/difad/, so those two are purely local outputs of the FW*EV0*..
! formula, not shared with rotgrid.f90/dadcoeft.f90's ECOD3/ECOD4),
! and eq_velocity_coeff0/eq_velocity_coeff2a/eq_velocity_coeff2b are
! used; fv1a/fv1b/fv2a/fv2b are unused placeholders. Naming matches
! dadcoeft.f90.
      double precision :: facd2(json), facd3(json), vesd2(json), &
           vesd3(json), am_2nd_deriv_coeff(json), am_3rd_deriv_coeff(json), &
           geometric_factor(json), velocity_coeff0(json), &
           velocity_coeff1a(json), velocity_coeff1b(json), &
           velocity_coeff2a(json), velocity_coeff2b(json), &
           eq_velocity_coeff0(json), eq_velocity_coeff1a(json), &
           eq_velocity_coeff1b(json), eq_velocity_coeff2a(json), &
           eq_velocity_coeff2b(json), shear_diffusion_coeff(json), &
           gsf_diffusion_coeff(json), shear_diffusion_coeff_eqgrid(json), &
           gsf_diffusion_coeff_eqgrid(json)
      common/difad3/ facd2, facd3, vesd2, vesd3, am_2nd_deriv_coeff, &
           am_3rd_deriv_coeff, geometric_factor, velocity_coeff0, &
           velocity_coeff1a, velocity_coeff1b, velocity_coeff2a, &
           velocity_coeff2b, eq_velocity_coeff0, eq_velocity_coeff1a, &
           eq_velocity_coeff1b, eq_velocity_coeff2a, eq_velocity_coeff2b, &
           shear_diffusion_coeff, gsf_diffusion_coeff, &
           shear_diffusion_coeff_eqgrid, gsf_diffusion_coeff_eqgrid

! Time change of theta
! common/oldrot2/: theta_new/theta_mean/theta_prev (THN/THM/THP),
! del_grad_diff_interface (DDELM), and es_relaxation_factor (FESTIME)
! are used here; tho/qwrst/wmst/qwrmst are unused placeholders. The
! precise physical meaning of the "theta" quantities beyond the
! algebra below is not confidently identified. Naming is local to this
! batch.
      double precision :: tho(json), theta_new(json), theta_mean(json), &
           del_grad_diff_interface(json), es_relaxation_factor(json), &
           theta_prev(json), qwrst(json), wmst(json), qwrmst(json)
      common/oldrot2/ tho, theta_new, theta_mean, del_grad_diff_interface, &
           es_relaxation_factor, theta_prev, qwrst, wmst, qwrmst

! DEFINITION TERMS FOR THE SECOND AND THIRD DERIVATIVE
! TERMS
! common/egridder/: only second_deriv_geom_factor_eqgrid/
! second_deriv_geom_factor (EQQCOD/QQCOD, used only under
! use_diffusion_advection_transport) are used here. Naming matches
! dadcoeft.f90.
      double precision :: second_deriv_geom_factor_eqgrid(json), &
           third_deriv_geom_factor_eqgrid(json), second_deriv_geom_factor(json), &
           third_deriv_geom_factor(json)
      common/egridder/ second_deriv_geom_factor_eqgrid, &
           third_deriv_geom_factor_eqgrid, second_deriv_geom_factor, &
           third_deriv_geom_factor

! SPEED OF LIGHT.
!       DATA CC/2.99792458D10/
      save

      integer :: i, i0, i1
      double precision :: dr, gmid, phistd, phis2, phi2, qw, &
           q1, q2, q0, qmu, qp, ddel, ddtest, ddtest2, qqq, v2, &
           qwrmx, fxx, qwrmx2, qwr, rmid, dlnwdr, dlnjmdr, fx, qwrmx_dyn, &
           wmin, wmid, dlnwdr0, dlnjmdr0, dell

! SKIPPING CONVECTIVE REGIONS, CHECK FOR INSTABILITY AGAINST EACH
! SECULAR MECHANISM.IF UNSTABLE,COMPUTE A CIRCULATION VELOCITY.
! THE PROCEDURE FOR CALCULATING STABILITY IS AS FOLLOWS:
!    STABILITY IS CALCULATED FOR THE I/I-1 INTERFACE.
!    EACH STABILITY CRITERION IS EXPRESSED AS FOLLOWS:
!     A COMBINATION OF VARIABLES(RHO,P,ETC.)STORED AT EACH SHELL
!     MULTIPLIED BY A GRADIENT BETWEEN THE SHELLS MUST BE LESS THAN
!     SOME CRITICAL NUMBER.

      any_transport_active = .false.
!  STORE THE OLD VELOCITY ESTIMATES FOR LATER USE.
      if (iteration.gt.1) then
         do 5 i = zone_min,zone_max
            mu_gradient_velocity_prev(i) = circ_vel%mu_gradient_velocity(i)
            circ_vel%gsf_circulation_velocity_prev(i) = circ_vel%gsf_circulation_velocity(i)
            circ_vel%es_circulation_velocity_prev(i) = circ_vel%es_circulation_velocity(i)
            circ_vel%secular_shear_velocity_prev(i) = circ_vel%secular_shear_velocity(i)
    5    continue
         if (use_diffusion_advection_transport) then
            do i = zone_min,zone_max
               es_advective_velocity_prev(i) = es_advective_velocity(i)
               es_diffusive_velocity_prev(i) = es_diffusive_velocity(i)
            end do
         end if
      end if
      do 10 i = 1,num_zones
         total_circulation_velocity(i) = 0.0d0
         circ_vel%es_circulation_velocity(i) = 0.0d0
         circ_vel%gsf_circulation_velocity(i) = 0.0d0
         circ_vel%secular_shear_velocity(i) = 0.0d0
         circ_vel%mu_gradient_velocity(i) = 0.0d0
   10 continue
!  MEAN MOLECULAR WEIGHT (AS WELL AS QUANTITIES WHICH DEPEND ON IT) AND
!  THE ANGULAR VELOCITY DISTRIBUTION CHANGE DURING A DIFFUSION TIMESTEP.
!  FIND THE NEW VALUES AT THE MIDPOINT IN RADIUS BETWEEN MASS SHELLS
!  WHERE STABILITY IS EVALUATED.  4 POINT LAGRANGIAN INTERPOLATION IS USED
!  FOR THE KINEMATIC VISCOSITY AND THERMAL DIFFUSIVITY.
!  TREATMENT OF FIRST INTERFACE.
      if (zone_min.lt.3) then
         kinematic_viscosity_interface(2)=exp(log(mix_phys%viscm(1))* &
              lagrange_interp_weights(1,2)+log(mix_phys%viscm(2))* &
              lagrange_interp_weights(2,2)+ &
              log(mix_phys%viscm(3))*lagrange_interp_weights(3,2)+log(mix_phys%viscm(4))* &
              lagrange_interp_weights(4,2))
         thermal_diffusivity_interface(2)=exp(log(mix_phys%thdifm(1))* &
              lagrange_interp_weights(1,2)+log(mix_phys%thdifm(2))* &
              lagrange_interp_weights(2,2) &
              +log(mix_phys%thdifm(3))*lagrange_interp_weights(3,2)+log(mix_phys%thdifm(4))* &
              lagrange_interp_weights(4,2))
         i0=3
      else
         i0=zone_min
      end if
!  TREATMENT OF LAST INTERFACE.
      if (zone_max.eq.num_zones) then
         kinematic_viscosity_interface(num_zones)= &
              exp(log(mix_phys%viscm(num_zones-3))*lagrange_interp_weights(1,num_zones)+ &
              log(mix_phys%viscm(num_zones-2))*lagrange_interp_weights(2,num_zones) &
              +log(mix_phys%viscm(num_zones-1))*lagrange_interp_weights(3,num_zones)+ &
              log(mix_phys%viscm(num_zones))*lagrange_interp_weights(4,num_zones))
         thermal_diffusivity_interface(num_zones)= &
              exp(log(mix_phys%thdifm(num_zones-3))* &
              lagrange_interp_weights(1,num_zones)+log(mix_phys%thdifm(num_zones-2))* &
              lagrange_interp_weights(2,num_zones)+log(mix_phys%thdifm(num_zones-1))* &
              lagrange_interp_weights(3,num_zones)+ &
              log(mix_phys%thdifm(num_zones))*lagrange_interp_weights(4,num_zones))
         i1=num_zones-1
      else
         i1=zone_max
      end if
!  GENERAL CASE.
      do 20 i = i0,i1
         kinematic_viscosity_interface(i)=exp(log(mix_phys%viscm(i-2))* &
              lagrange_interp_weights(1,i)+log(mix_phys%viscm(i-1)) &
              *lagrange_interp_weights(2,i)+log(mix_phys%viscm(i))* &
              lagrange_interp_weights(3,i)+log(mix_phys%viscm(i+1)) &
              *lagrange_interp_weights(4,i))
         thermal_diffusivity_interface(i)=exp(log(mix_phys%thdifm(i-2))* &
              lagrange_interp_weights(1,i)+log(mix_phys%thdifm(i-1)) &
              *lagrange_interp_weights(2,i)+log(mix_phys%thdifm(i))* &
              lagrange_interp_weights(3,i)+log(mix_phys%thdifm(i+1)) &
              *lagrange_interp_weights(4,i))
   20 continue
! USE LINEAR INTERPOLATION FOR OMEGA AND MU.
      do i = 2,num_zones
         mean_molecular_weight_interface(i) = 0.5d0*(mix_phys%amum(i)+mix_phys%amum(i-1))
         omega_interface(i) = 0.5d0*(omega(i)+omega(i-1))
      end do
! MHP 8/03 OMITTED OLD KM1974 MERIDIONAL CIRCULATION VELOCITY ESTIMATE.
! THE IES FLAG IS THEREFORE NO LONGER IN USE.
!
! KIPPENHAHN AND MOLLENHOF(1974) MERIDIONAL CIRCULATION VELOCITY.
!      IF(IES.EQ.2)THEN
!         DO 30 I = IMIN,IMAX
!  SKIP CONVECTIVE INTERFACES.
!            IF(LCZ(I).AND.LCZ(I-1))GOTO 30
!  EDDINGTON CIRCULATION.  VELOCITY DEFINED IN SR SETUPV.
!            DR = HRU(I)-HRU(I-1)
!            DV = ((OMEGA(I)*HRU(I))**2-(OMEGA(I-1)*HRU(I-1))**2)/DR
! QUADRUPOLE TERM, AS PER ZAHN (1992), ADDED.
!            GMID = 0.5D0*(GG(I)+GG(I-1))
! THE VARIATION IN GRAVITY ON A LEVEL SURFACE CAN BE EXPRESSED IN
! GENERAL AS 1/3 OMEGA**2 *D/DR((R**2 - PHI)/G), WHERE G IS THE
! LOCAL AVERAGE GRAVITY (~GM/R**2) AND PHI IS THE QUADRUPOLE MOMENT.
! THE CLASSICAL EXPRESSION NEGLECTS THE QUADRUPOLE AND THE SPATIAL
! DERIVATIVE OF THE MASS (PHISTD).  THE VELOCITY ESTIMATE FOR
! MERIDIONAL CIRCULATION IS THEREFORE MULTIPLED BY THE RATIO OF
! THE CORRECTED EXPRESSION TO THE CLASSICAL ONE (RAT).
!            PHISTD = 2.0D0*CC23*WM(I)**2*RM(I)/GMID
!            PHIS2 = CC13*WM(I)**2*(HRU(I)**2/GG(I)-
!     *              HRU(I-1)**2/GG(I-1))/DR
!            PHI2 = (QUAD(I-1)/GG(I-1)-QUAD(I)/GG(I))/DR
!            RAT(I) = (PHIS2+PHI2)/PHISTD
!            RAT(I) = PHIS2/PHISTD
!            VES(I) = ABS(FACT1(I)*(FACT2(I)*RAT(I)*WM(I)**2
!     *               +FACT6(I)*DV))
!   30    CONTINUE
!      ELSE
! ZAHN (1991) VELOCITY ESTIMATE.
      do 31 i = zone_min,zone_max
         if (am_transport_convective_flag(i).and. &
              am_transport_convective_flag(i-1)) goto 31
! ORIGINAL ESTIMATE,USING DG/G = W**2 R**3 / GM.
         circ_vel%es_circulation_velocity(i) = omega_interface(i)**2* &
              (es_velocity_coeff1(i)+omega_interface(i)**2*es_velocity_coeff2(i))
! QUADRUPOLE TERM ADDED, DG/G COMPUTED AS PER ZAHN 1992.
         dr = radius(i) - radius(i-1)
         gmid = 0.5d0*(local_gravity(i)+local_gravity(i-1))
         phistd = 2.0d0*cc23*omega_interface(i)**2*interface_radius(i)/gmid
         phis2 = cc13*omega_interface(i)**2*(radius(i)**2/local_gravity(i)- &
                 radius(i-1)**2/local_gravity(i-1))/dr
         phi2 = (quadrupole_moment(i-1)/local_gravity(i-1)- &
              quadrupole_moment(i)/local_gravity(i))/dr
         qw = omega_interface(i)*(omega(i)-omega(i-1))/dr
         circulation_correction_ratio(i) = (phis2+phi2)/phistd
         if (circulation_correction_ratio(i).lt.0.0d0) then
            write(*,303) i,circulation_correction_ratio(i),phistd,phis2,phi2, &
            quadrupole_moment(i),quadrupole_moment(i-1),local_gravity(i), &
            local_gravity(i-1),radius(i), &
            radius(i-1)
 303        format(i5,' RAT ',1pe10.3,' POT SPH,CYL,QUA', &
                 3e10.3/' QUA,G,R ',6e12.3)
         end if
         if (.not.use_diffusion_advection_transport) then
            circ_vel%es_circulation_velocity(i) = &
                 abs(circulation_correction_ratio(i)*circ_vel%es_circulation_velocity(i)- &
                 es_shear_coeff(i)*qw)
         else
! MHP 05/02 ADD FACTOR OF 1/5 HERE
!               VESA(I) = RAT(I)*VES(I)
!               VESD(I) = ABS(FES3(I)*WM(I)**2)
! MHP 06/02 ADDED D THETA/DT TERM
            theta_new(i) = es_relaxation_factor(i)*(theta_mean(i)*qw- &
                 theta_prev(i))/timestep
            es_advective_velocity(i) = 0.2d0* &
                 (circulation_correction_ratio(i)*circ_vel%es_circulation_velocity(i)+ &
                 theta_new(i))
            q1 = difad_shear_coeff1(i)*omega_interface(i)**2
            q2 = difad_shear_coeff2(i)*qw
            q0 = es_shear_coeff(i)*omega_interface(i)**2
!               VESD(I) = 0.2D0*ABS(FES3(I)*WM(I)**2)
!               VESD(I) = 0.2D0*ABS(Q0+Q1+Q2)
            es_diffusive_velocity(i) = 0.2d0*(q0+q1+q2)
! SECOND ORDER TERM
            vesd2(i) = 0.2d0*facd2(i)*omega_interface(i)
! THIRD ORDER TERM
            vesd3(i) = 0.2d0*facd3(i)*omega_interface(i)
!               VES(I) = ABS(RAT(I)*VES(I)+THN(I)-FES3(I)*QW)
            circ_vel%es_circulation_velocity(i) = &
                 circulation_correction_ratio(i)*circ_vel%es_circulation_velocity(i)+ &
                 theta_new(i)+ &
                 (es_shear_coeff(i)+difad_shear_coeff1(i))*qw+ &
                 difad_shear_coeff2(i)*((omega(i)-omega(i-1))/dr)**2
            second_deriv_geom_factor_eqgrid(i) = &
                 second_deriv_geom_factor(i)*omega_interface(i)
         end if
   31 continue
!      ENDIF
! INHIBITION CAUSED BY GRADIENTS IN MEAN MOLECULAR WEIGHT.
! INCLUDE ONLY THE ZAHN & MAEDER 1998 TREATMENT, WHERE (DEL - DEL AD) IS
! REPLACED BY (DEL + DEL MU - DEL AD).  RETAIN IMU FLAG FOR LATER USE.
!      IF(IMU.EQ.3)THEN
!  FICTITIOUS MU CURRENTS THAT OPPOSE ES CIRCULATION CALCULATED HERE.
!         DO 32 I = IMIN,IMAX
!            IF(LCZ(I).AND.LCZ(I-1))GOTO 32
! SQUARE ROOT OF TKH*,KIPPENHAHN,IAU#66,P.23,USING EQ.(12)FOR D ON P.25.
! NOTE FACTOR OF G IS SUBSUMED IN FACT4.
!            TKHS = SQRT(VISCMI(I)*CC/THDIFMI(I))
!            QMU = (AMUM(I)-AMUM(I-1))/(HRU(I)-HRU(I-1))
!            VMU(I) = ABS(FMU*FACT1(I)*FACT4(I)*AMUMI(I)*QMU/TKHS)
!            VMU2(I)=FMU*HGM(I)*ABS(QMU)/AMUMI(I)/WM(I)**2
!   32    CONTINUE
! LOCAL TIMESCALE ESTIMATE FOR MU INHIBITION.
! ASSUMES V/L = OMEGA.
!      ELSE IF(IMU.EQ.2)THEN
! MHP 06/02 REPLACE WITH THE ZAHN&MAEDER 1998 PRESCRIPTION
       do i = zone_min,zone_max
          qmu = log(mix_phys%amum(i))-log(mix_phys%amum(i-1))
          qp = ln10*(log_pressure(i)-log_pressure(i-1))
          ddel=del_grad_diff_interface(i)+ qmu/qp
          ddtest = max(del_grad_diff_interface(i),1.0d-3)
          ddtest2 = max(ddel,1.0d-3)
          qqq = ddtest/ddtest2
          circ_vel%es_circulation_velocity(i) = circ_vel%es_circulation_velocity(i)*qqq
          if (use_diffusion_advection_transport) then
             es_advective_velocity(i) = es_advective_velocity(i)*qqq
             es_diffusive_velocity(i) = es_diffusive_velocity(i)*qqq
             vesd2(i) = vesd2(i)*qqq
             vesd3(i) = vesd3(i)*qqq
          end if
          q1 = difad_shear_coeff1(i)*omega_interface(i)**2*qqq
          q0 = es_shear_coeff(i)*omega_interface(i)**2*qqq
          dr = radius(i) - radius(i-1)
          v2 = 0.2d0*(q0+q1)*(omega(i)-omega(i-1))/dr/omega_interface(i)
! ADD MU GRADIENTS TO VELOCITY ESTIMATES
          es_relaxation_factor(i) = es_relaxation_factor(i)*qqq
          velocity_coeff0(i) = velocity_coeff0(i)*qqq
!         WRITE(*,911)I,WM(I),VESA(I),V2,Q0,Q1,THN(I),VES(I),QQQ
! 911           FORMAT(I5,1P8E12.3)
       end do
!         DO 33 I = IMIN,IMAX
!            IF(LCZ(I).AND.LCZ(I-1))GOTO 33
!            DR = HRU(I) - HRU(I-1)
!            VMU(I)=FMU*HGM(I)*ABS(AMUM(I)-AMUM(I-1))/DR/AMUMI(I)
!     *             /WM(I)**2
!   33    CONTINUE
! ALTERNATE EXPRESSION : D = D0/(1+R*DEL MU/MU)
! THIS IS ACTUALLY SIMILAR TO THE ZM98 FORMULA, THEREFORE OBSOLETE.
! SHOULD HAVE PUBLISHED IT - OH WELL.
!      ELSE
!         DO 34 I = IMIN,IMAX
!            IF(LCZ(I).AND.LCZ(I-1))GOTO 34
!            DR = HRU(I) - HRU(I-1)
!            DDEL = DELAMI(I)/MAX(1.0D-6,DELAMI(I)-DELMI(I))
!            VMU(I)=FMU*RM(I)*DDEL*ABS(AMUM(I)-AMUM(I-1))/DR/AMUMI(I)
! 34         CONTINUE
!      ENDIF
!  GSF INSTABILITY.  VELOCITY DEFINED IN SR SETUPV.
! OMIT ALL BUT THE KIPPENHAHN 1980 ESTIMATE.
! IGSF IS RETAINED, AND THE VALUE DETERMINES WHAT IS PERMITTED TO
! INHIBIT IT.
! IGSF = 0: INHIBITED BY VISCOSITY AND MU GRADIENTS
! IGSF = 1: NEITHER VISCOSITY NOR MU GRADIENTS INHIBIT
! IGSF = 2: INHIBITED BY VISCOSITY ONLY
! NOTE THAT THE ABCD INSTABILITY - RELATED TO THE GSF INSTABILITY -
! WOULD NOT BE INHIBITED BY VISCOSITY/MU GRADIENTS, SO INCLUDING IT
! IS EQUIVALENT TO SETTING IGSF=1 (THE ORIGINAL K1980 ESTIMATE).
! THE INHIBITION FACTORS COME FROM SUBSEQUENT WORK BY SPRUIT.
!      IF(IGSF.EQ.0 .OR. IGSF.EQ.2 .OR. IGSF.EQ.1)THEN
! KIPPENHAHN (1980) ESTIMATE
      do 40 i = zone_min,zone_max
         if (am_transport_convective_flag(i).and. &
              am_transport_convective_flag(i-1)) goto 40
! MHP 8/93 STABILITY CONDITION ADDED, NEGLECTING THE EFFECTS OF
! MU GRADIENTS.
         if (gsf_inhibition_mode.eq.2 .or. gsf_inhibition_mode.eq.0) then
            qwrmx = 2.0d0*sqrt(kinematic_viscosity_interface(i)/ &
                 thermal_diffusivity_interface(i))*dlnomega_dlnr_max(i)
            if (abs(dlnomega_dlnr(i)).lt.qwrmx) then
               circ_vel%gsf_circulation_velocity(i) = 0.0d0
               goto 40
            else
              fxx = sqrt((abs(dlnomega_dlnr(i))-qwrmx)/qwrmx)
            end if
         else
            fxx = 1.0d0
         end if
         dr = radius(i)-radius(i-1)
         if (gsf_inhibition_mode.eq.0) then
            qwrmx=2.0d0*sqrt(interface_gravity_factor(i)* &
                 abs(mix_phys%amum(i)-mix_phys%amum(i-1)) &
                 /dr/mean_molecular_weight_interface(i))
            if (abs(dlnomega_dlnr(i)).lt.qwrmx) then
               circ_vel%gsf_circulation_velocity(i) = 0.0d0
               goto 40
            else
! MHP 05/02 ADDED TESTS TO AVOID DIVIDE BY ZERO
!                  FXX2 = SQRT((ABS(QWLNR(I))-QWRMX)/QWRMX)
               qwrmx2 = 1.0d-2*dlnomega_dlnr(i)
               if (abs(qwrmx).lt.abs(qwrmx2)) then
                  qwrmx2 = 1.0d0
               else
                  if (abs(qwrmx).lt.1.0d-20) then
                     qwrmx2 = 1.0d0
                  else
                     qwrmx2 = sqrt((abs(dlnomega_dlnr(i))-qwrmx)/qwrmx)
                  end if
               end if
               fxx = min(fxx,qwrmx2)
            end if
         end if
         rmid = 0.5d0*(radius(i)+radius(i-1))
         dlnwdr = abs(log(omega(i))-log(omega(i-1)))/dr
! GSF IS TRIGGERED BY D OMEGA/DZ NONZERO (I.E. ROTATION NOT
! ON CYLINDERS,WHICH IS TRUE IN GENERAL IN OUR MODELS), OR
! BY D LN(J/M)/DR < 0, WHICH SOMETIMES OCCURS.
! IF D LN(J/M)/DR < 0, CHECK SECOND CRITERIA AND USE THE LARGEST
! VELOCITY; OTHERWISE,USE ONLY D OMEGA/DZ NON-ZERO CRITERION.
         if (omega(i)*radius(i)**2.lt.omega(i-1)*radius(i-1)**2) then
            dlnjmdr = abs(2.0d0/rmid+(log(omega(i))- &
                      log(omega(i-1)))/dr)
            fx = max(2.0d0*dlnjmdr,0.25d0*dlnwdr)
! FACTOR OF R IN THE DENOMINATOR OCCURS BECAUSE ALL THE CIRCULATION
! VELOCITIES ARE LATER MULTIPLIED BY R (THE "LENGTH SCALE").
            circ_vel%gsf_circulation_velocity(i) = gsf_kippenhahn_coeff(i)*fx* &
                 thermal_diffusivity_interface(i)*omega_interface(i)**2/rmid
         else
            circ_vel%gsf_circulation_velocity(i)=0.25d0*gsf_kippenhahn_coeff(i)* &
                 thermal_diffusivity_interface(i)*dlnwdr* &
                 omega_interface(i)**2/rmid
         end if
         circ_vel%gsf_circulation_velocity(i) = abs(fxx*circ_vel%gsf_circulation_velocity(i))
   40 continue
! OMIT JAMES AND KAHN ESTIMATE
!      ELSE
! JAMES AND KAHN (1971) ESTIMATE.
! NOTE FACTOR OF OMEGA**4, WHICH IS PART OF THE MERIDIONAL
! CIRCULATION VELOCITY ESTIMATE.
!         DO 41 I = IMIN,IMAX
!            IF(LCZ(I).AND.LCZ(I-1))GOTO 41
! MHP 8/93 STABILITY CONDITION ADDED, NEGLECTING THE EFFECTS OF
! MU GRADIENTS.
!            QWRMX = 2.0D0*SQRT(VISCMI(I)/THDIFMI(I))*QWRMAX(I)
!            IF(ABS(QWLNR(I)).LT.QWRMX)THEN
!               VGSF(I) = 0.0D0
!               GOTO 41
!            ELSE
!               FXX = SQRT((ABS(QWLNR(I))-QWRMX)/QWRMX)
!            ENDIF
!            DR = HRU(I)-HRU(I-1)
!            RMID = 0.5D0*(HRU(I)+HRU(I-1))
!            DLNJMDR = ABS(2.0D0/RMID+(LOG(OMEGA(I))-
!     *                LOG(OMEGA(I-1)))/DR)
!            IF(IGSF.EQ.3)THEN
!               DWDR = abs(OMEGA(I)-OMEGA(I-1))/DR
!               VGSF(I) = 2.0D0*FGSFJ(I)*DWDR**2
!               VGSF(I) = ABS(VGSF(I))
!            ELSE IF(IES.EQ.2)THEN
!               VGSF(I) = 2.0D0*VES(I)*FACT3(I)*DLNJMDR**2
!            ELSE
!               VGSF(I) = 2.0D0*FGSFJ(I)*(DLNJMDR*WM(I))**2
!               VGSF(I) = ABS(FXX*VGSF(I))
!            ENDIF
!   41    CONTINUE
!      ENDIF
!  DIFFUSIVE AND DYNAMICAL SHEAR INSTABILITIES - REF. ENDAL&SOFIA PAPER II.
      do 60 i = zone_min,zone_max
!  CHECK FOR OPERATION OF DYNAMICAL SHEAR.
!  IF DYNAMICAL SHEAR IS OPERATING,SET SECULAR SHEAR VELOCITY TO MAXIMUM
!  VALUE AND COMPUTE (LARGE) DYNAMICAL SHEAR VELOCITY.
         if (abs(dlnomega_dlnr(i)).gt.dlnomega_dlnr_max(i)) then
            qwr = abs(dlnomega_dlnr(i))
            circ_vel%secular_shear_velocity(i)=8.0d0/4.5d1* &
                 thermal_diffusivity_interface(i)* &
                 (qwr/dlnomega_dlnr_max(i))**2/interface_radius(i)
!  CORRECT GSF VELOCITY AS WELL.
            dr = radius(i)-radius(i-1)
            rmid = 0.5d0*(radius(i)+radius(i-1))
            if (omega(i).lt.omega(i-1)) then
               wmin = omega(i) + dlnomega_dlnr_max(i)*ln10* &
                    (log_radius(i)-log_radius(i-1))
               wmid = omega(i) + dlnomega_dlnr_max(i)*ln10* &
                    (log_radius(i)-log10(rmid))
            else
               wmin = omega(i) - dlnomega_dlnr_max(i)*ln10* &
                    (log_radius(i)-log_radius(i-1))
               wmid = omega(i) - dlnomega_dlnr_max(i)*ln10* &
                    (log_radius(i)-log10(rmid))
            end if
            dlnwdr0 = abs(log(omega(i))-log(omega(i-1)))/dr
            dlnjmdr0 = abs(2.0d0/rmid+(log(omega(i))- &
                         log(omega(i-1)))/dr)
            dlnjmdr = abs(2.0d0/rmid+(log(omega(i))-log(wmin))/dr)
            dlnwdr = abs(log(omega(i))-log(wmin))/dr
! MHP 05/02 ADDED FIRST BRANCH FOR IGSF=0 - SOMEHOW
! MISSED IN CODE CHANGES EARLIER!
! kippenhahn (1980) estimate
! OMIT BRANCH FOR JAMES AND KAHN.
!            IF(IGSF.EQ.2 .OR. IGSF.EQ.1 .OR. IGSF.EQ.0)THEN
            if (omega(i)*radius(i)**2.lt.wmin*radius(i-1)**2) then
               fx = max(2.0d0*dlnjmdr,0.25d0*dlnwdr)
               qwrmx_dyn = max(2.0d0*dlnjmdr0,0.25d0*dlnwdr0)
               circ_vel%gsf_circulation_velocity(i) = circ_vel%gsf_circulation_velocity(i)* &
                    fx/qwrmx_dyn
            else
               circ_vel%gsf_circulation_velocity(i)=circ_vel%gsf_circulation_velocity(i)* &
                    dlnwdr/dlnwdr0
            end if
!            ELSE IF(IGSF.EQ.3)THEN
!               DWDR = abs(OMEGA(I)-WMIN)/DR
!               DWDR0 = ABS(OMEGA(I)-OMEGA(I-1))/DR
!               VGSF(I) = VGSF(I)*(DWDR/DWDR0)**2
!            ELSE IF(IES.EQ.2)THEN
!               VGSF(I) = VGSF(I)*(DLNJMDR/DLNJMDR0)**2
!            ELSE
!               WMID0 = 0.5D0*(OMEGA(I)+OMEGA(I-1))
!               VGSF(I) = VGSF(I)*(DLNJMDR*WMID/DLNJMDR0/WMID0)**2
!            ENDIF
            write(6,9911) i,omega(i),omega(i-1),wmin
 9911       format(1x,'DYNAMICAL SHEAR-SHELL',i5,1p,' WTOP',e11.3, &
                 ' WBOT',e11.3,' LIMIT',e11.3)
            goto 60
! *** END OF CHANGED SECTION
         end if
!   FIND MAXIMUM GRADIENT IN OMEGA ALLOWED BY SECULAR SHEAR.
!   THE RUN OF QWRMAX INPUT IS THAT ALLOWED BY THE DYNAMICAL SHEAR;
!   THE SECULAR SHEAR RICHARDSON # IS RELATED BY
!   RICHNO(SECULAR) = PRANDTL# * CRITICAL REYNOLDS#/8 * RICHNO(DYNAMICAL)
!   PRANDTL # = KINEMATIC VISCOSITY/THERMOMETRIC DIFFUSIVITY.
! MHP 3/92 SQUARE ROOT OF PR# NEEDED, NOT PR # - ERROR CORRECTED!
! THE VELOCITY ESTIMATE HERE IS FROM ZAHN 1991.
         qwrmx = sqrt(kinematic_viscosity_interface(i)/ &
              thermal_diffusivity_interface(i)*1.25d-1*critical_reynolds)* &
              dlnomega_dlnr_max(i)
         if (abs(dlnomega_dlnr(i)).gt.qwrmx) then
!  UNSTABLE; CHECK FOR MU GRADIENTS.
            if (abs((mix_phys%amum(i)-mix_phys%amum(i-1))/mean_molecular_weight_interface(i)) &
                 .lt.1.0d-10) then
               qwrmx2 = 0.0d0
               qwr = abs(dlnomega_dlnr(i)) - qwrmx
               circ_vel%secular_shear_velocity(i)=8.0d0/4.5d1* &
                    thermal_diffusivity_interface(i)* &
                    (qwr/dlnomega_dlnr_max(i))**2/interface_radius(i)
            else
!  CHECK FOR EFFECTS OF MU GRADIENT.
!  RICHNO = RHO/P*(-d lnT/d lnMU)*(del MU)/(del P)*(G/QWLNR)**2
!  WHERE -d lnT/d lnMU = (CON-1)/(1+3CON),CON=(a/3)T**4/P
!  GIVEN (1 - CON)P = CGAS*RHO*T/MU
!      FACT = (RHOM/PM)*QTMU*DMU/AMUMI/DP*HGM**2
               qwrmx2 = 2.0d0*sqrt(max(1.0d-20,mu_gradient_richardson_coeff(i)* &
                        abs((mix_phys%amum(i)-mix_phys%amum(i-1))/ &
                        mean_molecular_weight_interface(i))))
               if (abs(dlnomega_dlnr(i)).gt.qwrmx2) then
!  INTERFACE UNSTABLE WITH RESPECT TO BOTH CONDITIONS; CHOOSE THE
!  MAXIMUM GRADIENT IMPLIED BY THE SECOND CONDITION IF IT'S LARGER
!  THAN THE FIRST(I.E. IF A MU GRADIENT IS SLOWING J TRANSPORT).
                  qwrmx = max(qwrmx2,qwrmx)
                  qwr = abs(dlnomega_dlnr(i)) - qwrmx
            circ_vel%secular_shear_velocity(i)=8.0d0/4.5d1* &
                 thermal_diffusivity_interface(i)* &
                 (qwr/dlnomega_dlnr_max(i))**2/interface_radius(i)
               end if
            end if
         end if
   60 continue
!  NOW DETERMINE WHETHER OR NOT MU GRADIENTS ARE STEEP ENOUGH TO
!  INHIBIT TRANSPORT.  MULTIPLY THE RESULTING VELOCITY ESTIMATES
!  BY THE USER DEFINED PARAMETERS FES AND FGSF.
!  IMU=3 KIPPENHAHN AND MOLLENHOF(1974)METHOD;IMU=2 LOCAL DAMPING
!  FACTOR METHOD.
! AGAIN, OMIT OBSOLETE MU GRADIENT TREATMENTS.
!      IF(IMU.EQ.3)THEN
!         IF(IT.GT.1)THEN
!            DO I = 1,M
!               VMU(I) = 0.5D0*(VMU(I)+VMUP(I))
!            END DO
!         ENDIF
!         DO I = IMIN,IMAX
!            FM = 1.0D0+VMU2(I)
!            RMID = 0.5D0*(HRU(I)+HRU(I-1))
!            FCC = SQRT(FC*FESC)
!            VEST=MAX(0.0D0,FES*VES(I)-SQRT(FES*VES(I)*RMID)
!     *               *FCC*VMU(I))
!            IF(VEST.LE.0.0D0)THEN
!               VES(I) = ABS(FES*VES(I)/FM)
!            ELSE
!               VES(I) = VEST
!            ENDIF
!            FCC = SQRT(FC*FGSFC)
!            IF(VGSF(I).GT.0.0D0)THEN
!               VGSFT=MAX(0.0D0,FGSF*VGSF(I)-SQRT(FGSF*VGSF(I)*RMID)
!     *               *FCC*VMU(I))
!            ELSE
!               VGSFT = 0.0D0
!            ENDIF
!            IF(VGSFT.LE.0.0D0)THEN
!               VGSF(I) = ABS(FES*VGSF(I)/FM)
!            ELSE
!               VGSF(I) = VGSFT
!            ENDIF
!            VSS(I)=MAX(0.0D0,FSS*VSS(I))
!         END DO
!      ELSE
!            IF(IMU.NE.2)THEN
!               FM = 1.0D0+VMU(I)
!            ELSE
!               FM = 1.0D0
!            ENDIF
!            VES(I)=ABS(FES*VES(I)/FM)
! MHP 05/02 ONLY DO THIS IF MU GRADIENTS NOT
! ALREADY ACCOUNTED FOR
!            VGSF(I)=ABS(FGSF*VGSF(I)/FM)
!            IF(IGSF.NE.0)THEN
!               VGSF(I)=ABS(FGSF*VGSF(I)/FM)
!            ELSE
!               VGSF(I) = FGSF*VGSF(I)
!            ENDIF
! ALREADY INCLUDED - ONLY USE SCALE FACTOR
!            VSS(I)=ABS(FSS*VSS(I)/FM)
!            VSS(I)= FSS*VSS(I)
!         END DO
! MHP 8/03 MULTIPLY VELOCITY ESTIMATES BY USER PARAMETER
! SCALE FACTORS
      do i = zone_min,zone_max
         circ_vel%es_circulation_velocity(i)=abs(es_velocity_scale* &
              circ_vel%es_circulation_velocity(i))
         circ_vel%gsf_circulation_velocity(i) = gsf_velocity_scale* &
              circ_vel%gsf_circulation_velocity(i)
         circ_vel%secular_shear_velocity(i)= secular_shear_velocity_scale* &
              circ_vel%secular_shear_velocity(i)
      end do
! MHP 11/94
! REPEAT FOR DIF+AD
      if (use_diffusion_advection_transport) then
         do i = zone_min,zone_max
               es_diffusive_velocity(i)=abs(es_velocity_scale* &
                    es_diffusive_velocity(i))
               es_advective_velocity(i)=es_velocity_scale* &
                    es_advective_velocity(i)
         end do
      end if
!               IF(IMU.NE.2)THEN
!                  FM = 1.0D0+VMU(I)
!               ELSE
!                  FM = 1.0D0
!               ENDIF
!               FM = 1.0D0+VMU(I)
!               VESD(I)=ABS(FES*VESD(I)/FM)
!               VESA(I)=FES*VESA(I)/FM
!            END DO
!         ENDIF
!      ENDIF
      if (use_diffusion_advection_transport) then
! MHP 05/02
! CHANGED TO REFLECT THE DIFFERENT TREATMENT OF THE GSF INSTABILITY.
! THE ENDAL AND SOFIA VARIANT HAD V=THDIF*DLNWDR*W**2/R WHICH WAS
! THEN ADDED TO A DIFFUSION EQUATION D/DR(CON*V*R*DW/DR).
! IN THE ZAHN FORMULISM, WE ARE SOLVING AN EQUATION OF THE FORM
! D/DR(CON*V*W).  BECAUSE THE VELOCITY ESTIMATE ITSELF INCLUDES
! DW/DR THIS IS IN EFFECT A DIFFUSION TERM.  SO...WE NEED TO
! MULTIPLY THE ORIGINAL VELOCITY BY W AND DIVIDE BY DW/DR TO CAST
! IT AS AN ENTRY IN THE DIFFUSION EQUATION.
         do i = zone_min,zone_max
            es_diffusive_velocity(i) = es_diffusive_velocity(i)+ &
                 interface_radius(i)*(abs(circ_vel%gsf_circulation_velocity(i))+ &
                 abs(circ_vel%secular_shear_velocity(i)))
            shear_diffusion_coeff(i) = interface_radius(i)* &
                 abs(circ_vel%secular_shear_velocity(i))
            gsf_diffusion_coeff(i) = interface_radius(i)* &
                 abs(circ_vel%gsf_circulation_velocity(i))
!            IF(VGSF(I).GT.0.0D0)THEN
! D LN W/DR
!               DR = HRU(I) - HRU(I-1)
!               QLNWR = ABS(LOG(OMEGA(I))-LOG(OMEGA(I-1))/DR)
! CEILING SET BY DYNAMICAL SHEAR
!               QLNWRMAX = ABS(QWRMAX(I)/(WM(I)*RM(I)))
! TAKE THE SMALLER OF THE TWO
!               QLNWR = MIN(QLNWR,QLNWRMAX)
!               IF(QLNWR.GT.1.0D-32)VESD(I) = VESD(I)+ABS(VGSF(I)/QLNWR)
!            ENDIF
!            IF(VSS(I).GT.0.0D0)THEN
!               VESD(I) = VESD(I)+ABS(VSS(I)*RM(I))
!            ENDIF
         end do
      end if
!  AVERAGE PREVIOUS AND NEW VELOCITY ESTIMATES AFTER THE FIRST ITERATION.
      if (iteration.gt.1) then
         do 70 i = zone_min,zone_max
            circ_vel%gsf_circulation_velocity(i) = 0.5d0*(circ_vel%gsf_circulation_velocity(i) &
                 + circ_vel%gsf_circulation_velocity_prev(i))
            circ_vel%es_circulation_velocity(i) = 0.5d0*(circ_vel%es_circulation_velocity(i) &
                 + circ_vel%es_circulation_velocity_prev(i))
            circ_vel%secular_shear_velocity(i) = 0.5d0*(circ_vel%secular_shear_velocity(i) &
                 + circ_vel%secular_shear_velocity_prev(i))
   70    continue
! MHP 11/94
         if (use_diffusion_advection_transport) then
            do i = zone_min,zone_max
               es_advective_velocity(i) = 0.5d0*(es_advective_velocity(i) &
                    + es_advective_velocity_prev(i))
               es_diffusive_velocity(i) = 0.5d0*(es_diffusive_velocity(i) &
                    + es_diffusive_velocity_prev(i))
            end do
         end if
      end if
      do 80 i =zone_min,zone_max
         total_circulation_velocity(i) = circ_vel%gsf_circulation_velocity(i) + &
              circ_vel%es_circulation_velocity(i) + circ_vel%secular_shear_velocity(i)
         if (total_circulation_velocity(i).lt.1.0d-20) &
              total_circulation_velocity(i)=0.0d0
         if (total_circulation_velocity(i).gt.0.0d0) any_transport_active=.true.
   80 continue
! 9/93 MIXING WITHOUT TRANSPORT ADDED.
! ZERO OUT COEFFICIENTS IN CORE TO AVOID NUMERICAL PROBLEMS IN
! THE H-BURNING SHELL.
      if (no_am_transport_in_core) then
         do i = zone_min,zone_max
            dell = interface_luminosity(i)/total_luminosity
            if (dell.lt.9.9d-1) then
               total_circulation_velocity(i) = 0.0d0
               circ_vel%gsf_circulation_velocity(i) = 0.0d0
               circ_vel%es_circulation_velocity(i) = 0.0d0
               circ_vel%secular_shear_velocity(i) = 0.0d0
! MHP 11/94
               es_advective_velocity(i) = 0.0d0
               es_diffusive_velocity(i) = 0.0d0
            else
               goto 81
            end if
         end do
 81      continue
      end if

      return
end subroutine vcirc
