!----------------------------------------------------------------------
! dcoeft
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original dcoeft.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
!  DCOEFT SETS UP THE COEFFICIENTS FOR THE DIFFUSION DIFFERENCE EQUATION.
!
!  INPUT VARIABLES:
!  diffusion_coeff (ECOD) : DIFFUSION COEFFICIENTS AT THE EQUALLY SPACED
!     GRID POINTS.
!  grid_spacing (DR) : GRID SPACING.
!  timestep (DT) : TIMESTEP.
!  eq_moment_of_inertia (EI) : RUN OF MOMENTS OF INERTIA OF EQUALLY
!     SPACED GRID POINTS.
!  eq_angular_momentum (EJ) : RUN OF ANGULAR MOMENTUM OF EQUALLY SPACED
!     GRID POINTS (see rotgrid.f90's eq_angular_momentum).
!  eq_omega (EW) : RUN OF ANGULAR VELOCITY OF EQUALLY SPACED GRID POINTS.
!  num_eq_points (NTOT) : NUMBER OF EQUALLY SPACED GRID POINTS.
!  wind_loss_explicit (WIND1) : THE ANGULAR MOMENTUM LOSS FROM A
!     SURFACE C.Z. COMPUTED EXPLICITILY.
!  wind_loss_implicit (WIND2) : AS WIND1, BUT COMPUTED IMPLICITLY.
!  *NOTE: IF NOT APPLICABLE, WIND1 AND WIND2 ARE ZEROED OUT.
!
!  OUTPUT VARIABLES (via common/tridi/, see mixcom.f90/tridia.f90) :
!  THE ANGULAR VELOCITY OF SHELL I AT TIME N+1 (W(I,N+1)) IS A FUNCTION OF
!  W(I-1),W(I),AND W(I+1) AS DISCUSSED BELOW.  THIS CAN BE EXPRESSED
!  AS A TRIDIAGONAL MATRIX.
!  sub_diag(I) : CONTAINS ALL TERMS INVOLVING OMEGA(I-1).
!  diag(I) : CONTAINS ALL TERMS INVOLVING OMEGA(I).
!  super_diag(I) : CONTAINS ALL TERMS INVOLVING OMEGA(I+1)
!  rhs(I) : THE ANGULAR VELOCITY TERMS AT THE BEGINNING OF THE TIMESTEP.
!  *NOTE: ANGULAR MOMENTUM LOSS FROM A SURFACE C.Z. IS SUBTRACTED FROM
!     THE LAST ELEMENT OF ARRAY rhs IF APPLICABLE.
!
!  THE DIFFUSION EQUATION WE ARE SOLVING IS
!  dW/dT = 1/(4pi*RHO*R**2) 1/(I/M) d/dR(D*dW/dR)
!  WHERE D IS OUR DIFFUSION COEFFICIENT,W IS THE ANGULAR VELOCITY,
!  I IS THE MOMENT OF INERTIA,M IS MASS,T IS TIME,R IS RADIUS,
!  AND RHO IS DENSITY.
!  USING M = 4pi*RHO*R**2*DR THIS IS DIFFERENCED FOR SHELL I AS
!  W(I,N+1)-W(I,N)=(DT/DR)*(1/I)*(D(I+1/2)*(W(I+1,N+1)-W(I,N+1)) -
!  D(I-1/2)*(W(I,N+1)-W(I-1,N+1)))
!  WHERE INDEX N+1 REFERS TO VALUES AT THE END OF THE TIMESTEP,AND
!  INDEX N REFERS TO VALUES AT THE BEGINNING OF THE TIMESTEP.
!  THIS SYSTEM IS SUPPLEMENTED WITH APPROPRIATE BOUNDARY CONDITIONS.
!
! Simpler (single diffusion coefficient, no third/fourth-order terms)
! companion to dadcoeft.f90's 4-band angular-momentum-transport solve,
! and analogue of ccoeft.f90 (the composition-transport version of the
! same tridiagonal setup).
subroutine dcoeft(diffusion_coeff, grid_spacing, timestep, &
     eq_moment_of_inertia, eq_angular_momentum, eq_omega, num_eq_points, &
     wind_loss_explicit, wind_loss_implicit, fix_omega_at_surface)
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: diffusion_coeff(json)
      double precision, intent(in) :: grid_spacing, timestep
      double precision, intent(in) :: eq_moment_of_inertia(json), &
           eq_angular_momentum(json), eq_omega(json)
      integer, intent(in) :: num_eq_points
      double precision, intent(in) :: wind_loss_explicit, wind_loss_implicit

! common/tridi/: tridiagonal-solve work arrays (Thomas algorithm).
! sub_diag/diag/super_diag/rhs are filled in here; solution is
! solver-internal (not touched here). gamma_elim is used here purely
! as scratch storage for the surface wind-loss term (matching the
! original's reuse of GAMA(NTOT) for this purpose, before tridia.f90
! overwrites it as solver work space). Naming matches mixcom.f90/
! tridia.f90.
      double precision :: sub_diag(json), diag(json), super_diag(json), &
           rhs(json), solution(json), gamma_elim(json)
      common/tridi/ sub_diag, diag, super_diag, rhs, solution, gamma_elim

! common/difad/: am_advective_coeff/am_diffusive_coeff (originally
! ECOD3/ECOD4), used only under the diffusion+advection treatment.
! Naming matches dadcoeft.f90/rotgrid.f90.
      double precision :: am_advective_coeff(json), am_diffusive_coeff(json)
      common/difad/ am_advective_coeff, am_diffusive_coeff

! MHP 7/93
! common/varfc/: only use_diffusion_advection_transport (LDIFAD) is
! used here; vfc/lvfc are unused placeholders. Naming matches
! rotgrid.f90/vcirc.f90.
      double precision :: vfc(json)
      logical :: lvfc, use_diffusion_advection_transport
      common/varfc/ vfc, lvfc, use_diffusion_advection_transport

      logical, intent(in) :: fix_omega_at_surface
      save

      integer :: i
      double precision :: fact0, fact, facta

      fact0 = timestep/grid_spacing
      do 10 i = 1,num_eq_points
         rhs(i) = eq_omega(i)
   10 continue
      gamma_elim(num_eq_points) = -0.5d0*(wind_loss_explicit+ &
           wind_loss_implicit)*eq_moment_of_inertia(num_eq_points)
      rhs(num_eq_points) = rhs(num_eq_points)*(1.0d0+ &
           gamma_elim(num_eq_points)/eq_angular_momentum(num_eq_points))
      if (.not.use_diffusion_advection_transport) then
!  FIRST SHELL BOUNDARY CONDITIONS: NO FLOW BELOW THE BOUNDARY
!  I.E. THE ANGULAR MOMENTUM TRANSPORT AT THE FIRST SHELL DOES NOT
!  DEPEND ON THE SHELLS BELOW IT.
         fact = fact0/eq_moment_of_inertia(1)
         sub_diag(1) = 0.0d0
         diag(1) = 1.0d0 + fact*diffusion_coeff(2)
         super_diag(1) = -fact*diffusion_coeff(2)
         do 20 i = 2,num_eq_points-1
            fact = fact0/eq_moment_of_inertia(i)
            sub_diag(i) = -fact*diffusion_coeff(i)
            diag(i) = 1.0d0 + fact*(diffusion_coeff(i)+diffusion_coeff(i+1))
            super_diag(i) = -fact*diffusion_coeff(i+1)
   20    continue
!  LAST SHELL B.C. : SAME AS FIRST SHELL B.C.
         fact = fact0/eq_moment_of_inertia(num_eq_points)
         if (.not.fix_omega_at_surface) then
            sub_diag(num_eq_points) = -fact*diffusion_coeff(num_eq_points)
            diag(num_eq_points) = 1.0d0 + fact*diffusion_coeff(num_eq_points)
         else
            sub_diag(num_eq_points) = 0.0d0
            diag(num_eq_points) = 1.0d0
         end if
         super_diag(num_eq_points) = 0.0d0
      else
! IF LDIFAD=T, WE ARE SOLVING A COMBINED DIFFUSION/ADVECTION EQUATION.
! THIS ADDS A TERM D/DR(RHO*R**4*V*W) TO THE ORIGINAL D/DR(RHO*R**4*
! V*R*DW/DR) EQUATION.
!  FIRST SHELL BOUNDARY CONDITIONS: NO FLOW BELOW THE BOUNDARY
!  I.E. THE ANGULAR MOMENTUM TRANSPORT AT THE FIRST SHELL DOES NOT
!  DEPEND ON THE SHELLS BELOW IT.
         fact = fact0/eq_moment_of_inertia(1)
         facta = 0.5d0*timestep/eq_moment_of_inertia(1)
         sub_diag(1) = 0.0d0
         diag(1) = 1.0d0 + fact*am_diffusive_coeff(2) - &
              facta*am_advective_coeff(2)
         super_diag(1) = -fact*am_diffusive_coeff(2) - &
              facta*am_advective_coeff(2)
         do i = 2,num_eq_points-1
            fact = fact0/eq_moment_of_inertia(i)
            facta = 0.5d0*timestep/eq_moment_of_inertia(i)
            sub_diag(i) = -fact*am_diffusive_coeff(i) + &
                 facta*am_advective_coeff(i)
            diag(i) = 1.0d0 + fact*(am_diffusive_coeff(i)+ &
                 am_diffusive_coeff(i+1)) + &
                 facta*(am_advective_coeff(i)-am_advective_coeff(i+1))
            super_diag(i) = -fact*am_diffusive_coeff(i+1) - &
                 facta*am_advective_coeff(i+1)
         end do
!  LAST SHELL B.C. : SAME AS FIRST SHELL B.C.
         fact = fact0/eq_moment_of_inertia(num_eq_points)
         facta = 0.5d0*timestep/eq_moment_of_inertia(num_eq_points)
         if (.not.fix_omega_at_surface) then
            sub_diag(num_eq_points) = -fact*am_diffusive_coeff(num_eq_points) &
                 + facta*am_advective_coeff(num_eq_points)
            diag(num_eq_points) = 1.0d0 + fact*am_diffusive_coeff(num_eq_points) &
                 + facta*am_advective_coeff(num_eq_points)
         else
            sub_diag(num_eq_points) = 0.0d0
            diag(num_eq_points) = 1.0d0
         end if
         super_diag(num_eq_points) = 0.0d0
      end if

      return
end subroutine dcoeft
