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
!  OUTPUT VARIABLES (sub_diag/diag/super_diag/rhs dummy arguments) :
!  THE ANGULAR VELOCITY OF SHELL I AT TIME N+1 (W(I,N+1)) IS A FUNCTION OF
!  W(I-1),W(I),AND W(I+1) AS DISCUSSED BELOW.  THIS CAN BE EXPRESSED
!  AS A TRIDIAGONAL MATRIX.
!  sub_diag(I) : CONTAINS ALL TERMS INVOLVING OMEGA(I-1).
!  diag(I) : CONTAINS ALL TERMS INVOLVING OMEGA(I).
!  super_diag(I) : CONTAINS ALL TERMS INVOLVING OMEGA(I+1)
!  rhs(I) : THE ANGULAR VELOCITY TERMS AT THE BEGINNING OF THE TIMESTEP.
!  *NOTE: ANGULAR MOMENTUM LOSS FROM A SURFACE C.Z. IS SUBTRACTED FROM
!     THE LAST ELEMENT OF ARRAY rhs IF APPLICABLE, AND ALSO RETURNED
!     VIA surface_wind_loss_term FOR tridia.f90 TO SEED dj(n) WITH
!     (SEE tridia.f90'S HEADER NOTE -- THIS WAS PREVIOUSLY SMUGGLED
!     THROUGH common/tridi/'s gamma_elim(num_eq_points) SLOT).
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
! sub_diag/diag/super_diag/rhs/surface_wind_loss_term were originally
! shared with tridia.f90 (via seculr.f90, which calls both) through
! common/tridi/ (positional storage; surface_wind_loss_term was
! smuggled through the gamma_elim(num_eq_points) slot specifically so
! tridia.f90's first statement could read it before overwriting
! gamma_elim as pure solver scratch). Converted (2026, GUIDELINES.md)
! to explicit output arguments since this is real per-call data flow,
! not global configuration -- see tridia.f90's matching dj_n_seed
! input argument.
subroutine dcoeft(diffusion_coeff, grid_spacing, timestep, &
     eq_moment_of_inertia, eq_angular_momentum, eq_omega, num_eq_points, &
     wind_loss_explicit, wind_loss_implicit, fix_omega_at_surface, &
     sub_diag, diag, super_diag, rhs, surface_wind_loss_term)
      use star_info_lib, only: star
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: diffusion_coeff(json)
      double precision, intent(in) :: grid_spacing, timestep
      double precision, intent(in) :: eq_moment_of_inertia(json), &
           eq_angular_momentum(json), eq_omega(json)
      integer, intent(in) :: num_eq_points
      double precision, intent(in) :: wind_loss_explicit, wind_loss_implicit
      double precision, intent(out) :: sub_diag(json), diag(json), &
           super_diag(json), rhs(json)
      double precision, intent(out) :: surface_wind_loss_term



      logical, intent(in) :: fix_omega_at_surface
      integer :: i
      double precision :: fact0, fact, facta

      fact0 = timestep/grid_spacing
      do i = 1,num_eq_points
         rhs(i) = eq_omega(i)
      end do
      surface_wind_loss_term = -0.5d0*(wind_loss_explicit+ &
           wind_loss_implicit)*eq_moment_of_inertia(num_eq_points)
      rhs(num_eq_points) = rhs(num_eq_points)*(1.0d0+ &
           surface_wind_loss_term/eq_angular_momentum(num_eq_points))
      if (.not.use_diffusion_advection_transport) then
!  FIRST SHELL BOUNDARY CONDITIONS: NO FLOW BELOW THE BOUNDARY
!  I.E. THE ANGULAR MOMENTUM TRANSPORT AT THE FIRST SHELL DOES NOT
!  DEPEND ON THE SHELLS BELOW IT.
         fact = fact0/eq_moment_of_inertia(1)
         sub_diag(1) = 0.0d0
         diag(1) = 1.0d0 + fact*diffusion_coeff(2)
         super_diag(1) = -fact*diffusion_coeff(2)
         do i = 2,num_eq_points-1
            fact = fact0/eq_moment_of_inertia(i)
            sub_diag(i) = -fact*diffusion_coeff(i)
            diag(i) = 1.0d0 + fact*(diffusion_coeff(i)+diffusion_coeff(i+1))
            super_diag(i) = -fact*diffusion_coeff(i+1)
         end do
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
         diag(1) = 1.0d0 + fact*star%rot%am_diffusive_coeff(2) - &
              facta*star%rot%am_advective_coeff(2)
         super_diag(1) = -fact*star%rot%am_diffusive_coeff(2) - &
              facta*star%rot%am_advective_coeff(2)
         do i = 2,num_eq_points-1
            fact = fact0/eq_moment_of_inertia(i)
            facta = 0.5d0*timestep/eq_moment_of_inertia(i)
            sub_diag(i) = -fact*star%rot%am_diffusive_coeff(i) + &
                 facta*star%rot%am_advective_coeff(i)
            diag(i) = 1.0d0 + fact*(star%rot%am_diffusive_coeff(i)+ &
                 star%rot%am_diffusive_coeff(i+1)) + &
                 facta*(star%rot%am_advective_coeff(i)-star%rot%am_advective_coeff(i+1))
            super_diag(i) = -fact*star%rot%am_diffusive_coeff(i+1) - &
                 facta*star%rot%am_advective_coeff(i+1)
         end do
!  LAST SHELL B.C. : SAME AS FIRST SHELL B.C.
         fact = fact0/eq_moment_of_inertia(num_eq_points)
         facta = 0.5d0*timestep/eq_moment_of_inertia(num_eq_points)
         if (.not.fix_omega_at_surface) then
            sub_diag(num_eq_points) = -fact*star%rot%am_diffusive_coeff(num_eq_points) &
                 + facta*star%rot%am_advective_coeff(num_eq_points)
            diag(num_eq_points) = 1.0d0 + fact*star%rot%am_diffusive_coeff(num_eq_points) &
                 + facta*star%rot%am_advective_coeff(num_eq_points)
         else
            sub_diag(num_eq_points) = 0.0d0
            diag(num_eq_points) = 1.0d0
         end if
         super_diag(num_eq_points) = 0.0d0
      end if

      return
end subroutine dcoeft
