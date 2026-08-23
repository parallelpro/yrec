!----------------------------------------------------------------------
! ccoeft
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ccoeft.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
!  SET UP COEFFICIENTS OF DIFFUSION DIFFERENCE EQUATIONS FOR COMPOSITION
!  TRANSPORT DUE TO ROTATIONALLY INDUCED MIXING.
!
!  INPUT VARIABLES:
!  diffusion_coeff (COD) : RUN OF DIFFUSION COEFFICIENTS (COMPUTED IN
!     SR COEFFT).
!  grid_spacing (DR) : GRID SPACING IN RADIUS(CM).
!  timestep (DT) : DIFFUSION TIMESTEP (SEC).
!  eq_composition (ECOMP) : STARTING ARRAY OF COMPOSITION MASS FRACTION
!     AT THE EQUALLY SPACED GRID.
!  eq_mass (EM) : RUN OF MASSES OF THE EQUALLY SPACED GRID POINTS (GM).
!  num_eq_points (NTOT) : NUMBER OF POINTS IN THE EQUALLY SPACED GRID
!     USED TO SOLVE THE DIFFUSION EQUATION.
!  OUTPUT VARIABLES (sub_diag/diag/super_diag/rhs dummy arguments) :
!  THE ABUNDANCE OF SHELL I AT TIME N+1 (COMP(I,N+1)) IS A FUNCTION OF
!  COMP(I-1),COMP(I),AND COMP(I+1) AS DISCUSSED BELOW.  THIS CAN BE EXPRESSED
!  AS A TRIDIAGONAL MATRIX.  TERMS WHICH DEPEND ON COMP(I-1) ARE STORED IN
!  ARRAY sub_diag; THOSE WHICH DEPEND ON COMP(I) ARE STORED IN ARRAY diag;
!  THOSE WHICH DEPEND ON COMP(I+1) ARE STORED IN ARRAY super_diag; AND THE
!  STARTING RUN OF COMPOSITION IS STORED IN ARRAY rhs.
!
!  THE DIFFUSION EQUATION WE ARE SOLVING IS
!  dX/dT = 1/(4pi*RHO*R**2)d/dR(D*dX/dR)
!  WHERE D IS OUR DIFFUSION COEFFICIENT,X IS THE MASS FRACTION OF THE
!  SPECIES IN QUESTION,T IS TIME,R IS RADIUS, AND RHO IS DENSITY.
!  USING M = 4pi*RHO*R**2*DR THIS IS DIFFERENCED FOR SHELL I AS
!  X(I,N+1)-X(I,N)=(DT/DR)*(1/M)*(D(I+1/2)*(X(I+1,N+1)-X(I,N+1)) -
!  D(I-1/2)*(X(I,N+1)-X(I-1,N+1)))
!  WHERE M IS THE MASS OF SHELL , INDEX N+1 REFERS TO VALUES AT THE
!  END OF THE TIMESTEP, INDEX N REFERS TO VALUES AT THE BEGINNING OF
!  THE TIMESTEP.
!  THIS SYSTEM IS SUPPLEMENTED WITH APPROPRIATE BOUNDARY CONDITIONS.
!
! Composition-transport analogue of dcoeft.f90's angular-momentum-
! transport tridiagonal setup; called by mixcom.f90, which passes the
! output matrix straight on to ctridi.f90.
!
! sub_diag/diag/super_diag/rhs were originally shared with mixcom.f90/
! ctridi.f90 via common/tridi/ (positional storage); converted (2026,
! GUIDELINES.md) to explicit output arguments since this is real
! per-call data flow, not global configuration.
subroutine ccoeft(diffusion_coeff, grid_spacing, timestep, eq_composition, &
     eq_mass, num_eq_points, sub_diag, diag, super_diag, rhs)
      use const_lib
      use luout_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: diffusion_coeff(json)
      double precision, intent(in) :: grid_spacing, timestep
      double precision, intent(in) :: eq_composition(json), eq_mass(json)
      integer, intent(in) :: num_eq_points
      double precision, intent(out) :: sub_diag(json), diag(json), &
           super_diag(json), rhs(json)
      integer :: i
      double precision :: fact0, fact

!  STORE THE RUN OF COMPOSITION AT THE BEGINNING OF THE TIMESTEP IN ARRAY rhs.
      do i = 1,num_eq_points
         rhs(i) = eq_composition(i)
      end do
!  FIRST SHELL BOUNDARY CONDITIONS: NO FLOW BELOW THE BOUNDARY
!  I.E. THE ANGULAR MOMENTUM TRANSPORT AT THE FIRST SHELL DOES NOT
!  DEPEND ON THE SHELLS BELOW IT.(SAME USED FOR LAST SHELL B.C.)
      fact0 = timestep/grid_spacing
      fact = fact0/eq_mass(1)
      sub_diag(1) = 0.0d0
      diag(1) = 1.0d0 + fact*diffusion_coeff(2)
      super_diag(1) = -fact*diffusion_coeff(2)
      do i = 2,num_eq_points-1
         fact = fact0/eq_mass(i)
         sub_diag(i) = -fact*diffusion_coeff(i)
         diag(i) = 1.0d0 + fact*(diffusion_coeff(i)+diffusion_coeff(i+1))
         super_diag(i) = -fact*diffusion_coeff(i+1)
      end do
!  LAST SHELL BOUNDARY CONDITIONS.
      fact = fact0/eq_mass(num_eq_points)
      sub_diag(num_eq_points) = -fact*diffusion_coeff(num_eq_points)
      diag(num_eq_points) = 1.0d0 + fact*diffusion_coeff(num_eq_points)
      super_diag(num_eq_points) = 0.0d0

      return
end subroutine ccoeft
