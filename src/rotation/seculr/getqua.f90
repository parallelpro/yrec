!----------------------------------------------------------------------
! getqua
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original getqua.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Solves for the quadrupole moment of a rotating star, using the
! Sweet (1950) formulation as given in Zahn (1992). We are solving
! the equation
! 1/R D**2/DR2 (R*QUAD) -6*QUAD/R**2 - 4*PI*G*RHO*D(RHO)/DP * QUAD =
!  -4/3*PI*G*R**2/(GM/R**2) * D/DR (RHO*OMEGA**2)
! This can be expressed as D**2/DR2(QUAD) + 2/R DQUAD/DR -6*QUAD/R**2
! +4*PI*G*RHO*QUAD = 4*PI*G*RHO*R**2 *(OMEGA**2 - 2*OMEGA/G *DOMEGA/DR)
! using hydrostatic equilibrium and the ideal gas law
! D RHO/DR = DP/DR = -RHO*G.
! This defines a tridiagonal matrix system of the form
! A*QUAD(I-1)+B*QUAD(I)+C*QUAD(I+1) = D, which can be solved in
! the standard way for the run of QUAD(I) with appropriate boundary
! conditions. The central B.C. is that QUAD(0)=0; the surface B.C.
! is that there is no contribution to the quadrupole from outside the
! final shell, such that QUAD(I) varies as 1/R**4.
subroutine getqua(log_density, gravity, radius, angular_velocity, num_zones)

      use star_info_lib, only: star
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: log_density(json), gravity(json), &
           radius(json), angular_velocity(json)
      integer, intent(in) :: num_zones



! Tridiagonal-solve work arrays (Thomas algorithm) for this file's own
! self-contained inline solve below (not shared with any other file --
! was originally common/tridi/, which is genuinely shared elsewhere
! (ctridi.f90/tridia.f90/ccoeft.f90/dcoeft.f90/mixcom.f90), but this
! file never called those solvers, just reused the block's memory
! layout; converted (2026, GUIDELINES.md) to plain locals since there
! was never any real data flow to another file here).
      double precision :: sub_diag(json), diag(json), super_diag(json), &
           rhs(json), solution(json), gamma_elim(json)


      double precision :: density_omega2(json)
      save

      integer :: zone_index, matrix_row
      double precision :: four_pi_g
      double precision :: dr, weight_plus, weight_minus, inv_dr2, dr_inv_r
      double precision :: drho_dr
      double precision :: dr_plus, dr_minus, d_density_omega2
      double precision :: radius_plus, outer_boundary_ratio
      double precision :: pivot

      do zone_index = 1,num_zones
         density_omega2(zone_index) = exp(ln10*log_density(zone_index))* &
              angular_velocity(zone_index)**2
      end do
      four_pi_g = exp(ln10*(c4pil+cgl))
! CENTER NUMERICAL DERIVATIVES, USING THE SMALLER OF THE LEFT AND
! RIGHT HANDED INTERVALS.
      dr = min(radius(1),radius(2)-radius(1))
      weight_plus = dr/(radius(2)-radius(1))
      weight_minus = dr/radius(1)
      inv_dr2 = 1.0D0/dr**2
      dr_inv_r = 1.0D0/dr/radius(1)
! CENTRAL BOUNDARY CONDITION : QUADRUPOLE GOES TO ZERO.
      sub_diag(1) = 0.0D0
      drho_dr = (exp(ln10*log_density(2))-exp(ln10*log_density(1)))/dr
      diag(1) = -inv_dr2*(weight_plus+weight_minus)+dr_inv_r* &
           (weight_minus-weight_plus)-6.0D0/radius(1)**2 &
           -four_pi_g*drho_dr/gravity(1)
      super_diag(1) = weight_plus*(inv_dr2+dr_inv_r)
      rhs(1) = -four_pi_g*drho_dr*(radius(1)*angular_velocity(1))**2/3.0D0/ &
           gravity(1)
! GENERAL CASE : SECOND DERIVATIVE NUMERICALLY DIFFERENTIATED AS
!
      do zone_index = 2,num_zones-1
         dr_plus = radius(zone_index+1)-radius(zone_index)
         dr_minus = radius(zone_index)-radius(zone_index-1)
         dr = min(dr_minus,dr_plus)
         weight_plus = dr/dr_plus
         weight_minus = dr/dr_minus
         inv_dr2 = 1.0D0/dr**2
         dr_inv_r = 1.0D0/dr/radius(zone_index)
         drho_dr = 0.5D0*(weight_plus*exp(ln10*log_density(zone_index+1))+ &
              (weight_minus-weight_plus)* &
              exp(ln10*log_density(zone_index))-weight_minus* &
              exp(ln10*log_density(zone_index)))/dr
         sub_diag(zone_index) = weight_minus*(inv_dr2-dr_inv_r)
         diag(zone_index) = -inv_dr2*(weight_plus+weight_minus)+dr_inv_r* &
              (weight_minus-weight_plus)-6.0D0/radius(zone_index)**2 &
              -four_pi_g*drho_dr/gravity(zone_index)
         super_diag(zone_index) = weight_plus*(inv_dr2+dr_inv_r)
         d_density_omega2 = 0.5D0*(weight_plus*density_omega2(zone_index+1) + &
              (weight_minus-weight_plus)*density_omega2(zone_index) - &
              weight_minus*density_omega2(zone_index-1))/dr
         rhs(zone_index) = -d_density_omega2*four_pi_g* &
              radius(zone_index)**2/3.0D0/gravity(zone_index)
      end do
! SURFACE B.C. D PHI/DR + 3 PHI/R = 0
! THE SECOND AND THIRD TERMS IN THE EXPRESSION CANCEL IN THIS CASE.
! OUTSIDE THE STAR PHI(R) = PHI(RSTAR)*(RSTAR/R)**3.
      dr = radius(num_zones) - radius(num_zones-1)
      radius_plus = radius(num_zones)+dr
      outer_boundary_ratio = (radius(num_zones)/radius_plus)**3
      inv_dr2 = 1.0D0/dr**2
      sub_diag(num_zones) = inv_dr2
      diag(num_zones) = (outer_boundary_ratio-2.0D0)*inv_dr2-1.2D1/ &
           radius(num_zones)**2
      super_diag(num_zones) = 0.0D0
      rhs(num_zones) = 0.0D0
! SOLUTION OF SYSTEM FROM NUMERICAL RECIPES.
      if (diag(1).eq.0.0D0) stop
      pivot = diag(1)
      solution(1) = rhs(1)/pivot
      do matrix_row = 2,num_zones
         gamma_elim(matrix_row) = super_diag(matrix_row-1)/pivot
         pivot = diag(matrix_row) - sub_diag(matrix_row)*gamma_elim(matrix_row)
         if (pivot.eq.0.0D0) stop
         solution(matrix_row) = (rhs(matrix_row) - &
              sub_diag(matrix_row)*solution(matrix_row-1))/pivot
      end do
! BACKSUBSTITUTION.
      do matrix_row = num_zones-1,1,-1
         solution(matrix_row) = solution(matrix_row) - &
              gamma_elim(matrix_row+1)*solution(matrix_row+1)
      end do
      do zone_index = 1,num_zones
         star%rot%quadrupole_moment(zone_index) = solution(zone_index)
         star%rot%local_gravity(zone_index) = gravity(zone_index)
      end do

      return
end subroutine getqua
