!----------------------------------------------------------------------
! get_imp_diffco
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original get_imp_diffco.f; only variable names, source form, and
! comment style were updated. Validated against the Stage 0 regression
! suite (examples/run_standard_solar_model).
!
! Rebuilds the tridiagonal-matrix coefficients for the implicit
! (second-derivative) term of the microscopic-diffusion equation
! solved in microdiff_run.f90, correcting the zone-midpoint diffusion
! coefficient for the change in abundance found on the previous
! iteration. Also called by grsett.f (not part of this batch).
subroutine get_imp_diffco(alpha, diffusion_coeff_mid, delta_abundance_mid, &
     diffusion_coeff_deriv_mid, sub_diag, diag, super_diag, npt)

      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: alpha(json)
      double precision, intent(inout) :: diffusion_coeff_mid(json)
      double precision, intent(in) :: delta_abundance_mid(json)
      double precision, intent(in) :: diffusion_coeff_deriv_mid(json)
      double precision, intent(out) :: sub_diag(json), diag(json), &
           super_diag(json)
      integer, intent(in) :: npt
      integer :: i

! CORRECT DIFFUSION COEFFICEINTS FOR CHANGE IN X IN THE PREVIOUS ITERATION.
      do i=1,npt-1
         diffusion_coeff_mid(i)=diffusion_coeff_mid(i)+ &
              diffusion_coeff_deriv_mid(i)*delta_abundance_mid(i)
   10 continue
      end do
! NOW RECOMPUTE ELEMENTS OF THE TRIDIAGONAL MATRIX SYSTEM.
      sub_diag(1) = 0.0d0
      diag(1) = 1.0d0 + alpha(1)*diffusion_coeff_mid(1)
      super_diag(1) = -alpha(1)*diffusion_coeff_mid(1)
!   911 FORMAT(5X,I5,1P3E10.2)
      do i=2,npt-1
         sub_diag(i) = -alpha(i)*diffusion_coeff_mid(i-1)
         diag(i) = 1.0d0 + alpha(i)*(diffusion_coeff_mid(i-1)+ &
              diffusion_coeff_mid(i))
         super_diag(i) = -alpha(i)*diffusion_coeff_mid(i)
   20 continue
      end do
      sub_diag(npt) = -alpha(npt)*diffusion_coeff_mid(npt-1)
      diag(npt) = 1.0d0 + alpha(npt)*diffusion_coeff_mid(npt-1)
      super_diag(npt) = 0.0d0
      return
end subroutine get_imp_diffco
