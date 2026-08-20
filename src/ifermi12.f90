!----------------------------------------------------------------------
! ifermi12
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ifermi12.f; only variable names, source form, and comment
! style were updated.
!
! This routine applies a rational function expansion to get the
! inverse Fermi-Dirac integral of order 1/2 (i.e. the electron
! degeneracy parameter eta such that F_{1/2}(eta) = fermi_half_integral)
! when it is equal to fermi_half_integral.
! Maximum error is 4.19d-9. Reference: Antia, ApJS 84,101 1993.
double precision function ifermi12(fermi_half_integral)

      implicit none

! fermi_half_integral: value of the order-1/2 Fermi-Dirac integral
!                       whose inverse (the degeneracy parameter) is
!                       sought.
      double precision, intent(in) :: fermi_half_integral

!..declare
      integer :: term_idx, deg_num_small, deg_den_small, deg_num_large, &
           deg_den_large
      double precision :: fd_order, coef_num_small(12), coef_den_small(12), &
           coef_num_large(12), coef_den_large(12), numerator, denominator, &
           scaled_arg
!     unused in the original: z,drn
      save

!..load the coefficients of the expansion
      data  fd_order,deg_num_small,deg_den_small,deg_num_large,deg_den_large &
           /0.5d0, 4, 3, 6, 5/
      data  (coef_num_small(term_idx),term_idx=1,5)/ 1.999266880833d4,   5.702479099336d3, &
           6.610132843877d2,   3.818838129486d1, &
           1.0d0/
      data  (coef_den_small(term_idx),term_idx=1,4)/ 1.771804140488d4,  -2.014785161019d3, &
           9.130355392717d1,  -1.670718177489d0/
      data  (coef_num_large(term_idx),term_idx=1,7)/-1.277060388085d-2,  7.187946804945d-2, &
           -4.262314235106d-1,  4.997559426872d-1, &
           -1.285579118012d0,  -3.930805454272d-1, &
           1.0d0/
      data  (coef_den_large(term_idx),term_idx=1,6)/-9.745794806288d-3,  5.485432756838d-2, &
           -3.299466243260d-1,  4.077841975923d-1, &
           -1.145531476975d0,  -6.067091689181d-2/

      if (fermi_half_integral .lt. 4.0d0) then
       numerator  = fermi_half_integral + coef_num_small(deg_num_small)
       do term_idx=deg_num_small-1,1,-1
        numerator  = numerator*fermi_half_integral + coef_num_small(term_idx)
       enddo
       denominator = coef_den_small(deg_den_small+1)
       do term_idx=deg_den_small,1,-1
        denominator = denominator*fermi_half_integral + coef_den_small(term_idx)
       enddo
       ifermi12 = log(fermi_half_integral * numerator/denominator)

      else
       scaled_arg = 1.0d0/fermi_half_integral**(1.0d0/(1.0d0 + fd_order))
       numerator = scaled_arg + coef_num_large(deg_num_large)
       do term_idx=deg_num_large-1,1,-1
        numerator = numerator*scaled_arg + coef_num_large(term_idx)
       enddo
       denominator = coef_den_large(deg_den_large+1)
       do term_idx=deg_den_large,1,-1
        denominator = denominator*scaled_arg + coef_den_large(term_idx)
       enddo
       ifermi12 = numerator/(denominator*scaled_arg)
      end if
      return
end function ifermi12
