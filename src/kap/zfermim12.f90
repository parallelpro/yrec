!----------------------------------------------------------------------
! zfermim12
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original zfermim12.f; only variable names, source form, and comment
! style were updated.
!
!..this routine applies a rational function expansion to get the
!..fermi-dirac integral of order -1/2 evaluated at degeneracy_parameter.
!..maximum error is 1.23d-12. reference: antia apjs 84,101 1993.
double precision function zfermim12(degeneracy_parameter)

      implicit none

! degeneracy_parameter: electron degeneracy parameter eta at which the
!                        order -1/2 Fermi-Dirac integral is evaluated.
      double precision, intent(in) :: degeneracy_parameter

!..declare
      integer :: term_idx, deg_num_small, deg_den_small, deg_num_large, &
           deg_den_large
!       double precision fd_order,coef_num_small(12),coef_den_small(12),coef_num_large(12),coef_den_large(12),numerator,denominator,xx  ! KC 2025-05-31
      double precision :: coef_num_small(12), coef_den_small(12), &
           coef_num_large(12), coef_den_large(12), numerator, denominator, &
           scaled_arg
      save

!..load the coefficients of the expansion
!       data  fd_order,deg_num_small,deg_den_small,deg_num_large,deg_den_large /-0.5d0, 7, 7, 11, 11/  ! KC 2025-05-31
      data  deg_num_small,deg_den_small,deg_num_large,deg_den_large /7, 7, 11, 11/
      data  (coef_num_small(term_idx),term_idx=1,8)/ 1.71446374704454d7,    3.88148302324068d7, &
           3.16743385304962d7,    1.14587609192151d7, &
           1.83696370756153d6,    1.14980998186874d5, &
           1.98276889924768d3,    1.0d0/
      data  (coef_den_small(term_idx),term_idx=1,8)/ 9.67282587452899d6,    2.87386436731785d7, &
           3.26070130734158d7,    1.77657027846367d7, &
           4.81648022267831d6,    6.13709569333207d5, &
           3.13595854332114d4,    4.35061725080755d2/
      data (coef_num_large(term_idx),term_idx=1,12)/-4.46620341924942d-15, -1.58654991146236d-12, &
           -4.44467627042232d-10, -6.84738791621745d-8, &
           -6.64932238528105d-6,  -3.69976170193942d-4, &
           -1.12295393687006d-2,  -1.60926102124442d-1, &
           -8.52408612877447d-1,  -7.45519953763928d-1, &
           2.98435207466372d0,    1.0d0/
      data (coef_den_large(term_idx),term_idx=1,12)/-2.23310170962369d-15, -7.94193282071464d-13, &
           -2.22564376956228d-10, -3.43299431079845d-8, &
           -3.33919612678907d-6,  -1.86432212187088d-4, &
           -5.69764436880529d-3,  -8.34904593067194d-2, &
           -4.78770844009440d-1,  -4.99759250374148d-1, &
           1.86795964993052d0,    4.16485970495288d-1/

      if (degeneracy_parameter .lt. 2.0d0) then
       scaled_arg = exp(degeneracy_parameter)
       numerator = scaled_arg + coef_num_small(deg_num_small)
       do term_idx=deg_num_small-1,1,-1
        numerator = numerator*scaled_arg + coef_num_small(term_idx)
       enddo
       denominator = coef_den_small(deg_den_small+1)
       do term_idx=deg_den_small,1,-1
        denominator = denominator*scaled_arg + coef_den_small(term_idx)
       enddo
       zfermim12 = scaled_arg * numerator/denominator
!..
      else
       scaled_arg = 1.0d0/(degeneracy_parameter*degeneracy_parameter)
       numerator = scaled_arg + coef_num_large(deg_num_large)
       do term_idx=deg_num_large-1,1,-1
        numerator = numerator*scaled_arg + coef_num_large(term_idx)
       enddo
       denominator = coef_den_large(deg_den_large+1)
       do term_idx=deg_den_large,1,-1
        denominator = denominator*scaled_arg + coef_den_large(term_idx)
       enddo
       zfermim12 = sqrt(degeneracy_parameter)*numerator/denominator
      end if
      return
end function zfermim12
