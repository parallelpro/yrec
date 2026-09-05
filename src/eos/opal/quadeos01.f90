!----------------------------------------------------------------------
! quadeos01
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original quadeos01.f; only variable names, source form, and comment
! style were updated.
!
! Quadratic interpolation helper shared by the OPAL 2001 EOS table
! reader (esac01.f90/t6rinteos01.f90). See quad.f90 (the analogous
! OPAL 1995 helper) for a description of the recompute/cache scheme.
double precision function quadeos01(cache, recompute_flag, cache_slot, &
     eval_point, y1, y2, y3, x1, x2, x3)

      use opal_eos_lib
      implicit none

      type(opal_quad_cache), intent(inout) :: cache
      integer, intent(in) :: recompute_flag, cache_slot
      double precision, intent(in) :: eval_point, y1, y2, y3, x1, x2, x3

      double precision :: grid_x(3), grid_y(3)
      double precision :: coef2, coef1, coef0
      grid_x(1) = x1
      grid_x(2) = x2
      grid_x(3) = x3
      grid_y(1) = y1
      grid_y(2) = y2
      grid_y(3) = y3
      if (recompute_flag.eq.0) then
         cache%x12_inv(cache_slot) = 1.0d0/(grid_x(1) - grid_x(2))
         cache%x13_inv(cache_slot) = 1.0d0/(grid_x(1) - grid_x(3))
         cache%x23_inv(cache_slot) = 1.0d0/(grid_x(2) - grid_x(3))
         cache%x1_squared(cache_slot) = grid_x(1)*grid_x(1)
         cache%x1_plus_x2(cache_slot) = grid_x(1) + grid_x(2)
      end if
      coef2 = (grid_y(1) - grid_y(2))*cache%x12_inv(cache_slot)
      coef2 = coef2 - (grid_y(2) - grid_y(3))*cache%x23_inv(cache_slot)
      coef2 = coef2*cache%x13_inv(cache_slot)
      coef1 = (grid_y(1) - grid_y(2))*cache%x12_inv(cache_slot) - &
           cache%x1_plus_x2(cache_slot)*coef2
      coef0 = grid_y(1) - grid_x(1)*coef1 - cache%x1_squared(cache_slot)*coef2
      quadeos01 = coef0 + eval_point*(coef1 + eval_point*coef2)
      return
end function quadeos01
