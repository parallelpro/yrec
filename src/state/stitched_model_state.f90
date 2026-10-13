!----------------------------------------------------------------------
! stitched_model_state
!----------------------------------------------------------------------
! The record of the stitched model (core/diagnostics/stitched_model.f90):
! the whole converged star -- interior, envelope and atmosphere -- on one
! extended grid, which the star holds as s%stitched. build_stitched_model
! fills it; the profile and pulse writers, the seismic observables and
! the turnover time read it.
!   n_ext            the points of the extended grid
!   n_ie             the last point that is not in the atmosphere
!   region(i)        point i's region: 1 interior, 2 envelope, 3 atmosphere
!   index(i)         its index within that region
!   seismic(4,i)     brunt_N2, lamb_S2 (l = 1), gradL, gradr_div_grada
!   atm_height(i)    the height of atmosphere point i above the photosphere
!                    (cm; indexed like s%atm_points, 1 = outermost)
!   prof(col,i)      the profile columns (io/yrec_output's registry)
!   pulse(col,i)     the GYRE/FGONG pulse columns
! The arrays are allocated, and zeroed, by the first build; columns and
! points a build does not write keep their earlier values.
module stitched_model_state
      implicit none
      private

      type, public :: stitched_model_record
            integer :: n_ext = 0, n_ie = 0
            integer, allocatable :: region(:), index(:)
            double precision, allocatable :: seismic(:,:), atm_height(:), &
                 prof(:,:), pulse(:,:)
      end type stitched_model_record

end module stitched_model_state
