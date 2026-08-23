!------------------------    GROUP: SR_ALL   -------------------------------
!
!----------------------------------------------------------------------
! rtab
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original rtab.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Reads an unformatted opacity/EOS-style table of shape
! (num_t_points, num_rho_points, num_vars) from file_unit: for each of
! num_t_points temperature records, reads the record's density-point
! count (asserted equal for every record) and log10(T), then reads
! num_vars variables at each of its density points.
subroutine rtab(file_unit,max_t_points,max_rho_points,num_vars, &
     num_t_points,num_rho_points,log_t,table_data, ierr)
      use const_lib
      implicit none

      integer, intent(in) :: file_unit, max_t_points, max_rho_points, &
           num_vars, num_t_points
      integer, intent(out) :: num_rho_points

!     NT IS INPUT; NR,TL,TDVAR ARE OUTPUT
      double precision, intent(out) :: log_t(max_t_points), &
           table_data(max_t_points,max_rho_points,num_vars)
! --- locals ---
      integer :: t_idx, rho_idx, var_idx, rho_count_read

      integer, intent(out) :: ierr

      ierr = 0

      do t_idx = 1, num_t_points
!     READ(IR,1001) NRR,TL(N)
      read(file_unit     ) rho_count_read,log_t(t_idx)
      if(t_idx.eq.1) num_rho_points=rho_count_read
!     CHECK LIMITS AND NUMBER OF DENSITIY POINTS OF TABLE
      if ( num_t_points.gt.max_t_points .or.  num_rho_points.gt.max_rho_points ) then
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      if ( t_idx.gt.1 .and. rho_count_read.ne.num_rho_points ) then
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
!     END CHECK ....................................................
! KC 2025-05-30 fixed "Shared DO termination label"
!       DO 10 J = 1, NR
      do rho_idx = 1, num_rho_points
!     READ(IR,1002) (TDVAR(N,J,IV),IV=1,IVAR)
      read(file_unit     ) (table_data(t_idx,rho_idx,var_idx),var_idx=1,num_vars)
      end do
      end do
      return
! 1001  FORMAT(I5,E16.8)
! 1002  FORMAT(5E16.8)
! 1011  FORMAT(I5,1PE16.8)
! 1012  FORMAT(1P5E16.8)
! 9001  FORMAT(' ERROR IN RTAB. TOO SMALL PARAMETERS: ',
!      1 /' NT,NR,NTM,NRM = ',4I5)
! 9011  FORMAT(' ERROR IN RTAB. WRONG DENSITY CONSTRUCTION: ',
!      1 /' N,NR,NRR = ',3I5)
end subroutine rtab
