!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! MHDPX1
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!----------------------------------------------------------------------
! mhdpx1
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mhdpx1.f; only variable names, source form, and comment
! style were updated.
!
!     MHDST MUST BE CALLED IN MAIN.
!     INTERPOLATION IN TABLES WITH DIFFERENT X AND FIXED Z
!     TL < TLIM1:         LOWER PART OF ZAMS TABLES
!     TLIM1 < TL < TLIM2: UPPER PART OF ZAMS TABLES
!     TL > TLIM2:         CENTRE TABLES
!     TMINI,TMAXI:        TEMPERATURE INTERVAL COVERED
!                         BY THE TABLES
!     (TLIM1/TLIM2/TMINI/TMAXI are mhd_eos%zams_lower_upper_boundary_log10t,
!     zams_centre_boundary_log10t, table_log10t_min, table_log10t_max.)
!     An out-of-range TL returns silently with ierr = 0 and
!     mhd_eos%mhd_output untouched (historical behaviour, preserved).
subroutine mhdpx1(log10_pressure, log10_temperature, hydrogen_fraction, ierr)
      use mhd_eos_lib
      use luout_lib
      use numerics_lib
      implicit none
      integer, parameter :: ivarx = 25
      integer, parameter :: ndimt = 8

      double precision, intent(in) :: log10_pressure, log10_temperature, &
           hydrogen_fraction

      double precision :: table_vars(ndimt,ivarx), table_hfrac(ndimt)

!     QUANTITIES FOR INTERPOLATION IN X
      double precision :: cubic_vars(ivarx,4), cubic_x_nodes(4)
      integer :: cubic_table_index(4)
! LIR's "type" flag (an integer, despite the original's L-prefixed name).
      integer :: lir_type_flag
!     READ FROM APPROPRIATE TABLES
!     AND FILL ARRAYS VAROUT AND XC
!     CAN WE DO IT?
!
      integer :: i, itbl, iv, ixmin, num_vars, var_leading_dim, &
           num_points, interp_mode
      double precision :: x_grid_origin, x_grid_spacing, &
           var_at_x0, var_at_x1, var_at_x2

      integer, intent(out) :: ierr

      ierr = 0

      if (log10_temperature.lt.mhd_eos%table_log10t_min .or. &
           log10_temperature.gt.mhd_eos%table_log10t_max) then
          return
      end if
!     LOWER ZAMS TABLES
      if (log10_temperature.lt.mhd_eos%zams_lower_upper_boundary_log10t) then
         do i=1,3
         itbl = -i
         call mhdpx2(log10_pressure, log10_temperature, itbl, table_vars, table_hfrac, ndimt)
         end do
      else if (log10_temperature.lt.mhd_eos%zams_centre_boundary_log10t) then
!     UPPER ZAMS TABLES
         do i=1,3
         itbl = i
         call mhdpx2(log10_pressure, log10_temperature, itbl, table_vars, table_hfrac, ndimt)
         end do
      else
!     CENTRE TABLES
         do i=1,5
!        OFFSET (+3) TO ACCESS CENTER TABLES
         itbl = i + 3
         call mhdpx2(log10_pressure, log10_temperature, itbl, table_vars, table_hfrac, ndimt)
         end do
      end if
!     INTERPOLATION IN X
      if (log10_temperature.le.mhd_eos%zams_centre_boundary_log10t) then
!        QUADRATIC NEWTON (EQUIDISTANT XC'S)
         x_grid_origin = table_hfrac(1)
         x_grid_spacing  = table_hfrac(2) - table_hfrac(1)
         if (abs(table_hfrac(3)-table_hfrac(2)-x_grid_spacing).gt.1.d-4) then
            write(terminal_unit,*) 'ERROR (MHD): NON-EQUIDISTANT ZAMS TABLES.'
            write(terminal_unit,*) 'XC(1-3)= ',table_hfrac(1),table_hfrac(2),table_hfrac(3)
            write(run_log_unit,*) 'ERROR (MHD): NON-EQUIDISTANT ZAMS TABLES.'
            write(run_log_unit,*) 'XC(1-3)= ',table_hfrac(1),table_hfrac(2),table_hfrac(3)
            ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
            ! facades stop when their caller passes no ierr.
            ierr = 1
            return
         end if
         do iv=1,ivarx
         var_at_x0 = table_vars(1,iv)
         var_at_x1 = table_vars(2,iv)
         var_at_x2 = table_vars(3,iv)
         call quint(hydrogen_fraction, x_grid_origin, x_grid_spacing, &
              var_at_x0, var_at_x1, var_at_x2, mhd_eos%mhd_output(iv))
         end do
      else
!     CUBIC LAGRANGIAN (ARBITRARILY SPACED XC'S)
!     USE 4 HIGHEST XC(I)'S, IF POSSIBLE.
!     DISTINGUISH CASE WHERE THE XC(I), I=4..8,
!     INCREASE OR DECREASE.
         if (table_hfrac(4).lt.table_hfrac(8)) then
            ixmin = 5
            if (hydrogen_fraction.lt.table_hfrac(5)) ixmin=4
         else
            ixmin = 4
            if (hydrogen_fraction.lt.table_hfrac(7)) ixmin=5
         end if
         do i=1,4
            cubic_table_index(i)=ixmin+i-1
         end do
         do i =1,4
            cubic_x_nodes(i) = table_hfrac(cubic_table_index(i))
            do iv=1,ivarx
               cubic_vars(iv,i) = table_vars(cubic_table_index(i),iv)
            end do
         end do
         num_vars=ivarx
         var_leading_dim=ivarx
         num_points=4
         lir_type_flag=1
         interp_mode=1
         call lir(hydrogen_fraction, cubic_x_nodes, mhd_eos%mhd_output, cubic_vars, &
              num_vars, var_leading_dim, num_points, lir_type_flag, interp_mode)
      end if
      return
end subroutine mhdpx1
