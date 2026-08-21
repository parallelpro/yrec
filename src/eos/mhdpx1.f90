!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! MHDPX1
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!----------------------------------------------------------------------
! mhdpx1
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mhdpx1.f; only variable names, source form, and comment
! style were updated. common/luout/, common/mhdout/, and
! common/ccout2/ member names match those established in meqos.f90.
!
!     MHDST MUST BE CALLED IN MAIN.
!     INTERPOLATION IN TABLES WITH DIFFERENT X AND FIXED Z
!     TL < TLIM1:         LOWER PART OF ZAMS TABLES
!     TLIM1 < TL < TLIM2: UPPER PART OF ZAMS TABLES
!     TL > TLIM2:         CENTRE TABLES
!     TMINI,TMAXI:        TEMPERATURE INTERVAL COVERED
!                         BY THE TABLES
subroutine mhdpx1(log10_pressure, log10_temperature, hydrogen_fraction)
!
!     MHDST MUST BE CALLED IN MAIN.
!     INTERPOLATION IN TABLES WITH DIFFERENT X AND FIXED Z
      use numerics_lib
      implicit none
      integer, parameter :: ivarx = 25
      integer, parameter :: ndimt = 8

      double precision, intent(in) :: log10_pressure, log10_temperature, &
           hydrogen_fraction

      double precision :: table_vars(ndimt,ivarx), table_hfrac(ndimt)
! common/luout/: only short_file_unit and main_output_unit are used
! here. Naming matches meqos.f90.
      integer :: ilast, idebug, itrack, short_file_unit, imilne, imodpt, &
           istor, main_output_unit
      common/luout/ ilast, idebug, itrack, short_file_unit, imilne, &
           imodpt, istor, main_output_unit

      logical :: ldebug, lcorr, lmilne, ltrack, lstpch
      common/ccout2/ ldebug, lcorr, lmilne, ltrack, lstpch
! common/tttt/: temperature boundaries between the lower/upper ZAMS
! table regions and the centre-table region, and the overall
! temperature range covered by the loaded tables. Set by mhdst.
      double precision :: zams_lower_upper_boundary_log10t, &
           zams_centre_boundary_log10t, table_log10t_min, table_log10t_max
      common/tttt/zams_lower_upper_boundary_log10t, &
           zams_centre_boundary_log10t, table_log10t_min, table_log10t_max
!     QUANTITIES FOR INTERPOLATION IN X
      double precision :: cubic_vars(ivarx,4), cubic_x_nodes(4)
      integer :: cubic_table_index(4)
! LIR's "type" flag; the file's own IMPLICIT LOGICAL*4(L) rule is
! overridden below by an explicit INTEGER*4 declaration for L.
      integer :: lir_type_flag
!     OUTPUT
      double precision :: mhd_output(ivarx)
      common/mhdout/ mhd_output
      save
!     READ FROM APPROPRIATE TABLES
!     AND FILL ARRAYS VAROUT AND XC
!     CAN WE DO IT?
!
      integer :: i, itbl, iv, ixmin, num_vars, var_leading_dim, &
           num_points, interp_mode
      double precision :: x_grid_origin, x_grid_spacing, &
           var_at_x0, var_at_x1, var_at_x2

!     IRANGE = 1
      if (log10_temperature.lt.table_log10t_min .or. &
           log10_temperature.gt.table_log10t_max) then
!         IRANGE = 0
          go to 999
      end if
!     LOWER ZAMS TABLES
      if (log10_temperature.lt.zams_lower_upper_boundary_log10t) then
         do 10 i=1,3
         itbl = -i
         call mhdpx2(log10_pressure, log10_temperature, itbl, table_vars, table_hfrac, ndimt)
 10      continue
      else if (log10_temperature.lt.zams_centre_boundary_log10t) then
!     UPPER ZAMS TABLES
         do 20 i=1,3
         itbl = i
         call mhdpx2(log10_pressure, log10_temperature, itbl, table_vars, table_hfrac, ndimt)
 20      continue
      else
!     CENTRE TABLES
         do 30 i=1,5
!        OFFSET (+3) TO ACCESS CENTER TABLES
         itbl = i + 3
         call mhdpx2(log10_pressure, log10_temperature, itbl, table_vars, table_hfrac, ndimt)
 30      continue
      end if
!     INTERPOLATION IN X
      if (log10_temperature.le.zams_centre_boundary_log10t) then
!        QUADRATIC NEWTON (EQUIDISTANT XC'S)
         x_grid_origin = table_hfrac(1)
         x_grid_spacing  = table_hfrac(2) - table_hfrac(1)
         if (abs(table_hfrac(3)-table_hfrac(2)-x_grid_spacing).gt.1.d-4) then
            write(main_output_unit,*) 'ERROR (MHD): NON-EQUIDISTANT ZAMS TABLES.'
            write(main_output_unit,*) 'XC(1-3)= ',table_hfrac(1),table_hfrac(2),table_hfrac(3)
            write(short_file_unit,*) 'ERROR (MHD): NON-EQUIDISTANT ZAMS TABLES.'
            write(short_file_unit,*) 'XC(1-3)= ',table_hfrac(1),table_hfrac(2),table_hfrac(3)
            stop
         end if
!         IF( ( XC(3).GT.XC(1) .AND. (X.GT.XC(3) .OR. X.LT.XC(1)))
!    1                               .OR.
!    2       ( XC(3).LT.XC(1) .AND. (X.GT.XC(1) .OR. X.LT.XC(3))) )
!    3   THEN
!        END IF
         do 100 iv=1,ivarx
         var_at_x0 = table_vars(1,iv)
         var_at_x1 = table_vars(2,iv)
         var_at_x2 = table_vars(3,iv)
         call quint(hydrogen_fraction, x_grid_origin, x_grid_spacing, &
              var_at_x0, var_at_x1, var_at_x2, mhd_output(iv))
 100     continue
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
         do 200 i=1,4
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!  200     IX(I)=IXMIN+I-1
            cubic_table_index(i)=ixmin+i-1
 200     continue
         do 250 i =1,4
! KC 2025-05-30 fixed "Shared DO termination label"
!         XX(I) = XC(IX(I))
!         DO 250 IV=1,IVARX
!         VARI1(IV,I) = VAROUT(IX(I),IV)
            cubic_x_nodes(i) = table_hfrac(cubic_table_index(i))
            do 251 iv=1,ivarx
               cubic_vars(iv,i) = table_vars(cubic_table_index(i),iv)
 251        continue
 250     continue
         num_vars=ivarx
         var_leading_dim=ivarx
         num_points=4
         lir_type_flag=1
         interp_mode=1
         call lir(hydrogen_fraction, cubic_x_nodes, mhd_output, cubic_vars, &
              num_vars, var_leading_dim, num_points, lir_type_flag, interp_mode)
      end if
  999 return
!1000  FORMAT(' RESULTS FROM MHDPX2, ITBL,X = ',I6,1PE15.7/)
! 1001  FORMAT(12(/1X,1P5E15.6))
! 5001  FORMAT(1X,'******* WARNING: EXTRAPOLATION IN X (QUINT) ',
!      1          'PL,TL,X = ',3F12.6)
! 5011  FORMAT(1X,'******* WARNING: EXTRAPOLATION IN X (LIR) ',
!      1          'PL,TL,X = ',3F12.6)
! 9001  FORMAT(' ERROR IN MHDPX1. TL OUT OF RANGE. TL,TMINI,TMAXI=',
!      1 /1X,1P3E13.6)
end subroutine mhdpx1
