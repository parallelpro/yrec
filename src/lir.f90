!----------------------------------------------------------------------
! lir
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original lir.f; only variable names, source form, and comment style
! were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Generic table-lookup interpolation/extrapolation routine, entered
! either as LIR (cubic interpolation/extrapolation, unless num_points
! < 4) or as LIR1 (always linear interpolation/extrapolation).
!
! FOR A SUCH THAT target_z=table_z(A),  SETS result_y(I)=table_y(I,A), I=1,num_y
! table_z(N),table_y(I,N) MUST BE SUPPLIED FOR N=1,num_points AND I=1,num_y
! y_stride IS FIRST DIMENSION OF table_y
! interp_flag IS SET TO 1 FOR INTERPOLATION AND 0 FOR EXTRAPOLATION
! IF continue_search.LE.1, SCAN TO FIND THE table_z(N) WHICH IMMEDIATELY
!   BOUND target_z, STARTING AT N=1
! IF continue_search.GT.1, SCAN STARTS FROM VALUE OF N FROM PREVIOUS
!   CALL OF LIR
! NOTE
! MOST OF THE COMPUTATION IS PERFORMED IN SINGLE PRECISION
subroutine lir(target_z,table_z,result_y,table_y,num_y,y_stride, &
     num_points,continue_search,interp_flag)

      implicit none

      double precision, intent(in) :: target_z
      double precision, intent(in) :: table_z(*)
      double precision, intent(out) :: result_y(*)
      double precision, intent(in) :: table_y(*)
      integer, intent(in) :: num_y, y_stride, num_points, continue_search
      integer, intent(out) :: interp_flag

      double precision :: weight(4)
      integer :: search_idx
      data search_idx/-1/
      save

      integer :: linear_mode
      integer :: stride, stride_m1, y_strided, num_y_strided, table_end
      double precision :: diff
      integer :: pivot, closest
      double precision :: y1, y2, y3, y4, z1, z2, z3, z4, z12, z34
      integer :: y_idx, j, k, base_idx
      double precision :: yy

      linear_mode=0
      go to 1
      entry lir1(target_z,table_z,result_y,table_y,num_y,y_stride, &
           num_points,continue_search,interp_flag)
      linear_mode=1
    1 continue
      stride=1
! CHECK NT AND RESET IL IF NECESSARY
      if(num_points.lt.2) go to 101
      if(num_points.lt.4) linear_mode=1
! ADDRESSING CONSTANTS
      interp_flag=1
      stride_m1=stride-1
      y_strided=stride*y_stride
      num_y_strided=(num_y-1)*stride+1
      table_end=(num_points-1)*stride+1
      diff=table_z(table_end)-table_z(1)
! SET INDEX FOR START OF SEARCH
      search_idx=(search_idx-2)*stride+1
      if(continue_search.le.1.or.search_idx.lt.1) search_idx=1
! DETERMINE POSITION OF target_z WITHIN table_z
    2 if(search_idx.gt.table_end) go to 8
! KC 2025-05-30 fixed "Arithmetic IF statement"
!       IF(DIFF) 4,102,3
!     3 IF(ZI(N)-Z) 5,6,9
!     4 IF(ZI(N)-Z) 9,6,5
      if (diff .lt. 0.0) then
         goto 4
      else if (diff .eq. 0.0) then
         goto 102
      else
         goto 3
      end if
    3 if (table_z(search_idx) .lt. target_z) then
         goto 5
      else if (table_z(search_idx) .eq. target_z) then
         goto 6
      else
         goto 9
      end if
    4 if (table_z(search_idx) .lt. target_z) then
         goto 9
      else if (table_z(search_idx) .eq. target_z) then
         goto 6
      else
         goto 5
      end if
    5 search_idx=search_idx+stride
      go to 2
! SET Y WHEN Z LIES ON A MESH POINT
    6 base_idx=(search_idx-1)*y_stride
      do 7 y_idx=1,num_y_strided
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!       Y(I)=YI(I+J)
!     7 IF(Y(I).EQ.0.D0) Y(I+IR1)=0.D0
         result_y(y_idx)=table_y(y_idx+base_idx)
         if(result_y(y_idx).eq.0.d0) result_y(y_idx+stride_m1)=0.d0
    7 continue
      go to 30
! CONTROL WHEN Z DOES NOT LIE ON A MESH POINT
    8 interp_flag=0
    9 if(search_idx.le.1) interp_flag=0
      if(linear_mode.eq.1) go to 20
! CUBIC INTERPOLATION/EXTRAPOLATION
! PIVOTAL POINT (M) AND POINT (K) CLOSEST TO Z
!    10 M=N
      pivot=search_idx
      closest=3
      if(search_idx.gt.1+stride) go to 11
      pivot=1+stride+stride
      closest=search_idx
   11 if(search_idx.lt.table_end) go to 12
      pivot=table_end-stride
      closest=4
! WEIGHTING FACTORS
   12 y1=table_z(pivot-stride*2)
      y2=table_z(pivot-stride)
      y3=table_z(pivot)
      y4=table_z(pivot+stride)
      z1=target_z-y1
      z2=target_z-y2
      z3=target_z-y3
      z4=target_z-y4
!    13 Z12=Z1*Z2
      z12=z1*z2
      z34=z3*z4
!    14 A(1)=Z2*Z34/((Y1-Y2)*(Y1-Y3)*(Y1-Y4))
      weight(1)=z2*z34/((y1-y2)*(y1-y3)*(y1-y4))
      weight(2)=z1*z34/((y2-y1)*(y2-y3)*(y2-y4))
      weight(3)=z12*z4/((y3-y1)*(y3-y2)*(y3-y4))
      weight(4)=z12*z3/((y4-y1)*(y4-y2)*(y4-y3))
! CORRECT A(K)
!    15 DIFF=A(1)+A(2)+A(3)+A(4)
      diff=weight(1)+weight(2)+weight(3)+weight(4)
      weight(closest)=(1.d0+weight(closest))-diff
! COMPUTE Y
!    16 M=(M-1)/IR-3
      pivot=(pivot-1)/stride-3
      pivot=pivot*y_strided
      do 18 y_idx=1,num_y_strided
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!       K=I+M
!       YY=0.D0
!       DO 17 J=1,4
!       K=K+IRD
!       DIFF=YI(K)
!    17 YY=YY+A(J)*DIFF
!       Y(I)=YY
!    18 IF(Y(I).EQ.0.D0) Y(I+IR1)=0.D0
         k=y_idx+pivot
         yy=0.d0
         do 17 j=1,4
            k=k+y_strided
            diff=table_y(k)
            yy=yy+weight(j)*diff
   17    continue
         result_y(y_idx)=yy
         if(result_y(y_idx).eq.0.d0) result_y(y_idx+stride_m1)=0.d0
   18 continue
      go to 30
! LINEAR INTERPOLATION/EXTRAPOLATION
   20 if(search_idx.eq.1) search_idx=1+stride
      if(search_idx.gt.table_end) search_idx=table_end
      z1=table_z(search_idx)
      y1=(z1-target_z)/(z1-table_z(search_idx-stride))
      y2=1.0d0-y1
      base_idx=(search_idx-1)*y_stride
      pivot=base_idx-y_strided
      do 21 y_idx=1,num_y_strided,stride
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!       Y(I)=Y1*YI(I+M)+Y2*YI(I+J)
!    21 IF(Y(I).EQ.0.D0) Y(I+IR1)=0.D0
         result_y(y_idx)=y1*table_y(y_idx+pivot)+y2*table_y(y_idx+base_idx)
         if(result_y(y_idx).eq.0.d0) result_y(y_idx+stride_m1)=0.d0
   21 continue
! RESET N
   30 search_idx=(search_idx+stride-1)/stride
      return
! DIAGNOSTICS
  101 continue
      return
  102 continue
      return
!  1001 FORMAT(/1X,10('*'),5X,'THERE ARE FEWER THAN TWO DATA POINTS IN',
!      *      ' LIR     NT =',I4,5X,10('*')/)
!  1002 FORMAT(/1X,10('*'),5X,'EXTREME VALUES OF INDEPENDENT VARIABLE',
!      *      ' EQUAL IN LIR',5X,10('*')/16X,'ZI(1) =',1PE13.5,',   ',
!      *       'ZI(',I4,') =',1PE13.5/)
end subroutine lir
