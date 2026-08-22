!----------------------------------------------------------------------
! getopal95
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original getopal95.f; only variable names, source form, and comment
! style were updated.
!
! THIS IS THE INTERPOLATION FACILITY FOR THE LIVERMORE OPACITY TABLES
! USING CUBIC SPLINE INTERPLOATION SCHEME.   (FOR OPAL95)
! --------------------------------------------------------------
!      SUBROUTINE GETOPAL95(DL,TL,X,Z,O,OL,QOD,QOT) DRIVER- GET CAPPA
!      SUBROUTINE OP952D(O,OL,QOD,QOT) 2D INTERPOLATION IN (T,RHO)
!      SUBROUTINE OP953D(O,OL,QOD,QOT) 3D INTERPOLATION IN (X,T,RHO)
!      SUBROUTINE OP954D(O,OL,QOD,QOT) 4D INTERPOLATION IN (Z,X,T,RHO)
!      SUBROUTINE LL95TBL READ IN ALL TABLES
!      SUBROUTINE OP95ZTAB(Z) GENERATE TABLE AT FIXED Z
!      SUBROUTINE OP95XTAB(X) GENERATE TABLE AT FIXED X,Z
!
! 7/98 MHP DRIVER ROUTINE FOR OPACITY GIVEN RHO,T, X, AND Z.
! Determines whether a 2D (rho,T), 3D (X,rho,T), or 4D (Z,X,rho,T)
! interpolation is required and dispatches to op952d/op953d/op954d
! (not part of this batch) accordingly.
subroutine getopal95(log10_density, log10_temperature, hydrogen_fraction, &
     metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ierr)

      use opacity_table_lib
      use const_lib
      use numerics_lib
      implicit none
      integer, parameter :: num_t = 70
      integer, parameter :: num_d = 19
      integer, parameter :: num_x = 10
      integer, parameter :: num_z = 13
      integer, parameter :: num_xz = 126
! JVS Need this one too:
      integer, parameter :: json = 5000

      double precision, intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      integer, intent(out) :: ierr

! END JVS
      double precision :: interp_nodes(4), weight(4), dweight(4)
      integer :: i, k, j
      logical :: low_regime_flag, density_shifted
      integer :: x_table_index, density_base_index, temp_row_index, x_shift_base
      double precision :: extrap_logr, z_at_table

!     ENSURE THAT WE ARE WITHIN THE OPAL 95 TABLES.
!     COMPUTE LOG R = RHO/T6**3
      ierr = 0
      opacity_table%opal95_logr = log10_density - 3.0d0*(log10_temperature-6.0d0)
!     CHECK T
      if (log10_temperature.lt.opacity_table%opal95_grid_logt(1) .or. &
           log10_temperature.gt.opacity_table%opal95_grid_logt(num_t)) then
         write(*,5) log10_temperature
    5    format(' LOG T OF',f11.6,'OUT OF OPAL 95 TABLE RANGE')
! 2026 (ROADMAP.md stage 3): stop converted to ierr; the facade
! (kap_lib's kap_get) stops when its caller passes no ierr.
         ierr = 1
         return
      endif
!     CHECK TO SEE IF EXTRAPOLATION BELOW THE FIRST TABLE ELEMENT
!     IN DENSITY IS NEEDED.
      if (opacity_table%opal95_logr.lt.opacity_table%opal95_grid_logr(1)) then
         opacity_table%opal95_extrap_lo = .true.
      else
         opacity_table%opal95_extrap_lo = .false.
      endif
!     THIS SECTION IS REPLACED BY LINEAR EXTRAPOLATION IN R FROM THE
!     LAST TWO OPACITY VALUES.
!      IF(RL.LT.(R0GR(1)-4.0D0) .OR. RL.GT.(R0GR(NUMD)+1.0D0))THEN
!         WRITE(*,7)RL,TL
! 7       FORMAT(' LOG R OF',F11.6,' OUT OF OPAL 95 TABLE RANGE;TL'
!     *          ,F11.6)
!         STOP
!      ENDIF
!     GET INTERPOLATION FACTORS IN T
      if (log10_temperature.ge.opacity_table%opal95_grid_logt(opacity_table%opal95_index_t+2)) then
         do i = opacity_table%opal95_index_t+3,num_t-1
            if (log10_temperature.lt.opacity_table%opal95_grid_logt(i)) then
               opacity_table%opal95_index_t = i - 2
               exit
            endif
         end do
         if (i > (num_t-1)) then
         opacity_table%opal95_index_t = num_t - 3
         end if
   10    continue
      else
         do i = opacity_table%opal95_index_t+1,2,-1
            if (log10_temperature.gt.opacity_table%opal95_grid_logt(i)) then
               opacity_table%opal95_index_t = i - 1
               exit
            endif
         end do
         if (i < (2)) then
         opacity_table%opal95_index_t = 1
         end if
   20    continue
      endif
!     T INTERPOLATION FACTORS
      do i = 1,4
         interp_nodes(i) = opacity_table%opal95_grid_logt(opacity_table%opal95_index_t+i-1)
      end do
!     GET INTERPOLATION FACTORS IN T
      call interp(interp_nodes, weight, dweight, log10_temperature)
      do i = 1,4
         opacity_table%opal95_weight_t(i) = weight(i)
         opacity_table%opal95_dweight_t(i) = dweight(i)
      end do
!     GET INDICES IN RHO FOR EACH OF THE 4 T VALUES.
      if (opacity_table%opal95_logr.ge.opacity_table%opal95_grid_logr(opacity_table%opal95_index_rho(1)+2)) then
         do i = opacity_table%opal95_index_rho(1)+3,num_d-1
            if (opacity_table%opal95_logr.lt.opacity_table%opal95_grid_logr(i)) then
               opacity_table%opal95_index_rho(1) = i - 2
               exit
            endif
         end do
         if (i > (num_d-1)) then
         opacity_table%opal95_index_rho(1) = num_d - 3
         end if
   30    continue
      else
         do i = opacity_table%opal95_index_rho(1)+1,2,-1
            if (opacity_table%opal95_logr.gt.opacity_table%opal95_grid_logr(i)) then
               opacity_table%opal95_index_rho(1) = i - 1
               exit
            endif
         end do
         if (i < (2)) then
         opacity_table%opal95_index_rho(1) = 1
         end if
   40    continue
      endif
!     NOW CHECK IF WITHIN PORTION OF TABLE WITH DATA.
!     NOTE THAT A POINT CAN BE WITHIN THE TABLE AND STILL HAVE THE
!     4X4 MATRIX OF T,RHO AROUND IT OUTSIDE THE TABLE.
!     CHECK LOW RHO, LOW T REGIME; SINCE THE EMPTY REGION DEPENDS ON Z AT
!     LOW X, POSTPONE UNTIL YOU KNOW IF IT IS A 2D, 3D, OR 4D INTERPOLATION.
      if (hydrogen_fraction.lt.0.2d0 .and. opacity_table%opal95_index_t .le. 5) then
         low_regime_flag = .true.
         goto 60
      else
         low_regime_flag = .false.
      endif
!     CHECK HIGH RHO, HIGH T REGIME(BECAUSE OF THE TABLE GEOMETRY, THE
!     POINT WILL FALL OUTSIDE THE TABLE FIRST AT THE HIGHEST T).
      density_shifted = .false.
      opacity_table%opal95_extrap_hi = .false.
      if (opacity_table%opal95_fixed_z_opacity(opacity_table%opal95_index_x(1,1),opacity_table%opal95_index_t+3, &
           opacity_table%opal95_index_rho(1)+3) .ge. 9.9d0) then
         density_shifted = .true.
!        NOW THERE ARE TWO POSSIBLE SOLUTIONS: WE CAN USE A DIFFERENT SET
!        OF DENSITIES FOR EACH TEMPERATURE AND STAY WITHIN THE TABLE, OR
!        THE DESIRED DENSITY COULD BE OUTSIDE THE TABLE AT ONE OR MORE
!        TEMPERATURES.  IF THE DESIRED DENSITY IS OUTSIDE THE TABLE,
!        LINEAR EXTRAPOLATION IN R AT FIXED T IS USED.
         do k = 1,4
            if (opacity_table%opal95_fixed_z_opacity(opacity_table%opal95_index_x(1,1),opacity_table%opal95_index_t+k-1, &
                 opacity_table%opal95_index_rho(1)+1).ge.9.9d0) then
               opacity_table%opal95_extrap_hi = .true.
               opacity_table%opal95_extrap_hi_row(k) = .true.
            else
               opacity_table%opal95_extrap_hi_row(k) = .false.
            endif
         end do
!        FIND 4 NEAREST ELEMENTS WITHIN THE TABLE AT EACH OF THE 4 TEMPERATURES.
         x_table_index = opacity_table%opal95_index_x(1,1)
         density_base_index = opacity_table%opal95_index_rho(1)
         do i = 1,4
            temp_row_index = opacity_table%opal95_index_t+i-1
            if (opacity_table%opal95_fixed_z_opacity(x_table_index,temp_row_index, &
                 density_base_index+3).gt.9.9d0) then
               do j = density_base_index+2,1,-1
                  if (opacity_table%opal95_full_opacity(x_table_index,temp_row_index,j).le.9.9d0) then
                     opacity_table%opal95_index_rho(i) = j - 3
                     density_base_index = j - 3
                     exit
                  endif
               end do
   50          continue
            else
               opacity_table%opal95_index_rho(i) = density_base_index
            endif
         end do
      endif
   60 continue
!     GET INTERPOLATION FACTORS IN RHO FOR EACH OF THE 4 T VALUES.
!     ASSIGN ALL THE SAME VALUES IF GRID IN R IS SAME FOR ALL T.
      if (.not.density_shifted) then
         do i = 1,4
            interp_nodes(i) = opacity_table%opal95_grid_logr(opacity_table%opal95_index_rho(1)+i-1)
         end do
         if (.not.opacity_table%opal95_extrap_lo) then
            call interp(interp_nodes, weight, dweight, opacity_table%opal95_logr)
         else
            opacity_table%opal95_logr_lo_edge = opacity_table%opal95_grid_logr(1)
            call interp(interp_nodes, weight, dweight, opacity_table%opal95_logr_lo_edge)
         endif
         do j = 1,4
            opacity_table%opal95_index_rho(j) = opacity_table%opal95_index_rho(1)
            do i = 1,4
               opacity_table%opal95_weight_rho(j,i) = weight(i)
               opacity_table%opal95_dweight_rho(j,i) = dweight(i)
            end do
         end do
      else
         do j = 1,4
            do i = 1,4
               interp_nodes(i) = opacity_table%opal95_grid_logr(opacity_table%opal95_index_rho(j)+i-1)
            end do
!           IF EXTRAPOLATION IS BEING PERFORMED, ASSIGN R=REDGE
!           AND EXTRAPOLATE LINEARLY LATER
            if (.not.opacity_table%opal95_extrap_hi_row(j)) then
               call interp(interp_nodes, weight, dweight, opacity_table%opal95_logr)
            else
               extrap_logr = opacity_table%opal95_grid_logr(opacity_table%opal95_index_rho(j)+3)
               call interp(interp_nodes, weight, dweight, extrap_logr)
               opacity_table%opal95_logr_hi_edge(j) = extrap_logr
            endif
            do i = 1,4
               opacity_table%opal95_weight_rho(j,i) = weight(i)
               opacity_table%opal95_dweight_rho(j,i) = dweight(i)
            end do
         end do
      endif
!     DETERMINE WHETHER A 2D (RHO,T); 3D (X,RHO,T); OR 4D (Z,X,RHO,T)
!     INTERPOLATION IS NEEDED TO GET THE OPACITY.
!     JVS 04/11 force 4d interpolation for acoustic depth calculations
      if (.not. (compute_acoustic_depth .and. acoustic_depth_output)) then
      if (abs(metal_fraction-opacity_table%opal95_fixed_z)/max(opacity_table%opal95_fixed_z,1.0d-6).le.1.0d-4) then
         if (abs(hydrogen_fraction-opacity_table%opal95_surface_x).le.1.0d-4) then
!           2D INTERPOLATION IN SURFACE TABLE
            call op952d(opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
            continue
            return
         else
!           3D INTERPOLATION IN FIXED Z TABLE (X,T,RHO)
!           GET INTERPOLATION FACTORS IN X.
         if (hydrogen_fraction.gt.opacity_table%opal95_grid_x(opacity_table%opal95_index_x(1,1)+2)) then
            do i = opacity_table%opal95_index_x(1,1)+3,num_x-1
               if (hydrogen_fraction.lt.opacity_table%opal95_grid_x(i)) then
                  opacity_table%opal95_index_x(1,1) = i - 2
                  goto 70
               endif
            end do
            opacity_table%opal95_index_x(1,1) = num_x - 3
   70       continue
         else
            do i = opacity_table%opal95_index_x(1,1)+1,2,-1
               if (hydrogen_fraction.ge.opacity_table%opal95_grid_x(i)) then
                  opacity_table%opal95_index_x(1,1) = i - 1
                  goto 80
               endif
            end do
            opacity_table%opal95_index_x(1,1) = 1
   80       continue
         endif
!        CHECK FOR ABSENT TABLES AT HIGH X.
         if (opacity_table%opal95_fixed_z.ge.0.04d0) then
            if (hydrogen_fraction.ge.0.8d0) then
               if (opacity_table%opal95_fixed_z.lt.0.1d0) then
!                 NO TABLE 9 (X=0.95)
                  opacity_table%opal95_index_x(1,4) = num_x
                  opacity_table%opal95_index_x(1,3) = num_x - 2
                  opacity_table%opal95_index_x(1,2) = num_x - 3
                  opacity_table%opal95_index_x(1,1) = num_x - 4
               else
!                 NO TABLE 9 OR TABLE 10
                  opacity_table%opal95_index_x(1,4) = num_x - 2
                  opacity_table%opal95_index_x(1,3) = num_x - 3
                  opacity_table%opal95_index_x(1,2) = num_x - 4
                  opacity_table%opal95_index_x(1,1) = num_x - 5
               endif
            else
               x_shift_base = opacity_table%opal95_index_x(1,1) - 1
               do i = 1,4
                  opacity_table%opal95_index_x(1,i) = x_shift_base + i
               end do
            endif
         else
            x_shift_base = opacity_table%opal95_index_x(1,1) - 1
            do i = 1,4
               opacity_table%opal95_index_x(1,i) = x_shift_base + i
            end do
         endif
!        X INTERPOLATION FACTORS
         do i = 1,4
            x_table_index = opacity_table%opal95_index_x(1,i)
            interp_nodes(i) = opacity_table%opal95_grid_x(x_table_index)
         end do
         if (opacity_table%opal95_index_x(1,4).eq.num_x) then
            interp_nodes(4) = 1.0d0 - opacity_table%opal95_fixed_z
         endif
!        GET INTERPOLATION FACTORS IN X
         call intrp2(interp_nodes, weight, hydrogen_fraction)
         do j = 1,4
            opacity_table%opal95_weight_x(1,j) = weight(j)
         end do
         call op953d(opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
         continue
         return
      endif
      endif
      end if
  153 continue
!     4D INTERPOLATION IN Z,X,T,RHO
!     GET NEAREST TABLES IN Z.
      if (metal_fraction.gt.opacity_table%opal95_grid_z(opacity_table%opal95_index_z+2)) then
!mhp 7/12 corrected typo
!         DO I = JZ+3,NUMT-1
         do i = opacity_table%opal95_index_z+3,num_z-1
            if (metal_fraction.lt.opacity_table%opal95_grid_z(i)) then
               opacity_table%opal95_index_z = i - 2
               exit
            endif
         end do
         if (i > (num_z-1)) then
         opacity_table%opal95_index_z = num_z - 3
         end if
  110    continue
      else
         do i = opacity_table%opal95_index_z+1,2,-1
            if (metal_fraction.ge.opacity_table%opal95_grid_z(i)) then
               opacity_table%opal95_index_z = i - 1
               exit
            endif
         end do
         if (i < (2)) then
         opacity_table%opal95_index_z = 1
         end if
  120    continue
      endif
!     mhp 7/12 added failsafe
      opacity_table%opal95_index_z = min(opacity_table%opal95_index_z,num_z-3)
!     Z INTERPOLATION FACTORS
      do i = 1,4
         interp_nodes(i) = opacity_table%opal95_grid_z(opacity_table%opal95_index_z+i-1)
      end do
!     INTERPOLATION FACTORS IN Z
      call intrp2(interp_nodes, weight, metal_fraction)
      do i = 1,4
         opacity_table%opal95_weight_z(i) = weight(i)
      end do
!     GET NEAREST TABLES IN X AT EACH VALUE OF Z
      do k = 1,4
         z_at_table = opacity_table%opal95_grid_z(opacity_table%opal95_index_z+k-1)
         x_table_index = opacity_table%opal95_index_x(k,1)
         if (hydrogen_fraction.gt.opacity_table%opal95_grid_x(x_table_index+2)) then
            do i = opacity_table%opal95_index_x(k,1)+3,num_x-1
               if (hydrogen_fraction.lt.opacity_table%opal95_grid_x(i)) then
                  opacity_table%opal95_index_x(k,1) = i - 2
                  goto 75
               endif
            end do
            opacity_table%opal95_index_x(k,1) = num_x - 3
   75       continue
         else
            do i = opacity_table%opal95_index_x(k,1)+1,2,-1
               if (hydrogen_fraction.ge.opacity_table%opal95_grid_x(i)) then
                  opacity_table%opal95_index_x(k,1) = i - 1
                  goto 85
               endif
            end do
            opacity_table%opal95_index_x(k,1) = 1
   85       continue
         endif
!        CHECK FOR ABSENT TABLES AT HIGH X.
         if (z_at_table.ge.0.04d0) then
            if (hydrogen_fraction.ge.0.8d0) then
               if (z_at_table.lt.0.1d0) then
!                 NO TABLE 9 (X=0.95)
                  opacity_table%opal95_index_x(k,4) = num_x
                  opacity_table%opal95_index_x(k,3) = num_x - 2
                  opacity_table%opal95_index_x(k,2) = num_x - 3
                  opacity_table%opal95_index_x(k,1) = num_x - 4
               else
!                 NO TABLE 9 OR TABLE 10
                  opacity_table%opal95_index_x(k,4) = num_x - 2
                  opacity_table%opal95_index_x(k,3) = num_x - 3
                  opacity_table%opal95_index_x(k,2) = num_x - 4
                  opacity_table%opal95_index_x(k,1) = num_x - 5
               endif
            else
               x_shift_base = opacity_table%opal95_index_x(k,1) - 1
               do i = 1,4
                  opacity_table%opal95_index_x(k,i) = x_shift_base + i
               end do
            endif
         else
            x_shift_base = opacity_table%opal95_index_x(k,1) - 1
            do i = 1,4
               opacity_table%opal95_index_x(k,i) = x_shift_base + i
            end do
         endif
!        X INTERPOLATION FACTORS
         do i = 1,4
            x_table_index = opacity_table%opal95_index_x(k,i)
            interp_nodes(i) = opacity_table%opal95_grid_x(x_table_index)
         end do
         if (opacity_table%opal95_index_x(k,4).eq.num_x) then
            interp_nodes(4) = 1.0d0 - z_at_table
         endif
!        GET INTERPOLATION FACTORS IN X
         call intrp2(interp_nodes, weight, hydrogen_fraction)
         do j = 1,4
            opacity_table%opal95_weight_x(k,j) = weight(j)
         end do
      end do
      call op954d(opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
 9999 continue
      return
end subroutine getopal95
