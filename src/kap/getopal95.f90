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
     metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)

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

! FULL SET OF TABLES: OPACITY AS A FUNCTION OF Z AND X, T, RHO/T6**3
! TABLES ARE INCREMENTED IN SETS OF NZ*NX.  SO THE TABLES FOR THE
! THIRD METAL ABUNDANCE (3 X 10**-4)BEGIN AT TABLE 21 AND END AT TABLE 30.
! FOR THE HIGH VALUES OF Z, THE NUMBER OF X TABLES IS NOT THE SAME (I.E.
! X<0.9 FOR Z=0.1).
! FOR EACH COMPOSITION A FULL GRID IN (T,RHO/T6**3) IS RETAINED.
      double precision :: opal95_grid_logt(num_t), opal95_grid_x(num_x), &
           opal95_grid_logr(num_d), opal95_grid_z(num_z), &
           opal95_full_opacity(num_xz,num_t,num_d)
      integer :: opal95_num_x_at_z(num_z), opal95_table_start_index(num_z)
      common /llot95a/ opal95_grid_logt, opal95_grid_x, opal95_grid_logr, &
           opal95_grid_z, opal95_full_opacity, opal95_num_x_at_z, &
           opal95_table_start_index
! AS ABOVE FOR THE MODEL Z.
      double precision :: opal95_fixed_z_opacity(num_x,num_t,num_d), opal95_fixed_z
      common /llot95/opal95_fixed_z_opacity, opal95_fixed_z
! AS ABOVE FOR DESIRED SURFACE VALUE OF X.
      double precision :: opal95_surface_opacity(num_t,num_d), opal95_surface_x
      common /llot95e/opal95_surface_opacity, opal95_surface_x
! INDICES FOR INTERPOLATION IN Z,X,T, AND R
      integer :: opal95_index_z, opal95_index_x(4,4), opal95_index_t, &
           opal95_index_rho(4)
      common/op95indx/ opal95_index_z, opal95_index_x, opal95_index_t, &
           opal95_index_rho
! INTERPOLATION FACTORS FOR Z,X,T, AND R, AS WELL AS DERIVATIVE
! FACTORS FOR T AND RHO.
      double precision :: opal95_weight_z(4), opal95_weight_x(4,4), &
           opal95_weight_t(4), opal95_dweight_t(4), opal95_weight_rho(4,4), &
           opal95_dweight_rho(4,4)
      common/op95fact/ opal95_weight_z, opal95_weight_x, opal95_weight_t, &
           opal95_dweight_t, opal95_weight_rho, opal95_dweight_rho
      double precision :: opal95_logr, opal95_logr_lo_edge, opal95_logr_hi_edge(4)
      logical :: opal95_extrap_lo, opal95_extrap_hi, opal95_extrap_hi_row(4)
      common/op95ext/ opal95_logr, opal95_logr_lo_edge, opal95_logr_hi_edge, &
           opal95_extrap_lo, opal95_extrap_hi, opal95_extrap_hi_row
! JVS Acoustic depth common block
! KC 2025-05-30 reordered common block elements
!       COMMON/ACDPTH/TAUCZN,DELADJ(JSON),TAUHE, TNORM, TCZ, WHE, ICLCD,
! common/acdpth/: only compute_acoustic_depth (LADON) and
! acoustic_depth_output (LACOUT) are used here; the remaining members
! are unused placeholders preserving the shared storage layout, kept
! under their (lowercased) original cryptic names since their exact
! meaning is not confidently known from this file alone.
      double precision :: taucz_placeholder, deladj_placeholder(json), &
           tauhe_placeholder, tnorm_placeholder, tcz_placeholder, &
           whe_placeholder
      double precision :: acatmr_placeholder(json), acatmd_placeholder(json), &
           acatmp_placeholder(json), acatmt_placeholder(json), tatmos_placeholder
      double precision :: ageout_placeholder(5)
      logical :: lclcd_placeholder
      integer :: iclcd_placeholder, iacat_placeholder, ijlast_placeholder
      logical :: ljlast_placeholder, ljwrt_placeholder, compute_acoustic_depth, &
           laoly_placeholder
      integer :: ijvs_placeholder, ijent_placeholder, ijdel_placeholder
      logical :: acoustic_depth_output
      common/acdpth/taucz_placeholder,deladj_placeholder,tauhe_placeholder, &
           tnorm_placeholder, tcz_placeholder, whe_placeholder, &
           acatmr_placeholder, acatmd_placeholder, acatmp_placeholder, &
           acatmt_placeholder, tatmos_placeholder, &
           ageout_placeholder, lclcd_placeholder, iclcd_placeholder, &
           iacat_placeholder, ijlast_placeholder, ljlast_placeholder, &
           ljwrt_placeholder, compute_acoustic_depth, laoly_placeholder, &
           ijvs_placeholder, ijent_placeholder, ijdel_placeholder, &
           acoustic_depth_output
! END JVS
      double precision :: interp_nodes(4), weight(4), dweight(4)
      save

      integer :: i, k, j
      logical :: low_regime_flag, density_shifted
      integer :: x_table_index, density_base_index, temp_row_index, x_shift_base
      double precision :: extrap_logr, z_at_table

!     ENSURE THAT WE ARE WITHIN THE OPAL 95 TABLES.
!     COMPUTE LOG R = RHO/T6**3
      opal95_logr = log10_density - 3.0d0*(log10_temperature-6.0d0)
!     CHECK T
      if (log10_temperature.lt.opal95_grid_logt(1) .or. &
           log10_temperature.gt.opal95_grid_logt(num_t)) then
         write(*,5) log10_temperature
    5    format(' LOG T OF',f11.6,'OUT OF OPAL 95 TABLE RANGE')
         stop
      endif
!     CHECK TO SEE IF EXTRAPOLATION BELOW THE FIRST TABLE ELEMENT
!     IN DENSITY IS NEEDED.
      if (opal95_logr.lt.opal95_grid_logr(1)) then
         opal95_extrap_lo = .true.
      else
         opal95_extrap_lo = .false.
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
      if (log10_temperature.ge.opal95_grid_logt(opal95_index_t+2)) then
         do i = opal95_index_t+3,num_t-1
            if (log10_temperature.lt.opal95_grid_logt(i)) then
               opal95_index_t = i - 2
               goto 10
            endif
         end do
         opal95_index_t = num_t - 3
   10    continue
      else
         do i = opal95_index_t+1,2,-1
            if (log10_temperature.gt.opal95_grid_logt(i)) then
               opal95_index_t = i - 1
               goto 20
            endif
         end do
         opal95_index_t = 1
   20    continue
      endif
!     T INTERPOLATION FACTORS
      do i = 1,4
         interp_nodes(i) = opal95_grid_logt(opal95_index_t+i-1)
      end do
!     GET INTERPOLATION FACTORS IN T
      call interp(interp_nodes, weight, dweight, log10_temperature)
      do i = 1,4
         opal95_weight_t(i) = weight(i)
         opal95_dweight_t(i) = dweight(i)
      end do
!     GET INDICES IN RHO FOR EACH OF THE 4 T VALUES.
      if (opal95_logr.ge.opal95_grid_logr(opal95_index_rho(1)+2)) then
         do i = opal95_index_rho(1)+3,num_d-1
            if (opal95_logr.lt.opal95_grid_logr(i)) then
               opal95_index_rho(1) = i - 2
               goto 30
            endif
         end do
         opal95_index_rho(1) = num_d - 3
   30    continue
      else
         do i = opal95_index_rho(1)+1,2,-1
            if (opal95_logr.gt.opal95_grid_logr(i)) then
               opal95_index_rho(1) = i - 1
               goto 40
            endif
         end do
         opal95_index_rho(1) = 1
   40    continue
      endif
!     NOW CHECK IF WITHIN PORTION OF TABLE WITH DATA.
!     NOTE THAT A POINT CAN BE WITHIN THE TABLE AND STILL HAVE THE
!     4X4 MATRIX OF T,RHO AROUND IT OUTSIDE THE TABLE.
!     CHECK LOW RHO, LOW T REGIME; SINCE THE EMPTY REGION DEPENDS ON Z AT
!     LOW X, POSTPONE UNTIL YOU KNOW IF IT IS A 2D, 3D, OR 4D INTERPOLATION.
      if (hydrogen_fraction.lt.0.2d0 .and. opal95_index_t .le. 5) then
         low_regime_flag = .true.
         goto 60
      else
         low_regime_flag = .false.
      endif
!     CHECK HIGH RHO, HIGH T REGIME(BECAUSE OF THE TABLE GEOMETRY, THE
!     POINT WILL FALL OUTSIDE THE TABLE FIRST AT THE HIGHEST T).
      density_shifted = .false.
      opal95_extrap_hi = .false.
      if (opal95_fixed_z_opacity(opal95_index_x(1,1),opal95_index_t+3, &
           opal95_index_rho(1)+3) .ge. 9.9d0) then
         density_shifted = .true.
!        NOW THERE ARE TWO POSSIBLE SOLUTIONS: WE CAN USE A DIFFERENT SET
!        OF DENSITIES FOR EACH TEMPERATURE AND STAY WITHIN THE TABLE, OR
!        THE DESIRED DENSITY COULD BE OUTSIDE THE TABLE AT ONE OR MORE
!        TEMPERATURES.  IF THE DESIRED DENSITY IS OUTSIDE THE TABLE,
!        LINEAR EXTRAPOLATION IN R AT FIXED T IS USED.
         do k = 1,4
            if (opal95_fixed_z_opacity(opal95_index_x(1,1),opal95_index_t+k-1, &
                 opal95_index_rho(1)+1).ge.9.9d0) then
               opal95_extrap_hi = .true.
               opal95_extrap_hi_row(k) = .true.
            else
               opal95_extrap_hi_row(k) = .false.
            endif
         end do
!        FIND 4 NEAREST ELEMENTS WITHIN THE TABLE AT EACH OF THE 4 TEMPERATURES.
         x_table_index = opal95_index_x(1,1)
         density_base_index = opal95_index_rho(1)
         do i = 1,4
            temp_row_index = opal95_index_t+i-1
            if (opal95_fixed_z_opacity(x_table_index,temp_row_index, &
                 density_base_index+3).gt.9.9d0) then
               do j = density_base_index+2,1,-1
                  if (opal95_full_opacity(x_table_index,temp_row_index,j).le.9.9d0) then
                     opal95_index_rho(i) = j - 3
                     density_base_index = j - 3
                     goto 50
                  endif
               end do
   50          continue
            else
               opal95_index_rho(i) = density_base_index
            endif
         end do
      endif
   60 continue
!     GET INTERPOLATION FACTORS IN RHO FOR EACH OF THE 4 T VALUES.
!     ASSIGN ALL THE SAME VALUES IF GRID IN R IS SAME FOR ALL T.
      if (.not.density_shifted) then
         do i = 1,4
            interp_nodes(i) = opal95_grid_logr(opal95_index_rho(1)+i-1)
         end do
         if (.not.opal95_extrap_lo) then
            call interp(interp_nodes, weight, dweight, opal95_logr)
         else
            opal95_logr_lo_edge = opal95_grid_logr(1)
            call interp(interp_nodes, weight, dweight, opal95_logr_lo_edge)
         endif
         do j = 1,4
            opal95_index_rho(j) = opal95_index_rho(1)
            do i = 1,4
               opal95_weight_rho(j,i) = weight(i)
               opal95_dweight_rho(j,i) = dweight(i)
            end do
         end do
      else
         do j = 1,4
            do i = 1,4
               interp_nodes(i) = opal95_grid_logr(opal95_index_rho(j)+i-1)
            end do
!           IF EXTRAPOLATION IS BEING PERFORMED, ASSIGN R=REDGE
!           AND EXTRAPOLATE LINEARLY LATER
            if (.not.opal95_extrap_hi_row(j)) then
               call interp(interp_nodes, weight, dweight, opal95_logr)
            else
               extrap_logr = opal95_grid_logr(opal95_index_rho(j)+3)
               call interp(interp_nodes, weight, dweight, extrap_logr)
               opal95_logr_hi_edge(j) = extrap_logr
            endif
            do i = 1,4
               opal95_weight_rho(j,i) = weight(i)
               opal95_dweight_rho(j,i) = dweight(i)
            end do
         end do
      endif
!     DETERMINE WHETHER A 2D (RHO,T); 3D (X,RHO,T); OR 4D (Z,X,RHO,T)
!     INTERPOLATION IS NEEDED TO GET THE OPACITY.
!     JVS 04/11 force 4d interpolation for acoustic depth calculations
      if (compute_acoustic_depth .and. acoustic_depth_output) goto 153  ! If we're worrying about the acoustic depth, default to 4d interp
      if (abs(metal_fraction-opal95_fixed_z)/max(opal95_fixed_z,1.0d-6).le.1.0d-4) then
         if (abs(hydrogen_fraction-opal95_surface_x).le.1.0d-4) then
!           2D INTERPOLATION IN SURFACE TABLE
            call op952d(opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
            goto 9999
         else
!           3D INTERPOLATION IN FIXED Z TABLE (X,T,RHO)
!           GET INTERPOLATION FACTORS IN X.
         if (hydrogen_fraction.gt.opal95_grid_x(opal95_index_x(1,1)+2)) then
            do i = opal95_index_x(1,1)+3,num_x-1
               if (hydrogen_fraction.lt.opal95_grid_x(i)) then
                  opal95_index_x(1,1) = i - 2
                  goto 70
               endif
            end do
            opal95_index_x(1,1) = num_x - 3
   70       continue
         else
            do i = opal95_index_x(1,1)+1,2,-1
               if (hydrogen_fraction.ge.opal95_grid_x(i)) then
                  opal95_index_x(1,1) = i - 1
                  goto 80
               endif
            end do
            opal95_index_x(1,1) = 1
   80       continue
         endif
!        CHECK FOR ABSENT TABLES AT HIGH X.
         if (opal95_fixed_z.ge.0.04d0) then
            if (hydrogen_fraction.ge.0.8d0) then
               if (opal95_fixed_z.lt.0.1d0) then
!                 NO TABLE 9 (X=0.95)
                  opal95_index_x(1,4) = num_x
                  opal95_index_x(1,3) = num_x - 2
                  opal95_index_x(1,2) = num_x - 3
                  opal95_index_x(1,1) = num_x - 4
               else
!                 NO TABLE 9 OR TABLE 10
                  opal95_index_x(1,4) = num_x - 2
                  opal95_index_x(1,3) = num_x - 3
                  opal95_index_x(1,2) = num_x - 4
                  opal95_index_x(1,1) = num_x - 5
               endif
            else
               x_shift_base = opal95_index_x(1,1) - 1
               do i = 1,4
                  opal95_index_x(1,i) = x_shift_base + i
               end do
            endif
         else
            x_shift_base = opal95_index_x(1,1) - 1
            do i = 1,4
               opal95_index_x(1,i) = x_shift_base + i
            end do
         endif
!        X INTERPOLATION FACTORS
         do i = 1,4
            x_table_index = opal95_index_x(1,i)
            interp_nodes(i) = opal95_grid_x(x_table_index)
         end do
         if (opal95_index_x(1,4).eq.num_x) then
            interp_nodes(4) = 1.0d0 - opal95_fixed_z
         endif
!        GET INTERPOLATION FACTORS IN X
         call intrp2(interp_nodes, weight, hydrogen_fraction)
         do j = 1,4
            opal95_weight_x(1,j) = weight(j)
         end do
         call op953d(opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
         goto 9999
      endif
      endif
  153 continue
!     4D INTERPOLATION IN Z,X,T,RHO
!     GET NEAREST TABLES IN Z.
      if (metal_fraction.gt.opal95_grid_z(opal95_index_z+2)) then
!mhp 7/12 corrected typo
!         DO I = JZ+3,NUMT-1
         do i = opal95_index_z+3,num_z-1
            if (metal_fraction.lt.opal95_grid_z(i)) then
               opal95_index_z = i - 2
               goto 110
            endif
         end do
         opal95_index_z = num_z - 3
  110    continue
      else
         do i = opal95_index_z+1,2,-1
            if (metal_fraction.ge.opal95_grid_z(i)) then
               opal95_index_z = i - 1
               goto 120
            endif
         end do
         opal95_index_z = 1
  120    continue
      endif
!     mhp 7/12 added failsafe
      opal95_index_z = min(opal95_index_z,num_z-3)
!     Z INTERPOLATION FACTORS
      do i = 1,4
         interp_nodes(i) = opal95_grid_z(opal95_index_z+i-1)
      end do
!     INTERPOLATION FACTORS IN Z
      call intrp2(interp_nodes, weight, metal_fraction)
      do i = 1,4
         opal95_weight_z(i) = weight(i)
      end do
!     GET NEAREST TABLES IN X AT EACH VALUE OF Z
      do k = 1,4
         z_at_table = opal95_grid_z(opal95_index_z+k-1)
         x_table_index = opal95_index_x(k,1)
         if (hydrogen_fraction.gt.opal95_grid_x(x_table_index+2)) then
            do i = opal95_index_x(k,1)+3,num_x-1
               if (hydrogen_fraction.lt.opal95_grid_x(i)) then
                  opal95_index_x(k,1) = i - 2
                  goto 75
               endif
            end do
            opal95_index_x(k,1) = num_x - 3
   75       continue
         else
            do i = opal95_index_x(k,1)+1,2,-1
               if (hydrogen_fraction.ge.opal95_grid_x(i)) then
                  opal95_index_x(k,1) = i - 1
                  goto 85
               endif
            end do
            opal95_index_x(k,1) = 1
   85       continue
         endif
!        CHECK FOR ABSENT TABLES AT HIGH X.
         if (z_at_table.ge.0.04d0) then
            if (hydrogen_fraction.ge.0.8d0) then
               if (z_at_table.lt.0.1d0) then
!                 NO TABLE 9 (X=0.95)
                  opal95_index_x(k,4) = num_x
                  opal95_index_x(k,3) = num_x - 2
                  opal95_index_x(k,2) = num_x - 3
                  opal95_index_x(k,1) = num_x - 4
               else
!                 NO TABLE 9 OR TABLE 10
                  opal95_index_x(k,4) = num_x - 2
                  opal95_index_x(k,3) = num_x - 3
                  opal95_index_x(k,2) = num_x - 4
                  opal95_index_x(k,1) = num_x - 5
               endif
            else
               x_shift_base = opal95_index_x(k,1) - 1
               do i = 1,4
                  opal95_index_x(k,i) = x_shift_base + i
               end do
            endif
         else
            x_shift_base = opal95_index_x(k,1) - 1
            do i = 1,4
               opal95_index_x(k,i) = x_shift_base + i
            end do
         endif
!        X INTERPOLATION FACTORS
         do i = 1,4
            x_table_index = opal95_index_x(k,i)
            interp_nodes(i) = opal95_grid_x(x_table_index)
         end do
         if (opal95_index_x(k,4).eq.num_x) then
            interp_nodes(4) = 1.0d0 - z_at_table
         endif
!        GET INTERPOLATION FACTORS IN X
         call intrp2(interp_nodes, weight, hydrogen_fraction)
         do j = 1,4
            opal95_weight_x(k,j) = weight(j)
         end do
      end do
      call op954d(opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
 9999 continue
      return
end subroutine getopal95
