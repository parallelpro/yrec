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
! Companion routines in this directory:
!      ll95tbl               read in all tables
!      opal95_fixed_z_table  generate table at fixed Z
!      opal95_surface_table  generate table at fixed X,Z
!      opal95_interp2d/3d/4d 2D (T,RHO) / 3D (X,T,RHO) / 4D (Z,X,T,RHO)
!
! 7/98 MHP DRIVER ROUTINE FOR OPACITY GIVEN RHO,T, X, AND Z.
! Determines whether a 2D (rho,T), 3D (X,rho,T), or 4D (Z,X,rho,T)
! interpolation is required and dispatches to opal95_interp2d/3d/4d
! accordingly.
subroutine getopal95(log10_density, log10_temperature, hydrogen_fraction, &
     metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ierr)
      use opacity_table_lib
      use numerics_lib
      implicit none

      double precision, intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      double precision, intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      integer, intent(out) :: ierr

      double precision :: interp_nodes(4), weight(4), dweight(4)
      integer :: i, k, j
      logical :: low_regime_flag, density_shifted
      integer :: x_table_index, density_base_index, temp_row_index
      double precision :: extrap_logr, z_at_table

!     ENSURE THAT WE ARE WITHIN THE OPAL 95 TABLES.
!     COMPUTE LOG R = RHO/T6**3
      ierr = 0
      opacity_table%opal95_logr = log10_density - 3.0d0*(log10_temperature-6.0d0)
!     CHECK T
      if (log10_temperature.lt.opacity_table%opal95_grid_logt(1) .or. &
           log10_temperature.gt.opacity_table%opal95_grid_logt(n_opal95_t)) then
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
!     (No low-R range check: the region below the first table element
!     is handled by linear extrapolation in R from the last two opacity
!     values, opal95_extrap_lo above.)
!     GET INTERPOLATION FACTORS IN T
      call stencil4_locate_opal95(opacity_table%opal95_grid_logt, n_opal95_t, log10_temperature, opacity_table%opal95_index_t)
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
      call stencil4_locate_opal95(opacity_table%opal95_grid_logr, n_opal95_d, opacity_table%opal95_logr, opacity_table%opal95_index_rho(1))
!     NOW CHECK IF WITHIN PORTION OF TABLE WITH DATA.
!     NOTE THAT A POINT CAN BE WITHIN THE TABLE AND STILL HAVE THE
!     4X4 MATRIX OF T,RHO AROUND IT OUTSIDE THE TABLE.
!     CHECK LOW RHO, LOW T REGIME; SINCE THE EMPTY REGION DEPENDS ON Z AT
!     LOW X, POSTPONE UNTIL YOU KNOW IF IT IS A 2D, 3D, OR 4D INTERPOLATION.
      if (hydrogen_fraction.lt.0.2d0 .and. opacity_table%opal95_index_t .le. 5) then
         low_regime_flag = .true.
      else
         low_regime_flag = .false.
!     CHECK HIGH RHO, HIGH T REGIME(BECAUSE OF THE TABLE GEOMETRY, THE
!     POINT WILL FALL OUTSIDE THE TABLE FIRST AT THE HIGHEST T).
      density_shifted = .false.
      opacity_table%opal95_extrap_hi = .false.
      if (opacity_table%opal95_fixed_z_opacity(opacity_table%opal95_index_x(1,1),opacity_table%opal95_index_t+3, &
           opacity_table%opal95_index_rho(1)+3) .ge. opal95_missing_test) then
         density_shifted = .true.
!        NOW THERE ARE TWO POSSIBLE SOLUTIONS: WE CAN USE A DIFFERENT SET
!        OF DENSITIES FOR EACH TEMPERATURE AND STAY WITHIN THE TABLE, OR
!        THE DESIRED DENSITY COULD BE OUTSIDE THE TABLE AT ONE OR MORE
!        TEMPERATURES.  IF THE DESIRED DENSITY IS OUTSIDE THE TABLE,
!        LINEAR EXTRAPOLATION IN R AT FIXED T IS USED.
         do k = 1,4
            if (opacity_table%opal95_fixed_z_opacity(opacity_table%opal95_index_x(1,1),opacity_table%opal95_index_t+k-1, &
                 opacity_table%opal95_index_rho(1)+1).ge.opal95_missing_test) then
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
                 density_base_index+3).gt.opal95_missing_test) then
               do j = density_base_index+2,1,-1
                  if (opacity_table%opal95_full_opacity(x_table_index,temp_row_index,j).le.opal95_missing_test) then
                     opacity_table%opal95_index_rho(i) = j - 3
                     density_base_index = j - 3
                     exit
                  endif
               end do
            else
               opacity_table%opal95_index_rho(i) = density_base_index
            endif
         end do
      endif
      endif
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
!     (The retired acoustic-depth mode used to force 4D here.)
      if (abs(metal_fraction-opacity_table%opal95_fixed_z)/max(opacity_table%opal95_fixed_z,1.0d-6).le.opal95_composition_tol) then
         if (abs(hydrogen_fraction-opacity_table%opal95_surface_x).le.opal95_composition_tol) then
!           2D INTERPOLATION IN SURFACE TABLE
            call opal95_interp2d(opacity_table%opal95_index_t, opacity_table%opal95_index_rho, &
     opacity_table%opal95_weight_t, opacity_table%opal95_dweight_t, &
     opacity_table%opal95_weight_rho, opacity_table%opal95_dweight_rho, &
     opacity_table%opal95_logr, opacity_table%opal95_logr_lo_edge, &
     opacity_table%opal95_logr_hi_edge, opacity_table%opal95_extrap_lo, &
     opacity_table%opal95_extrap_hi, opacity_table%opal95_extrap_hi_row, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
            return
         else
!           3D INTERPOLATION IN FIXED Z TABLE (X,T,RHO)
!           GET INTERPOLATION FACTORS IN X.
         call opal95_x_stencil(1, opacity_table%opal95_fixed_z, hydrogen_fraction)
         call opal95_interp3d(opacity_table%opal95_index_x(1,:), &
     opacity_table%opal95_weight_x(1,:), opacity_table%opal95_index_t, opacity_table%opal95_index_rho, &
     opacity_table%opal95_weight_t, opacity_table%opal95_dweight_t, &
     opacity_table%opal95_weight_rho, opacity_table%opal95_dweight_rho, &
     opacity_table%opal95_logr, opacity_table%opal95_logr_lo_edge, &
     opacity_table%opal95_logr_hi_edge, opacity_table%opal95_extrap_lo, &
     opacity_table%opal95_extrap_hi, opacity_table%opal95_extrap_hi_row, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
         return
      endif
      endif
!     4D INTERPOLATION IN Z,X,T,RHO
!     GET NEAREST TABLES IN Z.
      if (metal_fraction.gt.opacity_table%opal95_grid_z(opacity_table%opal95_index_z+2)) then
!mhp 7/12 corrected typo
!         DO I = JZ+3,NUMT-1
         do i = opacity_table%opal95_index_z+3,n_opal95_z-1
            if (metal_fraction.lt.opacity_table%opal95_grid_z(i)) then
               opacity_table%opal95_index_z = i - 2
               exit
            endif
         end do
         if (i > (n_opal95_z-1)) then
         opacity_table%opal95_index_z = n_opal95_z - 3
         end if
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
      endif
!     mhp 7/12 added failsafe
      opacity_table%opal95_index_z = min(opacity_table%opal95_index_z,n_opal95_z-3)
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
         call opal95_x_stencil(k, z_at_table, hydrogen_fraction)
      end do
      call opal95_interp4d(opacity_table%opal95_index_z, opacity_table%opal95_weight_z, &
     opacity_table%opal95_index_x, opacity_table%opal95_weight_x, &
     opacity_table%opal95_index_t, opacity_table%opal95_index_rho, &
     opacity_table%opal95_weight_t, opacity_table%opal95_dweight_t, &
     opacity_table%opal95_weight_rho, opacity_table%opal95_dweight_rho, &
     opacity_table%opal95_logr, opacity_table%opal95_logr_lo_edge, &
     opacity_table%opal95_logr_hi_edge, opacity_table%opal95_extrap_lo, &
     opacity_table%opal95_extrap_hi, opacity_table%opal95_extrap_hi_row, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
      return

contains

! opal95_x_stencil
! For Z row k of the stencil (Z value z_here), pick the four X tables
! opacity_table%opal95_index_x(k,1:4) bracketing hydrogen_fraction and
! store their weights in opacity_table%opal95_weight_x(k,1:4). The
! warm-start search moves opal95_index_x(k,1) from its previous value;
! at Z >= 0.04 and X >= 0.8 the missing high-X tables (no X=0.95 for
! Z >= 0.04, no X=1-Z for Z >= 0.1) force the stencil to the top
! available tables. 2026 readability: extracted token-for-token from
! the 3D (k = 1, Z = opal95_fixed_z) and 4D (k = 1..4, Z of each row)
! branches above.
subroutine opal95_x_stencil(k, z_here, hydrogen_fraction)
      integer, intent(in) :: k
      double precision, intent(in) :: z_here, hydrogen_fraction
      integer :: i, j, x_table_index, x_shift_base
      double precision :: interp_nodes(4), weight(4)

      if (hydrogen_fraction.gt.opacity_table%opal95_grid_x(opacity_table%opal95_index_x(k,1)+2)) then
         do i = opacity_table%opal95_index_x(k,1)+3,n_opal95_x-1
            if (hydrogen_fraction.lt.opacity_table%opal95_grid_x(i)) then
               opacity_table%opal95_index_x(k,1) = i - 2
               exit
            endif
         end do
         if (i .gt. n_opal95_x-1) opacity_table%opal95_index_x(k,1) = n_opal95_x - 3
      else
         do i = opacity_table%opal95_index_x(k,1)+1,2,-1
            if (hydrogen_fraction.ge.opacity_table%opal95_grid_x(i)) then
               opacity_table%opal95_index_x(k,1) = i - 1
               exit
            endif
         end do
         if (i .lt. 2) opacity_table%opal95_index_x(k,1) = 1
      endif
!     CHECK FOR ABSENT TABLES AT HIGH X.
      if (z_here.ge.0.04d0) then
         if (hydrogen_fraction.ge.0.8d0) then
            if (z_here.lt.0.1d0) then
!              NO TABLE 9 (X=0.95)
               opacity_table%opal95_index_x(k,4) = n_opal95_x
               opacity_table%opal95_index_x(k,3) = n_opal95_x - 2
               opacity_table%opal95_index_x(k,2) = n_opal95_x - 3
               opacity_table%opal95_index_x(k,1) = n_opal95_x - 4
            else
!              NO TABLE 9 OR TABLE 10
               opacity_table%opal95_index_x(k,4) = n_opal95_x - 2
               opacity_table%opal95_index_x(k,3) = n_opal95_x - 3
               opacity_table%opal95_index_x(k,2) = n_opal95_x - 4
               opacity_table%opal95_index_x(k,1) = n_opal95_x - 5
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
!     X INTERPOLATION FACTORS
      do i = 1,4
         x_table_index = opacity_table%opal95_index_x(k,i)
         interp_nodes(i) = opacity_table%opal95_grid_x(x_table_index)
      end do
      if (opacity_table%opal95_index_x(k,4).eq.n_opal95_x) then
         interp_nodes(4) = 1.0d0 - z_here
      endif
!     GET INTERPOLATION FACTORS IN X
      call intrp2(interp_nodes, weight, hydrogen_fraction)
      do j = 1,4
         opacity_table%opal95_weight_x(k,j) = weight(j)
      end do
end subroutine opal95_x_stencil

end subroutine getopal95
