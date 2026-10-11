!----------------------------------------------------------------------
! alex06_opacity
!----------------------------------------------------------------------
! The low-temperature molecular opacities of Ferguson et al. (2005), the
! 'Alexander 2006' tables: one file of tables in X and Z, each in log T
! and log R.
!     read_alex06_tables       read the file
!     compute_alex06_opacity   the opacity and its log derivatives, from
!                              a table at the zone's (X, Z) built by
!                              build_alex06_composition_table and
!                              reused while the composition holds
! The tables live in opacity_tables (opacity_table%alex06_*).
module alex06_opacity
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: alex_composition_tol, n_alex06_d, n_alex06_t, &
           n_alex06_x, n_alex06_z, opacity_table, stencil4_locate
      use numerics_lib, only: interp, intrp2, lagrange4
      use physical_constants, only: ln10
      use math_lib
      implicit none
      private
      public :: read_alex06_tables, compute_alex06_opacity

contains

!----------------------------------------------------------------------
! read_alex06_tables (formerly readalex06)
!----------------------------------------------------------------------
! MHP 3/09  SUBROUTINE FOR READING IN ALEXANDER LOW T OPACITIES.
!
! ALEXANDER OPACITIES FERGUSON ET AL. 2005  APJ,623,585
! FOR EACH  X = 0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.8, 0.9, 1-Z
! INCREASING ORDER OF Z = 0.0, 0.00001, 0.00003, 0.0001, 0.0003,
!     (16)               0.001, 0.002, 0.004, 0.01, 0.02, 0.03,
!                        0.04, 0.05, 0.06, 0.08, 0.10
! DESCENDING ORDER OF LOGT 4.5-3.5 IN 0.05 DEX INCREMENTS
!     (85)                 3.5-2.9 IN 0.01 DEX INCREMENTS
!                          2.9-2.7 IN 0.05 DEX INCREMENTS
! ASCENDING ORDER OF  LOGR -8.0 -7.5 -7.0 -6.5 -6.0 -5.5 -5.0 -4.5 -4.0 -3.5
!     (19)                -3.0 -2.5 -2.0 -1.5 -1.0 -0.5  0.0  0.5  1.0
! WHERE R = RHO/T_6**3
! READ IN ALL Z VALUES AT FIXED X.  THEN GENERATE SURFACE X,Z TABLE.
! INTERPOLATION BETWEEN COMPOSITIONS IS 4-POINT LAGRANGIAN.
!
! MHP 8/25 Added file name to subroutine call
subroutine read_alex06_tables(alex06_table_path, ierr)
      integer, intent(out) :: ierr

! the table file's unit
      integer :: table_unit
      character(len=256), intent(in) :: alex06_table_path

!     LOCAL ARRAYS
      real(dp) :: row_logr_check(n_alex06_d), header_x, header_z, row_temp
! alex06_grid_x/alex06_grid_z/alex06_grid_logt/alex06_grid_logr/
! alex06_cached_x/alex06_cached_z/alex06_index_x/alex06_index_t/
! alex06_index_r defaults moved to opacity_tables.f90: DATA can no
! longer target them here now that they're use-associated.
      integer :: i, ii, j, jj, k, kk

      ierr = 0
      open(newunit=table_unit,file=alex06_table_path)
!     READ IN INITIAL X AND Z; ENSURE THAT THEY HAVE THE EXPECTED VALUES.
      do i = 1,n_alex06_x-1
         do ii = 1, n_alex06_z
!           INDEX FOR STORING TABLES : 15 SETS OF METAL ABUNDANCE WITH 7 SETS
!           OF HYDROGEN FOR EACH.  THESE ARE STORED IN THE ORDER (X,Z) OF
!           (0,0),(0,0.1),...(0.9,0),(0.9,0.1).  The final tables are defined
!           at X = 1-Z except for the 0.1 Z case (already read in at X=0.9).
!           HEADER INFORMATION: X AND Z; CHECK FOR CONSISTENCY
            read(table_unit,10) header_x, header_z
   10       format(30x,f9.6,7x,f9.6//)
            if (header_x /= opacity_table%alex06_grid_x(i) .or. header_z /= opacity_table%alex06_grid_z(ii)) then
               write(*,15) opacity_table%alex06_grid_x(i), header_x, opacity_table%alex06_grid_z(ii), header_z
   15          format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                    1x,'EXPECTED AND ACTUAL X,Z',4f7.2,' RUN STOPPED')
!               STOP
            endif
            read(table_unit,20) (row_logr_check(k),k=1,n_alex06_d)
   20       format(6x,19f7.3)
            do kk=1,16
               if (row_logr_check(kk) /= opacity_table%alex06_grid_logr(kk)) then
                  write(*,25) kk, opacity_table%alex06_grid_logr(kk), row_logr_check(kk)
   25             format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                       1x,'EXPECTED AND ACTUAL R',i3,2f7.3,' RUN STOPPED')
                  ! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
                  ierr = 1
                  return
               endif
            end do
            jj = (i-1)*n_alex06_z+ii
!           OPACITY INFORMATION AT EACH SHELL: CHECK FOR CONSISTENCY WITH T.
!           STORE IN A NUMXZ*NUMT*NUMR ARRAY.
            do j = n_alex06_t,1,-1
               read(table_unit,30) row_temp, (opacity_table%alex06_full_opacity(jj,j,k),k=1,n_alex06_d)
   30          format(f5.3,1x,19f7.3)
               if (row_temp /= opacity_table%alex06_grid_logt(j)) then
                  write(*,35) j, opacity_table%alex06_grid_logt(j), row_temp
   35             format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                       1x,'EXPECTED AND ACTUAL T',i3,2f7.3,' RUN STOPPED')
                  ! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
                  ierr = 1
                  return
               endif
            end do
         end do
      end do
!     FINAL SET OF X TABLES ARE DEFINED AT X = 1 - Z AND THERE IS NOT ONE FOR
!     Z = 0.1 (ALREADY READ IN AT X=0.9).
      do ii = 1, n_alex06_z-1
!        INDEX FOR STORING TABLES : 15 SETS OF METAL ABUNDANCE WITH 7 SETS
!        OF HYDROGEN FOR EACH.  THESE ARE STORED IN THE ORDER (X,Z) OF
!        (0,0),(0,0.1),...(0.9,0),(0.9,0.1).  The final tables are defined
!        at X = 1-Z except for the 0.1 Z case (already read in at X=0.9).
!        HEADER INFORMATION: X AND Z; CHECK FOR CONSISTENCY
         read(table_unit,10) header_x, header_z
         if (header_x /= (1.0d0 - opacity_table%alex06_grid_z(ii)) .or. &
              header_z /= opacity_table%alex06_grid_z(ii)) then
            write(*,15) 1.0d0 - opacity_table%alex06_grid_z(ii), header_x, opacity_table%alex06_grid_z(ii), header_z
            ! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
            ierr = 1
            return
         endif
         read(table_unit,20) (row_logr_check(k),k=1,n_alex06_d)
         do kk=1,16
            if (row_logr_check(kk) /= opacity_table%alex06_grid_logr(kk)) then
               write(*,25) kk, opacity_table%alex06_grid_logr(kk), row_logr_check(kk)
               ! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
               ierr = 1
               return
            endif
         end do
         jj = (n_alex06_x-1)*n_alex06_z+ii
!        OPACITY INFORMATION AT EACH SHELL: CHECK FOR CONSISTENCY WITH T.
!        STORE IN A NUMXZ*NUMT*NUMR ARRAY.
         do j = n_alex06_t,1,-1
            read(table_unit,30) row_temp, (opacity_table%alex06_full_opacity(jj,j,k),k=1,n_alex06_d)
            if (row_temp /= opacity_table%alex06_grid_logt(j)) then
               write(*,35) j, opacity_table%alex06_grid_logt(j), row_temp
               ! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
               ierr = 1
               return
            endif
         end do
      end do
      close(unit=table_unit)
!     INITIALIZE FIXED Z,X TABLE
      do i = 1,n_alex06_t
         do j = 1,n_alex06_d
            opacity_table%alex06_opacity(i,j)=opacity_table%alex06_full_opacity(1,i,j)
         end do
      end do
end subroutine read_alex06_tables

!----------------------------------------------------------------------
! build_alex06_composition_table (formerly alex06tab)
!----------------------------------------------------------------------
! Generates the fixed-(X,Z) table for the Alexander 2006 low-
! temperature opacities by 4-point Lagrangian interpolation, first in
! Z then in X, from the full set of tables read by read_alex06_tables.f90.
! The target X, Z are taken from opacity_table%alex06_cached_x/
! alex06_cached_z, which the caller (compute_alex06_opacity.f90) sets before
! calling this routine.
subroutine build_alex06_composition_table(ierr)

      integer, intent(out) :: ierr

      real(dp) :: interp_nodes(4), weight_z(4), weight_x(4)
      real(dp) :: opacity_by_x(4,n_alex06_t,n_alex06_d)
      real(dp) :: x_max, z_max, interp_target
      integer :: i, j, iz, k, kk, kk2, kk3, kk4

!     XE = DESIRED X; ZE = DESIRED Z
      ierr = 0
      x_max = 1.0d0 - opacity_table%alex06_cached_z
!     CHECK THAT THE REQUESTED COMPOSITION IS INSIDE TABLE BOUNDS
      if (opacity_table%alex06_cached_x < 0.0d0 .or. opacity_table%alex06_cached_x > x_max) then
         write(*,5) opacity_table%alex06_cached_x, opacity_table%alex06_cached_z
    5    format('ILLEGAL COMPOSITION (X,Z) = ',2f6.2,' IN ALEX06.RUN STOPPED')
! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib's kap_get).
         ierr = 1
         return
      endif
!     PERMIT EXTRAPOLATION IN Z BY UP TO 1 TABLE ELEMENT
      z_max = opacity_table%alex06_grid_z(n_alex06_z)+(opacity_table%alex06_grid_z(n_alex06_z)-opacity_table%alex06_grid_z(n_alex06_z-1))
      if (opacity_table%alex06_cached_z < 0.0d0 .or. opacity_table%alex06_cached_z > z_max) then
         write(*,5) opacity_table%alex06_cached_x, opacity_table%alex06_cached_z
! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib's kap_get).
         ierr = 1
         return
      endif
!     FIND 4 NEAREST TABLES IN Z.
      do i = 3,n_alex06_z-2
         if (opacity_table%alex06_cached_z <= opacity_table%alex06_grid_z(i)) then
            iz = i - 2
            exit
         endif
      end do
      if (i > (n_alex06_z-2)) then
      iz = n_alex06_z - 3
      end if
!     FIND 4 NEAREST TABLES IN X.
      do i = 3,n_alex06_x-2
         if (opacity_table%alex06_cached_x <= opacity_table%alex06_grid_x(i)) then
            opacity_table%alex06_index_x = i - 2
            exit
         endif
      end do
      if (i > (n_alex06_x-2)) then
!     NO TABLE FOR X > 0.9 IF Z =0.10 OR MORE
      if (opacity_table%alex06_cached_z >= 0.1d0) then
         opacity_table%alex06_index_x = n_alex06_x - 4
      else
         opacity_table%alex06_index_x = n_alex06_x - 3
      endif
      end if
!     INTERPOLATION FACTORS FOR Z
      do i = 1,4
         interp_nodes(i) = opacity_table%alex06_grid_z(iz+i-1)
      end do
      interp_target = opacity_table%alex06_cached_z
      call intrp2(interp_nodes, weight_z, interp_target)
!     THE DIFFERENCE IN THE NUMBER OF TABLES FOR THE Z=0.1 CASE REQUIRES SOME
!     CARE IN X INTERPOLATION.  FIRST 3 X CASES CAN BE TREATED NORMALLY.
!     (2026 R3: the four Z tables kk, kk+1, kk+2, kk+3 are consecutive
!     here, so the sum is lagrange4 over the contiguous slice.)
      do k = 1,3
         kk = n_alex06_z*(opacity_table%alex06_index_x+k-2)+iz
         do i = 1,n_alex06_t
            do j = 1,n_alex06_d
               opacity_by_x(k,i,j) = lagrange4(weight_z, opacity_table%alex06_full_opacity(kk:kk+3,i,j))
            end do
         end do
      end do
!     IF IN THE HIGH Z AND HIGH X DOMAIN THE TOP TABLE IS X = 1-Z (ENTRY NUMX) EXCEPT
!     FOR THE Z=0.10 CASE (WHERE THE X=0.9 CASE DOUBLES AS THE X=1-Z CASE).
      if (opacity_table%alex06_index_x == n_alex06_x-3.and.iz == n_alex06_z-3) then
!        USE DIFFERENT INDEXING FOR THE LAST TABLE
         kk = n_alex06_z*(opacity_table%alex06_index_x+2)+iz
         kk2 = kk + 1
         kk3 = kk2 + 1
         kk4 = n_alex06_z*(opacity_table%alex06_index_x+2)
      else
         kk = n_alex06_z*(opacity_table%alex06_index_x+2)+iz
         kk2 = kk + 1
         kk3 = kk2 + 1
         kk4 = kk3 + 1
      endif
!     (Inline rather than lagrange4: kk4 is not kk+3 in the Z=0.10 case.)
      do i = 1,n_alex06_t
         do j = 1,n_alex06_d
            opacity_by_x(4,i,j) = weight_z(1)*opacity_table%alex06_full_opacity(kk,i,j)+ &
                 weight_z(2)*opacity_table%alex06_full_opacity(kk2,i,j) + &
                 weight_z(3)*opacity_table%alex06_full_opacity(kk3,i,j) + &
                 weight_z(4)*opacity_table%alex06_full_opacity(kk4,i,j)
         end do
      end do
!     NOW DO X INTERPOLATION
!     INTERPOLATION FACTORS FOR X
      do i = 1,3
         interp_nodes(i) = opacity_table%alex06_grid_x(opacity_table%alex06_index_x+i-1)
      end do
      if (iz == n_alex06_z-3) then
         interp_nodes(4) = 1.0d0-opacity_table%alex06_cached_z
      else
         interp_nodes(4) = opacity_table%alex06_grid_x(opacity_table%alex06_index_x+3)
      endif
      interp_target = opacity_table%alex06_cached_x
      call intrp2(interp_nodes, weight_x, interp_target)
      do i = 1,n_alex06_t
         do j = 1,n_alex06_d
            opacity_table%alex06_opacity(i,j) = lagrange4(weight_x, opacity_by_x(1:4,i,j))
         end do
      end do
end subroutine build_alex06_composition_table

!----------------------------------------------------------------------
! compute_alex06_opacity (formerly getalex06)
!----------------------------------------------------------------------
! MHP 3/09 Computes Alexander 2006 (Ferguson et al. 2005) low-
! temperature opacities by 4-point Lagrangian interpolation in log T
! and log R = rho/T6**3, reloading the fixed-(X,Z) table via
! build_alex06_composition_table.f90 whenever the requested composition changes.
subroutine compute_alex06_opacity(log10_density, log10_temperature, hydrogen_fraction, &
     metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ierr)


      real(dp), intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      integer, intent(out) :: ierr

      real(dp) :: interp_nodes(4), weight_t(4), dweight_t(4), &
           weight_r(4), dweight_r(4), opacity_row(4), dlnkap_dlnr_row(4)
      real(dp) :: delta_z, delta_x, logr, saved_r
      logical :: extrapolate_linear
      integer :: i, ii

      ierr = 0
      delta_z = abs(metal_fraction-opacity_table%alex06_cached_z)
      delta_x = abs(hydrogen_fraction-opacity_table%alex06_cached_x)
!     ENSURE THAT OPACITY TABLE HAS THE SAME X,Z VALUE AS THE POINT
      if (delta_z > alex_composition_tol .or. delta_x > alex_composition_tol) then
         opacity_table%alex06_cached_z = metal_fraction
         opacity_table%alex06_cached_x = hydrogen_fraction
         call build_alex06_composition_table(ierr)
         if (ierr /= 0) return
      endif
!     COMPUTE R FOR GRID POINT
      logr = log10_density - 3.0d0*(log10_temperature-6.0d0)
!     Insure that index IT is within the required array bounds
      if (opacity_table%alex06_index_t < 1) opacity_table%alex06_index_t=1
      if ((opacity_table%alex06_index_t+2) > n_alex06_t) opacity_table%alex06_index_t=n_alex06_t-2
!     FIND NEAREST GRID POINTS IN T.
      call stencil4_locate(opacity_table%alex06_grid_logt, n_alex06_t, log10_temperature, opacity_table%alex06_index_t)
!     INTERPOLATION FACTORS IN LOG T
      do i = 1,4
         interp_nodes(i) = opacity_table%alex06_grid_logt(opacity_table%alex06_index_t+i-1)
      end do
      call interp(interp_nodes, weight_t, dweight_t, log10_temperature)
!     Insure that index ID is within the required array bounds
      if (opacity_table%alex06_index_r < 1) opacity_table%alex06_index_r=1
      if ((opacity_table%alex06_index_r+2) > n_alex06_d) opacity_table%alex06_index_r=n_alex06_d-2
!     FIND NEAREST GRID POINTS IN R = RHO/T6**3
      call stencil4_locate(opacity_table%alex06_grid_logr, n_alex06_d, logr, opacity_table%alex06_index_r)
!     INTERPOLATION FACTORS IN LOG R
      if (logr > opacity_table%alex06_grid_logr(n_alex06_d)) then
         extrapolate_linear = .true.
         saved_r = logr
         logr = opacity_table%alex06_grid_logr(n_alex06_d)
      else
         extrapolate_linear = .false.
      endif
      do i = 1,4
         interp_nodes(i) = opacity_table%alex06_grid_logr(opacity_table%alex06_index_r+i-1)
      end do
      call interp(interp_nodes, weight_r, dweight_r, logr)
!     INTERPOLATE IN LOG R AT FIXED T
      do i = 1,4
         ii = opacity_table%alex06_index_t + i - 1
         opacity_row(i) = weight_r(1)*opacity_table%alex06_opacity(ii,opacity_table%alex06_index_r) + &
              weight_r(2)*opacity_table%alex06_opacity(ii,opacity_table%alex06_index_r+1) + &
              weight_r(3)*opacity_table%alex06_opacity(ii,opacity_table%alex06_index_r+2) + &
              weight_r(4)*opacity_table%alex06_opacity(ii,opacity_table%alex06_index_r+3)
         dlnkap_dlnr_row(i) = dweight_r(1)*opacity_table%alex06_opacity(ii,opacity_table%alex06_index_r) + &
              dweight_r(2)*opacity_table%alex06_opacity(ii,opacity_table%alex06_index_r+1) + &
              dweight_r(3)*opacity_table%alex06_opacity(ii,opacity_table%alex06_index_r+2) + &
              dweight_r(4)*opacity_table%alex06_opacity(ii,opacity_table%alex06_index_r+3)
      end do
      if (extrapolate_linear) then
         do i = 1,4
            ii = opacity_table%alex06_index_t+i-1
            opacity_row(i) = opacity_row(i)+(saved_r-logr)* &
                 (opacity_table%alex06_opacity(ii,n_alex06_d)-opacity_table%alex06_opacity(ii,n_alex06_d-1))/ &
                 (opacity_table%alex06_grid_logr(n_alex06_d)-opacity_table%alex06_grid_logr(n_alex06_d-1))
         end do
         logr = saved_r
      endif
!     INTERPOLATE IN T
      log10_opacity = lagrange4(weight_t, opacity_row)
!     D LN CAPPA/D LN T AT FIXED R
      dlnkap_dlnt = lagrange4(dweight_t, opacity_row)
!     INTERPOLATE IN D LN CAPPA/ D LN R AT FIXED T
      dlnkap_dlnrho = lagrange4(weight_t, dlnkap_dlnr_row)
!     CORRECT FROM DERIVATIVE AT FIXED R TO DERIVATIVE AT FIXED RHO
      dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
      opacity = exp(ln10*log10_opacity)
end subroutine compute_alex06_opacity

end module alex06_opacity
