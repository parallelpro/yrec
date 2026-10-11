!----------------------------------------------------------------------
! opal95_opacity
!----------------------------------------------------------------------
! The atomic opacities in the OPAL 1995 table layout (the OPAL and OP
! tables YREC ships, e.g. GS98.OP17): 126 tables in X and Z, each in
! log T and log R = log rho - 3 log T6.
!     read_opal95_tables          read a table file, then build the
!                                 fixed-Z table at the run's Z
!     build_opal95_surface_table  the fixed-X, fixed-Z table at the
!                                 envelope composition (refreshed when
!                                 the envelope X changes)
!     compute_opal95_opacity      the opacity and its log derivatives:
!                                 on the surface table near the
!                                 envelope composition, else by 3D
!                                 (X, T, R) or 4D (Z, X, T, R)
!                                 Lagrange interpolation
! The tables and the interpolation state live in opacity_tables
! (opacity_table%opal95_*).
module opal95_opacity
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables
      use numerics_lib, only: interp, intrp2, lagrange4
      use math_lib
      implicit none
      private
      public :: read_opal95_tables, build_opal95_surface_table, compute_opal95_opacity

contains

!----------------------------------------------------------------------
! read_opal95_tables (formerly ll95tbl)
!----------------------------------------------------------------------
! MHP 7/98 MODIFIED TO READ IN ALL METAL ABUNDACES FOR OPACITY TABLES.
! Reads the full set of OPAL95 opacity tables (all tabulated Z and X),
! skipping the header, then builds the fixed-Z table via
! build_opal95_fixed_z_table. Any premature end of file (fewer than n_opal95_tables
! tables, or a truncated table body) is reported and returns ierr = 1.
subroutine read_opal95_tables(opal95_table_path, ierr)
      use io_units
      integer, intent(out) :: ierr

! the table file's unit
      integer :: table_unit
      character(len=256), intent(in) :: opal95_table_path

! opal95_table_start_index/opal95_num_x_at_z/opal95_grid_x/
! opal95_grid_z/opal95_index_z/opal95_index_x/opal95_index_t/
! opal95_index_rho defaults live in opacity_tables.f90.
! LOCAL VECTOR, USED TO SKIP LONG HEADER.
      character(len=132) :: header_line

      integer :: fmt_start, iz, ix, nn, n, i, j, nmax
      integer :: read_status
      real(dp) :: xx, zz, xxt, target_z

!     READ IN OPACITY TABLES, SKIPPING HEADER.  HEY, IT WORKS.
      ierr = 0
      open(newunit=table_unit,file=opal95_table_path,status='OLD',access='SEQUENTIAL')
      fmt_start=1

      do
         read(table_unit,1,end=9999) header_line
      if (header_line(fmt_start:fmt_start+4).eq.'TABLE') exit
      end do
    1 format(a)
      read(header_line,'(36X,F7.4,11X,F7.4)') xx,zz
      iz = 1
      ix = 1
      nn = 1
!     ENTRY POINT FOR GETTING NEW TABLES.
      table_loop: do
      n = opacity_table%opal95_table_start_index(iz)+ix
      if (ix.lt.n_opal95_x) then
         xxt = opacity_table%opal95_grid_x(ix)
      else
         xxt = 1.0d0 - opacity_table%opal95_grid_z(iz)
      endif
      if (opacity_table%opal95_grid_z(iz).ne.zz .or. xxt.ne.xx) then
       write(run_log_unit,*)' OPAL95: Z ERROR INCOMPATIBLE TABLE'
         ! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
         ierr = 1
         return
      endif

!     READ IN HEADER INFO: GRID IN RHO/T6**3
      read(table_unit,20,iostat=read_status) (opacity_table%opal95_grid_logr(i),i=1,n_opal95_d)
      if (read_status .lt. 0) then
         write(*,*) 'read_opal95_tables: OPAL95 opacity table file ended after ', nn-1, &
              ' of ', n_opal95_tables, ' tables'
         close(table_unit)
         ierr = 1
         return
      end if
      if (read_status .gt. 0) then
         write(*,*) 'read_opal95_tables: malformed OPAL95 opacity table (header read failed)'
         ierr = 1
         return
      end if
   20 format(///,4x,19f7.1,/)
!     READ IN FULL TABLE: LOG CAPPA AS A FUNCTION OF LOG T AND
!     LOG R = RHO/T6**3
      do i = 1,57
         read(table_unit,30,end=9999) opacity_table%opal95_grid_logt(i), &
              (opacity_table%opal95_full_opacity(n,i,j),j=1,n_opal95_d)
      end do
   30 format(f4.2,19f7.3)
! MHP 12/97 NOW TREAT CORNER WITHOUT DATA.
      read(table_unit,31,end=9999) opacity_table%opal95_grid_logt(58), &
           (opacity_table%opal95_full_opacity(n,58,j),j=1,n_opal95_d-1)
   31 format(f4.2,18f7.3)

      opacity_table%opal95_full_opacity(n,58,n_opal95_d) = opal95_missing_opacity
      do i = 59,60
         read(table_unit,32,end=9999) opacity_table%opal95_grid_logt(i), &
              (opacity_table%opal95_full_opacity(n,i,j),j=1,n_opal95_d-2)
   32    format(f4.2,17f7.3)
         do j = n_opal95_d-1,n_opal95_d
            opacity_table%opal95_full_opacity(n,i,j) = opal95_missing_opacity
         end do
      end do
      do i = 61,64
         read(table_unit,33,end=9999) opacity_table%opal95_grid_logt(i), &
              (opacity_table%opal95_full_opacity(n,i,j),j=1,n_opal95_d-3)
   33    format(f4.2,16f7.3)
         do j = n_opal95_d-2,n_opal95_d
            opacity_table%opal95_full_opacity(n,i,j) = opal95_missing_opacity
         end do
      end do
      do i = 65,69
         read(table_unit,34,end=9999) opacity_table%opal95_grid_logt(i), &
              (opacity_table%opal95_full_opacity(n,i,j),j=1,n_opal95_d-4)
   34    format(f4.2,15f7.3)
         do j = n_opal95_d-3,n_opal95_d
            opacity_table%opal95_full_opacity(n,i,j) = opal95_missing_opacity
         end do
      end do
      i = n_opal95_t
      read(table_unit,35,end=9999) opacity_table%opal95_grid_logt(i), &
           (opacity_table%opal95_full_opacity(n,i,j),j=1,n_opal95_d-5)
   35 format(f4.2,14f7.3)
      do j = n_opal95_d-4,n_opal95_d
         opacity_table%opal95_full_opacity(n,i,j) = opal95_missing_opacity
      end do

!     EXIT IF CORRECT NUMBER OF TABLES READ IN.
      if (nn.lt.n_opal95_tables) then
      nn = nn + 1
!     NEED TO ACCOUNT FOR FEWER X VALUES AT HIGHER Z.
      if (ix.le.8) then
         nmax = n_opal95_z
      else if (ix.eq.9) then
         nmax = n_opal95_z - 3
      else
         nmax = n_opal95_z - 1
      endif
!     READ IN NEXT METAL ABUNDANCE (AT FIXED X) OR READ IN FIRST METAL
!     ABUNDANCE AT NEXT X.
      if (iz.lt.nmax) then
         iz = iz + 1
      else
         iz = 1
         ix = ix + 1
      endif

!     RETURN TO READ IN NEXT TABLE.
      read(table_unit,900) xx,zz
  900 format(/36x,f7.4,11x,f7.4)
      cycle table_loop
      end if
      exit table_loop
      end do table_loop

      close(table_unit)
!     NOW GENERATE A TABLE AT A FIXED VALUE OF Z.
!     NOTE THAT FOR METAL DIFFUSION A 4-D INTERPOLATION (IN X,Z,T,RHO)
!     IS PERFORMED RATHER THAN A LINEAR INTERPOLATION BETWEEN TWO FIXED
!     Z TABLES.
      target_z = kap_settings%opal95_single_table_z

      call build_opal95_fixed_z_table(target_z, ierr)
      return
!     PREMATURE END OF FILE INSIDE A TABLE BODY (end=9999 above).
 9999 write(*,*) 'read_opal95_tables: OPAL95 opacity table file ended inside table ', nn
      close(table_unit)
      ierr = 1
      return
end subroutine read_opal95_tables

!----------------------------------------------------------------------
! build_opal95_fixed_z_table (formerly opal95_fixed_z_table)
!----------------------------------------------------------------------
! MHP 7/98 COMPUTE OPACITY TABLE AT FIXED Z FROM THE FULL OP95 SET.
! Builds the fixed-Z (X, T, rho) opacity table by cubic interpolation
! across the 4 nearest Z tables in the full OPAL95 table set
! (opacity_table%opal95_full_opacity, filled by read_opal95_tables).
subroutine build_opal95_fixed_z_table(metal_fraction, ierr)

      integer, intent(out) :: ierr

      real(dp), intent(in) :: metal_fraction

      real(dp) :: table_z_nodes(4), z_weight(4), z_weight_hix(4)
      integer :: i, j, k, z_table_index, z_table_index2
      integer :: table1_index, table2_index, table3_index, table4_index

      ierr = 0
      opacity_table%opal95_fixed_z = metal_fraction
!  FIND 4 NEAREST TABLES IN Z TO DESIRED VALUE.
      do i = 3,n_opal95_z-1
         if (opacity_table%opal95_grid_z(i).ge.metal_fraction) then
            z_table_index = i - 2
            exit
         endif
      end do
      if (i > (n_opal95_z-1)) then
! DESIRED Z > 0.1D0; STOP
      write(*,5)metal_fraction
    5 format(1x,'DESIRED Z',f10.6,'OUTSIDE OP95 TABLE RANGE'/ &
           ' RUN STOPPED')
      ! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
      ierr = 1
      return
      end if
!  FIND INTERPOLATION FACTORS IN Z.
      do i = 1,4
         table_z_nodes(i) = opacity_table%opal95_grid_z(z_table_index+i-1)
      end do
      call intrp2(table_z_nodes, z_weight, metal_fraction)
!  FOR X TABLES 1-8, ALL VALUES OF Z PRESENT; PERFORM INTERPOLATION.
      do i = 1,8
! INDICES FOR 4 DESIRED COMPOSITIONS.
         table1_index = opacity_table%opal95_table_start_index(z_table_index)+i
         table2_index = opacity_table%opal95_table_start_index(z_table_index+1)+i
         table3_index = opacity_table%opal95_table_start_index(z_table_index+2)+i
         table4_index = opacity_table%opal95_table_start_index(z_table_index+3)+i
         do j = 1,n_opal95_t
            do k = 1,n_opal95_d
               opacity_table%opal95_fixed_z_opacity(i,j,k) = z_weight(1)*opacity_table%opal95_full_opacity(table1_index,j,k)+ &
               z_weight(2)*opacity_table%opal95_full_opacity(table2_index,j,k)+z_weight(3)*opacity_table%opal95_full_opacity(table3_index,j,k)+ &
               z_weight(4)*opacity_table%opal95_full_opacity(table4_index,j,k)
            end do
         end do
      end do
!  FOR X TABLES 9 AND 10, HIGH VALUES OF Z ARE NOT PRESENT
!  OMIT TABLE 9 (X=0.95) IF DESIRED Z > 0.04
      if (metal_fraction.lt.0.04d0) then
!  OTHERWISE, CHECK TO ENSURE THAT THE 4 Z TABLES USED HAVE Z < 0.04
!  ADJUST INTERPOLATION FACTORS IF NEEDED
      if (opacity_table%opal95_grid_z(z_table_index+3).gt.0.04d0) then
         z_table_index2 = n_opal95_z - 6
         do i = 1,4
            table_z_nodes(i) = opacity_table%opal95_grid_z(z_table_index2+i-1)
         end do
         call intrp2(table_z_nodes, z_weight_hix, metal_fraction)
      else
         z_table_index2 = z_table_index
         do i = 1,4
            z_weight_hix(i) = z_weight(i)
         end do
      endif
      table1_index = opacity_table%opal95_table_start_index(z_table_index2)+9
      table2_index = opacity_table%opal95_table_start_index(z_table_index2+1)+9
      table3_index = opacity_table%opal95_table_start_index(z_table_index2+2)+9
      table4_index = opacity_table%opal95_table_start_index(z_table_index2+3)+9
      do j = 1,n_opal95_t
         do k = 1,n_opal95_d
            opacity_table%opal95_fixed_z_opacity(9,j,k) = z_weight_hix(1)*opacity_table%opal95_full_opacity(table1_index,j,k)+ &
            z_weight_hix(2)*opacity_table%opal95_full_opacity(table2_index,j,k)+z_weight_hix(3)*opacity_table%opal95_full_opacity(table3_index,j,k)+ &
            z_weight_hix(4)*opacity_table%opal95_full_opacity(table4_index,j,k)
         end do
      end do
      end if
!  TABLE 10 (X = 1-Z) IS ALWAYS BUILT; ITS Z STENCIL MUST NOT INCLUDE
!  THE Z=0.1 TABLE, SO PICK z_table_index2/z_weight_hix AGAIN HERE.
!  (Z > grid_z(n_opal95_z-1) IS REJECTED ABOVE, SO THIS TEST IS ALWAYS
!  TRUE TODAY; WITHOUT IT THE TABLE-9 STENCIL WOULD BE REUSED.)
      if (metal_fraction.lt.0.1d0) then
!  CHECK TO ENSURE THAT Z=0.1 TABLE IS NOT ONE OF THE 4 TABLES; ADJUST
!  INTERPOLATION FACTORS IF NEEDED.
      if (opacity_table%opal95_grid_z(z_table_index+3).ge.0.1d0) then
         z_table_index2 = n_opal95_z - 4
         do i = 1,4
            table_z_nodes(i) = opacity_table%opal95_grid_z(z_table_index2+i-1)
         end do
         call intrp2(table_z_nodes, z_weight_hix, metal_fraction)
      else
         z_table_index2 = z_table_index
         do i = 1,4
            z_weight_hix(i) = z_weight(i)
         end do
      endif
      end if
      table1_index = opacity_table%opal95_table_start_index(z_table_index2)+10
      table2_index = opacity_table%opal95_table_start_index(z_table_index2+1)+10
      table3_index = opacity_table%opal95_table_start_index(z_table_index2+2)+10
      table4_index = opacity_table%opal95_table_start_index(z_table_index2+3)+10
      do j = 1,n_opal95_t
         do k = 1,n_opal95_d
            opacity_table%opal95_fixed_z_opacity(10,j,k) = z_weight_hix(1)*opacity_table%opal95_full_opacity(table1_index,j,k)+ &
            z_weight_hix(2)*opacity_table%opal95_full_opacity(table2_index,j,k)+z_weight_hix(3)*opacity_table%opal95_full_opacity(table3_index,j,k)+ &
            z_weight_hix(4)*opacity_table%opal95_full_opacity(table4_index,j,k)
         end do
      end do
      return
end subroutine build_opal95_fixed_z_table

!----------------------------------------------------------------------
! build_opal95_surface_table (formerly opal95_surface_table)
!----------------------------------------------------------------------
! MHP 7/98 COMPUTE OPACITY TABLE AT FIXED X FROM THE OP95 TABLE AT
! THE MODEL Z. Builds the surface (T, rho) opacity table for the
! given hydrogen abundance by cubic interpolation across the 4
! nearest X tables in the fixed-Z table set
! (opacity_table%opal95_fixed_z_opacity, filled by
! build_opal95_fixed_z_table).
subroutine build_opal95_surface_table(hydrogen_fraction)


      real(dp), intent(in) :: hydrogen_fraction

      real(dp) :: table_x_nodes(4), x_weight(4)
      integer :: i, j, k, x_table_index
      integer :: table1_index, table2_index, table3_index, table4_index
! x_table_index normally names the first of four consecutive X tables;
! this value means "tables 6, 7, 8 and 10" (table 9 absent).
      integer, parameter :: x_stencil_6_7_8_10 = 11

      opacity_table%opal95_surface_x = hydrogen_fraction
!  FIND 4 NEAREST TABLES IN X TO DESIRED VALUE.
      if (hydrogen_fraction.le.0.8d0 .or. &
           opacity_table%opal95_fixed_z.le.0.04d0) then
! ALL FOUR X TABLES PRESENT: EITHER X IS BELOW THE HIGH-X GAP, OR
! THE HIGH X TABLES ARE PRESENT AT THIS LOW Z.
! (2026 readability: the two branches were identical and were merged.)
         do i = 3,n_opal95_x-1
            if (opacity_table%opal95_grid_x(i).ge.hydrogen_fraction) then
               x_table_index = i - 2
               exit
            endif
         end do
         if (i > (n_opal95_x-1)) then
         x_table_index = n_opal95_x - 3
         end if
         do i = 1,4
            table_x_nodes(i) = opacity_table%opal95_grid_x(x_table_index+i-1)
         end do
      else if (opacity_table%opal95_fixed_z.ge.0.1d0) then
! USE TABLES 6-9
         x_table_index = 6
         do i = 1,4
            table_x_nodes(i) = opacity_table%opal95_grid_x(x_table_index+i-1)
         end do
      else
! IF Z IS BETWEEN 0.04 AND 0.1, TABLES 1-8 AND 10 EXIST.
! SINCE WE HAVE ALREADY DETERMINED THAT X > 0.8, USE 6-8 AND 10.
         x_table_index = x_stencil_6_7_8_10
         table_x_nodes(1) = opacity_table%opal95_grid_x(6)
         table_x_nodes(2) = opacity_table%opal95_grid_x(7)
         table_x_nodes(3) = opacity_table%opal95_grid_x(8)
         table_x_nodes(4) = opacity_table%opal95_grid_x(10)
      endif
!
!  FIND INTERPOLATION FACTORS IN X.
      call intrp2(table_x_nodes, x_weight, hydrogen_fraction)
! INDICES FOR 4 DESIRED COMPOSITIONS.
      if (x_table_index.ne.x_stencil_6_7_8_10) then
         table1_index = x_table_index
         table2_index = x_table_index + 1
         table3_index = x_table_index + 2
         table4_index = x_table_index + 3
      else
         table1_index = 6
         table2_index = 7
         table3_index = 8
         table4_index = 10
      endif
      do j = 1,n_opal95_t
         do k = 1,n_opal95_d
            opacity_table%opal95_surface_opacity(j,k) = x_weight(1)*opacity_table%opal95_fixed_z_opacity(table1_index,j,k) + &
            x_weight(2)*opacity_table%opal95_fixed_z_opacity(table2_index,j,k) + x_weight(3)*opacity_table%opal95_fixed_z_opacity(table3_index,j,k) + &
            x_weight(4)*opacity_table%opal95_fixed_z_opacity(table4_index,j,k)
         end do
      end do
      return
end subroutine build_opal95_surface_table

!----------------------------------------------------------------------
! compute_opal95_opacity (formerly getopal95)
!----------------------------------------------------------------------
! THIS IS THE INTERPOLATION FACILITY FOR THE LIVERMORE OPACITY TABLES
! USING CUBIC SPLINE INTERPLOATION SCHEME.   (FOR OPAL95)
! --------------------------------------------------------------
! Companion routines in this directory:
!      read_opal95_tables               read in all tables
!      build_opal95_fixed_z_table  generate table at fixed Z
!      build_opal95_surface_table  generate table at fixed X,Z
!      interpolate_opal95_2d/3d/4d 2D (T,RHO) / 3D (X,T,RHO) / 4D (Z,X,T,RHO)
!
! 7/98 MHP DRIVER ROUTINE FOR OPACITY GIVEN RHO,T, X, AND Z.
! Determines whether a 2D (rho,T), 3D (X,rho,T), or 4D (Z,X,rho,T)
! interpolation is required and dispatches to interpolate_opal95_2d/3d/4d
! accordingly.
subroutine compute_opal95_opacity(log10_density, log10_temperature, hydrogen_fraction, &
     metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ierr)

      real(dp), intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      integer, intent(out) :: ierr

      real(dp) :: interp_nodes(4), weight(4), dweight(4)
      integer :: i, k, j
      logical :: low_regime_flag, density_shifted
      integer :: x_table_index, density_base_index, temp_row_index
      real(dp) :: extrap_logr, z_at_table

!     ENSURE THAT WE ARE WITHIN THE OPAL 95 TABLES.
!     COMPUTE LOG R = RHO/T6**3
      ierr = 0
      opacity_table%opal95_logr = log10_density - 3.0d0*(log10_temperature-6.0d0)
!     CHECK T
      if (log10_temperature.lt.opacity_table%opal95_grid_logt(1) .or. &
           log10_temperature.gt.opacity_table%opal95_grid_logt(n_opal95_t)) then
         write(*,5) log10_temperature
    5    format(' LOG T OF',f11.6,'OUT OF OPAL 95 TABLE RANGE')
! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the facade
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
            call interpolate_opal95_2d(opacity_table%opal95_index_t, opacity_table%opal95_index_rho, &
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
         call locate_opal95_x_stencil(1, opacity_table%opal95_fixed_z, hydrogen_fraction)
         call interpolate_opal95_3d(opacity_table%opal95_index_x(1,:), &
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
         call locate_opal95_x_stencil(k, z_at_table, hydrogen_fraction)
      end do
      call interpolate_opal95_4d(opacity_table%opal95_index_z, opacity_table%opal95_weight_z, &
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

! locate_opal95_x_stencil
! For Z row k of the stencil (Z value z_here), pick the four X tables
! opacity_table%opal95_index_x(k,1:4) bracketing hydrogen_fraction and
! store their weights in opacity_table%opal95_weight_x(k,1:4). The
! warm-start search moves opal95_index_x(k,1) from its previous value;
! at Z >= 0.04 and X >= 0.8 the missing high-X tables (no X=0.95 for
! Z >= 0.04, no X=1-Z for Z >= 0.1) force the stencil to the top
! available tables. 2026 readability: extracted token-for-token from
! the 3D (k = 1, Z = opal95_fixed_z) and 4D (k = 1..4, Z of each row)
! branches above.
subroutine locate_opal95_x_stencil(k, z_here, hydrogen_fraction)
      integer, intent(in) :: k
      real(dp), intent(in) :: z_here, hydrogen_fraction
      integer :: i, j, x_table_index, x_shift_base
      real(dp) :: interp_nodes(4), weight(4)

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
end subroutine locate_opal95_x_stencil

end subroutine compute_opal95_opacity

!----------------------------------------------------------------------
! interpolate_opal95_2d (formerly opal95_interp2d)
!----------------------------------------------------------------------
! MHP 7/98 GET OPAL95 OPACITY FROM TABLE AT SURFACE X,Z
! 2D interpolation (T, rho) in the surface-X table produced by
! build_opal95_surface_table. Called from compute_opal95_opacity when both X and Z match the model
! surface composition exactly.
! 2026 wave 3 (R5): the T/rho stencil (indices, weights, log R
! extrapolation state) that compute_opal95_opacity computes per call used to reach
! this routine as hidden module state (opacity_table%opal95_index_*,
! opal95_weight_*, opal95_logr*, opal95_extrap_*); it now arrives as
! explicit intent(in) arguments, passed straight from those members by
! compute_opal95_opacity. The opacity tables themselves stay in opacity_table.
subroutine interpolate_opal95_2d(index_t, index_rho, weight_t, dweight_t, &
     weight_rho, dweight_rho, logr, logr_lo_edge, logr_hi_edge, &
     extrap_lo, extrap_hi, extrap_hi_row, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)


      integer, intent(in) :: index_t, index_rho(4)
      real(dp), intent(in) :: weight_t(4), dweight_t(4), &
           weight_rho(4,4), dweight_rho(4,4)
      real(dp), intent(in) :: logr, logr_lo_edge, logr_hi_edge(4)
      logical, intent(in) :: extrap_lo, extrap_hi, extrap_hi_row(4)
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

      real(dp) :: logcappa_at_t(4), dlogcappa_dlogr_at_t(4)
      integer :: i, temp_index, density_index
      real(dp) :: delta_logr

! FIND OPACITY AT EACH OF THE 4 NEARBY VALUES OF T.
      do i = 1,4
         temp_index = index_t + i - 1
! FOR EACH T, INTERPOLATE IN 4 VALUES OF R
         density_index = index_rho(i)
! LOG CAPPA AT DESIRED RHO FOR EACH OF THE 4 DESIRED T.
         logcappa_at_t(i) = weight_rho(i,1)*opacity_table%opal95_surface_opacity(temp_index,density_index) &
              + weight_rho(i,2)*opacity_table%opal95_surface_opacity(temp_index,density_index+1) &
              + weight_rho(i,3)*opacity_table%opal95_surface_opacity(temp_index,density_index+2) &
              + weight_rho(i,4)*opacity_table%opal95_surface_opacity(temp_index,density_index+3)
! D LOG CAPPA/D LOG R FOR EACH OF THE 4 DESIRED T.
         dlogcappa_dlogr_at_t(i) = dweight_rho(i,1)*opacity_table%opal95_surface_opacity(temp_index,density_index) &
              + dweight_rho(i,2)*opacity_table%opal95_surface_opacity(temp_index,density_index+1) &
              + dweight_rho(i,3)*opacity_table%opal95_surface_opacity(temp_index,density_index+2) &
              + dweight_rho(i,4)*opacity_table%opal95_surface_opacity(temp_index,density_index+3)
      end do
! CHECK ON WHETHER RHO IS OUTSIDE OF THE TABLE AND NEEDS EXTRAPOLATION
      if (extrap_lo) then
! FACTOR IN R
         delta_logr = logr - logr_lo_edge
! CORRECT CAPPA BY USING THE DERIVATIVE AT THE BOUNDARY.
         do i = 1,4
            logcappa_at_t(i) = logcappa_at_t(i) + delta_logr*dlogcappa_dlogr_at_t(i)
         end do
      else if (extrap_hi) then
         do i = 1,4
            if (extrap_hi_row(i)) then
               delta_logr = logr - logr_hi_edge(i)
               logcappa_at_t(i) = logcappa_at_t(i) + delta_logr*dlogcappa_dlogr_at_t(i)
            endif
         end do
      endif
! INTERPOLATE FOR LOG CAPPA IN T.
      log10_opacity = lagrange4(weight_t, logcappa_at_t)
! D LOG CAPPA/D LOG T
      dlnkap_dlnt = lagrange4(dweight_t, logcappa_at_t)
! D LOG CAPPA/D LOG R
      dlnkap_dlnrho = lagrange4(weight_t, dlogcappa_dlogr_at_t)
! CORRECT FROM DERIVATE AT FIXED R TO DERIVATIVE AT FIXED RHO.
      dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
      opacity = pow(1.0d1, log10_opacity)
      return
end subroutine interpolate_opal95_2d

!----------------------------------------------------------------------
! interpolate_opal95_3d (formerly opal95_interp3d)
!----------------------------------------------------------------------
! MHP 7/98 GET OPAL95 OPACITY FROM TABLE AT FIXED Z
! 3D interpolation (X, T, rho) in the fixed-Z table produced by
! build_opal95_fixed_z_table. Called from compute_opal95_opacity when Z matches the model metal
! fraction but X does not match the surface composition.
! 2026 wave 3 (R5): the T/rho stencil (indices, weights, log R
! extrapolation state) that compute_opal95_opacity computes per call used to reach
! this routine as hidden module state (opacity_table%opal95_index_*,
! opal95_weight_*, opal95_logr*, opal95_extrap_*); it now arrives as
! explicit intent(in) arguments, passed straight from those members by
! compute_opal95_opacity. The opacity tables themselves stay in opacity_table.
subroutine interpolate_opal95_3d(index_x, weight_x, index_t, index_rho, &
     weight_t, dweight_t, weight_rho, dweight_rho, &
     logr, logr_lo_edge, logr_hi_edge, extrap_lo, extrap_hi, extrap_hi_row, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)


! X stencil: row 1 of compute_opal95_opacity's opal95_index_x/opal95_weight_x.
      integer, intent(in) :: index_x(4)
      real(dp), intent(in) :: weight_x(4)
      integer, intent(in) :: index_t, index_rho(4)
      real(dp), intent(in) :: weight_t(4), dweight_t(4), &
           weight_rho(4,4), dweight_rho(4,4)
      real(dp), intent(in) :: logr, logr_lo_edge, logr_hi_edge(4)
      logical, intent(in) :: extrap_lo, extrap_hi, extrap_hi_row(4)
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

      real(dp) :: logcappa_at_t(4), dlogcappa_dlogr_at_t(4)
      real(dp) :: logcappa_at_x(4), dlogcappa_dlogt_at_x(4), &
           dlogcappa_dlogr_at_x(4)
      integer :: i, j, x_table_index, temp_index, density_index
      real(dp) :: delta_logr

! FIND OPACITY AT EACH OF THE 4 NEARBY VALUES OF X
      do j = 1,4
         x_table_index = index_x(j)
! FOR EACH X, GET CAPPA FOR 4 VALUES OF T
         do i = 1,4
            temp_index = index_t + i - 1
! FOR EACH T, GET CAPPA FOR 4 VALUES OF R
            density_index = index_rho(i)
! LOG CAPPA AT DESIRED RHO FOR EACH OF THE 4 DESIRED T.
            logcappa_at_t(i) = weight_rho(i,1)*opacity_table%opal95_fixed_z_opacity(x_table_index,temp_index,density_index) &
                 + weight_rho(i,2)*opacity_table%opal95_fixed_z_opacity(x_table_index,temp_index,density_index+1) &
                 + weight_rho(i,3)*opacity_table%opal95_fixed_z_opacity(x_table_index,temp_index,density_index+2) &
                 + weight_rho(i,4)*opacity_table%opal95_fixed_z_opacity(x_table_index,temp_index,density_index+3)
! D LOG CAPPA/D LOG R FOR EACH OF THE 4 DESIRED T.
            dlogcappa_dlogr_at_t(i) = dweight_rho(i,1)*opacity_table%opal95_fixed_z_opacity(x_table_index,temp_index,density_index) &
                 + dweight_rho(i,2)*opacity_table%opal95_fixed_z_opacity(x_table_index,temp_index,density_index+1) &
                 + dweight_rho(i,3)*opacity_table%opal95_fixed_z_opacity(x_table_index,temp_index,density_index+2) &
                 + dweight_rho(i,4)*opacity_table%opal95_fixed_z_opacity(x_table_index,temp_index,density_index+3)
         end do
! CHECK ON WHETHER RHO IS OUTSIDE OF THE TABLE AND NEEDS EXTRAPOLATION
         if (extrap_lo) then
! FACTOR IN R
            delta_logr = logr - logr_lo_edge
! CORRECT CAPPA BY USING THE DERIVATIVE AT THE BOUNDARY.
            do i = 1,4
               logcappa_at_t(i) = logcappa_at_t(i) + delta_logr*dlogcappa_dlogr_at_t(i)
            end do
         else if (extrap_hi) then
            do i = 1,4
               if (extrap_hi_row(i)) then
                  delta_logr = logr - logr_hi_edge(i)
                  logcappa_at_t(i) = logcappa_at_t(i) + delta_logr*dlogcappa_dlogr_at_t(i)
               endif
            end do
         endif
! INTERPOLATE FOR LOG CAPPA IN T.
         logcappa_at_x(j) = lagrange4(weight_t, logcappa_at_t)
! D LOG CAPPA/D LOG T
         dlogcappa_dlogt_at_x(j) = lagrange4(dweight_t, logcappa_at_t)
! D LOG CAPPA/D LOG R
         dlogcappa_dlogr_at_x(j) = lagrange4(weight_t, dlogcappa_dlogr_at_t)
      end do
! INTERPOLATE FOR LOG CAPPA IN X.
      log10_opacity = lagrange4(weight_x, logcappa_at_x)
! INTERPOLATE FOR QOT IN X.
      dlnkap_dlnt = lagrange4(weight_x, dlogcappa_dlogt_at_x)
! INTERPOLATE FOR QOD IN X.
      dlnkap_dlnrho = lagrange4(weight_x, dlogcappa_dlogr_at_x)
! CORRECT FROM DERIVATE AT FIXED R TO DERIVATIVE AT FIXED RHO.
      dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
      opacity = pow(1.0d1, log10_opacity)
      return
end subroutine interpolate_opal95_3d

!----------------------------------------------------------------------
! interpolate_opal95_4d (formerly opal95_interp4d)
!----------------------------------------------------------------------
! MHP 7/98 GET OPAL95 OPACITY FROM FULL SET OF OPACITY TABLES.
! 4D interpolation (Z, X, T, rho) across the full OPAL95 table set.
! Called from compute_opal95_opacity when neither X nor Z match a precomputed
! reduced table.
! 2026 wave 3 (R5): the T/rho stencil (indices, weights, log R
! extrapolation state) that compute_opal95_opacity computes per call used to reach
! this routine as hidden module state (opacity_table%opal95_index_*,
! opal95_weight_*, opal95_logr*, opal95_extrap_*); it now arrives as
! explicit intent(in) arguments, passed straight from those members by
! compute_opal95_opacity. The opacity tables themselves stay in opacity_table.
subroutine interpolate_opal95_4d(index_z, weight_z, index_x, weight_x, &
     index_t, index_rho, weight_t, dweight_t, weight_rho, dweight_rho, &
     logr, logr_lo_edge, logr_hi_edge, extrap_lo, extrap_hi, extrap_hi_row, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)


! Z stencil, and the X stencil (index_x(k,:), weight_x(k,:)) at each
! of its four Z rows.
      integer, intent(in) :: index_z, index_x(4,4)
      real(dp), intent(in) :: weight_z(4), weight_x(4,4)
      integer, intent(in) :: index_t, index_rho(4)
      real(dp), intent(in) :: weight_t(4), dweight_t(4), &
           weight_rho(4,4), dweight_rho(4,4)
      real(dp), intent(in) :: logr, logr_lo_edge, logr_hi_edge(4)
      logical, intent(in) :: extrap_lo, extrap_hi, extrap_hi_row(4)
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

      real(dp) :: logcappa_at_t(4), dlogcappa_dlogr_at_t(4)
      real(dp) :: logcappa_at_x(4), dlogcappa_dlogt_at_x(4), &
           dlogcappa_dlogr_at_x(4)
      real(dp) :: logcappa_at_z(4), dlogcappa_dlogt_at_z(4), &
           dlogcappa_dlogr_at_z(4)
      integer :: i, j, k, z_table_offset, xz_table_index, temp_index, &
           density_index
      real(dp) :: delta_logr

! FIND OPACITY AT EACH OF THE 4 NEARBY VALUES OF Z
      do k = 1,4
         z_table_offset = opacity_table%opal95_table_start_index(index_z+k-1)
! FIND OPACITY AT EACH OF THE 4 NEARBY VALUES OF X
         do j = 1,4
            xz_table_index = index_x(k,j) + z_table_offset
! FOR EACH X, GET CAPPA FOR 4 VALUES OF T
            do i = 1,4
               temp_index = index_t + i - 1
! FOR EACH T, GET CAPPA FOR 4 VALUES OF R
               density_index = index_rho(i)
! LOG CAPPA AT DESIRED RHO FOR EACH OF THE 4 DESIRED T.
               logcappa_at_t(i) = weight_rho(i,1)*opacity_table%opal95_full_opacity(xz_table_index,temp_index,density_index) + &
                    weight_rho(i,2)*opacity_table%opal95_full_opacity(xz_table_index,temp_index,density_index+1) &
                    + weight_rho(i,3)*opacity_table%opal95_full_opacity(xz_table_index,temp_index,density_index+2) &
                    + weight_rho(i,4)*opacity_table%opal95_full_opacity(xz_table_index,temp_index,density_index+3)
! D LOG CAPPA/D LOG R FOR EACH OF THE 4 DESIRED T.
               dlogcappa_dlogr_at_t(i) = dweight_rho(i,1)*opacity_table%opal95_full_opacity(xz_table_index,temp_index,density_index) &
                    + dweight_rho(i,2)*opacity_table%opal95_full_opacity(xz_table_index,temp_index,density_index+1) &
                    + dweight_rho(i,3)*opacity_table%opal95_full_opacity(xz_table_index,temp_index,density_index+2) &
                    + dweight_rho(i,4)*opacity_table%opal95_full_opacity(xz_table_index,temp_index,density_index+3)
            end do
! CHECK ON WHETHER RHO IS OUTSIDE OF THE TABLE AND NEEDS EXTRAPOLATION
            if (extrap_lo) then
! FACTOR IN R
               delta_logr = logr - logr_lo_edge
! CORRECT CAPPA BY USING THE DERIVATIVE AT THE BOUNDARY.
               do i = 1,4
                  logcappa_at_t(i) = logcappa_at_t(i) + delta_logr*dlogcappa_dlogr_at_t(i)
               end do
            else if (extrap_hi) then
               do i = 1,4
                  if (extrap_hi_row(i)) then
                     delta_logr = logr - logr_hi_edge(i)
                     logcappa_at_t(i) = logcappa_at_t(i) + delta_logr*dlogcappa_dlogr_at_t(i)
                  endif
               end do
            endif
! INTERPOLATE FOR LOG CAPPA IN T. (Inline rather than lagrange4: the
! original accumulates three terms, stores, then adds the fourth.)
            logcappa_at_x(j) = weight_t(1)*logcappa_at_t(1) + weight_t(2)*logcappa_at_t(2) &
                 + weight_t(3)*logcappa_at_t(3)
            logcappa_at_x(j) = logcappa_at_x(j) + weight_t(4)*logcappa_at_t(4)
! D LOG CAPPA/D LOG T
            dlogcappa_dlogt_at_x(j) = dweight_t(1)*logcappa_at_t(1) + dweight_t(2)*logcappa_at_t(2) &
                 + dweight_t(3)*logcappa_at_t(3)
            dlogcappa_dlogt_at_x(j) = dlogcappa_dlogt_at_x(j) + dweight_t(4)*logcappa_at_t(4)
! D LOG CAPPA/D LOG R
            dlogcappa_dlogr_at_x(j) = weight_t(1)*dlogcappa_dlogr_at_t(1) + weight_t(2)*dlogcappa_dlogr_at_t(2) &
                 + weight_t(3)*dlogcappa_dlogr_at_t(3)
            dlogcappa_dlogr_at_x(j) = dlogcappa_dlogr_at_x(j) + weight_t(4)*dlogcappa_dlogr_at_t(4)
         end do
! INTERPOLATE FOR LOG CAPPA IN X. (Inline rather than lagrange4: the
! weights are the strided row weight_x(k,1:4).)
         logcappa_at_z(k) = weight_x(k,1)*logcappa_at_x(1) + weight_x(k,2)*logcappa_at_x(2) + &
              weight_x(k,3)*logcappa_at_x(3) + &
              weight_x(k,4)*logcappa_at_x(4)
! INTERPOLATE FOR QOT IN X.
         dlogcappa_dlogt_at_z(k) = weight_x(k,1)*dlogcappa_dlogt_at_x(1) + weight_x(k,2)*dlogcappa_dlogt_at_x(2) + &
              weight_x(k,3)*dlogcappa_dlogt_at_x(3) + &
              weight_x(k,4)*dlogcappa_dlogt_at_x(4)
! INTERPOLATE FOR QOD IN X.
         dlogcappa_dlogr_at_z(k) = weight_x(k,1)*dlogcappa_dlogr_at_x(1) + weight_x(k,2)*dlogcappa_dlogr_at_x(2) + &
              weight_x(k,3)*dlogcappa_dlogr_at_x(3) + &
              weight_x(k,4)*dlogcappa_dlogr_at_x(4)
      end do
! INTERPOLATE FOR LOG CAPPA IN Z.
      log10_opacity = lagrange4(weight_z, logcappa_at_z)
! INTERPOLATE FOR QOT IN Z.
      dlnkap_dlnt = lagrange4(weight_z, dlogcappa_dlogt_at_z)
! INTERPOLATE FOR QOD IN Z.
      dlnkap_dlnrho = lagrange4(weight_z, dlogcappa_dlogr_at_z)
! CORRECT FROM DERIVATE AT FIXED R TO DERIVATIVE AT FIXED RHO.
      dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
      opacity = pow(1.0d1, log10_opacity)
      return
end subroutine interpolate_opal95_4d

end module opal95_opacity
