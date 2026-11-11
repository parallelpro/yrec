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
! The tables and the lookup's cursors are an opal95_table_set
! (opacity_tables), which kap_lib holds and passes in; each lookup's
! stencil is a local opal95_stencil.
! Thresholds and limits are in opacity_tables and below; the sources and
! the code's history are at the end of the file.
module opal95_opacity
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: n_opal95_d, n_opal95_t, n_opal95_tables, &
           n_opal95_x, n_opal95_z, opal95_composition_tol, &
           opal95_missing_opacity, opal95_missing_test, opal95_table_set, &
           stencil4_locate_opal95, kap_settings_t, kap_table_paths_t, &
           opacity_table_state
      use opacity_source, only: kap_point, kap_result, kap_source
      use numerics_lib, only: interp, lagrange4
      use kap_messages, only: kap_log, kap_report, kap_terminal, message_length, &
           unwritten_line
      use opal95_reader, only: read_opal95_tables
      use math_lib
      implicit none
      private
      public :: opal95_source
      public :: build_opal95_fixed_z_table, build_opal95_surface_table, compute_opal95_opacity

! The tables the set lacks: X = 0.95 (table 9) from Z = 0.04 up, and
! X = 1 - Z (table 10) at Z = 0.1; a stencil above X = 0.8 takes the
! X tables that exist there
      real(dp), parameter :: no_x095_from_z = 0.04d0, no_x1mz_from_z = 0.1d0, &
           high_x = 0.8d0
! The low-X, low-T corner the lookup does not check (see
! compute_opal95_opacity): X below low_corner_x with the T stencil's
! first row at or below low_corner_t_row
      real(dp), parameter :: low_corner_x = 0.2d0
      integer, parameter :: low_corner_t_row = 5

! A lookup's stencil (compute_opal95_opacity's local): the Lagrange
! weights in Z, in X for each Z row (weight_x(:,k) for row k), in T (with
! the T derivative's) and
! in log R for each T row (with the log R derivative's); log R; and where
! log R is extrapolated from: below the table (extrap_lo, from its first
! column, logr_lo_edge) or, on the T rows of extrap_hi_row, past the
! row's last filled column (logr_hi_edge). The stencil's first indices
! are the table set's cursors (opal95_table_set%index_*).
      type :: opal95_stencil
         real(dp) :: weight_z(4), weight_x(4,4), weight_t(4), dweight_t(4), &
              weight_rho(4,4), dweight_rho(4,4)
         real(dp) :: logr, logr_lo_edge, logr_hi_edge(4)
         logical :: extrap_lo, extrap_hi, extrap_hi_row(4)
      end type opal95_stencil

contains


!----------------------------------------------------------------------
! build_opal95_fixed_z_table (formerly opal95_fixed_z_table)
!----------------------------------------------------------------------
! Builds the fixed-Z table (in X, T and log R) by 4-point Lagrange
! interpolation in Z of the full set (tbl%full_opacity), for each X table
! from the Z tables it has: X tables 1-8 at every Z, table 9 (X = 0.95)
! below Z = 0.04 only, table 10 (X = 1 - Z) below Z = 0.1. A Z above the
! top node (0.1) is reported and returns ierr = 1.
subroutine build_opal95_fixed_z_table(tbl, metal_fraction, ierr)
      type(opal95_table_set), intent(inout) :: tbl

      integer, intent(out) :: ierr

      real(dp), intent(in) :: metal_fraction

      real(dp) :: table_z_nodes(4), z_weight(4), z_weight_hix(4)
      integer :: i, j, k, z_table_index, z_table_index2
      integer :: table1_index, table2_index, table3_index, table4_index
      character(len=message_length) :: lines(2)

      ierr = 0
      tbl%fixed_z = metal_fraction
!  Find 4 nearest tables in Z to desired value.
! (the last stencil, the top four nodes, for Z above the next-to-last
! node; only a Z above the top node is refused)
      do i = 3,n_opal95_z-1
         if (tbl%grid_z(i) >= metal_fraction) then
            z_table_index = i - 2
            exit
         end if
      end do
      if (i > (n_opal95_z-1) .and. metal_fraction <= tbl%grid_z(n_opal95_z)) then
         z_table_index = n_opal95_z - 3
      else if (i > (n_opal95_z-1)) then
      lines = unwritten_line
      write(lines,5)metal_fraction
      call kap_terminal(lines)
    5 format(1x,'DESIRED Z',f10.6,'OUTSIDE OP95 TABLE RANGE'/ &
           ' RUN STOPPED')
      ierr = 1
      return
      end if
!  Find interpolation factors in Z.
      do i = 1,4
         table_z_nodes(i) = tbl%grid_z(z_table_index+i-1)
      end do
      call interp(table_z_nodes, z_weight, x_eval=metal_fraction)
!  For X tables 1-8, all values of Z present; perform interpolation.
      do i = 1,8
! Indices for 4 desired compositions.
         table1_index = tbl%table_start_index(z_table_index)+i
         table2_index = tbl%table_start_index(z_table_index+1)+i
         table3_index = tbl%table_start_index(z_table_index+2)+i
         table4_index = tbl%table_start_index(z_table_index+3)+i
         do j = 1,n_opal95_t
            do k = 1,n_opal95_d
               tbl%fixed_z_opacity(i,j,k) = z_weight(1)*tbl%full_opacity(table1_index,j,k)+ &
               z_weight(2)*tbl%full_opacity(table2_index,j,k)+z_weight(3)*tbl%full_opacity(table3_index,j,k)+ &
               z_weight(4)*tbl%full_opacity(table4_index,j,k)
            end do
         end do
      end do
!  Table 9 (X = 0.95) exists below Z = 0.04 only: built below it, from
!  a stencil of Z nodes below 0.04
      if (metal_fraction < no_x095_from_z) then
      if (tbl%grid_z(z_table_index+3) > no_x095_from_z) then
         z_table_index2 = n_opal95_z - 6
         do i = 1,4
            table_z_nodes(i) = tbl%grid_z(z_table_index2+i-1)
         end do
         call interp(table_z_nodes, z_weight_hix, x_eval=metal_fraction)
      else
         z_table_index2 = z_table_index
         do i = 1,4
            z_weight_hix(i) = z_weight(i)
         end do
      end if
      table1_index = tbl%table_start_index(z_table_index2)+9
      table2_index = tbl%table_start_index(z_table_index2+1)+9
      table3_index = tbl%table_start_index(z_table_index2+2)+9
      table4_index = tbl%table_start_index(z_table_index2+3)+9
      do j = 1,n_opal95_t
         do k = 1,n_opal95_d
            tbl%fixed_z_opacity(9,j,k) = z_weight_hix(1)*tbl%full_opacity(table1_index,j,k)+ &
            z_weight_hix(2)*tbl%full_opacity(table2_index,j,k)+z_weight_hix(3)*tbl%full_opacity(table3_index,j,k)+ &
            z_weight_hix(4)*tbl%full_opacity(table4_index,j,k)
         end do
      end do
      end if
!  Table 10 (X = 1-Z) is always built. There is no X = 1-Z table at
!  Z = 0.1, so a stencil that reaches Z = 0.1 is replaced by the four
!  nodes below it (extrapolating for Z above 0.08).
      if (tbl%grid_z(z_table_index+3) >= no_x1mz_from_z) then
         z_table_index2 = n_opal95_z - 4
         do i = 1,4
            table_z_nodes(i) = tbl%grid_z(z_table_index2+i-1)
         end do
         call interp(table_z_nodes, z_weight_hix, x_eval=metal_fraction)
      else
         z_table_index2 = z_table_index
         do i = 1,4
            z_weight_hix(i) = z_weight(i)
         end do
      end if
      table1_index = tbl%table_start_index(z_table_index2)+10
      table2_index = tbl%table_start_index(z_table_index2+1)+10
      table3_index = tbl%table_start_index(z_table_index2+2)+10
      table4_index = tbl%table_start_index(z_table_index2+3)+10
      do j = 1,n_opal95_t
         do k = 1,n_opal95_d
            tbl%fixed_z_opacity(10,j,k) = z_weight_hix(1)*tbl%full_opacity(table1_index,j,k)+ &
            z_weight_hix(2)*tbl%full_opacity(table2_index,j,k)+z_weight_hix(3)*tbl%full_opacity(table3_index,j,k)+ &
            z_weight_hix(4)*tbl%full_opacity(table4_index,j,k)
         end do
      end do
end subroutine build_opal95_fixed_z_table

!----------------------------------------------------------------------
! build_opal95_surface_table (formerly opal95_surface_table)
!----------------------------------------------------------------------
! Builds the surface table (in T and log R) at hydrogen_fraction by
! 4-point Lagrange interpolation in X of the fixed-Z table set
! (tbl%fixed_z_opacity), from the X tables that exist at its Z.
subroutine build_opal95_surface_table(tbl, hydrogen_fraction)
      type(opal95_table_set), intent(inout) :: tbl


      real(dp), intent(in) :: hydrogen_fraction

      real(dp) :: table_x_nodes(4), x_weight(4)
      integer :: i, j, k, x_table_index
      integer :: table1_index, table2_index, table3_index, table4_index
! x_table_index normally names the first of four consecutive X tables;
! this value means "tables 6, 7, 8 and 10" (table 9 absent).
      integer, parameter :: x_stencil_6_7_8_10 = 11

      tbl%surface_x = hydrogen_fraction
!  Find 4 nearest tables in X to desired value.
      if (hydrogen_fraction <= high_x .or. &
           tbl%fixed_z <= no_x095_from_z) then
! All four X tables present: Either X is below the high-X gap, or
! the high X tables are present at this low Z.
         do i = 3,n_opal95_x-1
            if (tbl%grid_x(i) >= hydrogen_fraction) then
               x_table_index = i - 2
               exit
            end if
         end do
         if (i > (n_opal95_x-1)) then
         x_table_index = n_opal95_x - 3
         end if
         do i = 1,4
            table_x_nodes(i) = tbl%grid_x(x_table_index+i-1)
         end do
      else if (tbl%fixed_z >= no_x1mz_from_z) then
! Use tables 6-9
         x_table_index = 6
         do i = 1,4
            table_x_nodes(i) = tbl%grid_x(x_table_index+i-1)
         end do
      else
! If Z is between 0.04 and 0.1, tables 1-8 and 10 exist.
! Since we have already determined that X > 0.8, use 6-8 and 10.
         x_table_index = x_stencil_6_7_8_10
         table_x_nodes(1) = tbl%grid_x(6)
         table_x_nodes(2) = tbl%grid_x(7)
         table_x_nodes(3) = tbl%grid_x(8)
         table_x_nodes(4) = tbl%grid_x(10)
      end if
!
!  Find interpolation factors in X.
      call interp(table_x_nodes, x_weight, x_eval=hydrogen_fraction)
! Indices for 4 desired compositions.
      if (x_table_index /= x_stencil_6_7_8_10) then
         table1_index = x_table_index
         table2_index = x_table_index + 1
         table3_index = x_table_index + 2
         table4_index = x_table_index + 3
      else
         table1_index = 6
         table2_index = 7
         table3_index = 8
         table4_index = 10
      end if
      do j = 1,n_opal95_t
         do k = 1,n_opal95_d
            tbl%surface_opacity(j,k) = x_weight(1)*tbl%fixed_z_opacity(table1_index,j,k) + &
            x_weight(2)*tbl%fixed_z_opacity(table2_index,j,k) + x_weight(3)*tbl%fixed_z_opacity(table3_index,j,k) + &
            x_weight(4)*tbl%fixed_z_opacity(table4_index,j,k)
         end do
      end do
end subroutine build_opal95_surface_table

!----------------------------------------------------------------------
! compute_opal95_opacity (formerly getopal95)
!----------------------------------------------------------------------
! The opacity at (log rho, log T, X, Z) and its log derivatives, by
! 4-point Lagrange interpolation: it finds the stencil in T and log R
! (its first indices are tbl's warm-start cursors), checks the stencil
! against the table's empty high-R corner, and interpolates
!   2D (T, log R) in the surface table, at the surface X and the fixed Z
!     (each within opal95_composition_tol);
!   3D (X, T, log R) in the fixed-Z table, at the fixed Z and another X;
!   4D (Z, X, T, log R) in the full set, at any other Z.
! A log T off the table is reported and returns ierr = 1.
subroutine compute_opal95_opacity(tbl, log10_density, log10_temperature, hydrogen_fraction, &
     metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ierr)
      type(opal95_table_set), intent(inout) :: tbl

      real(dp), intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      integer, intent(out) :: ierr

      real(dp) :: interp_nodes(4), weight(4), dweight(4)
      integer :: i, k, j
      logical :: density_shifted
      integer :: x_table_index, density_base_index, temp_row_index
      real(dp) :: extrap_logr, z_at_table
! the call's stencil: its weights and where log R is extrapolated from
! (its first indices are tbl's cursors)
      type(opal95_stencil) :: st
      character(len=message_length) :: lines(1)

!     log R = log rho - 3 log T6, and log T within the table
      ierr = 0
      st%logr = log10_density - 3.0d0*(log10_temperature-6.0d0)
      if (log10_temperature < tbl%grid_logt(1) .or. &
           log10_temperature > tbl%grid_logt(n_opal95_t)) then
         lines = unwritten_line
         write(lines,5) log10_temperature
    5    format(' LOG T OF',f11.6,'OUT OF OPAL 95 TABLE RANGE')
         call kap_terminal(lines)
         ierr = 1
         return
      end if
!     Check to see if extrapolation below the first table element
!     in density is needed.
      if (st%logr < tbl%grid_logr(1)) then
         st%extrap_lo = .true.
      else
         st%extrap_lo = .false.
      end if
!     (No low-R range check: below the table's first column the opacity
!     is extrapolated linearly in log R, extrapolate_in_logr.)
!     get interpolation factors in T
      call stencil4_locate_opal95(tbl%grid_logt, n_opal95_t, log10_temperature, tbl%index_t)
!     T interpolation factors
      do i = 1,4
         interp_nodes(i) = tbl%grid_logt(tbl%index_t+i-1)
      end do
!     Get interpolation factors in T
      call interp(interp_nodes, weight, dweight, log10_temperature)
      do i = 1,4
         st%weight_t(i) = weight(i)
         st%dweight_t(i) = dweight(i)
      end do
!     Get indices in rho for each of the 4 T values.
      call stencil4_locate_opal95(tbl%grid_logr, n_opal95_d, st%logr, tbl%index_rho(1))
!     Now check if within portion of table with data.
!     Note that a point can be within the table and still have the
!     4x4 matrix of T,rho around it outside the table.
!     The low-rho, low-T corner at low X (X < 0.2, the stencil's first T
!     row at or below the fifth) is not checked: its empty region depends
!     on Z, and the check the F77 code postponed to the 2D/3D/4D choice
!     was never written, so the stencil is used as found, without the
!     high-rho shift or extrapolation.
      density_shifted = .false.
      st%extrap_hi = .false.
      st%extrap_hi_row = .false.
      if (.not. (hydrogen_fraction < low_corner_x .and. tbl%index_t <= low_corner_t_row)) then
!     Check high rho, high T regime(because of the table geometry, the
!     point will fall outside the table first at the highest T).
      if (tbl%fixed_z_opacity(tbl%index_x(1,1),tbl%index_t+3, &
           tbl%index_rho(1)+3) >= opal95_missing_test) then
         density_shifted = .true.
!        Now there are two possible solutions: We can use a different set
!        of densities for each temperature and stay within the table, or
!        the desired density could be outside the table at one or more
!        temperatures.  If the desired density is outside the table,
!        linear extrapolation in R at fixed T is used.
         do k = 1,4
            if (tbl%fixed_z_opacity(tbl%index_x(1,1),tbl%index_t+k-1, &
                 tbl%index_rho(1)+1) >= opal95_missing_test) then
               st%extrap_hi = .true.
               st%extrap_hi_row(k) = .true.
            else
               st%extrap_hi_row(k) = .false.
            end if
         end do
!        Find 4 nearest elements within the table at each of the 4 temperatures.
         x_table_index = tbl%index_x(1,1)
         density_base_index = tbl%index_rho(1)
         do i = 1,4
            temp_row_index = tbl%index_t+i-1
            if (tbl%fixed_z_opacity(x_table_index,temp_row_index, &
                 density_base_index+3) > opal95_missing_test) then
               do j = density_base_index+2,1,-1
                  if (tbl%full_opacity(x_table_index,temp_row_index,j) <= opal95_missing_test) then
                     tbl%index_rho(i) = j - 3
                     density_base_index = j - 3
                     exit
                  end if
               end do
            else
               tbl%index_rho(i) = density_base_index
            end if
         end do
      end if
      end if
!     Get interpolation factors in rho for each of the 4 T values.
!     Assign all the same values if grid in R is same for all T.
      if (.not.density_shifted) then
         do i = 1,4
            interp_nodes(i) = tbl%grid_logr(tbl%index_rho(1)+i-1)
         end do
         if (.not.st%extrap_lo) then
            call interp(interp_nodes, weight, dweight, st%logr)
         else
            st%logr_lo_edge = tbl%grid_logr(1)
            call interp(interp_nodes, weight, dweight, st%logr_lo_edge)
         end if
         do j = 1,4
            tbl%index_rho(j) = tbl%index_rho(1)
            do i = 1,4
               st%weight_rho(j,i) = weight(i)
               st%dweight_rho(j,i) = dweight(i)
            end do
         end do
      else
         do j = 1,4
            do i = 1,4
               interp_nodes(i) = tbl%grid_logr(tbl%index_rho(j)+i-1)
            end do
!           If extrapolation is being performed, assign R=redge
!           and extrapolate linearly later
            if (.not.st%extrap_hi_row(j)) then
               call interp(interp_nodes, weight, dweight, st%logr)
            else
               extrap_logr = tbl%grid_logr(tbl%index_rho(j)+3)
               call interp(interp_nodes, weight, dweight, extrap_logr)
               st%logr_hi_edge(j) = extrap_logr
            end if
            do i = 1,4
               st%weight_rho(j,i) = weight(i)
               st%dweight_rho(j,i) = dweight(i)
            end do
         end do
      end if
!     Determine whether a 2d (rho,T); 3d (X,rho,T); or 4d (Z,X,rho,T)
!     interpolation is needed to get the opacity.
      if (abs(metal_fraction-tbl%fixed_z)/max(tbl%fixed_z,1.0d-6) <= opal95_composition_tol) then
         if (abs(hydrogen_fraction-tbl%surface_x) <= opal95_composition_tol) then
!           2d interpolation in surface table
            call interpolate_opal95_2d(tbl, st, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
            return
         else
!           3d interpolation in fixed Z table (X,T,rho)
!           get interpolation factors in X.
         call locate_opal95_x_stencil(1, tbl%fixed_z, hydrogen_fraction)
         call interpolate_opal95_3d(tbl, st, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
         return
      end if
      end if
!     4d interpolation in Z,X,T,rho
!     get nearest tables in Z.
      if (metal_fraction > tbl%grid_z(tbl%index_z+2)) then
         do i = tbl%index_z+3,n_opal95_z-1
            if (metal_fraction < tbl%grid_z(i)) then
               tbl%index_z = i - 2
               exit
            end if
         end do
         if (i > (n_opal95_z-1)) then
         tbl%index_z = n_opal95_z - 3
         end if
      else
         do i = tbl%index_z+1,2,-1
            if (metal_fraction >= tbl%grid_z(i)) then
               tbl%index_z = i - 1
               exit
            end if
         end do
         if (i < (2)) then
         tbl%index_z = 1
         end if
      end if
!     (the last stencil at most)
      tbl%index_z = min(tbl%index_z,n_opal95_z-3)
!     Z interpolation factors
      do i = 1,4
         interp_nodes(i) = tbl%grid_z(tbl%index_z+i-1)
      end do
!     Interpolation factors in Z
      call interp(interp_nodes, weight, x_eval=metal_fraction)
      do i = 1,4
         st%weight_z(i) = weight(i)
      end do
!     Get nearest tables in X at each value of Z
      do k = 1,4
         z_at_table = tbl%grid_z(tbl%index_z+k-1)
         call locate_opal95_x_stencil(k, z_at_table, hydrogen_fraction)
      end do
      call interpolate_opal95_4d(tbl, st, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
      return

contains

! locate_opal95_x_stencil
! For Z row k of the stencil (Z value z_here), pick the four X tables
! tbl%index_x(k,1:4) bracketing hydrogen_fraction and
! store their weights in st%weight_x(1:4,k). The
! warm-start search moves tbl%index_x(k,1) from its previous value;
! at Z >= 0.04 and X >= 0.8 the missing high-X tables (no X=0.95 for
! Z >= 0.04, no X=1-Z for Z >= 0.1) force the stencil to the top
! available tables. The 3D path calls it for k = 1 at the fixed Z, the
! 4D path for each of its Z rows.
subroutine locate_opal95_x_stencil(k, z_here, hydrogen_fraction)
      integer, intent(in) :: k
      real(dp), intent(in) :: z_here, hydrogen_fraction
      integer :: i, j, x_table_index, x_shift_base
      real(dp) :: interp_nodes(4), weight(4)

      if (hydrogen_fraction > tbl%grid_x(tbl%index_x(k,1)+2)) then
         do i = tbl%index_x(k,1)+3,n_opal95_x-1
            if (hydrogen_fraction < tbl%grid_x(i)) then
               tbl%index_x(k,1) = i - 2
               exit
            end if
         end do
         if (i > n_opal95_x-1) tbl%index_x(k,1) = n_opal95_x - 3
      else
         do i = tbl%index_x(k,1)+1,2,-1
            if (hydrogen_fraction >= tbl%grid_x(i)) then
               tbl%index_x(k,1) = i - 1
               exit
            end if
         end do
         if (i < 2) tbl%index_x(k,1) = 1
      end if
!     Check for absent tables at high X.
      if (z_here >= no_x095_from_z) then
         if (hydrogen_fraction >= high_x) then
            if (z_here < no_x1mz_from_z) then
!              no table 9 (X = 0.95)
               tbl%index_x(k,4) = n_opal95_x
               tbl%index_x(k,3) = n_opal95_x - 2
               tbl%index_x(k,2) = n_opal95_x - 3
               tbl%index_x(k,1) = n_opal95_x - 4
            else
!              No table 9 or table 10
               tbl%index_x(k,4) = n_opal95_x - 2
               tbl%index_x(k,3) = n_opal95_x - 3
               tbl%index_x(k,2) = n_opal95_x - 4
               tbl%index_x(k,1) = n_opal95_x - 5
            end if
         else
            x_shift_base = tbl%index_x(k,1) - 1
            do i = 1,4
               tbl%index_x(k,i) = x_shift_base + i
            end do
         end if
      else
         x_shift_base = tbl%index_x(k,1) - 1
         do i = 1,4
            tbl%index_x(k,i) = x_shift_base + i
         end do
      end if
!     X interpolation factors
      do i = 1,4
         x_table_index = tbl%index_x(k,i)
         interp_nodes(i) = tbl%grid_x(x_table_index)
      end do
      if (tbl%index_x(k,4) == n_opal95_x) then
         interp_nodes(4) = 1.0d0 - z_here
      end if
!     Get interpolation factors in X
      call interp(interp_nodes, weight, x_eval=hydrogen_fraction)
      do j = 1,4
         st%weight_x(j,k) = weight(j)
      end do
end subroutine locate_opal95_x_stencil

end subroutine compute_opal95_opacity

!----------------------------------------------------------------------
! interpolate_opal95_2d, _3d, _4d: the opacity from the stencil st
!----------------------------------------------------------------------
! Lagrange interpolation, four points in each dimension, in the surface
! table (2D: T, log R), the fixed-Z table (3D: X, T, log R) or the full
! set (4D: Z, X, T, log R). Each T row of the stencil gives log kappa
! and d log kappa/d log R at the point's log R (opal95_rows_at_logr,
! extrapolated by extrapolate_in_logr where log R is off the table);
! these are interpolated in T, then X, then Z, and the T derivative
! at fixed R becomes the derivative at fixed rho:
!   dln kappa/dln T|rho = dln kappa/dln T|R - 3 dln kappa/dln R.
pure subroutine interpolate_opal95_2d(tbl, st, opacity, log10_opacity, dlnkap_dlnrho, &
     dlnkap_dlnt)
      type(opal95_table_set), intent(in) :: tbl
      type(opal95_stencil), intent(in) :: st
      real(dp), intent(out) :: opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt
      real(dp) :: logcappa_at_t(4), dlogcappa_dlogr_at_t(4)

      call opal95_rows_at_logr(tbl%surface_opacity, tbl%index_t, tbl%index_rho, st, &
           logcappa_at_t, dlogcappa_dlogr_at_t)
      log10_opacity = lagrange4(st%weight_t, logcappa_at_t)
      dlnkap_dlnt = lagrange4(st%dweight_t, logcappa_at_t)
      dlnkap_dlnrho = lagrange4(st%weight_t, dlogcappa_dlogr_at_t)
      dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
      opacity = pow(1.0d1, log10_opacity)
end subroutine interpolate_opal95_2d

pure subroutine interpolate_opal95_3d(tbl, st, opacity, log10_opacity, dlnkap_dlnrho, &
     dlnkap_dlnt)
      type(opal95_table_set), intent(in) :: tbl
      type(opal95_stencil), intent(in) :: st
      real(dp), intent(out) :: opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt
      real(dp) :: logcappa_at_t(4), dlogcappa_dlogr_at_t(4)
      real(dp) :: logcappa_at_x(4), dlogcappa_dlogt_at_x(4), dlogcappa_dlogr_at_x(4)
      integer :: j

! the X stencil is the first Z row's (the fixed-Z table's Z)
      do j = 1,4
         call opal95_rows_at_logr(tbl%fixed_z_opacity(tbl%index_x(1,j),:,:), tbl%index_t, &
              tbl%index_rho, st, logcappa_at_t, dlogcappa_dlogr_at_t)
         logcappa_at_x(j) = lagrange4(st%weight_t, logcappa_at_t)
         dlogcappa_dlogt_at_x(j) = lagrange4(st%dweight_t, logcappa_at_t)
         dlogcappa_dlogr_at_x(j) = lagrange4(st%weight_t, dlogcappa_dlogr_at_t)
      end do
      log10_opacity = lagrange4(st%weight_x(:,1), logcappa_at_x)
      dlnkap_dlnt = lagrange4(st%weight_x(:,1), dlogcappa_dlogt_at_x)
      dlnkap_dlnrho = lagrange4(st%weight_x(:,1), dlogcappa_dlogr_at_x)
      dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
      opacity = pow(1.0d1, log10_opacity)
end subroutine interpolate_opal95_3d

pure subroutine interpolate_opal95_4d(tbl, st, opacity, log10_opacity, dlnkap_dlnrho, &
     dlnkap_dlnt)
      type(opal95_table_set), intent(in) :: tbl
      type(opal95_stencil), intent(in) :: st
      real(dp), intent(out) :: opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt
      real(dp) :: logcappa_at_t(4), dlogcappa_dlogr_at_t(4)
      real(dp) :: logcappa_at_x(4), dlogcappa_dlogt_at_x(4), dlogcappa_dlogr_at_x(4)
      real(dp) :: logcappa_at_z(4), dlogcappa_dlogt_at_z(4), dlogcappa_dlogr_at_z(4)
      integer :: j, k, xz_slot

! each Z row of the stencil has its own X stencil (index_x(k,:)), the
! slots of the full set being table_start_index(Z index) + X index
      do k = 1,4
         do j = 1,4
            xz_slot = tbl%index_x(k,j) + tbl%table_start_index(tbl%index_z+k-1)
            call opal95_rows_at_logr(tbl%full_opacity(xz_slot,:,:), tbl%index_t, &
                 tbl%index_rho, st, logcappa_at_t, dlogcappa_dlogr_at_t)
            logcappa_at_x(j) = lagrange4(st%weight_t, logcappa_at_t)
            dlogcappa_dlogt_at_x(j) = lagrange4(st%dweight_t, logcappa_at_t)
            dlogcappa_dlogr_at_x(j) = lagrange4(st%weight_t, dlogcappa_dlogr_at_t)
         end do
         logcappa_at_z(k) = lagrange4(st%weight_x(:,k), logcappa_at_x)
         dlogcappa_dlogt_at_z(k) = lagrange4(st%weight_x(:,k), dlogcappa_dlogt_at_x)
         dlogcappa_dlogr_at_z(k) = lagrange4(st%weight_x(:,k), dlogcappa_dlogr_at_x)
      end do
      log10_opacity = lagrange4(st%weight_z, logcappa_at_z)
      dlnkap_dlnt = lagrange4(st%weight_z, dlogcappa_dlogt_at_z)
      dlnkap_dlnrho = lagrange4(st%weight_z, dlogcappa_dlogr_at_z)
      dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
      opacity = pow(1.0d1, log10_opacity)
end subroutine interpolate_opal95_4d

! log kappa and d log kappa/d log R at the point's log R on each of the
! stencil's four T rows of one table (table(T index, log R index)), each
! row with its own log R stencil index_rho(i), then extrapolated where
! log R is off the table (extrapolate_in_logr)
pure subroutine opal95_rows_at_logr(table, index_t, index_rho, st, logcappa_at_t, &
     dlogcappa_dlogr_at_t)
      real(dp), intent(in) :: table(:,:)
      integer, intent(in) :: index_t, index_rho(:)
      type(opal95_stencil), intent(in) :: st
      real(dp), intent(out) :: logcappa_at_t(:), dlogcappa_dlogr_at_t(:)
      integer :: i, temp_index, density_index

      do i = 1,4
         temp_index = index_t + i - 1
         density_index = index_rho(i)
         logcappa_at_t(i) = st%weight_rho(i,1)*table(temp_index,density_index) &
              + st%weight_rho(i,2)*table(temp_index,density_index+1) &
              + st%weight_rho(i,3)*table(temp_index,density_index+2) &
              + st%weight_rho(i,4)*table(temp_index,density_index+3)
         dlogcappa_dlogr_at_t(i) = st%dweight_rho(i,1)*table(temp_index,density_index) &
              + st%dweight_rho(i,2)*table(temp_index,density_index+1) &
              + st%dweight_rho(i,3)*table(temp_index,density_index+2) &
              + st%dweight_rho(i,4)*table(temp_index,density_index+3)
      end do
      call extrapolate_in_logr(st, logcappa_at_t, dlogcappa_dlogr_at_t)
end subroutine opal95_rows_at_logr

! Where log R is off the table, each T row's log kappa, evaluated at the
! edge (the stencil's weights were taken there), is extrapolated
! linearly to the point's log R with the row's derivative at the edge:
! below the table's first column on every row, or past a row's last
! filled column at high R on the rows where that is so.
pure subroutine extrapolate_in_logr(st, logcappa_at_t, dlogcappa_dlogr_at_t)
      type(opal95_stencil), intent(in) :: st
      real(dp), intent(inout) :: logcappa_at_t(:)
      real(dp), intent(in) :: dlogcappa_dlogr_at_t(:)
      real(dp) :: delta_logr
      integer :: i

      if (st%extrap_lo) then
         delta_logr = st%logr - st%logr_lo_edge
         do i = 1,4
            logcappa_at_t(i) = logcappa_at_t(i) + delta_logr*dlogcappa_dlogr_at_t(i)
         end do
      else if (st%extrap_hi) then
         do i = 1,4
            if (st%extrap_hi_row(i)) then
               delta_logr = st%logr - st%logr_hi_edge(i)
               logcappa_at_t(i) = logcappa_at_t(i) + delta_logr*dlogcappa_dlogr_at_t(i)
            end if
         end do
      end if
end subroutine extrapolate_in_logr

!----------------------------------------------------------------------
! The family's descriptor (opacity_source) and its procedures
!----------------------------------------------------------------------
! An interior family that interpolates in Z itself (the 4D path), so it
! needs no second table at another Z; in a helium-burning core it is
! read at its top Z.
function opal95_source(settings) result(source)
      type(kap_settings_t), intent(in) :: settings
      type(kap_source) :: source

      source%name = 'OPAL 1995'
      source%selected = settings%use_opal95_tables
      source%covers_any_z = .true.
      source%load => load_opal95
      source%refresh => refresh_opal95
      source%lookup => lookup_opal95
      source%core_lookup => core_lookup_opal95
end function opal95_source

subroutine load_opal95(settings, paths, envelope_hydrogen_fraction, tables, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(kap_table_paths_t), intent(in) :: paths
      real(dp), intent(in) :: envelope_hydrogen_fraction
      type(opacity_table_state), intent(inout) :: tables
      integer, intent(out) :: ierr

      call read_opal95_tables(tables%opal95, paths%opal95, ierr)
      if (ierr /= 0) return
      call build_opal95_fixed_z_table(tables%opal95, settings%opal95_single_table_z, ierr)
      if (ierr /= 0) return
      call build_opal95_surface_table(tables%opal95, envelope_hydrogen_fraction)
end subroutine load_opal95

subroutine refresh_opal95(settings, hydrogen_fraction, tables, ierr)
      type(kap_settings_t), intent(in) :: settings
      real(dp), intent(in) :: hydrogen_fraction
      type(opacity_table_state), intent(inout) :: tables
      integer, intent(out) :: ierr

      ierr = 0
      call build_opal95_surface_table(tables%opal95, hydrogen_fraction)
end subroutine refresh_opal95

subroutine lookup_opal95(settings, tables, point, res, in_table, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(opacity_table_state), intent(inout) :: tables
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      logical, intent(out) :: in_table
      integer, intent(out) :: ierr

      in_table = .true.
      call compute_opal95_opacity(tables%opal95, point%log10_density, &
           point%log10_temperature, point%hydrogen_fraction, point%metal_fraction, &
           res%opacity, res%log10_opacity, res%dlnkap_dlnrho, res%dlnkap_dlnt, ierr)
end subroutine lookup_opal95

! in a helium-burning core: the tables at their top Z
subroutine core_lookup_opal95(settings, tables, point, res, table_z, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(opacity_table_state), intent(inout) :: tables
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      real(dp), intent(out) :: table_z
      integer, intent(out) :: ierr

      table_z = tables%opal95%grid_z(n_opal95_z)
      call compute_opal95_opacity(tables%opal95, point%log10_density, &
           point%log10_temperature, point%hydrogen_fraction, table_z, &
           res%opacity, res%log10_opacity, res%dlnkap_dlnrho, res%dlnkap_dlnt, ierr)
end subroutine core_lookup_opal95

!----------------------------------------------------------------------
! Sources
!----------------------------------------------------------------------
! The tables: the OPAL opacities (Iglesias & Rogers 1996, ApJ 464, 943)
! and the OP ones (Badnell et al. 2005, MNRAS 360, 458) in the OPAL 1995
! file layout, as YREC ships them in input/opacity/atomic/legacy/opal95/.
! The code: MHP's F77 routines LL95TBL, OP95ZTAB,
! GETOPAL95 and their surface-table and 2D/3D/4D companions (7/98; the
! corner without data, 12/97; the Z search's typo and failsafe, 7/12),
! converted to Fortran 90 in 2026. Since then: the reads checked (L10,
! M220); the fixed-Z table reaches Z = 0.10 (M224); the low corner's
! flags reset at each call (M223); the tables and cursors in a table set
! passed in, the per-call stencil local (KA-12); the three interpolators
! sharing one row helper and one log R extrapolation.

end module opal95_opacity
