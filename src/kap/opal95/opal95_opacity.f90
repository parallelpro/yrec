!----------------------------------------------------------------------
! opal95_opacity
!----------------------------------------------------------------------
! The atomic opacities in the OPAL 1995 table layout (the OPAL and OP
! tables YREC ships, e.g. GS98.OP17): 126 tables in X and Z, each in
! log T and log R = log rho - 3 log T6.
!     build_opal95_fixed_z_table  the tables at the run's Z
!                                 (opal95_single_table_z), from the full
!                                 set opal95_reader read
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
      use opacity_tables, only: n_opal95_d, n_opal95_t, &
           n_opal95_x, n_opal95_z, opal95_composition_tol, &
           opal95_missing_test, opal95_table_set, &
           stencil4_locate_opal95, kap_settings_t, kap_table_paths_t, &
           opacity_table_state
      use opacity_source, only: kap_point, kap_result, kap_source
      use numerics_lib, only: interp, lagrange4
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
! The table a lookup interpolates in (locate_opal95_composition): the
! surface table (2D), the fixed-Z table (3D) or the full set (4D); and
! the floor of the fixed Z in the relative test for "at the fixed Z", so
! that a fixed Z of 0 has a tolerance
      integer, parameter :: path_surface_table = 2, path_fixed_z_table = 3, path_full_set = 4
      real(dp), parameter :: fixed_z_floor = 1.0d-6

! A lookup's stencil (compute_opal95_opacity's local): the Lagrange
! weights in Z, in X for each Z row (weight_x(:,k) for row k), in T (with
! the T derivative's) and in log R for each T row (weight_rho(row,:),
! with the log R derivative's); log R; and where
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
! build_opal95_fixed_z_table
!----------------------------------------------------------------------
! Builds the fixed-Z table (in X, T and log R) by 4-point Lagrange
! interpolation in Z of the full set (tbl%full_opacity), for each X table
! from the Z tables it has: X tables 1-8 at every Z, table 9 (X = 0.95)
! below Z = 0.04 only, table 10 (X = 1 - Z) below Z = 0.1. A Z above the
! top node (0.1) returns ierr = 1 with the reason in message.
subroutine build_opal95_fixed_z_table(tbl, metal_fraction, message, ierr)
      type(opal95_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: metal_fraction
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

! the first of the four Z tables around Z, and their weights; the same
! for an X table that lacks the upper Z tables
      integer :: i, ix, first_z, first_z_high_x
      real(dp) :: z_weight(4), z_weight_high_x(4)

      ierr = 0
      message = ''
      tbl%fixed_z = metal_fraction
!  Find 4 nearest tables in Z to desired value.
! (the last stencil, the top four nodes, for Z above the next-to-last
! node; only a Z above the top node is refused)
      do i = 3,n_opal95_z-1
         if (tbl%grid_z(i) >= metal_fraction) then
            first_z = i - 2
            exit
         end if
      end do
      if (i > (n_opal95_z-1) .and. metal_fraction <= tbl%grid_z(n_opal95_z)) then
         first_z = n_opal95_z - 3
      else if (i > (n_opal95_z-1)) then
         write(message,'(a,f10.6,a,f5.2,a)') 'opal95_single_table_z =', metal_fraction, &
              ' is above the OPAL 1995 tables'' top Z (', tbl%grid_z(n_opal95_z), ')'
         ierr = 1
         return
      end if
      call interp(tbl%grid_z(first_z:first_z+3), z_weight, x_eval=metal_fraction)
!  For X tables 1-8, all values of Z present; perform interpolation.
      do ix = 1,8
         call interpolate_x_table_in_z(tbl, ix, first_z, z_weight)
      end do
!  Table 9 (X = 0.95) exists below Z = 0.04 only: built below it, from
!  a stencil of Z nodes below 0.04
      if (metal_fraction < no_x095_from_z) then
         if (tbl%grid_z(first_z+3) > no_x095_from_z) then
            first_z_high_x = n_opal95_z - 6
            call interp(tbl%grid_z(first_z_high_x:first_z_high_x+3), z_weight_high_x, &
                 x_eval=metal_fraction)
         else
            first_z_high_x = first_z
            z_weight_high_x = z_weight
         end if
         call interpolate_x_table_in_z(tbl, 9, first_z_high_x, z_weight_high_x)
      end if
!  Table 10 (X = 1-Z) is always built. There is no X = 1-Z table at
!  Z = 0.1, so a stencil that reaches Z = 0.1 is replaced by the four
!  nodes below it (extrapolating for Z above 0.08).
      if (tbl%grid_z(first_z+3) >= no_x1mz_from_z) then
         first_z_high_x = n_opal95_z - 4
         call interp(tbl%grid_z(first_z_high_x:first_z_high_x+3), z_weight_high_x, &
              x_eval=metal_fraction)
      else
         first_z_high_x = first_z
         z_weight_high_x = z_weight
      end if
      call interpolate_x_table_in_z(tbl, 10, first_z_high_x, z_weight_high_x)
end subroutine build_opal95_fixed_z_table

! The fixed-Z table's X table ix: the full set's tables of that X at the
! four Z from Z index first_z, combined with z_weight (the slot of X
! index ix at Z index iz is table_start_index(iz) + ix).
subroutine interpolate_x_table_in_z(tbl, ix, first_z, z_weight)
      type(opal95_table_set), intent(inout) :: tbl
      integer, intent(in) :: ix, first_z
      real(dp), intent(in) :: z_weight(4)
      integer :: slot(4), j, k

      slot = tbl%table_start_index(first_z:first_z+3) + ix
      do j = 1,n_opal95_t
         do k = 1,n_opal95_d
            tbl%fixed_z_opacity(ix,j,k) = z_weight(1)*tbl%full_opacity(slot(1),j,k)+ &
                 z_weight(2)*tbl%full_opacity(slot(2),j,k)+z_weight(3)*tbl%full_opacity(slot(3),j,k)+ &
                 z_weight(4)*tbl%full_opacity(slot(4),j,k)
         end do
      end do
end subroutine interpolate_x_table_in_z

!----------------------------------------------------------------------
! build_opal95_surface_table
!----------------------------------------------------------------------
! Builds the surface table (in T and log R) at hydrogen_fraction by
! 4-point Lagrange interpolation in X of the fixed-Z table set
! (tbl%fixed_z_opacity), from the X tables that exist at its Z, each at
! its own X (opal95_x_tables_at_z: the 3D and 4D paths' rule).
subroutine build_opal95_surface_table(tbl, hydrogen_fraction)
      type(opal95_table_set), intent(inout) :: tbl


      real(dp), intent(in) :: hydrogen_fraction

      real(dp) :: table_x_nodes(4), x_weight(4)
      integer :: i, j, k, first_table, table_index(4)

      tbl%surface_x = hydrogen_fraction
!  The first of the four consecutive X tables around X: the first whose
!  third node reaches X, the top four above the next-to-last node.
      first_table = n_opal95_x - 3
      do i = 3,n_opal95_x-1
         if (tbl%grid_x(i) >= hydrogen_fraction) then
            first_table = i - 2
            exit
         end if
      end do
!  The four tables that exist at the fixed Z, and their X.
      call opal95_x_tables_at_z(tbl%grid_x, tbl%fixed_z, hydrogen_fraction, first_table, &
           table_index, table_x_nodes)
      tbl%surface_x_tables = table_index
!
!  Find interpolation factors in X.
      call interp(table_x_nodes, x_weight, x_eval=hydrogen_fraction)
      do j = 1,n_opal95_t
         do k = 1,n_opal95_d
            tbl%surface_opacity(j,k) = x_weight(1)*tbl%fixed_z_opacity(table_index(1),j,k) + &
                 x_weight(2)*tbl%fixed_z_opacity(table_index(2),j,k) + &
                 x_weight(3)*tbl%fixed_z_opacity(table_index(3),j,k) + &
                 x_weight(4)*tbl%fixed_z_opacity(table_index(4),j,k)
         end do
      end do
end subroutine build_opal95_surface_table

!----------------------------------------------------------------------
! compute_opal95_opacity
!----------------------------------------------------------------------
! The opacity at (log rho, log T, X, Z) and its log derivatives, by
! 4-point Lagrange interpolation, in five steps on the call's stencil st
! (its first indices are tbl's warm-start cursors):
!   1. locate_opal95_t_and_logr: the stencil in T, with its weights, and
!      the first column of the stencil in log R;
!   2. locate_opal95_composition: which table the point is interpolated
!      in, and its stencils in X and Z:
!        2D (T, log R) in the surface table, at the surface X and the
!           fixed Z (each within opal95_composition_tol);
!        3D (X, T, log R) in the fixed-Z table, at the fixed Z and
!           another X;
!        4D (Z, X, T, log R) in the full set, at any other Z;
!   3. shift_rows_off_empty_corner: where the stencil reaches the tables'
!      empty high-R corner, each T row's columns moved down to the row's
!      last filled ones, and the rows the point lies beyond marked for
!      extrapolation;
!   4. opal95_logr_weights: each T row's weights in log R;
!   5. the interpolation on the chosen path (interpolate_opal95_2d, _3d,
!      _4d).
! A log T off the table returns ierr = 1 with the reason in message.
! Below the table's first log R column the opacity is extrapolated
! linearly in log R (extrapolate_in_logr), without a limit.
subroutine compute_opal95_opacity(tbl, log10_density, log10_temperature, hydrogen_fraction, &
     metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, message, ierr)
      type(opal95_table_set), intent(inout) :: tbl

      real(dp), intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

! the call's stencil: its weights and where log R is extrapolated from
! (its first indices are tbl's cursors)
      type(opal95_stencil) :: st
      integer :: path, corner_x_tables(4)
      logical :: rows_shifted

      ierr = 0
      message = ''
      if (log10_temperature < tbl%grid_logt(1) .or. &
           log10_temperature > tbl%grid_logt(n_opal95_t)) then
         write(message,'(a,f9.6,a,2(f5.2,a))') 'log T =', log10_temperature, &
              ' is outside the OPAL 1995 tables (log T', tbl%grid_logt(1), ' to', &
              tbl%grid_logt(n_opal95_t), ')'
         ierr = 1
         return
      end if
      call locate_opal95_t_and_logr(tbl, log10_density, log10_temperature, st)
      call locate_opal95_composition(tbl, hydrogen_fraction, metal_fraction, st, path, &
           corner_x_tables)
      call shift_rows_off_empty_corner(tbl, hydrogen_fraction, corner_x_tables, st, rows_shifted)
      call opal95_logr_weights(tbl, rows_shifted, st)
      select case (path)
      case (path_surface_table)
         call interpolate_opal95_2d(tbl, st, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
      case (path_fixed_z_table)
         call interpolate_opal95_3d(tbl, st, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
      case default
         call interpolate_opal95_4d(tbl, st, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)
      end select
end subroutine compute_opal95_opacity

! Step 1: log R = log rho - 3 log T6, and whether it is below the
! table's first column (extrap_lo); the stencil in T (tbl%index_t) with
! its weights and their T derivatives; and the first column of the
! stencil in log R (tbl%index_rho(1)), which step 3 may move row by row.
subroutine locate_opal95_t_and_logr(tbl, log10_density, log10_temperature, st)
      type(opal95_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: log10_density, log10_temperature
      type(opal95_stencil), intent(inout) :: st
      real(dp) :: interp_nodes(4)

      st%logr = log10_density - 3.0d0*(log10_temperature-6.0d0)
      st%extrap_lo = st%logr < tbl%grid_logr(1)
      call stencil4_locate_opal95(tbl%grid_logt, n_opal95_t, log10_temperature, tbl%index_t)
      interp_nodes = tbl%grid_logt(tbl%index_t:tbl%index_t+3)
      call interp(interp_nodes, st%weight_t, st%dweight_t, log10_temperature)
      call stencil4_locate_opal95(tbl%grid_logr, n_opal95_d, st%logr, tbl%index_rho(1))
end subroutine locate_opal95_t_and_logr

! Step 2: the table the point is interpolated in (path: the surface
! table, the fixed-Z table or the full set), with the point's stencil in
! X (and, on the 4D path, in Z, one X stencil per Z row) and its
! weights; and x_tables, the X tables of the fixed-Z table around the
! point's X, which step 3 reads: the tables the surface table was built
! from, the 3D path's stencil, or the 4D path's at its first Z row. The
! 2D path locates no X stencil, so it leaves the X cursors alone.
subroutine locate_opal95_composition(tbl, hydrogen_fraction, metal_fraction, st, path, x_tables)
      type(opal95_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      type(opal95_stencil), intent(inout) :: st
      integer, intent(out) :: path, x_tables(4)
      logical :: at_fixed_z, at_surface_x
      integer :: k

      at_fixed_z = abs(metal_fraction-tbl%fixed_z)/max(tbl%fixed_z,fixed_z_floor) <= &
           opal95_composition_tol
      at_surface_x = at_fixed_z .and. abs(hydrogen_fraction-tbl%surface_x) <= opal95_composition_tol
      if (at_surface_x) then
         path = path_surface_table
         x_tables = tbl%surface_x_tables
      else if (at_fixed_z) then
         path = path_fixed_z_table
         call locate_opal95_x_stencil(tbl, 1, tbl%fixed_z, hydrogen_fraction, st)
         x_tables = tbl%index_x(1,:)
      else
         path = path_full_set
         call locate_opal95_z_stencil(tbl, metal_fraction, st)
         do k = 1,4
            call locate_opal95_x_stencil(tbl, k, tbl%grid_z(tbl%index_z+k-1), hydrogen_fraction, st)
         end do
         x_tables = tbl%index_x(1,:)
      end if
end subroutine locate_opal95_composition

! The 4D path's stencil in Z: the warm-start search moves tbl%index_z
! from its previous value to the first of the four Z tables around
! metal_fraction (the last four above the next-to-last node), and
! st%weight_z takes their weights.
subroutine locate_opal95_z_stencil(tbl, metal_fraction, st)
      type(opal95_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: metal_fraction
      type(opal95_stencil), intent(inout) :: st
      real(dp) :: interp_nodes(4)
      integer :: i

      if (metal_fraction > tbl%grid_z(tbl%index_z+2)) then
         do i = tbl%index_z+3,n_opal95_z-1
            if (metal_fraction < tbl%grid_z(i)) then
               tbl%index_z = i - 2
               exit
            end if
         end do
         if (i > n_opal95_z-1) tbl%index_z = n_opal95_z - 3
      else
         do i = tbl%index_z+1,2,-1
            if (metal_fraction >= tbl%grid_z(i)) then
               tbl%index_z = i - 1
               exit
            end if
         end do
         if (i < 2) tbl%index_z = 1
      end if
!     (the last stencil at most)
      tbl%index_z = min(tbl%index_z,n_opal95_z-3)
      interp_nodes = tbl%grid_z(tbl%index_z:tbl%index_z+3)
      call interp(interp_nodes, st%weight_z, x_eval=metal_fraction)
end subroutine locate_opal95_z_stencil

! The stencil in X of Z row k (at metal fraction z_here): the warm-start
! search moves tbl%index_x(k,1) from its previous value to the first of
! the four consecutive X tables around hydrogen_fraction;
! opal95_x_tables_at_z then takes the tables that exist at this Z
! (tbl%index_x(k,1:4)), and st%weight_x(1:4,k) takes their weights. The
! 3D path calls it for k = 1 at the fixed Z, the 4D path for each of its
! Z rows.
subroutine locate_opal95_x_stencil(tbl, k, z_here, hydrogen_fraction, st)
      type(opal95_table_set), intent(inout) :: tbl
      integer, intent(in) :: k
      real(dp), intent(in) :: z_here, hydrogen_fraction
      type(opal95_stencil), intent(inout) :: st
      integer :: i, tables_at_z(4)
      real(dp) :: interp_nodes(4)

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
      call opal95_x_tables_at_z(tbl%grid_x, z_here, hydrogen_fraction, tbl%index_x(k,1), &
           tables_at_z, interp_nodes)
      tbl%index_x(k,:) = tables_at_z
      call interp(interp_nodes, st%weight_x(:,k), x_eval=hydrogen_fraction)
end subroutine locate_opal95_x_stencil

! Step 3: the tables have no data in their high-R, high-T corner (each T
! row from log T 7.2 up ends before the last log R column), and a point
! inside the table can still have part of its 4 x 4 stencil in T and
! log R there. The corner is read in the fixed-Z table, in each of the
! X tables x_tables: a cell empty in any of them counts as empty.
! Where the stencil's top cell (its last T row, its last column) is
! empty, rows_shifted is set and each T row gets its own columns
! (tbl%index_rho(i)): the last four filled in every X table. A row whose
! second column is empty too has the point beyond its last filled column
! by more than a column: it is marked (st%extrap_hi_row) to be evaluated
! at that column and extrapolated linearly in log R (extrapolate_in_logr).
! The low-R, low-T corner at low X (X below low_corner_x, the stencil's
! first T row at or below low_corner_t_row) is not checked: its empty
! region depends on Z, and the check the F77 code postponed was never
! written, so there the stencil is used as found.
subroutine shift_rows_off_empty_corner(tbl, hydrogen_fraction, x_tables, st, rows_shifted)
      type(opal95_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: hydrogen_fraction
      integer, intent(in) :: x_tables(4)
      type(opal95_stencil), intent(inout) :: st
      logical, intent(out) :: rows_shifted
      integer :: i, j, k, first_column, t_row

      rows_shifted = .false.
      st%extrap_hi = .false.
      st%extrap_hi_row = .false.
      if (hydrogen_fraction < low_corner_x .and. tbl%index_t <= low_corner_t_row) return
!     the stencil's top cell (the point leaves the data first at the
!     highest T)
      if (.not. any(tbl%fixed_z_opacity(x_tables,tbl%index_t+3, &
           tbl%index_rho(1)+3) >= opal95_missing_test)) return
      rows_shifted = .true.
!     the rows the point lies beyond by more than a column
      do k = 1,4
         if (any(tbl%fixed_z_opacity(x_tables,tbl%index_t+k-1, &
              tbl%index_rho(1)+1) >= opal95_missing_test)) then
            st%extrap_hi = .true.
            st%extrap_hi_row(k) = .true.
         end if
      end do
!     each row's four columns: those found for the row before it, or, if
!     their last cell is empty on this row, the row's last four filled
      first_column = tbl%index_rho(1)
      do i = 1,4
         t_row = tbl%index_t+i-1
         if (any(tbl%fixed_z_opacity(x_tables,t_row,first_column+3) > opal95_missing_test)) then
            do j = first_column+2,1,-1
               if (all(tbl%fixed_z_opacity(x_tables,t_row,j) <= opal95_missing_test)) then
                  tbl%index_rho(i) = j - 3
                  first_column = j - 3
                  exit
               end if
            end do
         else
            tbl%index_rho(i) = first_column
         end if
      end do
end subroutine shift_rows_off_empty_corner

! Step 4: each T row's weights in log R (st%weight_rho(row,:)) and their
! log R derivatives, at the point's log R on the row's four columns.
! Without a shift every row has the first row's columns, and below the
! table's first column the weights are taken at that column
! (st%logr_lo_edge). With a shift each row has its own columns, and a
! row marked for extrapolation takes its weights at its last column
! (st%logr_hi_edge(row)). extrapolate_in_logr carries each such row from
! its edge to the point.
subroutine opal95_logr_weights(tbl, rows_shifted, st)
      type(opal95_table_set), intent(inout) :: tbl
      logical, intent(in) :: rows_shifted
      type(opal95_stencil), intent(inout) :: st
      real(dp) :: interp_nodes(4), weight(4), dweight(4)
      integer :: j

      if (.not. rows_shifted) then
         interp_nodes = tbl%grid_logr(tbl%index_rho(1):tbl%index_rho(1)+3)
         if (.not.st%extrap_lo) then
            call interp(interp_nodes, weight, dweight, st%logr)
         else
            st%logr_lo_edge = tbl%grid_logr(1)
            call interp(interp_nodes, weight, dweight, st%logr_lo_edge)
         end if
         do j = 1,4
            tbl%index_rho(j) = tbl%index_rho(1)
            st%weight_rho(j,:) = weight
            st%dweight_rho(j,:) = dweight
         end do
      else
         do j = 1,4
            interp_nodes = tbl%grid_logr(tbl%index_rho(j):tbl%index_rho(j)+3)
            if (.not.st%extrap_hi_row(j)) then
               call interp(interp_nodes, weight, dweight, st%logr)
            else
               st%logr_hi_edge(j) = tbl%grid_logr(tbl%index_rho(j)+3)
               call interp(interp_nodes, weight, dweight, st%logr_hi_edge(j))
            end if
            st%weight_rho(j,:) = weight
            st%dweight_rho(j,:) = dweight
         end do
      end if
end subroutine opal95_logr_weights

!----------------------------------------------------------------------
! opal95_x_tables_at_z: the X stencil at one Z
!----------------------------------------------------------------------
! The four X tables an interpolation in X uses at metal fraction z_here
! (their indices in table_index) and the X each is at (nodes), given
! first_table, the first of the four consecutive tables around
! hydrogen_fraction that a search of grid_x found. The set lacks the
! X = 0.95 table (9) from Z = 0.04 up and the X = 1 - Z table (10) at
! Z = 0.1, so there a stencil at X >= 0.8 takes the top tables that
! exist: 6, 7, 8 and 10, or 5 to 8 at Z = 0.1. The tenth table is at
! X = 1 - Z, not at grid_x's 1. The surface table's build and the 3D
! and 4D paths' stencils (locate_opal95_x_stencil) share this rule.
pure subroutine opal95_x_tables_at_z(grid_x, z_here, hydrogen_fraction, first_table, &
     table_index, nodes)
      real(dp), intent(in) :: grid_x(:), z_here, hydrogen_fraction
      integer, intent(in) :: first_table
      integer, intent(out) :: table_index(4)
      real(dp), intent(out) :: nodes(4)
      integer :: i

      if (z_here >= no_x095_from_z .and. hydrogen_fraction >= high_x) then
         if (z_here < no_x1mz_from_z) then
!           no table 9 (X = 0.95)
            table_index = [n_opal95_x - 4, n_opal95_x - 3, n_opal95_x - 2, n_opal95_x]
         else
!           no table 9 or table 10
            table_index = [n_opal95_x - 5, n_opal95_x - 4, n_opal95_x - 3, n_opal95_x - 2]
         end if
      else
         table_index = [(first_table + i - 1, i = 1, 4)]
      end if
      nodes = grid_x(table_index)
      if (table_index(4) == n_opal95_x) nodes(4) = 1.0d0 - z_here
end subroutine opal95_x_tables_at_z

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

subroutine load_opal95(settings, paths, envelope_hydrogen_fraction, tables, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(kap_table_paths_t), intent(in) :: paths
      real(dp), intent(in) :: envelope_hydrogen_fraction
      type(opacity_table_state), intent(inout) :: tables
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

! (the reader reports its own errors)
      message = ''
      call read_opal95_tables(tables%opal95, paths%opal95, ierr)
      if (ierr /= 0) return
      call build_opal95_fixed_z_table(tables%opal95, settings%opal95_single_table_z, message, ierr)
      if (ierr /= 0) return
      call build_opal95_surface_table(tables%opal95, envelope_hydrogen_fraction)
end subroutine load_opal95

subroutine refresh_opal95(settings, hydrogen_fraction, tables, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      real(dp), intent(in) :: hydrogen_fraction
      type(opacity_table_state), intent(inout) :: tables
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

      ierr = 0
      message = ''
      call build_opal95_surface_table(tables%opal95, hydrogen_fraction)
end subroutine refresh_opal95

subroutine lookup_opal95(settings, tables, point, res, in_table, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(opacity_table_state), intent(inout) :: tables
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      logical, intent(out) :: in_table
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

      in_table = .true.
      call compute_opal95_opacity(tables%opal95, point%log10_density, &
           point%log10_temperature, point%hydrogen_fraction, point%metal_fraction, &
           res%opacity, res%log10_opacity, res%dlnkap_dlnrho, res%dlnkap_dlnt, message, ierr)
end subroutine lookup_opal95

! in a helium-burning core: the tables at their top Z
subroutine core_lookup_opal95(settings, tables, point, res, table_z, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(opacity_table_state), intent(inout) :: tables
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      real(dp), intent(out) :: table_z
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

      table_z = tables%opal95%grid_z(n_opal95_z)
      call compute_opal95_opacity(tables%opal95, point%log10_density, &
           point%log10_temperature, point%hydrogen_fraction, table_z, &
           res%opacity, res%log10_opacity, res%dlnkap_dlnrho, res%dlnkap_dlnt, message, ierr)
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
! sharing one row helper and one log R extrapolation; the high-R shift's
! search in the fixed-Z table, where the F77 code searched the full set's
! Z = 0 slot (round 1B, M277), and in the X tables of the point's own
! stencil, where it read the X index the previous lookup had left
! (round 3, M298); the surface table's X tables those of the lookups'
! stencil (M339).

end module opal95_opacity
