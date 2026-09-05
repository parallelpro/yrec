!----------------------------------------------------------------------
! opal92_table_prep
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ylloc.f; only variable names, source form, and comment
! style were updated.
!
! DBG 5/94 Modified to include ZRAMP and ZDIFF stuff (second opacity
! table). Builds the cubic-spline coefficients (in density) for the
! OPAL92 opacity tables, and for the second (different-Z) OPAL92
! table set when use_two_z_tables is set. Note the two strides for
! the same logical row: index1 uses the number of temperature rows
! actually read, while log10_opacity is addressed with the fixed
! n_opal92_t stride that read_opal92_tables stores with.
!
! 2026 wave 3 (R5): the verbatim second-Z loop nest (on
! opacity_table%opal92(2)) became the internal prep_table below,
! called once per set in the original order; the loop's working
! scalars stay host-associated. The stdout diagnostic no longer says
! "second".
subroutine opal92_table_prep(ierr)
      use star_info_lib, only: star
      use opacity_table_lib
      use numerics_lib
      implicit none
      integer, intent(out) :: ierr
! MHP 10/02 made array dimensions consistent
      integer, parameter :: np = 100

      double precision :: spline_work(4,np), density_nodes(n_opal92_d)
      integer :: ix, it, index1, jd, id, index2, j, i
      double precision :: chkd, chko

      ierr = 0

      call prep_table(opacity_table%opal92(1), ierr)
      if (ierr /= 0) return
!
! DBG 5/94 ZRAMP stuff
      if (star%use_two_z_tables) then
         call prep_table(opacity_table%opal92(2), ierr)
      end if

      return

contains

! prep_table: density-spline coefficients for every (X,T) row of tbl.
subroutine prep_table(tbl, ierr)
      type(opal92_table_set), intent(inout) :: tbl
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success

      do ix = 1,tbl%num_x
      do it = 1,tbl%num_temps
       index1 = it + (ix-1)*tbl%num_temps
       jd = 0
        do id = 1,n_opal92_d
          chkd = tbl%grid_logr(id)
          chko = tbl%log10_opacity(it+n_opal92_t*(ix-1),id)
!>>>> CHECK THE EMPTY REGION
          if (chko.le.opal92_missing_opacity) cycle
          if (jd.le.0) then
             tbl%density_start_index(index1) = id
             if (id.ne.1) then
                write(*,*) 'opal92_table_prep: opal92 table density grid does not start at index 1'
                ierr = 1
                return
             end if
          endif
          jd = jd + 1
          density_nodes(jd) = chkd
!>>>> CHECK THE OPACITY VALUE IN THE TABLE
          spline_work(1,jd) = chko
        end do
        tbl%density_count(index1) = jd
        if (jd.le.1) cycle
        call ysplin(density_nodes, spline_work, jd)
        do j = 1,jd
        do i = 1,4
         index2 = i + (j-1)*4
         tbl%spline_coeffs(index1,index2) = spline_work(i,j)
        end do
        end do
      end do
      end do
end subroutine prep_table

end subroutine opal92_table_prep
