!----------------------------------------------------------------------
! build_kurucz_splines
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ykoeff.f; only variable names, source form, and comment
! style were updated.
!
! YCK 3/91. Builds the cubic-spline coefficients (in density) for
! one Kurucz90 opacity table set tbl. Reads the tbl%grid_logt/
! log10_rho/log10_opacity grids filled by read_kurucz_tables.f90 and
! writes the per-row spline coefficients, density start index and
! count used by kurucz.f90.
!
! 2026 wave 3 (R5): was one routine that did this for
! opacity_table%kurucz_* and then (if use_two_z_tables) repeated the
! identical loop for kurucz2_*; read_kurucz_tables now calls it once
! per table set in the same order. The two stdout diagnostics that
! said "kurucz table" / "second kurucz table" are now generic.
subroutine build_kurucz_splines(tbl, ierr)
      use opacity_table_lib
      use numerics_lib
      use math_lib
      implicit none
      type(kurucz_table_set), intent(inout) :: tbl
      integer, intent(out) :: ierr
! MHP 10/02 made array dimensions consistent
      integer, parameter :: np = 100

      double precision :: spline_work(4,np), density_nodes(kurucz_max_num_densities)
      integer :: it, jd, id, index2, j, i
      double precision :: chkd, chko

      ierr = 0
      do it = 1,tbl%num_temps
         jd = 0
         do id = 1,kurucz_max_num_densities
            chkd = tbl%log10_rho(it,id)
            chko = tbl%log10_opacity(it,id)
            if (chko.le.0.0d0) cycle
            if (jd.le.0) tbl%density_start_index(it) = id
            jd = jd + 1
            density_nodes(jd) = chkd
            spline_work(1,jd) = log10(chko)
         end do
         tbl%density_count(it) = jd
         if (tbl%density_start_index(it).ne.1) then
            write(*,*) 'build_kurucz_splines: kurucz table density grid does not start at index 1'
            ierr = 1
            return
         end if
         if (tbl%density_count(it).lt.25) then
            write(*,*) 'build_kurucz_splines: kurucz table has fewer than 25 density points'
            ierr = 1
            return
         end if
         call ysplin(density_nodes, spline_work, jd)
         do j = 1,jd
            do i = 1,4
               index2 = i + (j-1)*4
               tbl%spline_coeffs(it,index2) = spline_work(i,j)
            end do
         end do
      end do

      return
end subroutine build_kurucz_splines
