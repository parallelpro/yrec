!----------------------------------------------------------------------
! build_kurucz_splines
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ykoeff.f; only variable names, source form, and comment
! style were updated.
!
! YCK 3/91. Builds the cubic-spline coefficients (in density) for
! the Kurucz90 opacity tables, and for the second (different-Z)
! Kurucz table set when use_two_z_tables is set. Reads the
! opacity_table%kurucz_*/kurucz2_* grids filled by
! read_kurucz_tables.f90 and writes the per-row spline coefficients,
! density start index and count used by kurucz.f90/kurucz2.f90.
subroutine build_kurucz_splines(ierr)
      use star_info_lib, only: star

      use opacity_table_lib
      use numerics_lib
      use math_lib
      implicit none
      integer, intent(out) :: ierr
! MHP 10/02 made array dimensions consistent
      integer, parameter :: np = 100

      double precision :: spline_work(4,np), density_nodes(kurucz_max_num_densities)
      integer :: it, jd, id, index2, j, i
      double precision :: chkd, chko

      ierr = 0
      do it = 1,opacity_table%kurucz_num_temps
         jd = 0
         do id = 1,kurucz_max_num_densities
            chkd = opacity_table%kurucz_log10_rho(it,id)
            chko = opacity_table%kurucz_log10_opacity(it,id)
            if (chko.le.0.0d0) cycle
            if (jd.le.0) opacity_table%kurucz_density_start_index(it) = id
            jd = jd + 1
            density_nodes(jd) = chkd
            spline_work(1,jd) = log10(chko)
         end do
         opacity_table%kurucz_density_count(it) = jd
         if (opacity_table%kurucz_density_start_index(it).ne.1) then
            write(*,*) 'build_kurucz_splines: kurucz table density grid does not start at index 1'
            ierr = 1
            return
         end if
         if (opacity_table%kurucz_density_count(it).lt.25) then
            write(*,*) 'build_kurucz_splines: kurucz table has fewer than 25 density points'
            ierr = 1
            return
         end if
         call ysplin(density_nodes, spline_work, jd)
         do j = 1,jd
            do i = 1,4
               index2 = i + (j-1)*4
               opacity_table%kurucz_spline_coeffs(it,index2) = spline_work(i,j)
            end do
         end do
      end do
! DBG 12/95 second Z table
      if (star%use_two_z_tables) then
         do it = 1,opacity_table%kurucz2_num_temps
            jd = 0
            do id = 1,kurucz_max_num_densities
               chkd = opacity_table%kurucz2_log10_rho(it,id)
               chko = opacity_table%kurucz2_log10_opacity(it,id)
               if (chko.le.0.0d0) cycle
               if (jd.le.0) opacity_table%kurucz2_density_start_index(it) = id
               jd = jd + 1
               density_nodes(jd) = chkd
               spline_work(1,jd) = log10(chko)
            end do
            opacity_table%kurucz2_density_count(it) = jd
            if (opacity_table%kurucz2_density_start_index(it).ne.1) then
               write(*,*) 'build_kurucz_splines: second kurucz table density grid does not start at index 1'
               ierr = 1
               return
            end if
            if (opacity_table%kurucz2_density_count(it).lt.25) then
               write(*,*) 'build_kurucz_splines: second kurucz table has fewer than 25 density points'
               ierr = 1
               return
            end if
            call ysplin(density_nodes, spline_work, jd)
            do j = 1,jd
               do i = 1,4
                  index2 = i + (j-1)*4
                  opacity_table%kurucz2_spline_coeffs(it,index2) = spline_work(i,j)
               end do
            end do
         end do
      end if

      return
end subroutine build_kurucz_splines
