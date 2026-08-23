!----------------------------------------------------------------------
! ykoeff
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ykoeff.f; only variable names, source form, and comment
! style were updated.
!
! YCK 3/91. Builds the cubic-spline coefficients (in density) for
! the Kurucz90 opacity tables, and for the second (different-Z)
! Kurucz table set when use_two_z_tables is set. common/krz/,
! common/gkrz/, common/intpl2/ (and their "2" siblings) match the
! names established in setkrz.f90/kurucz.f90/kurucz2.f90.
subroutine ykoeff

      use opacity_table_lib
      use const_lib
      use numerics_lib
      implicit none
      integer, parameter :: num_t = 60
      integer, parameter :: num_d = 50
      integer, parameter :: num_x = 1
      integer, parameter :: num_xt = num_t*num_x
      integer, parameter :: num_4d = 4*num_d
! MHP 10/02 made array dimensions consistent
      integer, parameter :: np = 100

      double precision :: spline_work(4,np), density_nodes(num_d)
      integer :: it, index1, jd, ids, idf, id, index2, j, i
      double precision :: chkd, chko

      do it = 1,opacity_table%kurucz_num_temps
         index1 = it
         jd = 0
         ids = 1
         idf = num_d
         do id = ids,idf
            chkd = opacity_table%kurucz_log10_rho(it,id)
            chko = opacity_table%kurucz_log10_opacity(it,id)
            if (chko.le.0.0d0) cycle
            if (jd.le.0) opacity_table%kurucz_density_start_index(index1) = id
            jd = jd + 1
            density_nodes(jd) = chkd
            spline_work(1,jd) = dlog10(chko)
         end do
         opacity_table%kurucz_density_count(index1) = jd
         if (opacity_table%kurucz_density_start_index(index1).ne.1) stop ' ERROR KURUCZ OPACITY NDS'
         if (opacity_table%kurucz_density_count(index1).lt.25) stop ' ERROR KURUCZ OPACITY NDD'
         if (jd.le.1) cycle
         call ysplin(density_nodes, spline_work, jd)
         do j = 1,jd
            do i = 1,4
               index2 = i + (j-1)*4
               opacity_table%kurucz_spline_coeffs(index1,index2) = spline_work(i,j)
            end do
         end do
      end do
!
!
!
! DBG 12/95 second Z table
      if (use_two_z_tables) then
         do it = 1,opacity_table%kurucz2_num_temps
            index1 = it
            jd = 0
            ids = 1
            idf = num_d
            do id = ids,idf
               chkd = opacity_table%kurucz2_log10_rho(it,id)
               chko = opacity_table%kurucz2_log10_opacity(it,id)
               if (chko.le.0.0d0) cycle
               if (jd.le.0) opacity_table%kurucz2_density_start_index(index1) = id
               jd = jd + 1
               density_nodes(jd) = chkd
               spline_work(1,jd) = dlog10(chko)
            end do
            opacity_table%kurucz2_density_count(index1) = jd
            if (opacity_table%kurucz2_density_start_index(index1).ne.1) stop ' NDS2'
            if (opacity_table%kurucz2_density_count(index1).lt.25) stop ' NDD2'
            if (jd.le.1) cycle
            call ysplin(density_nodes, spline_work, jd)
            do j = 1,jd
               do i = 1,4
                  index2 = i + (j-1)*4
                  opacity_table%kurucz2_spline_coeffs(index1,index2) = spline_work(i,j)
               end do
            end do
         end do
      end if

      return
end subroutine ykoeff
