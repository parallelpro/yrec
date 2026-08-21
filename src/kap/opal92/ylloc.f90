!----------------------------------------------------------------------
! ylloc
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ylloc.f; only variable names, source form, and comment
! style were updated. common/gllot/, common/llot/, and
! common/lintpl/ member names match those established in
! ll4th.f90/setllo.f90.
!
! DBG 5/94 Modified to include ZRAMP and ZDIFF stuff (second opacity
! table). Builds the cubic-spline coefficients (in density) for the
! OPAL92 opacity tables, and for the second (different-Z) OPAL92
! table set when use_two_z_tables is set.
subroutine ylloc
      use opacity_table_lib
      use const_lib
      use numerics_lib
      implicit none
      integer, parameter :: num_t = 50
      integer, parameter :: num_d = 17
      integer, parameter :: num_x = 3
      integer, parameter :: num_xt = num_t*num_x
      integer, parameter :: num_4d = 4*num_d
! MHP 10/02 made array dimensions consistent
      integer, parameter :: np = 100

      double precision :: spline_work(4,np), density_nodes(num_d)
      save

      integer :: ix, it, index1, jd, ids, idf, id, index2, j, i
      double precision :: chkd, chko

      do 1 ix = 1,opacity_table%opal92_num_x
      do 2 it = 1,opacity_table%opal92_num_temps
       index1 = it + (ix-1)*opacity_table%opal92_num_temps
       jd = 0
        ids = 1
        idf = num_d
        do 3 id = ids,idf
          chkd = opacity_table%opal92_grid_logr(id)
          chko = opacity_table%opal92_log10_opacity(it+num_t*(ix-1),id)
!>>>> CHECK THE EMPTY REGION
          if (chko.le.-9.999d0) go to 3
          if (jd.le.0) then
             opacity_table%opal92_density_start_index(index1) = id
             if (id.ne.1) stop ' CHECK NDS '
          endif
          jd = jd + 1
          density_nodes(jd) = chkd
!>>>> CHECK THE OPACITY VALUE IN THE TABLE
          spline_work(1,jd) = chko
    3   continue
        opacity_table%opal92_density_count(index1) = jd
        if (jd.le.1) go to 2
        call ysplin(density_nodes, spline_work, jd)
        do 100 j = 1,jd
        do 200 i = 1,4
         index2 = i + (j-1)*4
         opacity_table%opal92_spline_coeffs(index1,index2) = spline_work(i,j)
  200   continue
  100   continue
    2 continue
    1 continue
!
! DBG 5/94 ZRAMP stuff
      if (use_two_z_tables) then
         do ix = 1,opacity_table%opal92_num_x_z2
            do it = 1,opacity_table%opal92_num_temps_z2
               index1 = it + (ix-1)*opacity_table%opal92_num_temps_z2
               jd = 0
               ids = 1
               idf = num_d
               do id = ids,idf
                  chkd = opacity_table%opal92_grid_logr_z2(id)
                  chko = opacity_table%opal92_log10_opacity_z2(it+num_t*(ix-1),id)
                  if (chko.le.-9.999d0) go to 503
                  if (jd.le.0) then
                     opacity_table%opal92_density_start_index_z2(index1) = id
                     if (id.ne.1) stop ' CHECK NDS2 '
                  endif
                  jd = jd + 1
                  density_nodes(jd) = chkd
                  spline_work(1,jd) = chko
  503             continue
               end do
               opacity_table%opal92_density_count_z2(index1) = jd
               if (jd.le.1) go to 502
               call ysplin(density_nodes, spline_work, jd)
               do j = 1,jd
                  do i = 1,4
                     index2 = i + (j-1)*4
                     opacity_table%opal92_spline_coeffs_z2(index1,index2) = spline_work(i,j)
                  end do
               end do
  502          continue
            end do
         end do
      end if

      return
end subroutine ylloc
