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

      use ysplin_mod
      implicit none
      integer, parameter :: num_t = 60
      integer, parameter :: num_d = 50
      integer, parameter :: num_x = 1
      integer, parameter :: num_xt = num_t*num_x
      integer, parameter :: num_4d = 4*num_d
! MHP 10/02 made array dimensions consistent
      integer, parameter :: np = 100

      double precision :: spline_work(4,np), density_nodes(num_d)
!      DIMENSION C(4,NUMD),XD(NUMD)
      double precision :: kurucz_grid_logt(num_t)
      common /gkrz/kurucz_grid_logt
      double precision :: kurucz_log10_opacity(num_xt,num_d), kurucz_log10_rho(num_xt,num_d)
      integer :: kurucz_num_temps
      common /krz/kurucz_log10_opacity, kurucz_log10_rho, kurucz_num_temps
      double precision :: kurucz_spline_coeffs(num_xt,num_4d)
      integer :: kurucz_density_start_index(num_xt), kurucz_density_count(num_xt)
      common /intpl2/ kurucz_spline_coeffs, kurucz_density_start_index, &
           kurucz_density_count
      double precision :: kurucz2_grid_logt(num_t)
      common /gkrz2/kurucz2_grid_logt
      double precision :: kurucz2_log10_opacity(num_xt,num_d), kurucz2_log10_rho(num_xt,num_d)
      integer :: kurucz2_num_temps
      common /krz2/kurucz2_log10_opacity, kurucz2_log10_rho, kurucz2_num_temps
      double precision :: kurucz2_spline_coeffs(num_xt,num_4d)
      integer :: kurucz2_density_start_index(num_xt), kurucz2_density_count(num_xt)
      common /intpl22/ kurucz2_spline_coeffs, kurucz2_density_start_index, &
           kurucz2_density_count
! OPACITY COMMON BLOCKS - modified 3/09
      double precision :: laol_table_z1, laol_table_z2, opal_table_z1, &
           opal_table_z2, opal95_single_table_z, alex_table_z1, &
           kurucz_table_z1, kurucz_table_z2, molecular_opacity_logt_min, &
           molecular_opacity_logt_max
      logical :: use_alex06_tables, use_laol89_tables, use_opal92_tables, &
           use_opal95_tables, use_kurucz90_tables, use_alex95_tables, &
           use_two_z_tables
      common /newopac/ laol_table_z1, laol_table_z2, opal_table_z1, &
           opal_table_z2, opal95_single_table_z, alex_table_z1, &
           kurucz_table_z1, kurucz_table_z2, molecular_opacity_logt_min, &
           molecular_opacity_logt_max, use_alex06_tables, &
           use_laol89_tables, use_opal92_tables, use_opal95_tables, &
           use_kurucz90_tables, use_alex95_tables, use_two_z_tables
      save

      integer :: it, index1, jd, ids, idf, id, index2, j, i
      double precision :: chkd, chko

      do 102 it = 1,kurucz_num_temps
         index1 = it
         jd = 0
         ids = 1
         idf = num_d
         do 103 id = ids,idf
            chkd = kurucz_log10_rho(it,id)
            chko = kurucz_log10_opacity(it,id)
            if (chko.le.0.0d0) go to 103
            if (jd.le.0) kurucz_density_start_index(index1) = id
            jd = jd + 1
            density_nodes(jd) = chkd
            spline_work(1,jd) = dlog10(chko)
 103     continue
         kurucz_density_count(index1) = jd
         if (kurucz_density_start_index(index1).ne.1) stop ' ERROR KURUCZ OPACITY NDS'
         if (kurucz_density_count(index1).lt.25) stop ' ERROR KURUCZ OPACITY NDD'
         if (jd.le.1) go to 102
         call ysplin(density_nodes, spline_work, jd)
         do j = 1,jd
            do i = 1,4
               index2 = i + (j-1)*4
               kurucz_spline_coeffs(index1,index2) = spline_work(i,j)
            end do
         end do
 102  continue
!
!
!
! DBG 12/95 second Z table
      if (use_two_z_tables) then
         do 202 it = 1,kurucz2_num_temps
            index1 = it
            jd = 0
            ids = 1
            idf = num_d
            do 203 id = ids,idf
               chkd = kurucz2_log10_rho(it,id)
               chko = kurucz2_log10_opacity(it,id)
               if (chko.le.0.0d0) go to 203
               if (jd.le.0) kurucz2_density_start_index(index1) = id
               jd = jd + 1
               density_nodes(jd) = chkd
               spline_work(1,jd) = dlog10(chko)
 203        continue
            kurucz2_density_count(index1) = jd
            if (kurucz2_density_start_index(index1).ne.1) stop ' NDS2'
            if (kurucz2_density_count(index1).lt.25) stop ' NDD2'
            if (jd.le.1) go to 202
            call ysplin(density_nodes, spline_work, jd)
            do j = 1,jd
               do i = 1,4
                  index2 = i + (j-1)*4
                  kurucz2_spline_coeffs(index1,index2) = spline_work(i,j)
               end do
            end do
 202     continue
      end if

      return
end subroutine ykoeff
