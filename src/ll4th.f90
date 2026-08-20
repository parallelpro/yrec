!----------------------------------------------------------------------
! ll4th
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ll4th.f; only variable names, source form, and comment
! style were updated.
!
! THIS IS THE 4TH TABLE OF THE LAWRENCE LIVERMORE OPACITY TABLES.
! IT CONTAINS THE OPACITIES FOR THE SURFACE ABUNDANCES OF X AND Z, AND
! IS USED TO AVOID MANY INTERPOLATIONS TO THE SAME X.
subroutine ll4th(hydrogen_fraction)

      implicit none
      integer, parameter :: num_t = 50
      integer, parameter :: num_d = 17
      integer, parameter :: num_x = 3
      integer, parameter :: num_xt = num_t*num_x
      integer, parameter :: num_4d = 4*num_d
! MHP 10/02 made array dimensions consistent
      integer, parameter :: np = 100

      double precision, intent(in) :: hydrogen_fraction

      double precision :: coeff(4,np)
!     OUTPUT OF LL4TH
      double precision :: opal92_surface_x, opal92_surface_z, &
           opal92_surface_spline_coeffs(num_t,num_4d)
      integer :: opal92_surface_x_index
      common /llot4/opal92_surface_x, opal92_surface_z, &
           opal92_surface_spline_coeffs, opal92_surface_x_index
      double precision :: opal92_surface_x_z2, opal92_surface_z_z2, &
           opal92_surface_spline_coeffs_z2(num_t,num_4d)
      integer :: opal92_surface_x_index_z2
      common /llot42/opal92_surface_x_z2, opal92_surface_z_z2, &
           opal92_surface_spline_coeffs_z2, opal92_surface_x_index_z2
!     GRIDS
      double precision :: opal92_grid_logt(num_t), opal92_grid_x(num_x), &
           opal92_grid_logr(num_d)
      common /gllot/ opal92_grid_logt, opal92_grid_x, opal92_grid_logr
      double precision :: opal92_grid_logt_z2(num_t), opal92_grid_x_z2(num_x), &
           opal92_grid_logr_z2(num_d)
      common /gllot2/ opal92_grid_logt_z2, opal92_grid_x_z2, opal92_grid_logr_z2
!     LL OPACITY TABLES
      double precision :: opal92_log10_opacity(num_xt, num_d)
      integer :: opal92_num_x, opal92_num_temps
      common /llot/ opal92_log10_opacity, opal92_num_x, opal92_num_temps
      double precision :: opal92_log10_opacity_z2(num_xt, num_d)
      integer :: opal92_num_x_z2, opal92_num_temps_z2
      common /llot2/ opal92_log10_opacity_z2, opal92_num_x_z2, opal92_num_temps_z2
!     THE SPLINE COEFFICIENTS AND THE SHAPE OF THE LL TABLE.
      double precision :: opal92_spline_coeffs(num_xt, num_4d)
      integer :: opal92_density_start_index(num_xt), opal92_density_count(num_xt)
      common /lintpl/ opal92_spline_coeffs, opal92_density_start_index, &
           opal92_density_count
      double precision :: opal92_spline_coeffs_z2(num_xt, num_4d)
      integer :: opal92_density_start_index_z2(num_xt), opal92_density_count_z2(num_xt)
      common /lintpl2/ opal92_spline_coeffs_z2, opal92_density_start_index_z2, &
           opal92_density_count_z2
! OPACITY COMMON BLOCKS - modified 3/09
! common/newopac/: only use_two_z_tables is used here.
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

      integer :: x_index, temp_index, rho_search_index, im2, im3, &
           row_index, coeff_base_index, density_start, density_end, j, i
      double precision :: temp6, density_rhot3, opacity0, log10_opacity0, &
           dlnkap_dlnrho0, dlnkap_dlnt0, opacity1, log10_opacity1, &
           dlnkap_dlnrho1, dlnkap_dlnt1, x_fraction_within, log10_opacity_final, &
           opacity_final

!     KEEP THE COMPOSITION OF THE 4TH TABLE.
      opal92_surface_x=hydrogen_fraction
      opal92_surface_z=opal_table_z1
      call findex(opal92_grid_x, num_x, hydrogen_fraction, x_index)
      if (x_index.lt.0) x_index=-x_index
      if (x_index.ge.3) x_index=2
      if (x_index.le.0) stop ' ERROR IN X GRID'
      opal92_surface_x_index=x_index
      do im2=1,opal92_num_temps
         x_index=opal92_surface_x_index
         temp6=opal92_grid_logt(im2)
         temp_index=im2
         row_index=im2+(x_index-1)*opal92_num_temps
         density_start=opal92_density_start_index(row_index)
         if (density_start.ne.1) stop ' LL4TH NDSS '
         density_end=density_start+opal92_density_count(row_index)-1
         do im3=density_start,density_end
            density_rhot3=opal92_grid_logr(im3)
            call yllo2d(temp6,density_rhot3,x_index,temp_index,rho_search_index, &
                 opacity0,log10_opacity0,dlnkap_dlnrho0,dlnkap_dlnt0)
            call yllo2d(temp6,density_rhot3,x_index+1,temp_index,rho_search_index, &
                 opacity1,log10_opacity1,dlnkap_dlnrho1,dlnkap_dlnt1)
            x_fraction_within=(hydrogen_fraction-opal92_grid_x(x_index))/ &
                 (opal92_grid_x(x_index+1)-opal92_grid_x(x_index))
            log10_opacity_final=(log10_opacity1-log10_opacity0)*x_fraction_within+log10_opacity0
            opacity_final=10.0d0**log10_opacity_final
!     CONVERSION FROM THE DERIVATIVE WITH CONSTANT RHOT3 TO CONSTANT RHO
            coeff(1,im3)=log10_opacity_final
         end do
         call ysplin(opal92_grid_logr,coeff,density_end)
         do j=1,density_end
            do i=1,4
               coeff_base_index=i+(j-1)*4
               opal92_surface_spline_coeffs(im2,coeff_base_index)=coeff(i,j)
            end do
         end do
      end do



      if (use_two_z_tables) then
         opal92_surface_x_z2=hydrogen_fraction
         opal92_surface_z_z2=opal_table_z2
         call findex(opal92_grid_x_z2, num_x, hydrogen_fraction, x_index)
         if (x_index.lt.0) x_index=-x_index
         if (x_index.ge.3) x_index=2
         if (x_index.le.0) stop ' ERROR IN X GRID'
         opal92_surface_x_index_z2=x_index
         do im2=1,opal92_num_temps_z2
            x_index=opal92_surface_x_index_z2
            temp6=opal92_grid_logt_z2(im2)
            temp_index=im2
            row_index=im2+(x_index-1)*opal92_num_temps_z2
            density_start=opal92_density_start_index_z2(row_index)
            if (density_start.ne.1) stop ' LL4TH Z1 CHECK NDSS '
            density_end=density_start+opal92_density_count_z2(row_index)-1
            do im3=density_start,density_end
               density_rhot3=opal92_grid_logr_z2(im3)
               call yllo2d2(temp6,density_rhot3,x_index,temp_index,rho_search_index, &
                    opacity0,log10_opacity0,dlnkap_dlnrho0,dlnkap_dlnt0)
               call yllo2d2(temp6,density_rhot3,x_index+1,temp_index,rho_search_index, &
                    opacity1,log10_opacity1,dlnkap_dlnrho1,dlnkap_dlnt1)
               x_fraction_within=(hydrogen_fraction-opal92_grid_x_z2(x_index))/ &
                    (opal92_grid_x_z2(x_index+1)-opal92_grid_x_z2(x_index))
               log10_opacity_final=(log10_opacity1-log10_opacity0)*x_fraction_within+log10_opacity0
               opacity_final=10.0d0**log10_opacity_final
! CONVERSION FROM THE DERIVATIVE WITH CONSTANT RHOT3 TO CONSTANT RHO
               coeff(1,im3)=log10_opacity_final
            end do
            call ysplin(opal92_grid_logr_z2,coeff,density_end)
            do j=1,density_end
               do i=1,4
                  coeff_base_index=i+(j-1)*4
                  opal92_surface_spline_coeffs_z2(im2,coeff_base_index)=coeff(i,j)
               end do
            end do
         end do
      end if



      return
end subroutine ll4th
