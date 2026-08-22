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

      double precision, intent(in) :: hydrogen_fraction

      double precision :: coeff(4,np)
      integer :: x_index, temp_index, rho_search_index, im2, im3, &
           row_index, coeff_base_index, density_start, density_end, j, i
      double precision :: temp6, density_rhot3, opacity0, log10_opacity0, &
           dlnkap_dlnrho0, dlnkap_dlnt0, opacity1, log10_opacity1, &
           dlnkap_dlnrho1, dlnkap_dlnt1, x_fraction_within, log10_opacity_final, &
           opacity_final

!     KEEP THE COMPOSITION OF THE 4TH TABLE.
      opacity_table%opal92_surface_x=hydrogen_fraction
      opacity_table%opal92_surface_z=opal_table_z1
      call findex(opacity_table%opal92_grid_x, num_x, hydrogen_fraction, x_index)
      if (x_index.lt.0) x_index=-x_index
      if (x_index.ge.3) x_index=2
      if (x_index.le.0) stop ' ERROR IN X GRID'
      opacity_table%opal92_surface_x_index=x_index
      do im2=1,opacity_table%opal92_num_temps
         x_index=opacity_table%opal92_surface_x_index
         temp6=opacity_table%opal92_grid_logt(im2)
         temp_index=im2
         row_index=im2+(x_index-1)*opacity_table%opal92_num_temps
         density_start=opacity_table%opal92_density_start_index(row_index)
         if (density_start.ne.1) stop ' LL4TH NDSS '
         density_end=density_start+opacity_table%opal92_density_count(row_index)-1
         do im3=density_start,density_end
            density_rhot3=opacity_table%opal92_grid_logr(im3)
            call yllo2d(temp6,density_rhot3,x_index,temp_index,rho_search_index, &
                 opacity0,log10_opacity0,dlnkap_dlnrho0,dlnkap_dlnt0)
            call yllo2d(temp6,density_rhot3,x_index+1,temp_index,rho_search_index, &
                 opacity1,log10_opacity1,dlnkap_dlnrho1,dlnkap_dlnt1)
            x_fraction_within=(hydrogen_fraction-opacity_table%opal92_grid_x(x_index))/ &
                 (opacity_table%opal92_grid_x(x_index+1)-opacity_table%opal92_grid_x(x_index))
            log10_opacity_final=(log10_opacity1-log10_opacity0)*x_fraction_within+log10_opacity0
            opacity_final=10.0d0**log10_opacity_final
!     CONVERSION FROM THE DERIVATIVE WITH CONSTANT RHOT3 TO CONSTANT RHO
            coeff(1,im3)=log10_opacity_final
         end do
         call ysplin(opacity_table%opal92_grid_logr,coeff,density_end)
         do j=1,density_end
            do i=1,4
               coeff_base_index=i+(j-1)*4
               opacity_table%opal92_surface_spline_coeffs(im2,coeff_base_index)=coeff(i,j)
            end do
         end do
      end do



      if (use_two_z_tables) then
         opacity_table%opal92_surface_x_z2=hydrogen_fraction
         opacity_table%opal92_surface_z_z2=opal_table_z2
         call findex(opacity_table%opal92_grid_x_z2, num_x, hydrogen_fraction, x_index)
         if (x_index.lt.0) x_index=-x_index
         if (x_index.ge.3) x_index=2
         if (x_index.le.0) stop ' ERROR IN X GRID'
         opacity_table%opal92_surface_x_index_z2=x_index
         do im2=1,opacity_table%opal92_num_temps_z2
            x_index=opacity_table%opal92_surface_x_index_z2
            temp6=opacity_table%opal92_grid_logt_z2(im2)
            temp_index=im2
            row_index=im2+(x_index-1)*opacity_table%opal92_num_temps_z2
            density_start=opacity_table%opal92_density_start_index_z2(row_index)
            if (density_start.ne.1) stop ' LL4TH Z1 CHECK NDSS '
            density_end=density_start+opacity_table%opal92_density_count_z2(row_index)-1
            do im3=density_start,density_end
               density_rhot3=opacity_table%opal92_grid_logr_z2(im3)
               call yllo2d2(temp6,density_rhot3,x_index,temp_index,rho_search_index, &
                    opacity0,log10_opacity0,dlnkap_dlnrho0,dlnkap_dlnt0)
               call yllo2d2(temp6,density_rhot3,x_index+1,temp_index,rho_search_index, &
                    opacity1,log10_opacity1,dlnkap_dlnrho1,dlnkap_dlnt1)
               x_fraction_within=(hydrogen_fraction-opacity_table%opal92_grid_x_z2(x_index))/ &
                    (opacity_table%opal92_grid_x_z2(x_index+1)-opacity_table%opal92_grid_x_z2(x_index))
               log10_opacity_final=(log10_opacity1-log10_opacity0)*x_fraction_within+log10_opacity0
               opacity_final=10.0d0**log10_opacity_final
! CONVERSION FROM THE DERIVATIVE WITH CONSTANT RHOT3 TO CONSTANT RHO
               coeff(1,im3)=log10_opacity_final
            end do
            call ysplin(opacity_table%opal92_grid_logr_z2,coeff,density_end)
            do j=1,density_end
               do i=1,4
                  coeff_base_index=i+(j-1)*4
                  opacity_table%opal92_surface_spline_coeffs_z2(im2,coeff_base_index)=coeff(i,j)
               end do
            end do
         end do
      end if



      return
end subroutine ll4th
