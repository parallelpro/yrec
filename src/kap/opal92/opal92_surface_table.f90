!----------------------------------------------------------------------
! opal92_surface_table
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ll4th.f; only variable names, source form, and comment
! style were updated.
!
! THIS IS THE 4TH TABLE OF THE LAWRENCE LIVERMORE OPACITY TABLES.
! IT CONTAINS THE OPACITIES FOR THE SURFACE ABUNDANCES OF X AND Z, AND
! IS USED TO AVOID MANY INTERPOLATIONS TO THE SAME X. ONLY THE
! LOG OPACITY IS KEPT; THE DERIVATIVES RETURNED BY opal92_interp2d
! ARE NOT USED HERE.
!
! 2026 wave 3 (R5): the verbatim second-Z half (on opacity_table%
! opal92(2), Z = opal_table_z2, through the former
! opal92_interp2d_z2) became the internal build_surface_table below,
! called once per set in the original order. Working scalars stay
! host-associated (rho_search_index in particular is the warm-start
! cursor handed to opal92_interp2d and carries over from the first
! set into the second exactly as before). The stdout diagnostics no
! longer say "second-Z".
subroutine opal92_surface_table(hydrogen_fraction, ierr)
      use star_info_lib, only: star

      use opacity_table_lib
      use numerics_lib
      use math_lib
      implicit none
! MHP 10/02 made array dimensions consistent
      integer, parameter :: np = 100

      double precision, intent(in) :: hydrogen_fraction
      integer, intent(out) :: ierr

      double precision :: coeff(4,np)
      integer :: x_index, temp_index, rho_search_index, im2, im3, &
           row_index, coeff_base_index, density_start, density_end, j, i
      double precision :: temp6, density_rhot3, opacity0, log10_opacity0, &
           dlnkap_dlnrho0, dlnkap_dlnt0, opacity1, log10_opacity1, &
           dlnkap_dlnrho1, dlnkap_dlnt1, x_fraction_within, log10_opacity_final

      ierr = 0
      call build_surface_table(opacity_table%opal92(1), star%ctrl%opal_table_z1, ierr)
      if (ierr /= 0) return

      if (star%use_two_z_tables) then
         call build_surface_table(opacity_table%opal92(2), star%ctrl%opal_table_z2, ierr)
      end if

      return

contains

! build_surface_table: build tbl's surface-X spline table at Z = table_z.
subroutine build_surface_table(tbl, table_z, ierr)
      type(opal92_table_set), intent(inout) :: tbl
      double precision, intent(in) :: table_z
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success

!     KEEP THE COMPOSITION OF THE 4TH TABLE.
      tbl%surface_x=hydrogen_fraction
      tbl%surface_z=table_z
      call findex(tbl%grid_x, n_opal92_x, hydrogen_fraction, x_index)
      if (x_index.lt.0) x_index=-x_index
      if (x_index.ge.3) x_index=2
      if (x_index.le.0) then
         write(*,*) 'opal92_surface_table: surface X outside table X grid'
         ierr = 1
         return
      end if
      tbl%surface_x_index=x_index
      do im2=1,tbl%num_temps
         x_index=tbl%surface_x_index
         temp6=tbl%grid_logt(im2)
         temp_index=im2
         row_index=im2+(x_index-1)*tbl%num_temps
         density_start=tbl%density_start_index(row_index)
         if (density_start.ne.1) then
            write(*,*) 'opal92_surface_table: surface-table density grid does not start at index 1'
            ierr = 1
            return
         end if
         density_end=density_start+tbl%density_count(row_index)-1
         do im3=density_start,density_end
            density_rhot3=tbl%grid_logr(im3)
            call opal92_interp2d(tbl,temp6,density_rhot3,x_index,temp_index,rho_search_index, &
                 opacity0,log10_opacity0,dlnkap_dlnrho0,dlnkap_dlnt0, ierr)
            if (ierr /= 0) return
            call opal92_interp2d(tbl,temp6,density_rhot3,x_index+1,temp_index,rho_search_index, &
                 opacity1,log10_opacity1,dlnkap_dlnrho1,dlnkap_dlnt1, ierr)
            if (ierr /= 0) return
            x_fraction_within=(hydrogen_fraction-tbl%grid_x(x_index))/ &
                 (tbl%grid_x(x_index+1)-tbl%grid_x(x_index))
            log10_opacity_final=(log10_opacity1-log10_opacity0)*x_fraction_within+log10_opacity0
            coeff(1,im3)=log10_opacity_final
         end do
         call ysplin(tbl%grid_logr,coeff,density_end)
         do j=1,density_end
            do i=1,4
               coeff_base_index=i+(j-1)*4
               tbl%surface_spline_coeffs(im2,coeff_base_index)=coeff(i,j)
            end do
         end do
      end do
end subroutine build_surface_table

end subroutine opal92_surface_table
