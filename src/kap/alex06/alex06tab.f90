!----------------------------------------------------------------------
! alex06tab
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original alex06tab.f; only variable names, source form, and comment
! style were updated.
!
! Generates the fixed-(X,Z) table for the Alexander 2006 low-
! temperature opacities by 4-point Lagrangian interpolation, first in
! Z then in X, from the full set of tables read by readalex06.f90.
! The target X, Z are taken from opacity_table%alex06_cached_x/
! alex06_cached_z, which the caller (getalex06.f90) sets before
! calling this routine.
subroutine alex06tab(ierr)

      use opacity_table_lib
      use numerics_lib
      implicit none
      integer, intent(out) :: ierr

      double precision :: interp_nodes(4), weight_z(4), weight_x(4)
      double precision :: opacity_by_x(4,n_alex06_t,n_alex06_d)
      double precision :: x_max, z_max, interp_target
      integer :: i, j, iz, k, kk, kk2, kk3, kk4

!     XE = DESIRED X; ZE = DESIRED Z
      ierr = 0
      x_max = 1.0d0 - opacity_table%alex06_cached_z
!     CHECK THAT THE REQUESTED COMPOSITION IS INSIDE TABLE BOUNDS
      if (opacity_table%alex06_cached_x.lt.0.0d0 .or. opacity_table%alex06_cached_x .gt. x_max) then
         write(*,5) opacity_table%alex06_cached_x, opacity_table%alex06_cached_z
    5    format('ILLEGAL COMPOSITION (X,Z) = ',2f6.2,' IN ALEX06.RUN STOPPED')
! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib's kap_get).
         ierr = 1
         return
      endif
!     PERMIT EXTRAPOLATION IN Z BY UP TO 1 TABLE ELEMENT
      z_max = opacity_table%alex06_grid_z(n_alex06_z)+(opacity_table%alex06_grid_z(n_alex06_z)-opacity_table%alex06_grid_z(n_alex06_z-1))
      if (opacity_table%alex06_cached_z.lt.0.0d0 .or. opacity_table%alex06_cached_z .gt. z_max) then
         write(*,5) opacity_table%alex06_cached_x, opacity_table%alex06_cached_z
! 2026 (ROADMAP.md stage 3): stop -> ierr (see kap_lib's kap_get).
         ierr = 1
         return
      endif
!     FIND 4 NEAREST TABLES IN Z.
      do i = 3,n_alex06_z-2
         if (opacity_table%alex06_cached_z.le.opacity_table%alex06_grid_z(i)) then
            iz = i - 2
            exit
         endif
      end do
      if (i > (n_alex06_z-2)) then
      iz = n_alex06_z - 3
      end if
!     FIND 4 NEAREST TABLES IN X.
      do i = 3,n_alex06_x-2
         if (opacity_table%alex06_cached_x.le.opacity_table%alex06_grid_x(i)) then
            opacity_table%alex06_index_x = i - 2
            exit
         endif
      end do
      if (i > (n_alex06_x-2)) then
!     NO TABLE FOR X > 0.9 IF Z =0.10 OR MORE
      if (opacity_table%alex06_cached_z.ge.0.1d0) then
         opacity_table%alex06_index_x = n_alex06_x - 4
      else
         opacity_table%alex06_index_x = n_alex06_x - 3
      endif
      end if
!     INTERPOLATION FACTORS FOR Z
      do i = 1,4
         interp_nodes(i) = opacity_table%alex06_grid_z(iz+i-1)
      end do
      interp_target = opacity_table%alex06_cached_z
      call intrp2(interp_nodes, weight_z, interp_target)
!     THE DIFFERENCE IN THE NUMBER OF TABLES FOR THE Z=0.1 CASE REQUIRES SOME
!     CARE IN X INTERPOLATION.  FIRST 3 X CASES CAN BE TREATED NORMALLY.
!     (2026 R3: the four Z tables kk, kk+1, kk+2, kk+3 are consecutive
!     here, so the sum is lagrange4 over the contiguous slice.)
      do k = 1,3
         kk = n_alex06_z*(opacity_table%alex06_index_x+k-2)+iz
         do i = 1,n_alex06_t
            do j = 1,n_alex06_d
               opacity_by_x(k,i,j) = lagrange4(weight_z, opacity_table%alex06_full_opacity(kk:kk+3,i,j))
            end do
         end do
      end do
!     IF IN THE HIGH Z AND HIGH X DOMAIN THE TOP TABLE IS X = 1-Z (ENTRY NUMX) EXCEPT
!     FOR THE Z=0.10 CASE (WHERE THE X=0.9 CASE DOUBLES AS THE X=1-Z CASE).
      if (opacity_table%alex06_index_x.eq.n_alex06_x-3.and.iz.eq.n_alex06_z-3) then
!        USE DIFFERENT INDEXING FOR THE LAST TABLE
         kk = n_alex06_z*(opacity_table%alex06_index_x+2)+iz
         kk2 = kk + 1
         kk3 = kk2 + 1
         kk4 = n_alex06_z*(opacity_table%alex06_index_x+2)
      else
         kk = n_alex06_z*(opacity_table%alex06_index_x+2)+iz
         kk2 = kk + 1
         kk3 = kk2 + 1
         kk4 = kk3 + 1
      endif
!     (Inline rather than lagrange4: kk4 is not kk+3 in the Z=0.10 case.)
      do i = 1,n_alex06_t
         do j = 1,n_alex06_d
            opacity_by_x(4,i,j) = weight_z(1)*opacity_table%alex06_full_opacity(kk,i,j)+ &
                 weight_z(2)*opacity_table%alex06_full_opacity(kk2,i,j) + &
                 weight_z(3)*opacity_table%alex06_full_opacity(kk3,i,j) + &
                 weight_z(4)*opacity_table%alex06_full_opacity(kk4,i,j)
         end do
      end do
!     NOW DO X INTERPOLATION
!     INTERPOLATION FACTORS FOR X
      do i = 1,3
         interp_nodes(i) = opacity_table%alex06_grid_x(opacity_table%alex06_index_x+i-1)
      end do
      if (iz.eq.n_alex06_z-3) then
         interp_nodes(4) = 1.0d0-opacity_table%alex06_cached_z
      else
         interp_nodes(4) = opacity_table%alex06_grid_x(opacity_table%alex06_index_x+3)
      endif
      interp_target = opacity_table%alex06_cached_x
      call intrp2(interp_nodes, weight_x, interp_target)
      do i = 1,n_alex06_t
         do j = 1,n_alex06_d
            opacity_table%alex06_opacity(i,j) = lagrange4(weight_x, opacity_by_x(1:4,i,j))
         end do
      end do
      return
end subroutine alex06tab
