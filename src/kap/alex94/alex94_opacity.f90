!----------------------------------------------------------------------
! alex94_opacity
!----------------------------------------------------------------------
! The low-temperature molecular opacities of Alexander & Ferguson
! (1994/95): seven files, one per X, each with tables in Z, log T and
! log R.
!     read_alex94_tables          read the seven files, then build the
!                                 fixed-Z tables at the run's Z
!     build_alex94_surface_table  the fixed-X table at the envelope X
!     compute_alex94_opacity      the opacity and its log derivatives
! The tables live in opacity_tables (opacity_table%alex94_*).
module alex94_opacity
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables
      use numerics_lib, only: interp, intrp2, lagrange4
      use physical_constants, only: ln10
      use io_units, only: run_log_unit
      use math_lib
      implicit none
      private
      public :: read_alex94_tables, build_alex94_surface_table, compute_alex94_opacity

contains

!----------------------------------------------------------------------
! read_alex94_tables
!----------------------------------------------------------------------
! alxtbl
! ALEXANDER OPACITIES   APJ,437,879,1994
! IN EACH 7 FILE FOR X = 0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.8
! DESCENDING ORDER OF Z = 0.0, 0.00001, 0.00003, 0.0001, 0.0003,
!     (15)               0.001, 0.002, 0.004, 0.01, 0.02, 0.03,
!                        0.04, 0.06, 0.08, 0.10
! DESCENDING ORDER OF LOGT  4.10 4.05 4.00 3.95 3.90 3.85 3.80 3.75
!     (23)                 3.70 3.65 3.60 3.55 3.50 3.45 3.40 3.35
!                          3.30 3.25 3.20 3.15 3.10 3.05 3.00
! ASCENDING ORDER OF  LOGR -7.0 -6.5 -6.0 -5.5 -5.0 -4.5 -4.0 -3.5
!     (17)                -3.0 -2.5 -2.0 -1.5 -1.0 -0.5  0.0  0.5  1.0
! WHERE R = RHO/T_6**3
! READ IN ALL Z VALUES AT FIXED X AND GENERATE OPACITY ENTRIES AT
! EITHER ONE OR TWO FIXED VALUES OF Z.  THEN GENERATE SURFACE X TABLE
! IF MODEL HAS NO DIFFUSION.  INTERPOLATION BETWEEN COMPOSITIONS IS
! 4-POINT LAGRANGIAN.
!
! MHP 8/25 Added file name to subroutine call
subroutine read_alex94_tables(alex95_table_paths, ierr)
      integer, intent(out) :: ierr

! the table file's unit
      integer :: table_unit
      character(len=256), intent(in) :: alex95_table_paths(7)

! --- local arrays ---
      real(dp) :: row_opacity_temp(n_alex94_d)
! alex94_grid_x/alex94_grid_z/alex94_grid_logt/alex94_grid_logr/
! alex94_cached_x/alex94_cached_z/alex94_index_x/alex94_index_t/
! alex94_index_r defaults moved to opacity_tables.f90: DATA can no
! longer target them here now that they're use-associated.
      integer :: table_index, i, ii, j, k, row_density_count
      real(dp) :: header_x, header_z, row_temp, row_logr0
      real(dp) :: target_z

      do table_index=1,n_alex94_x
!        OPEN EACH OF THE TABLES AT FIXED X WITH A RANGE OF Z.
         open(newunit=table_unit,file=alex95_table_paths(table_index), &
              status='OLD', form='FORMATTED')
!        READ IN INITIAL X AND Z; ENSURE THAT THEY HAVE THE EXPECTED VALUES.
         do i = n_alex94_z,1,-1
!           INDEX FOR STORING TABLES : 15 SETS OF METAL ABUNDANCE WITH 7 SETS
!           OF HYDROGEN FOR EACH.  THESE ARE STORED IN THE ORDER (X,Z) OF
!           (0,0),(0.1,0),...(0.8,0),(0,0.00001),...
            ii = table_index + 7*(i-1)
!           HEADER INFORMATION: X AND Z
            read(table_unit,10) header_x, header_z
   10       format(18x,f6.2,2x,f7.2)
            if (header_x.ne.opacity_table%alex94_grid_x(table_index) .or. &
                 header_z.ne.opacity_table%alex94_grid_z(i)) then
               write(run_log_unit,20) header_x, opacity_table%alex94_grid_x(table_index), &
                    header_z, opacity_table%alex94_grid_z(i)
   20          format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                    1x,'EXPECTED AND ACTUAL X,Z',4f7.2,' RUN STOPPED')
               ! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
               ierr = 1
               return
            endif
!           OPACITY INFORMATION AT EACH SHELL: CHECK FOR CONSISTENCY WITH T,
!           STARTING R.  STORE IN A NUMZ*NUMT*NUMR ARRAY.
            do j = n_alex94_t,1,-1
               read(table_unit,30) opacity_table%alex94_index_t, row_density_count, &
                    row_temp, row_logr0, (row_opacity_temp(k),k=1,n_alex94_d)
               do k = 1, n_alex94_d
                  opacity_table%alex94_full_opacity(ii,j,k) = row_opacity_temp(k)
               end do
   30          format(i2,i3,f6.3,f5.1,8f8.3/9f8.3)
               if (row_density_count.ne.n_alex94_d .or. row_temp.ne.opacity_table%alex94_grid_logt(j) &
                    .or. row_logr0.ne.opacity_table%alex94_grid_logr(1)) then
                  write(run_log_unit,40) row_density_count, row_temp, &
                       opacity_table%alex94_grid_logt(j), row_logr0, opacity_table%alex94_grid_logr(1)
   40             format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                       1x,'EXPECTED AND ACTUAL T,RHO',i3,4f7.2, &
                       ' RUN STOPPED')
                  ! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
                  ierr = 1
                  return
               endif
            end do
         end do
         close(unit=table_unit)
      end do
      ierr = 0
      target_z = kap_settings%alex_table_z1
      call build_alex94_fixed_z_table(target_z)
      return
end subroutine read_alex94_tables

!----------------------------------------------------------------------
! build_alex94_fixed_z_table (formerly alex94_fixed_z_table)
!----------------------------------------------------------------------
! alxztab
! Generates the fixed-Z table (one table per X, interpolated in Z)
! from the full set of Alexander 1994 low-temperature opacity tables,
! by 4-point Lagrangian interpolation in Z.
subroutine build_alex94_fixed_z_table(metal_fraction)


      real(dp), intent(in) :: metal_fraction

      real(dp) :: interp_z(4), weight_z(4)
      integer :: i, jj, kk, idz, ii1, ii2, ii3, ii4

!     LOCATE FOUR NEAREST TABLES IN Z
      do i=3,n_alex94_z-1
         if (metal_fraction.lt.opacity_table%alex94_grid_z(i)) then
            idz = i - 2
            exit
         endif
      end do
      if (i > (n_alex94_z-1)) then
      idz = n_alex94_z - 3
      end if
      do i = 1,4
         interp_z(i) = opacity_table%alex94_grid_z(idz+i-1)
      end do
!     GET INTERPOLATION FACTORS FOR Z.
      call intrp2(interp_z, weight_z, metal_fraction)
!     INTERPOLATE IN Z AT FIXED X.
      do i = 1,n_alex94_x
!        INDICES FOR TABLES: SETS OF 7 X AT FIXED Z, STARTING FROM
!        Z= 0 AND GOING TO Z = 0.1.
         ii1 = i + n_alex94_x*(idz-1)
         ii2 = ii1 + n_alex94_x
         ii3 = ii2 + n_alex94_x
         ii4 = ii3 + n_alex94_x
         do jj = 1,n_alex94_t
            do kk = 1,n_alex94_d
               opacity_table%alex94_opacity(i,jj,kk) = weight_z(1)*opacity_table%alex94_full_opacity(ii1,jj,kk)+ &
                    weight_z(2)*opacity_table%alex94_full_opacity(ii2,jj,kk) + &
                    weight_z(3)*opacity_table%alex94_full_opacity(ii3,jj,kk) + &
                    weight_z(4)*opacity_table%alex94_full_opacity(ii4,jj,kk)
            end do
         end do
      end do
      opacity_table%alex94_cached_z = metal_fraction
      return
end subroutine build_alex94_fixed_z_table

!----------------------------------------------------------------------
! build_alex94_surface_table (formerly alex94_surface_table)
!----------------------------------------------------------------------
! alx8th
! Generates the surface-X (envelope-abundance) table for the
! Alexander 1994 low-temperature opacity tables by 4-point Lagrangian
! interpolation in X, storing it in table slot 8 of
! opacity_table%alex94_full_opacity.
subroutine build_alex94_surface_table(hydrogen_fraction)


      real(dp), intent(in) :: hydrogen_fraction

      real(dp) :: interp_x(4), weight_x(4)
      integer :: i, j

      opacity_table%alex94_cached_x = hydrogen_fraction
!     FIND 4 NEAREST TABLES IN X.
      opacity_table%alex94_index_x = alex_x_stencil_start(hydrogen_fraction)
      do i = 1,4
         interp_x(i) = opacity_table%alex94_grid_x(opacity_table%alex94_index_x+i-1)
      end do
      call intrp2(interp_x, weight_x, hydrogen_fraction)
      do i = 1,n_alex94_t
         do j = 1,n_alex94_d
            opacity_table%alex94_opacity(8,i,j) = weight_x(1)*opacity_table%alex94_opacity(opacity_table%alex94_index_x,i,j)+ &
                 weight_x(2)*opacity_table%alex94_opacity(opacity_table%alex94_index_x+1,i,j) + &
                 weight_x(3)*opacity_table%alex94_opacity(opacity_table%alex94_index_x+2,i,j) + &
                 weight_x(4)*opacity_table%alex94_opacity(opacity_table%alex94_index_x+3,i,j)
         end do
      end do
      return
end subroutine build_alex94_surface_table

!----------------------------------------------------------------------
! compute_alex94_opacity (formerly alex94_interp3d)
!----------------------------------------------------------------------
! yalo3d
! MHP 4/98 Computes Alexander 1994 low-temperature opacities by
! 4-point Lagrangian interpolation in log T and log R = rho/T6**3,
! and (away from the envelope X, Z) also in X.
subroutine compute_alex94_opacity(log10_density, log10_temperature, hydrogen_fraction, &
     metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt)


      real(dp), intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

      real(dp) :: interp_nodes(4), weight_t(4), dweight_t(4), &
           weight_r(4), dweight_r(4), weight_x(4), opacity_row(4), &
           dlnkap_dlnr_row(4), opacity_x(4), dlnkap_dlnt_x(4), &
           dlnkap_dlnr_x(4)
      real(dp) :: delta_z, logr, saved_r
      logical :: extrapolate_linear
      integer :: i, j, jj, ii

      delta_z = abs(metal_fraction-opacity_table%alex94_cached_z)
!     ENSURE THAT OPACITY TABLE HAS THE SAME Z VALUE AS THE ENVELOPE.
      if (delta_z.gt.alex_composition_tol) then
         call build_alex94_fixed_z_table(metal_fraction)
         call build_alex94_surface_table(hydrogen_fraction)
      endif
!     COMPUTE R FOR GRID POINT
      logr = log10_density - 3.0d0*(log10_temperature-6.0d0)
!     FIND NEAREST GRID POINTS IN T.
!     Insure that index IT is within the required array bounds  llp  8/19/08
      if (opacity_table%alex94_index_t .lt. 1) opacity_table%alex94_index_t=1
      if ((opacity_table%alex94_index_t+2) .gt. n_alex94_t) opacity_table%alex94_index_t=n_alex94_t-2
      call stencil4_locate(opacity_table%alex94_grid_logt, n_alex94_t, log10_temperature, opacity_table%alex94_index_t)
!     INTERPOLATION FACTORS IN LOG T
      do i = 1,4
         interp_nodes(i) = opacity_table%alex94_grid_logt(opacity_table%alex94_index_t+i-1)
      end do
      call interp(interp_nodes, weight_t, dweight_t, log10_temperature)
!     FIND NEAREST GRID POINTS IN R = RHO/T6**3
      call stencil4_locate(opacity_table%alex94_grid_logr, n_alex94_d, logr, opacity_table%alex94_index_r)
!     INTERPOLATION FACTORS IN LOG R
      if (logr.gt.opacity_table%alex94_grid_logr(n_alex94_d).and. &
           abs(hydrogen_fraction-opacity_table%alex94_cached_x).lt.alex_composition_tol) then
         extrapolate_linear = .true.
         saved_r = logr
         logr = opacity_table%alex94_grid_logr(n_alex94_d)
      else
         extrapolate_linear = .false.
      endif
      do i = 1,4
         interp_nodes(i) = opacity_table%alex94_grid_logr(opacity_table%alex94_index_r+i-1)
      end do
      call interp(interp_nodes, weight_r, dweight_r, logr)
!     NOW EITHER INTERPOLATE IN SURFACE X TABLE OR CALCULATE OPACITY AT
!     4 DIFFERENT VALUES OF X AND INTERPOLATE IN X.
      if (abs(hydrogen_fraction-opacity_table%alex94_cached_x).lt.alex_composition_tol) then
!        SURFACE ABUNDANCE TABLE
!        INTERPOLATE IN LOG R AT FIXED T
         do i = 1,4
            ii = opacity_table%alex94_index_t+i - 1
            opacity_row(i) = weight_r(1)*opacity_table%alex94_opacity(8,ii,opacity_table%alex94_index_r) + &
                 weight_r(2)*opacity_table%alex94_opacity(8,ii,opacity_table%alex94_index_r+1) + &
                 weight_r(3)*opacity_table%alex94_opacity(8,ii,opacity_table%alex94_index_r+2) + &
                 weight_r(4)*opacity_table%alex94_opacity(8,ii,opacity_table%alex94_index_r+3)
            dlnkap_dlnr_row(i) = dweight_r(1)*opacity_table%alex94_opacity(8,ii,opacity_table%alex94_index_r) + &
                 dweight_r(2)*opacity_table%alex94_opacity(8,ii,opacity_table%alex94_index_r+1) + &
                 dweight_r(3)*opacity_table%alex94_opacity(8,ii,opacity_table%alex94_index_r+2) + &
                 dweight_r(4)*opacity_table%alex94_opacity(8,ii,opacity_table%alex94_index_r+3)
         end do
         if (extrapolate_linear) then
            do i = 1,4
               ii = opacity_table%alex94_index_t+i-1
               opacity_row(i) = opacity_row(i)+(saved_r-logr)* &
                    (opacity_table%alex94_opacity(8,ii,n_alex94_d)-opacity_table%alex94_opacity(8,ii,n_alex94_d-1))/ &
                    (opacity_table%alex94_grid_logr(n_alex94_d)-opacity_table%alex94_grid_logr(n_alex94_d-1))
            end do
            logr = saved_r
         endif
!        INTERPOLATE IN T
         log10_opacity = lagrange4(weight_t, opacity_row)
!        D LN CAPPA/D LN T AT FIXED R
         dlnkap_dlnt = lagrange4(dweight_t, opacity_row)
!        INTERPOLATE IN D LN CAPPA/ D LN R AT FIXED T
         dlnkap_dlnrho = lagrange4(weight_t, dlnkap_dlnr_row)
!        CORRECT FROM DERIVATIVE AT FIXED R TO DERIVATIVE AT FIXED RHO
         dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
         opacity = exp(ln10*log10_opacity)
      else
!        FIND 4 NEAREST TABLES IN X.
         opacity_table%alex94_index_x = alex_x_stencil_start(hydrogen_fraction)
         do i = 1,4
            interp_nodes(i) = opacity_table%alex94_grid_x(opacity_table%alex94_index_x+i-1)
         end do
         call intrp2(interp_nodes, weight_x, hydrogen_fraction)
!        INTERPOLATE IN LOG R AT FIXED T
         do j = 1,4
            jj = opacity_table%alex94_index_x+j-1
            do i = 1,4
               ii = opacity_table%alex94_index_t+i - 1
               opacity_row(i) = weight_r(1)*opacity_table%alex94_opacity(jj,ii,opacity_table%alex94_index_r) + &
                    weight_r(2)*opacity_table%alex94_opacity(jj,ii,opacity_table%alex94_index_r+1) + &
                    weight_r(3)*opacity_table%alex94_opacity(jj,ii,opacity_table%alex94_index_r+2) + &
                    weight_r(4)*opacity_table%alex94_opacity(jj,ii,opacity_table%alex94_index_r+3)
               dlnkap_dlnr_row(i) = dweight_r(1)*opacity_table%alex94_opacity(jj,ii,opacity_table%alex94_index_r) + &
                    dweight_r(2)*opacity_table%alex94_opacity(jj,ii,opacity_table%alex94_index_r+1) + &
                    dweight_r(3)*opacity_table%alex94_opacity(jj,ii,opacity_table%alex94_index_r+2) + &
                    dweight_r(4)*opacity_table%alex94_opacity(jj,ii,opacity_table%alex94_index_r+3)
            end do
!           INTERPOLATE IN T
            opacity_x(j) = lagrange4(weight_t, opacity_row)
!           D LN CAPPA/D LN T AT FIXED R
            dlnkap_dlnt_x(j) = lagrange4(dweight_t, opacity_row)
!           INTERPOLATE IN D LN CAPPA/ D LN R AT FIXED T
            dlnkap_dlnr_x(j) = lagrange4(weight_t, dlnkap_dlnr_row)
         end do
!        INTERPOLATE IN X
         log10_opacity = lagrange4(weight_x, opacity_x)
         dlnkap_dlnt = lagrange4(weight_x, dlnkap_dlnt_x)
         dlnkap_dlnrho = lagrange4(weight_x, dlnkap_dlnr_x)
!        CORRECT FROM DERIVATIVE AT FIXED R TO DERIVATIVE AT FIXED RHO
         dlnkap_dlnt = dlnkap_dlnt - 3.0d0*dlnkap_dlnrho
         opacity = exp(ln10*log10_opacity)
      endif
      return
end subroutine compute_alex94_opacity

end module alex94_opacity
