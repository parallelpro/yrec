!----------------------------------------------------------------------
! laol_opacity
!----------------------------------------------------------------------
! The Los Alamos (LAOL89) atomic opacity tables, one per Z (a second at
! another Z when the run blends two), and the LAOL pure-Z table (X = 0,
! Y + Z = 1) that helium-burning cores blend toward.
!     read_laol_tables, build_laol_splines    load and spline the
!                                             LAOL89 tables
!     compute_laol_opacity                    the opacity and its log
!                                             derivatives from them
!     read_pure_z_table, build_pure_z_splines,
!     compute_pure_z_opacity                  the same for the pure-Z
!                                             table
! Both interpolate with cubic splines in log rho and four-point windows
! in T and X, extrapolating up to laol_extrapolation_tolerance; the
! tables live in opacity_tables (opacity_table%laol(1:2), %laol_pure_z).
module laol_opacity
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: kap_settings, laol_opacity_cap, &
           laol_table_set, n_laol_rho, n_laol_t, n_laol_x, opacity_table
      use numerics_lib, only: cspline, locate, splint
      use io_units, only: run_log_unit
      use math_lib
      implicit none
      private
      public :: read_laol_tables, build_laol_splines, compute_laol_opacity, read_pure_z_table, build_pure_z_splines, compute_pure_z_opacity

contains

!----------------------------------------------------------------------
! read_laol_tables (formerly rdlaol)
!----------------------------------------------------------------------
! DBG 4/94 Modified to read in second table for ZRAMP core stuff.
! Reads the LAOL89 opacity table(s) (and, if a second Z table is
! requested, a second LAOL89 table) into opacity_table%laol(1)/(2).
!
! 2026 wave 3 (R5): the two verbatim read sequences became the
! internal read_laol_table below, called once per table; the second
! table's header mixture was read into work_array2(6),(9),(8),(11),
! (1),(3),(2),(5),(10),(4), which is exactly (ix_c, ix_n, ix_o, ix_ne,
! ix_na, ix_mg, ix_al, ix_si, ix_ar, ix_fe), so one read list serves
! both. The run-log message "SECOND OPACITY ARRAY TOO LARGE." is
! kept verbatim for table 2.
subroutine read_laol_tables(laol_work_array, laol_debye_huckel_z, laol_table_path, &
     laol_table2_path, ierr)
      use composition_layout, only: ix_na, ix_al, ix_mg, ix_fe, ix_si, ix_c, &
           ix_o, ix_n, ix_ar, ix_ne

      integer, intent(out) :: ierr
! the 18-element metal mixture from the table header (see the read
! below); the caller hands it to the eos domain (eos_set_debye_huckel_z)
      real(dp), intent(out) :: laol_debye_huckel_z(18)
      real(dp), intent(inout) :: laol_work_array(12)
      character(len=256), intent(in) :: laol_table_path, laol_table2_path

!     THE work_array2/zdh2 TWINS OF laol_work_array/laol_debye_huckel_z
!     are read-list targets for the second table that are never used
!     AFTERWARDS.
      real(dp) :: work_array2(12), zdh2(18)
      integer :: ii, ix, ir, it

      ierr = 0
      opacity_table%laol(2)%table_number = 2
      call read_laol_table(laol_table_path, opacity_table%laol(1), &
           laol_work_array, laol_debye_huckel_z, ierr)
      if (ierr /= 0) return

! DBG 4/94 New stuff follows
      if (kap_settings%use_two_z_tables) then
         call read_laol_table(laol_table2_path, opacity_table%laol(2), &
              work_array2, zdh2, ierr)
         if (ierr /= 0) return
      end if
      return

contains

! read_laol_table: read one LAOL89 table file at `path` into tbl. The
! header's Cox & Stewart metal mixture goes to work_array (the
! caller's laol_work_array for table 1, a discarded local for table
! 2) and its 18-element Debye-Huckel mixture to debye_huckel_z.
subroutine read_laol_table(path, tbl, work_array, debye_huckel_z, ierr)
      character(len=256), intent(in) :: path
      type(laol_table_set), intent(inout) :: tbl
      real(dp), intent(inout) :: work_array(12)
      real(dp), intent(out) :: debye_huckel_z(18)
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success
! runtime-allocated unit for the table (formerly the fixed iolaol =
! 61 for table 1, io_units's fixed iolaol2 = 63 for table 2)
      integer :: laol_unit
!     Zlot/zhit are read-list targets that are never used afterwards.
      real(dp) :: zlot, zhit

      open(newunit=laol_unit,file=path, form='FORMATTED', &
           status='OLD')
!     Read in array sizes
      read(laol_unit,100) tbl%num_x,tbl%num_rho,tbl%num_t
  100 format(/,18x,i2,9x,i3,14x,i3)
      if (tbl%num_x > n_laol_x-1.or.tbl%num_rho > n_laol_rho.or.tbl%num_t > n_laol_t) then
         if (tbl%table_number == 2) then
         write(run_log_unit,*)' SECOND OPACITY ARRAY TOO LARGE.'
         else
         write(run_log_unit,*)' OPACITY ARRAY TOO LARGE.'
         end if
         ! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
         ierr = 1
         return
      end if
!     Read in relative abundances by weigth of the metals c,n,o,ne,na,
!     mg,si,ar,fe (the cox&stewart mix).  Mixture is for the zlot part
!     of the opacity table with LAOL relative abundances scaled
!     TO SUM TO ZLOT.
      read(laol_unit,120) zlot, work_array(ix_c), work_array(ix_n), &
           work_array(ix_o), work_array(ix_ne), work_array(ix_na), &
           work_array(ix_mg), work_array(ix_al), work_array(ix_si), &
           work_array(ix_ar), work_array(ix_fe)
  120 format(47x,f8.5,//,1p6e12.5,/,1p4e12.5)
!     Read zhit the metal abundance for the hit part of table.
! DBG 7/92 need relative abundances of metals for debye-huckel
!     correction. 18 elements, c,n,o,ne,na,mg,al,si,p,
!     S,Cl,Ar,Ca,Ti,Cr,Mn,Fe,Ni scaled to sum to ZHIT
      read(laol_unit,130) zhit, debye_huckel_z
! 130 FORMAT(54X,F8.5,///////)
  130 format(54x,f8.5,/////,1p6e12.5,/,1p6e12.5,/1p6e12.5)
!     Read in h mass fractions of table
      read(laol_unit,140) (tbl%grid_x(ii),ii=1,tbl%num_x)
  140 format(/,(1p6e12.5))
!     Read in density grid of table
      read(laol_unit,150) (tbl%grid_rho(ii),ii=1,tbl%num_rho)
  150 format(/,(1p6e12.5))
!     Read in temperature grid of table
      read(laol_unit,160) (tbl%grid_t(ii),ii=1,tbl%num_t)
  160 format(/,(1p6e12.5))
!     Read in opacities
      read(laol_unit,170)
  170 format(1x)
      do ix=1,tbl%num_x
         do ir=1,tbl%num_rho
            read(laol_unit,200)
  200       format(1x)
            read(laol_unit,210) (tbl%opacity(ix,ir,it),it=1,tbl%num_t)
  210       format(1p6e12.5)
         end do
      end do
      close(laol_unit)
end subroutine read_laol_table

end subroutine read_laol_tables

!----------------------------------------------------------------------
! build_laol_splines (formerly sulaol)
!----------------------------------------------------------------------
! DBG 4/94 Modified to do ZRAMP stuff.
! Builds the natural cubic spline coefficients (in log rho, at each
! tabulated X and T) for the LAOL89 opacity table(s) read by
! read_laol_tables.f90, for use by compute_laol_opacity.f90.
!
! 2026 wave 3 (R5): the spline loop, identical for both tables, is
! the internal spline_laol_table below; the log10 conversion of each
! table's T grid stays here because the second table's loop bound
! was (and is) the first table's num_t.
subroutine build_laol_splines

! MHP 8/25 Removed unused variables
      real(dp) :: row_log10_opacity(n_laol_rho), row_log_rho(n_laol_rho), &
           row_d2opacity(n_laol_rho)
      integer :: it, ix, ir, num_valid_rho

      do it=1, opacity_table%laol(1)%num_t
         opacity_table%laol(1)%grid_t(it) = log10(opacity_table%laol(1)%grid_t(it))
      end do
      call spline_laol_table(opacity_table%laol(1))
! DBG 4/94 Do SPLINE on second opacity table if ZRAMP
      if (kap_settings%use_two_z_tables) then
! NOTE: the log10 conversion of the second table's T grid runs over
! the FIRST table's num_t, as the original did (an R6 decision
! whether that is a latent bug when the two tables differ in T
! extent); it is therefore kept here rather than in the helper.
       do it=1, opacity_table%laol(1)%num_t
          opacity_table%laol(2)%grid_t(it) = log10(opacity_table%laol(2)%grid_t(it))
       end do
         call spline_laol_table(opacity_table%laol(2))
      end if
      return

contains

! spline_laol_table: the natural cubic spline in log rho of every
! (X, T) row of tbl -- token-identical in both former halves.
subroutine spline_laol_table(tbl)
      type(laol_table_set), intent(inout) :: tbl

      do ix=1, tbl%num_x
         do it=1, tbl%num_t
            num_valid_rho=0
            do ir=1, tbl%num_rho
                tbl%slaol_opacity(ix,ir,it) = 0.0d0
                tbl%slaol_log_rho(ix,ir,it) = 0.0d0
                tbl%slaol_d2opacity(ix,ir,it) = 0.0d0
                if (tbl%opacity(ix,ir,it) /= 0.0d0) then
                   num_valid_rho = num_valid_rho+1
                   row_log10_opacity(num_valid_rho) = log10(tbl%opacity(ix,ir,it))
                   row_log_rho(num_valid_rho) = log10(tbl%grid_rho(ir))
                end if
            end do
            if (num_valid_rho >= 4) then
               tbl%slaol_num_points(ix,it)=num_valid_rho
               call cspline(row_log_rho, row_log10_opacity, num_valid_rho, &
                    1.0d30, 1.0d30, row_d2opacity)
               do ir=1,num_valid_rho
                   tbl%slaol_opacity(ix,ir,it) = row_log10_opacity(ir)
                   tbl%slaol_log_rho(ix,ir,it) = row_log_rho(ir)
                   tbl%slaol_d2opacity(ix,ir,it) = row_d2opacity(ir)
               end do
            else
               tbl%slaol_num_points(ix,it) = 0
            end if
         end do
      end do
end subroutine spline_laol_table

end subroutine build_laol_splines

!----------------------------------------------------------------------
! compute_laol_opacity (formerly gtlaol)
!----------------------------------------------------------------------
! Get LAOL89 opacity: cubic-spline interpolate in density, then in
! temperature, then linearly interpolate (or, within tollaol of the
! table edge, linearly extrapolate) in hydrogen abundance.
!
! 2026 wave 3 (R5): tbl is the table to interpolate in --
! opacity_table%laol(1) (the former compute_laol_opacity.f90) or laol(2) (the
! former gtlaol2.f90, deleted: a member-for-member rename of this
! file whose two run-log diagnostics said "OPACITY TABLE #2"; those
! are kept verbatim under `if (tbl%table_number == 2)`).
subroutine compute_laol_opacity(tbl, log10_density, log10_temperature, hydrogen_fraction, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ierr)
      integer :: jerr_gate
      type(laol_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      integer, intent(out) :: ierr

      real(dp) :: row_log10_opacity(n_laol_rho), row_log_rho(n_laol_rho), &
           row_d2opacity(n_laol_rho), dlnkap_dlnrho_by_t(n_laol_t), dlnkap_dlnrho_by_x(4)
      real(dp) :: logt_interp_opacity(n_laol_t), logt_values(n_laol_t), &
           logt_d2opacity(n_laol_t), dlnkap_dlnt_by_x(4)
      real(dp) :: opacity_by_x(4), x_values(4)
      real(dp) :: log_extrap_tolerance
      integer :: t_locate_guess, t_index, x_grid_index, x_grid_index_hi, &
           num_valid_x, x_loop_index, num_valid_t, t_range_lo, t_range_hi, &
           num_valid_rho, rho_loop_index, spline_index_lo, spline_index_hi
      real(dp) :: log10_opacity_value, slope

!     1. Cubic spline interpolate in density
!     2. Cubic spline interpolate in temperature
!     3. Linear interpolate in composition
!     if within tollaol of edge then linear extrapolate
!
!     Tollaol permits some extraploation beyond table edge.
      ierr = 0
      log_extrap_tolerance = log(kap_settings%laol_extrapolation_tolerance)
      call locate(tbl%grid_t, tbl%num_t, log10_temperature, t_locate_guess)
      call locate(tbl%grid_x, tbl%num_x, hydrogen_fraction, x_grid_index)
      if (x_grid_index == tbl%num_x) then
          x_grid_index = tbl%num_x-1
      end if
      if (x_grid_index == 0) then
          x_grid_index = 1
      end if
      num_valid_x = 0
      x_grid_index_hi = x_grid_index + 1
      do x_loop_index=x_grid_index, x_grid_index_hi
         num_valid_t = 0
!        Get range of four tt surrounding T
         call xrng4(t_locate_guess, tbl%num_t, t_range_lo, t_range_hi)
         do t_index=t_range_lo, t_range_hi
            num_valid_rho = tbl%slaol_num_points(x_loop_index,t_index)
            if (num_valid_rho >= 4) then
               do rho_loop_index=1, num_valid_rho
                  row_log10_opacity(rho_loop_index) = &
                       tbl%slaol_opacity(x_loop_index,rho_loop_index,t_index)
                  row_log_rho(rho_loop_index) = &
                       tbl%slaol_log_rho(x_loop_index,rho_loop_index,t_index)
                  row_d2opacity(rho_loop_index) = &
                       tbl%slaol_d2opacity(x_loop_index,rho_loop_index,t_index)
               end do
               if (log10_density > row_log_rho(1) .and. &
                    log10_density < row_log_rho(num_valid_rho)) then
                  call splint(row_log_rho, row_log10_opacity, num_valid_rho, &
                       row_d2opacity, log10_density, log10_opacity_value, &
                       spline_index_lo, spline_index_hi, jerr_gate)
                  ! 2026 numerics-gate opt-in: interpolation failure returns via
                  ! ierr (diagnostic printed at the gate) instead of stopping.
                  if (jerr_gate /= 0) then
                     ierr = jerr_gate
                     return
                  end if
                  num_valid_t = num_valid_t+1
                  logt_interp_opacity(num_valid_t) = log10_opacity_value
                  logt_values(num_valid_t) = tbl%grid_t(t_index)
                  dlnkap_dlnrho_by_t(num_valid_t) = &
                       (row_log10_opacity(spline_index_hi)-row_log10_opacity(spline_index_lo))/ &
                       (row_log_rho(spline_index_hi)-row_log_rho(spline_index_lo))
               else if (log10_density > row_log_rho(1)-log_extrap_tolerance .and. &
                    log10_density <= row_log_rho(1)) then
!                 Put linear extrapolation routines here
                  slope = (row_log10_opacity(2)-row_log10_opacity(1))/ &
                       (row_log_rho(2)-row_log_rho(1))
                  log10_opacity_value = row_log10_opacity(1)+slope*(log10_density-row_log_rho(1))
                  num_valid_t = num_valid_t+1
                  logt_interp_opacity(num_valid_t) = log10_opacity_value
                  logt_values(num_valid_t) = tbl%grid_t(t_index)
                  dlnkap_dlnrho_by_t(num_valid_t) = slope
               else if (log10_density >= row_log_rho(num_valid_rho) .and. &
                    log10_density < row_log_rho(num_valid_rho)+log_extrap_tolerance) then
!                 Put linear extrapolation routines here
                  slope = (row_log10_opacity(num_valid_rho-1)-row_log10_opacity(num_valid_rho))/ &
                       (row_log_rho(num_valid_rho-1)-row_log_rho(num_valid_rho))
                  log10_opacity_value = row_log10_opacity(num_valid_rho)+ &
                       slope*(log10_density-row_log_rho(num_valid_rho))
                  num_valid_t = num_valid_t+1
                  logt_interp_opacity(num_valid_t) = log10_opacity_value
                  logt_values(num_valid_t) = tbl%grid_t(t_index)
                  dlnkap_dlnrho_by_t(num_valid_t) = slope
               end if
            else
               if (tbl%table_number == 2) then
               write(run_log_unit,220) log10_density, log10_temperature
  220          format(' OUTSIDE OPACITY TABLE #2, IN DENSITY.  ', &
                    'LOG(RHO)=',1pe12.3, ' LOG(T)=', 1pe12.3)
               else
               write(run_log_unit,120) log10_density, log10_temperature
  120          format(' OUTSIDE OPACITY TABLE, IN DENSITY.  ', &
                    'LOG(RHO)=',1pe12.3, ' LOG(T)=', 1pe12.3)
               end if
! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib's kap_get).
               ierr = 1
               return
            end if
         end do
         if (num_valid_t >= 4) then
            call cspline(logt_values, logt_interp_opacity, num_valid_t, &
                 1.0d30, 1.0d30, logt_d2opacity)
            call splint(logt_values, logt_interp_opacity, num_valid_t, &
                 logt_d2opacity, log10_temperature, log10_opacity_value, &
                 spline_index_lo, spline_index_hi, jerr_gate)
            if (jerr_gate /= 0) then
               ierr = jerr_gate
               return
            end if
            num_valid_x = num_valid_x + 1
            opacity_by_x(num_valid_x) = log10_opacity_value
            slope=(dlnkap_dlnrho_by_t(spline_index_hi)-dlnkap_dlnrho_by_t(spline_index_lo))/ &
                 (logt_values(spline_index_hi)-logt_values(spline_index_lo))
            dlnkap_dlnrho_by_x(num_valid_x)=dlnkap_dlnrho_by_t(spline_index_lo)+ &
                 slope*(log10_temperature-logt_values(spline_index_lo))
            dlnkap_dlnt_by_x(num_valid_x) = &
                 (logt_interp_opacity(spline_index_hi)-logt_interp_opacity(spline_index_lo))/ &
                 (logt_values(spline_index_hi)-logt_values(spline_index_lo))
            x_values(num_valid_x) = tbl%grid_x(x_loop_index)
         else
            if (tbl%table_number == 2) then
            write(run_log_unit,221) log10_density, log10_temperature
  221       format(' OUTSIDE OPACITY TABLE #2, IN TEMPERATURE.  ', &
                 'LOG(RHO)=',1pe12.3, ' LOG(T)=', 1pe12.3)
            else
            write(run_log_unit,121) log10_density, log10_temperature
  121       format(' OUTSIDE OPACITY TABLE, IN TEMPERATURE.  ', &
                 'LOG(RHO)=',1pe12.3, ' LOG(T)=', 1pe12.3)
            end if
! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib's kap_get).
            ierr = 1
            return
         end if
      end do
!     Both X columns succeeded (each iteration above either adds one or
!     RETURNS WITH ierr), SO num_valid_x IS 2 HERE.
      slope = (opacity_by_x(2)-opacity_by_x(1))/(x_values(2)-x_values(1))
      log10_opacity_value = opacity_by_x(1)+slope*(hydrogen_fraction-x_values(1))
      slope = (dlnkap_dlnt_by_x(2)-dlnkap_dlnt_by_x(1))/(x_values(2)-x_values(1))
      dlnkap_dlnt = dlnkap_dlnt_by_x(1) + slope*(hydrogen_fraction-x_values(1))
      slope = (dlnkap_dlnrho_by_x(2)-dlnkap_dlnrho_by_x(1))/(x_values(2)-x_values(1))
      dlnkap_dlnrho = dlnkap_dlnrho_by_x(1)+slope*(hydrogen_fraction-x_values(1))
      if (log10_opacity_value > 35) then
         opacity = laol_opacity_cap
         log10_opacity = 35.0d0
      else
         opacity = exp10(log10_opacity_value)
         log10_opacity = log10_opacity_value
      end if
end subroutine compute_laol_opacity

!----------------------------------------------------------------------
! read_pure_z_table (formerly rdzlaol)
!----------------------------------------------------------------------
! Reads the pure-Z (Z=1) LAOL89 opacity table.
! MHP 10/02 vector v not used
subroutine read_pure_z_table(pure_z_table_path, ierr)

      integer, intent(out) :: ierr
      character(len=256), intent(in) :: pure_z_table_path
      integer :: n, i, ii, ir, it
! runtime-allocated unit for the pure-Z table (formerly the fixed
! iopurez = 62)
      integer :: pure_z_unit
      real(dp) :: dummy(104)

      ierr = 0
      open(newunit=pure_z_unit, file=pure_z_table_path,form='FORMATTED', &
              status='OLD')
!     Read in array sizes
      read(pure_z_unit,100) n,opacity_table%zlaol_num_rho,opacity_table%zlaol_num_t
  100 format(/,18x,i2,9x,i3,14x,i3)
      if (n /= 1.or.opacity_table%zlaol_num_rho > n_laol_rho.or.opacity_table%zlaol_num_t > n_laol_t) then
         write(run_log_unit,*)' Z OPACITY INPUT ERROR.'
         ! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib facades).
         ierr = 1
         return
      end if
      read(pure_z_unit,120) (dummy(i),i=1,11)
  120 format(47x,f8.5,//,1p6e12.5,/,1p4e12.5)
      read(pure_z_unit,131) dummy(1)
  131 format(54x,f8.5,///////)
!     Read in h mass fractions of table
      read(pure_z_unit,140) dummy(1)
  140 format(/,(1p6e12.5))
!     Read in density grid of table
      read(pure_z_unit,150) (opacity_table%zlaol_logrho_grid(ii),ii=1,opacity_table%zlaol_num_rho)
  150 format(/,(1p6e12.5))
!     Read in temperature grid of table
      read(pure_z_unit,160) (opacity_table%zlaol_logt_grid(ii),ii=1,opacity_table%zlaol_num_t)
  160 format(/,(1p6e12.5))
!     Read in pure Z opacities
      read(pure_z_unit,170)
  170 format(1x)
      do ir=1,opacity_table%zlaol_num_rho
         read(pure_z_unit,200)
  200    format(1x)
         read(pure_z_unit,210) (opacity_table%zlaol_opacity(ir,it),it=1,opacity_table%zlaol_num_t)
  210    format(1p6e12.5)
         end do
      close(pure_z_unit)

end subroutine read_pure_z_table

!----------------------------------------------------------------------
! build_pure_z_splines (formerly zsulaol)
!----------------------------------------------------------------------
! DBG 4/94 Modified to do ZRAMP stuff.
! Builds the natural cubic spline coefficients (in log rho, at each
! tabulated T) for the pure-Z LAOL89 opacity table read by
! read_pure_z_table.f90, for use by compute_pure_z_opacity.f90.
subroutine build_pure_z_splines

      real(dp) :: row_log10_opacity(n_laol_rho), row_log_rho(n_laol_rho), &
           row_d2opacity(n_laol_rho)
      integer :: it, ir, num_valid_rho

      do it=1, opacity_table%zlaol_num_t
         opacity_table%zlaol_logt_grid(it) = log10(opacity_table%zlaol_logt_grid(it))
      end do
        do it=1, opacity_table%zlaol_num_t
            num_valid_rho=0
            do ir=1, opacity_table%zlaol_num_rho
                  opacity_table%zslaol_opacity(ir,it) = 0.0d0
                  opacity_table%zslaol_log_rho(ir,it) = 0.0d0
                  opacity_table%zslaol_d2opacity(ir,it) = 0.0d0
                  if (opacity_table%zlaol_opacity(ir,it) /= 0.0d0) then
                     num_valid_rho = num_valid_rho+1
                     row_log10_opacity(num_valid_rho) = log10(opacity_table%zlaol_opacity(ir,it))
                     row_log_rho(num_valid_rho) = log10(opacity_table%zlaol_logrho_grid(ir))
                  end if
            end do
            if (num_valid_rho >= 4) then
               opacity_table%zslaol_num_points(it)=num_valid_rho
               call cspline(row_log_rho, row_log10_opacity, num_valid_rho, &
                    1.0d30, 1.0d30, row_d2opacity)
               do ir=1,num_valid_rho
                     opacity_table%zslaol_opacity(ir,it) = row_log10_opacity(ir)
                     opacity_table%zslaol_log_rho(ir,it) = row_log_rho(ir)
                     opacity_table%zslaol_d2opacity(ir,it) = row_d2opacity(ir)
               end do
            else
               opacity_table%zslaol_num_points(it) = 0
            end if
        end do
end subroutine build_pure_z_splines

!----------------------------------------------------------------------
! compute_pure_z_opacity (formerly gtpurz)
!----------------------------------------------------------------------
! Get LAOL89 opacity for the pure-Z (Z=1) table: cubic-spline
! interpolate in density, then in temperature (no composition
! interpolation needed since the table is pure Z).
subroutine compute_pure_z_opacity(log10_density, log10_temperature, opacity, &
     log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ierr)
      integer :: jerr_gate
      real(dp), intent(in) :: log10_density, log10_temperature
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      integer, intent(out) :: ierr

      real(dp) :: row_log10_opacity(n_laol_rho), row_log_rho(n_laol_rho), &
           row_d2opacity(n_laol_rho), dlnkap_dlnrho_by_t(n_laol_t)
      real(dp) :: logt_interp_opacity(n_laol_t), logt_values(n_laol_t), &
           logt_d2opacity(n_laol_t)
      real(dp) :: log_extrap_tolerance
      integer :: t_locate_guess, t_index, t_range_lo, t_range_hi, &
           num_valid_t, num_valid_rho, rho_loop_index, spline_index_lo, &
           spline_index_hi
      real(dp) :: log10_opacity_value, slope

!     1. Cubic spline interpolate in density
!     2. Cubic spline interpolate in temperature
!     if within tollaol of edge then linear extrapolate
!
!     Tollaol permits some extraploation beyond table edge.
      ierr = 0
      log_extrap_tolerance = log(kap_settings%laol_extrapolation_tolerance)
      call locate(opacity_table%zlaol_logt_grid, opacity_table%zlaol_num_t, log10_temperature, t_locate_guess)
      num_valid_t = 0
!     Get range of four tt surrounding T
      call xrng4(t_locate_guess, opacity_table%zlaol_num_t, t_range_lo, t_range_hi)
      do t_index=t_range_lo, t_range_hi
         num_valid_rho = opacity_table%zslaol_num_points(t_index)
         if (num_valid_rho >= 4) then
            do rho_loop_index=1, num_valid_rho
               row_log10_opacity(rho_loop_index) = opacity_table%zslaol_opacity(rho_loop_index,t_index)
               row_log_rho(rho_loop_index) = opacity_table%zslaol_log_rho(rho_loop_index,t_index)
               row_d2opacity(rho_loop_index) = opacity_table%zslaol_d2opacity(rho_loop_index,t_index)
            end do
            if (log10_density > row_log_rho(1) .and. &
                 log10_density < row_log_rho(num_valid_rho)) then
               call splint(row_log_rho, row_log10_opacity, num_valid_rho, &
                    row_d2opacity, log10_density, log10_opacity_value, &
                    spline_index_lo, spline_index_hi, jerr_gate)
               ! 2026 numerics-gate opt-in: interpolation failure returns via
               ! ierr (diagnostic printed at the gate) instead of stopping.
               if (jerr_gate /= 0) then
                  ierr = jerr_gate
                  return
               end if
               num_valid_t = num_valid_t+1
               logt_interp_opacity(num_valid_t) = log10_opacity_value
               logt_values(num_valid_t) = opacity_table%zlaol_logt_grid(t_index)
               dlnkap_dlnrho_by_t(num_valid_t) = &
                    (row_log10_opacity(spline_index_hi)-row_log10_opacity(spline_index_lo))/ &
                    (row_log_rho(spline_index_hi)-row_log_rho(spline_index_lo))
            else if (log10_density > row_log_rho(1)-log_extrap_tolerance .and. &
                 log10_density <= row_log_rho(1)) then
!              Put linear extrapolation routines here
               slope = (row_log10_opacity(2)-row_log10_opacity(1))/ &
                    (row_log_rho(2)-row_log_rho(1))
               log10_opacity_value = row_log10_opacity(1)+slope*(log10_density-row_log_rho(1))
               num_valid_t = num_valid_t+1
               logt_interp_opacity(num_valid_t) = log10_opacity_value
               logt_values(num_valid_t) = opacity_table%zlaol_logt_grid(t_index)
               dlnkap_dlnrho_by_t(num_valid_t) = slope
            else if (log10_density >= row_log_rho(num_valid_rho) .and. &
                 log10_density < row_log_rho(num_valid_rho)+log_extrap_tolerance) then
!              Put linear extrapolation routines here
               slope = (row_log10_opacity(num_valid_rho-1)-row_log10_opacity(num_valid_rho))/ &
                    (row_log_rho(num_valid_rho-1)-row_log_rho(num_valid_rho))
               log10_opacity_value = row_log10_opacity(num_valid_rho)+ &
                    slope*(log10_density-row_log_rho(num_valid_rho))
               num_valid_t = num_valid_t+1
               logt_interp_opacity(num_valid_t) = log10_opacity_value
               logt_values(num_valid_t) = opacity_table%zlaol_logt_grid(t_index)
               dlnkap_dlnrho_by_t(num_valid_t) = slope
            end if
         else
            write(run_log_unit,120) log10_density, log10_temperature
  120       format(' OUTSIDE Z OPACITY TABLE, IN DENSITY.  ', &
                 'LOG(RHO)=',1pe12.3, ' LOG(T)=', 1pe12.3)
! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib's kap_get).
            ierr = 1
            return
         end if
      end do
      if (num_valid_t >= 4) then
         call cspline(logt_values, logt_interp_opacity, num_valid_t, &
              1.0d30, 1.0d30, logt_d2opacity)
         call splint(logt_values, logt_interp_opacity, num_valid_t, &
              logt_d2opacity, log10_temperature, log10_opacity_value, &
              spline_index_lo, spline_index_hi, jerr_gate)
         if (jerr_gate /= 0) then
            ierr = jerr_gate
            return
         end if
         slope=(dlnkap_dlnrho_by_t(spline_index_hi)-dlnkap_dlnrho_by_t(spline_index_lo))/ &
              (logt_values(spline_index_hi)-logt_values(spline_index_lo))
         dlnkap_dlnrho = dlnkap_dlnrho_by_t(spline_index_lo)+ &
              slope*(log10_temperature-logt_values(spline_index_lo))
         dlnkap_dlnt = (logt_interp_opacity(spline_index_hi)-logt_interp_opacity(spline_index_lo))/ &
              (logt_values(spline_index_hi)-logt_values(spline_index_lo))
      else
         write(run_log_unit,121) log10_density, log10_temperature
  121    format(' OUTSIDE Z OPACITY TABLE, IN TEMPERATURE.  ', &
              'LOG(RHO)=',1pe12.3, ' LOG(T)=', 1pe12.3)
! 2026 (docs/ROADMAP.md stage 3): stop -> ierr (see kap_lib's kap_get).
         ierr = 1
         return
      end if
      if (log10_opacity_value > 35) then
         opacity = laol_opacity_cap
         log10_opacity = 35.0d0
      else
         opacity = exp10(log10_opacity_value)
         log10_opacity = log10_opacity_value
      end if
end subroutine compute_pure_z_opacity

end module laol_opacity
