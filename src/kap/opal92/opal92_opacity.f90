!----------------------------------------------------------------------
! opal92_opacity
!----------------------------------------------------------------------
! The atomic opacities of the 1992 OPAL (Livermore) tables, one file per
! Z with tables in X, log T and log R; a second file at another Z when
! the run blends two Z tables.
!     read_opal92_tables          read the table file(s) and spline them
!     build_opal92_surface_table  the fixed-X table at the envelope
!                                 composition, per Z table
!     compute_opal92_opacity      the opacity and its log derivatives by
!                                 interpolation in X, T and R
! The tables live in opacity_tables (opacity_table%opal92(1:2)).
module opal92_opacity
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: kap_settings, n_opal92_d, n_opal92_t, &
           n_opal92_x, opacity_table, opal92_missing_opacity, &
           opal92_table_set, opal92_x_match_tol
      use numerics_lib, only: findex, intpol, ysplin
      use math_lib
      implicit none
      private
      public :: read_opal92_tables, build_opal92_surface_table, compute_opal92_opacity

contains

!----------------------------------------------------------------------
! read_opal92_tables
!----------------------------------------------------------------------
! DBG 5/94 Modified to read in second opacity table at different Z.
! Reads the OPAL92 (Lawrence Livermore) opacity table(s) and, via
! prepare_opal92_tables, builds the spline interpolation coefficients used
! later by build_opal92_surface_table.f90 and the OPAL92 lookup routine
! compute_opal92_opacity.
!
! 2026 wave 3 (R5): the verbatim second-Z read loop (into
! opacity_table%opal92(2)) became the internal read_opal92_table
! below, called once per set in the original order. The two OPEN
! statements stay here as they were (the first table on the fixed
! star%ctrl%opal92_table_unit, the second on a NEWUNIT); the READ,
! FORMAT and CLOSE statements moved into the helper verbatim, and the
! loop's working scalars stay host-associated.
subroutine read_opal92_tables(opal92_table_path, opal92_table2_path, ierr)
! runtime-allocated unit for the second (Z-ramp) OPAL92 table
! (formerly io_units's fixed ioopal2 = 64)
      integer :: opal92b_unit

      character(len=256), intent(in) :: opal92_table_path, opal92_table2_path
      integer, intent(out) :: ierr
! the table file's unit
      integer :: opal92a_unit
!     local_grid_z IS A READ-LIST TARGET THAT IS NEVER USED AFTERWARDS.
      real(dp) :: local_grid_z(n_opal92_x)
      integer :: i, k, density_index, num_temps_read
      real(dp) :: grid_temp_k

!     OPEN TABLE
      ierr = 0
      open(newunit=opal92a_unit,file=opal92_table_path)
      call read_opal92_table(opal92a_unit, opacity_table%opal92(1), ierr)
      if (ierr /= 0) return

! DBG 5/94 Second Opacity Table read here
      if (kap_settings%use_two_z_tables) then
         open(newunit=opal92b_unit,file=opal92_table2_path)
         call read_opal92_table(opal92b_unit, opacity_table%opal92(2), ierr)
         if (ierr /= 0) return
      end if
!
      call prepare_opal92_tables(ierr)

      return

contains

! read_opal92_table: read the already-opened table on `unit` into tbl
! and close the unit.
subroutine read_opal92_table(unit, tbl, ierr)
      integer, intent(in) :: unit
      type(opal92_table_set), intent(inout) :: tbl
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success

      do i=1,n_opal92_x
!        READ GRID POINT FOR ABUNDANCE
!        READ NUMBER OF GRIDS FOR DENSITY, AND TEMPERATURE
        read(unit,190,end=97) tbl%grid_x(i), local_grid_z(i)
  190   format(33x,f7.4,2x,f7.4)
         read(unit,'()')
!        READ  LOG(DENSITY/TEMPERATURE**3)
            read(unit, 200) (tbl%grid_logr(density_index), density_index=1, n_opal92_d)
  200   format (6x, 17f7.1)
!        READ GRID VALUES FOR TEMPERATURE, AND OPACITY TABLE
         do k=1, n_opal92_t
         read(unit,196,end=93) grid_temp_k, &
              (tbl%log10_opacity(k+(i-1)*n_opal92_t,density_index),density_index=1,n_opal92_d)
         tbl%grid_logt(k)=log10(grid_temp_k)
         end do
   93    num_temps_read=k-1
  196    format(18f7.3)
!
      end do
!     CLOSE THE TABLE WE HAVE READ
   97 close(unit,err=99)
      tbl%num_temps = num_temps_read
      tbl%num_x=i-1
      return
   99 continue
      write(*,*) 'read_opal92_tables: error closing table file'
      ierr = 1
end subroutine read_opal92_table

end subroutine read_opal92_tables

!----------------------------------------------------------------------
! prepare_opal92_tables (formerly opal92_table_prep)
!----------------------------------------------------------------------
! DBG 5/94 Modified to include ZRAMP and ZDIFF stuff (second opacity
! table). Builds the cubic-spline coefficients (in density) for the
! OPAL92 opacity tables, and for the second (different-Z) OPAL92
! table set when use_two_z_tables is set. Note the two strides for
! the same logical row: index1 uses the number of temperature rows
! actually read, while log10_opacity is addressed with the fixed
! n_opal92_t stride that read_opal92_tables stores with.
!
! 2026 wave 3 (R5): the verbatim second-Z loop nest (on
! opacity_table%opal92(2)) became the internal prep_table below,
! called once per set in the original order; the loop's working
! scalars stay host-associated. The stdout diagnostic no longer says
! "second".
subroutine prepare_opal92_tables(ierr)
      integer, intent(out) :: ierr
! MHP 10/02 made array dimensions consistent
      integer, parameter :: np = 100

      real(dp) :: spline_work(4,np), density_nodes(n_opal92_d)
      integer :: ix, it, index1, jd, id, index2, j, i
      real(dp) :: chkd, chko

      ierr = 0

      call prep_table(opacity_table%opal92(1), ierr)
      if (ierr /= 0) return
!
! DBG 5/94 ZRAMP stuff
      if (kap_settings%use_two_z_tables) then
         call prep_table(opacity_table%opal92(2), ierr)
      end if

      return

contains

! prep_table: density-spline coefficients for every (X,T) row of tbl.
subroutine prep_table(tbl, ierr)
      type(opal92_table_set), intent(inout) :: tbl
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success

      do ix = 1,tbl%num_x
      do it = 1,tbl%num_temps
       index1 = it + (ix-1)*tbl%num_temps
       jd = 0
        do id = 1,n_opal92_d
          chkd = tbl%grid_logr(id)
          chko = tbl%log10_opacity(it+n_opal92_t*(ix-1),id)
!>>>> CHECK THE EMPTY REGION
          if (chko <= opal92_missing_opacity) cycle
          if (jd <= 0) then
             tbl%density_start_index(index1) = id
             if (id /= 1) then
                write(*,*) 'prepare_opal92_tables: opal92 table density grid does not start at index 1'
                ierr = 1
                return
             end if
          endif
          jd = jd + 1
          density_nodes(jd) = chkd
!>>>> CHECK THE OPACITY VALUE IN THE TABLE
          spline_work(1,jd) = chko
        end do
        tbl%density_count(index1) = jd
        if (jd <= 1) cycle
        call ysplin(density_nodes, spline_work, jd)
        do j = 1,jd
        do i = 1,4
         index2 = i + (j-1)*4
         tbl%spline_coeffs(index1,index2) = spline_work(i,j)
        end do
        end do
      end do
      end do
end subroutine prep_table

end subroutine prepare_opal92_tables

!----------------------------------------------------------------------
! build_opal92_surface_table (formerly opal92_surface_table)
!----------------------------------------------------------------------
! THIS IS THE 4TH TABLE OF THE LAWRENCE LIVERMORE OPACITY TABLES.
! IT CONTAINS THE OPACITIES FOR THE SURFACE ABUNDANCES OF X AND Z, AND
! IS USED TO AVOID MANY INTERPOLATIONS TO THE SAME X. ONLY THE
! LOG OPACITY IS KEPT; THE DERIVATIVES RETURNED BY interpolate_opal92_2d
! ARE NOT USED HERE.
!
! 2026 wave 3 (R5): the verbatim second-Z half (on opacity_table%
! opal92(2), Z = opal_table_z2, through the former
! opal92_interp2d_z2) became the internal build_surface_table below,
! called once per set in the original order. Working scalars stay
! host-associated (rho_search_index in particular is the warm-start
! cursor handed to interpolate_opal92_2d and carries over from the first
! set into the second exactly as before). The stdout diagnostics no
! longer say "second-Z".
subroutine build_opal92_surface_table(hydrogen_fraction, ierr)
! MHP 10/02 made array dimensions consistent
      integer, parameter :: np = 100

      real(dp), intent(in) :: hydrogen_fraction
      integer, intent(out) :: ierr

      real(dp) :: coeff(4,np)
      integer :: x_index, temp_index, rho_search_index, im2, im3, &
           row_index, coeff_base_index, density_start, density_end, j, i
      real(dp) :: temp6, density_rhot3, opacity0, log10_opacity0, &
           dlnkap_dlnrho0, dlnkap_dlnt0, opacity1, log10_opacity1, &
           dlnkap_dlnrho1, dlnkap_dlnt1, x_fraction_within, log10_opacity_final

      ierr = 0
      call build_surface_table(opacity_table%opal92(1), kap_settings%opal_table_z1, ierr)
      if (ierr /= 0) return

      if (kap_settings%use_two_z_tables) then
         call build_surface_table(opacity_table%opal92(2), kap_settings%opal_table_z2, ierr)
      end if

      return

contains

! build_surface_table: build tbl's surface-X spline table at Z = table_z.
subroutine build_surface_table(tbl, table_z, ierr)
      type(opal92_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: table_z
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success

!     KEEP THE COMPOSITION OF THE 4TH TABLE.
      tbl%surface_x=hydrogen_fraction
      tbl%surface_z=table_z
      call findex(tbl%grid_x, n_opal92_x, hydrogen_fraction, x_index)
      if (x_index < 0) x_index=-x_index
      if (x_index >= 3) x_index=2
      if (x_index <= 0) then
         write(*,*) 'build_opal92_surface_table: surface X outside table X grid'
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
         if (density_start /= 1) then
            write(*,*) 'build_opal92_surface_table: surface-table density grid does not start at index 1'
            ierr = 1
            return
         end if
         density_end=density_start+tbl%density_count(row_index)-1
         do im3=density_start,density_end
            density_rhot3=tbl%grid_logr(im3)
            call interpolate_opal92_2d(tbl,temp6,density_rhot3,x_index,temp_index,rho_search_index, &
                 opacity0,log10_opacity0,dlnkap_dlnrho0,dlnkap_dlnt0, ierr)
            if (ierr /= 0) return
            call interpolate_opal92_2d(tbl,temp6,density_rhot3,x_index+1,temp_index,rho_search_index, &
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

end subroutine build_opal92_surface_table

!----------------------------------------------------------------------
! compute_opal92_opacity (formerly opal92_interp3d)
!----------------------------------------------------------------------
! THIS IS THE INTERPOLATION FACILITY FOR THE LIVERMORE (OPAL92)
! OPACITY TABLES USING CUBIC SPLINE INTERPOLATION. Companion routines
! in this directory: read_opal92_tables (read), prepare_opal92_tables
! (spline coefficients), build_opal92_surface_table (surface X,Z table),
! interpolate_opal92_2d (2D in T,rho).
!
! 3D interpolation (X, T, rho) in the OPAL92 opacity tables. Finds
! the nearest X table(s), converts to the table's internal (T6,
! rho/T6**3) coordinates, calls interpolate_opal92_2d for the 2D
! interpolation, and linearly extrapolates/interpolates in X if
! needed.
!
! 2026 wave 3 (R5): tbl is the table set (opacity_table%opal92(1) or
! (2)); the former opal92_interp3d_z2.f90 was a member-for-member
! rename of this file (its two stdout messages named it; its grdnt
! line was only continued differently) and is deleted.
subroutine compute_opal92_opacity(tbl, log10_density, log10_temperature, hydrogen_fraction, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ierr)

      type(opal92_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction
! 2026 ierr campaign: failures return via ierr (kap_eval gates).
      integer, intent(out) :: ierr
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

! tbl%index_x/index_t/index_rho (the warm-start cursors) default to
! 1 in opacity_tables.f90 (opal92_table_set).
      logical :: single_x_table
      real(dp) :: t6, rhot3
      integer :: im1
      real(dp) :: o0, ol0, qod0, qot0, o1, ol1, qod1, qot1
      real(dp) :: qodi, qoti, grdnt

      ierr = 0
      single_x_table = .false.
! INDEPENDENT PARAMETER IN LIVERMORE OPACITY TABLE;
! T6 = LOG10(T/10E6)
      t6 = log10_temperature - 6.0d0
! RHOT3 = LOG10(RHO/(T/10E6)**3)
      rhot3 = log10_density - 3.0d0*t6

      if (abs(tbl%surface_x-hydrogen_fraction) <= opal92_x_match_tol) then
         tbl%index_x = 4
         single_x_table = .true.
      endif
      if (.not. single_x_table) then
      do im1 = 1,n_opal92_x
         if (abs(tbl%grid_x(im1)-hydrogen_fraction) <= opal92_x_match_tol) then
            tbl%index_x = im1
            single_x_table = .true.
            exit
         endif
      end do
      end if
      if (.not. single_x_table) then
      call findex(tbl%grid_x, n_opal92_x, hydrogen_fraction, tbl%index_x)
      if (tbl%index_x < 0) tbl%index_x = -tbl%index_x
      if (tbl%index_x <= 1.and.rhot3 > -1.0d0) tbl%index_x = 2
      if (tbl%index_x >= 3) tbl%index_x = 2
      if (tbl%index_x <= 0) then
         write(*,*) 'compute_opal92_opacity: error in X grid'
         ierr = 1
         return
      end if
      end if
      call findex(tbl%grid_logt, tbl%num_temps, t6, tbl%index_t)
      if (tbl%index_t < 0.and.tbl%grid_logt(tbl%num_temps) == t6) tbl%index_t = -tbl%index_t
      if (tbl%index_t < 0) then
         write(*,*) 'compute_opal92_opacity: T out of table'
         ierr = 1
         return
      end if
      call interpolate_opal92_2d(tbl, t6, rhot3, tbl%index_x, tbl%index_t, tbl%index_rho, o0, ol0, qod0, qot0, ierr)
      if (ierr /= 0) return
      if (single_x_table) then
! USE ONLY ONE X TABLE
         log10_opacity = ol0
         opacity = o0
         qodi = qod0
         qoti = qot0
      else
! LINEAR EXTRAPOLATION IN X
         call interpolate_opal92_2d(tbl, t6, rhot3, tbl%index_x+1, tbl%index_t, tbl%index_rho, o1, ol1, qod1, qot1, ierr)
      if (ierr /= 0) return
         grdnt = (hydrogen_fraction-tbl%grid_x(tbl%index_x))/(tbl%grid_x(tbl%index_x+1)-tbl%grid_x(tbl%index_x))
         log10_opacity = (ol1-ol0)*grdnt + ol0
         qodi = (qod1-qod0)*grdnt + qod0
         qoti = (qot1-qot0)*grdnt + qot0
         opacity = exp10(log10_opacity)
      endif
! CONVERSION FROM THE DERIVATIVE WITH CONSTANT RHOT3 TO CONSTANT RHO
      dlnkap_dlnrho = qodi
      dlnkap_dlnt = qoti - 3.0d0*dlnkap_dlnrho

end subroutine compute_opal92_opacity

!----------------------------------------------------------------------
! interpolate_opal92_2d (formerly opal92_interp2d)
!----------------------------------------------------------------------
! TWO DIMENSIONAL INTERPOLATION FOR OPACITY AND THE DERIVATIVES
! (OPAL92 tables). abund_index, temp_index, and dens_index are the
! nearest grid point of abundance, temperature, and density found by
! compute_opal92_opacity. opacity is opacity, log10_opacity is log(opacity),
! dlnkap_dlnrho is the partial derivative of opacity wrt density,
! dlnkap_dlnt is the partial derivative of opacity wrt temperature.
!
! 2026 wave 3 (R5): tbl is the table set (opacity_table%opal92(1) or
! (2)); the former opal92_interp2d_z2.f90 was a member-for-member
! rename of this file (its two stdout messages named it) and is
! deleted.
subroutine interpolate_opal92_2d(tbl, temperature, density, abund_index, temp_index, &
     dens_index, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, ierr)


      type(opal92_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: temperature, density
! 2026 ierr campaign: interpolation/extrapolation failures return
! via ierr (kap_eval gates); the historical stops are gone.
      integer, intent(out) :: ierr
      integer, intent(in) :: abund_index, temp_index
      integer, intent(inout) :: dens_index
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt

      real(dp) :: xt(n_opal92_t), yto(n_opal92_t)
      real(dp) :: aqod(n_opal92_t)
      integer :: jt, it, its, itf
      integer :: mm1, index1, ndss, ndf, knot, index2
      real(dp) :: dx, c1, c2, c3, c4, ol0, qodi
      real(dp) :: ol00, unused_ddensity_dtemp

      ierr = 0
! FOR SIX GRID POINTS OF TEMPERATURE
      its = temp_index - 2
      if (its <= 0) its = 1
      itf = temp_index + 3
      if (itf > tbl%num_temps) itf = tbl%num_temps
      if (abund_index < 4) then
         mm1 = abund_index
      else
         mm1 = tbl%surface_x_index
      endif
      jt = 0
      do it = its,itf
         index1 = it + (mm1-1)*tbl%num_temps
         ndss = tbl%density_start_index(index1)
         if (ndss /= 1) then
            write(*,*) 'interpolate_opal92_2d: CHECK NDSS'
            ierr = 1
            return
         end if
         ndf = ndss + tbl%density_count(index1) - 1
         call findex(tbl%grid_logr, ndf, density, dens_index)
         if (dens_index < 0) then
! OUT SIDE THEN  LINEAR EXTRAPOLATION
            dens_index = -dens_index
            knot = dens_index - ndss + 1
            index2 = 4*(knot-1)
            dx = density - tbl%grid_logr(dens_index)
            if (abund_index < 4) then
               c1 = tbl%spline_coeffs(index1,index2+1)
               c2 = tbl%spline_coeffs(index1,index2+2)
            else
               c1 = tbl%surface_spline_coeffs(it,index2+1)
               c2 = tbl%surface_spline_coeffs(it,index2+2)
            endif
            ol0 = c2*dx + c1
            qodi = c2
         else
! IN SIDE THEN  SPLINE INTERPOLATION
            knot = dens_index - ndss + 1
            index2 = 4*(knot-1)
            dx = density - tbl%grid_logr(dens_index)
            if (abund_index < 4) then
               c1 = tbl%spline_coeffs(index1,index2+1)
               c2 = tbl%spline_coeffs(index1,index2+2)
               c3 = tbl%spline_coeffs(index1,index2+3)
               c4 = tbl%spline_coeffs(index1,index2+4)
            else
               c1 = tbl%surface_spline_coeffs(it,index2+1)
               c2 = tbl%surface_spline_coeffs(it,index2+2)
               c3 = tbl%surface_spline_coeffs(it,index2+3)
               c4 = tbl%surface_spline_coeffs(it,index2+4)
            endif
! INTERPOLATION FOR OPACITY(OL) IN THE ENTRY D AND THE EACH T-GRID
! ESTIMATES THE PARTIAL DERIVATIVE OF OL WRT D
! EVALUATES THE INTERPOLATION VALUE IN THE SUB-RANGE WE DETERMINED.
            ol0 = ((c4*dx+c3)*dx+c2)*dx + c1
            qodi = (3.0d0*c4*dx+2.0d0*c3)*dx + c2
         endif
!
         jt = jt + 1
         xt(jt) = tbl%grid_logt(it)
         yto(jt) = ol0
         aqod(jt) = qodi
      end do
      if (xt(1) > temperature.or.xt(jt) < temperature) then
         write(*,*) 'interpolate_opal92_2d: extrapolation fails'
         ierr = 1
         return
      end if
!! INTERPOLATION FOR THE OPACITY IN THE ENTRY T AND D.
!! GET THE PARTIAL DERIVATIVE OF OL WRT T.
      call intpol(xt, yto, jt, temperature, ol00, dlnkap_dlnt, ierr)
      if (ierr /= 0) return
      log10_opacity = ol00
      opacity = exp10(log10_opacity)
! QOTF = D LN(O)/D LN(T)
!! FIND THE PARTIAL DERIVATIVE VALUE OF OL WRT D IN THE GIVEN T AND D
      call intpol(xt, aqod, jt, temperature, dlnkap_dlnrho, unused_ddensity_dtemp, ierr)
      if (ierr /= 0) return
! QODF = D LN(O)/D LN(D)

end subroutine interpolate_opal92_2d

end module opal92_opacity
