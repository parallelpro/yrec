!----------------------------------------------------------------------
! mhd_eos
!----------------------------------------------------------------------
! The Mihalas-Hummer-Daeppen equation of state, from tables in
! (log P, log T) for several compositions.
!     read_mhd_tables  open the 3 ZAMS-type and up to 5 centre-type
!                      table files and read them (at eos_init)
!     compute_mhd_eos  density and the thermodynamic derivatives at
!                      (log T, log P, X, Z)
! Under these, Daeppen's routines under their original names: mhdpx,
! mhdpx1 and mhdpx2 interpolate the tables (mhd_table_reader reads
! them, with mhdst, mhdst1, rabu and rtab). The tables live in
! mhd_eos_tables (mhd).
module mhd_eos
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use eos_results, only: eos_result
      use mhd_eos_tables, only: i_mhd_chi_rho, i_mhd_chi_t, &
           i_mhd_dgrad_ad_dlog10_rho, i_mhd_dgrad_ad_dlog10_t, &
           i_mhd_dgrad_ad_dx, i_mhd_dlog10cp_dx, i_mhd_dlog10p_dx, &
           i_mhd_dlog10u_dx, i_mhd_dlogcp_dlogrho, i_mhd_dlogcp_dlogt, &
           i_mhd_eta, i_mhd_grad_ad, i_mhd_ion_frac_1, i_mhd_log10_cp, &
           i_mhd_log10_p, i_mhd_log10_pgas, i_mhd_log10_rho, i_mhd_log10_u, &
           i_mhd_placeholder, imhd_centre_1, imhd_centre_2, imhd_centre_3, &
           imhd_centre_4, imhd_centre_5, imhd_zams_a, imhd_zams_b, &
           imhd_zams_c, mhd, mhd_ivarc, mhd_ivarx, mhd_n_centre, mhd_n_zams, &
           mhd_nchem0, mhd_nr1m, mhd_nr2m, mhd_nrxm, mhd_nt1m, mhd_nt2m, &
           mhd_ntxm, imhd_family_zams_lower, imhd_family_zams_upper, imhd_family_centre
      use eos_messages, only: eos_log, eos_report, message_length, unwritten_line
      use numerics_lib, only: lir, quint
      use eos_mixture, only: compute_mean_weights
      use mhd_table_reader, only: read_mhd_tables
      use eos_configuration, only: eos_settings, eos_table_paths_t
      use eos_table_source, only: eos_table_source_t, stage_whole_eos, always_in_window, &
           whole_weight_inside, adopt_whole_result, blend_table_result
      use math_lib
      implicit none
      private
      public :: compute_mhd_eos, interpolate_mhd_table, mhd_table_source

contains

!----------------------------------------------------------------------
! mhd_table_source
!----------------------------------------------------------------------
! MHD as a table source (eos_table_source): selected by use_mhd_eos, its
! eight tables read at eos_init, and the whole EOS where it is selected
! (eos_get takes compute_mhd_eos's result alone; a point outside its
! tables is ierr).
function mhd_table_source() result(src)
      type(eos_table_source_t) :: src
      src%name = 'MHD'
      src%stage = stage_whole_eos
      src%selected => mhd_selected
      src%read_tables => read_mhd_source
      src%in_window => always_in_window
      src%lookup => mhd_lookup
      src%edge_weight => whole_weight_inside
      src%adopt => adopt_whole_result
      src%blend => blend_table_result
end function mhd_table_source

logical function mhd_selected()
      mhd_selected = eos_settings%use_mhd_eos
end function mhd_selected

subroutine read_mhd_source(paths, ierr)
      type(eos_table_paths_t), intent(in) :: paths
      integer, intent(out) :: ierr
      call read_mhd_tables(paths%mhd_zams, paths%mhd_centre, ierr)
end subroutine read_mhd_source

subroutine mhd_lookup(log10_temperature, log10_pressure, hydrogen_fraction, metal_fraction, &
     table_r, covers, ierr)
      real(dp), intent(in) :: log10_temperature, log10_pressure, hydrogen_fraction, metal_fraction
      type(eos_result), intent(inout) :: table_r
      logical, intent(out) :: covers
      integer, intent(out) :: ierr
      call compute_mhd_eos(log10_temperature, log10_pressure, hydrogen_fraction, metal_fraction, &
           table_r, ierr)
      covers = ierr == 0
end subroutine mhd_lookup

!----------------------------------------------------------------------
! compute_mhd_eos (formerly meqos)
!----------------------------------------------------------------------
! Mihalas, Hummer, Dappen equation of state. Implemented by YC,
! March 7, 1990. Given log10(T), log10(P), X, and Z, calls the MHD
! table interpolator (mhdpx) and converts its output
! (mhd%mhd_output) into the same thermodynamic derivative set
! used by the rest of the EOS machinery. A point outside the tables
! (log T outside their range, or log P above a row's highest) returns
! ierr = 1, with a message to the terminal and the run log (eos_report).
subroutine compute_mhd_eos(log10_temperature, log10_pressure, &
     hydrogen_fraction, metal_fraction, r, ierr)


      real(dp), intent(in) :: log10_temperature, log10_pressure
      type(eos_result), intent(inout) :: r
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction

      integer, intent(out) :: ierr

! how far R/mu from the table's P_gas may be from the mean weights' (a
! table result inconsistent with its own P: a stale or empty row)
      real(dp), parameter :: gas_constant_tolerance = 5.0d-7
! log10 e to the original's nine digits (clni is the exact 1/ln 10):
! turns the tables' d/dlog10 rho and d/dlog10 T of grad_ad into d/dln
      real(dp), parameter :: log10_e = 0.434294481d0
      real(dp) :: mhdpx_r10
      real(dp) :: chi_rho, chi_t, log10_specific_heat_cp
      real(dp) :: specific_gas_constant_check
      integer :: ion_idx
      character(len=message_length) :: eos_lines(4)

      ierr = 0

      r%temperature = exp10(log10_temperature)
      r%pressure = exp10(log10_pressure)
      call mhdpx(log10_pressure, log10_temperature, hydrogen_fraction, &
           mhdpx_r10, ierr)
      if (ierr /= 0) then
         eos_lines = unwritten_line
         write(eos_lines,'(a,2f10.5)') ' MHD EOS: (log T, log P) outside the tables or a bad table:', &
              log10_temperature, log10_pressure
         call eos_report(eos_lines)
         return
      end if
      r%log10_density = mhd%mhd_output(i_mhd_log10_rho)
      r%density = exp10(r%log10_density)
      r%dlnrho_dlnp = 1.0d0/mhd%mhd_output(i_mhd_chi_rho)
      chi_rho = mhd%mhd_output(i_mhd_chi_rho)
      chi_t = mhd%mhd_output(i_mhd_chi_t)
      r%dlnrho_dlnt = -chi_t/chi_rho
      log10_specific_heat_cp = mhd%mhd_output(i_mhd_log10_cp)
      r%specific_heat_cp = exp10(log10_specific_heat_cp)
      r%specific_heat_cp_dp = r%dlnrho_dlnp*mhd%mhd_output(i_mhd_dlogcp_dlogrho)
      r%specific_heat_cp_dt = mhd%mhd_output(i_mhd_dlogcp_dlogt) + mhd%mhd_output(i_mhd_dlogcp_dlogrho)*r%dlnrho_dlnt
      r%adiabatic_gradient = mhd%mhd_output(i_mhd_grad_ad)
      r%adiabatic_gradient_dp = r%dlnrho_dlnp*mhd%mhd_output(i_mhd_dgrad_ad_dlog10_rho)*log10_e/r%adiabatic_gradient
      r%adiabatic_gradient_dt = (mhd%mhd_output(i_mhd_dgrad_ad_dlog10_t) + mhd%mhd_output(i_mhd_dgrad_ad_dlog10_rho)*r%dlnrho_dlnt)* &
           log10_e/r%adiabatic_gradient
      r%dlnrho_dlnt_dp = r%dlnrho_dlnt*(r%adiabatic_gradient_dp - 1.0d0 + &
           r%dlnrho_dlnp + r%specific_heat_cp_dp)
      r%dlnrho_dlnt_dt = r%dlnrho_dlnt*(r%adiabatic_gradient_dt + 1.0d0 + &
           r%dlnrho_dlnt + r%specific_heat_cp_dt)
      r%beta = exp10((mhd%mhd_output(i_mhd_log10_pgas) - mhd%mhd_output(i_mhd_log10_p)))
      r%beta_inverse = 1.0d0/r%beta
      r%beta14 = 4.0d0*r%beta_inverse*(1.0d0 - r%beta)
      do ion_idx = 1, 3
      r%ion_fraction(ion_idx) = mhd%mhd_output(ion_idx+i_mhd_ion_frac_1-1)
      end do
      r%electron_degeneracy_parameter = mhd%mhd_output(i_mhd_eta)
      r%specific_gas_constant = exp10(mhd%mhd_output(i_mhd_log10_pgas))/(r%density*r%temperature)
      call compute_mean_weights(r%temperature, r%pressure, r%density, hydrogen_fraction, &
           metal_fraction, eos_settings%use_diffusion_z, specific_gas_constant_check, &
           r%ion_mean_weight_inverse, r%electron_mean_weight_inverse, r%beta)
      if (abs((r%specific_gas_constant-specific_gas_constant_check)/ &
           r%specific_gas_constant) > gas_constant_tolerance) then
         eos_lines = unwritten_line
         write(eos_lines(1),'(a,2es24.15)') 'MHD EOS: the gas constant of the mean weights ' // &
              'is not the table''s P_gas/(rho T) (table, mean weights):', r%specific_gas_constant, &
              specific_gas_constant_check
         write(eos_lines(2),'(a)') '     check the composition'
         call eos_report(eos_lines)
          ierr = 1
          return
      end if
end subroutine compute_mhd_eos


!----------------------------------------------------------------------
! mhdpx
!----------------------------------------------------------------------
!
!     MHD interpolation package
!     pressure as argument interpolation in X, Z fixed
!
!     Input   PGL = log10 (pressure [dyn/cm2]) (Daeppen's header says
!             the gas pressure; compute_mhd_eos passes the total, an open
!             question: eos/README.md section 2.2)
!             TL  = log10 (temperature [K])
!             XC  = hydrogen abundance (by mass)
!
!     Output  RL  = log10 (density [g/cm3])        : Argument
!             OTHER (SEE SEPARATE INSTRUCTIONS)    : mhd%mhd_output
!
!     Error   ierr = 1 signals PGL,TL outside the domain of tables
!             (log T outside the tables' range, or log P above the
!             highest of a row the interpolation uses), or
!             non-equidistant ZAMS tables; otherwise ierr = 0.
!
! Top-level entry point to the MHD (Mihalas, Hummer, Dappen)
! equation-of-state table interpolator. S/R MHDST (via read_mhd_tables) must
! be called in main to load the tables first. Dispatches to the
! variable-X interpolation in mhdpx1, then returns log10(density).
subroutine mhdpx(log10_pressure, log10_temperature, hydrogen_fraction, &
     log10_density, ierr)

      real(dp), intent(in) :: log10_pressure, log10_temperature, &
           hydrogen_fraction
      real(dp), intent(out) :: log10_density
      integer, intent(out) :: ierr

      ierr = 0

      call mhdpx1(log10_pressure, log10_temperature, hydrogen_fraction, ierr)
      if (ierr /= 0) return
      log10_density = mhd%mhd_output(i_mhd_log10_rho)
end subroutine mhdpx

!----------------------------------------------------------------------
! mhdpx1
!----------------------------------------------------------------------
!     The tables (read by mhdst at eos_init), with different X and fixed Z:
!     TL < TLIM1:         Lower part of ZAMS tables
!     TLIM1 < TL < TLIM2: Upper part of ZAMS tables
!     TL > TLIM2:         Centre tables
!     TMINI,TMAXI:        Temperature interval covered
!                         by the tables
!     (TLIM1/TLIM2/TMINI/TMAXI are mhd%zams_lower_upper_boundary_log10t,
!     zams_centre_boundary_log10t, table_log10t_min, table_log10t_max.)
!     A point outside the tables returns ierr = 1 (mhd%mhd_output is
!     then undefined). The X interpolation follows the family the T
!     selects: quadratic across the three ZAMS compositions, cubic
!     across four of the five centre ones.
subroutine mhdpx1(log10_pressure, log10_temperature, hydrogen_fraction, ierr)

      real(dp), intent(in) :: log10_pressure, log10_temperature, &
           hydrogen_fraction

      integer, intent(out) :: ierr
      integer, parameter :: ivarx = mhd_ivarx
      integer, parameter :: ndimt = 8
! how far the ZAMS tables' X may be from equidistant for the quadratic
      real(dp), parameter :: x_spacing_tolerance = 1.0d-4

      real(dp) :: table_vars(ndimt,ivarx), table_hfrac(ndimt)
      logical :: in_range, use_centre_tables

!     Quantities for interpolation in X
      real(dp) :: cubic_vars(ivarx,4), cubic_x_nodes(4)
      integer :: cubic_table_index(4)
! the X interpolation's lir settings (interp_mode: lir's flag, unused)
      integer :: i, family, n_tables, iv, ixmin, num_vars, var_leading_dim, &
           num_points, interp_mode
      real(dp) :: x_grid_origin, x_grid_spacing, &
           var_at_x0, var_at_x1, var_at_x2
      character(len=message_length) :: eos_lines(4)

      ierr = 0

      if (log10_temperature < mhd%table_log10t_min .or. &
           log10_temperature > mhd%table_log10t_max) then
         ierr = 1
         return
      end if
      table_vars = 0.0d0
      table_hfrac = 0.0d0
      use_centre_tables = log10_temperature >= mhd%zams_centre_boundary_log10t
! the family the temperature selects, each of its compositions in turn
      if (log10_temperature < mhd%zams_lower_upper_boundary_log10t) then
         family = imhd_family_zams_lower
         n_tables = mhd_n_zams
      else if (.not. use_centre_tables) then
         family = imhd_family_zams_upper
         n_tables = mhd_n_zams
      else
         family = imhd_family_centre
         n_tables = mhd_n_centre
      end if
      do i = 1, n_tables
         call mhdpx2(log10_pressure, log10_temperature, family, i, table_vars, table_hfrac, &
              ndimt, in_range, ierr)
         if (ierr /= 0) return
         if (.not. in_range) exit
      end do
      if (.not. in_range) then
         ierr = 1
         return
      end if
!     Interpolation in X
      if (.not. use_centre_tables) then
!        Quadratic Newton (equidistant XC's)
         x_grid_origin = table_hfrac(1)
         x_grid_spacing  = table_hfrac(2) - table_hfrac(1)
         if (abs(table_hfrac(3)-table_hfrac(2)-x_grid_spacing) > x_spacing_tolerance) then
            eos_lines = unwritten_line
            write(eos_lines(1),'(a,3es14.6)') 'MHD EOS: the ZAMS tables'' X are not equidistant: X =', &
                 table_hfrac(1), table_hfrac(2), table_hfrac(3)
            call eos_report(eos_lines)
            ierr = 1
            return
         end if
         do iv=1,ivarx
         var_at_x0 = table_vars(1,iv)
         var_at_x1 = table_vars(2,iv)
         var_at_x2 = table_vars(3,iv)
         call quint(hydrogen_fraction, x_grid_origin, x_grid_spacing, &
              var_at_x0, var_at_x1, var_at_x2, mhd%mhd_output(iv))
         end do
      else
!     Cubic Lagrangian (arbitrarily spaced XC's)
!     use 4 highest XC(I)'s, if possible.
!     Distinguish case where the XC(I), I=4..8,
!     increase or decrease.
         if (table_hfrac(4) < table_hfrac(8)) then
            ixmin = 5
            if (hydrogen_fraction < table_hfrac(5)) ixmin=4
         else
            ixmin = 4
            if (hydrogen_fraction < table_hfrac(7)) ixmin=5
         end if
         do i=1,4
            cubic_table_index(i)=ixmin+i-1
         end do
         do i =1,4
            cubic_x_nodes(i) = table_hfrac(cubic_table_index(i))
            do iv=1,ivarx
               cubic_vars(iv,i) = table_vars(cubic_table_index(i),iv)
            end do
         end do
         num_vars=ivarx
         var_leading_dim=ivarx
         num_points=4
         call lir(hydrogen_fraction, cubic_x_nodes, mhd%mhd_output, cubic_vars, &
              num_vars, var_leading_dim, num_points, interp_mode, ierr)
      end if
end subroutine mhdpx1

!----------------------------------------------------------------------
! mhdpx2
!----------------------------------------------------------------------
! Interpolates table ic of a family (imhd_family_zams_lower or _upper,
! compositions A, B, C; imhd_family_centre, compositions 1..5) at
! (log10 P, log10 T) and stores the result as row ic (ZAMS) or
! mhd_n_zams + ic (centre) of table_vars, with that table's hydrogen
! mass fraction in table_hfrac; in_range is false when log P is above
! the highest of a row the interpolation uses, and ierr lir's status
! (numerics_status). The tables are the composition-indexed arrays of
! mhd (mhd_eos_tables).
subroutine mhdpx2(log10_pressure, log10_temperature, family, ic, &
     table_vars, table_hfrac, ndimt, in_range, ierr)
      integer, parameter :: ivarx = mhd_ivarx
      integer, intent(in) :: family, ic, ndimt
      real(dp), intent(in) :: log10_pressure, log10_temperature
      real(dp), intent(inout) :: table_vars(:,:), &
           table_hfrac(:)
      logical, intent(out) :: in_range
      integer, intent(out) :: ierr
      integer, parameter :: ivarc = mhd_ivarc
      integer, parameter :: nt1m = mhd_nt1m
      integer, parameter :: nt2m = mhd_nt2m
      integer, parameter :: ntxm = mhd_ntxm
      integer, parameter :: nr1m = mhd_nr1m
      integer, parameter :: nr2m = mhd_nr2m
      integer, parameter :: nrxm = mhd_nrxm
      real(dp) :: intpt_work1(ivarx,4), intpt_work2(ivarx,4), &
           intpt_y(ivarx), interpolated_vars(ivarx)

      integer :: i

      interpolated_vars = 0.0d0
      in_range = .true.
      ierr = 0
      select case (family)
      case (imhd_family_zams_lower)
         call interpolate_and_store(mhd%zams_lower_table(:,:,:,ic), &
              nt1m, nr1m, ivarc, mhd%zams_lower_log10t, mhd%zams_lower_num_t, &
              mhd%zams_lower_num_r, ic, mhd%zams_mass_fraction(1,ic))
      case (imhd_family_zams_upper)
         call interpolate_and_store(mhd%zams_upper_table(:,:,:,ic), &
              nt2m, nr2m, ivarc, mhd%zams_upper_log10t, mhd%zams_upper_num_t, &
              mhd%zams_upper_num_r, ic, mhd%zams_mass_fraction(1,ic))
      case (imhd_family_centre)
         call interpolate_and_store(mhd%centre_table(:,:,:,ic), &
              ntxm, nrxm, ivarx, mhd%centre_log10t, mhd%centre_num_t, &
              mhd%centre_num_r, mhd_n_zams + ic, mhd%centre_mass_fraction(1,ic))
      end select

contains

! Interpolate one table at (log10 P, log10 T) and store its num_vars
! results as row `row` of table_vars, with the table's hydrogen mass
! fraction in table_hfrac(row). The work arrays and the (P,T) point are
! host-associated from mhdpx2.
subroutine interpolate_and_store(table_data, table_dim_t, table_dim_r, &
     num_vars, table_log10t, num_t, num_r, row, hydrogen_mass_fraction)
      integer, intent(in) :: table_dim_t, table_dim_r, num_vars, num_t, num_r, row
      real(dp), intent(in) :: table_data(:,:,:)
      real(dp), intent(in) :: table_log10t(:)
      real(dp), intent(in) :: hydrogen_mass_fraction

      call interpolate_mhd_table(log10_pressure, log10_temperature, table_data, &
           table_dim_t, table_dim_r, num_vars, table_log10t, num_t, &
           num_r, intpt_work1(1:num_vars,:), intpt_work2(1:num_vars,:), intpt_y(1:num_vars), &
           interpolated_vars(1:num_vars), in_range, ierr)
      do i=1,num_vars
         table_vars(row,i) = interpolated_vars(i)
      end do
      table_hfrac(row)    = hydrogen_mass_fraction
end subroutine interpolate_and_store

end subroutine mhdpx2

!----------------------------------------------------------------------
! interpolate_mhd_table (formerly intpt, numerics_lib)
!----------------------------------------------------------------------
! Bicubic interpolation in one MHD table at (log10_pressure,
! log10_temperature): table_data(table_dim_t,table_dim_r,num_vars) on
! the table_log10t rows, with column 2 holding log P. Two passes of
! numerics' 4-point lir: at each of the 4 bracketing T rows, in log P
! (work1 -> work2); then across the 4 rows in log T (work2 ->
! interp_vars). Above the highest log P of a row it returns with
! in_range false and interp_vars unset; below the lowest, and outside
! the T rows, the stencils clamp to the table's edge.
subroutine interpolate_mhd_table(log10_pressure, log10_temperature, table_data, &
     table_dim_t, table_dim_r, num_vars, table_log10t, num_t, num_r, &
     work1, work2, y_work, interp_vars, in_range, ierr)
      integer, intent(in) :: table_dim_t, table_dim_r, num_vars
      real(dp), intent(in) :: log10_pressure, log10_temperature
      real(dp), intent(in) :: table_data(:,:,:)
      real(dp), intent(in) :: table_log10t(:)
      integer, intent(in) :: num_t, num_r
! work arrays of num_vars rows (lir reads them so: num_vars is its
! leading dimension)
      real(dp), intent(inout) :: work1(:,:), work2(:,:)
      real(dp), intent(inout) :: y_work(:)
      real(dp), intent(out) :: interp_vars(:)
      logical, intent(out) :: in_range
! lir's status (numerics_status): nonzero for coincident table nodes
      integer, intent(out) :: ierr

      integer :: r_lo_guess(4), r_indices(4,4), t_indices(4)
      real(dp) :: x_nodes(4)

      integer :: n, i, m, j, t_col, iv, t_idx, r_idx, t_idx_max, r_idx_max
      real(dp) :: p_max
      integer :: lir_num_vars, lir_leading_dim, lir_num_points, lir_interp_mode

! Bracket the temperature: the original's GOTO 101 out of this scan
! is an EXIT, not a RETURN (a RETURN left every output unset for any
! in-range temperature).  t_indices(1) starts at 1 so a temperature
! below the first table row clamps to the bottom stencil instead of
! reading a stale/zero index.
      in_range = .true.
      ierr = 0
      t_indices(1)=1
      do n=1,num_t
         if (table_log10t(n) >= log10_temperature) exit
         t_indices(1)=n
      end do
      if (t_indices(1) >= 2) t_indices(1)=t_indices(1)-1
      t_idx_max=num_t-3
      if (t_indices(1) > t_idx_max) t_indices(1)=t_idx_max
      do i=2,4
         t_indices(i)=t_indices(1)+i-1
      end do
      do i=1,4
         r_lo_guess(i)=1
         t_idx=t_indices(i)
         p_max=table_data(t_idx,num_r,2)
         if (log10_pressure > p_max) then
            in_range = .false.
            return
         end if
         do m=1,num_r
            if (table_data(t_idx,m,2) >= log10_pressure) exit
            r_lo_guess(i)=m
         end do
         if (r_lo_guess(i) >= 2) r_lo_guess(i)=r_lo_guess(i)-1
         r_idx_max=num_r-3
         if (r_lo_guess(i) > r_idx_max) r_lo_guess(i)=r_idx_max
      end do
      do i=1,4
         do j=1,4
            r_indices(j,i)=r_lo_guess(i)+j-1
         end do
      end do
! lir call settings, shared by both interpolation passes
      lir_num_vars=num_vars
      lir_leading_dim=num_vars
      lir_num_points=4
      do t_col=1,4
         t_idx=t_indices(t_col)
         do i=1,4
            r_idx=r_indices(i,t_col)
            x_nodes(i)=table_data(t_idx,r_idx,2)
         end do
         do i=1,4
            r_idx=r_indices(i,t_col)
            do iv=1,num_vars
               work1(iv,i)=table_data(t_idx,r_idx,iv)
            end do
         end do

         call lir(log10_pressure, x_nodes, y_work, work1, lir_num_vars, &
              lir_leading_dim, lir_num_points, lir_interp_mode, ierr)
         if (ierr /= 0) return
         do iv=1,num_vars
            work2(iv,t_col) = y_work(iv)
         end do

      end do
      do i=1,4
         t_idx=t_indices(i)
         x_nodes(i) = table_log10t(t_idx)
      end do

      call lir(log10_temperature, x_nodes, interp_vars, work2, lir_num_vars, &
           lir_leading_dim, lir_num_points, lir_interp_mode, ierr)

      return
end subroutine interpolate_mhd_table
!======================================================================
! Sources: the MHD EOS of Hummer & Mihalas (1988), Mihalas, Daeppen &
! Hummer (1988) and Daeppen et al. (1988); the interpolation package and
! its names (mhdpx, mhdpx1, mhdpx2, intpt) are Daeppen's, wrapped for
! YREC by YC (MEQOS, March 1990).

end module mhd_eos
