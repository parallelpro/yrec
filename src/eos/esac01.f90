!----------------------------------------------------------------------
! esac01
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original esac01.f; only variable names, source form, and comment
! style were updated.
!
! OPAL 2001 EOS analogue of esac.f90 (see there for the general
! description of the interpolation scheme). Differences from the
! 1995 version are preserved verbatim, including: iorder capped at 9
! (not 10, since the 2001 tables only carry 9 thermodynamic
! variables); the extra iqu/ipu pre-checks that gate whether the
! "3x4"/"4x4" boundary sums are even accumulated; and the boundary
! reset checking both l3==nr and k3==nt (1995's esac.f only checks
! l3==nr).
!
!   eos_output(1) pressure, in megabars (1e12 dyne/cm**2)
!   eos_output(2) energy, in 1e12 erg/gm (zero at T6=0)
!   eos_output(3) entropy, in units of energy/T6
!   eos_output(4) specific heat, dE/dT6 at constant volume
!   eos_output(5) dlogP/dlogRho at constant T6 (Cox & Guili eq 9.82)
!   eos_output(6) dlogP/dlogT6 at constant Rho (Cox & Guili eq 9.81)
!   eos_output(7) gamma1 (Cox & Guili eq 9.88)
!   eos_output(8) gamma2/(gamma2-1) (Cox & Guili eq 9.88)
!   eos_output(9) gamma3-1 (Cox & Guili eq 9.88)
subroutine esac01(hydrogen_fraction, t6_temperature, density, &
     deriv_order, rad_flag, *)

      use luout_lib
      implicit none

      double precision, intent(in) :: hydrogen_fraction, t6_temperature, &
           density
      integer, intent(in) :: deriv_order, rad_flag

      integer, parameter :: mx = 5, mv = 10, nr = 169, nt = 191

! common/lreadco/: table_loaded_flag is shared (by COMMON block name)
! with the OPAL 1995/2006 readers -- see esac.f90's header.
      integer :: table_loaded_flag
      common/lreadco/ table_loaded_flag

! common/eeeos/: both members used here.
!   x_interp_workspace holds the per-variable table slice used for
!     the X (hydrogen_fraction) interpolation.
!   x_grid_copy is an auxiliary copy of x_grid, set from it below.
      double precision :: x_interp_workspace(mx,nt,nr), x_grid_copy(mx)
      common/eeeos/ x_interp_workspace, x_grid_copy

! common/aaeos/: not used directly here (consumed by t6rinteos01.f90).
      double precision :: rho_interp_hi(4), rho_interp_lo(4), xxh
      common/aaeos/ rho_interp_hi, rho_interp_lo, xxh

! common/aeos/: the main OPAL 2001 EOS table and its interpolation
! grids/scratch arrays.
!   eos_table(x,var,t6,rho): the thermodynamic-variable table. arg 1
!     selects x (x_grid), arg 2 selects the variable (matching
!     eos_output's default ordering), arg 3 selects T6 (t6_grid /
!     t6_list), arg 4 selects density (density_grid).
!   t6_list(density_row, t6_row): T6 grid, varies per density row
!     (the table is homogeneous in T6 per row, but the number of
!     valid rows shrinks with density -- see t6_index_lo).
!   density_grid: tabulated density values, g/cc.
!   t6_grid: auxiliary 1-D copy of t6_list's first column.
!   x_grid: tabulated hydrogen-fraction values (0.0, 0.2, ..., 0.8).
!   x_grid_spacing_inv/t6_grid_spacing_inv/density_grid_spacing_inv:
!     1/(grid(i)-grid(i-1)) for each respective grid.
      double precision :: eos_table(mx,mv,nt,nr), t6_list(nr,nt), &
           density_grid(nr), t6_grid(nt), x_interp_result(nt,nr), &
           x_interp_result_alt(nt,nr), x_grid_spacing_inv(mx), &
           t6_grid_spacing_inv(nt), density_grid_spacing_inv(nr)
      integer :: x_loop_index, x_index_lo
      double precision :: x_grid(mx)
      common/aeos/ eos_table, t6_list, density_grid, t6_grid, &
           x_interp_result, x_interp_result_alt, x_grid_spacing_inv, &
           t6_grid_spacing_inv, density_grid_spacing_inv, x_loop_index, &
           x_index_lo, x_grid

! common/beos/: z_table is the table's Z (metal fraction), used for
! consistency checks. eos_var_order defines the order thermodynamic
! variables are stored in eos_table/eos_output; eos_index_inverse is
! its inverse mapping. t6_index_lo(density_row) is the index in
! t6_grid of the lowest available temperature at that density row.
      double precision :: z_table(mx)
      integer :: eos_index_inverse(10), eos_var_order(10), t6_index_lo(nr)
      common/beos/ z_table, eos_index_inverse, eos_var_order, t6_index_lo

! common/bbeos/: density index window (l1..l4), t6 index window
! (k1..k4), and the interpolation-order flags t6_interp_order
! (original ip) and density_interp_order (original iq).
      integer :: density_index_1, density_index_2, density_index_3, &
           density_index_4, t6_index_1, t6_index_2, t6_index_3, &
           t6_index_4, t6_interp_order, density_interp_order
      common/bbeos/ density_index_1, density_index_2, density_index_3, &
           density_index_4, t6_index_1, t6_index_2, t6_index_3, &
           t6_index_4, t6_interp_order, density_interp_order

! common/eeos/: eos_output is the caller-facing interpolated result.
      double precision :: esact, eos_output(mv)
      common/eeos/ esact, eos_output


      double precision :: species_mass_fraction(7)
      double precision :: molar_gas_constant_mbcc
      character(len=1) :: blank_line
      character(len=15) :: routine_id

      integer :: fill_idx

! NOTE: this literal has no D0 suffix in the original (esac01.f), unlike
! esac.f90's analogous 83.1446304d0 -- preserved verbatim: it is parsed
! as a single-precision constant and then widened to double precision,
! which is NOT bit-identical to the correctly-rounded double value.
      data molar_gas_constant_mbcc/83.14511/
! NOTE: no D-suffix in the original (esac01.f) -- preserved verbatim:
! these are parsed as single-precision constants and then widened to
! double precision (0.2/0.4/0.6/0.8 are not exactly representable in
! binary floating point, so this is NOT bit-identical to the
! correctly-rounded double values).
      data (x_grid(fill_idx), fill_idx=1,mx) /0.0, 0.2, 0.4, 0.6, 0.8/
      data (eos_var_order(fill_idx), fill_idx=1,10) /1,2,3,4,5,6,7,8,9,10/
      data (t6_index_lo(fill_idx), fill_idx=1,nr) &
           /92*191, 190, 189, 188, 187, 186, 185, 184, 174, 4*134, 3*133, &
           2*132, 98, 92, 2*85, 2*77, 71, 3*63, 2*59, 53, 51, 2*46, 9*44, &
           3*38, 6*33, 16*29, 27, 26, 25, 23, 22, 20, 19, 18, 17, 16/
      data routine_id/"OPALEOS/ESAC01:"/
      save

! --- locals ---
      double precision :: hydrogen_fraction_copy, density_copy   ! xxi/ri: assigned but not used further
      double precision :: t6_value, density_value   ! working copies (slt/slr)
      double precision :: table_metal_fraction   ! z: set once, save'd across calls
      integer :: species_idx, index_idx
      integer :: lo_idx, hi_idx, mid_idx, result_idx
      integer :: x_index_2, x_index_3, x_index_4, x_index_hi
      integer :: t6_order_hi, density_order_hi   ! ipu/iqu
      double precision :: table_sum_2x2, table_sum_3x3, table_sum_3x4, &
           table_sum_4x4
      integer :: eos_var_idx, x_print_idx
      integer :: density_scan_idx, t6_scan_idx
      integer :: recompute_flag, cache_slot
      double precision :: x_interp_weight
      double precision :: pressure_scale
      double precision :: total_moles, ground_state_energy, &
           metal_mole_fraction, mean_molecular_weight
      double precision, external :: quadeos01, gmass01

      blank_line = ' '
      if (deriv_order.gt.9) then
         write (short_file_unit,'(A, " iorder cannot exceed 9")') routine_id
      end if
      if ((rad_flag.ne.0) .and. (rad_flag.ne.1)) then
         write (short_file_unit,'(A, " Irad must be 0 or 1")') routine_id
         stop
      end if

      hydrogen_fraction_copy = hydrogen_fraction
      density_copy = density
!
      t6_value = t6_temperature
      density_value = density
!
      if (table_loaded_flag.ne.12345678) then
         table_loaded_flag = 12345678
         do index_idx = 1, 10
            do species_idx = 1, 10
               if (eos_var_order(index_idx).eq.species_idx) &
                    eos_index_inverse(index_idx) = species_idx
            end do
         end do
         do species_idx = 1, mx
            x_grid_copy(species_idx) = x_grid(species_idx)
         end do
!
! ..... read the data files
         call readcoeos01
         table_metal_fraction = z_table(1)

         if (table_metal_fraction+hydrogen_fraction-1.0d-6.gt.1) go to 61
      end if
!
!
! ..... Determine T6,rho grid points to use in the
!       interpolation.
      if ((t6_value.gt.t6_grid(1)) .or. (t6_value.lt.t6_grid(nt))) go to 62
      if ((density_value.lt.density_grid(1)) .or. &
           (density_value.gt.density_grid(nr))) go to 62
!
!
!
      lo_idx = 2
      hi_idx = mx
    8 if (hi_idx-lo_idx.gt.1) then
         mid_idx = (hi_idx+lo_idx)/2
         if (hydrogen_fraction.le.x_grid(mid_idx)+1.0d-7) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
         go to 8
      end if
      result_idx = hi_idx
      x_index_lo = result_idx - 2
      x_index_2 = result_idx - 1
      x_index_3 = result_idx
      x_index_4 = result_idx + 1
      x_index_hi = x_index_4
      if (hydrogen_fraction.lt.1.0d-6) then
         x_index_lo = 1
         x_index_2 = 1
         x_index_3 = 1
         x_index_4 = 2
         x_index_hi = 1
      end if
      if ((hydrogen_fraction.le.x_grid(2)+1.0d-7) .or. &
           (hydrogen_fraction.ge.x_grid(mx-2)-1.0d-7)) x_index_hi = x_index_3
!
      lo_idx = 2
      hi_idx = nr
   12 if (hi_idx-lo_idx.gt.1) then
         mid_idx = (hi_idx+lo_idx)/2
         if (density_value.eq.density_grid(mid_idx)) then
            hi_idx = mid_idx
            go to 13
         end if
         if (density_value.le.density_grid(mid_idx)) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
         go to 12
      end if
   13 result_idx = hi_idx
      density_index_1 = result_idx - 2
      density_index_2 = result_idx - 1
      density_index_3 = result_idx
      density_index_4 = density_index_3 + 1
      density_order_hi = 3
      if (density_index_4.gt.nr) density_order_hi = 2
!
      lo_idx = nt
      hi_idx = 2
   11 if (lo_idx-hi_idx.gt.1) then
         mid_idx = (hi_idx+lo_idx)/2
         if (t6_temperature.eq.t6_list(1,mid_idx)) then
            lo_idx = mid_idx
            go to 14
         end if
         if (t6_temperature.le.t6_list(1,mid_idx)) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
         go to 11
      end if
   14 result_idx = lo_idx
      t6_index_1 = result_idx - 2
      t6_index_2 = result_idx - 1
      t6_index_3 = result_idx
      t6_index_4 = t6_index_3 + 1
      t6_order_hi = 3
      if (t6_index_4.gt.nt) t6_order_hi = 2
      if (t6_index_3.eq.0) then
         write (short_file_unit,'(A, " ihi,ilo,imd",3I5)') routine_id
      end if

!     check to determine if interpolation indices fall within
!     table boundaries.  choose largest allowed size.
      table_sum_2x2 = 0.0d0
      table_sum_3x3 = 0.0d0
      table_sum_3x4 = 0.0d0
      table_sum_4x4 = 0.0d0
      do x_loop_index = x_index_lo, x_index_hi
         do density_scan_idx = density_index_1, density_index_1+1
            do t6_scan_idx = t6_index_1, t6_index_1+1
               table_sum_2x2 = table_sum_2x2 + &
                    eos_table(x_loop_index,1,t6_scan_idx,density_scan_idx)
            end do
         end do
         do density_scan_idx = density_index_1, density_index_1+2
            do t6_scan_idx = t6_index_1, t6_index_1+2
               table_sum_3x3 = table_sum_3x3 + &
                    eos_table(x_loop_index,1,t6_scan_idx,density_scan_idx)
            end do
         end do
         if (t6_order_hi.eq.3) then
            do density_scan_idx = density_index_1, density_index_1+2
               do t6_scan_idx = t6_index_1, t6_index_1+t6_order_hi
                  table_sum_3x4 = table_sum_3x4 + &
                       eos_table(x_loop_index,1,t6_scan_idx,density_scan_idx)
               end do
            end do
         else
            table_sum_3x4 = 2.0d+30
         end if
         if (density_order_hi.eq.3) then
            do density_scan_idx = density_index_1, density_index_1+3
               do t6_scan_idx = t6_index_1, t6_index_1+t6_order_hi
                  table_sum_4x4 = table_sum_4x4 + &
                       eos_table(x_loop_index,1,t6_scan_idx,density_scan_idx)
               end do
            end do
         else
            table_sum_4x4 = 2.0d+30
         end if
      end do
      density_interp_order = 2
      t6_interp_order = 2
      if (table_sum_3x3.gt.1.0d+30) then
         if (table_sum_2x2.lt.1.0d+25) then
            t6_index_1 = t6_index_3 - 3
            t6_index_2 = t6_index_1 + 1
            t6_index_3 = t6_index_2 + 1
            density_index_1 = density_index_3 - 3
            density_index_2 = density_index_1 + 1
            density_index_3 = density_index_2 + 1
            go to 15
         else
            go to 65
         end if
      end if
      if (table_sum_3x4.lt.1.0d+30) t6_interp_order = 3
      if (table_sum_4x4.lt.1.0d+30) density_interp_order = 3

      if (t6_temperature.ge.t6_list(1,2)+1.0d-7) t6_interp_order = 2
      if (density_value.le.density_grid(2)+1.0d-15) density_interp_order = 2

      if ((density_index_3.eq.nr) .or. (t6_index_3.eq.nt)) then
         density_interp_order = 2
         t6_interp_order = 2
      end if

   15 continue
      do 124 eos_var_idx = 1, deriv_order
      do 123 x_loop_index = x_index_lo, x_index_hi

      recompute_flag = 0

!__________
      do density_scan_idx = density_index_1, density_index_1+density_interp_order
         do t6_scan_idx = t6_index_1, t6_index_1+t6_interp_order
            x_interp_workspace(x_loop_index,t6_scan_idx,density_scan_idx) = &
                 eos_table(x_loop_index,eos_var_idx,t6_scan_idx,density_scan_idx)
            recompute_flag = 1
         end do
      end do
  123 continue
      if ((z_table(x_index_2).ne.z_table(x_index_lo)) .or. &
           (z_table(x_index_3).ne.z_table(x_index_lo))) then
         write(short_file_unit,'(A,"Z does not match Z in OEOS01 files you are" &
              &," using")') routine_id
         stop
      end if
      if (table_metal_fraction.ne.z_table(x_index_lo)) go to 66
      recompute_flag = 0
      cache_slot = 1
      do 45 density_scan_idx = density_index_1, density_index_1+density_interp_order
         do t6_scan_idx = t6_index_1, t6_index_1+t6_interp_order
            if (x_index_hi.eq.1) then
               x_interp_result(t6_scan_idx,density_scan_idx) = &
                    x_interp_workspace(x_index_lo,t6_scan_idx,density_scan_idx)
               go to 46
            end if
            x_interp_result(t6_scan_idx,density_scan_idx) = &
                 quadeos01(recompute_flag, cache_slot, hydrogen_fraction, &
                 x_interp_workspace(x_index_lo,t6_scan_idx,density_scan_idx), &
                 x_interp_workspace(x_index_2,t6_scan_idx,density_scan_idx), &
                 x_interp_workspace(x_index_3,t6_scan_idx,density_scan_idx), &
                 x_grid_copy(x_index_lo), x_grid_copy(x_index_2), &
                 x_grid_copy(x_index_3))
            if (x_interp_result(t6_scan_idx,density_scan_idx).gt.1.0d+20) then
               write(short_file_unit,'(A," problem it ir,l3,k3,iq,ip=", 6I5)') &
                    routine_id, t6_scan_idx, density_scan_idx, density_index_3, &
                    t6_index_3, density_interp_order, t6_interp_order
               write(short_file_unit,'(3E12.4)') &
                    (x_interp_workspace(x_print_idx,t6_scan_idx,density_scan_idx), &
                    x_print_idx=x_index_lo,x_index_lo+2)
            end if
            recompute_flag = 1
   46       continue
         end do
   45 continue

      if (x_index_4.eq.x_index_hi) then  ! interpolate between quadratics
      recompute_flag = 0
      cache_slot = 1
      x_interp_weight = (x_grid_copy(x_index_3) - hydrogen_fraction)* &
           x_grid_spacing_inv(x_index_3)
      do 47 density_scan_idx = density_index_1, density_index_1+density_interp_order
         do t6_scan_idx = t6_index_1, t6_index_1+t6_interp_order
            x_interp_result_alt(t6_scan_idx,density_scan_idx) = &
                 quadeos01(recompute_flag, cache_slot, hydrogen_fraction, &
                 x_interp_workspace(x_index_2,t6_scan_idx,density_scan_idx), &
                 x_interp_workspace(x_index_3,t6_scan_idx,density_scan_idx), &
                 x_interp_workspace(x_index_4,t6_scan_idx,density_scan_idx), &
                 x_grid_copy(x_index_2), x_grid_copy(x_index_3), &
                 x_grid_copy(x_index_4))
            if (x_interp_result(t6_scan_idx,density_scan_idx).gt.1.0d+20) then
               write(short_file_unit,'(A," problem it ir,l3,k3,iq,ip=", 6I5)') &
                    routine_id, t6_scan_idx, density_scan_idx, density_index_3, &
                    t6_index_3, density_interp_order, t6_interp_order
               write(short_file_unit,'(3E12.4)') &
                    (x_interp_workspace(x_print_idx,t6_scan_idx,density_scan_idx), &
                    x_print_idx=x_index_2,x_index_2+2)
            end if
            x_interp_result(t6_scan_idx,density_scan_idx) = &
                 x_interp_result(t6_scan_idx,density_scan_idx)*x_interp_weight + &
                 x_interp_result_alt(t6_scan_idx,density_scan_idx)* &
                 (1.0d0 - x_interp_weight)
            recompute_flag = 1
         end do
   47 continue


      end if

      recompute_flag = 0
!
! ..... completed X interpolation. Now interpolate T6 and rho on a
!       4x4 grid. (t6a(i),i=i1,i1+3),rho(j),j=j1,j1+3)).Procedure
!       mixes overlapping quadratics to obtain smoothed derivatives.
!
!
      call t6rinteos01(density_value, t6_value)
      eos_output(eos_var_idx) = esact
  124 continue

      pressure_scale = t6_temperature*density
      eos_output(eos_index_inverse(1)) = eos_output(eos_index_inverse(1))* &
           pressure_scale   ! interpolated in p/po
      eos_output(eos_index_inverse(2)) = eos_output(eos_index_inverse(2))* &
           t6_temperature   ! interpolated in E/T6
      mean_molecular_weight = gmass01(hydrogen_fraction, table_metal_fraction, &
           total_moles, ground_state_energy, metal_mole_fraction, &
           species_mass_fraction)
      if (rad_flag.eq.1) then
         call radsub01(t6_temperature, density, total_moles, &
              mean_molecular_weight)
      else
         eos_output(eos_index_inverse(4)) = eos_output(eos_index_inverse(4))* &
              total_moles*molar_gas_constant_mbcc/mean_molecular_weight
      end if
      return

   61 write(short_file_unit,*) routine_id, "Mass fractions exceed unity (61)"
      write(short_file_unit,*) 'Z, XH', table_metal_fraction, hydrogen_fraction
      stop
   62 continue
      write(short_file_unit,*) routine_id, " T6 or rho outside of table range (62)"
      write(short_file_unit,*) "t6, t6a(1),t6a(nt):", t6_value, t6_grid(1), &
           t6_grid(nt)
      write(short_file_unit,*) "slr,r,rho(1),rho(nr):", density_value, &
           density, density_grid(1), density_grid(nr)
      return 1
   65 continue
      write(short_file_unit,*) routine_id, "T6/rho in empty region of OPAL 2001 EOS", &
           " table (65)"
      write(short_file_unit,'("xh,t6,r=", 3E12.4)') hydrogen_fraction, &
           t6_temperature, density
      return 1
   66 write(short_file_unit,*) routine_id, " Z does not match Z in OPAL 2001 EOS files", &
           " you are using (66)"
      write(short_file_unit,'("mf,zz(mf)=",I5,E12.4)') x_index_lo, &
           z_table(x_index_lo)
      write(short_file_unit,'("  iq,ip,k3,l3,xh,t6,r,z= ",4I5,4E12.4)') &
           t6_interp_order, density_interp_order, t6_index_3, density_index_3, &
           hydrogen_fraction, t6_temperature, density, table_metal_fraction
      stop

end subroutine esac01
