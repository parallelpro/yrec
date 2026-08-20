!----------------------------------------------------------------------
! esac06
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original esac06.f; only variable names, source form, and comment
! style were updated.
!
! OPAL 2006 EOS analogue of esac01.f90 (see there and esac.f90 for
! the general description of the interpolation scheme); the
! structure is otherwise identical to esac01.f90 (same iorder cap of
! 9, same iqu/ipu pre-checks, same l3==nr .or. k3==nt boundary reset).
!
!   eos_output(1) pressure, in megabars (1e12 dyne/cm**2)
!   eos_output(2) energy, in 1e12 erg/gm (zero at T6=0)
!   eos_output(3) entropy, in units of energy/T6
!   eos_output(4) dE/dRho at constant T6
!   eos_output(5) specific heat, dE/dT6 at constant volume
!   eos_output(6) dlogP/dlogRho at constant T6 (Cox & Guili eq 9.82)
!   eos_output(7) dlogP/dlogT6 at constant Rho (Cox & Guili eq 9.81)
!   eos_output(8) gamma1 (Cox & Guili eq 9.88)
!   eos_output(9) gamma2/(gamma2-1) (Cox & Guili eq 9.88)
subroutine esac06(hydrogen_fraction, t6_temperature, density, &
     deriv_order, rad_flag, *)

      implicit none

      double precision, intent(in) :: hydrogen_fraction, t6_temperature, &
           density
      integer, intent(in) :: deriv_order, rad_flag

      integer, parameter :: mx = 5, mv = 10, nr = 169, nt = 197

! common/lreadco/: table_loaded_flag is shared (by COMMON block name)
! with the OPAL 1995/2001 readers -- see esac.f90's header.
      integer :: table_loaded_flag
      common/lreadco/ table_loaded_flag

! common/eeeos06/: both members used here (see esac01.f90's analogous
! common/eeeos/ for the description).
      double precision :: x_interp_workspace(mx,nt,nr), x_grid_copy(mx)
      common/eeeos06/ x_interp_workspace, x_grid_copy

! common/aaeos06/: not used directly here (consumed by t6rinteos06.f90).
      double precision :: rho_interp_hi(4), rho_interp_lo(4), xxh
      common/aaeos06/ rho_interp_hi, rho_interp_lo, xxh

! common/aeos06/: the main OPAL 2006 EOS table and its interpolation
! grids/scratch arrays. See esac01.f90's analogous common/aeos/ for
! the full member-by-member description.
      double precision :: eos_table(mx,mv,nt,nr), t6_list(nr,nt), &
           density_grid(nr), t6_grid(nt), x_interp_result(nt,nr), &
           x_interp_result_alt(nt,nr), x_grid_spacing_inv(mx), &
           t6_grid_spacing_inv(nt), density_grid_spacing_inv(nr)
      integer :: x_loop_index, x_index_lo
      double precision :: x_grid(mx)
      common/aeos06/ eos_table, t6_list, density_grid, t6_grid, &
           x_interp_result, x_interp_result_alt, x_grid_spacing_inv, &
           t6_grid_spacing_inv, density_grid_spacing_inv, x_loop_index, &
           x_index_lo, x_grid

! common/beos06/: z_table, eos_var_order, eos_index_inverse, and
! t6_index_lo as in common/beos/ (esac01.f90); density_index_edge
! (original nra) is new here -- see rhoofp06.f90.
      double precision :: z_table(mx)
      integer :: eos_index_inverse(10), eos_var_order(10), &
           t6_index_lo(nr), density_index_edge(nt)
      common/beos06/ z_table, eos_index_inverse, eos_var_order, &
           t6_index_lo, density_index_edge

! common/bbeos06/: density index window (l1..l4), t6 index window
! (k1..k4), and the interpolation-order flags t6_interp_order
! (original ip) and density_interp_order (original iq).
      integer :: density_index_1, density_index_2, density_index_3, &
           density_index_4, t6_index_1, t6_index_2, t6_index_3, &
           t6_index_4, t6_interp_order, density_interp_order
      common/bbeos06/ density_index_1, density_index_2, density_index_3, &
           density_index_4, t6_index_1, t6_index_2, t6_index_3, &
           t6_index_4, t6_interp_order, density_interp_order

! common/eeos06/: eos_output is the caller-facing interpolated result.
      double precision :: esact, eos_output(mv)
      common/eeos06/ esact, eos_output

! common/luout/: only short_file_unit is used here.
      integer :: ilast, idebug, itrack, short_file_unit, imilne, imodpt, &
           istor, main_output_unit
      common/luout/ ilast, idebug, itrack, short_file_unit, imilne, &
           imodpt, istor, main_output_unit

      double precision :: species_mass_fraction(7)
      double precision :: molar_gas_constant_mbcc
      character(len=1) :: blank_line
      character(len=15) :: routine_id

      integer :: fill_idx

! NOTE: this literal has no D0 suffix in the original (esac06.f) --
! preserved verbatim: it is parsed as a single-precision constant and
! then widened to double precision, which is NOT bit-identical to the
! correctly-rounded double value.
      data molar_gas_constant_mbcc/83.14511/
      data routine_id/"OPALEOS/ESAC06:"/
! NOTE: no D-suffix in the original (esac06.f) -- preserved verbatim:
! these are parsed as single-precision constants and then widened to
! double precision (0.2/0.4/0.6/0.8 are not exactly representable in
! binary floating point, so this is NOT bit-identical to the
! correctly-rounded double values).
      data (x_grid(fill_idx), fill_idx=1,mx) /0.0, 0.2, 0.4, 0.6, 0.8/
      data (eos_var_order(fill_idx), fill_idx=1,10) /1,2,3,4,5,6,7,8,9,10/
      data (t6_index_lo(fill_idx), fill_idx=1,nr) &
           /87*197, 7*191, 190, 2*189, 185, 179, 170, 2*149, 133, 125, 123, &
           122, 120, 115, 113, 107, 102, 2*80, 72, 68, 66, 64, 62, 56, 54, &
           52, 51, 2*50, 49, 47, 2*45, 43, 42, 28*40, 39, 37, 36, 35, 34, &
           32, 31, 30, 29, 27, 26/
      data (density_index_edge(fill_idx), fill_idx=1,nt) &
           /26*169, 168, 2*167, 166, 165, 164, 2*163, 162, 161, 160, 2*159, &
           158, 2*130, 129, 2*128, 2*126, 2*125, 124, 122, 121, 2*120, &
           2*119, 6*118, 2*117, 2*116, 2*115, 4*114, 8*113, 22*111, 5*110, &
           6*109, 2*108, 5*107, 2*106, 105, 2*104, 8*103, 16*102, 21*100, &
           9*99, 6*98, 4*97, 95, 94, 6*87/

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
      double precision, external :: quadeos06, gmass06

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
         call readcoeos06
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
         write (short_file_unit,'(A, " ihi,ilo,imd",3I5)') routine_id, &
              hi_idx, lo_idx, mid_idx
      end if
!     endif

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
         write(short_file_unit,'(A,"Z does not match Z in OEOS06 files you are" &
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
                 quadeos06(recompute_flag, cache_slot, hydrogen_fraction, &
                 x_interp_workspace(x_index_lo,t6_scan_idx,density_scan_idx), &
                 x_interp_workspace(x_index_2,t6_scan_idx,density_scan_idx), &
                 x_interp_workspace(x_index_3,t6_scan_idx,density_scan_idx), &
                 x_grid_copy(x_index_lo), x_grid_copy(x_index_2), &
                 x_grid_copy(x_index_3))
            if (x_interp_result(t6_scan_idx,density_scan_idx).gt.1.0d+20) then
               write(short_file_unit,'(A," problem it,ir,l3,k3,iq,ip=", 6I5)') &
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
                 quadeos06(recompute_flag, cache_slot, hydrogen_fraction, &
                 x_interp_workspace(x_index_2,t6_scan_idx,density_scan_idx), &
                 x_interp_workspace(x_index_3,t6_scan_idx,density_scan_idx), &
                 x_interp_workspace(x_index_4,t6_scan_idx,density_scan_idx), &
                 x_grid_copy(x_index_2), x_grid_copy(x_index_3), &
                 x_grid_copy(x_index_4))
            if (x_interp_result(t6_scan_idx,density_scan_idx).gt.1.0d+20) then
               write(short_file_unit,'(A," problem it,ir,l3,k3,iq,ip=", 6I5)') &
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
      call t6rinteos06(density_value, t6_value)
      eos_output(eos_var_idx) = esact
  124 continue
      pressure_scale = t6_temperature*density
      eos_output(eos_index_inverse(1)) = eos_output(eos_index_inverse(1))* &
           pressure_scale   ! interpolated in p/po
      eos_output(eos_index_inverse(2)) = eos_output(eos_index_inverse(2))* &
           t6_temperature   ! interpolated in E/T6
      mean_molecular_weight = gmass06(hydrogen_fraction, table_metal_fraction, &
           total_moles, ground_state_energy, metal_mole_fraction, &
           species_mass_fraction)
      if (rad_flag.eq.1) then
         call radsub06(rad_flag, t6_temperature, density, total_moles, &
              mean_molecular_weight)
      else
         eos_output(eos_index_inverse(5)) = eos_output(eos_index_inverse(5))* &
              total_moles*molar_gas_constant_mbcc/mean_molecular_weight
      end if
      return

   61 write(short_file_unit,*) routine_id, "Mass fractions exceed unity (61)"
      write(short_file_unit,*) 'Z, XH', table_metal_fraction, hydrogen_fraction
      stop
   62 write(short_file_unit,*) routine_id, " T6/LogR outside of table range (62)"
      write(short_file_unit,*) "slt, t6a(1),t6a(nt):", t6_value, t6_grid(1), &
           t6_grid(nt)
      write(short_file_unit,*) "slr, rho(1), rho(nr):", density_value, &
           density_grid(1), density_grid(nr)
      return 1
   65 write(short_file_unit,*) routine_id, "T6/log rho in empty region of table (65)"
      write(short_file_unit,'("xh,t6,r=", 3E12.4)') hydrogen_fraction, &
           t6_temperature, density
      return 1
   66 write(short_file_unit,*) routine_id, " Z does not match Z in EOSdata* files you ", &
           "are using (66)"
      write(short_file_unit,'("mf,zz(mf)=",I5,E12.4)') x_index_lo, &
           z_table(x_index_lo)
      write(short_file_unit,'("  iq,ip,k3,l3,xh,t6,r,z= ",4I5,4E12.4)') &
           t6_interp_order, density_interp_order, t6_index_3, density_index_3, &
           hydrogen_fraction, t6_temperature, density, table_metal_fraction
      stop

end subroutine esac06
