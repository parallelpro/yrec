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
!   opal_eos%eos_output_06(1) pressure, in megabars (1e12 dyne/cm**2)
!   opal_eos%eos_output_06(2) energy, in 1e12 erg/gm (zero at T6=0)
!   opal_eos%eos_output_06(3) entropy, in units of energy/T6
!   opal_eos%eos_output_06(4) dE/dRho at constant T6
!   opal_eos%eos_output_06(5) specific heat, dE/dT6 at constant volume
!   opal_eos%eos_output_06(6) dlogP/dlogRho at constant T6 (Cox & Guili eq 9.82)
!   opal_eos%eos_output_06(7) dlogP/dlogT6 at constant Rho (Cox & Guili eq 9.81)
!   opal_eos%eos_output_06(8) gamma1 (Cox & Guili eq 9.88)
!   opal_eos%eos_output_06(9) gamma2/(gamma2-1) (Cox & Guili eq 9.88)
subroutine esac06(hydrogen_fraction, t6_temperature, density, &
     deriv_order, rad_flag, ierr, *)

      use opal_eos_lib
      use luout_lib
      implicit none

      double precision, intent(in) :: hydrogen_fraction, t6_temperature, &
           density
      integer, intent(in) :: deriv_order, rad_flag

      integer, parameter :: mx = 5, mv = 10, nr = 169, nt = 197









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
! x_grid_06/eos_var_order_06/t6_index_lo_06/density_index_edge_06
! defaults moved to opal_eos_lib.f90: DATA can no longer target them
! here now that they're use-associated.
! --- locals ---
      integer :: x_loop_index_06
      double precision :: hydrogen_fraction_copy, density_copy   ! xxi/ri: assigned but not used further
      double precision :: t6_value, density_value   ! working copies (slt/slr)
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

      integer, intent(out) :: ierr

      ierr = 0

      blank_line = ' '
      if (deriv_order.gt.9) then
         write (short_file_unit,'(A, " iorder cannot exceed 9")') routine_id
      end if
      if ((rad_flag.ne.0) .and. (rad_flag.ne.1)) then
         write (short_file_unit,'(A, " Irad must be 0 or 1")') routine_id
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      hydrogen_fraction_copy = hydrogen_fraction
      density_copy = density
!
      t6_value = t6_temperature
      density_value = density
!
      if (opal_eos%table_loaded_flag.ne.12345678) then
         opal_eos%table_loaded_flag = 12345678
         do index_idx = 1, 10
            do species_idx = 1, 10
               if (opal_eos%eos_var_order_06(index_idx).eq.species_idx) &
                    opal_eos%eos_index_inverse_06(index_idx) = species_idx
            end do
         end do
         do species_idx = 1, mx
            opal_eos%x_grid_copy_06(species_idx) = opal_eos%x_grid_06(species_idx)
         end do
!
! ..... read the data files
         call readcoeos06(ierr)
         if (ierr /= 0) return
         opal_eos%table_metal_fraction_06 = opal_eos%z_table_06(1)

         if (opal_eos%table_metal_fraction_06+hydrogen_fraction-1.0d-6.gt.1) then
            write(short_file_unit,*) routine_id, "Mass fractions exceed unity (61)"
            write(short_file_unit,*) 'Z, XH', opal_eos%table_metal_fraction_06, hydrogen_fraction
            
            
            ierr = 1
            return
         end if
      end if
!
!
! ..... Determine T6,rho grid points to use in the
!       interpolation.
      if ((t6_value.gt.opal_eos%t6_grid_06(1)) .or. (t6_value.lt.opal_eos%t6_grid_06(nt))) then
         write(short_file_unit,*) routine_id, " T6/LogR outside of table range (62)"
         write(short_file_unit,*) "slt, t6a(1),t6a(nt):", t6_value, opal_eos%t6_grid_06(1), &
         opal_eos%t6_grid_06(nt)
         write(short_file_unit,*) "slr, rho(1), rho(nr):", density_value, &
         opal_eos%density_grid_06(1), opal_eos%density_grid_06(nr)
         return 1
      end if
      if ((density_value.lt.opal_eos%density_grid_06(1)) .or. (density_value.gt.opal_eos%density_grid_06(nr))) then
         write(short_file_unit,*) routine_id, " T6/LogR outside of table range (62)"
         write(short_file_unit,*) "slt, t6a(1),t6a(nt):", t6_value, opal_eos%t6_grid_06(1), &
         opal_eos%t6_grid_06(nt)
         write(short_file_unit,*) "slr, rho(1), rho(nr):", density_value, &
         opal_eos%density_grid_06(1), opal_eos%density_grid_06(nr)
         return 1
      end if
!
!
!
      lo_idx = 2
      hi_idx = mx
    8 if (hi_idx-lo_idx.gt.1) then
         mid_idx = (hi_idx+lo_idx)/2
         if (hydrogen_fraction.le.opal_eos%x_grid_06(mid_idx)+1.0d-7) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
         go to 8
      end if
      result_idx = hi_idx
      opal_eos%x_index_lo_06 = result_idx - 2
      x_index_2 = result_idx - 1
      x_index_3 = result_idx
      x_index_4 = result_idx + 1
      x_index_hi = x_index_4
      if (hydrogen_fraction.lt.1.0d-6) then
         opal_eos%x_index_lo_06 = 1
         x_index_2 = 1
         x_index_3 = 1
         x_index_4 = 2
         x_index_hi = 1
      end if
      if ((hydrogen_fraction.le.opal_eos%x_grid_06(2)+1.0d-7) .or. &
           (hydrogen_fraction.ge.opal_eos%x_grid_06(mx-2)-1.0d-7)) x_index_hi = x_index_3
!
      lo_idx = 2
      hi_idx = nr
   12 if (hi_idx-lo_idx.gt.1) then
         mid_idx = (hi_idx+lo_idx)/2
         if (density_value.eq.opal_eos%density_grid_06(mid_idx)) then
            hi_idx = mid_idx
            go to 13
         end if
         if (density_value.le.opal_eos%density_grid_06(mid_idx)) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
         go to 12
      end if
   13 result_idx = hi_idx
      opal_eos%density_index_1_06 = result_idx - 2
      opal_eos%density_index_2_06 = result_idx - 1
      opal_eos%density_index_3_06 = result_idx
      opal_eos%density_index_4_06 = opal_eos%density_index_3_06 + 1
      density_order_hi = 3
      if (opal_eos%density_index_4_06.gt.nr) density_order_hi = 2
!
      lo_idx = nt
      hi_idx = 2
   11 if (lo_idx-hi_idx.gt.1) then
         mid_idx = (hi_idx+lo_idx)/2
         if (t6_temperature.eq.opal_eos%t6_list_06(1,mid_idx)) then
            lo_idx = mid_idx
            go to 14
         end if
         if (t6_temperature.le.opal_eos%t6_list_06(1,mid_idx)) then
            hi_idx = mid_idx
         else
            lo_idx = mid_idx
         end if
         go to 11
      end if
   14 result_idx = lo_idx
      opal_eos%t6_index_1_06 = result_idx - 2
      opal_eos%t6_index_2_06 = result_idx - 1
      opal_eos%t6_index_3_06 = result_idx
      opal_eos%t6_index_4_06 = opal_eos%t6_index_3_06 + 1
      t6_order_hi = 3
      if (opal_eos%t6_index_4_06.gt.nt) t6_order_hi = 2
      if (opal_eos%t6_index_3_06.eq.0) then
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
      do x_loop_index_06 = opal_eos%x_index_lo_06, x_index_hi
         do density_scan_idx = opal_eos%density_index_1_06, opal_eos%density_index_1_06+1
            do t6_scan_idx = opal_eos%t6_index_1_06, opal_eos%t6_index_1_06+1
               table_sum_2x2 = table_sum_2x2 + &
                    opal_eos%eos_table_06(x_loop_index_06,1,t6_scan_idx,density_scan_idx)
            end do
         end do
         do density_scan_idx = opal_eos%density_index_1_06, opal_eos%density_index_1_06+2
            do t6_scan_idx = opal_eos%t6_index_1_06, opal_eos%t6_index_1_06+2
               table_sum_3x3 = table_sum_3x3 + &
                    opal_eos%eos_table_06(x_loop_index_06,1,t6_scan_idx,density_scan_idx)
            end do
         end do
         if (t6_order_hi.eq.3) then
            do density_scan_idx = opal_eos%density_index_1_06, opal_eos%density_index_1_06+2
               do t6_scan_idx = opal_eos%t6_index_1_06, opal_eos%t6_index_1_06+t6_order_hi
                  table_sum_3x4 = table_sum_3x4 + &
                       opal_eos%eos_table_06(x_loop_index_06,1,t6_scan_idx,density_scan_idx)
               end do
            end do
         else
            table_sum_3x4 = 2.0d+30
         end if
         if (density_order_hi.eq.3) then
            do density_scan_idx = opal_eos%density_index_1_06, opal_eos%density_index_1_06+3
               do t6_scan_idx = opal_eos%t6_index_1_06, opal_eos%t6_index_1_06+t6_order_hi
                  table_sum_4x4 = table_sum_4x4 + &
                       opal_eos%eos_table_06(x_loop_index_06,1,t6_scan_idx,density_scan_idx)
               end do
            end do
         else
            table_sum_4x4 = 2.0d+30
         end if
      end do
      opal_eos%density_interp_order_06 = 2
      opal_eos%t6_interp_order_06 = 2
      if (table_sum_3x3.gt.1.0d+30) then
         if (table_sum_2x2.lt.1.0d+25) then
            opal_eos%t6_index_1_06 = opal_eos%t6_index_3_06 - 3
            opal_eos%t6_index_2_06 = opal_eos%t6_index_1_06 + 1
            opal_eos%t6_index_3_06 = opal_eos%t6_index_2_06 + 1
            opal_eos%density_index_1_06 = opal_eos%density_index_3_06 - 3
            opal_eos%density_index_2_06 = opal_eos%density_index_1_06 + 1
            opal_eos%density_index_3_06 = opal_eos%density_index_2_06 + 1
            go to 15
         else
            write(short_file_unit,*) routine_id, "T6/log rho in empty region of table (65)"
            write(short_file_unit,'("xh,t6,r=", 3E12.4)') hydrogen_fraction, &
            t6_temperature, density
            return 1
         end if
      end if
      if (table_sum_3x4.lt.1.0d+30) opal_eos%t6_interp_order_06 = 3
      if (table_sum_4x4.lt.1.0d+30) opal_eos%density_interp_order_06 = 3

      if (t6_temperature.ge.opal_eos%t6_list_06(1,2)+1.0d-7) opal_eos%t6_interp_order_06 = 2
      if (density_value.le.opal_eos%density_grid_06(2)+1.0d-15) opal_eos%density_interp_order_06 = 2

      if ((opal_eos%density_index_3_06.eq.nr) .or. (opal_eos%t6_index_3_06.eq.nt)) then
         opal_eos%density_interp_order_06 = 2
         opal_eos%t6_interp_order_06 = 2
      end if

   15 continue
      do eos_var_idx = 1, deriv_order
      do x_loop_index_06 = opal_eos%x_index_lo_06, x_index_hi

      recompute_flag = 0
!__________
      do density_scan_idx = opal_eos%density_index_1_06, opal_eos%density_index_1_06+opal_eos%density_interp_order_06
         do t6_scan_idx = opal_eos%t6_index_1_06, opal_eos%t6_index_1_06+opal_eos%t6_interp_order_06
            opal_eos%x_interp_workspace_06(x_loop_index_06,t6_scan_idx,density_scan_idx) = &
                 opal_eos%eos_table_06(x_loop_index_06,eos_var_idx,t6_scan_idx,density_scan_idx)
            recompute_flag = 1
         end do
      end do
  123 continue
      end do
      if ((opal_eos%z_table_06(x_index_2).ne.opal_eos%z_table_06(opal_eos%x_index_lo_06)) .or. &
           (opal_eos%z_table_06(x_index_3).ne.opal_eos%z_table_06(opal_eos%x_index_lo_06))) then
         write(short_file_unit,'(A,"Z does not match Z in OEOS06 files you are" &
              &," using")') routine_id
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      if (opal_eos%table_metal_fraction_06.ne.opal_eos%z_table_06(opal_eos%x_index_lo_06)) go to 66
      recompute_flag = 0
      cache_slot = 1
      do density_scan_idx = opal_eos%density_index_1_06, opal_eos%density_index_1_06+opal_eos%density_interp_order_06
         do t6_scan_idx = opal_eos%t6_index_1_06, opal_eos%t6_index_1_06+opal_eos%t6_interp_order_06
            if (x_index_hi.eq.1) then
               opal_eos%x_interp_result_06(t6_scan_idx,density_scan_idx) = &
                    opal_eos%x_interp_workspace_06(opal_eos%x_index_lo_06,t6_scan_idx,density_scan_idx)
               cycle
            end if
            opal_eos%x_interp_result_06(t6_scan_idx,density_scan_idx) = &
                 quadeos06(recompute_flag, cache_slot, hydrogen_fraction, &
                 opal_eos%x_interp_workspace_06(opal_eos%x_index_lo_06,t6_scan_idx,density_scan_idx), &
                 opal_eos%x_interp_workspace_06(x_index_2,t6_scan_idx,density_scan_idx), &
                 opal_eos%x_interp_workspace_06(x_index_3,t6_scan_idx,density_scan_idx), &
                 opal_eos%x_grid_copy_06(opal_eos%x_index_lo_06), opal_eos%x_grid_copy_06(x_index_2), &
                 opal_eos%x_grid_copy_06(x_index_3))
            if (opal_eos%x_interp_result_06(t6_scan_idx,density_scan_idx).gt.1.0d+20) then
               write(short_file_unit,'(A," problem it,ir,l3,k3,iq,ip=", 6I5)') &
                    routine_id, t6_scan_idx, density_scan_idx, opal_eos%density_index_3_06, &
                    opal_eos%t6_index_3_06, opal_eos%density_interp_order_06, opal_eos%t6_interp_order_06
               write(short_file_unit,'(3E12.4)') &
                    (opal_eos%x_interp_workspace_06(x_print_idx,t6_scan_idx,density_scan_idx), &
                    x_print_idx=opal_eos%x_index_lo_06,opal_eos%x_index_lo_06+2)
            end if
            recompute_flag = 1
   46       continue
         end do
   45 continue
      end do

      if (x_index_4.eq.x_index_hi) then  ! interpolate between quadratics
      recompute_flag = 0
      cache_slot = 1
      x_interp_weight = (opal_eos%x_grid_copy_06(x_index_3) - hydrogen_fraction)* &
           opal_eos%x_grid_spacing_inv_06(x_index_3)
      do density_scan_idx = opal_eos%density_index_1_06, opal_eos%density_index_1_06+opal_eos%density_interp_order_06
         do t6_scan_idx = opal_eos%t6_index_1_06, opal_eos%t6_index_1_06+opal_eos%t6_interp_order_06
            opal_eos%x_interp_result_alt_06(t6_scan_idx,density_scan_idx) = &
                 quadeos06(recompute_flag, cache_slot, hydrogen_fraction, &
                 opal_eos%x_interp_workspace_06(x_index_2,t6_scan_idx,density_scan_idx), &
                 opal_eos%x_interp_workspace_06(x_index_3,t6_scan_idx,density_scan_idx), &
                 opal_eos%x_interp_workspace_06(x_index_4,t6_scan_idx,density_scan_idx), &
                 opal_eos%x_grid_copy_06(x_index_2), opal_eos%x_grid_copy_06(x_index_3), &
                 opal_eos%x_grid_copy_06(x_index_4))
            if (opal_eos%x_interp_result_06(t6_scan_idx,density_scan_idx).gt.1.0d+20) then
               write(short_file_unit,'(A," problem it,ir,l3,k3,iq,ip=", 6I5)') &
                    routine_id, t6_scan_idx, density_scan_idx, opal_eos%density_index_3_06, &
                    opal_eos%t6_index_3_06, opal_eos%density_interp_order_06, opal_eos%t6_interp_order_06
               write(short_file_unit,'(3E12.4)') &
                    (opal_eos%x_interp_workspace_06(x_print_idx,t6_scan_idx,density_scan_idx), &
                    x_print_idx=x_index_2,x_index_2+2)
            end if
            opal_eos%x_interp_result_06(t6_scan_idx,density_scan_idx) = &
                 opal_eos%x_interp_result_06(t6_scan_idx,density_scan_idx)*x_interp_weight + &
                 opal_eos%x_interp_result_alt_06(t6_scan_idx,density_scan_idx)* &
                 (1.0d0 - x_interp_weight)
            recompute_flag = 1
         end do
   47 continue
      end do


      end if

      recompute_flag = 0
!
! ..... completed X interpolation. Now interpolate T6 and rho on a
!       4x4 grid. (t6a(i),i=i1,i1+3),rho(j),j=j1,j1+3)).Procedure
!       mixes overlapping quadratics to obtain smoothed derivatives.
!
!
      call t6rinteos06(density_value, t6_value, ierr)
      if (ierr /= 0) return
      opal_eos%eos_output_06(eos_var_idx) = opal_eos%esact_06
  124 continue
      end do
      pressure_scale = t6_temperature*density
      opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(1)) = opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(1))* &
           pressure_scale   ! interpolated in p/po
      opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(2)) = opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(2))* &
           t6_temperature   ! interpolated in E/T6
      mean_molecular_weight = gmass06(hydrogen_fraction, opal_eos%table_metal_fraction_06, &
           total_moles, ground_state_energy, metal_mole_fraction, &
           species_mass_fraction)
      if (rad_flag.eq.1) then
         call radsub06(rad_flag, t6_temperature, density, total_moles, &
              mean_molecular_weight)
      else
         opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(5)) = opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(5))* &
              total_moles*molar_gas_constant_mbcc/mean_molecular_weight
      end if
      return

   61 write(short_file_unit,*) routine_id, "Mass fractions exceed unity (61)"
      write(short_file_unit,*) 'Z, XH', opal_eos%table_metal_fraction_06, hydrogen_fraction
      ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
      ! facades stop when their caller passes no ierr.
      ierr = 1
      return
   62 write(short_file_unit,*) routine_id, " T6/LogR outside of table range (62)"
      write(short_file_unit,*) "slt, t6a(1),t6a(nt):", t6_value, opal_eos%t6_grid_06(1), &
           opal_eos%t6_grid_06(nt)
      write(short_file_unit,*) "slr, rho(1), rho(nr):", density_value, &
           opal_eos%density_grid_06(1), opal_eos%density_grid_06(nr)
      return 1
   65 write(short_file_unit,*) routine_id, "T6/log rho in empty region of table (65)"
      write(short_file_unit,'("xh,t6,r=", 3E12.4)') hydrogen_fraction, &
           t6_temperature, density
      return 1
   66 write(short_file_unit,*) routine_id, " Z does not match Z in EOSdata* files you ", &
           "are using (66)"
      write(short_file_unit,'("mf,zz(mf)=",I5,E12.4)') opal_eos%x_index_lo_06, &
           opal_eos%z_table_06(opal_eos%x_index_lo_06)
      write(short_file_unit,'("  iq,ip,k3,l3,xh,t6,r,z= ",4I5,4E12.4)') &
           opal_eos%t6_interp_order_06, opal_eos%density_interp_order_06, opal_eos%t6_index_3_06, opal_eos%density_index_3_06, &
           hydrogen_fraction, t6_temperature, density, opal_eos%table_metal_fraction_06
      ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
      ! facades stop when their caller passes no ierr.
      ierr = 1
      return

end subroutine esac06
