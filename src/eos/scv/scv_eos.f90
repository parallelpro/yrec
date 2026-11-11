!----------------------------------------------------------------------
! scv_eos
!----------------------------------------------------------------------
! The Saumon-Chabrier-Van Horn equation of state for cool, dense
! H/He gas, from the pure-H, pure-He and pure-metal tables.
!     build_scv_envelope_table  combine the three tables at the envelope
!                               X and Z (eos_set_mixture, whenever the
!                               envelope composition changes)
!     compute_scv_envelope_eos  density, derivatives and ionization
!                               fractions at (log T, P): from the
!                               envelope table at the envelope
!                               composition, else from compute_scv_eos
!     compute_scv_eos           the same for any (X, Z), from the three
!                               tables directly
! The tables live in scv_eos_tables.
module scv_eos
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use eos_results, only: eos_result
      use scv_eos_tables, only: scv_cell, last_scv_cell, iscv_dlnrho_dlnp, iscv_dlnrho_dlnt, &
           iscv_dlns_dlnt, iscv_frac_atom_ion1, iscv_frac_neutral, &
           iscv_log10_p, iscv_log10_rho, iscv_log10_s, iscvenv_cp, &
           iscvenv_dlnrho_dlnp, iscvenv_dlnrho_dlnt, iscvenv_du_dt, &
           iscvenv_log10_p, iscvenv_log10_rho, iscvz_dlnrho_dlnp, &
           iscvz_dlnrho_dlnt, iscvz_log10_du_dt, iscvz_log10_rho, nptsx, &
           scv_nt, smix, tablenv, tablex, tabley, tablez, &
           tlogx, tablenv_hydrogen_fraction, tablenv_metal_fraction, scv_tables_loaded
      use physical_constants, only: cc13, ln10, radiation_constant_over_3, &
           boltzmann_constant_cgs, hydrogen_atom_mass_cgs, helium_atom_mass_cgs
      use numerics_lib, only: inter3, interp, lagrange4, natural_spline, splint, &
           numerics_status_text
      use eos_messages, only: eos_log, eos_report, message_length, unwritten_line
      use eos_mixture, only: mean_weights_from_density
      use scv_table_reader, only: read_scv_tables
      use eos_configuration, only: eos_settings, eos_table_paths_t
      use eos_table_source, only: eos_table_source_t, stage_in_place_of_saha, always_in_window, &
           whole_weight_inside, adopt_whole_result, blend_table_result
      use math_lib
      implicit none
      private
      public :: compute_scv_envelope_eos, build_scv_envelope_table, scv_table_source

! the lowest log10 P_gas (cgs) the tables cover; how far a zone's X or
! Z may be from the envelope table's for the table to serve it
      real(dp), parameter :: scv_min_log10_gas_pressure = 4.0d0, &
           scv_composition_tolerance = 1.0d-5

contains

!----------------------------------------------------------------------
! scv_table_source
!----------------------------------------------------------------------
! SCV as a table source (eos_table_source): selected by use_scv_eos, its
! three tables read at eos_init, consulted below the Saha cutoff in place
! of the Saha solution (compute_scv_envelope_eos, from a copy of it), and
! taken whole wherever its tables cover the point: it has no edge weight
! yet (M88).
function scv_table_source() result(src)
      type(eos_table_source_t) :: src
      src%name = 'SCV'
      src%stage = stage_in_place_of_saha
      src%selected => scv_selected
      src%read_tables => read_scv_source
      src%in_window => always_in_window
      src%lookup => scv_lookup
      src%edge_weight => whole_weight_inside
      src%adopt => adopt_whole_result
      src%blend => blend_table_result
end function scv_table_source

logical function scv_selected()
      scv_selected = eos_settings%use_scv_eos
end function scv_selected

subroutine read_scv_source(paths, ierr)
      type(eos_table_paths_t), intent(in) :: paths
      integer, intent(out) :: ierr
      call read_scv_tables(trim(paths%scv_h), trim(paths%scv_he), trim(paths%scv_z), ierr)
end subroutine read_scv_source

! the lookup takes P from table_r (the copy of the Saha result), not from
! log10_pressure
subroutine scv_lookup(log10_temperature, log10_pressure, hydrogen_fraction, metal_fraction, &
     table_r, covers, ierr)
      real(dp), intent(in) :: log10_temperature, log10_pressure, hydrogen_fraction, metal_fraction
      type(eos_result), intent(inout) :: table_r
      logical, intent(out) :: covers
      integer, intent(out) :: ierr
      call compute_scv_envelope_eos(log10_temperature, hydrogen_fraction, metal_fraction, &
           table_r, covers, ierr)
end subroutine scv_lookup

!----------------------------------------------------------------------
! compute_scv_envelope_eos (formerly eqscve)
!----------------------------------------------------------------------
! The SCV EOS at (log T, P in r) from the envelope table, the table
! built for the envelope composition: a zone whose X or Z differs from
! the table's by more than scv_composition_tolerance takes
! compute_scv_eos (the three tables, any composition) instead. With the
! cell around (log T, log P_gas) (find_scv_cell; covers false outside
! the table, r untouched), a 4-point natural spline in log P_gas and then
! in log T of the table's five columns, smoothed across the next cell
! within tol_temp_smooth / tol_pressure_smooth of a node; from them the
! density, the mean weights it implies, dln rho/dln T and dln rho/dln P
! with radiation, cp and grad_ad; and the ionization fractions from the
! pure H and He tables (scv_ionization_fractions). A spline failure is
! ierr (with a message).
subroutine compute_scv_envelope_eos(log10_temperature, hydrogen_fraction, &
     metal_fraction, r, covers, ierr)
      real(dp), intent(in) :: log10_temperature
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      type(eos_result), intent(inout) :: r
      logical, intent(out) :: covers
      integer, intent(out) :: ierr
      integer :: jerr_gate

      integer, parameter :: nts = scv_nt

! Replace 4-pt Lagrange interpolate with 4-pt natural spline
      real(dp) :: press_interp_weights(4), temp_interp_weights(4)
      real(dp) :: interp_nodes(4), ytab_work(4), &
           spline_second_deriv(4), interp_result(5), temp_grid(4,5), &
           interp_result_temp_shift(5), interp_result_press_shift(5), &
           interp_result_both_shift(5), interp_result_temp_smoothed(5), &
           interp_result_press_smoothed(5)
      real(dp), parameter :: tol_pressure_smooth = 0.08d0, &
           tol_temp_smooth = 0.032d0

! --- locals ---
      real(dp) :: beta_complement, log10_gas_pressure
      real(dp) :: temp_dist_above, temp_dist_below, &
           temp_smooth_weight
      real(dp) :: press_dist_above, press_dist_below, &
           press_smooth_weight
      real(dp) :: spline_value
      real(dp) :: radiation_pressure
      real(dp) :: change_t1, change_t2
      real(dp) :: dlnrho_dlnp_gas, dlnrho_dlnt_gas, dlnp_dlnt_gas, &
           dlnp_dlnt, du_dt, &
           radiation_energy_density_term
      logical :: temp_needs_smoothing, press_needs_smoothing, t_centred, p_centred
      integer :: temp_smooth_direction, press_smooth_direction
      integer :: i, ii, j, k
! the cell's first T row and first P_gas column (find_scv_cell)
      integer :: idtt, idp
      ierr = 0
      if (abs(hydrogen_fraction-tablenv_hydrogen_fraction) > scv_composition_tolerance &
           .or. abs(metal_fraction-tablenv_metal_fraction) > scv_composition_tolerance) then
         call compute_scv_eos(log10_temperature, hydrogen_fraction, metal_fraction, &
              r, covers)
         return
      end if
! find gas pressure (which is the quantity which is tabulated).
      beta_complement = 1.0d0 - r%beta
      log10_gas_pressure = log10(r%beta*r%pressure)
! the table's cell around the point (idtt, idp), if it covers it
      call find_scv_cell(log10_temperature, log10_gas_pressure, last_scv_cell, covers, &
           t_centred, p_centred)
      if (.not. covers) return
      idtt = last_scv_cell%first_row
      idp = last_scv_cell%first_column
      temp_needs_smoothing = .false.
      if (t_centred) then
         temp_dist_above = log10_temperature - &
              tlogx(idtt+1)
         temp_dist_below = tlogx(idtt+2) - &
              log10_temperature
         if (temp_dist_above < tol_temp_smooth) then
            if (idtt > 1) then
               temp_needs_smoothing = .true.
               temp_smooth_direction = -1
               temp_smooth_weight = (temp_dist_above+tol_temp_smooth)/ &
                    (2.0d0*tol_temp_smooth)
            end if
         else if (temp_dist_below < tol_temp_smooth) then
            if (idtt < nts-3) then
               temp_needs_smoothing = .true.
               temp_smooth_direction = 1
! the weight is 0.5 at the node, falling to 0 at tol_temp_smooth (the
! mirror of the downward branch)
               temp_smooth_weight = 0.5d0*(1.0d0 - temp_dist_below/tol_temp_smooth)
            end if
         end if
      end if
      press_needs_smoothing = .false.
      if (p_centred) then
         press_dist_above = log10_gas_pressure - &
              tablenv(idtt,idp+1,iscvenv_log10_p)
         press_dist_below = tablenv(idtt, &
              idp+2,iscvenv_log10_p) - log10_gas_pressure
         if (press_dist_above < tol_pressure_smooth) then
            if (idp > 1) then
               press_needs_smoothing = .true.
               press_smooth_direction = -1
               press_smooth_weight = (press_dist_above+tol_pressure_smooth)/ &
                    (2.0d0*tol_pressure_smooth)
            end if
         else if (press_dist_below < tol_pressure_smooth) then
            if (idp < nptsx(idtt)-3) then
               press_needs_smoothing = .true.
               press_smooth_direction = 1
! (as in T)
               press_smooth_weight = 0.5d0*(1.0d0 - press_dist_below/tol_pressure_smooth)
            end if
         end if
      end if
! include radiation pressure in the equation of state.
      radiation_pressure = beta_complement*r%pressure
      call spline_cell(0, 0, interp_result)
      if (ierr /= 0) return
! perform interpolation at adjacent temperature table points
      if (temp_needs_smoothing) then
      call spline_cell(temp_smooth_direction, 0, interp_result_temp_shift)
      if (ierr /= 0) return
      if (temp_smooth_direction == -1) then
         do j = 1,5
            interp_result_temp_smoothed(j) = interp_result_temp_shift(j) + &
                 temp_smooth_weight*(interp_result(j) - &
                 interp_result_temp_shift(j))
         end do
      else
         do j = 1,5
            interp_result_temp_smoothed(j) = interp_result(j) + &
                 temp_smooth_weight*(interp_result_temp_shift(j) - &
                 interp_result(j))
         end do
      end if
      end if
! perform interpolation at adjacent pressure table points
      if (press_needs_smoothing) then
      call spline_cell(0, press_smooth_direction, interp_result_press_shift)
      if (ierr /= 0) return
      if (press_smooth_direction == -1) then
         do j = 1,5
            interp_result_press_smoothed(j) = interp_result_press_shift(j) + &
                 press_smooth_weight*(interp_result(j) - &
                 interp_result_press_shift(j))
         end do
      else
         do j = 1,5
            interp_result_press_smoothed(j) = interp_result(j) + &
                 press_smooth_weight*(interp_result_press_shift(j) - &
                 interp_result(j))
         end do
      end if
      end if
! perform interpolation at adjacent t+p table points
      if (press_needs_smoothing .and. temp_needs_smoothing) then
      call spline_cell(temp_smooth_direction, press_smooth_direction, &
           interp_result_both_shift)
      if (ierr /= 0) return
      end if
      if (temp_needs_smoothing) then
! add changes for both t and p interpolation
         if (press_needs_smoothing) then
            do j = 1,5
               if (temp_smooth_direction == -1) then
! interpolate in t at fixed p
                  change_t1 = interp_result_temp_shift(j) + &
                       temp_smooth_weight*(interp_result(j) - &
                       interp_result_temp_shift(j))
! interpolate in t at different p
                  change_t2 = interp_result_both_shift(j) + &
                       temp_smooth_weight*(interp_result_press_shift(j) - &
                       interp_result_both_shift(j))
               else
                  change_t1 = interp_result(j) + &
                       temp_smooth_weight*(interp_result_temp_shift(j) - &
                       interp_result(j))
                  change_t2 = interp_result_press_shift(j) + &
                       temp_smooth_weight*(interp_result_both_shift(j) - &
                       interp_result_press_shift(j))
               end if
               if (press_smooth_direction == -1) then
                  interp_result(j) = change_t2+press_smooth_weight* &
                       (change_t1-change_t2)
               else
                  interp_result(j) = change_t1+press_smooth_weight* &
                       (change_t2-change_t1)
               end if
            end do
          else
! add t interpolation changes only
             do j = 1,5
                interp_result(j) = interp_result_temp_smoothed(j)
             end do
          end if
       else if (press_needs_smoothing) then
! add p interpolation changes only
          do j = 1,5
             interp_result(j) = interp_result_press_smoothed(j)
          end do
       end if
!  density, and the mean weights it implies (R/mu = beta P/(rho T))
      r%log10_density = interp_result(1)
      r%density = exp(ln10*r%log10_density)
      call mean_weights_from_density(r)
! d ln rho/ d ln p = qdp
! for gas pressure qdp(tot) = qdp(x)*x*rho/rho(x)+qdp(y)*y*rho/rho(y)
! for radiation pressure qdp = 0, so qdp(tot) = qdp(gas)*p/pgas
      dlnrho_dlnp_gas = interp_result(3)
      r%dlnrho_dlnp = dlnrho_dlnp_gas/r%beta
! d ln rho/ d ln t = qdt (note : d ln p/ d ln t = qpt)
! for gas pressure, correct as per qdp
! for radiation pressure, use qdt = -qdp*qpt.  correct qpt for
! radiation pressure and use the corrected qdp, qpt to get qdt.
      dlnrho_dlnt_gas = interp_result(2)
      dlnp_dlnt_gas = -dlnrho_dlnt_gas/dlnrho_dlnp_gas
      dlnp_dlnt = dlnp_dlnt_gas*r%beta + 4.0d0*beta_complement
      r%dlnrho_dlnt = -dlnp_dlnt*r%dlnrho_dlnp
!  cp = du/dt + p*(d ln rho/d ln t)**2/rho/t/(d ln rho/ d ln p), with
!  du/dt from the table (its gas-only cp, interp_result(4), is not
!  needed) corrected for radiation pressure, as are the other terms.
      du_dt = interp_result(5)
! correct du/dt for radiation
      radiation_energy_density_term = 3.0d0*radiation_pressure/r%density
      du_dt = du_dt + 4.0d0*radiation_energy_density_term/r%temperature
! correct cp for radiation pressure
      r%specific_heat_cp = du_dt - r%pressure*r%dlnrho_dlnt*dlnp_dlnt/r%density/ &
           r%temperature
! adiabatic temperature gradient
      r%adiabatic_gradient = -r%pressure*r%dlnrho_dlnt/r%density/r%temperature/ &
           r%specific_heat_cp

! the ionization fractions, from the pure H and He tables at the cell
      do k = 1,4
         interp_nodes(k) = tlogx(k+idtt-1)
      end do
      call interp(interp_nodes, temp_interp_weights, x_eval=log10_temperature)
      do k = 1,4
         interp_nodes(k) = tablenv(idtt, &
              idp+k-1,iscvenv_log10_p)
      end do
      call interp(interp_nodes, press_interp_weights, x_eval=log10_gas_pressure)
      call scv_ionization_fractions(idtt, idp, temp_interp_weights, press_interp_weights, &
           r%ion_fraction)

contains

! Bicubic-spline interpolation of the five envelope-table variables
! (columns 2-6 of tablenv) at (log10_gas_pressure, log10_temperature)
! on the 4x4 cell whose corner is (idtt + t_shift, idp + p_shift):
! first a pressure spline along each of the four temperature rows,
! then a temperature spline through the four results. t_shift/p_shift
! are 0 for the home cell and +-1 for the smoothing cells. Table,
! indices, work arrays and ierr are host-associated from compute_scv_envelope_eos; on a
! spline failure ierr is set (with a message) and the caller must return.
subroutine spline_cell(t_shift, p_shift, result)
      integer, intent(in) :: t_shift, p_shift
      real(dp), intent(out) :: result(:)

! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = idtt+i-1+t_shift
         do k = 1,4
            interp_nodes(k) = tablenv(idtt,idp+k-1+p_shift,iscvenv_log10_p)
         end do
         do j = 1,5
            do k = 1,4
               ytab_work(k) = tablenv(ii,idp+k-1+p_shift,j+1)
            end do
            call natural_spline(interp_nodes, ytab_work, 4, spline_second_deriv)
            call splint(interp_nodes, ytab_work, 4, spline_second_deriv, &
                 log10_gas_pressure, spline_value, jerr_gate)
            if (jerr_gate /= 0) then
               call report_spline_failure(jerr_gate)
               return
            end if
            temp_grid(i,j) = spline_value
         end do
      end do
! interpolate in temperature
      do k = 1,4
         interp_nodes(k) = tlogx(k+idtt-1+t_shift)
      end do
      do j = 1,5
         do k = 1,4
            ytab_work(k) = temp_grid(k,j)
         end do
         call natural_spline(interp_nodes, ytab_work, 4, spline_second_deriv)
         call splint(interp_nodes, ytab_work, 4, spline_second_deriv, &
              log10_temperature, spline_value, jerr_gate)
         if (jerr_gate /= 0) then
            call report_spline_failure(jerr_gate)
            return
         end if
         result(j) = spline_value
      end do
end subroutine spline_cell

! splint's failure (numerics' status in status): ierr, and the message
! with the point and the cell
subroutine report_spline_failure(status)
      integer, intent(in) :: status
      character(len=message_length) :: eos_lines(4)
      ierr = status
      eos_lines = unwritten_line
      write(eos_lines,'(a,a,a,2f10.5,a,2i4)') 'SCV envelope table: the spline failed (', &
           numerics_status_text(status), ') at log T, log P_gas', log10_temperature, &
           log10_gas_pressure, '; cell', idtt, idp
      call eos_report(eos_lines)
end subroutine report_spline_failure

end subroutine compute_scv_envelope_eos

!----------------------------------------------------------------------
! compute_scv_eos (formerly eqscvg)
!----------------------------------------------------------------------
! The SCV EOS at (log T, P in r) for any (X, Z), from the pure H, He and
! metal tables: at the cell around (log T, log P_gas) (find_scv_cell;
! covers false outside the tables, r untouched), 4-point Lagrange
! interpolation of each table in log P_gas and then log T, combined by
! the additive-volume rule (with the entropy of mixing) into the
! density, the mean weights it implies, dln rho/dln T and dln rho/dln P
! with radiation, cp and grad_ad; and the ionization fractions
! (scv_ionization_fractions).
subroutine compute_scv_eos(log10_temperature, hydrogen_fraction, &
     metal_fraction, r, covers)
      real(dp), intent(in) :: log10_temperature
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      type(eos_result), intent(inout) :: r
      logical, intent(out) :: covers
      real(dp) :: interp_nodes(4), temp_interp_weights(4), &
           temp_interp_weight_derivs(4), press_interp_weights(4), temp_work(4,3)

! --- locals ---
      real(dp) :: beta_complement, log10_gas_pressure, &
           helium_fraction, radiation_pressure, gas_pressure
      real(dp) :: log10_density_pure_h, density_pure_h, &
           log10_density_pure_he, density_pure_he, log10_density_pure_z, &
           density_pure_z
      real(dp) :: volume_frac_h, volume_frac_he, volume_frac_z
      real(dp) :: dlnrho_dlnp_pure_h, dlnrho_dlnp_pure_he, &
           dlnrho_dlnp_pure_z, dlnrho_dlnp_gas
      real(dp) :: dlnrho_dlnt_pure_h, dlnrho_dlnt_pure_he, &
           dlnrho_dlnt_pure_z, dlnrho_dlnt_gas, dlnp_dlnt_gas, dlnp_dlnt
      real(dp) :: dlns_dlnt_p, specific_heat_cp_gas, du_dt, &
           radiation_energy_density_term
      real(dp) :: log10_entropy_pure_h, entropy_pure_h, &
           log10_entropy_pure_he, entropy_pure_he, entropy_of_mixing, &
           dlnsmix_dlnt
      real(dp) :: dlns_dlnt_pure_h, dlns_dlnt_pure_he, &
           log10_du_dt_pure_z, du_dt_pure_z
      logical :: t_centred, p_centred
      integer :: i, ii, k
! the cell's first T row and first P_gas column (find_scv_cell)
      integer :: idtt, idp
      beta_complement = 1.0d0 - r%beta
      log10_gas_pressure = log10(r%beta*r%pressure)
! the tables' cell around the point (idtt, idp), if they cover it
      call find_scv_cell(log10_temperature, log10_gas_pressure, last_scv_cell, covers, &
           t_centred, p_centred)
      if (.not. covers) return
      idtt = last_scv_cell%first_row
      idp = last_scv_cell%first_column
      do k = 1,4
         interp_nodes(k) = tlogx(k+idtt-1)
      end do
      call interp(interp_nodes, temp_interp_weights, &
           temp_interp_weight_derivs, log10_temperature)
      do k = 1,4
         interp_nodes(k) = tablex(idtt,idp+k-1,iscv_log10_p)
      end do
      call interp(interp_nodes, press_interp_weights, x_eval=log10_gas_pressure)
      helium_fraction = 1.0d0 - hydrogen_fraction - metal_fraction
! include radiation pressure in the equation of state.
      radiation_pressure = beta_complement*r%pressure
      gas_pressure = r%beta*r%pressure
! density for x=1, y=1, z=1
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = idtt+i-1
         temp_work(i,1) = lagrange4(press_interp_weights, tablex(ii,idp:idp+3,iscv_log10_rho))
         temp_work(i,2) = lagrange4(press_interp_weights, tabley(ii,idp:idp+3,iscv_log10_rho))
         temp_work(i,3) = lagrange4(press_interp_weights, tablez(ii,idp:idp+3,iscvz_log10_rho))
      end do
! interpolate in temperature
      log10_density_pure_h = lagrange4(temp_interp_weights, temp_work(1:4,1))
      density_pure_h = exp(ln10*log10_density_pure_h)
      log10_density_pure_he = lagrange4(temp_interp_weights, temp_work(1:4,2))
      density_pure_he = exp(ln10*log10_density_pure_he)
      log10_density_pure_z = lagrange4(temp_interp_weights, temp_work(1:4,3))
      density_pure_z = exp(ln10*log10_density_pure_z)
! density
      r%density = 1.0d0/(hydrogen_fraction/density_pure_h + &
           helium_fraction/density_pure_he + metal_fraction/density_pure_z)
      r%log10_density = log10(r%density)
! the mean weights it implies (R/mu = beta P/(rho T))
      call mean_weights_from_density(r)
      volume_frac_h = hydrogen_fraction*r%density/density_pure_h
      volume_frac_he = helium_fraction*r%density/density_pure_he
      volume_frac_z = metal_fraction*r%density/density_pure_z
! d ln rho/ d ln p = qdp
! for gas pressure qdp(tot) = qdp(x)*x*rho/rho(x)+qdp(y)*y*rho/rho(y)
! for radiation pressure qdp = 0, so qdp(tot) = qdp(gas)*p/pgas
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = idtt+i-1
         temp_work(i,1) = lagrange4(press_interp_weights, tablex(ii,idp:idp+3,iscv_dlnrho_dlnp))
         temp_work(i,2) = lagrange4(press_interp_weights, tabley(ii,idp:idp+3,iscv_dlnrho_dlnp))
         temp_work(i,3) = lagrange4(press_interp_weights, tablez(ii,idp:idp+3,iscvz_dlnrho_dlnp))
      end do
! interpolate in temperature
      dlnrho_dlnp_pure_h = lagrange4(temp_interp_weights, temp_work(1:4,1))
      dlnrho_dlnp_pure_he = lagrange4(temp_interp_weights, temp_work(1:4,2))
      dlnrho_dlnp_pure_z = lagrange4(temp_interp_weights, temp_work(1:4,3))
      dlnrho_dlnp_gas = dlnrho_dlnp_pure_h*volume_frac_h + &
           dlnrho_dlnp_pure_he*volume_frac_he + dlnrho_dlnp_pure_z*volume_frac_z
      r%dlnrho_dlnp = dlnrho_dlnp_gas/r%beta
! d ln rho/ d ln t = qdt (note : d ln p/ d ln t = qpt)
! for gas pressure, correct as per qdp
! for radiation pressure, use qdt = -qdp*qpt.  correct qpt for
! radiation pressure and use the corrected qdp, qpt to get qdt.
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = idtt+i-1
         temp_work(i,1) = lagrange4(press_interp_weights, tablex(ii,idp:idp+3,iscv_dlnrho_dlnt))
         temp_work(i,2) = lagrange4(press_interp_weights, tabley(ii,idp:idp+3,iscv_dlnrho_dlnt))
         temp_work(i,3) = lagrange4(press_interp_weights, tablez(ii,idp:idp+3,iscvz_dlnrho_dlnt))
      end do
! interpolate in temperature
      dlnrho_dlnt_pure_h = lagrange4(temp_interp_weights, temp_work(1:4,1))
      dlnrho_dlnt_pure_he = lagrange4(temp_interp_weights, temp_work(1:4,2))
      dlnrho_dlnt_pure_z = lagrange4(temp_interp_weights, temp_work(1:4,3))
      dlnrho_dlnt_gas = dlnrho_dlnt_pure_h*volume_frac_h + &
           dlnrho_dlnt_pure_he*volume_frac_he + dlnrho_dlnt_pure_z*volume_frac_z
      dlnp_dlnt_gas = -dlnrho_dlnt_gas/dlnrho_dlnp_gas
      dlnp_dlnt = dlnp_dlnt_gas*r%beta + 4.0d0*beta_complement
      r%dlnrho_dlnt = -dlnp_dlnt*r%dlnrho_dlnp
! cp = s*(d ln s/ d ln t)|p is tabulated. use
! cp = du/dt + p*(d ln rho/d ln t)**2/rho/t/(d ln rho/ d ln p)
! to find du/dt.  then include the effects of radiation pressure
! on du/dt and the other terms and get a corrected cp.
! cp (gas pressure only).
! note that unlike the envelope table, we need to interpolate in
! entropy s, d ln s/d lnt, and entropy of mixing, then calculate
! du/dt
! entropy and entropy of mixing
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = idtt+i-1
         temp_work(i,1) = lagrange4(press_interp_weights, tablex(ii,idp:idp+3,iscv_log10_s))
         temp_work(i,2) = lagrange4(press_interp_weights, tabley(ii,idp:idp+3,iscv_log10_s))
         temp_work(i,3) = lagrange4(press_interp_weights, smix(ii,idp:idp+3))
      end do
! interpolate in temperature
      log10_entropy_pure_h = lagrange4(temp_interp_weights, temp_work(1:4,1))
      entropy_pure_h = exp(ln10*log10_entropy_pure_h)
      log10_entropy_pure_he = lagrange4(temp_interp_weights, temp_work(1:4,2))
      entropy_pure_he = exp(ln10*log10_entropy_pure_he)
      entropy_of_mixing = lagrange4(temp_interp_weights, temp_work(1:4,3))
      dlnsmix_dlnt = temp_interp_weight_derivs(1)*log10(temp_work(1,3)) + &
           temp_interp_weight_derivs(2)*log10(temp_work(2,3)) &
           + temp_interp_weight_derivs(3)*log10(temp_work(3,3)) &
           + temp_interp_weight_derivs(4)*log10(temp_work(4,3))
! d ln s/ d ln t (x and y) and du/dt (z)
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = idtt+i-1
         temp_work(i,1) = lagrange4(press_interp_weights, tablex(ii,idp:idp+3,iscv_dlns_dlnt))
         temp_work(i,2) = lagrange4(press_interp_weights, tabley(ii,idp:idp+3,iscv_dlns_dlnt))
         temp_work(i,3) = lagrange4(press_interp_weights, tablez(ii,idp:idp+3,iscvz_log10_du_dt))
      end do
! interpolate in temperature
      dlns_dlnt_pure_h = lagrange4(temp_interp_weights, temp_work(1:4,1))
      dlns_dlnt_pure_he = lagrange4(temp_interp_weights, temp_work(1:4,2))
      log10_du_dt_pure_z = lagrange4(temp_interp_weights, temp_work(1:4,3))
      du_dt_pure_z = exp(ln10*log10_du_dt_pure_z)
      dlns_dlnt_p = hydrogen_fraction*entropy_pure_h*dlns_dlnt_pure_h + &
           helium_fraction*entropy_pure_he*dlns_dlnt_pure_he + &
           entropy_of_mixing*dlnsmix_dlnt
      specific_heat_cp_gas = dlns_dlnt_p + metal_fraction*(du_dt_pure_z + &
           gas_pressure*dlnrho_dlnt_pure_z**2/dlnrho_dlnp_pure_z/ &
           density_pure_z/r%temperature)
! now find du/dt from the original table.
      du_dt = specific_heat_cp_gas - gas_pressure*dlnrho_dlnt_gas**2/ &
           dlnrho_dlnp_gas/r%density/r%temperature
! correct du/dt for radiation
      radiation_energy_density_term = 3.0d0*radiation_pressure/r%density
      du_dt = du_dt + 4.0d0*radiation_energy_density_term/r%temperature
! correct cp for radiation pressure
      r%specific_heat_cp = du_dt + r%pressure*r%dlnrho_dlnt**2/r%dlnrho_dlnp/ &
           r%density/r%temperature
! adiabatic temperature gradient
      r%adiabatic_gradient = -r%pressure*r%dlnrho_dlnt/r%density/r%temperature/ &
           r%specific_heat_cp

! the ionization fractions, from the pure H and He tables at the cell
      call scv_ionization_fractions(idtt, idp, temp_interp_weights, press_interp_weights, &
           r%ion_fraction)

end subroutine compute_scv_eos

!----------------------------------------------------------------------
! find_scv_cell
!----------------------------------------------------------------------
! The 4 x 4 cell of the SCV tables around (log T, log P_gas): its first
! row idtt and first column idp, on the H table's grid (tlogx and the log
! P_gas column, which the He, Z and envelope tables share), searched from
! the cell given (the last lookup's), which becomes this one; a column
! below the table's first is clamped to it. covers is false outside the
! tables: log T beyond their rows, log P_gas below
! scv_min_log10_gas_pressure or above the row's highest (the cell is then
! left as it was found). t_centred and p_centred tell whether the cell
! is centred on the point (not clamped at the table's edge), where the
! envelope lookup smooths across cells.
subroutine find_scv_cell(log10_temperature, log10_gas_pressure, cell, covers, t_centred, &
     p_centred)
      real(dp), intent(in) :: log10_temperature, log10_gas_pressure
      type(scv_cell), intent(inout) :: cell
      logical, intent(out) :: covers, t_centred, p_centred
      integer, parameter :: nts = scv_nt
      integer :: i, ii, j, jj, jjj, idtt, idp

      idtt = cell%first_row
      idp = cell%first_column
      covers = .false.
      t_centred = .false.
      p_centred = .false.
      if (log10_temperature < tlogx(1) .or. log10_temperature > tlogx(nts) .or. &
           log10_gas_pressure < scv_min_log10_gas_pressure) return
! the rows: down or up from the last cell
      if (log10_temperature < tlogx(idtt+1)) then
         do i = idtt, 1, -1
            if (log10_temperature > tlogx(i)) then
               ii = i - 1
               exit
            end if
         end do
         if (i < 1) ii = 1
      else
         do i = idtt+2, nts
            if (log10_temperature < tlogx(i)) then
               ii = i - 2
               exit
            end if
         end do
         if (i > nts) ii = nts - 3
      end if
      idtt = max(1,ii)
      idtt = min(nts-3,idtt)
      t_centred = idtt == ii
      cell%first_row = idtt
! the columns, at the cell's first row (the one with the fewest points)
      jjj = min(nptsx(idtt)-3, idp)
      if (log10_gas_pressure < tablex(idtt,jjj+1,iscv_log10_p)) then
         do j = jjj, 1, -1
            if (log10_gas_pressure > tablex(idtt,j,iscv_log10_p)) then
               jj = j - 1
               exit
            end if
         end do
         if (j < 1) jj = 1
         idp = max(1,jj)
         idp = min(nptsx(idtt)-3, idp)
      else
         do j = jjj+2, nptsx(idtt)
            if (log10_gas_pressure < tablex(idtt,j,iscv_log10_p)) then
               jj = j - 2
               exit
            end if
         end do
         if (j > nptsx(idtt)) return
         idp = max(1, min(nptsx(idtt)-3, jj))
      end if
      covers = .true.
      p_centred = idp == jj
      cell%first_column = idp
end subroutine find_scv_cell

!----------------------------------------------------------------------
! scv_ionization_fractions
!----------------------------------------------------------------------
! The H+, He+ and He++ fractions of the nuclei at the cell (idtt, idp),
! with the 4-point weights in T and P: the pure H table's particle
! fractions of H2 and H (columns 2 and 3) and the pure He table's of
! neutral He and He+, interpolated in P at the cell's four rows and then
! in T; the H+, electron and He++ fractions follow from the
! conservation of particles and of charge.
subroutine scv_ionization_fractions(idtt, idp, temp_interp_weights, press_interp_weights, &
     ion_fraction)
      integer, intent(in) :: idtt, idp
      real(dp), intent(in) :: temp_interp_weights(:), press_interp_weights(:)
      real(dp), intent(out) :: ion_fraction(:)
      real(dp) :: temp_grid(4,2)
      real(dp) :: xtf_h2, xtf_he, xtf_h1, xtf_hep, xtf_h_e, xtf_hp, xtf_he_e, xtf_hepp
      integer :: i, ii

! the fractions of all particles (electrons included) that are H2 and
! neutral He (column 2 of the H and He tables)
      do i = 1,4
         ii = idtt+i-1
         temp_grid(i,1) = lagrange4(press_interp_weights, tablex(ii,idp:idp+3,iscv_frac_neutral))
         temp_grid(i,2) = lagrange4(press_interp_weights, tabley(ii,idp:idp+3,iscv_frac_neutral))
      end do
      xtf_h2 = lagrange4(temp_interp_weights, temp_grid(1:4,1))
      xtf_he = lagrange4(temp_interp_weights, temp_grid(1:4,2))
! that are neutral H and He+ (column 3)
      do i = 1,4
         ii = idtt+i-1
         temp_grid(i,1) = lagrange4(press_interp_weights, tablex(ii,idp:idp+3,iscv_frac_atom_ion1))
         temp_grid(i,2) = lagrange4(press_interp_weights, tabley(ii,idp:idp+3,iscv_frac_atom_ion1))
      end do
      xtf_h1 = lagrange4(temp_interp_weights, temp_grid(1:4,1))
      xtf_hep = lagrange4(temp_interp_weights, temp_grid(1:4,2))
! the hydrogen-related electrons and H+, the helium-related electrons
! and He++, by the conservation of particles and charge
      xtf_h_e = .5d0*(1d0 - xtf_h2 - xtf_h1)
      xtf_hp = xtf_h_e
      xtf_he_e = (1d0/3d0)*(2d0 - 2d0*xtf_he - xtf_hep)
      xtf_hepp = 1d0 - xtf_he - xtf_hep - xtf_he_e
! the fractions of the H nuclei that are H+, and of the He nuclei that
! are He+ and He++
      ion_fraction(1) = xtf_hp / (xtf_h2 + xtf_h1 + xtf_hp)
      ion_fraction(2) = xtf_hep / (xtf_he + xtf_hep + xtf_hepp)
      ion_fraction(3) = xtf_hepp / (xtf_he + xtf_hep + xtf_hepp)
end subroutine scv_ionization_fractions

!----------------------------------------------------------------------
! build_scv_envelope_table (formerly setscv)
!----------------------------------------------------------------------
! Builds the envelope-mixture (X/Y/Z) SCVH-style equation-of-state
! table tablenv from the separately tabulated pure-hydrogen (tablex),
! pure-helium (tabley), and pure-metal (tablez) tables: first the
! entropy of mixing (smix) on the native H/He table grid, then the
! three species combined by the additive-volume rule at the envelope
! composition (X = envelope_hydrogen_fraction, Z =
! envelope_metal_fraction) into tablenv's six columns, which
! compute_scv_envelope_eos interpolates. Nothing to do before
! read_scv_tables has read the tables, or when the table is built for
! this composition already. The lookups' warm-start cell (last_scv_cell)
! is not touched.
subroutine build_scv_envelope_table(envelope_hydrogen_fraction, &
     envelope_metal_fraction)
! the envelope composition the table is built for
      real(dp), intent(in) :: envelope_hydrogen_fraction, &
           envelope_metal_fraction

      integer, parameter :: nts = scv_nt

! interp_x: the three abscissae handed to inter3, and its weights and
! their derivatives
      real(dp) :: interp_x(3), t_interp_weight(3), t_interp_dweight(3)
! --- locals ---
      integer :: t_idx, p_idx, k_idx, t_row
! log10 T and log10 P of the table point, ln T, and T
      real(dp) :: log_t_work, log_p_work, ln_t_work, temp_value
! helium_fraction_hhe: Y = 1 - X of the pure H/He mixture used for the
! entropy of mixing; helium_fraction_local: Y = 1 - X - Z of the
! envelope mixture.
      real(dp) :: helium_fraction_hhe, one_minus_y_local
      real(dp) :: helium_fraction_local, &
           hydrogen_fraction_local, metal_fraction_local
      real(dp) :: density_h, density_he, density_z, density_mix
      real(dp) :: hydrogen_number_density, helium_number_density, &
           hydrogen_electron_density, helium_electron_density, &
           total_electron_density
      real(dp) :: mixing_beta, mixing_gamma, hydrogen_xd, &
           helium_yd, mix_de, mix_dei, total_number_density
      real(dp) :: pressure_value, radiation_pressure, gas_pressure, &
           rho_mix, log_rho_mix
      real(dp) :: dlnrho_dlnp_gas, dlnrho_dlnt_gas
      real(dp) :: entropy_h, entropy_he, entropy_total, &
           dsmix_dlnt, dlns_dlnt, cp_gas, du_dt

      if (.not. scv_tables_loaded) return
      if (envelope_hydrogen_fraction == tablenv_hydrogen_fraction .and. &
           envelope_metal_fraction == tablenv_metal_fraction) return
!  Read in equation of state tables for hydrogen and helium
      do t_idx = 1, nts
         do p_idx = 1, nptsx(t_idx)
! Compute the entropy of mixing.
! Number density of hydrogen and helium
            helium_fraction_hhe = 1.0d0-envelope_hydrogen_fraction
            one_minus_y_local = 1.0d0 - helium_fraction_hhe
            density_h = exp(ln10*tablex(t_idx,p_idx,iscv_log10_rho))
            density_he = exp(ln10*tabley(t_idx,p_idx,iscv_log10_rho))
            density_mix = 1.0d0/(one_minus_y_local/density_h + &
                 helium_fraction_hhe/density_he)
            hydrogen_number_density = 2.0d0*one_minus_y_local*density_mix/ &
                 hydrogen_atom_mass_cgs/(1.0d0+3.0d0*tablex(t_idx,p_idx,iscv_frac_neutral)+ &
                 tablex(t_idx,p_idx,iscv_frac_atom_ion1))
            helium_number_density = 3.0d0*helium_fraction_hhe*density_mix/ &
                 helium_atom_mass_cgs/(1.0d0 + 2.0d0*tabley(t_idx,p_idx,iscv_frac_neutral) + &
                 tabley(t_idx,p_idx,iscv_frac_atom_ion1))
! Hydrogen and helium electron number densities
            hydrogen_electron_density = 0.5d0*hydrogen_number_density* &
                 (1.0d0 - tablex(t_idx,p_idx,iscv_frac_neutral) - tablex(t_idx,p_idx,iscv_frac_atom_ion1))
            helium_electron_density = cc13*helium_number_density* &
                 (2.0d0 - 2.0d0*tabley(t_idx,p_idx,iscv_frac_neutral) - tabley(t_idx,p_idx,iscv_frac_atom_ion1))
            total_electron_density = hydrogen_electron_density + &
                 helium_electron_density
            hydrogen_electron_density = max(hydrogen_electron_density,1.0d0)
            helium_electron_density = max(helium_electron_density,1.0d0)
            total_electron_density = max(total_electron_density,1.0d0)
            if (hydrogen_number_density == 0.0d0 .or. &
                 helium_number_density == 0.0d0) then
               smix(t_idx,p_idx) = 0.0d0
            else
               mixing_beta = (hydrogen_atom_mass_cgs/helium_atom_mass_cgs)* &
                    (helium_fraction_hhe/one_minus_y_local)
               mixing_gamma = 1.5d0*(1.0d0 + tablex(t_idx,p_idx,iscv_frac_atom_ion1)+ 3.0d0* &
                    tablex(t_idx,p_idx,iscv_frac_neutral))/(1.0d0+2.0d0*tabley(t_idx,p_idx,iscv_frac_neutral)+ &
                    tabley(t_idx,p_idx,iscv_frac_atom_ion1))
               hydrogen_xd = 1.0d0 - tablex(t_idx,p_idx,iscv_frac_neutral) - tablex(t_idx,p_idx,iscv_frac_atom_ion1)
               helium_yd = 2.0d0 - 2.0d0*tabley(t_idx,p_idx,iscv_frac_neutral) - &
                    tabley(t_idx,p_idx,iscv_frac_atom_ion1)
               if (hydrogen_xd > 0.0d0) then
                 mix_de = 1.5d0*helium_yd*mixing_beta*mixing_gamma/hydrogen_xd
               else
                 mix_de = 0.0d0
               end if
               if (helium_yd > 0.0d0) then
                  mix_dei = hydrogen_xd/(1.5d0*helium_yd*mixing_beta*mixing_gamma)
               else
                  mix_dei = 0.0d0
               end if
               total_number_density = hydrogen_number_density + &
                    helium_number_density
               smix(t_idx,p_idx) = one_minus_y_local/hydrogen_atom_mass_cgs*2.0d0/ &
                    (1.0d0+tablex(t_idx,p_idx,iscv_frac_atom_ion1)+3.0d0*tablex(t_idx,p_idx,iscv_frac_neutral))* &
                    (log(1.0d0+mixing_beta*mixing_gamma)- &
                    hydrogen_electron_density/total_number_density* &
                    log(1.0d0 + mix_de)+mixing_beta*mixing_gamma* &
                    log(1.0d0+1.0d0/mixing_beta/mixing_gamma)- &
                    helium_electron_density/total_number_density* &
                    log(1.0d0+mix_dei))
               smix(t_idx,p_idx) = boltzmann_constant_cgs*smix(t_idx,p_idx)
            end if
         end do
      end do
!  Now compute equation of state variables for the surface mixture.
      do t_idx=1,nts
         log_t_work = tlogx(t_idx)
         temp_value = (exp10(log_t_work))
! the three-row stencil in T for d smix/d ln T (its first row t_row)
        if (t_idx == 1) then
           t_row = 1
        else if (t_idx == nts) then
           t_row = t_idx - 2
        else
           t_row = t_idx - 1
        end if
        do k_idx = 1,3
           interp_x(k_idx) = ln10*tlogx(k_idx+t_row-1)
        end do
        ln_t_work = ln10*log_t_work
        call inter3(interp_x,t_interp_weight,t_interp_dweight,ln_t_work)
         do p_idx=1,nptsx(t_idx)
            log_p_work = tablex(t_idx,p_idx,iscv_log10_p)
            hydrogen_fraction_local = envelope_hydrogen_fraction
            metal_fraction_local = envelope_metal_fraction
            helium_fraction_local = 1.0d0 - hydrogen_fraction_local - &
                 metal_fraction_local
            pressure_value = (exp10(log_p_work))
! Include radiation pressure in the equation of state.
            radiation_pressure = radiation_constant_over_3*(temp_value**2)**2
            gas_pressure = pressure_value
            pressure_value = (pressure_value + radiation_pressure)
            log_p_work = log10(pressure_value)
!  Density : Add inversely : 1/rho = X/rho(X) + Y/rho(Y) + Z/rho(Z)
            density_h = exp(ln10*tablex(t_idx,p_idx,iscv_log10_rho))
            density_he = exp(ln10*tabley(t_idx,p_idx,iscv_log10_rho))
            density_z = exp(ln10*tablez(t_idx,p_idx,iscvz_log10_rho))
            rho_mix = 1.0d0/(hydrogen_fraction_local/density_h + &
                 helium_fraction_local/density_he + &
                 metal_fraction_local/density_z)
            log_rho_mix = log10(rho_mix)
! D ln rho/ d ln P = qdp
! for gas pressure qdp(tot) = qdp(X)*X*rho/rho(X)+qdp(Y)*Y*rho/rho(Y)
! for radiation pressure qdp = 0, so qdp(tot) = qdp(gas)*P/pgas
            dlnrho_dlnp_gas = tablex(t_idx,p_idx,iscv_dlnrho_dlnp)*hydrogen_fraction_local* &
                 (rho_mix/density_h)+tabley(t_idx,p_idx,iscv_dlnrho_dlnp)* &
                 helium_fraction_local*(rho_mix/density_he) &
                 + tablez(t_idx,p_idx,iscvz_dlnrho_dlnp)*metal_fraction_local* &
                 (rho_mix/density_z)
! D ln rho/ d ln T = qdt (note : D ln P/ d ln T = qpt)
! for gas pressure, correct as per qdp
! for radiation pressure, use qdt = qdp*qpt.  Correct qpt for
! radiation pressure and use the corrected qdp, qpt to get qdt.
            dlnrho_dlnt_gas = hydrogen_fraction_local*(rho_mix/density_h)* &
                 tablex(t_idx,p_idx,iscv_dlnrho_dlnt)+helium_fraction_local* &
                 (rho_mix/density_he)*tabley(t_idx,p_idx,iscv_dlnrho_dlnt) &
                 +metal_fraction_local*(rho_mix/density_z)* &
                 tablez(t_idx,p_idx,iscvz_dlnrho_dlnt)
            entropy_h = exp(ln10*tablex(t_idx,p_idx,iscv_log10_s))
            entropy_he = exp(ln10*tabley(t_idx,p_idx,iscv_log10_s))
! Entropy.  Obeys additive volume rule, but also need to include the
! entropy of mixing.
            entropy_total = hydrogen_fraction_local*entropy_h + &
                 helium_fraction_local*entropy_he + smix(t_idx,p_idx)
            dsmix_dlnt = t_interp_dweight(1)*smix(t_row,p_idx)+ &
                 t_interp_dweight(2)*smix(t_row+1,p_idx)+ &
                 t_interp_dweight(3)*smix(t_row+2,p_idx)
! D LN S/ D LN T
            dlns_dlnt = (hydrogen_fraction_local*entropy_h* &
                 tablex(t_idx,p_idx,iscv_dlns_dlnt) + helium_fraction_local*entropy_he* &
                 tabley(t_idx,p_idx,iscv_dlns_dlnt) + dsmix_dlnt)/entropy_total
!  Cp = s*(d ln s/ d ln T)|P is tabulated. Use
!  cp = du/dt + P*(d ln rho/d ln T)**2/rho/T/(d ln rho/ d ln P)
!  to include the effects of radiation pressure.
! Cp (gas pressure only).
            cp_gas = entropy_total*dlns_dlnt + metal_fraction_local* &
                 (exp(ln10*tablez(t_idx,p_idx,iscvz_log10_du_dt)) + &
                 gas_pressure*tablez(t_idx,p_idx,iscvz_dlnrho_dlnt)**2/ &
                 tablez(t_idx,p_idx,iscvz_dlnrho_dlnp)/density_z/temp_value)
! Now find du/dt from the original table.
            du_dt = cp_gas - gas_pressure*dlnrho_dlnt_gas**2/ &
                 dlnrho_dlnp_gas/rho_mix/temp_value
! Now store the relevant variables in a table for the envelope
! MIXTURE.
            tablenv(t_idx,p_idx,iscvenv_log10_p) = tablex(t_idx,p_idx,iscv_log10_p)
            tablenv(t_idx,p_idx,iscvenv_log10_rho) = log_rho_mix
            tablenv(t_idx,p_idx,iscvenv_dlnrho_dlnt) = dlnrho_dlnt_gas
            tablenv(t_idx,p_idx,iscvenv_dlnrho_dlnp) = dlnrho_dlnp_gas
            tablenv(t_idx,p_idx,iscvenv_cp) = cp_gas
            tablenv(t_idx,p_idx,iscvenv_du_dt) = du_dt
         end do
      end do
      tablenv_hydrogen_fraction = envelope_hydrogen_fraction
      tablenv_metal_fraction = envelope_metal_fraction

end subroutine build_scv_envelope_table

!======================================================================
! Sources: Saumon, Chabrier & Van Horn (1995) for the H and He tables;
! the metal table was computed with the Prather EOS of the old Yale code.
! The lookups were EQSCVE (compute_scv_envelope_eos) and EQSCVG
! (compute_scv_eos), the table build SETSCV. The upward smoothing
! weights, 0.5 at the node and falling to 0 at the tolerance, were made
! the mirror of the downward ones in 2026 (the inherited 0.5 dist/tol
! was reversed, a jump at 60% of every cell). The search start became a
! defined cell (scv_cell) and the column search was clamped to the table
! in 2026 (M107).

end module scv_eos
