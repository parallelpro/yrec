!----------------------------------------------------------------------
! scv_eos
!----------------------------------------------------------------------
! The Saumon-Chabrier-Van Horn equation of state for cool, dense
! H/He gas, from the pure-H, pure-He and pure-metal tables.
!     build_scv_envelope_table  combine the three tables at the envelope
!                               X and Z (at eos_init and when the
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
      use scv_eos_tables, only: idp, idtt, iscv_dlnrho_dlnp, iscv_dlnrho_dlnt, &
           iscv_dlns_dlnt, iscv_du_dt, iscv_frac_atom_ion1, iscv_frac_neutral, &
           iscv_log10_p, iscv_log10_rho, iscv_log10_s, iscvenv_cp, &
           iscvenv_dlncp_dlnt, iscvenv_dlnp_dlnrho, iscvenv_dlnrho_dlnp, &
           iscvenv_dlnrho_dlnt, iscvenv_dqdt_dlnp, iscvenv_dqdt_dlnt, &
           iscvenv_dqut_dlnp, iscvenv_dqut_dlnt, iscvenv_du_dt, &
           iscvenv_log10_p, iscvenv_log10_rho, iscvz_dlnrho_dlnp, &
           iscvz_dlnrho_dlnt, iscvz_log10_du_dt, iscvz_log10_rho, nptsx, &
           scv_nt, scv_weighted_sum4, smix, tablenv, tablex, tabley, tablez, &
           tlogx
      use physical_constants, only: cc13, ln10, radiation_constant_over_3
      use numerics_lib, only: inter3, interp, kspline, ksplint
      use io_units, only: run_log_unit
      use eos_mixture, only: eos_mix
      use math_lib
      implicit none
      private
      public :: compute_scv_envelope_eos, build_scv_envelope_table, read_scv_tables

contains

!----------------------------------------------------------------------
! compute_scv_envelope_eos (formerly eqscve)
!----------------------------------------------------------------------
! Interpolates in the surface (envelope) abundance table of the
! Saumon-Chabrier-Van Horn (SCV) equation of state -- i.e. the table
! built specifically for the envelope composition (X, Z). If the
! requested composition differs from the envelope composition, falls
! back to compute_scv_eos (the general-abundance SCV table). Locates the SCV
! table cell in (log10 T, log10 Pgas), interpolates density and its
! logarithmic derivatives, the specific heat, adiabatic gradient, and
! ionization fractions, with an additional smoothing pass across
! adjacent table cells (the 4-pt Lagrange interpolation of the
! original version was replaced with a 4-pt natural spline).
subroutine compute_scv_envelope_eos(log10_temperature, hydrogen_fraction, &
     metal_fraction, r, valid_table_point, ierr)

      integer, intent(out) :: ierr
      integer :: jerr_gate

      integer, parameter :: nts = scv_nt

      type(eos_result), intent(inout) :: r
real(dp), intent(in) :: log10_temperature
real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      logical, intent(out) :: valid_table_point

! Replace 4-pt Lagrange interpolate with 4-pt natural spline
      real(dp) :: press_interp_weights(4), &
           press_interp_weight_derivs(4), temp_interp_weights(4), &
           temp_interp_weight_derivs(4)
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
           dlnp_dlnt, specific_heat_cp_gas, du_dt, &
           radiation_energy_density_term
      real(dp) :: xtf_h2, xtf_he, xtf_h1, xtf_hep, xtf_h_e, xtf_hp, &
           xtf_he_e, xtf_hepp, xhp, xhep, xhepp
      logical :: temp_needs_smoothing, press_needs_smoothing
      integer :: temp_smooth_direction, press_smooth_direction
      integer :: i, ii, j, jj, jjj, k
      ierr = 0
      if (abs(hydrogen_fraction-eos_mix%envelope_hydrogen_fraction) > 1.0d-5 &
           .or. abs(metal_fraction-eos_mix%envelope_metal_fraction) > 1.0d-5) then
         call compute_scv_eos(log10_temperature, hydrogen_fraction, metal_fraction, &
              r, valid_table_point)
         return
      end if
! find gas pressure (which is the quantity which is tabulated).
      beta_complement = 1.0d0 - r%beta
      log10_gas_pressure = log10(r%beta*r%pressure)
!  check if the point is within the table
      if (log10_temperature < tlogx(1) .or. &
           log10_temperature > tlogx(nts) .or. &
           log10_gas_pressure < 4.0d0) then
         valid_table_point = .false.   ! Error exit - no valid table entry
         return
      end if
!  find nearest points in temperature.
      if (log10_temperature < tlogx(idtt+1)) then
! search down to find nearest 4 table elements
         do i = idtt, 1, -1
            if (log10_temperature > tlogx(i)) then
               ii = i - 1
               exit
            end if
         end do
         if (i < (1)) then
         ii = 1
         end if
      else
! search up for nearest 4 table elements
         do i = idtt+2, nts
            if (log10_temperature < tlogx(i)) then
               ii = i - 2
               exit
            end if
         end do
         if (i > nts) then
         ii = nts - 3
         end if
      end if
      idtt = max(1,ii)
      idtt = min(nts-3,idtt)
      temp_needs_smoothing = .false.
      if (idtt == ii) then
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
! 2026 (bugsweep Batch 2): weight must be 0.5 AT the node and fall
! to 0 at dist = tol (mirror of the downward branch above); the
! inherited 0.5*dist/tol was reversed, giving a 0.5|S_n+1 - S_n|
! jump at 60% of every cell instead of a continuous blend.
               temp_smooth_weight = 0.5d0*(1.0d0 - temp_dist_below/tol_temp_smooth)
            end if
         end if
      end if
!  find nearest points in pressure.
      jjj = min(nptsx(idtt)-3, idp)
      if (log10_gas_pressure < tablenv(idtt,jjj+1,iscvenv_log10_p)) then
! search down to find nearest 4 table elements
         do j = jjj, 1, -1
            if (log10_gas_pressure > tablenv(idtt,j,iscvenv_log10_p)) then
               jj = j - 1
               exit
            end if
         end do
         if (j < (1)) then
         jj = 1
         end if
         idp = max(1,jj)
         idp = min(nptsx(idtt)-3, &
              idp)
      else
! search up for nearest 4 table elements.  note search is done at lowest
! temperature point (with the minimum range in p).
         do j = jjj+2, nptsx(idtt)
            if (log10_gas_pressure < tablenv(idtt,j,iscvenv_log10_p)) then
               jj = j - 2
               exit
            end if
         end do
         if (j > (nptsx(idtt))) then
! point is outside table; return.
         valid_table_point = .false.   ! Error exit - no valid table entry
         return
         end if
         idp = min(nptsx(idtt)-3, jj)
      end if
      valid_table_point = .true.
      press_needs_smoothing = .false.
      if (idp == jj) then
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
! 2026 (bugsweep Batch 2): same reversal as the T branch above.
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
!  density
      r%log10_density = interp_result(1)
      r%density = exp(ln10*r%log10_density)
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
!  cp = s*(d ln s/ d ln t)|p is tabulated. use
!  cp = du/dt + p*(d ln rho/d ln t)**2/rho/t/(d ln rho/ d ln p)
!  to find du/dt.  then include the effects of radiation pressure
!  on du/dt and the otehr terms and get a corrected cp.
! cp (gas pressure only).
      specific_heat_cp_gas = interp_result(4)
! now find du/dt from the original table.
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

! Get fractions of total particles (including electrons), as follows:
!   XTF_H2  the fraction that is neutral hydrogen molecules, and
!   XTF_He  the fraction thst is neutral helium).
!   These are in column 2 respectively of the SCV hydrogen and
!   Helium tables.
! interpolate in pressure at 4 different temperature points.

      do k = 1,4
         interp_nodes(k) = tlogx(k+idtt-1)
      end do
      call interp(interp_nodes, temp_interp_weights, &
           temp_interp_weight_derivs, log10_temperature)
      do k = 1,4
         interp_nodes(k) = tablenv(idtt, &
              idp+k-1,iscvenv_log10_p)
      end do
      call interp(interp_nodes, press_interp_weights, &
           press_interp_weight_derivs, log10_gas_pressure)

      do i = 1,4
         ii = idtt+i-1
         temp_grid(i,1) = scv_weighted_sum4(press_interp_weights, tablex(ii,idp:idp+3,iscv_frac_neutral))
         temp_grid(i,2) = scv_weighted_sum4(press_interp_weights, tabley(ii,idp:idp+3,iscv_frac_neutral))
      end do
! interpolate in temperature
      xtf_h2 = scv_weighted_sum4(temp_interp_weights, temp_grid(1:4,1))
      xtf_he = scv_weighted_sum4(temp_interp_weights, temp_grid(1:4,2))

! Get more fractions of total particles (including electrons), as follows:
!   XTF_H1  the fraction that is neutral hydrogen atoms, and
!   XTF_HeP  the fraction thst is singly ionized helium).
!   These are in column 3 respectively of the SCV hydrogen and
!   Helium tables.
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = idtt+i-1
         temp_grid(i,1) = scv_weighted_sum4(press_interp_weights, tablex(ii,idp:idp+3,iscv_frac_atom_ion1))
         temp_grid(i,2) = scv_weighted_sum4(press_interp_weights, tabley(ii,idp:idp+3,iscv_frac_atom_ion1))
      end do
! interpolate in temperature
      xtf_h1 = scv_weighted_sum4(temp_interp_weights, temp_grid(1:4,1))
      xtf_hep = scv_weighted_sum4(temp_interp_weights, temp_grid(1:4,2))

! At tis time we can calculate the ramaining particles fractions using
! conservation of particles and coservation of charge.  The variable
! names are:
!  XTF_H_e  the hycrogen related fraction that is electrons
!  XTF_HP   the hydrogen ralated fraction that is H+ ions
!  XTF_He_e the helium related fraction that is electrons
!  XTF_HePP the helium related fraction that is doubly ionized helium (He++)
! Particle and charge conservation yields:
      xtf_h_e = .5d0*(1d0 - xtf_h2 - xtf_h1)
      xtf_hp = xtf_h_e
      xtf_he_e = (1d0/3d0)*(2d0 - 2d0*xtf_he - xtf_hep)
      xtf_hepp = 1d0 - xtf_he - xtf_hep - xtf_he_e

! We are seeking the fraction XHP of hydrogen nuclei that are singly
! ionized, the fraction XHeP of helium huclei that are singly ionized
! and XHePP, the helium fraction that is doubly ionized.
      xhp = xtf_hp / (xtf_h2 + xtf_h1 + xtf_hp)
      xhep = xtf_hep / (xtf_he + xtf_hep + xtf_hepp)
      xhepp = xtf_hepp / (xtf_he + xtf_hep + xtf_hepp)

! In the same order, these are the desired elements of array ion_fraction(3)
      r%ion_fraction(1) = xhp
      r%ion_fraction(2) = xhep
      r%ion_fraction(3) = xhepp

      return  ! Normal exit. Valid table entry. valid_table_point is true.

contains

! Bicubic-spline interpolation of the five envelope-table variables
! (columns 2-6 of tablenv) at (log10_gas_pressure, log10_temperature)
! on the 4x4 cell whose corner is (idtt + t_shift, idp + p_shift):
! first a pressure spline along each of the four temperature rows,
! then a temperature spline through the four results. t_shift/p_shift
! are 0 for the home cell and +-1 for the smoothing cells. Table,
! indices, work arrays and ierr are host-associated from compute_scv_envelope_eos; on a
! spline failure ierr is set and the caller must return.
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
            call kspline(interp_nodes, ytab_work, spline_second_deriv)
            call ksplint(interp_nodes, ytab_work, spline_second_deriv, &
                 log10_gas_pressure, spline_value, jerr_gate)
            ! 2026 numerics-gate opt-in (see numerics_lib): interpolation
            ! failure returns via ierr instead of stopping.
            if (jerr_gate /= 0) then
               ierr = jerr_gate
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
         call kspline(interp_nodes, ytab_work, spline_second_deriv)
         call ksplint(interp_nodes, ytab_work, spline_second_deriv, &
              log10_temperature, spline_value, jerr_gate)
         if (jerr_gate /= 0) then
            ierr = jerr_gate
            return
         end if
         result(j) = spline_value
      end do
end subroutine spline_cell

end subroutine compute_scv_envelope_eos

!----------------------------------------------------------------------
! compute_scv_eos (formerly eqscvg)
!----------------------------------------------------------------------
! Interpolates in the general-abundance tables of the Saumon-
! Chabrier-Van Horn (SCV) equation of state. Given log10(T), the
! total pressure P, and the composition (X, Z), locates the SCV
! table cell in (log10 T, log10 Pgas) and interpolates density and
! its logarithmic derivatives, plus the specific heat, adiabatic
! gradient, and ionization fractions.
subroutine compute_scv_eos(log10_temperature, hydrogen_fraction, &
     metal_fraction, r, valid_table_point)


      integer, parameter :: nts = scv_nt

      type(eos_result), intent(inout) :: r
real(dp), intent(in) :: log10_temperature
real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      logical, intent(out) :: valid_table_point

      real(dp) :: interp_nodes(4), temp_interp_weights(4), &
           temp_interp_weight_derivs(4), press_interp_weights(4), &
           press_interp_weight_derivs(4), temp_work(4,3)

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
      real(dp) :: xtf_h2, xtf_he, xtf_h1, xtf_hep, xtf_h_e, xtf_hp, &
           xtf_he_e, xtf_hepp, xhp, xhep, xhepp
      integer :: i, ii, j, jj, jjj, k
      beta_complement = 1.0d0 - r%beta
      log10_gas_pressure = log10(r%beta*r%pressure)
! check if the point is within the table
      if (log10_temperature < tlogx(1) .or. &
           log10_temperature > tlogx(nts) .or. &
           log10_gas_pressure < 4.0d0) then
         valid_table_point = .false.
         return
      end if
! find nearest points in temperature.
      if (log10_temperature < tlogx(idtt+1)) then
! search down to find nearest 4 table elements
         do i = idtt, 1, -1
            if (log10_temperature > tlogx(i)) then
               ii = i - 1
               exit
            end if
         end do
         if (i < (1)) then
         ii = 1
         end if
         idtt = max(1,ii)
         idtt = min(nts-3,idtt)
      else
! search up for nearest 4 table elements
         do i = idtt+2, nts
            if (log10_temperature < tlogx(i)) then
               ii = i - 2
               exit
            end if
         end do
         if (i > nts) then
         ii = nts - 3
         end if
         idtt = max(1,ii)
         idtt = min(nts-3,idtt)
      end if
! find nearest points in pressure.
      jjj = min(nptsx(idtt)-3, idp)
      if (log10_gas_pressure < tablenv(idtt,jjj+1,iscvenv_log10_p)) then
! search down to find nearest 4 table elements
         do j = jjj, 1, -1
            if (log10_gas_pressure > tablenv(idtt,j,iscvenv_log10_p)) then
               jj = j - 1
               exit
            end if
         end do
         if (j < (1)) then
         jj = 1
         end if
         idp = max(1,jj)
         idp = min(nptsx(idtt)-3, &
              idp)
      else
! search up for nearest 4 table elements.  note search is done at lowest
! temperature point (with the minimum range in p).
         do j = jjj+2, nptsx(idtt)
            if (log10_gas_pressure < tablenv(idtt,j,iscvenv_log10_p)) then
               jj = j - 2
               exit
            end if
         end do
         if (j > (nptsx(idtt))) then
! point is outside table; return.
         valid_table_point = .false.
         return
         end if
         idp = min(nptsx(idtt)-3, jj)
      end if
      valid_table_point = .true.
      do k = 1,4
         interp_nodes(k) = tlogx(k+idtt-1)
      end do
      call interp(interp_nodes, temp_interp_weights, &
           temp_interp_weight_derivs, log10_temperature)
      do k = 1,4
         interp_nodes(k) = tablenv(idtt, &
              idp+k-1, iscvenv_log10_p)
      end do
      call interp(interp_nodes, press_interp_weights, &
           press_interp_weight_derivs, log10_gas_pressure)
      helium_fraction = 1.0d0 - hydrogen_fraction - metal_fraction
! include radiation pressure in the equation of state.
      radiation_pressure = beta_complement*r%pressure
      gas_pressure = r%beta*r%pressure
! density for x=1, y=1, z=1
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = idtt+i-1
         temp_work(i,1) = scv_weighted_sum4(press_interp_weights, tablex(ii,idp:idp+3,iscv_log10_rho))
         temp_work(i,2) = scv_weighted_sum4(press_interp_weights, tabley(ii,idp:idp+3,iscv_log10_rho))
         temp_work(i,3) = scv_weighted_sum4(press_interp_weights, tablez(ii,idp:idp+3,iscvz_log10_rho))
      end do
! interpolate in temperature
      log10_density_pure_h = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,1))
      density_pure_h = exp(ln10*log10_density_pure_h)
      log10_density_pure_he = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,2))
      density_pure_he = exp(ln10*log10_density_pure_he)
      log10_density_pure_z = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,3))
      density_pure_z = exp(ln10*log10_density_pure_z)
! density
      r%density = 1.0d0/(hydrogen_fraction/density_pure_h + &
           helium_fraction/density_pure_he + metal_fraction/density_pure_z)
      r%log10_density = log10(r%density)
      volume_frac_h = hydrogen_fraction*r%density/density_pure_h
      volume_frac_he = helium_fraction*r%density/density_pure_he
      volume_frac_z = metal_fraction*r%density/density_pure_z
! d ln rho/ d ln p = qdp
! for gas pressure qdp(tot) = qdp(x)*x*rho/rho(x)+qdp(y)*y*rho/rho(y)
! for radiation pressure qdp = 0, so qdp(tot) = qdp(gas)*p/pgas
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = idtt+i-1
         temp_work(i,1) = scv_weighted_sum4(press_interp_weights, tablex(ii,idp:idp+3,iscv_dlnrho_dlnp))
         temp_work(i,2) = scv_weighted_sum4(press_interp_weights, tabley(ii,idp:idp+3,iscv_dlnrho_dlnp))
         temp_work(i,3) = scv_weighted_sum4(press_interp_weights, tablez(ii,idp:idp+3,iscvz_dlnrho_dlnp))
      end do
! interpolate in temperature
      dlnrho_dlnp_pure_h = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,1))
      dlnrho_dlnp_pure_he = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,2))
      dlnrho_dlnp_pure_z = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,3))
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
         temp_work(i,1) = scv_weighted_sum4(press_interp_weights, tablex(ii,idp:idp+3,iscv_dlnrho_dlnt))
         temp_work(i,2) = scv_weighted_sum4(press_interp_weights, tabley(ii,idp:idp+3,iscv_dlnrho_dlnt))
         temp_work(i,3) = scv_weighted_sum4(press_interp_weights, tablez(ii,idp:idp+3,iscvz_dlnrho_dlnt))
      end do
! interpolate in temperature
      dlnrho_dlnt_pure_h = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,1))
      dlnrho_dlnt_pure_he = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,2))
      dlnrho_dlnt_pure_z = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,3))
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
         temp_work(i,1) = scv_weighted_sum4(press_interp_weights, tablex(ii,idp:idp+3,iscv_log10_s))
         temp_work(i,2) = scv_weighted_sum4(press_interp_weights, tabley(ii,idp:idp+3,iscv_log10_s))
         temp_work(i,3) = scv_weighted_sum4(press_interp_weights, smix(ii,idp:idp+3))
      end do
! interpolate in temperature
      log10_entropy_pure_h = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,1))
      entropy_pure_h = exp(ln10*log10_entropy_pure_h)
      log10_entropy_pure_he = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,2))
      entropy_pure_he = exp(ln10*log10_entropy_pure_he)
      entropy_of_mixing = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,3))
      dlnsmix_dlnt = temp_interp_weight_derivs(1)*log10(temp_work(1,3)) + &
           temp_interp_weight_derivs(2)*log10(temp_work(2,3)) &
           + temp_interp_weight_derivs(3)*log10(temp_work(3,3)) &
           + temp_interp_weight_derivs(4)*log10(temp_work(4,3))
! d ln s/ d ln t (x and y) and du/dt (z)
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = idtt+i-1
         temp_work(i,1) = scv_weighted_sum4(press_interp_weights, tablex(ii,idp:idp+3,iscv_dlns_dlnt))
         temp_work(i,2) = scv_weighted_sum4(press_interp_weights, tabley(ii,idp:idp+3,iscv_dlns_dlnt))
         temp_work(i,3) = scv_weighted_sum4(press_interp_weights, tablez(ii,idp:idp+3,iscvz_log10_du_dt))
      end do
! interpolate in temperature
      dlns_dlnt_pure_h = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,1))
      dlns_dlnt_pure_he = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,2))
      log10_du_dt_pure_z = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,3))
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

! Get fractions of total particles (including electrons), as follows:
!   XTF_H2  the fraction that is neutral hydrogen molecules, and
!   XTF_He  the fraction thst is neutral helium).
!   These are in column 2 respectively of the SCV hydrogen and
!   Helium tables.
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = idtt+i-1
         temp_work(i,1) = scv_weighted_sum4(press_interp_weights, tablex(ii,idp:idp+3,iscv_frac_neutral))
         temp_work(i,2) = scv_weighted_sum4(press_interp_weights, tabley(ii,idp:idp+3,iscv_frac_neutral))
      end do
! interpolate in temperature
      xtf_h2 = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,1))
      xtf_he = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,2))

! Get more fractions of total particles (including electrons), as follows:
!   XTF_H1  the fraction that is neutral hydrogen atoms, and
!   XTF_HeP  the fraction thst is singly ionized helium).
!   These are in column 3 respectively of the SCV hydrogen and
!   Helium tables.
! interpolate in pressure at 4 different temperature points.
      do i = 1,4
         ii = idtt+i-1
         temp_work(i,1) = scv_weighted_sum4(press_interp_weights, tablex(ii,idp:idp+3,iscv_frac_atom_ion1))
         temp_work(i,2) = scv_weighted_sum4(press_interp_weights, tabley(ii,idp:idp+3,iscv_frac_atom_ion1))
      end do
! interpolate in temperature
      xtf_h1 = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,1))
      xtf_hep = scv_weighted_sum4(temp_interp_weights, temp_work(1:4,2))

! At tis time we can calculate the ramaining particles fractions using
! conservation of particles and coservation of charge.  The variable
! names are:
!  XTF_H_e  the hycrogen related fraction that is electrons
!  XTF_HP   the hydrogen ralated fraction that is H+ ions
!  XTF_He_e the helium related fraction that is electrons
!  XTF_HePP the helium related fraction that is doubly ionized helium (He++)
! Particle and charge conservation yields:
      xtf_h_e = .5d0*(1d0 - xtf_h2 - xtf_h1)
      xtf_hp = xtf_h_e
      xtf_he_e = (1d0/3d0)*(2d0 - 2d0*xtf_he - xtf_hep)
      xtf_hepp = 1d0 - xtf_he - xtf_hep - xtf_he_e

! We are seeking the fraction XHP of hydrogen nuclei that are singly
! ionized, the fraction XHeP of helium huclei that are singly ionized
! and XHePP, the helium fraction that is doubly ionized.
      xhp = xtf_hp / (xtf_h2 + xtf_h1 + xtf_hp)
      xhep = xtf_hep / (xtf_he + xtf_hep + xtf_hepp)
      xhepp = xtf_hepp / (xtf_he + xtf_hep + xtf_hepp)

! In the same order, these are the desired elements of array ion_fraction(3)
      r%ion_fraction(1) = xhp
      r%ion_fraction(2) = xhep
      r%ion_fraction(3) = xhepp

      return    ! Primay exit. Success. valid_table_point is true.

end subroutine compute_scv_eos

!----------------------------------------------------------------------
! build_scv_envelope_table (formerly setscv)
!----------------------------------------------------------------------
! Builds the envelope-mixture (X/Y/Z) SCVH-style equation-of-state
! table tablenv from the separately tabulated pure-hydrogen (tablex),
! pure-helium (tabley), and pure-metal (tablez) tables: first
! completes du/dt and the entropy-of-mixing correction (smix) on the
! native H/He table grid, then combines the three species by the
! additive-volume rule at the envelope composition (X=envelope_hydrogen_
! fraction, Z=envelope_metal_fraction) to fill tablenv columns 1-6,
! then numerically differentiates those columns (in log T and log P)
! to fill tablenv columns 7-12.
subroutine build_scv_envelope_table(envelope_hydrogen_fraction, &
     envelope_metal_fraction)

      integer, parameter :: nts = scv_nt
! the envelope composition the table is built for
      real(dp), intent(in) :: envelope_hydrogen_fraction, &
           envelope_metal_fraction

! interp_x: the three abscissae handed to inter3; interp_f: the three
! tabulated values combined with the weights inter3 returns.
      real(dp) :: interp_x(3), interp_f(3), t_interp_weight(3), &
           t_interp_dweight(3), p_interp_weight(3), p_interp_dweight(3)
      real(dp), parameter :: hydrogen_atom_mass = 1.67357d-24
      real(dp), parameter :: helium_atom_mass = 6.646442d-24
      real(dp), parameter :: boltzmann_constant = 1.380658d-16
! --- locals ---
      integer :: t_idx, p_idx, k_idx, ii, jj
! log_t_work/log_p_work are reused as scratch across the file: log10
! scale in the first two loops below, natural-log scale (CLN*...) in
! the final derivative loop -- this mirrors the original TL/PL reuse.
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
      real(dp) :: dqdt_dlnp, dqdt_dlnt, dlncp_dlnt, &
           dqut_dlnt, dqut_dlnp, dlnp_dlnrho, ln_rho_value

!  Read in equation of state tables for hydrogen and helium
      do t_idx = 1, nts
         do p_idx = 1, nptsx(t_idx)
!  Compute du/dt - needed for specific heat derivatives
!  cp = T (ds/dt)|P = s(d ln s/d ln T)|P (eq 1)
! and cp = du/dt + P*(d ln rho/d ln T)**2/(rho*T*d ln rho/d ln P) (eq 2)
! du/dt can be inferred from eq 1 and the second term in eq 2, all of
! which are contained in the tables.
            tablex(t_idx,p_idx,iscv_du_dt) = &
                 exp(ln10*tablex(t_idx,p_idx,iscv_log10_s))*tablex(t_idx,p_idx,iscv_dlns_dlnt) - &
                 exp(ln10*(tablex(t_idx,p_idx,iscv_log10_p)-tablex(t_idx,p_idx,iscv_log10_rho)- &
                 tlogx(t_idx)))*tablex(t_idx,p_idx,iscv_dlnrho_dlnt)**2/tablex(t_idx,p_idx,iscv_dlnrho_dlnp)
            tabley(t_idx,p_idx,iscv_du_dt) = &
                 exp(ln10*tabley(t_idx,p_idx,iscv_log10_s))*tabley(t_idx,p_idx,iscv_dlns_dlnt) - &
                 exp(ln10*(tabley(t_idx,p_idx,iscv_log10_p)-tabley(t_idx,p_idx,iscv_log10_rho)- &
                 tlogx(t_idx)))*tabley(t_idx,p_idx,iscv_dlnrho_dlnt)**2/tabley(t_idx,p_idx,iscv_dlnrho_dlnp)
! Compute the entropy of mixing.
! Number density of hydrogen and helium
            helium_fraction_hhe = 1.0d0-envelope_hydrogen_fraction
            one_minus_y_local = 1.0d0 - helium_fraction_hhe
            density_h = exp(ln10*tablex(t_idx,p_idx,iscv_log10_rho))
            density_he = exp(ln10*tabley(t_idx,p_idx,iscv_log10_rho))
            density_mix = 1.0d0/(one_minus_y_local/density_h + &
                 helium_fraction_hhe/density_he)
            hydrogen_number_density = 2.0d0*one_minus_y_local*density_mix/ &
                 hydrogen_atom_mass/(1.0d0+3.0d0*tablex(t_idx,p_idx,iscv_frac_neutral)+ &
                 tablex(t_idx,p_idx,iscv_frac_atom_ion1))
            helium_number_density = 3.0d0*helium_fraction_hhe*density_mix/ &
                 helium_atom_mass/(1.0d0 + 2.0d0*tabley(t_idx,p_idx,iscv_frac_neutral) + &
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
               mixing_beta = (hydrogen_atom_mass/helium_atom_mass)* &
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
               smix(t_idx,p_idx) = one_minus_y_local/hydrogen_atom_mass*2.0d0/ &
                    (1.0d0+tablex(t_idx,p_idx,iscv_frac_atom_ion1)+3.0d0*tablex(t_idx,p_idx,iscv_frac_neutral))* &
                    (log(1.0d0+mixing_beta*mixing_gamma)- &
                    hydrogen_electron_density/total_number_density* &
                    log(1.0d0 + mix_de)+mixing_beta*mixing_gamma* &
                    log(1.0d0+1.0d0/mixing_beta/mixing_gamma)- &
                    helium_electron_density/total_number_density* &
                    log(1.0d0+mix_dei))
               smix(t_idx,p_idx) = boltzmann_constant*smix(t_idx,p_idx)
            end if
         end do
      end do
!  Now compute equation of state variables for the surface mixture.
      do t_idx=1,nts
         log_t_work = tlogx(t_idx)
         temp_value = (exp10(log_t_work))
! Temperature interpolation factors
        if (t_idx == 1) then
           idtt = 1
        else if (t_idx == nts) then
           idtt = t_idx - 2
        else
           idtt = t_idx - 1
        end if
        do k_idx = 1,3
           interp_x(k_idx) = ln10*tlogx(k_idx+idtt-1)
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
            dsmix_dlnt = t_interp_dweight(1)*smix(idtt,p_idx)+ &
                 t_interp_dweight(2)*smix(idtt+1,p_idx)+ &
                 t_interp_dweight(3)*smix(idtt+2,p_idx)
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
! Now find the following derivatives numerically :
! D QDT/ D LN T
! D QDT/ D LN P
! D LN CP / D LN P
! D LN CP / D LN T
! D LN DEL(AD) / D LN P
! D LN DEL(AD) / D LN T
      do t_idx = 1, nts
! Temperature interpolation factors
        if (t_idx == 1) then
           idtt = 1
        else if (t_idx == nts) then
           idtt = t_idx - 2
        else
           idtt = t_idx - 1
        end if
        do k_idx = 1,3
           interp_x(k_idx) = ln10*tlogx(k_idx+idtt-1)
        end do
        log_t_work = ln10*tlogx(t_idx)
        call inter3(interp_x,t_interp_weight,t_interp_dweight,log_t_work)
        temp_value = exp(log_t_work)
        do p_idx = 1, nptsx(t_idx)
! Pressure interpolation factors
           if (p_idx == 1) then
              idp = 1
           else if (p_idx == nptsx(t_idx)) then
              idp = p_idx - 2
           else
              idp = p_idx - 1
           end if
           do k_idx = 1,3
              interp_x(k_idx) = ln10*tablenv(t_idx,idp+k_idx-1,iscvenv_log10_p)
           end do
           log_p_work = ln10*tablenv(t_idx,p_idx,iscvenv_log10_p)
           call inter3(interp_x,p_interp_weight,p_interp_dweight,log_p_work)
! Derivatives of d ln rho/ d ln T
           dqdt_dlnp = tablenv(t_idx,idp,iscvenv_dlnrho_dlnt)*p_interp_dweight(1)+ &
                tablenv(t_idx,idp+1,iscvenv_dlnrho_dlnt)*p_interp_dweight(2) &
                +tablenv(t_idx,idp+2,iscvenv_dlnrho_dlnt)*p_interp_dweight(3)
           dqdt_dlnt = tablenv(idtt,p_idx,iscvenv_dlnrho_dlnt)*t_interp_dweight(1)+ &
                tablenv(idtt+1,p_idx,iscvenv_dlnrho_dlnt)*t_interp_dweight(2) &
                +tablenv(idtt+2,p_idx,iscvenv_dlnrho_dlnt)*t_interp_dweight(3)
! DERIVATIVES OF LN CP (only the T derivative is stored; the P
! derivative was never used)
           do k_idx = 1,3
              ii = idtt+k_idx-1
            jj = min(nptsx(ii),p_idx)
            interp_f(k_idx) = log(tablenv(ii,jj,iscvenv_cp))
           end do
           dlncp_dlnt = interp_f(1)*t_interp_dweight(1)+ &
                interp_f(2)*t_interp_dweight(2)+interp_f(3)*t_interp_dweight(3)
! Derivatives of du/dt
           do k_idx = 1,3
              ii = idtt+k_idx-1
            jj = min(nptsx(ii),p_idx)
            interp_f(k_idx) = tablenv(ii,jj,iscvenv_du_dt)
              interp_f(k_idx) = tablenv(idtt+k_idx-1,p_idx,iscvenv_du_dt)
           end do
           dqut_dlnt = interp_f(1)*t_interp_dweight(1)+ &
                interp_f(2)*t_interp_dweight(2)+interp_f(3)*t_interp_dweight(3)
           do k_idx = 1,3
              interp_f(k_idx) = tablenv(t_idx,idp+k_idx-1,iscvenv_du_dt)
           end do
           dqut_dlnp = interp_f(1)*p_interp_dweight(1)+ &
                interp_f(2)*p_interp_dweight(2)+interp_f(3)*p_interp_dweight(3)
! D/d ln rho (d ln P/d ln rho)
           do k_idx = 1,3
              interp_x(k_idx) = ln10*tablenv(t_idx,idp+k_idx-1,iscvenv_log10_rho)
           end do
           ln_rho_value = ln10*tablenv(t_idx,p_idx,iscvenv_log10_rho)
           call inter3(interp_x,p_interp_weight,p_interp_dweight,ln_rho_value)
           dlnp_dlnrho = p_interp_dweight(1)/tablenv(t_idx,idp,iscvenv_dlnrho_dlnp)+ &
                p_interp_dweight(2)/tablenv(t_idx,idp+1,iscvenv_dlnrho_dlnp) &
                +p_interp_dweight(3)/tablenv(t_idx,idp+2,iscvenv_dlnrho_dlnp)
           tablenv(t_idx,p_idx,iscvenv_dqdt_dlnp) = dqdt_dlnp
           tablenv(t_idx,p_idx,iscvenv_dqdt_dlnt) = dqdt_dlnt
           tablenv(t_idx,p_idx,iscvenv_dlnp_dlnrho) = dlnp_dlnrho
           tablenv(t_idx,p_idx,iscvenv_dlncp_dlnt) = dlncp_dlnt
           tablenv(t_idx,p_idx,iscvenv_dqut_dlnp) = dqut_dlnp
           tablenv(t_idx,p_idx,iscvenv_dqut_dlnt) = dqut_dlnt
        end do
      end do

end subroutine build_scv_envelope_table

!----------------------------------------------------------------------
! read_scv_tables
!----------------------------------------------------------------------
! Reads the SCV tables (columns 1-11 of the pure-H and pure-He tables,
! 1-13 of the pure-metal one, computed with the Prather EOS) into
! tablex, tabley and tablez. The He and Z tables share the H table's
! T grid (tlogx) and P(T) point counts (nptsx).
subroutine read_scv_tables(scv_h_table_path, scv_he_table_path, scv_z_table_path, &
     ierr)
      character(len=*), intent(in) :: scv_h_table_path, scv_he_table_path, &
           scv_z_table_path
      integer, intent(out) :: ierr
      integer :: scv_h_unit, scv_he_unit, scv_z_unit
      integer :: t_idx, p_idx, col_idx
      real(dp) :: scvhe_dummy_val, scvz_dummy_val
      integer :: scvhe_dummy_npts, scvz_dummy_npts

      call open_scv_table(scv_h_table_path, scv_h_unit, ierr)
      if (ierr /= 0) return
      call open_scv_table(scv_he_table_path, scv_he_unit, ierr)
      if (ierr /= 0) return
      call open_scv_table(scv_z_table_path, scv_z_unit, ierr)
      if (ierr /= 0) return
!  Read in equation of state tables for hydrogen and helium
      do t_idx = 1, scv_nt
         read(scv_h_unit,1) tlogx(t_idx),nptsx(t_idx)
         read(scv_he_unit,1) scvhe_dummy_val,scvhe_dummy_npts
         read(scv_z_unit,1) scvz_dummy_val,scvz_dummy_npts
    1       format(f5.2,i4)
! The He and Z tables share the H table's T grid and P(T) point
! counts; their copies are read into dummies.
         do p_idx = 1, nptsx(t_idx)
            read(scv_h_unit,2) (tablex(t_idx,p_idx,col_idx),col_idx=1,11)
            read(scv_he_unit,2) (tabley(t_idx,p_idx,col_idx),col_idx=1,11)
    2          format(f6.2,1p2e13.5,0p,8f9.4)
         end do
!  Read in metal equation of state table; computed using the Prather
! equation of state in the old Yale code.
         do p_idx = nptsx(t_idx),1,-1
            read(scv_z_unit,3)(tablez(t_idx,p_idx,col_idx),col_idx=1,13)
 3             format(f6.2,12f9.4)
         end do
      end do
      close(scv_h_unit)
      close(scv_he_unit)
      close(scv_z_unit)
end subroutine read_scv_tables

subroutine open_scv_table(path, table_unit, ierr)
      character(len=*), intent(in) :: path
      integer, intent(out) :: table_unit, ierr
      open(newunit=table_unit, file=path, status='OLD', iostat=ierr)
      if (ierr /= 0) then
         write(run_log_unit,*) 'eos_init: cannot open the SCV table ', trim(path)
         ierr = 1
      end if
end subroutine open_scv_table

end module scv_eos
