!----------------------------------------------------------------------
! opacity_blending
!----------------------------------------------------------------------
! The formulas kap_eval combines opacities with, each written once (pure
! kernels, star-blind):
!   blend_in_z            linear interpolation in Z between two tables'
!                         values (log10 kappa or a log derivative)
!   blend_results_in_z    a whole kap_result blended so: log10 kappa,
!                         then kappa = 10**(log10 kappa), then the two
!                         derivatives (the two-Z families and the
!                         helium-core blend toward the pure-Z table)
!   ramp_to_molecular     the ramp from the molecular opacity (below
!                         logt_min) to the interior one (above logt_max),
!                         linear in kappa, with the blend's derivatives
!   hubbard_lampe_applies, hubbard_lampe_opacity
!                         the conductive opacity of Hubbard & Lampe's
!                         fit, where the Potekhin table gives none
!   add_conduction        the radiative and conductive opacities added
!                         harmonically, 1/kappa = 1/kappa_rad + 1/kappa_cond
module opacity_blending
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_source, only: kap_result
      use physical_constants, only: ln10
      use math_lib
      implicit none
      private
      public :: blend_in_z, blend_results_in_z, ramp_to_molecular, &
           hubbard_lampe_applies, hubbard_lampe_opacity, add_conduction

! Hubbard & Lampe (1969)'s conductive opacity, as the F77 code's fit:
!   log10 kappa_c = log10(1 - 0.6 X) + c0 - (c_rho + c_rho2 log rho) log rho
!                   + (c_t + c_t2 log T) log T + c_rho_t log rho log T,
! used only where log T >= hubbard_lampe_logt_min and
! log rho >= 2 log T + hubbard_lampe_logrho_offset (degenerate enough
! for conduction to matter)
      real(dp), parameter :: hubbard_lampe_c0 = -14.6196d0, &
           hubbard_lampe_c_rho = 3.5853d0, hubbard_lampe_c_rho2 = 0.1386d0, &
           hubbard_lampe_c_t = 5.1324d0, hubbard_lampe_c_t2 = 0.3219d0, &
           hubbard_lampe_c_rho_t = 0.3901d0
      real(dp), parameter :: hubbard_lampe_logt_min = 4.2d0, &
           hubbard_lampe_logrho_offset = -13.0d0

contains

! value_z1, a quantity from the table at table_z1, replaced by its linear
! interpolation at metal_fraction between it and value_z2, the same
! quantity from the table at table_z2 (slope, then base + (Z - Z2) slope)
pure subroutine blend_in_z(value_z1, value_z2, metal_fraction, table_z1, table_z2)
      real(dp), intent(inout) :: value_z1
      real(dp), intent(in) :: value_z2, metal_fraction, table_z1, table_z2
      real(dp) :: slope

      slope = (value_z1 - value_z2) / &
           (table_z1 - table_z2)
      value_z1 = value_z2 + &
           (metal_fraction - table_z2)*slope
end subroutine blend_in_z

! res, from the table at table_z1, replaced by its blend at
! metal_fraction with res_z2 from the table at table_z2
pure subroutine blend_results_in_z(res, res_z2, metal_fraction, table_z1, table_z2)
      type(kap_result), intent(inout) :: res
      type(kap_result), intent(in) :: res_z2
      real(dp), intent(in) :: metal_fraction, table_z1, table_z2

      call blend_in_z(res%log10_opacity, res_z2%log10_opacity, &
           metal_fraction, table_z1, table_z2)
      res%opacity = exp10(res%log10_opacity)
      call blend_in_z(res%dlnkap_dlnrho, res_z2%dlnkap_dlnrho, &
           metal_fraction, table_z1, table_z2)
      call blend_in_z(res%dlnkap_dlnt, res_z2%dlnkap_dlnt, &
           metal_fraction, table_z1, table_z2)
end subroutine blend_results_in_z

! res, the interior opacity at log10_temperature, replaced by the ramp's
! value: kappa = w kappa_int + (1 - w) kappa_mol with w = (log T -
! logt_min) / (logt_max - logt_min) from logt_min up, and the molecular
! opacity below it. The log derivatives are kappa's own:
!   dln kappa/dln rho = (w kappa_int d_int + (1 - w) kappa_mol d_mol)/kappa
!   dln kappa/dln T   = (w kappa_int d_int + (1 - w) kappa_mol d_mol
!                        + (dw/dln T)(kappa_int - kappa_mol))/kappa,
! d_int and d_mol the interior and molecular derivatives, and
! dw/dln T = 1/(ln 10 (logt_max - logt_min)).
pure subroutine ramp_to_molecular(log10_temperature, logt_min, logt_max, molecular, res)
      real(dp), intent(in) :: log10_temperature, logt_min, logt_max
      type(kap_result), intent(in) :: molecular
      type(kap_result), intent(inout) :: res
      real(dp) :: ramp_weight, dweight_dlnt, interior_part, molecular_part
      type(kap_result) :: interior

      if (log10_temperature >= logt_min) then
         interior = res
         ramp_weight = (log10_temperature - logt_min) / (logt_max - logt_min)
         dweight_dlnt = 1.0d0 / (ln10*(logt_max - logt_min))
         interior_part = ramp_weight*interior%opacity
         molecular_part = (1.0d0 - ramp_weight)*molecular%opacity
         res%opacity = interior_part + molecular_part
         res%log10_opacity = log10(res%opacity)
         res%dlnkap_dlnrho = (interior_part*interior%dlnkap_dlnrho + &
              molecular_part*molecular%dlnkap_dlnrho) / res%opacity
         res%dlnkap_dlnt = (interior_part*interior%dlnkap_dlnt + &
              molecular_part*molecular%dlnkap_dlnt + &
              dweight_dlnt*(interior%opacity - molecular%opacity)) / res%opacity
      else
         res = molecular
      end if
end subroutine ramp_to_molecular

pure logical function hubbard_lampe_applies(log10_density, log10_temperature)
      real(dp), intent(in) :: log10_density, log10_temperature

      hubbard_lampe_applies = .not. (log10_temperature < hubbard_lampe_logt_min .or. &
           log10_density < (2.0d0*log10_temperature + hubbard_lampe_logrho_offset))
end function hubbard_lampe_applies

pure function hubbard_lampe_opacity(log10_density, log10_temperature, &
     hydrogen_fraction) result(res)
      real(dp), intent(in) :: log10_density, log10_temperature, hydrogen_fraction
      type(kap_result) :: res

      res%log10_opacity = &
           log10(1.0d0 - 0.6d0*hydrogen_fraction) + hubbard_lampe_c0 - &
           (hubbard_lampe_c_rho + hubbard_lampe_c_rho2*log10_density)*log10_density + &
           (hubbard_lampe_c_t - hubbard_lampe_c_t2*log10_temperature)*log10_temperature + &
           hubbard_lampe_c_rho_t*log10_density*log10_temperature
      res%opacity = exp10(res%log10_opacity)
      res%dlnkap_dlnrho = hubbard_lampe_c_rho_t*log10_temperature - &
           2.0d0*hubbard_lampe_c_rho2*log10_density - hubbard_lampe_c_rho
      res%dlnkap_dlnt = hubbard_lampe_c_rho_t*log10_density - &
           2.0d0*hubbard_lampe_c_t2*log10_temperature + hubbard_lampe_c_t
end function hubbard_lampe_opacity

! res, the radiative opacity, replaced by its harmonic sum with the
! conductive one, 1/kappa = 1/kappa_rad + 1/kappa_cond, and the sum's
! log derivatives
pure subroutine add_conduction(res, conductive)
      type(kap_result), intent(inout) :: res
      type(kap_result), intent(in) :: conductive
      type(kap_result) :: radiative

      radiative = res
      res%opacity = radiative%opacity*conductive%opacity / &
           (radiative%opacity + conductive%opacity)
      res%log10_opacity = log10(res%opacity)
      res%dlnkap_dlnrho = (radiative%dlnkap_dlnrho + conductive%dlnkap_dlnrho - &
           (radiative%opacity*radiative%dlnkap_dlnrho + &
           conductive%opacity*conductive%dlnkap_dlnrho) / &
           (radiative%opacity + conductive%opacity))
      res%dlnkap_dlnt = (radiative%dlnkap_dlnt + conductive%dlnkap_dlnt - &
           (radiative%opacity*radiative%dlnkap_dlnt + &
           conductive%opacity*conductive%dlnkap_dlnt) / &
           (radiative%opacity + conductive%opacity))
end subroutine add_conduction

end module opacity_blending
