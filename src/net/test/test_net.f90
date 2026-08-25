!----------------------------------------------------------------------
! test_net
!----------------------------------------------------------------------
! Standalone net-domain test program (2026), companion to
! eos/kap/atm's (see eos/test/test_eos.f90 for the initialization
! conventions). Exercises the composition-in/rates-out core of the
! domain -- the pieces that are pure functions of (logT, logRho,
! composition) once the cross-section scales are set:
!   * setup/remap.f90's cross-section scale computation (from the
!     newcross namelist defaults),
!   * rates (the 13 reaction rates + branching fractions),
!   * sneut (Itoh et al. 1996 neutrino losses),
!   * azbar (abar/zbar bookkeeping),
!   * deutrate (deuterium burning, results in star%light_burn).
! engeb is deliberately NOT exercised: it is the per-zone energy
! generation driver coupled to the full model state (star%prev, flux
! diagnostics, timestep) -- a full-model concern covered by Stage-0.
! Results print for byte-comparison against expected_test_net.out.
program test_net
      use const_lib
      use luout_lib
      use star_info_lib, only: star, i_h2
      use net_lib
      use scv_eos_lib, only: use_scv_eos
      use opacity_table_lib, only: use_pure_z_table
      use burn_lib
      implicit none

      character(len=256) :: yrec_input
      character(len=256) :: fermi_path, dummy_path
      character(len=256) :: dummy_paths7(7)
      double precision :: laol_work(12)
      integer :: setups_ierr, ipt, i

      integer, parameter :: npts = 4
      double precision :: grid_logt(npts), grid_logd(npts)
      data grid_logt /6.90d0, 7.10d0, 7.30d0, 7.60d0/
      data grid_logd /1.0d0,  1.8d0,  2.2d0,  3.0d0/

      integer, parameter :: json = 5000
      double precision, dimension(json) :: rpp, r33, r34, rc12, rc13, &
           rn14, ro16, rc13a, rz9, rc12a, rn14a, r3a, rz13, fc12a, fbe7e
      double precision :: snu, dsnudt, dsnudd, dsnuda, dsnudz
      double precision :: xmass(3), aion(3), zion(3), ymass(3), &
           abar, zbar

      call get_environment_variable("YREC_INPUT", yrec_input)
      if (len_trim(yrec_input) == 0) then
         write(*,'(a)') "test_net: FAIL (YREC_INPUT not set)"
         stop 1
      end if
      fermi_path = trim(yrec_input)//"/eos/yale/FERMI.TAB"
      dummy_path = ""
      dummy_paths7 = ""
      laol_work = 0.0d0

! unit numbers, per core/parmin.f90
      short_file_unit = 20
      star%ctrl%fermi_unit = 15
      open(short_file_unit, file="test_net.short", status="replace")

! everything gated off; only the Fermi table is a hard requirement
      star%ctrl%use_mhd_eos = .false.
      use_scv_eos = .false.
      star%ctrl%use_opal95_eos = .false.
      star%ctrl%use_opal2001_eos = .false.
      star%ctrl%use_opal2006_eos = .false.
      star%ctrl%use_opal95_tables = .false.
      star%ctrl%use_opal92_tables = .false.
      star%ctrl%use_laol89_tables = .false.
      star%ctrl%use_alex06_tables = .false.
      star%ctrl%use_alex95_tables = .false.
      star%ctrl%use_kurucz90_tables = .false.
      star%use_two_z_tables = .false.
      use_pure_z_table = .false.
      star%ctrl%use_conductive_opacity = .false.
      star%job%atm_choice = 0

      call setups(laol_work, dummy_path, dummy_path, dummy_path, &
           fermi_path, dummy_path, dummy_path, dummy_path, dummy_path, &
           dummy_path, dummy_path, dummy_path, dummy_path, dummy_path, &
           dummy_path, dummy_path, dummy_path, dummy_path, dummy_path, &
           dummy_path, dummy_path, dummy_path, dummy_path, dummy_path, &
           dummy_paths7, setups_ierr)
      if (setups_ierr /= 0) then
         write(*,'(a)') "test_net: FAIL (setups error)"
         stop 1
      end if

! Nuclear cross-sections: parmin keeps the S-factor namelist locals
! with Solar Fusion II defaults (Adelberger et al. 2011; see
! core/parmin.f90's DATA statements around line 1177) and
! copy-assigns them into controls. This test bypasses parmin, so it
! sets the same SFII values explicitly, with the new-rates path on.
      use_new_nuclear_rates = .true.
      s0_pp = 4.01d-22
      s0_he3he3 = 5.21d3
      s0_he3he4 = 5.6d-1
      s0_p_c12 = 1.34d0
      s0_p_c13 = 7.6d0
      s0_p_n14 = 1.66d0
      s0_p_o16 = 1.06d1
      s0_be7_electron = 1.7709d-10
      s0_be7_p = 0.0208d0
      s0_n15_p_c12_branch = 7.3d4
      s0_n15_p_o16_branch = 3.6d1
      s0p_pp = 4.49d-24
      s0p_he3he3 = -4.9d0
      s0p_he3he4 = -3.6d-4
      s0p_p_c12 = 2.6d-3
      s0p_p_c13 = -7.83d-3
      s0p_p_n14 = -3.3d-3
      s0p_p_o16 = -5.4d-2
      s0p_be7_p = -3.12d-5
      s0pp_p_c12 = 8.3d-5
      s0pp_p_c13 = 7.29d-4
      s0pp_p_o16 = 0.0d0
      s0pp_be7_p = -2.288d-7

! cross-section scales via the real remap
      call remap
      write(*,'(a)') "# test_net: cross_section_scale from remap " // &
           "(newcross defaults)"
      write(*,'(4(1pe20.12))') (star%cross_section_scale(i), i = 1, 16)

! reaction rates over a (logT, logRho) grid at a solar-ish mixture
      write(*,'(a)') "# test_net: rates, X=0.70 Y3=3e-5 C12=3.5e-3 " // &
           "C13=4e-5 N14=1e-3 O16=8e-3 O18=2e-5"
      do ipt = 1, npts
         call rates(grid_logd(ipt), grid_logt(ipt), 0.70d0, 0.28d0, &
              3.0d-5, 3.5d-3, 4.0d-5, 1.0d-3, 8.0d-3, 2.0d-5, ipt, &
              rpp, r33, r34, rc12, rc13, rn14, ro16, rc13a, rz9, &
              rc12a, rn14a, r3a, rz13, fc12a, fbe7e)
         write(*,'(a,i2)') "rates ", ipt
         write(*,'(4(1pe20.12))') rpp(ipt), r33(ipt), r34(ipt), rc12(ipt)
         write(*,'(4(1pe20.12))') rc13(ipt), rn14(ipt), ro16(ipt), rc13a(ipt)
         write(*,'(4(1pe20.12))') rc12a(ipt), rn14a(ipt), r3a(ipt), fbe7e(ipt)
      end do

! Itoh neutrino losses (pure; first call initializes its own tables)
      write(*,'(a)') "# test_net: sneut (Itoh et al. 1996)"
      do ipt = 1, npts
         call sneut(10.0d0**grid_logt(ipt), 10.0d0**grid_logd(ipt), &
              1.30d0, 1.10d0, snu, dsnudt, dsnudd, dsnuda, dsnudz)
         write(*,'(a,i2,3(1pe20.12))') "snu ", ipt, snu, dsnudt, dsnudd
      end do

! azbar bookkeeping on a 3-species H/He4/C12 mixture
      xmass = (/0.70d0, 0.28d0, 0.02d0/)
      aion  = (/1.0d0, 4.0d0, 12.0d0/)
      zion  = (/1.0d0, 2.0d0, 6.0d0/)
      call azbar(xmass, aion, zion, 3, ymass, abar, zbar)
      write(*,'(a)') "# test_net: azbar H/He4/C12"
      write(*,'(5(1pe20.12))') abar, zbar, ymass(1), ymass(2), ymass(3)

! deuterium burning rate (result lands in star%light_burn)
      write(*,'(a)') "# test_net: deutrate"
      star%xa(i_h2,1) = 2.0d-5
      call deutrate(1.5d0, 6.2d0, 0.70d0, 1, 1)
      write(*,'(1pe20.12)') star%light_burn%deuterium_burning_rate(1)

      close(short_file_unit)
      write(*,'(a)') "test_net: done"
end program test_net
