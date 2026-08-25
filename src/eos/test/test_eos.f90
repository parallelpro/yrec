!----------------------------------------------------------------------
! test_eos
!----------------------------------------------------------------------
! Standalone eos-domain test program (2026, phase three -- ROADMAP.md
! stage 2), following MESA's per-module test convention: boot only
! what the eos domain needs, evaluate the public facade entries over
! a fixed grid of points, and print the results in a stable format
! for byte-comparison against the checked-in expected output
! (expected_test_eos.out, regenerated deliberately and reviewed like
! any other diff).
!
! Initialization notes -- everything here replicates what the full
! code does at startup, with the source of each piece cited:
!  * unit numbers: core/parmin.f90's assignments (short_file_unit=20,
!    fermi_unit=15, scv units 72-74, iopale=49).
!  * flags/cutoffs: explicit assignments below (TSCUT=6.0 matching
!    the reference solar-model namelists; all kap table flags false
!    so setups' kap_init is a no-op; atm_choice=0 so atm_init is a
!    no-op).
!  * constants + eos table loads: the real setup/setups.f90 +
!    eos_lib's eos_init, with real Fermi/SCV table paths and dummy
!    paths for everything gated off.
!  * env_comp (the envelope-composition state eqstat2 consumes):
!    core/starin.f90's exact mixture algorithm (its statements at the
!    "COMPUTE SURFACE MIX VALUES" block), seeded from const_lib's
!    vnew defaults (the G&N93 12-species mixture) and starin's
!    atomic_weight data, for X=0.70, Z=0.016492 (matching the shipped
!    EOSOPAL06Z0.016492 table).
!  * the OPAL-2006 EOS file is opened on iopale exactly as
!    core/parmin.f90 does; the table itself is lazily read by esac06
!    on first use, as in production.
!
! The MHD path is exercised only if YREC_MHD_TABLES is set to a
! directory holding the 8 MHD table files -- the repository ships no
! MHD tables, so by default that section reports SKIPPED. This is
! the infrastructure half of closing ROADMAP.md's LMHD coverage gap;
! the data half needs tables from outside the repo.
program test_eos
      use const_lib
      use luout_lib
      use star_info_lib, only: star
      use eos_lib
      use opacity_table_lib
      use yale_eos_lib
      use scv_eos_lib
      implicit none

      character(len=256) :: yrec_input, mhd_dir
      character(len=256) :: fermi_path, scvh_path, scvhe_path, scvz_path, &
           opal06_path, dummy_path
      character(len=256) :: dummy_paths7(7)
      double precision :: laol_work(12)

! starin.f90's atomic_weight data (12-species set, same order as vnew)
      double precision :: atomic_weight(12)
      data atomic_weight /23.0d0,26.99d0,24.32d0,55.86d0,28.1d0,12.015d0, &
           1.008d0,16.0d0,14.01d0,39.96d0,20.19d0,4.004d0/

      double precision :: w(12), wsum, scale
      integer :: setups_ierr, i, ipt, eos_ierr

! eos_get argument set
      double precision :: logt, t, logp, p, logd, d
      double precision :: x_frac, z_frac
      double precision :: beta, betai, beta14, fxion(3), rmu, amu, emu, eta
      double precision :: qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp
      logical :: lderiv, latmo
      integer :: ksaha

! eos_get_gamma1 argument set
      double precision :: t6, rho, pres, gamma1, grad_ad

! test point grids
      integer, parameter :: npts = 8
      double precision :: grid_logt(npts), grid_logp(npts)
      data grid_logt /3.70d0, 3.90d0, 4.50d0, 5.50d0, 6.30d0, 6.90d0, &
           7.10d0, 7.161d0/
      data grid_logp /4.5d0, 6.0d0, 8.0d0, 11.5d0, 14.0d0, 16.0d0, &
           16.9d0, 17.35d0/

      integer, parameter :: ng1 = 4
      double precision :: g1_t6(ng1), g1_rho(ng1), g1_p(ng1)
      data g1_t6  /0.006d0, 0.10d0,  2.0d0,   12.0d0/
      data g1_rho /1.0d-7,  1.0d-4,  1.0d-1,  80.0d0/
      data g1_p   /3.16d4,  1.0d8,   3.16d13, 1.6d17/

      call get_environment_variable("YREC_INPUT", yrec_input)
      if (len_trim(yrec_input) == 0) then
         write(*,'(a)') "test_eos: FAIL (YREC_INPUT not set)"
         stop 1
      end if
      fermi_path  = trim(yrec_input)//"/eos/yale/FERMI.TAB"
      scvh_path   = trim(yrec_input)//"/eos/scv/h_tab_i.dat"
      scvhe_path  = trim(yrec_input)//"/eos/scv/he_tab_i.dat"
      scvz_path   = trim(yrec_input)//"/eos/scv/z_tab_i.dat"
      opal06_path = trim(yrec_input)//"/eos/opal2006/EOSOPAL06Z0.016492"
      dummy_path = ""
      dummy_paths7 = ""
      laol_work = 0.0d0

! unit numbers, per core/parmin.f90
      short_file_unit = 20
      star%ctrl%fermi_unit = 15
      star%ctrl%scv_h_unit = 72
      star%ctrl%scv_he_unit = 73
      star%ctrl%scv_z_unit = 74
      star%ctrl%iopale = 49
      open(short_file_unit, file="test_eos.short", status="replace")

! eos configuration
      star%ctrl%use_mhd_eos = .false.
      use_scv_eos = .true.
      star%ctrl%use_opal95_eos = .false.
      star%ctrl%use_opal2001_eos = .false.
      star%ctrl%use_opal2006_eos = .true.
      star%job%use_diffusion_z = .false.
      star%ctrl%use_numerical_derivatives = .false.
      use_debye_huckel_correction = .false.
      star%ctrl%saha_log10t_cutoff = 6.0d0

! kap/atm gated off so setups' kap_init/atm_init are no-ops
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

! envelope composition state, per starin.f90's mixture algorithm
      star%xnew = 0.70d0
      star%znew = 0.016492d0
      star%envelope_hydrogen_fraction = star%xnew
      star%envelope_metal_fraction = star%znew
      do i = 1, 12
         w(i) = star%ctrl%vnew(i)
      end do
      wsum = w(1)+w(2)+w(3)+w(4)+w(5)+w(6)+w(8)+w(9)+w(10)+w(11)
      star%zenvm = star%envelope_metal_fraction* &
           (wsum - w(6) - w(8) - w(9))/wsum
      scale = star%envelope_metal_fraction/wsum
      w(7) = star%envelope_hydrogen_fraction/scale
      w(12) = (1.0d0 - star%envelope_hydrogen_fraction - &
           star%envelope_metal_fraction)/scale
      wsum = 0.0d0
      do i = 1, 12
         w(i) = scale*w(i)/atomic_weight(i)
         wsum = wsum + w(i)
      end do
      star%amuenv = wsum
      scale = 1.0d0/star%amuenv
      do i = 1, 12
         star%fxenv(i) = w(i)*scale
      end do
! push the mixture to the eos domain, as starin does (physics-purity)
      call eos_set_mixture(star%envelope_hydrogen_fraction, &
           star%envelope_metal_fraction, star%amuenv, &
           star%fxenv)

! constants + Fermi/SCV table loads (real setups + eos_init inside it)
      call setups(laol_work, dummy_path, dummy_path, dummy_path, &
           fermi_path, dummy_path, dummy_path, dummy_path, dummy_path, &
           dummy_path, dummy_path, dummy_path, dummy_path, dummy_path, &
           dummy_path, dummy_path, dummy_path, dummy_path, dummy_path, &
           dummy_path, dummy_path, scvh_path, scvhe_path, scvz_path, &
           dummy_paths7, setups_ierr)
      if (setups_ierr /= 0) then
         write(*,'(a)') "test_eos: FAIL (setups error)"
         stop 1
      end if

! OPAL-2006 EOS file, opened as core/parmin.f90 does (esac06 reads it
! lazily on first use)
      open(star%ctrl%iopale, file=opal06_path, status='OLD')

      x_frac = 0.70d0
      z_frac = 0.016492d0

      write(*,'(a)') "# test_eos: eos_get over (logT, logP) grid, " // &
           "X=0.70 Z=0.016492, SCV+OPAL2006, no MHD, no DH"
      do ipt = 1, npts
         logt = grid_logt(ipt)
         logp = grid_logp(ipt)
         logd = 0.0d0
         d = 0.0d0
         ksaha = 0
         lderiv = .true.
         latmo = .false.
         fxion = 0.0d0
         call eos_get(logt, t, logp, p, logd, d, x_frac, z_frac, &
              beta, betai, beta14, fxion, rmu, amu, emu, eta, &
              qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp, &
              lderiv, latmo, ksaha)
         write(*,'(a,i2)') "point ", ipt
         write(*,'(4(1pe24.15))') logt, logp, logd, beta
         write(*,'(4(1pe24.15))') qdt, qdp, qcp, dela
         write(*,'(4(1pe24.15))') rmu, amu, emu, eta
         write(*,'(3(1pe24.15))') fxion(1), fxion(2), fxion(3)
         write(*,'(4(1pe24.15))') qdtt, qdtp, qat, qap
         write(*,'(2(1pe24.15),i4)') qcpt, qcpp, ksaha
      end do

      write(*,'(a)') "# test_eos: eos_get_gamma1, OPAL2006 branch"
      do ipt = 1, ng1
         t6 = g1_t6(ipt)
         rho = g1_rho(ipt)
         pres = g1_p(ipt)
         ksaha = 0
         call eos_get_gamma1(x_frac, z_frac, t6, rho, pres, gamma1, &
              grad_ad, ksaha)
         write(*,'(a,i2,2(1pe24.15))') "g1 ", ipt, gamma1, grad_ad
      end do

! Only the two hotter points here: at the two low-T points the
! Yale/SCV branch's Cox & Guili gamma1 reconstruction consumes
! quantities eqstat2 leaves undefined under the branch's historical
! lderiv=.false. setting (inherited verbatim from calcad's original
! formulation; NaN/uninitialized in practice) -- a pre-existing
! property of that code path, unreachable in production configs
! (which set LOPALE06), and not meaningfully pinnable.
      write(*,'(a)') "# test_eos: eos_get_gamma1, Yale/SCV branch " // &
           "(defined regime only)"
      star%ctrl%use_opal2006_eos = .false.
      do ipt = 3, ng1
         t6 = g1_t6(ipt)
         rho = g1_rho(ipt)
         pres = g1_p(ipt)
         ksaha = 0
         call eos_get_gamma1(x_frac, z_frac, t6, rho, pres, gamma1, &
              grad_ad, ksaha)
         write(*,'(a,i2,2(1pe24.15))') "g1 ", ipt, gamma1, grad_ad
      end do
      star%ctrl%use_opal2006_eos = .true.

! MHD path: infrastructure for closing the LMHD coverage gap; runs
! only when tables are supplied from outside the repository.
      call get_environment_variable("YREC_MHD_TABLES", mhd_dir)
      if (len_trim(mhd_dir) == 0) then
         write(*,'(a)') "# test_eos: MHD path SKIPPED " // &
              "(YREC_MHD_TABLES not set; no MHD tables ship with YREC)"
      else
         star%ctrl%use_mhd_eos = .true.
         call eos_init(fermi_path, scvh_path, scvhe_path, scvz_path, &
              trim(mhd_dir)//"/zams_a", trim(mhd_dir)//"/zams_b", &
              trim(mhd_dir)//"/zams_c", trim(mhd_dir)//"/centre1", &
              trim(mhd_dir)//"/centre2", trim(mhd_dir)//"/centre3", &
              trim(mhd_dir)//"/centre4", trim(mhd_dir)//"/centre5")
         write(*,'(a)') "# test_eos: eos_get over the same grid, MHD"
         do ipt = 1, npts
            logt = grid_logt(ipt)
            logp = grid_logp(ipt)
            logd = 0.0d0
            ksaha = 0
            lderiv = .true.
            latmo = .false.
            fxion = 0.0d0
            call eos_get(logt, t, logp, p, logd, d, x_frac, z_frac, &
                 beta, betai, beta14, fxion, rmu, amu, emu, eta, &
                 qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp, &
                 lderiv, latmo, ksaha)
            write(*,'(a,i2)') "mhd point ", ipt
            write(*,'(4(1pe24.15))') logt, logp, logd, beta
            write(*,'(4(1pe24.15))') qdt, qdp, qcp, dela
         end do
         star%ctrl%use_mhd_eos = .false.
      end if

! Error paths (2026, ROADMAP.md stage 3): in the eos domain the error
! conditions behind the converted stops are corrupt-table /
! misconfiguration cases (out-of-table lookups stay *recoverable* by
! design here -- alternate returns and sentinels -- unlike kap), so
! with healthy tables the facade cannot be driven to ierr /= 0. The
! first check asserts the success path threads ierr = 0 end to end
! through eqstat/eqstat2/oeqos06/esac06. The second drives esac06 (a
! domain internal -- this test is in-domain, white-box by design)
! with an invalid rad_flag, the same class of internal error the
! converted stops guarded, and asserts ierr = 1 with no crash; its
! diagnostic goes to short_file_unit as always.
      write(*,'(a)') "# test_eos: error paths via ierr"
      logt = 6.30d0
      logp = 14.0d0
      logd = 0.0d0
      d = 0.0d0
      ksaha = 0
      lderiv = .true.
      latmo = .false.
      fxion = 0.0d0
      call eos_get(logt, t, logp, p, logd, d, x_frac, z_frac, &
           beta, betai, beta14, fxion, rmu, amu, emu, eta, &
           qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp, &
           lderiv, latmo, ksaha, ierr=eos_ierr)
      write(*,'(a,i4)') "err facade-success ierr = ", eos_ierr
      call esac06(0.70d0, 2.0d0, 0.1d0, 9, 2, eos_ierr, *300)
  300 continue
      write(*,'(a,i4)') "err bad-rad-flag   ierr = ", eos_ierr

      close(short_file_unit)
      write(*,'(a)') "test_eos: done"
end program test_eos
