!----------------------------------------------------------------------
! test_structure
!----------------------------------------------------------------------
! The structure equations' partials against finite differences.
! Loads a deck's starting model with the deck's own microphysics (the
! controls, the tables, the network and the convection settings, as
! run_yrec does before its first step), then at a handful of shells
! compares every partial evaluate_shell_equations returns with central
! differences of its own terms qp, qt, qr, ql in ln P, ln T, ln R and L.
! The change of log P and log T over the step moves with log P and
! log T, as in the Newton iteration, so the eps_grav partials are
! checked too (the step is set to 1e6 yr, with changes of 5e-4 in log P
! and 2e-4 in log T).
!
! Usage: test_structure DECK.inlist [TOLERANCE]   (run with YREC_INPUT,
! YREC_START and YREC_OUTPUT set, in a directory where the deck may write
! output). TOLERANCE (default 3e-3) bounds the relative difference of an
! entry marked ok.
! Prints, after the line 'test_structure report', per shell the 16
! entries of dq/dx: the analytic value, the finite difference (step 1e-4)
! and their relative difference |f - a|/max(|a|, |f|, 1e-6 |q|) at the
! best of the steps 1e-3, 1e-4, 1e-5, with a status:
!   ok     within the tolerance: 3e-3 with the OPAL tables, whose own
!          derivatives dln rho/dln P, dln rho/dln T agree with finite
!          differences of rho(P, T) to about 2e-3 (they enter qr, the
!          radiative qt and ql); 1e-4 with the Yale EOS, which is
!          accurate to that level and so also checks the eps_grav
!          partials where chi departs from -1 (the 1/chi of 18869257
!          moves them by 2e-3 there);
!   known  a documented approximation of the Jacobian:
!          (qt, L) in a convective shell: the convective gradient's
!            dependence on L is not linearized (qt_dl = 0);
!          (qt, ln P) and (qt, ln T) in a convective shell: dln grad/dln P
!            and dln grad/dln T from the convection routine, which near
!            the adiabat rest on the EOS's derivatives of grad_ad, known
!            to be inaccurate (docs/ROADMAP.md, the eos follow-up);
!   FAIL   anything else; the program then stops with status 1.
program test_structure
      use controls_reader, only: read_controls
      use run_setup, only: set_up_constants_and_tables
      use burn_step, only: configure_network
      use convection_lib, only: set_convection_settings
      use microphysics_settings, only: convection_settings_from_star
      use rotation_configuration, only: set_rotation_settings
      use rotation_settings_setup, only: rotation_settings_from_star
      use starting_model, only: read_starting_model
      use star_info_lib, only: star, control_nml_override, physics_nml_override
      use shell_equations, only: henyey_shell_terms, shell_state, &
           structure_settings, shell_diagnostics, evaluate_shell_equations
      use eos_lib, only: eos_result
      implicit none

      character(len=512) :: deck
      character(len=32) :: tol_arg
      integer :: ierr, nfail
      double precision :: tolerance
      type(structure_settings) :: cfg

      if (command_argument_count() < 1 .or. command_argument_count() > 2) then
         write(*,'(a)') 'usage: test_structure DECK.inlist [TOLERANCE]'
         stop 2
      end if
      call get_command_argument(1, deck)
      tolerance = 3.0d-3
      if (command_argument_count() == 2) then
         call get_command_argument(2, tol_arg)
         read(tol_arg,*) tolerance
      end if
      call load_deck_model(deck)
      call set_settings()
      nfail = 0
      call check_selected_shells()
      if (nfail > 0) then
         write(*,'(a,i4,a)') 'test_structure: ', nfail, ' FAIL'
         stop 1
      end if
      write(*,'(a)') 'test_structure: done'

contains

! ---------------------------------------------------------------
! The deck's controls, tables, network and convection settings, and its
! first card's starting model, as run_yrec sets them up.
subroutine load_deck_model(path)
      character(len=*), intent(in) :: path
      ierr = 0
      control_nml_override = path
      physics_nml_override = path
      call read_controls(ierr)
      if (ierr /= 0) call die('read_controls')
      call set_up_constants_and_tables(ierr)
      if (ierr /= 0) call die('set_up_constants_and_tables')
      call configure_network()
      star%job%nk = 1
      star%job%initial_envelope_x = star%job%initial_x_array(1)
      star%job%initial_envelope_z = star%job%initial_z_array(1)
      star%mixing_length_alpha = star%job%mixing_length_array(1)
      call set_convection_settings(convection_settings_from_star())
      call set_rotation_settings(rotation_settings_from_star())
      star%job%change_envelope_mass_flag = star%job%has_senv0_array(1)
      star%job%requested_envelope_mass = star%job%senv0_array(1)
      star%total_angular_momentum = 0.0d0
      star%total_rotational_ke = 0.0d0
      call read_starting_model(star%timestep_yr, star%dt, star%chosen_dt, star%trial_sign_flag, &
           star%ikut_flag, star%istore_flag, star%model_diverged_flag, &
           star%recompute_envelope_triangle, star%job%nk, star%dlnrho_dlnp, star%dlnrho_dlnt, &
           star%total_angular_momentum, star%total_rotational_ke, &
           star%convective_velocity, star%job%mixture_weights, ierr)
      if (ierr /= 0) call die('read_starting_model')
end subroutine load_deck_model

subroutine die(what)
      character(len=*), intent(in) :: what
      write(*,'(a)') 'test_structure: FAIL ('//what//')'
      stop 1
end subroutine die

! ---------------------------------------------------------------
! The shared settings, as compute_henyey_coefficients sets them, with a
! 1e6 yr step.
subroutine set_settings()
      cfg%dt = 1.0d6*3.1558d7
      cfg%burning_log_temperature_cutoff = star%ctrl%nuclear_logT_cutoffs(1)
      cfg%log_teff = star%log_Teff
      cfg%turnover_timescale = star%convective_turnover_timescale
      cfg%ln10_over_lsun = star%ln10_over_lsun
      cfg%rotation_active = star%job%rotation_active
end subroutine set_settings

! ---------------------------------------------------------------
! Shell im's state as compute_henyey_coefficients gathers it, with the
! test's changes over the step.
function shell_at(im) result(sh)
      integer, intent(in) :: im
      type(shell_state) :: sh
      sh%log_mass = star%log_mass(im)
      sh%mass = star%m(im)
      sh%shell_mass = star%dm(im)
      sh%log_temperature = star%logT(im)
      sh%log_pressure = star%logP(im)
      sh%log_radius = star%logR(im)
      sh%luminosity_lsun = star%luminosity_lsun(im)
      sh%log_density = star%logRho(im)
      sh%xa = star%xa(:,im)
      sh%fp = star%fp_rot(im)
      sh%ft = star%ft_rot(im)
      sh%log_temperature_delta = 2.0d-4
      sh%log_pressure_delta = 5.0d-4
      sh%rotational_ke = 0.0d0
      sh%rotational_ke_start = 0.0d0
      sh%in_overshoot = im >= star%iov1 .and. im <= star%iov2
      sh%in_surface_cz = im >= star%jcz
end function shell_at

! ---------------------------------------------------------------
! The four terms at a shell state, from a fresh copy of the eos record
! (every evaluation starts the eos from the same point).
subroutine evaluate(sh, eos0, q, terms, diag)
      type(shell_state), intent(in) :: sh
      type(eos_result), intent(in) :: eos0
      double precision, intent(out) :: q(4)
      type(henyey_shell_terms), intent(out) :: terms
      type(shell_diagnostics), intent(out) :: diag
      type(eos_result) :: eos
      integer :: saha_state
      eos = eos0
      saha_state = 0
      call evaluate_shell_equations(sh, cfg, eos, saha_state, terms, diag, ierr)
      if (ierr /= 0) call die('evaluate_shell_equations')
      q = [terms%qp, terms%qt, terms%qr, terms%ql]
end subroutine evaluate

! ---------------------------------------------------------------
! Pick the shells: the centre, the strongest burning, a radiative shell
! near the middle in mass, the base of the outer convection zone, the
! most superadiabatic shell, and the outer point.
subroutine check_selected_shells()
      type(shell_state) :: sh
      type(eos_result) :: eos0
      type(henyey_shell_terms) :: terms
      type(shell_diagnostics) :: diag
      double precision :: q(4), best_eps, mid, mid_rad, x
      integer :: im, i_burn, i_mid, i_rad, i_czbase, nz, k
      integer :: picked(6)
      character(len=24) :: labels(6)
      logical :: conv(5000)
      nz = star%nz
      i_burn = 0; i_mid = 1; i_rad = 0; i_czbase = 0
      best_eps = 0.0d0; mid = 1.0d99; mid_rad = 1.0d99
      do im = 1, nz
         sh = shell_at(im)
         eos0 = eos_result()
         eos0%log10_density = star%logRho(im)
         call evaluate(sh, eos0, q, terms, diag)
         conv(im) = diag%convective
         if (diag%burned) then
            if (diag%burn%total_energy_gen_rate > best_eps) then
               best_eps = diag%burn%total_energy_gen_rate; i_burn = im
            end if
         end if
         x = abs(star%m(im)/star%m(nz) - 0.5d0)
         if (x < mid) then
            mid = x; i_mid = im
         end if
         if (.not.diag%convective .and. x < mid_rad) then
            mid_rad = x; i_rad = im
         end if
      end do
      do im = nz, 1, -1
         if (.not.conv(im)) exit
         i_czbase = im
      end do
      if (i_rad > 0) i_mid = i_rad
      picked = [1, i_burn, i_mid, i_czbase, max(nz - 3, 1), nz]
      labels = [character(len=24) :: 'centre', 'peak burning', 'mid mass', &
           'convection zone base', 'outer convection zone', 'outer point']
      write(*,'(a)') 'test_structure report'
      write(*,'(a,1pe9.1)') 'tolerance', tolerance
      do k = 1, 6
         if (picked(k) <= 0) cycle
         if (any(picked(1:k-1) == picked(k))) cycle
         call check_shell(trim(labels(k)), picked(k))
      end do
end subroutine check_selected_shells

! ---------------------------------------------------------------
! The 16 partials at shell im: analytic (zero where the Jacobian has no
! entry) against central differences with steps h in ln x.
subroutine check_shell(label, im)
      character(len=*), intent(in) :: label
      integer, intent(in) :: im
      type(shell_state) :: sh0, sp, sm
      type(eos_result) :: eos0, eos_base
      type(henyey_shell_terms) :: t0, tdum
      type(shell_diagnostics) :: d0, ddum
      double precision :: q0(4), qp(4), qm(4), analytic(4,4), fd(4,4,3), h
      double precision, parameter :: ln10 = 2.302585092994046d0
      double precision, parameter :: steps(3) = [1.0d-3, 1.0d-4, 1.0d-5]
      character(len=2), parameter :: qname(4) = ['qp','qt','qr','ql']
      character(len=4), parameter :: xname(4) = ['lnP ','lnT ','lnR ','L   ']
      integer :: k, v, s
      double precision :: best
      character(len=8) :: status
      sh0 = shell_at(im)
      eos_base = eos_result()
      eos_base%log10_density = star%logRho(im)
! Every evaluation starts the eos from the same point.
      eos0 = eos_base
      call evaluate(sh0, eos0, q0, t0, d0)
      analytic(:,1) = [t0%qp_dp, t0%qt_dp, t0%qr_dp, t0%ql_dp]
      analytic(:,2) = [0.0d0,    t0%qt_dt, t0%qr_dt, t0%ql_dt]
      analytic(:,3) = [t0%qp_dr, t0%qt_dr, t0%qr_dr, 0.0d0]
      analytic(:,4) = [0.0d0,    t0%qt_dl*ln10, 0.0d0, 0.0d0]
      do s = 1, 3
         h = steps(s)
         do v = 1, 4
            sp = sh0; sm = sh0
            select case (v)
            case (1)
               sp%log_pressure = sh0%log_pressure + h/ln10
               sm%log_pressure = sh0%log_pressure - h/ln10
               sp%log_pressure_delta = sh0%log_pressure_delta + h/ln10
               sm%log_pressure_delta = sh0%log_pressure_delta - h/ln10
            case (2)
               sp%log_temperature = sh0%log_temperature + h/ln10
               sm%log_temperature = sh0%log_temperature - h/ln10
               sp%log_temperature_delta = sh0%log_temperature_delta + h/ln10
               sm%log_temperature_delta = sh0%log_temperature_delta - h/ln10
            case (3)
               sp%log_radius = sh0%log_radius + h/ln10
               sm%log_radius = sh0%log_radius - h/ln10
            case (4)
               sp%luminosity_lsun = sh0%luminosity_lsun*(1.0d0 + h)
               sm%luminosity_lsun = sh0%luminosity_lsun*(1.0d0 - h)
            end select
            call evaluate(sp, eos0, qp, tdum, ddum)
            call evaluate(sm, eos0, qm, tdum, ddum)
            fd(:,v,s) = (qp - qm)/(2.0d0*h)
            if (v == 4) fd(:,v,s) = fd(:,v,s)/sh0%luminosity_lsun
         end do
      end do
      write(*,'(a,a,a,i5,a,l2,a,l2,a,f8.4)') '# shell ', label, ' (', im, ') convective', &
           d0%convective, ' burning', d0%burned, ' log T', sh0%log_temperature
      write(*,'(a,4(1pe12.4))') '  q = ', q0
      write(*,'(a)') '  term/var     analytic        f.d.   rel.diff  status'
      do k = 1, 4
         do v = 1, 4
            best = minval([(rel(analytic(k,v), fd(k,v,s), q0(k)), s=1,3)])
            if (best <= tolerance) then
               status = 'ok'
            else if (d0%convective .and. k == 2 .and. v /= 3) then
               status = 'known'
            else
               status = 'FAIL'
               nfail = nfail + 1
            end if
            write(*,'(2x,a,a,a,1pe12.4,1pe12.4,1pe9.1,2x,a)') qname(k), '/', xname(v), &
                 analytic(k,v), fd(k,v,2), best, trim(status)
         end do
      end do
end subroutine check_shell

! The relative difference of an analytic partial a and a finite
! difference f of a term with value q: |f - a|/max(|a|, |f|, 1e-6 |q|).
double precision function rel(a, f, q)
      double precision, intent(in) :: a, f, q
      rel = abs(f - a)/max(abs(a), abs(f), 1.0d-6*abs(q), tiny(1.0d0))
end function rel

end program test_structure
