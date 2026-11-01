!----------------------------------------------------------------------
! stitched_model
!----------------------------------------------------------------------
! The one place that assembles the whole converged star -- interior
! (center -> fitting point) + envelope (fitting point -> photosphere) +
! atmosphere (photosphere -> tau ~ 0) -- and computes every quantity
! the profile and pulse writers need, as plain arrays over the extended
! grid. build_stitched_model re-integrates the envelope and atmosphere
! at the converged (Teff, L) with the fixed output step sizes, then
! fills:
!   s%stitched%prof(1:n_prof_cols, 1:n_ext) -- the profile columns, at
!     their named indices ip_<name>; column ip_zone is left 0 (the
!     writer numbers the zones)
!   s%stitched%pulse(1:n_pulse_cols, 1:n_ext) -- the GYRE/FGONG pulse set
!     (its first nz points, the interior's, are also the interior GYRE
!     file every pulse_gyre_interval models)
! The profile column table (profile_columns: each column's name and
! kind) is here, next to the values (interior_value, envelope_value,
! atmosphere_value). The io writers only read these arrays; nothing in
! the output path integrates or stitches anything. compute_observables
! reads them too (the turnover time, the seismic observables).
module stitched_model
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use physical_constants
      use math_lib
      use star_info_lib, only: json, star_info
      use output_columns, only: output_column, real_column, integer_column
      implicit none
      private
      public :: build_stitched_model, profile_columns, n_prof_cols, n_pulse_cols
      public :: ip_zone, ip_mass, ip_logR, ip_logT, ip_logRho, ip_logP, &
           ip_mixing_type, ip_gradr, ip_gradT, ip_grada, ip_conv_vel, &
           ip_brunt_N2, ip_csound
      public :: ipul_r, ipul_m, ipul_L, ipul_P, ipul_T, ipul_rho, &
           ipul_grad, ipul_N2, ipul_gamma1, ipul_grada, ipul_delta, &
           ipul_kap, ipul_kap_kap_T, ipul_kap_kap_rho, ipul_eps, &
           ipul_eps_eps_T, ipul_eps_eps_rho, ipul_omega, ipul_cp, &
           ipul_mu_e_inv, ipul_h1, ipul_z, ipul_species_base, &
           ipul_eps_grav, ipul_gyre_last

! The profile columns' named indices, in the order of the column table.
      integer, parameter :: &
           ip_zone = 1, ip_mass = 2, ip_logR = 3, ip_logT = 4, &
           ip_logRho = 5, ip_logP = 6, ip_luminosity = 7, ip_dm = 8, &
           ip_mixing_type = 9, ip_gamma1 = 10, ip_opacity = 11, ip_gradr = 12, &
           ip_gradT = 13, ip_grada = 14, ip_conv_vel = 15, ip_beta = 16, &
           ip_eta = 17, ip_grav = 18, ip_eps_ppI = 19, ip_eps_ppII = 20, &
           ip_eps_ppIII = 21, ip_eps_cno = 22, ip_eps_he3 = 23, ip_eps_nuc = 24, &
           ip_eps_neu = 25, ip_eps_grav = 26, ip_h1 = 27, ip_h2 = 28, &
           ip_he3 = 29, ip_he4 = 30, ip_li6 = 31, ip_li7 = 32, &
           ip_be9 = 33, ip_c12 = 34, ip_c13 = 35, ip_n14 = 36, &
           ip_n15 = 37, ip_o16 = 38, ip_o17 = 39, ip_o18 = 40, &
           ip_z = 41, ip_omega = 42, ip_j_rot = 43, ip_i_rot = 44, &
           ip_fp_rot = 45, ip_ft_rot = 46, ip_v_es = 47, ip_v_gsf = 48, &
           ip_v_ss = 49, ip_cp = 50, ip_delta = 51, ip_mu = 52, &
           ip_mu_e_inv = 53, ip_brunt_N2 = 54, ip_lamb_S2 = 55, ip_gradL = 56, &
           ip_gradr_div_grada = 57, ip_csound = 58, ip_D_omega = 59, ip_D_mix_rot = 60, &
           ip_eps_extra = 61, ip_v_sh = 62
      integer, parameter :: n_prof_cols = ip_v_sh
! The four seismic columns, which compute_seismic_columns fills over the
! whole extended grid (s%stitched%seismic(1:4,:), in this order).
      integer, parameter :: ip_seismic_first = ip_brunt_N2, &
           ip_seismic_last = ip_gradr_div_grada
! The regions of the extended grid (s%stitched%region).
      integer, parameter :: region_interior = 1, region_envelope = 2, &
           region_atmosphere = 3
      integer, parameter :: n_pulse_cols = 35
! The pulse columns' named indices. Columns 1..18 are the GYRE
! schema-101 set, in file order (ipul_gyre_last marks the end of that
! block); 19..34 the FGONG extras. 23..33 are the n_pulse_species
! FGONG-ordered species at ipul_species_base + k, k = 1..n_pulse_species
! (see species_slot).
      integer, parameter :: ipul_r = 1, ipul_m = 2, ipul_L = 3, &
           ipul_P = 4, ipul_T = 5, ipul_rho = 6, ipul_grad = 7, &
           ipul_N2 = 8, ipul_gamma1 = 9, ipul_grada = 10, &
           ipul_delta = 11, ipul_kap = 12, ipul_kap_kap_T = 13, &
           ipul_kap_kap_rho = 14, ipul_eps = 15, ipul_eps_eps_T = 16, &
           ipul_eps_eps_rho = 17, ipul_omega = 18, ipul_cp = 19, &
           ipul_mu_e_inv = 20, ipul_h1 = 21, ipul_z = 22, &
           ipul_species_base = 22, ipul_eps_grav = 34
      integer, parameter :: ipul_gyre_last = ipul_omega, n_pulse_species = 11
! The extended model: interior (center -> fitting point) + envelope
! (fitting point -> photosphere) + atmosphere (photosphere -> tau~0),
! assembled inward-to-outward. Profiles and pulse files both cover the
! full star: truncating at the fitting point would drop the
! superadiabatic layer and photosphere, which dominate p-mode
! frequencies.
      integer, parameter :: max_ext = 3*json
! The extended model itself, its index map, its seismic columns and the
! atmosphere heights live in the star's record s%stitched
! (state/stitched_model_state.f90), sized max_ext.

contains

! Regenerate the envelope/atmosphere structures for the CONVERGED
! model and build the inward-to-outward index map. The envelope the
! solver last integrated belongs to some trial (Teff, L) -- or was
! never integrated at all, when the envelope triangle interpolated --
! so a fresh atm_get at the converged values is required, with the
! fixed output step sizes. ierr = 1 when a profile column has no value
! in a region (fill_stitched_arrays).
subroutine build_stitched_model(s, ierr)
      use star_info_lib, only: i_h1, i_metals
      use envelope_integration, only: atm_get, envint_step_config, fixed_envint_step, &
           surface_log10_radius_cm
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr
      integer :: j, i, jerr
      type(envint_step_config) :: atm_steps, env_steps
      real(dp) :: b, gl, rl, ateffl, plim
      integer :: saha_iterations
      logical :: lprt, lsbc0

! The record's arrays, allocated and zeroed on the first build.
      if (.not. allocated(s%stitched%prof)) then
         allocate(s%stitched%region(max_ext), s%stitched%index(max_ext), &
              s%stitched%seismic(4,max_ext), s%stitched%atm_height(max_ext), &
              s%stitched%prof(n_prof_cols,max_ext), s%stitched%pulse(n_pulse_cols,max_ext))
         s%stitched%region = 0
         s%stitched%index = 0
         s%stitched%seismic = 0.0d0
         s%stitched%atm_height = 0.0d0
         s%stitched%prof = 0.0d0
         s%stitched%pulse = 0.0d0
      end if
      ierr = 0
! interior always present
      s%stitched%n_ext = 0
      do j = 1, s%nz
         s%stitched%n_ext = s%stitched%n_ext + 1
         s%stitched%region(s%stitched%n_ext) = region_interior
         s%stitched%index(s%stitched%n_ext) = j
      end do
      if (.not. s%physics%calc_envelope_flag) then
         s%stitched%n_ie = s%stitched%n_ext
! interior-only: materialize what we have. (Historical quirk kept
! byte-for-byte: the seismic columns are NOT recomputed on this
! branch -- they hold whatever the last full build left.)
         call fill_stitched_arrays(s, ierr)
         return
      end if

! ---- re-integrate at the converged model ----
! with the fixed output step sizes, for this atm_get call only
      atm_steps = fixed_envint_step(s%physics%atm_step_size)
      env_steps = fixed_envint_step(s%physics%envelope_step_size)

      saha_iterations = 0
      lprt = .false.
      lsbc0 = .false.
      b = exp(ln10*s%log_L_lsun)
      rl = surface_log10_radius_cm(s, s%log_L_lsun, s%log_Teff)
      gl = cgl + s%log_total_mass - rl - rl
      plim = s%logP(s%nz)
      if (s%convective_flag(s%nz) .and. s%physics%spot_filling_factor /= 0.0d0 &
          .and. s%physics%spot_temp_contrast /= 1.0d0) then
         ateffl = s%log_Teff - 0.25d0*log10(s%physics%spot_filling_factor* &
              pow(s%physics%spot_temp_contrast, 4.0d0) + 1.0d0 - s%physics%spot_filling_factor)
      else
         ateffl = s%log_Teff
      end if
      jerr = 0
! The stitch is the one writer of star%photosphere_logP: atm_get sets
! it from this integration at the converged model, and no other path
! recomputes it. The stitch runs every step, so photosphere_logP is
! always current.
      call atm_get(s, b, s%fp_rot(s%nz), s%ft_rot(s%nz), gl, &
           s%log_total_mass, lprt, lsbc0, plim, rl, ateffl, &
           s%xa(i_h1,s%nz), s%xa(i_metals,s%nz), s%xa(:,s%nz), saha_iterations, jerr, &
           atm_steps=atm_steps, env_steps=env_steps)

      if (jerr /= 0) then
         s%stitched%n_ie = s%stitched%n_ext
         call fill_stitched_arrays(s, ierr)
         return
      end if

! envelope: env_points runs fitting point -> photosphere (envint
! inverts it), so it appends directly. Its innermost point repeats
! the fitting point (to re-integration roundoff, so its radius can
! even land marginally BELOW the last interior point's) -- skip any
! envelope point not strictly outside the interior, keeping the
! extended radius strictly monotonic for GYRE.
      do i = 1, s%env_points%num_env_points
         if (s%stitched%n_ext >= max_ext) exit
         if (s%env_points%env_log10_radius(i) <= s%logR(s%nz)) cycle
         s%stitched%n_ext = s%stitched%n_ext + 1
         s%stitched%region(s%stitched%n_ext) = region_envelope
         s%stitched%index(s%stitched%n_ext) = i
      end do
      s%stitched%n_ie = s%stitched%n_ext
! atmosphere: atm_points runs outward-in, so walk it in reverse.
! atmo_delta_depth(k) is the PER-STEP geometric length between tau
! points k-1 and k (envint's Cox p590 quadrature), not a height, so
! accumulate it into atm_height: point i sits sum(delta(i+1:num))
! above the deepest point (tau = 2/3, the photosphere) (Sources, 1).
! delta(1), the one value polluted by the stale-prev_tau quirk envint
! notes, never enters any height. The
! deepest point (height 0) repeats the envelope's photosphere radius
! and is skipped, same as the fitting-point repeat above.
      if (s%atm_points%num_atm_points > 0) then
         s%stitched%atm_height(s%atm_points%num_atm_points) = 0.0d0
         do i = s%atm_points%num_atm_points - 1, 1, -1
            s%stitched%atm_height(i) = s%stitched%atm_height(i+1) + &
                 s%atm_points%atmo_delta_depth(i+1)
         end do
      end if
      do i = s%atm_points%num_atm_points, 1, -1
         if (s%stitched%n_ext >= max_ext) exit
         if (s%stitched%atm_height(i) <= 0.0d0) cycle
         s%stitched%n_ext = s%stitched%n_ext + 1
         s%stitched%region(s%stitched%n_ext) = region_atmosphere
         s%stitched%index(s%stitched%n_ext) = i
      end do
      call compute_seismic_columns(s)
      call fill_stitched_arrays(s, ierr)
end subroutine build_stitched_model

! ---------------------------------------------------------------
! Fill s%stitched%seismic(1:4,:) over the assembled extended grid:
!   1 brunt_N2 = g*[(1/Gamma1) dlnP/dr - dlnRho/dr]  (centered
!     differences in r -- the same derivative content as FGONG's A4)
!   2 lamb_S2  = l(l+1) c^2 / r^2 with l=1: 2*Gamma1*P/(rho r^2)
!   3 gradL    = grada + (1/delta)*dln(mu)/dlnP  (Ledoux gradient
!     with the ideal-fully-ionized phi=1 approximation; the mu
!     gradient is a centered difference in lnP)
!   4 gradr_div_grada
! Endpoints copy their neighbor's derivative-based values.
subroutine compute_seismic_columns(s)
      use math_lib
      type(star_info), intent(inout) :: s
      integer :: j
      real(dp) :: r(max_ext), lnp(max_ext), lnrho(max_ext), &
           lnmu(max_ext), g1(max_ext), grada_(max_ext), delta_(max_ext), &
           gradr_(max_ext), mass_g(max_ext)
      real(dp) :: dr, dlnp_dr, dlnrho_dr, grav, dlnmu_dlnp, p, rho
      logical :: mu_ok(max_ext)

      mu_ok = .true.
      do j = 1, s%stitched%n_ext
         r(j)      = exp(ln10*ext_profile_value(s, ip_logR, j))
         lnp(j)    = ln10*ext_profile_value(s, ip_logP, j)
         lnrho(j)  = ln10*ext_profile_value(s, ip_logRho, j)
         g1(j)     = ext_profile_value(s, ip_gamma1, j)
         grada_(j) = ext_profile_value(s, ip_grada, j)
         gradr_(j) = ext_profile_value(s, ip_gradr, j)
         delta_(j) = ext_profile_value(s, ip_delta, j)
         mass_g(j) = ext_profile_value(s, ip_mass, j)*s%solar_mass_cgs
         lnmu(j)   = ext_profile_value(s, ip_mu, j)
         if (lnmu(j) > 0.0d0) then
            lnmu(j) = log(lnmu(j))
         else
            lnmu(j) = 0.0d0
            mu_ok(j) = .false.
         end if
      end do
      do j = 1, s%stitched%n_ext
         p   = exp(lnp(j))
         rho = exp(lnrho(j))
! lamb_S2 (pointwise)
         if (r(j) > 0.0d0 .and. rho > 0.0d0 .and. g1(j) > 0.0d0) then
            s%stitched%seismic(2,j) = 2.0d0*g1(j)*p/(rho*r(j)*r(j))
         else
            s%stitched%seismic(2,j) = 0.0d0
         end if
! gradr_div_grada (pointwise)
         if (grada_(j) /= 0.0d0) then
            s%stitched%seismic(4,j) = gradr_(j)/grada_(j)
         else
            s%stitched%seismic(4,j) = 0.0d0
         end if
      end do
      do j = 2, s%stitched%n_ext-1
         dr = r(j+1) - r(j-1)
         if (dr /= 0.0d0 .and. r(j) > 0.0d0 .and. g1(j) > 0.0d0) then
            dlnp_dr   = (lnp(j+1) - lnp(j-1))/dr
            dlnrho_dr = (lnrho(j+1) - lnrho(j-1))/dr
            grav = local_gravity_cgs(mass_g(j), r(j))
            s%stitched%seismic(1,j) = grav*(dlnp_dr/g1(j) - dlnrho_dr)
         else
            s%stitched%seismic(1,j) = 0.0d0
         end if
! the mu gradient's term, dropped where mu is unavailable (Sources, 2)
         if (lnp(j+1) /= lnp(j-1) .and. delta_(j) /= 0.0d0 .and. &
              mu_ok(j-1) .and. mu_ok(j+1)) then
            dlnmu_dlnp = (lnmu(j+1) - lnmu(j-1))/(lnp(j+1) - lnp(j-1))
            s%stitched%seismic(3,j) = grada_(j) + dlnmu_dlnp/delta_(j)
         else
            s%stitched%seismic(3,j) = grada_(j)
         end if
      end do
      if (s%stitched%n_ext >= 2) then
         s%stitched%seismic(1,1) = s%stitched%seismic(1,2)
         s%stitched%seismic(3,1) = s%stitched%seismic(3,2)
         s%stitched%seismic(1,s%stitched%n_ext) = s%stitched%seismic(1,s%stitched%n_ext-1)
         s%stitched%seismic(3,s%stitched%n_ext) = s%stitched%seismic(3,s%stitched%n_ext-1)
      end if
end subroutine compute_seismic_columns

! ---------------------------------------------------------------
! The profile column table: each column's name and kind, at its named
! index (profile_output writes the selected ones).
function profile_columns() result(cols)
      type(output_column) :: cols(n_prof_cols)

      cols(ip_zone)            = output_column('zone', integer_column)
      cols(ip_mass)            = output_column('mass', real_column)
      cols(ip_logR)            = output_column('logR', real_column)
      cols(ip_logT)            = output_column('logT', real_column)
      cols(ip_logRho)          = output_column('logRho', real_column)
      cols(ip_logP)            = output_column('logP', real_column)
      cols(ip_luminosity)      = output_column('luminosity', real_column)
      cols(ip_dm)              = output_column('dm', real_column)
      cols(ip_mixing_type)     = output_column('mixing_type', integer_column)
      cols(ip_gamma1)          = output_column('gamma1', real_column)
      cols(ip_opacity)         = output_column('opacity', real_column)
      cols(ip_gradr)           = output_column('gradr', real_column)
      cols(ip_gradT)           = output_column('gradT', real_column)
      cols(ip_grada)           = output_column('grada', real_column)
      cols(ip_conv_vel)        = output_column('conv_vel', real_column)
      cols(ip_beta)            = output_column('beta', real_column)
      cols(ip_eta)             = output_column('eta', real_column)
      cols(ip_grav)            = output_column('grav', real_column)
      cols(ip_eps_ppI)         = output_column('eps_ppI', real_column)
      cols(ip_eps_ppII)        = output_column('eps_ppII', real_column)
      cols(ip_eps_ppIII)       = output_column('eps_ppIII', real_column)
      cols(ip_eps_cno)         = output_column('eps_cno', real_column)
      cols(ip_eps_he3)         = output_column('eps_he3', real_column)
      cols(ip_eps_nuc)         = output_column('eps_nuc', real_column)
      cols(ip_eps_neu)         = output_column('eps_neu', real_column)
      cols(ip_eps_grav)        = output_column('eps_grav', real_column)
      cols(ip_h1)              = output_column('h1', real_column)
      cols(ip_h2)              = output_column('h2', real_column)
      cols(ip_he3)             = output_column('he3', real_column)
      cols(ip_he4)             = output_column('he4', real_column)
      cols(ip_li6)             = output_column('li6', real_column)
      cols(ip_li7)             = output_column('li7', real_column)
      cols(ip_be9)             = output_column('be9', real_column)
      cols(ip_c12)             = output_column('c12', real_column)
      cols(ip_c13)             = output_column('c13', real_column)
      cols(ip_n14)             = output_column('n14', real_column)
      cols(ip_n15)             = output_column('n15', real_column)
      cols(ip_o16)             = output_column('o16', real_column)
      cols(ip_o17)             = output_column('o17', real_column)
      cols(ip_o18)             = output_column('o18', real_column)
      cols(ip_z)               = output_column('z', real_column)
      cols(ip_omega)           = output_column('omega', real_column)
      cols(ip_j_rot)           = output_column('j_rot', real_column)
      cols(ip_i_rot)           = output_column('i_rot', real_column)
      cols(ip_fp_rot)          = output_column('fp_rot', real_column)
      cols(ip_ft_rot)          = output_column('ft_rot', real_column)
      cols(ip_v_es)            = output_column('v_es', real_column)
      cols(ip_v_gsf)           = output_column('v_gsf', real_column)
      cols(ip_v_ss)            = output_column('v_ss', real_column)
      cols(ip_cp)              = output_column('cp', real_column)
      cols(ip_delta)           = output_column('delta', real_column)
      cols(ip_mu)              = output_column('mu', real_column)
      cols(ip_mu_e_inv)        = output_column('mu_e_inv', real_column)
      cols(ip_brunt_N2)        = output_column('brunt_N2', real_column)
      cols(ip_lamb_S2)         = output_column('lamb_S2', real_column)
      cols(ip_gradL)           = output_column('gradL', real_column)
      cols(ip_gradr_div_grada) = output_column('gradr_div_grada', real_column)
      cols(ip_csound)          = output_column('csound', real_column)
      cols(ip_D_omega)         = output_column('D_omega', real_column)
      cols(ip_D_mix_rot)       = output_column('D_mix_rot', real_column)
      cols(ip_eps_extra)       = output_column('eps_extra', real_column)
      cols(ip_v_sh)            = output_column('v_sh', real_column)
end function profile_columns

! ---------------------------------------------------------------
! Profile column icol at extended point j: a seismic column from
! compute_seismic_columns' pass, otherwise its value in the point's
! region. found is false when that region has no case for the column
! (a column added without one), which fill_stitched_arrays reports.
subroutine extended_point_value(s, icol, j, value, found)
      type(star_info), intent(in) :: s
      integer, intent(in) :: icol, j
      real(dp), intent(out) :: value
      logical, intent(out) :: found
      integer :: i

      if (icol >= ip_seismic_first .and. icol <= ip_seismic_last) then
         value = s%stitched%seismic(icol-ip_seismic_first+1, j)
         found = .true.
         return
      end if
      i = s%stitched%index(j)
      select case (s%stitched%region(j))
      case (region_interior)
         call interior_value(s, icol, i, value, found)
      case (region_envelope)
         call envelope_value(s, icol, i, value, found)
      case (region_atmosphere)
         call atmosphere_value(s, icol, i, value, found)
      case default
         value = 0.0d0
         found = .false.
      end select
end subroutine extended_point_value

! Profile column icol at extended point j, for compute_seismic_columns,
! which reads only columns every region has (fill_stitched_arrays
! reports any that is missing).
real(dp) function ext_profile_value(s, icol, j)
      type(star_info), intent(in) :: s
      integer, intent(in) :: icol, j
      logical :: found
      call extended_point_value(s, icol, j, ext_profile_value, found)
end function ext_profile_value

! ---------------------------------------------------------------
! Column kernels: the expressions the region values, compute_seismic_columns
! and build_pulse_points evaluate at more than one site, each written
! once; the sites differ only in which array they read (the interior's
! star%, the envelope's env_points, the atmosphere's atm_points). Two
! are not shared, because their sites differ in operand order (and so
! in the last bit): the interior gravity of column ip_grav,
! exp(ln10*(cgl - 2*logR))*m, versus local_gravity_cgs; and the N^2 of
! compute_seismic_columns, grav*(dlnp_dr/g1 - dlnrho_dr) with the two
! derivatives divided by dr first, versus build_pulse_points' second
! pass, which divides the bracket by dr last.

! Local gravity G*m/r^2 [cgs] from mass [g] and radius [cm].
real(dp) function local_gravity_cgs(mass_g, radius_cm)
      real(dp), intent(in) :: mass_g, radius_cm
      local_gravity_cgs = exp(ln10*cgl)*mass_g/(radius_cm*radius_cm)
end function local_gravity_cgs

! Sound speed sqrt(Gamma1*P/rho) [cm/s] from Gamma1, log10 P, log10 rho
! (profile column 58).
real(dp) function sound_speed_cgs(gamma1, log10_p, log10_rho)
      real(dp), intent(in) :: gamma1, log10_p, log10_rho
      sound_speed_cgs = sqrt(gamma1*exp(ln10*(log10_p - log10_rho)))
end function sound_speed_cgs

! Profile column 9: 1 in a convective shell, 0 otherwise.
real(dp) function convective_indicator(is_convective)
      logical, intent(in) :: is_convective
      if (is_convective) then
         convective_indicator = 1.0d0
      else
         convective_indicator = 0.0d0
      end if
end function convective_indicator

! Mass [g] enclosed at envelope point i (env_log10_mass is relative to
! the total mass).
real(dp) function envelope_mass_g(s, i)
      type(star_info), intent(in) :: s
      integer, intent(in) :: i
      envelope_mass_g = exp(ln10*(s%env_points%env_log10_mass(i) + s%log_total_mass))
end function envelope_mass_g

! Radius [cm] of atmosphere point i: the photosphere (the envelope's
! outermost point) plus its accumulated height (see build_stitched_model).
real(dp) function atmosphere_radius_cm(s, i)
      type(star_info), intent(in) :: s
      integer, intent(in) :: i
      atmosphere_radius_cm = exp(ln10* &
           s%env_points%env_log10_radius(s%env_points%num_env_points)) + s%stitched%atm_height(i)
end function atmosphere_radius_cm

! mu = R*rho*T/(beta*P): the ideal-gas mean molecular weight from the
! gas pressure (compute_blended_eos's own specific-gas-constant identity). Zero
! when the inputs cannot support it.
real(dp) function ideal_gas_mu(log10_p, log10_t, log10_rho, beta)
      use math_lib
      use physical_constants, only: gas_constant
      real(dp), intent(in) :: log10_p, log10_t, log10_rho, beta
      if (beta > 0.0d0) then
         ideal_gas_mu = gas_constant*exp(ln10*(log10_rho + log10_t - log10_p))/beta
      else
         ideal_gas_mu = 0.0d0
      end if
end function ideal_gas_mu

! ---------------------------------------------------------------
! Profile column icol at interior zone k (1 = center .. nz = surface):
! the converged model's values. found is false for a column without a
! case.
subroutine interior_value(s, icol, k, value, found)
      use math_lib
      use star_info_lib, only: i_be9, i_c12, i_c13, i_eps_cno, i_eps_grav, &
           i_eps_3alpha, i_eps_neu, i_eps_pp1, i_eps_pp2, i_eps_pp3, i_h1, &
           i_h2, i_he3, i_he4, i_li6, i_li7, i_metals, i_n14, i_n15, i_o16, &
           i_o17, i_o18
      use star_info_lib, only: i_mech_eddington_sweet, i_mech_gsf, i_mech_shear, &
           i_mech_solberg_hoiland
      type(star_info), intent(in) :: s
      integer, intent(in) :: icol, k
      real(dp), intent(out) :: value
      logical, intent(out) :: found

      found = .true.
      select case (icol)
      case (ip_zone);         value = 0.0d0   ! the writer numbers the zones
      case (ip_mass);         value = s%m(k)/s%solar_mass_cgs
      case (ip_logR);         value = s%logR(k)
      case (ip_logT);         value = s%logT(k)
      case (ip_logRho);       value = s%logRho(k)
      case (ip_logP);         value = s%logP(k)
      case (ip_luminosity);   value = s%luminosity_lsun(k)
      case (ip_dm);           value = s%dm(k)/s%solar_mass_cgs
      case (ip_mixing_type);  value = convective_indicator(s%convective_flag(k))
      case (ip_gamma1);       value = s%adiabatic_index_gamma1(k)
      case (ip_opacity);      value = s%opacity_zone(k)
      case (ip_gradr);        value = s%gradr(k)
      case (ip_gradT);        value = s%gradT(k)
      case (ip_grada);        value = s%grada(k)
      case (ip_conv_vel);     value = s%conv_vel(k)
      case (ip_beta);         value = s%beta(k)
      case (ip_eta);          value = s%eta(k)
      case (ip_grav);         value = exp(ln10*(cgl - 2.0d0*s%logR(k)))*s%m(k)
! eps_ppI .. eps_he3 are fractions of eps_nuc, one per chain with a
! slot (eps_he3 is the triple-alpha share); the alpha captures have none,
! so where they burn the five add up to 1 minus their share: in a 1 Msun
! core-helium-burning model that share is 4% at the centre at Yc = 0.76,
! 18% at Yc = 0.30 and 99% at Yc = 0.002, and nearly all of eps_nuc
! toward the convective core's edge, where 12C(alpha,gamma) outruns the
! triple alpha. eps_nuc, eps_neu and eps_grav are erg/g/s.
      case (ip_eps_ppI);      value = s%eps_channels(i_eps_pp1,k)
      case (ip_eps_ppII);     value = s%eps_channels(i_eps_pp2,k)
      case (ip_eps_ppIII);    value = s%eps_channels(i_eps_pp3,k)
      case (ip_eps_cno);      value = s%eps_channels(i_eps_cno,k)
      case (ip_eps_he3);      value = s%eps_channels(i_eps_3alpha,k)
      case (ip_eps_nuc);      value = s%eps_total(k)
      case (ip_eps_neu);      value = s%eps_channels(i_eps_neu,k)
      case (ip_eps_grav);     value = s%eps_channels(i_eps_grav,k)
      case (ip_h1);           value = s%xa(i_h1,k)
      case (ip_h2);           value = s%xa(i_h2,k)
      case (ip_he3);          value = s%xa(i_he3,k)
      case (ip_he4);          value = s%xa(i_he4,k)
      case (ip_li6);          value = s%xa(i_li6,k)
      case (ip_li7);          value = s%xa(i_li7,k)
      case (ip_be9);          value = s%xa(i_be9,k)
      case (ip_c12);          value = s%xa(i_c12,k)
      case (ip_c13);          value = s%xa(i_c13,k)
      case (ip_n14);          value = s%xa(i_n14,k)
      case (ip_n15);          value = s%xa(i_n15,k)
      case (ip_o16);          value = s%xa(i_o16,k)
      case (ip_o17);          value = s%xa(i_o17,k)
      case (ip_o18);          value = s%xa(i_o18,k)
      case (ip_z);            value = s%xa(i_metals,k)
      case (ip_omega);        value = s%omega(k)
      case (ip_j_rot);        value = s%j_rot(k)
      case (ip_i_rot);        value = s%i_rot(k)
      case (ip_fp_rot);       value = s%fp_rot(k)
      case (ip_ft_rot);       value = s%ft_rot(k)
      case (ip_v_es);         value = s%mechanism_velocity(i_mech_eddington_sweet,k)
      case (ip_v_gsf);        value = s%mechanism_velocity(i_mech_gsf,k)
      case (ip_v_ss);         value = s%mechanism_velocity(i_mech_shear,k)
      case (ip_cp);           value = s%pulse_specific_heat(k)
      case (ip_delta);        value = -s%pulse_dlnrho_dlnt(k)
      case (ip_mu);           value = s%pulse_mean_molecular_weight(k)
! (the array holds 1/mu_e, eos's electron_mean_weight_inverse)
      case (ip_mu_e_inv);     value = s%pulse_electron_mean_weight_inverse(k)
      case (ip_csound);       value = sound_speed_cgs(s%adiabatic_index_gamma1(k), &
                                   s%logP(k), s%logRho(k))
      case (ip_D_omega);      value = s%am_diffusion_coeff(k)
      case (ip_D_mix_rot);    value = s%mixing_diffusion_coeff(k)
      case (ip_eps_extra);    value = s%eps_extra(k)
      case (ip_v_sh);         value = s%mechanism_velocity(i_mech_solberg_hoiland,k)
      case default
         value = 0.0d0
         found = .false.
      end select
end subroutine interior_value

! ---------------------------------------------------------------
! Profile column icol at envelope point i: what the envelope
! integration stores (env_points), the surface omega, and zero for the
! columns it does not compute (burning, the composition beyond X and Z,
! the rotation's internals and transport). found is false for a column
! without a case.
subroutine envelope_value(s, icol, i, value, found)
      type(star_info), intent(in) :: s
      integer, intent(in) :: icol, i
      real(dp), intent(out) :: value
      logical, intent(out) :: found

      found = .true.
      select case (icol)
      case (ip_mass);         value = envelope_mass_g(s, i)/s%solar_mass_cgs
      case (ip_logR);         value = s%env_points%env_log10_radius(i)
      case (ip_logT);         value = s%env_points%env_log10_temperature(i)
      case (ip_logRho);       value = s%env_points%env_log10_density(i)
      case (ip_logP);         value = s%env_points%env_log10_pressure(i)
      case (ip_luminosity);   value = s%env_points%env_luminosity(i)
      case (ip_mixing_type);  value = convective_indicator(s%env_points%env_convective_flag(i))
      case (ip_gamma1);       value = s%env_points%env_gamma1(i)
      case (ip_opacity);      value = s%env_points%env_opacity(i)
! env_gradients holds (1) the radiative, (2) the adiabatic and (3) the
! actual gradient (envint's current_gradients).
      case (ip_gradr);        value = s%env_points%env_gradients(1,i)
      case (ip_gradT);        value = s%env_points%env_gradients(3,i)
      case (ip_grada);        value = s%env_points%env_gradients(2,i)
      case (ip_conv_vel);     value = s%env_points%env_convective_velocity(i)
      case (ip_beta);         value = s%env_points%env_beta(i)
      case (ip_h1);           value = s%env_points%env_hydrogen_fraction(i)
      case (ip_z);            value = s%env_points%env_metal_fraction(i)
      case (ip_omega);        value = s%omega(s%nz)
      case (ip_cp);           value = s%env_points%env_specific_heat_cp(i)
      case (ip_delta);        value = -s%env_points%env_dlnrho_dlnt(i)
! mu from the gas pressure, R/mu = beta*P/(rho*T) (compute_blended_eos's
! specific gas constant), continuous across the fitting point
      case (ip_mu);           value = ideal_gas_mu( &
                                   s%env_points%env_log10_pressure(i), &
                                   s%env_points%env_log10_temperature(i), &
                                   s%env_points%env_log10_density(i), s%env_points%env_beta(i))
      case (ip_csound);       value = sound_speed_cgs(s%env_points%env_gamma1(i), &
                                   s%env_points%env_log10_pressure(i), &
                                   s%env_points%env_log10_density(i))
      case (ip_zone, ip_dm, ip_eta, ip_grav, &
           ip_eps_ppI, ip_eps_ppII, ip_eps_ppIII, ip_eps_cno, ip_eps_he3, &
           ip_eps_nuc, ip_eps_neu, ip_eps_grav, ip_eps_extra, &
           ip_h2, ip_he3, ip_he4, ip_li6, ip_li7, ip_be9, ip_c12, ip_c13, &
           ip_n14, ip_n15, ip_o16, ip_o17, ip_o18, &
           ip_j_rot, ip_i_rot, ip_fp_rot, ip_ft_rot, ip_v_es, ip_v_gsf, &
           ip_v_ss, ip_v_sh, ip_mu_e_inv, ip_D_omega, ip_D_mix_rot)
         value = 0.0d0
      case default
         value = 0.0d0
         found = .false.
      end select
end subroutine envelope_value

! ---------------------------------------------------------------
! Profile column icol at atmosphere point i: what the atmosphere
! integration stores (atm_points), the surface's mass, luminosity,
! composition and omega, and zero for the columns it does not compute
! (the envelope's, and the mixing type and convective velocity). found
! is false for a column without a case.
subroutine atmosphere_value(s, icol, i, value, found)
      use star_info_lib, only: i_h1, i_metals
      type(star_info), intent(in) :: s
      integer, intent(in) :: icol, i
      real(dp), intent(out) :: value
      logical, intent(out) :: found

      found = .true.
      select case (icol)
      case (ip_mass);         value = s%star_mass
      case (ip_logR);         value = log10(atmosphere_radius_cm(s, i))
      case (ip_logT);         value = s%atm_points%atmo_log10_temperature(i)
      case (ip_logRho);       value = s%atm_points%atmo_log10_density(i)
      case (ip_logP);         value = s%atm_points%atmo_log10_pressure(i)
      case (ip_luminosity);   value = exp(ln10*s%log_L_lsun)
      case (ip_gamma1);       value = s%atm_points%atmo_gamma1(i)
      case (ip_opacity);      value = s%atm_points%atmo_opacity(i)
! atmo_gradients holds (1) the radiative, (2) the actual and (3) the
! adiabatic gradient.
      case (ip_gradr);        value = s%atm_points%atmo_gradients(1,i)
      case (ip_gradT);        value = s%atm_points%atmo_gradients(2,i)
      case (ip_grada);        value = s%atm_points%atmo_gradients(3,i)
      case (ip_beta);         value = s%atm_points%atmo_beta(i)
      case (ip_h1);           value = s%xa(i_h1,s%nz)
      case (ip_z);            value = s%xa(i_metals,s%nz)
      case (ip_omega);        value = s%omega(s%nz)
      case (ip_cp);           value = s%atm_points%atmo_specific_heat_cp(i)
      case (ip_delta);        value = -s%atm_points%atmo_dlnrho_dlnt(i)
      case (ip_mu);           value = ideal_gas_mu( &
                                   s%atm_points%atmo_log10_pressure(i), &
                                   s%atm_points%atmo_log10_temperature(i), &
                                   s%atm_points%atmo_log10_density(i), s%atm_points%atmo_beta(i))
      case (ip_csound);       value = sound_speed_cgs(s%atm_points%atmo_gamma1(i), &
                                   s%atm_points%atmo_log10_pressure(i), &
                                   s%atm_points%atmo_log10_density(i))
      case (ip_zone, ip_dm, ip_eta, ip_grav, ip_mixing_type, ip_conv_vel, &
           ip_eps_ppI, ip_eps_ppII, ip_eps_ppIII, ip_eps_cno, ip_eps_he3, &
           ip_eps_nuc, ip_eps_neu, ip_eps_grav, ip_eps_extra, &
           ip_h2, ip_he3, ip_he4, ip_li6, ip_li7, ip_be9, ip_c12, ip_c13, &
           ip_n14, ip_n15, ip_o16, ip_o17, ip_o18, &
           ip_j_rot, ip_i_rot, ip_fp_rot, ip_ft_rot, ip_v_es, ip_v_gsf, &
           ip_v_ss, ip_v_sh, ip_mu_e_inv, ip_D_omega, ip_D_mix_rot)
         value = 0.0d0
      case default
         value = 0.0d0
         found = .false.
      end select
end subroutine atmosphere_value

! Assemble the per-point pulse data for the FULL extended model
! (interior + envelope + atmosphere; build_stitched_model must have run).
! Columns 1-18 are the GYRE schema-101 set; 19-34 the extras FGONG
! needs. Envelope/atmosphere points carry what those integrations
! store, the opacity derivatives (kap_kap_T, kap_kap_rho) among them
! (Sources, 3). The energy generation derivative columns are zero there
! (no burning above the fitting point), and the composition above the
! fitting point is the surface composition.
subroutine build_pulse_points(s, pts)
      use math_lib
      use star_info_lib, only: i_eps_grav, i_h1, i_metals
      type(star_info), intent(inout) :: s
      real(dp), intent(out) :: pts(n_pulse_cols, s%stitched%n_ext)
      integer :: j, i, k
      real(dp) :: r, m, P, T, rho, delta, nab, nab_ad, grav

      do j = 1, s%stitched%n_ext
         i = s%stitched%index(j)
         pts(:,j) = 0.0d0
         select case (s%stitched%region(j))
         case (region_interior)
            call gather_interior_point()
         case (region_envelope)
            call gather_envelope_point()
         case default   ! region_atmosphere
            call gather_atmosphere_point()
         end select
         call store_common_columns()
      end do
      call difference_brunt_n2()

contains

! ---------------------------------------------------------------
! Interior point i: its structure, and the pulse columns the converged
! interior model holds.
subroutine gather_interior_point()
      r = exp(ln10*s%logR(i))
      m = s%m(i)
      P = exp(ln10*s%logP(i))
      T = exp(ln10*s%logT(i))
      rho = exp(ln10*s%logRho(i))
      delta = -s%pulse_dlnrho_dlnt(i)
      nab = s%gradT(i)
      nab_ad = s%grada(i)
      pts(ipul_L,j) = s%luminosity_lsun(i)*s%solar_luminosity_cgs
      pts(ipul_gamma1,j) = s%adiabatic_index_gamma1(i)
      pts(ipul_kap,j) = s%opacity_zone(i)
! GSM/GYRE convention (MESA pulse_gyre.f90): kap_kap_T = kap*dlnkap/dlnT
! (the absolute derivative dkap/dlnT), NOT the bare log-derivative;
! likewise eps_eps_T = eps*dlneps/dlnT = deps/dlnT. Only nonadiabatic
! GYRE runs read these columns. pulse_dlnkap_* are true log-derivatives,
! so they are multiplied by kap here. pulse_deps_* already hold the
! absolute derivatives of eps_total, the eps written beside them
! (nuclear_energy's deps_nuc_*: eps_nuc's, without the thermal neutrino
! loss, as MESA writes eps_nuc with d_epsnuc_dlnT and d_epsnuc_dlnd), so
! they are written as they are (Sources, 4).
      pts(ipul_kap_kap_T,j) = s%opacity_zone(i)*s%pulse_dlnkap_dlnt(i)
      pts(ipul_kap_kap_rho,j) = s%opacity_zone(i)*s%pulse_dlnkap_dlnrho(i)
      pts(ipul_eps,j) = s%eps_total(i)
      pts(ipul_eps_eps_T,j) = s%pulse_deps_dlnt(i)
      pts(ipul_eps_eps_rho,j) = s%pulse_deps_dlnrho(i)
      pts(ipul_omega,j) = s%omega(i)
      pts(ipul_cp,j) = s%pulse_specific_heat(i)
! (the array holds 1/mu_e, which FGONG's var(14) wants)
      pts(ipul_mu_e_inv,j) = s%pulse_electron_mean_weight_inverse(i)
      pts(ipul_h1,j) = s%xa(i_h1,i)
      pts(ipul_z,j) = s%xa(i_metals,i)
      do k = 1, n_pulse_species
         pts(ipul_species_base+k,j) = s%xa(species_slot(k),i)
      end do
      pts(ipul_eps_grav,j) = s%eps_channels(i_eps_grav,i)
end subroutine gather_interior_point

! ---------------------------------------------------------------
! Envelope point i: its structure and the columns env_points holds; the
! rotation and composition are the surface's.
subroutine gather_envelope_point()
      r = exp(ln10*s%env_points%env_log10_radius(i))
      m = envelope_mass_g(s, i)
      P = exp(ln10*s%env_points%env_log10_pressure(i))
      T = exp(ln10*s%env_points%env_log10_temperature(i))
      rho = exp(ln10*s%env_points%env_log10_density(i))
      delta = -s%env_points%env_dlnrho_dlnt(i)
! env_gradients holds (1) the radiative, (2) the adiabatic and (3) the
! actual gradient (envelope_derivs' current_gradients), unlike
! atmo_gradients, whose order is (radiative, actual, adiabatic)
! (Sources, 5).
      nab = s%env_points%env_gradients(3,i)
      nab_ad = s%env_points%env_gradients(2,i)
      pts(ipul_L,j) = s%env_points%env_luminosity(i)*s%solar_luminosity_cgs
      pts(ipul_gamma1,j) = s%env_points%env_gamma1(i)
      pts(ipul_kap,j) = s%env_points%env_opacity(i)
! absolute derivatives kap*dlnkap/dln{T,rho}, the GSM/GYRE convention
! used for the interior above
      pts(ipul_kap_kap_T,j) = s%env_points%env_opacity(i)*s%env_points%env_dlnkap_dlnt(i)
      pts(ipul_kap_kap_rho,j) = s%env_points%env_opacity(i)*s%env_points%env_dlnkap_dlnrho(i)
      pts(ipul_omega,j) = s%omega(s%nz)
      pts(ipul_cp,j) = s%env_points%env_specific_heat_cp(i)
      pts(ipul_h1,j) = s%env_points%env_hydrogen_fraction(i)
      pts(ipul_z,j) = s%env_points%env_metal_fraction(i)
      do k = 1, n_pulse_species
         pts(ipul_species_base+k,j) = s%xa(species_slot(k),s%nz)
      end do
end subroutine gather_envelope_point

! ---------------------------------------------------------------
! Atmosphere point i: its structure and the columns atm_points holds;
! the luminosity, rotation and composition are the surface's.
subroutine gather_atmosphere_point()
      r = atmosphere_radius_cm(s, i)
      m = exp(ln10*s%log_total_mass)
      P = exp(ln10*s%atm_points%atmo_log10_pressure(i))
      T = exp(ln10*s%atm_points%atmo_log10_temperature(i))
      rho = exp(ln10*s%atm_points%atmo_log10_density(i))
      delta = -s%atm_points%atmo_dlnrho_dlnt(i)
      nab = s%atm_points%atmo_gradients(2,i)
      nab_ad = s%atm_points%atmo_gradients(3,i)
      pts(ipul_L,j) = exp(ln10*s%log_L_lsun)*s%solar_luminosity_cgs
      pts(ipul_gamma1,j) = s%atm_points%atmo_gamma1(i)
      pts(ipul_kap,j) = s%atm_points%atmo_opacity(i)
      pts(ipul_kap_kap_T,j) = s%atm_points%atmo_opacity(i)*s%atm_points%atmo_dlnkap_dlnt(i)
      pts(ipul_kap_kap_rho,j) = s%atm_points%atmo_opacity(i)*s%atm_points%atmo_dlnkap_dlnrho(i)
      pts(ipul_omega,j) = s%omega(s%nz)
      pts(ipul_cp,j) = s%atm_points%atmo_specific_heat_cp(i)
      pts(ipul_h1,j) = s%xa(i_h1,s%nz)
      pts(ipul_z,j) = s%xa(i_metals,s%nz)
      do k = 1, n_pulse_species
         pts(ipul_species_base+k,j) = s%xa(species_slot(k),s%nz)
      end do
end subroutine gather_atmosphere_point

! ---------------------------------------------------------------
! The columns every region fills alike: r, m, P, T, rho, the gradients,
! delta, and the pointwise thermal N^2.
subroutine store_common_columns()
      pts(ipul_r,j) = r
      pts(ipul_m,j) = m
      pts(ipul_P,j) = P
      pts(ipul_T,j) = T
      pts(ipul_rho,j) = rho
      pts(ipul_grad,j) = nab
      pts(ipul_grada,j) = nab_ad
      pts(ipul_delta,j) = delta
      if (r > 0.0d0) then
         grav = local_gravity_cgs(m, r)
         pts(ipul_N2,j) = grav*grav*(rho/P)*delta*(nab_ad - nab)
      end if
end subroutine store_common_columns

! ---------------------------------------------------------------
! The second pass over N^2: the gradient form by centered differences
! outside convective zones (see below).
subroutine difference_brunt_n2()
! Brunt-Vaisala N^2 (column 8), second pass: overwrite the pointwise
! thermal-only value computed above with the exact gradient form
!     N^2 = g * [ (1/Gamma_1) dlnP/dr - dlnRho/dr ]
! (the same derivative content as FGONG's A4 and profile column 54,
! centered differences in r). The thermal form g^2 rho delta
! (nab_ad - nab)/P is exact ONLY for homogeneous composition: it
! omits the Ledoux mu-gradient term, which dominates N^2 in the
! composition-gradient layers above a retreating core / around the
! H-burning shell -- on a 1.2 Msun subgiant it understated the
! g-cavity buoyancy by ~30% (12% in the period spacing DeltaPi_1)
! and shifted l=1 mixed-mode frequencies by up to ~16 uHz. The
! actual density gradient carries the composition term for free.
! CONVECTIVE zones keep the first-pass thermal value: there the
! mixture is homogeneous by mixing (mu-gradient zero), so the thermal
! form is already exact -- and correctly, smoothly negative -- while
! the centered difference of a near-adiabatic stratification is
! cancellation noise of random sign. The switch is the sign of the
! thermal value itself (< 0 means Schwarzschild-unstable). The
! positive N^2 spike the difference produces AT a convective-zone
! base (the mu-discontinuity's buoyancy interface) lands on the
! radiative side and is kept -- it is physics, not noise.
! Endpoints copy their neighbor, matching compute_seismic_columns.
      do j = 2, s%stitched%n_ext - 1
         if (pts(ipul_N2,j) >= 0.0d0 .and. &
              pts(ipul_r,j) > 0.0d0 .and. pts(ipul_r,j+1) > pts(ipul_r,j-1) .and. &
              pts(ipul_gamma1,j) > 0.0d0) then
            grav = local_gravity_cgs(pts(ipul_m,j), pts(ipul_r,j))
            pts(ipul_N2,j) = grav*( &
                 (log(pts(ipul_P,j+1)) - log(pts(ipul_P,j-1)))/pts(ipul_gamma1,j) - &
                 (log(pts(ipul_rho,j+1)) - log(pts(ipul_rho,j-1))) ) / &
                 (pts(ipul_r,j+1) - pts(ipul_r,j-1))
         end if
      end do
      if (s%stitched%n_ext >= 2) then
         pts(ipul_N2,1) = pts(ipul_N2,2)
         pts(ipul_N2,s%stitched%n_ext) = pts(ipul_N2,s%stitched%n_ext-1)
      end if
end subroutine difference_brunt_n2

end subroutine build_pulse_points

! star%xa slot for pulse column ipul_species_base + k (k = 1 ..
! n_pulse_species), in FGONG species order (column 34 is eps_grav, not a
! species; be9 has no FGONG column).
integer function species_slot(k)
      use star_info_lib, only: i_he3, i_c12, i_c13, i_n14, i_o16, i_h2, i_he4, i_li7, i_n15, i_o17, i_o18
      integer, intent(in) :: k
      integer, parameter :: slots(n_pulse_species) = [i_he3, i_c12, i_c13, i_n14, &
           i_o16, i_h2, i_he4, i_li7, i_n15, i_o17, i_o18]
      species_slot = slots(k)
end function species_slot

! ---------------------------------------------------------------
! The writers' arrays from the assembled grid: every profile column at
! every point, then the pulse points. A column without a value in a
! point's region is an error (ierr = 1), reported once.
subroutine fill_stitched_arrays(s, ierr)
      use run_messages, only: report_error
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr
      type(output_column) :: cols(n_prof_cols)
      integer :: j, icol
      real(dp) :: value
      logical :: found

      ierr = 0
      do j = 1, s%stitched%n_ext
         do icol = 1, n_prof_cols
            call extended_point_value(s, icol, j, value, found)
            if (.not. found) then
               cols = profile_columns()
               call report_error('the profile column ' // trim(cols(icol)%name) // &
                    ' has no value in the ' // &
                    trim(region_name(s%stitched%region(j))) // &
                    ' (core/diagnostics/stitched_model.f90)')
               ierr = 1
               return
            end if
            s%stitched%prof(icol, j) = value
         end do
      end do
      call build_pulse_points(s, s%stitched%pulse)
end subroutine fill_stitched_arrays

! The name of a region of the extended grid, for messages.
function region_name(region) result(name)
      integer, intent(in) :: region
      character(len=16) :: name
      select case (region)
      case (region_interior);   name = 'interior'
      case (region_envelope);   name = 'envelope'
      case (region_atmosphere); name = 'atmosphere'
      case default;             name = 'unknown region'
      end select
end function region_name

! ---------------------------------------------------------------
! Sources (the fixes behind some of the lines above)
! 1. The atmosphere heights: the legacy .store stitch added
!    atmo_delta_depth's raw per-step value to the photosphere radius,
!    which left the atmosphere radii non-monotonic (every atmosphere
!    point hugged r_phot); fixed in 2026, not preserved.
! 2. Column mu (ip_mu) is the real mu everywhere since 2026 (it was R/mu
!    in the interior, a wrong-signed composition term, and 0 in the
!    envelope, a log(1e-30) spike at the junction).
! 3. The envelope and atmosphere integrations always computed
!    dlnkap/dlnT and dlnkap/dlnrho but did not keep them before 2026,
!    so kap_kap_T and kap_kap_rho were zero above the fitting point,
!    which removed the kappa mechanism of the hydrogen ionisation zone
!    and the atmosphere's radiative damping for a nonadiabatic
!    pulsation code (found with GYRE and nap on the standard solar
!    model).
! 4. 2ed398ca multiplied pulse_deps_* by eps as well (misled by their
!    former pulse_dlneps_* names): wrong by a factor eps until the
!    2026-10 audit. Until 2026-10-29 they were the Henyey's derivatives,
!    which include the thermal neutrino loss that eps_total leaves out,
!    so above log T = 7.5 (where the loss is computed) eps and its
!    derivatives described two different eps.
! 5. The envelope's two gradient slots were swapped here before 2026
!    (as in the profile columns), which flipped the thermal N^2 across
!    the envelope of every pulse file.
end module stitched_model
