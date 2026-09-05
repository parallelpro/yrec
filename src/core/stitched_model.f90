!----------------------------------------------------------------------
! stitched_model_lib
!----------------------------------------------------------------------
! Added 2026 (stitched-model restructure): the one place that
! assembles the FULL converged star -- interior (center -> fitting
! point) + envelope (fitting point -> photosphere) + atmosphere
! (photosphere -> tau~0) -- and materializes every quantity the
! profile and pulse writers need, as plain arrays over the extended
! grid. build_stitched_model re-integrates the envelope/atmosphere
! at the converged (Teff, L) with the fixed output step sizes
! (stitch's recipe, formerly io/yrec_output.f90's build_extended),
! then fills:
!   stx_prof(1:n_prof_cols, 1:n_ext)  -- the profile column values
!     (column meanings = io/yrec_output.f90's profile column
!     registry; column 1, the zone number, is left 0 -- the writer
!     numbers zones itself)
!   stx_pulse(1:n_pulse_cols, 1:n_ext) -- the GYRE/FGONG pulse set
! The io writers are pure readers of these arrays; nothing in the
! output path integrates or stitches anything. compute_observables
! themes may read them too (profile-based observables).
module stitched_model_lib
      use phys_const_lib
      use math_lib
      use star_info_lib, only: json
      implicit none
      private
      public :: build_stitched_model, n_ext, n_ie, stx_prof, &
           stx_pulse, n_prof_cols, n_pulse_cols, &
           ip_mass, ip_logR, ip_logT, ip_logRho, ip_logP, ip_conv, &
           ip_gradr, ip_gradT, ip_grada, ip_conv_vel, ip_brunt_N2, &
           ip_csound
      public :: ipul_r, ipul_m, ipul_L, ipul_P, ipul_T, ipul_rho, &
           ipul_grad, ipul_N2, ipul_gamma1, ipul_grada, ipul_delta, &
           ipul_kap, ipul_kap_kap_T, ipul_kap_kap_rho, ipul_eps, &
           ipul_eps_eps_T, ipul_eps_eps_rho, ipul_omega, ipul_cp, &
           ipul_mu_e_inv, ipul_h1, ipul_z, ipul_species_base, &
           ipul_eps_grav, ipul_gyre_last

      integer, parameter :: n_prof_cols = 60
      integer, parameter :: n_pulse_cols = 35
! Named pulse-column indices (2026 audit, section 8: the pts(35,*)
! map was bare numbers in every writer -- reverse-engineered once
! per bug fix). Columns 1..18 are the GYRE schema-101 set, in file
! order (ipul_gyre_last marks the end of that block); 19..34 the
! FGONG extras. 23..33 are the 11 FGONG-ordered species at
! ipul_species_base + k, k = 1..11 (see species_slot).
      integer, parameter :: ipul_r = 1, ipul_m = 2, ipul_L = 3, &
           ipul_P = 4, ipul_T = 5, ipul_rho = 6, ipul_grad = 7, &
           ipul_N2 = 8, ipul_gamma1 = 9, ipul_grada = 10, &
           ipul_delta = 11, ipul_kap = 12, ipul_kap_kap_T = 13, &
           ipul_kap_kap_rho = 14, ipul_eps = 15, ipul_eps_eps_T = 16, &
           ipul_eps_eps_rho = 17, ipul_omega = 18, ipul_cp = 19, &
           ipul_mu_e_inv = 20, ipul_h1 = 21, ipul_z = 22, &
           ipul_species_base = 22, ipul_eps_grav = 34
      integer, parameter :: ipul_gyre_last = ipul_omega
! Named indices for the stx_prof columns physics consumers read
! (column meanings = the profile column registry in io/yrec_output).
      integer, parameter :: ip_mass = 2, ip_logR = 3, ip_logT = 4, &
           ip_logRho = 5, ip_logP = 6, ip_conv = 9, ip_gradr = 12, &
           ip_gradT = 13, ip_grada = 14, ip_conv_vel = 15, &
           ip_brunt_N2 = 54, ip_csound = 58
! further columns read (not exported) by compute_seismic_columns:
! gamma1, delta = -dlnrho/dlnT, mean molecular weight; and the four
! seismic columns 54-57 (brunt_N2, lamb_S2, gradL, gradr_div_grada)
! that it fills.
      integer, parameter :: ip_gamma1 = 10, ip_delta = 51, ip_mu = 52
      integer, parameter :: ip_seismic_first = ip_brunt_N2, &
           ip_seismic_last = ip_brunt_N2 + 3
! The extended model: interior (center -> fitting point) + envelope
! (fitting point -> photosphere) + atmosphere (photosphere -> tau~0),
! assembled inward-to-outward, exactly the regions io/write_stitched_profile.f90
! splices for the legacy .store format. Profiles and pulse files both
! cover the full star: truncating at the fitting point would drop the
! superadiabatic layer and photosphere, which dominate p-mode
! frequencies.
      integer, parameter :: max_ext = 3*json
      integer :: n_ext = 0
      integer :: ext_region(max_ext)   ! 1 interior, 2 envelope, 3 atmosphere
      integer :: ext_index(max_ext)    ! index within that region
! Seismic profile columns (54-57), precomputed over the extended grid
! by compute_seismic_columns (called from build_stitched_model): brunt_N2,
! lamb_S2 (l=1), gradL, gradr_div_grada.
      double precision :: ext_seismic(4,max_ext)
! Geometric height of each atmosphere point above the photosphere,
! indexed like atmo_struct (1 = outermost): the cumulative sum of
! envint's per-step lengths, built by build_stitched_model.
      double precision :: atm_height(max_ext)
      double precision :: stx_prof(n_prof_cols, max_ext)
      double precision :: stx_pulse(n_pulse_cols, max_ext)
! Index of the last NON-atmosphere point: the turnover calculation
! walks regions 1-2 only (interior + envelope), matching the
! pre-restructure combined-array assembly.
      integer :: n_ie = 0

contains

! Regenerate the envelope/atmosphere structures for the CONVERGED
! model and build the inward-to-outward index map. The envelope the
! solver last integrated belongs to some trial (Teff, L) -- or was
! never integrated at all, when the envelope triangle interpolated --
! so a fresh atm_get at the converged values is required (the same recipe
! io/write_stitched_profile.f90 uses for the legacy .store profiles, with the
! same fixed output step sizes).
subroutine build_stitched_model
      use star_info_lib, only: star, i_h1, i_metals
      use envstruct_lib
      use atmstruct_lib
      use envint_lib, only: atm_get, envint_step_config, fixed_envint_step
      integer :: j, i, jerr
      type(envint_step_config) :: atm_steps, env_steps
      double precision :: b, gl, rl, ateffl, plim
      integer :: ksaha
      logical :: lprt, lsbc0

! interior always present
      n_ext = 0
      do j = 1, star%nz
         n_ext = n_ext + 1
         ext_region(n_ext) = 1
         ext_index(n_ext) = j
      end do
      if (.not. star%job%calc_envelope_flag) then
         n_ie = n_ext
! interior-only: materialize what we have. (Historical quirk kept
! byte-for-byte: the seismic columns are NOT recomputed on this
! branch -- they hold whatever the last full build left.)
         call fill_stitched_arrays
         return
      end if

! ---- re-integrate at the converged model (stitch's recipe) ----
! fixed output step sizes for this atm_get call only (2026 W2: passed
! in instead of overwriting and restoring star%job%{atm,env}_step_*)
      atm_steps = fixed_envint_step(star%ctrl%atm_step_size)
      env_steps = fixed_envint_step(star%ctrl%envelope_step_size)

      ksaha = 0
      lprt = .false.
      lsbc0 = .false.
      b = exp(ln10*star%log_L)
      rl = 0.5d0*(star%log_L + star%log10_solar_luminosity - 4.0d0*star%log_Teff &
           - c4pil - csigl)
      gl = cgl + star%log_total_mass - rl - rl
      plim = star%logP(star%nz)
      if (star%convective_flag(star%nz) .and. star%ctrl%spot_filling_factor /= 0.0d0 &
          .and. star%ctrl%spot_temp_contrast /= 1.0d0) then
         ateffl = star%log_Teff - 0.25d0*log10(star%ctrl%spot_filling_factor* &
              pow(star%ctrl%spot_temp_contrast, 4.0d0) + 1.0d0 - star%ctrl%spot_filling_factor)
      else
         ateffl = star%log_Teff
      end if
      jerr = 0
! The stitch is the ONE writer of star%pphot: atm_get sets it from
! this integration at the converged model, and no other path
! recomputes it (the former gettau/wrtout own atm_get calls are
! gone). The stitch runs every step, so pphot is always current.
      call atm_get(b, star%fp_rot(star%nz), star%ft_rot(star%nz), gl, &
           star%log_total_mass, lprt, lsbc0, plim, rl, ateffl, &
           star%xa(i_h1,star%nz), star%xa(i_metals,star%nz), ksaha, jerr, &
           atm_steps=atm_steps, env_steps=env_steps)

      if (jerr /= 0) then
         n_ie = n_ext
         call fill_stitched_arrays
         return
      end if

! envelope: env_struct runs fitting point -> photosphere (envint
! inverts it), so it appends directly. Its innermost point repeats
! the fitting point (to re-integration roundoff, so its radius can
! even land marginally BELOW the last interior point's) -- skip any
! envelope point not strictly outside the interior, keeping the
! extended radius strictly monotonic for GYRE.
      do i = 1, env_struct%num_env_points
         if (n_ext >= max_ext) exit
         if (env_struct%env_log10_radius(i) <= star%logR(star%nz)) cycle
         n_ext = n_ext + 1
         ext_region(n_ext) = 2
         ext_index(n_ext) = i
      end do
      n_ie = n_ext
! atmosphere: atmo_struct runs outward-in, so walk it in reverse.
! atmo_delta_depth(k) is the PER-STEP geometric length between tau
! points k-1 and k (envint's Cox p590 quadrature), not a height, so
! accumulate it into atm_height: point i sits sum(delta(i+1:num))
! above the deepest point (tau = 2/3, the photosphere). The legacy
! .store stitch added the raw per-step value to the photosphere
! radius -- a longstanding bug that left the atmosphere radii
! non-monotonic (and made every atmosphere point hug r_phot); it is
! fixed here, not preserved. delta(1), the one value polluted by the
! stale-prev_tau quirk envint notes, never enters any height. The
! deepest point (height 0) repeats the envelope's photosphere radius
! and is skipped, same as the fitting-point repeat above.
      if (atmo_struct%num_atm_points > 0) then
         atm_height(atmo_struct%num_atm_points) = 0.0d0
         do i = atmo_struct%num_atm_points - 1, 1, -1
            atm_height(i) = atm_height(i+1) + &
                 atmo_struct%atmo_delta_depth(i+1)
         end do
      end if
      do i = atmo_struct%num_atm_points, 1, -1
         if (n_ext >= max_ext) exit
         if (atm_height(i) <= 0.0d0) cycle
         n_ext = n_ext + 1
         ext_region(n_ext) = 3
         ext_index(n_ext) = i
      end do
      call compute_seismic_columns
      call fill_stitched_arrays
end subroutine build_stitched_model

! ---------------------------------------------------------------
! Fill ext_seismic(1:4,:) over the assembled extended grid:
!   1 brunt_N2 = g*[(1/Gamma1) dlnP/dr - dlnRho/dr]  (centered
!     differences in r -- the same derivative content as FGONG's A4)
!   2 lamb_S2  = l(l+1) c^2 / r^2 with l=1: 2*Gamma1*P/(rho r^2)
!   3 gradL    = grada + (1/delta)*dln(mu)/dlnP  (Ledoux gradient
!     with the ideal-fully-ionized phi=1 approximation; the mu
!     gradient is a centered difference in lnP)
!   4 gradr_div_grada
! Endpoints copy their neighbor's derivative-based values.
subroutine compute_seismic_columns
      use math_lib
      use star_info_lib, only: star
      integer :: j
      double precision :: r(max_ext), lnp(max_ext), lnrho(max_ext), &
           lnmu(max_ext), g1(max_ext), grada_(max_ext), delta_(max_ext), &
           gradr_(max_ext), mass_g(max_ext)
      double precision :: dr, dlnp_dr, dlnrho_dr, grav, dlnmu_dlnp, p, rho
      logical :: mu_ok(max_ext)

      mu_ok = .true.
      do j = 1, n_ext
         r(j)      = exp(ln10*ext_profile_value(ip_logR, j))
         lnp(j)    = ln10*ext_profile_value(ip_logP, j)
         lnrho(j)  = ln10*ext_profile_value(ip_logRho, j)
         g1(j)     = ext_profile_value(ip_gamma1, j)
         grada_(j) = ext_profile_value(ip_grada, j)
         gradr_(j) = ext_profile_value(ip_gradr, j)
         delta_(j) = ext_profile_value(ip_delta, j)
         mass_g(j) = ext_profile_value(ip_mass, j)*star%solar_mass_cgs
         lnmu(j)   = ext_profile_value(ip_mu, j)
         if (lnmu(j) > 0.0d0) then
            lnmu(j) = log(lnmu(j))
         else
            lnmu(j) = 0.0d0
            mu_ok(j) = .false.
         end if
      end do
      do j = 1, n_ext
         p   = exp(lnp(j))
         rho = exp(lnrho(j))
! lamb_S2 (pointwise)
         if (r(j) > 0.0d0 .and. rho > 0.0d0 .and. g1(j) > 0.0d0) then
            ext_seismic(2,j) = 2.0d0*g1(j)*p/(rho*r(j)*r(j))
         else
            ext_seismic(2,j) = 0.0d0
         end if
! gradr_div_grada (pointwise)
         if (grada_(j) /= 0.0d0) then
            ext_seismic(4,j) = gradr_(j)/grada_(j)
         else
            ext_seismic(4,j) = 0.0d0
         end if
      end do
      do j = 2, n_ext-1
         dr = r(j+1) - r(j-1)
         if (dr /= 0.0d0 .and. r(j) > 0.0d0 .and. g1(j) > 0.0d0) then
            dlnp_dr   = (lnp(j+1) - lnp(j-1))/dr
            dlnrho_dr = (lnrho(j+1) - lnrho(j-1))/dr
            grav = exp(ln10*cgl)*mass_g(j)/(r(j)*r(j))
            ext_seismic(1,j) = grav*(dlnp_dr/g1(j) - dlnrho_dr)
         else
            ext_seismic(1,j) = 0.0d0
         end if
! 2026 (bugsweep sec-11): column 52 is now the real mu everywhere
! (it was R/mu in the interior -- wrong-signed composition term --
! and 0 in the envelope, a log(1e-30) spike at the junction); where
! mu is still unavailable the term is dropped rather than faked.
         if (lnp(j+1) /= lnp(j-1) .and. delta_(j) /= 0.0d0 .and. &
              mu_ok(j-1) .and. mu_ok(j+1)) then
            dlnmu_dlnp = (lnmu(j+1) - lnmu(j-1))/(lnp(j+1) - lnp(j-1))
            ext_seismic(3,j) = grada_(j) + dlnmu_dlnp/delta_(j)
         else
            ext_seismic(3,j) = grada_(j)
         end if
      end do
      if (n_ext >= 2) then
         ext_seismic(1,1) = ext_seismic(1,2)
         ext_seismic(3,1) = ext_seismic(3,2)
         ext_seismic(1,n_ext) = ext_seismic(1,n_ext-1)
         ext_seismic(3,n_ext) = ext_seismic(3,n_ext-1)
      end if
end subroutine compute_seismic_columns

! Profile column value at extended point j. Envelope/atmosphere points
! carry what those integrations store; quantities they do not track
! (per-species abundances beyond X/Z, burning terms, rotation
! internals) are zero, as io/write_stitched_profile.f90 also writes them.
double precision function ext_profile_value(icol, j)
      use math_lib
      use star_info_lib, only: star, i_h1, i_metals
      use envstruct_lib
      use atmstruct_lib
      integer, intent(in) :: icol, j
      integer :: i

      if (icol >= ip_seismic_first .and. icol <= ip_seismic_last) then
         ext_profile_value = ext_seismic(icol-ip_seismic_first+1, j)
         return
      end if
      i = ext_index(j)
      select case (ext_region(j))
      case (1)
         ext_profile_value = profile_value(icol, i)
      case (2)
         select case (icol)
         case (2);  ext_profile_value = exp(ln10*(env_struct%env_log10_mass(i) &
                       + star%log_total_mass))/star%solar_mass_cgs
         case (3);  ext_profile_value = env_struct%env_log10_radius(i)
         case (4);  ext_profile_value = env_struct%env_log10_temperature(i)
         case (5);  ext_profile_value = env_struct%env_log10_density(i)
         case (6);  ext_profile_value = env_struct%env_log10_pressure(i)
         case (7);  ext_profile_value = env_struct%env_luminosity(i)
         case (9)
            if (env_struct%env_convective_flag(i)) then
               ext_profile_value = 1.0d0
            else
               ext_profile_value = 0.0d0
            end if
         case (10); ext_profile_value = env_struct%env_gamma1(i)
         case (11); ext_profile_value = env_struct%env_opacity(i)
         case (12); ext_profile_value = env_struct%env_gradients(1,i)
! env_gradients order is (1) radiative, (2) adiabatic, (3) actual
! (set from current_gradients in envint) -- profile columns are
! 13 = gradT (actual), 14 = grada. The pre-2026 mapping had these
! two swapped in the envelope region (found when the turnover
! walker started reading the stitched columns).
         case (13); ext_profile_value = env_struct%env_gradients(3,i)
         case (14); ext_profile_value = env_struct%env_gradients(2,i)
         case (15); ext_profile_value = env_struct%env_convective_velocity(i)
         case (16); ext_profile_value = env_struct%env_beta(i)
         case (27); ext_profile_value = env_struct%env_hydrogen_fraction(i)
         case (41); ext_profile_value = env_struct%env_metal_fraction(i)
         case (42); ext_profile_value = star%omega(star%nz)
         case (50); ext_profile_value = env_struct%env_specific_heat_cp(i)
         case (51); ext_profile_value = -env_struct%env_dlnrho_dlnt(i)
! 2026 (bugsweep sec-11): mu from the gas-pressure ideal-gas relation
! R/mu = beta*P/(rho*T) -- the same identity eqstat uses for its
! specific gas constant -- so the Ledoux column below sees a
! continuous mu across the interior/envelope junction instead of 0.
         case (52); ext_profile_value = ideal_gas_mu( &
              env_struct%env_log10_pressure(i), &
              env_struct%env_log10_temperature(i), &
              env_struct%env_log10_density(i), env_struct%env_beta(i))
         case (58); ext_profile_value = sqrt(env_struct%env_gamma1(i)* &
              exp(ln10*(env_struct%env_log10_pressure(i) - &
              env_struct%env_log10_density(i))))
         case default; ext_profile_value = 0.0d0
         end select
      case (3)
         select case (icol)
         case (2);  ext_profile_value = star%star_mass
         case (3);  ext_profile_value = log10(exp(ln10* &
                       env_struct%env_log10_radius(env_struct%num_env_points)) &
                       + atm_height(i))
         case (4);  ext_profile_value = atmo_struct%atmo_log10_temperature(i)
         case (5);  ext_profile_value = atmo_struct%atmo_log10_density(i)
         case (6);  ext_profile_value = atmo_struct%atmo_log10_pressure(i)
         case (7);  ext_profile_value = exp(ln10*star%log_L)
         case (10); ext_profile_value = atmo_struct%atmo_gamma1(i)
         case (11); ext_profile_value = atmo_struct%atmo_opacity(i)
         case (12); ext_profile_value = atmo_struct%atmo_gradients(1,i)
         case (13); ext_profile_value = atmo_struct%atmo_gradients(2,i)
         case (14); ext_profile_value = atmo_struct%atmo_gradients(3,i)
         case (16); ext_profile_value = atmo_struct%atmo_beta(i)
         case (27); ext_profile_value = star%xa(i_h1,star%nz)
         case (41); ext_profile_value = star%xa(i_metals,star%nz)
         case (42); ext_profile_value = star%omega(star%nz)
         case (50); ext_profile_value = atmo_struct%atmo_specific_heat_cp(i)
         case (51); ext_profile_value = -atmo_struct%atmo_dlnrho_dlnt(i)
         case (52); ext_profile_value = ideal_gas_mu( &
              atmo_struct%atmo_log10_pressure(i), &
              atmo_struct%atmo_log10_temperature(i), &
              atmo_struct%atmo_log10_density(i), atmo_struct%atmo_beta(i))
         case (58); ext_profile_value = sqrt(atmo_struct%atmo_gamma1(i)* &
              exp(ln10*(atmo_struct%atmo_log10_pressure(i) - &
              atmo_struct%atmo_log10_density(i))))
         case default; ext_profile_value = 0.0d0
         end select
      case default
         ext_profile_value = 0.0d0
      end select
end function ext_profile_value

! mu = R*rho*T/(beta*P): the ideal-gas mean molecular weight from the
! gas pressure (eqstat's own specific-gas-constant identity). Zero
! when the inputs cannot support it.
double precision function ideal_gas_mu(log10_p, log10_t, log10_rho, beta)
      use math_lib
      use phys_const_lib, only: gas_constant
      double precision, intent(in) :: log10_p, log10_t, log10_rho, beta
      if (beta > 0.0d0) then
         ideal_gas_mu = gas_constant*exp(ln10*(log10_rho + log10_t - log10_p))/beta
      else
         ideal_gas_mu = 0.0d0
      end if
end function ideal_gas_mu

! Per-zone value for profile column icol at YREC zone index k
! (1 = center .. nz = surface). Sources match putstore's per-shell
! block and the pulse arrays coefft fills every model.
double precision function profile_value(icol, k)
      use math_lib
      use star_info_lib, only: star, i_be9, i_c12, i_c13, i_eps_cno, i_eps_grav, i_eps_he3, i_eps_neu, i_eps_pp1, i_eps_pp2, i_eps_pp3, i_h1, i_h2, i_he3, i_he4, i_li6, i_li7, i_metals, i_n14, i_n15, i_o16, i_o17, i_o18
      integer, intent(in) :: icol, k

      select case (icol)
      case (1);  profile_value = 0.0d0   ! zone number set by the writer
      case (2);  profile_value = star%m(k)/star%solar_mass_cgs
      case (3);  profile_value = star%logR(k)
      case (4);  profile_value = star%logT(k)
      case (5);  profile_value = star%logRho(k)
      case (6);  profile_value = star%logP(k)
      case (7);  profile_value = star%luminosity_lsun(k)
      case (8);  profile_value = star%dm(k)/star%solar_mass_cgs
      case (9)
         if (star%convective_flag(k)) then
            profile_value = 1.0d0
         else
            profile_value = 0.0d0
         end if
      case (10); profile_value = star%adiabatic_index_gamma1(k)
      case (11); profile_value = star%opacity_zone(k)
      case (12); profile_value = star%gradr(k)
      case (13); profile_value = star%gradT(k)
      case (14); profile_value = star%grada(k)
      case (15); profile_value = star%conv_vel(k)
      case (16); profile_value = star%beta(k)
      case (17); profile_value = star%eta(k)
      case (18)
         profile_value = exp(ln10*(cgl - 2.0d0*star%logR(k)))*star%m(k)
      case (19); profile_value = star%eps_channels(i_eps_pp1,k)
      case (20); profile_value = star%eps_channels(i_eps_pp2,k)
      case (21); profile_value = star%eps_channels(i_eps_pp3,k)
      case (22); profile_value = star%eps_channels(i_eps_cno,k)
      case (23); profile_value = star%eps_channels(i_eps_he3,k)
      case (24); profile_value = star%eps_total(k)
      case (25); profile_value = star%eps_channels(i_eps_neu,k)
      case (26); profile_value = star%eps_channels(i_eps_grav,k)
      case (27); profile_value = star%xa(i_h1,k)
      case (28); profile_value = star%xa(i_h2,k)
      case (29); profile_value = star%xa(i_he3,k)
      case (30); profile_value = star%xa(i_he4,k)
      case (31); profile_value = star%xa(i_li6,k)
      case (32); profile_value = star%xa(i_li7,k)
      case (33); profile_value = star%xa(i_be9,k)
      case (34); profile_value = star%xa(i_c12,k)
      case (35); profile_value = star%xa(i_c13,k)
      case (36); profile_value = star%xa(i_n14,k)
      case (37); profile_value = star%xa(i_n15,k)
      case (38); profile_value = star%xa(i_o16,k)
      case (39); profile_value = star%xa(i_o17,k)
      case (40); profile_value = star%xa(i_o18,k)
      case (41); profile_value = star%xa(i_metals,k)
      case (42); profile_value = star%omega(k)
      case (43); profile_value = star%j_rot(k)
      case (44); profile_value = star%i_rot(k)
      case (45); profile_value = star%fp_rot(k)
      case (46); profile_value = star%ft_rot(k)
      case (47); profile_value = star%es_circulation_velocity(k)
      case (48); profile_value = star%gsf_circulation_velocity(k)
      case (49); profile_value = star%secular_shear_velocity(k)
      case (50); profile_value = star%pulse_specific_heat(k)
      case (51); profile_value = -star%pulse_dlnrho_dlnt(k)
      case (52); profile_value = star%pulse_mean_molecular_weight(k)
! 2026 (bugsweep sec-11): the array already holds 1/mu_e (eos
! electron_mean_weight_inverse); it used to be inverted again here.
      case (53); profile_value = star%pulse_electron_mean_weight_inverse(k)
      case (58)
! sound speed sqrt(Gamma1*P/rho) [cm/s]
         profile_value = sqrt(star%adiabatic_index_gamma1(k)* &
              exp(ln10*(star%logP(k) - star%logRho(k))))
      case (59); profile_value = star%am_diffusion_coeff(k)
      case (60); profile_value = star%mixing_diffusion_coeff(k)
      case default
         profile_value = 0.0d0
      end select
end function profile_value

! Assemble the per-point pulse data for the FULL extended model
! (interior + envelope + atmosphere; build_stitched_model must have run).
! Columns 1-18 are the GYRE schema-101 set; 19-34 the extras FGONG
! needs. Envelope/atmosphere points carry what those integrations
! store; opacity/energy derivative columns they do not track are
! zero there (they matter only for nonadiabatic work), and the
! composition above the fitting point is the surface composition, as
! io/write_stitched_profile.f90 also writes.
subroutine build_pulse_points(pts)
      use math_lib
      use star_info_lib, only: star, i_eps_grav, i_h1, i_metals
      use envstruct_lib
      use atmstruct_lib
      double precision, intent(out) :: pts(n_pulse_cols, n_ext)
      integer :: j, i, k
      double precision :: r, m, P, T, rho, delta, nab, nab_ad, grav

      do j = 1, n_ext
         i = ext_index(j)
         pts(:,j) = 0.0d0
         select case (ext_region(j))
         case (1)
            r = exp(ln10*star%logR(i))
            m = star%m(i)
            P = exp(ln10*star%logP(i))
            T = exp(ln10*star%logT(i))
            rho = exp(ln10*star%logRho(i))
            delta = -star%pulse_dlnrho_dlnt(i)
            nab = star%gradT(i)
            nab_ad = star%grada(i)
            pts(ipul_L,j) = star%luminosity_lsun(i)*star%solar_luminosity_cgs
            pts(ipul_gamma1,j) = star%adiabatic_index_gamma1(i)
            pts(ipul_kap,j) = star%opacity_zone(i)
! GSM/GYRE convention (MESA pulse_gyre.f90): kap_kap_T = kap*dlnkap/dlnT
! (the absolute derivative dkap/dlnT), NOT the bare log-derivative;
! likewise eps_eps_T = eps*dlneps/dlnT. Bare log-derivatives were
! written here originally -- wrong by factors kap and eps (only
! nonadiabatic GYRE runs read these columns).
            pts(ipul_kap_kap_T,j) = star%opacity_zone(i)*star%pulse_dlnkap_dlnt(i)
            pts(ipul_kap_kap_rho,j) = star%opacity_zone(i)*star%pulse_dlnkap_dlnrho(i)
            pts(ipul_eps,j) = star%eps_total(i)
            pts(ipul_eps_eps_T,j) = star%eps_total(i)*star%pulse_dlneps_dlnt(i)
            pts(ipul_eps_eps_rho,j) = star%eps_total(i)*star%pulse_dlneps_dlnrho(i)
            pts(ipul_omega,j) = star%omega(i)
            pts(ipul_cp,j) = star%pulse_specific_heat(i)
! 2026 (bugsweep sec-11): 1/mu_e is what the array holds -- FGONG
! var(14) wants exactly that (it used to be inverted a second time).
            pts(ipul_mu_e_inv,j) = star%pulse_electron_mean_weight_inverse(i)
            pts(ipul_h1,j) = star%xa(i_h1,i)
            pts(ipul_z,j) = star%xa(i_metals,i)
            do k = 1, 11
               pts(ipul_species_base+k,j) = star%xa(species_slot(k),i)
            end do
            pts(ipul_eps_grav,j) = star%eps_channels(i_eps_grav,i)
         case (2)
            r = exp(ln10*env_struct%env_log10_radius(i))
            m = exp(ln10*(env_struct%env_log10_mass(i) + &
                 star%log_total_mass))
            P = exp(ln10*env_struct%env_log10_pressure(i))
            T = exp(ln10*env_struct%env_log10_temperature(i))
            rho = exp(ln10*env_struct%env_log10_density(i))
            delta = -env_struct%env_dlnrho_dlnt(i)
! 2026 sweep fix: env_gradients order is (1) radiative, (2) ADIABATIC,
! (3) ACTUAL (see the profile-column mapping above and
! envelope_derivs' current_gradients stores) -- unlike atmo_gradients,
! whose order is (rad, actual, adiabatic). This block had the two
! slots swapped (the pulse-side twin of the pre-2026 profile-column
! swap), flipping the thermal N^2 sign across the envelope region of
! every pulse file and swapping the grad/grad_ad pulse columns there.
            nab = env_struct%env_gradients(3,i)
            nab_ad = env_struct%env_gradients(2,i)
            pts(ipul_L,j) = env_struct%env_luminosity(i)*star%solar_luminosity_cgs
            pts(ipul_gamma1,j) = env_struct%env_gamma1(i)
            pts(ipul_kap,j) = env_struct%env_opacity(i)
            pts(ipul_omega,j) = star%omega(star%nz)
            pts(ipul_cp,j) = env_struct%env_specific_heat_cp(i)
            pts(ipul_h1,j) = env_struct%env_hydrogen_fraction(i)
            pts(ipul_z,j) = env_struct%env_metal_fraction(i)
            do k = 1, 11
               pts(ipul_species_base+k,j) = star%xa(species_slot(k),star%nz)
            end do
         case default   ! atmosphere
            r = exp(ln10* &
                 env_struct%env_log10_radius(env_struct%num_env_points)) &
                 + atm_height(i)
            m = exp(ln10*star%log_total_mass)
            P = exp(ln10*atmo_struct%atmo_log10_pressure(i))
            T = exp(ln10*atmo_struct%atmo_log10_temperature(i))
            rho = exp(ln10*atmo_struct%atmo_log10_density(i))
            delta = -atmo_struct%atmo_dlnrho_dlnt(i)
            nab = atmo_struct%atmo_gradients(2,i)
            nab_ad = atmo_struct%atmo_gradients(3,i)
            pts(ipul_L,j) = exp(ln10*star%log_L)*star%solar_luminosity_cgs
            pts(ipul_gamma1,j) = atmo_struct%atmo_gamma1(i)
            pts(ipul_kap,j) = atmo_struct%atmo_opacity(i)
            pts(ipul_omega,j) = star%omega(star%nz)
            pts(ipul_cp,j) = atmo_struct%atmo_specific_heat_cp(i)
            pts(ipul_h1,j) = star%xa(i_h1,star%nz)
            pts(ipul_z,j) = star%xa(i_metals,star%nz)
            do k = 1, 11
               pts(ipul_species_base+k,j) = star%xa(species_slot(k),star%nz)
            end do
         end select
         pts(ipul_r,j) = r
         pts(ipul_m,j) = m
         pts(ipul_P,j) = P
         pts(ipul_T,j) = T
         pts(ipul_rho,j) = rho
         pts(ipul_grad,j) = nab
         pts(ipul_grada,j) = nab_ad
         pts(ipul_delta,j) = delta
         if (r > 0.0d0) then
            grav = exp(ln10*cgl)*m/(r*r)
            pts(ipul_N2,j) = grav*grav*(rho/P)*delta*(nab_ad - nab)
         end if
      end do

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
      do j = 2, n_ext - 1
         if (pts(ipul_N2,j) >= 0.0d0 .and. &
              pts(ipul_r,j) > 0.0d0 .and. pts(ipul_r,j+1) > pts(ipul_r,j-1) .and. &
              pts(ipul_gamma1,j) > 0.0d0) then
            grav = exp(ln10*cgl)*pts(ipul_m,j)/(pts(ipul_r,j)*pts(ipul_r,j))
            pts(ipul_N2,j) = grav*( &
                 (log(pts(ipul_P,j+1)) - log(pts(ipul_P,j-1)))/pts(ipul_gamma1,j) - &
                 (log(pts(ipul_rho,j+1)) - log(pts(ipul_rho,j-1))) ) / &
                 (pts(ipul_r,j+1) - pts(ipul_r,j-1))
         end if
      end do
      if (n_ext >= 2) then
         pts(ipul_N2,1) = pts(ipul_N2,2)
         pts(ipul_N2,n_ext) = pts(ipul_N2,n_ext-1)
      end if
end subroutine build_pulse_points

! star%xa slot for pulse column 22+k (k = 1..11), in FGONG species
! order (column 34 is eps_grav, not a species; be9 has no FGONG
! column).
integer function species_slot(k)
      use star_info_lib, only: i_he3, i_c12, i_c13, i_n14, i_o16, i_h2, i_he4, i_li7, i_n15, i_o17, i_o18
      integer, intent(in) :: k
      integer, parameter :: slots(11) = [i_he3, i_c12, i_c13, i_n14, &
           i_o16, i_h2, i_he4, i_li7, i_n15, i_o17, i_o18]
      species_slot = slots(k)
end function species_slot

! ---------------------------------------------------------------
! Materialize the writer-facing arrays from the assembled grid.
subroutine fill_stitched_arrays
      integer :: j, icol
      do j = 1, n_ext
         do icol = 1, n_prof_cols
            stx_prof(icol, j) = ext_profile_value(icol, j)
         end do
      end do
      call build_pulse_points(stx_pulse)
end subroutine fill_stitched_arrays

end module stitched_model_lib
