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
      implicit none
      private
      public :: build_stitched_model, n_ext, n_ie, stx_prof, &
           stx_pulse, n_prof_cols, n_pulse_cols, &
           ip_mass, ip_logR, ip_logT, ip_logRho, ip_logP, ip_conv, &
           ip_gradr, ip_gradT, ip_grada, ip_conv_vel, ip_brunt_N2, &
           ip_csound

      integer, parameter :: n_prof_cols = 60
      integer, parameter :: n_pulse_cols = 35
! Named indices for the stx_prof columns physics consumers read
! (column meanings = the profile column registry in io/yrec_output).
      integer, parameter :: ip_mass = 2, ip_logR = 3, ip_logT = 4, &
           ip_logRho = 5, ip_logP = 6, ip_conv = 9, ip_gradr = 12, &
           ip_gradT = 13, ip_grada = 14, ip_conv_vel = 15, &
           ip_brunt_N2 = 54, ip_csound = 58
! The extended model: interior (center -> fitting point) + envelope
! (fitting point -> photosphere) + atmosphere (photosphere -> tau~0),
! assembled inward-to-outward, exactly the regions io/write_stitched_profile.f90
! splices for the legacy .store format. Profiles and pulse files both
! cover the full star: truncating at the fitting point would drop the
! superadiabatic layer and photosphere, which dominate p-mode
! frequencies.
      integer, parameter :: max_ext = 3*5000
      integer :: n_ext = 0
      integer :: ext_region(max_ext)   ! 1 interior, 2 envelope, 3 atmosphere
      integer :: ext_index(max_ext)    ! index within that region
! Seismic profile columns (54-57), precomputed over the extended grid
! by compute_seismic_columns (called from build_extended): brunt_N2,
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
      use envint_lib, only: atm_get
      integer :: j, i, jerr
      double precision :: atm_beg0, atm_min0, atm_max0
      double precision :: env_beg0, env_min0, env_max0
      double precision :: b, gl, rl, ateffl, plim, dum1(4), dum2(3), &
           dum3(3), dum4(3)
      integer :: ixx, idum, katm, kenv, ksaha
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
      atm_beg0 = star%job%atm_step_begin
      atm_min0 = star%job%atm_step_min
      atm_max0 = star%job%atm_step_max
      env_beg0 = star%job%env_step_begin
      env_min0 = star%job%env_step_min
      env_max0 = star%job%env_step_max
      star%job%atm_step_begin = star%ctrl%atm_step_size
      star%job%atm_step_min = star%ctrl%atm_step_size
      star%job%atm_step_max = star%ctrl%atm_step_size
      star%job%env_step_begin = star%ctrl%envelope_step_size
      star%job%env_step_min = star%ctrl%envelope_step_size
      star%job%env_step_max = star%ctrl%envelope_step_size

      idum = 0
      ixx = 0
      katm = 0
      kenv = 0
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
           star%log_total_mass, ixx, lprt, lsbc0, plim, rl, ateffl, &
           star%xa(i_h1,star%nz), star%xa(i_metals,star%nz), dum1, idum, katm, &
           kenv, ksaha, dum2, dum3, dum4, jerr)

      star%job%atm_step_begin = atm_beg0
      star%job%atm_step_min = atm_min0
      star%job%atm_step_max = atm_max0
      star%job%env_step_begin = env_beg0
      star%job%env_step_min = env_min0
      star%job%env_step_max = env_max0
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

      do j = 1, n_ext
         r(j)      = exp(ln10*ext_profile_value(3, j))
         lnp(j)    = ln10*ext_profile_value(6, j)
         lnrho(j)  = ln10*ext_profile_value(5, j)
         g1(j)     = ext_profile_value(10, j)
         grada_(j) = ext_profile_value(14, j)
         gradr_(j) = ext_profile_value(12, j)
         delta_(j) = ext_profile_value(51, j)
         mass_g(j) = ext_profile_value(2, j)*star%solar_mass_cgs
         lnmu(j)   = log(max(ext_profile_value(52, j), 1.0d-30))
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
         if (lnp(j+1) /= lnp(j-1) .and. delta_(j) /= 0.0d0) then
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

      if (icol >= 54 .and. icol <= 57) then
         ext_profile_value = ext_seismic(icol-53, j)
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
         case (27); ext_profile_value = env_struct%env_hydrogen_fraction(i)
         case (41); ext_profile_value = env_struct%env_metal_fraction(i)
         case (42); ext_profile_value = star%omega(star%nz)
         case (50); ext_profile_value = env_struct%env_specific_heat_cp(i)
         case (51); ext_profile_value = -env_struct%env_dlnrho_dlnt(i)
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
         case (58); ext_profile_value = sqrt(atmo_struct%atmo_gamma1(i)* &
              exp(ln10*(atmo_struct%atmo_log10_pressure(i) - &
              atmo_struct%atmo_log10_density(i))))
         case default; ext_profile_value = 0.0d0
         end select
      case default
         ext_profile_value = 0.0d0
      end select
end function ext_profile_value

! Per-zone value for profile column icol at YREC zone index k
! (1 = center .. nz = surface). Sources match putstore's per-shell
! block and the pulse arrays coefft fills every model.
double precision function profile_value(icol, k)
      use math_lib
      use star_info_lib, only: star, i_be9, i_c12, i_c13, i_eps_cno, i_eps_grav, i_eps_he3, i_eps_neu, i_eps_pp1, i_eps_pp2, i_eps_pp3, i_grad_actual, i_grad_ad, i_grad_rad, i_h1, i_h2, i_he3, i_he4, i_li6, i_li7, i_metals, i_n14, i_n15, i_o16, i_o17, i_o18
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
      case (53)
         if (star%pulse_electron_mean_molecular_weight(k) > 0.0d0) then
            profile_value = &
                 1.0d0/star%pulse_electron_mean_molecular_weight(k)
         else
            profile_value = 0.0d0
         end if
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
! (interior + envelope + atmosphere; build_extended must have run).
! Columns 1-18 are the GYRE schema-101 set; 19-34 the extras FGONG
! needs. Envelope/atmosphere points carry what those integrations
! store; opacity/energy derivative columns they do not track are
! zero there (they matter only for nonadiabatic work), and the
! composition above the fitting point is the surface composition, as
! io/write_stitched_profile.f90 also writes.
subroutine build_pulse_points(pts)
      use math_lib
      use star_info_lib, only: star, i_eps_grav, i_grad_actual, i_grad_ad, i_h1, i_metals
      use envstruct_lib
      use atmstruct_lib
      double precision, intent(out) :: pts(35, n_ext)
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
            pts(3,j) = star%luminosity_lsun(i)*star%solar_luminosity_cgs
            pts(9,j) = star%adiabatic_index_gamma1(i)
            pts(12,j) = star%opacity_zone(i)
! GSM/GYRE convention (MESA pulse_gyre.f90): kap_kap_T = kap*dlnkap/dlnT
! (the absolute derivative dkap/dlnT), NOT the bare log-derivative;
! likewise eps_eps_T = eps*dlneps/dlnT. Bare log-derivatives were
! written here originally -- wrong by factors kap and eps (only
! nonadiabatic GYRE runs read these columns).
            pts(13,j) = star%opacity_zone(i)*star%pulse_dlnkap_dlnt(i)
            pts(14,j) = star%opacity_zone(i)*star%pulse_dlnkap_dlnrho(i)
            pts(15,j) = star%eps_total(i)
            pts(16,j) = star%eps_total(i)*star%pulse_dlneps_dlnt(i)
            pts(17,j) = star%eps_total(i)*star%pulse_dlneps_dlnrho(i)
            pts(18,j) = star%omega(i)
            pts(19,j) = star%pulse_specific_heat(i)
            if (star%pulse_electron_mean_molecular_weight(i) &
                 > 0.0d0) pts(20,j) = &
                 1.0d0/star%pulse_electron_mean_molecular_weight(i)
            pts(21,j) = star%xa(i_h1,i)
            pts(22,j) = star%xa(i_metals,i)
            do k = 1, 11
               pts(22+k,j) = star%xa(species_slot(k),i)
            end do
            pts(34,j) = star%eps_channels(i_eps_grav,i)
         case (2)
            r = exp(ln10*env_struct%env_log10_radius(i))
            m = exp(ln10*(env_struct%env_log10_mass(i) + &
                 star%log_total_mass))
            P = exp(ln10*env_struct%env_log10_pressure(i))
            T = exp(ln10*env_struct%env_log10_temperature(i))
            rho = exp(ln10*env_struct%env_log10_density(i))
            delta = -env_struct%env_dlnrho_dlnt(i)
            nab = env_struct%env_gradients(2,i)
            nab_ad = env_struct%env_gradients(3,i)
            pts(3,j) = env_struct%env_luminosity(i)*star%solar_luminosity_cgs
            pts(9,j) = env_struct%env_gamma1(i)
            pts(12,j) = env_struct%env_opacity(i)
            pts(18,j) = star%omega(star%nz)
            pts(19,j) = env_struct%env_specific_heat_cp(i)
            pts(21,j) = env_struct%env_hydrogen_fraction(i)
            pts(22,j) = env_struct%env_metal_fraction(i)
            do k = 1, 11
               pts(22+k,j) = star%xa(species_slot(k),star%nz)
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
            pts(3,j) = exp(ln10*star%log_L)*star%solar_luminosity_cgs
            pts(9,j) = atmo_struct%atmo_gamma1(i)
            pts(12,j) = atmo_struct%atmo_opacity(i)
            pts(18,j) = star%omega(star%nz)
            pts(19,j) = atmo_struct%atmo_specific_heat_cp(i)
            pts(21,j) = star%xa(i_h1,star%nz)
            pts(22,j) = star%xa(i_metals,star%nz)
            do k = 1, 11
               pts(22+k,j) = star%xa(species_slot(k),star%nz)
            end do
         end select
         pts(1,j) = r
         pts(2,j) = m
         pts(4,j) = P
         pts(5,j) = T
         pts(6,j) = rho
         pts(7,j) = nab
         pts(10,j) = nab_ad
         pts(11,j) = delta
         if (r > 0.0d0) then
            grav = exp(ln10*cgl)*m/(r*r)
            pts(8,j) = grav*grav*(rho/P)*delta*(nab_ad - nab)
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
         if (pts(8,j) >= 0.0d0 .and. &
              pts(1,j) > 0.0d0 .and. pts(1,j+1) > pts(1,j-1) .and. &
              pts(9,j) > 0.0d0) then
            grav = exp(ln10*cgl)*pts(2,j)/(pts(1,j)*pts(1,j))
            pts(8,j) = grav*( &
                 (log(pts(4,j+1)) - log(pts(4,j-1)))/pts(9,j) - &
                 (log(pts(6,j+1)) - log(pts(6,j-1))) ) / &
                 (pts(1,j+1) - pts(1,j-1))
         end if
      end do
      if (n_ext >= 2) then
         pts(8,1) = pts(8,2)
         pts(8,n_ext) = pts(8,n_ext-1)
      end if
end subroutine build_pulse_points

! star%xa slot for pulse column 22+k (k = 1..11), in FGONG species
! order (column 34 is eps_grav, not a species; be9 has no FGONG
! column).
integer function species_slot(k)
      use star_info_lib, only: star, i_he3, i_c12, i_c13, i_n14, i_o16, i_h2, i_he4, i_li7, i_n15, i_o17, i_o18
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
