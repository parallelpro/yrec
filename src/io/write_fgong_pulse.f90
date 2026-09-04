!----------------------------------------------------------------------
! write_fgong_pulse
!----------------------------------------------------------------------
! New (2026, MESA-style output): writes the converged model in FGONG
! format (ivers 300: ICONST=15 globals, IVAR=40 columns, 5 values per
! line in 1P5E16.9). Column definitions and layout ported from
! mesa-26.04.1/star/private/pulse_fgong.f90 rather than reinvented;
! quantities YREC does not track per zone (the central d2P/d2rho
! terms, dlnGamma1 partials, X_Be7, X_Ne20, r_X) are written as zero,
! as MESA itself does for several of them. Points run surface to
! center (FGONG convention).
!
! Takes the extended point set assembled by yrec_output's
! build_pulse_points -- interior + envelope + atmosphere, center (1)
! to outermost atmosphere point (n) -- so the model reaches the top
! of the atmosphere rather than stopping at the fitting point. The
! global M/R/L refer to the photosphere; atmosphere points sit above
! R_star, so their var(18) = R_star - r is negative (height above
! the photosphere), matching MESA's add_atmosphere convention.
! Column schema of pts(35,n) is documented at build_pulse_points.
subroutine write_fgong_pulse(n, pts, mstar_g, rstar_cm, lstar_cgs, &
     pulse_path)
      use star_info_lib, only: star
      use phys_const_lib
      use math_lib
      use stitched_model_lib, only: n_pulse_cols, ipul_r, ipul_m, &
           ipul_L, ipul_P, ipul_T, ipul_rho, ipul_grad, ipul_N2, &
           ipul_gamma1, ipul_grada, ipul_delta, ipul_kap, ipul_eps, &
           ipul_cp, ipul_mu_e_inv, ipul_h1, ipul_z, &
           ipul_species_base, ipul_eps_grav
      implicit none

      integer, intent(in) :: n
      double precision, intent(in) :: pts(n_pulse_cols, n)
      double precision, intent(in) :: mstar_g, rstar_cm, lstar_cgs
      character(len=*), intent(in) :: pulse_path

      integer, parameter :: iconst = 15, ivar = 40, ivers = 1300
      double precision :: glob(iconst), var(ivar)
      integer :: u, j, k, i
      double precision :: radius_cm, mass_g, grav

      open(newunit=u, file=pulse_path, status='REPLACE', form='FORMATTED')

      glob = 0.0d0
      glob(1) = mstar_g
      glob(2) = rstar_cm
      glob(3) = lstar_cgs
      glob(4) = pts(ipul_z, n)
      glob(5) = pts(ipul_h1, n)
! glob(6) (mixing length alpha), glob(11)/glob(12) (central second
! derivatives) not tracked here -- zero, as tools tolerate.
      glob(13) = star%dage*1.0d9
      glob(14) = exp10(star%log_Teff)
      glob(15) = exp(ln10*cgl)

! MESA-style 4-line comment header (FGONG spec: lines 1-4 free text),
! followed by the standard `nn iconst ivar ivers` record. ivers 1300
! selects the wide E26.18E3 layout, exactly as MESA writes it.
      write(u,'(a)') 'FGONG file'
      write(u,'(a)') 'Created by YREC'
      write(u,'(a,i8,a,1pe16.9,a,0pf10.6)') 'model_number ', &
           star%model_number, '  star_age_yr ', star%dage*1.0d9, &
           '  log_Teff ', star%log_Teff
      write(u,'(a)') ''
      write(u,'(4I10)') n, iconst, ivar, ivers
      write(u,'(1P,5(X,E26.18E3))') (glob(i), i = 1, iconst)

      do j = 1, n
         k = n - j + 1
         radius_cm = pts(ipul_r,k)
         mass_g = pts(ipul_m,k)

         var = 0.0d0
         var(1) = radius_cm
         if (mass_g > 0.0d0) then
            var(2) = log(mass_g/mstar_g)
         else
            var(2) = -1.0d99
         end if
         var(3) = pts(ipul_T,k)
         var(4) = pts(ipul_P,k)
         var(5) = pts(ipul_rho,k)
         var(6) = pts(ipul_h1,k)
         var(7) = pts(ipul_L,k)
         var(8) = pts(ipul_kap,k)
         var(9) = pts(ipul_eps,k)
         var(10) = pts(ipul_gamma1,k)
         var(11) = pts(ipul_grada,k)
         var(12) = pts(ipul_delta,k)
         var(13) = pts(ipul_cp,k)
         var(14) = pts(ipul_mu_e_inv,k)
         if (radius_cm > 0.0d0) then
            grav = exp(ln10*cgl)*mass_g/(radius_cm*radius_cm)
            var(15) = pts(ipul_N2,k)*radius_cm/grav
         end if
! var(16) r_X: not tracked -> 0
         var(17) = pts(ipul_z,k)
         var(18) = rstar_cm - radius_cm
         var(19) = pts(ipul_eps_grav,k)
! var(20) unused
         var(21) = pts(ipul_species_base+1,k)
         var(22) = pts(ipul_species_base+2,k)
         var(23) = pts(ipul_species_base+3,k)
         var(24) = pts(ipul_species_base+4,k)
         var(25) = pts(ipul_species_base+5,k)
! var(26)-var(28) dlnGamma1 partials: not tracked -> 0
         var(29) = pts(ipul_species_base+6,k)
         var(30) = pts(ipul_species_base+7,k)
         var(31) = pts(ipul_species_base+8,k)
! var(32) X_Be7: YREC tracks Be9, not Be7 -> 0
         var(33) = pts(ipul_species_base+9,k)
         var(34) = pts(ipul_species_base+10,k)
         var(35) = pts(ipul_species_base+11,k)
! var(36) X_Ne20 and var(37)-var(40): not tracked -> 0

         write(u,'(1P,5(X,E26.18E3))') (var(i), i = 1, ivar)
      end do

      close(u)
      return
end subroutine write_fgong_pulse
