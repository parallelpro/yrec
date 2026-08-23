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
! center (FGONG convention); YREC stores center (1) .. surface (nz).
!
! A PURE READER over star_info: the pulse arrays (coefft fills them
! unconditionally each model) and the diagnostics arrays supply
! kappa, eps, Gamma1, delta, cp, 1/mu_e and the gradients; the
! Brunt-Vaisala A* uses the same N^2 formula as write_gyre_pulse.
subroutine write_fgong_pulse(num_shells, model_number, pulse_path)
      use star_info_lib, only: star
      use const_lib
      implicit none

      integer, intent(in) :: num_shells, model_number
      character(len=*), intent(in) :: pulse_path

      integer, parameter :: iconst = 15, ivar = 40, ivers = 300
      double precision :: glob(iconst), var(ivar)
      integer :: u, j, k, i
      double precision :: r_outer, m_outer, radius_cm, mass_g, &
           pressure_cgs, temperature_k, density_cgs, delta, grav, &
           brunt_n2, nabla_ad, nabla

      open(newunit=u, file=pulse_path, status='REPLACE', form='FORMATTED')

      m_outer = star%m(num_shells)
      r_outer = exp(ln10*star%logR(num_shells))

      glob = 0.0d0
      glob(1) = m_outer
      glob(2) = r_outer
      glob(3) = star%luminosity_lsun(num_shells)*solar_luminosity_cgs
      glob(4) = star%xa(3,num_shells)
      glob(5) = star%xa(1,num_shells)
! glob(6) (mixing length alpha), glob(11)/glob(12) (central second
! derivatives) not tracked here -- zero, as tools tolerate.
      glob(13) = star%run%dage*1.0d9
      glob(14) = 10.0d0**star%log_Teff
      glob(15) = exp(ln10*cgl)

      write(u,'(a)') 'YREC FGONG output (ported layout: MESA pulse_fgong)'
      write(u,'(a,i8)') 'model_number ', model_number
      write(u,'(a,1pe16.9,a,0pf10.6)') 'star_age_yr ', star%run%dage*1.0d9, &
           '  log_Teff ', star%log_Teff
      write(u,'(a)') ''
      write(u,'(4I10)') num_shells, iconst, ivar, ivers
      write(u,'(1P5E16.9,x)') (glob(i), i = 1, iconst)

      do j = 1, num_shells
         k = num_shells - j + 1
         radius_cm = exp(ln10*star%logR(k))
         mass_g = star%m(k)
         pressure_cgs = exp(ln10*star%logP(k))
         temperature_k = exp(ln10*star%logT(k))
         density_cgs = exp(ln10*star%logRho(k))
         delta = -star%pulse%pulse_dlnrho_dlnt(k)
         nabla_ad = star%diag%del_grad(3,k)
         nabla = star%diag%del_grad(2,k)

         var = 0.0d0
         var(1) = radius_cm
         if (mass_g > 0.0d0) then
            var(2) = log(mass_g/m_outer)
         else
            var(2) = -1.0d99
         end if
         var(3) = temperature_k
         var(4) = pressure_cgs
         var(5) = density_cgs
         var(6) = star%xa(1,k)
         var(7) = star%luminosity_lsun(k)*solar_luminosity_cgs
         var(8) = star%diag%so(k)
         var(9) = star%diag%sesum(k)
         var(10) = star%run%adiabatic_index_gamma1(k)
         var(11) = nabla_ad
         var(12) = delta
         var(13) = star%pulse%pulse_specific_heat(k)
         if (star%pulse%pulse_electron_mean_molecular_weight(k) > 0.0d0) then
            var(14) = 1.0d0/star%pulse%pulse_electron_mean_molecular_weight(k)
         end if
         if (radius_cm > 0.0d0) then
            grav = exp(ln10*cgl)*mass_g/(radius_cm*radius_cm)
            brunt_n2 = grav*grav*(density_cgs/pressure_cgs)*delta* &
                 (nabla_ad - nabla)
            var(15) = brunt_n2*radius_cm/grav
         end if
! var(16) r_X: not tracked -> 0
         var(17) = star%xa(3,k)
         var(18) = r_outer - radius_cm
         var(19) = star%diag%seg(7,k)
! var(20) unused
         var(21) = star%xa(4,k)
         var(22) = star%xa(5,k)
         var(23) = star%xa(6,k)
         var(24) = star%xa(7,k)
         var(25) = star%xa(9,k)
! var(26)-var(28) dlnGamma1 partials: not tracked -> 0
         var(29) = star%xa(12,k)
         var(30) = star%xa(2,k)
         var(31) = star%xa(14,k)
! var(32) X_Be7: YREC tracks Be9, not Be7 -> 0
         var(33) = star%xa(8,k)
         var(34) = star%xa(10,k)
         var(35) = star%xa(11,k)
! var(36) X_Ne20 and var(37)-var(40): not tracked -> 0

         write(u,'(1P5E16.9,x)') (var(i), i = 1, ivar)
      end do

      close(u)
      return
end subroutine write_fgong_pulse
