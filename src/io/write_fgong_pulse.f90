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
      use const_lib
      implicit none

      integer, intent(in) :: n
      double precision, intent(in) :: pts(35, n)
      double precision, intent(in) :: mstar_g, rstar_cm, lstar_cgs
      character(len=*), intent(in) :: pulse_path

      integer, parameter :: iconst = 15, ivar = 40, ivers = 300
      double precision :: glob(iconst), var(ivar)
      integer :: u, j, k, i
      double precision :: radius_cm, mass_g, grav

      open(newunit=u, file=pulse_path, status='REPLACE', form='FORMATTED')

      glob = 0.0d0
      glob(1) = mstar_g
      glob(2) = rstar_cm
      glob(3) = lstar_cgs
      glob(4) = pts(22, n)
      glob(5) = pts(21, n)
! glob(6) (mixing length alpha), glob(11)/glob(12) (central second
! derivatives) not tracked here -- zero, as tools tolerate.
      glob(13) = star%dage*1.0d9
      glob(14) = 10.0d0**star%log_Teff
      glob(15) = exp(ln10*cgl)

      write(u,'(a)') 'YREC FGONG output (ported layout: MESA pulse_fgong)'
      write(u,'(a,i8)') 'model_number ', star%model_number
      write(u,'(a,1pe16.9,a,0pf10.6)') 'star_age_yr ', star%dage*1.0d9, &
           '  log_Teff ', star%log_Teff
      write(u,'(a)') ''
      write(u,'(4I10)') n, iconst, ivar, ivers
      write(u,'(1P5E16.9,x)') (glob(i), i = 1, iconst)

      do j = 1, n
         k = n - j + 1
         radius_cm = pts(1,k)
         mass_g = pts(2,k)

         var = 0.0d0
         var(1) = radius_cm
         if (mass_g > 0.0d0) then
            var(2) = log(mass_g/mstar_g)
         else
            var(2) = -1.0d99
         end if
         var(3) = pts(5,k)
         var(4) = pts(4,k)
         var(5) = pts(6,k)
         var(6) = pts(21,k)
         var(7) = pts(3,k)
         var(8) = pts(12,k)
         var(9) = pts(15,k)
         var(10) = pts(9,k)
         var(11) = pts(10,k)
         var(12) = pts(11,k)
         var(13) = pts(19,k)
         var(14) = pts(20,k)
         if (radius_cm > 0.0d0) then
            grav = exp(ln10*cgl)*mass_g/(radius_cm*radius_cm)
            var(15) = pts(8,k)*radius_cm/grav
         end if
! var(16) r_X: not tracked -> 0
         var(17) = pts(22,k)
         var(18) = rstar_cm - radius_cm
         var(19) = pts(34,k)
! var(20) unused
         var(21) = pts(23,k)
         var(22) = pts(24,k)
         var(23) = pts(25,k)
         var(24) = pts(26,k)
         var(25) = pts(27,k)
! var(26)-var(28) dlnGamma1 partials: not tracked -> 0
         var(29) = pts(28,k)
         var(30) = pts(29,k)
         var(31) = pts(30,k)
! var(32) X_Be7: YREC tracks Be9, not Be7 -> 0
         var(33) = pts(31,k)
         var(34) = pts(32,k)
         var(35) = pts(33,k)
! var(36) X_Ne20 and var(37)-var(40): not tracked -> 0

         write(u,'(1P5E16.9,x)') (var(i), i = 1, ivar)
      end do

      close(u)
      return
end subroutine write_fgong_pulse
