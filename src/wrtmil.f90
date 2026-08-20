!----------------------------------------------------------------------
! wrtmil
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original wrtmil.f; only variable names, source form, and comment
! style were updated.
!
! Computes the Milne invariants U and V (and the polytropic-index-like
! quantity N+1) shell by shell and writes them to the .milne logical
! unit, printed every nprtpt points (plus the first and last points).
subroutine wrtmil(hcomp, hd, hl, hp, hr, hs1, m, model)

      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: hcomp(15,json), hd(json), hl(json), &
           hp(json), hr(json), hs1(json)
      integer, intent(in) :: m, model

! common/luout/: only milne_file_unit is used here. Naming matches
! getopac.f90 (short_file_unit there is the ishort slot; here we need
! the imilne slot).
      integer :: ilast, idebug, itrack, ishort, milne_file_unit, imodpt, &
           istor, iowr
      common/luout/ ilast, idebug, itrack, ishort, milne_file_unit, &
           imodpt, istor, iowr

! common/ccout1/: only print_point_interval (NPRTPT) is used here.
      integer :: npenv, nprtmod, print_point_interval, npoint
      common/ccout1/ npenv, nprtmod, print_point_interval, npoint

! common/comp/: only envelope_hydrogen_fraction/envelope_metal_fraction
! are used here. Naming matches getopac.f90.
      double precision :: envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, amuenv, fxenv(12), xnew, znew, stotal, senv
      common/comp/ envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, amuenv, fxenv, xnew, znew, stotal, senv

! common/const/: only solar_luminosity_cgs/solar_mass_cgs are used
! here. Naming matches vcirc.f90.
      double precision :: solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude
      common/const/ solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude

! common/const1/: only ln10/c4pi are used here. Naming matches
! eqburn.f90.
      double precision :: ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi
      common/const1/ ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi

! common/const2/: only cgl is used here. Naming matches meqos.f90.
      double precision :: gas_constant, radiation_constant_over_3, ca3l, &
           csig, csigl, cgl, cmkh, cmkhn
      common/const2/ gas_constant, radiation_constant_over_3, ca3l, csig, &
           csigl, cgl, cmkh, cmkhn

! common/scrtch/: only sesum/seg are used here. Naming matches
! microdiff_setup.f90.
      double precision :: sesum(json), seg(7,json), sbeta(json), seta(json)
      logical :: locons(json)
      double precision :: so(json), del_grad(3,json), sfxion(3,json), &
           svel(json), scp(json)
      common/scrtch/ sesum, seg, sbeta, seta, locons, so, del_grad, &
           sfxion, svel, scp

      save

! --- locals ---
      double precision :: np1
      double precision :: smtot, d, p, r, u, v, w
      integer :: i, ibeg, iend

!  FIND THE MILNE INVARIANTS U AND V, ALONG WITH THE INDEX N+1.
!  WRITE THEM OUT TO LOGICAL UNIT IMILNE.
!  HEADER
      smtot = dexp(ln10*stotal)/solar_mass_cgs
      write(milne_file_unit,5) model,envelope_hydrogen_fraction, &
           envelope_metal_fraction,smtot
    5 format(10X,'MODEL',I5,'  XENV =',1PD10.3,'  ZENV =',D10.3, &
           '  MASS(SOLAR UNITS) =',D10.3)
!  CALCULATE FOR FIRST POINT(ALWAYS DONE)
      d = dexp(ln10*hd(1))
      p = dexp(ln10*hp(1))
      r = dexp(ln10*hr(1))
      u = c4pi*d*r**3/hs1(1)
      v = dexp(ln10*cgl)*hs1(1)*d/(p*r)
      w = u*hs1(1)*(sesum(1)+seg(7,1))/(hl(1)*solar_luminosity_cgs)
      np1 = 1.0d0/del_grad(2,1)
      write(milne_file_unit,10)1,hs1(1),r,p,d,hcomp(1,1),hl(1),u,v,w,np1
   10 format(1X,I4,10(1PE11.3))
!  PRINT OUT EVERY NPRTPT POINTS;LAST POINT ALWAYS PRINTED.
      iend = 1
      if(print_point_interval.le.m) then
       ibeg = max(2,print_point_interval)
       iend = m - mod(m,print_point_interval)
       do 20 i = ibeg,iend,print_point_interval
          d = dexp(ln10*hd(i))
          p = dexp(ln10*hp(i))
          r = dexp(ln10*hr(i))
          u = c4pi*d*r**3/hs1(i)
          v = dexp(ln10*cgl)*hs1(i)*d/(p*r)
          w = u*hs1(i)*(sesum(i)+seg(7,i))/(hl(i)*solar_luminosity_cgs)
          np1 = 1.0d0/del_grad(2,i)
          write(milne_file_unit,10)i,hs1(i),r,p,d,hcomp(1,i),hl(i), &
                            u,v,w,np1
   20    continue
      endif
      if(iend.lt.m) then
!  PRINT OUT LAST POINT IF NPRTPT DOESNT DIVIDE EVENLY INTO M.
       d = dexp(ln10*hd(m))
       p = dexp(ln10*hp(m))
       r = dexp(ln10*hr(m))
       u = c4pi*d*r**3/hs1(m)
       v = dexp(ln10*cgl)*hs1(m)*d/(p*r)
       w = u*hs1(m)*(sesum(m)+seg(7,m))/(hl(m)*solar_luminosity_cgs)
       np1 = 1.0d0/del_grad(2,m)
       write(milne_file_unit,10)m,hs1(m),r,p,d,hcomp(1,m),hl(m), &
                         u,v,w,np1
      endif
      close(milne_file_unit)
      return
end subroutine wrtmil
