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

      use star_info_lib, only: star
      use star_info_lib, only: star
      use luout_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: hcomp(15,json), hd(json), hl(json), &
           hp(json), hr(json), hs1(json)
      integer, intent(in) :: m, model








      save

! --- locals ---
      double precision :: np1
      double precision :: smtot, d, p, r, u, v, w
      integer :: i, ibeg, iend

!  FIND THE MILNE INVARIANTS U AND V, ALONG WITH THE INDEX N+1.
!  WRITE THEM OUT TO LOGICAL UNIT IMILNE.
!  HEADER
      smtot = dexp(ln10*star%env_comp%stotal)/solar_mass_cgs
      write(imilne,5) model,star%env_comp%envelope_hydrogen_fraction, &
           star%env_comp%envelope_metal_fraction,smtot
    5 format(10X,'MODEL',I5,'  XENV =',1PD10.3,'  ZENV =',D10.3, &
           '  MASS(SOLAR UNITS) =',D10.3)
!  CALCULATE FOR FIRST POINT(ALWAYS DONE)
      d = dexp(ln10*hd(1))
      p = dexp(ln10*hp(1))
      r = dexp(ln10*hr(1))
      u = c4pi*d*r**3/hs1(1)
      v = dexp(ln10*cgl)*hs1(1)*d/(p*r)
      w = u*hs1(1)*(star%diag%sesum(1)+star%diag%seg(7,1))/(hl(1)*solar_luminosity_cgs)
      np1 = 1.0d0/star%diag%del_grad(2,1)
      write(imilne,10)1,hs1(1),r,p,d,hcomp(1,1),hl(1),u,v,w,np1
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
          w = u*hs1(i)*(star%diag%sesum(i)+star%diag%seg(7,i))/(hl(i)*solar_luminosity_cgs)
          np1 = 1.0d0/star%diag%del_grad(2,i)
          write(imilne,10)i,hs1(i),r,p,d,hcomp(1,i),hl(i), &
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
       w = u*hs1(m)*(star%diag%sesum(m)+star%diag%seg(7,m))/(hl(m)*solar_luminosity_cgs)
       np1 = 1.0d0/star%diag%del_grad(2,m)
       write(imilne,10)m,hs1(m),r,p,d,hcomp(1,m),hl(m), &
                         u,v,w,np1
      endif
      close(imilne)
      return
end subroutine wrtmil
