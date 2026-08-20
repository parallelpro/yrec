!----------------------------------------------------------------------
! lirate88
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original lirate88.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! This routine computes the burning rates for Li6, Li7, Be9 at the
! beginning and end of a timestep and stores them in common blocks
! newrat and oldrat.
!
! Burning rates from Caughlin and Fowler (1988).
subroutine lirate88(composition, log_density, log_temperature, num_zones, &
     use_current_model)
      implicit none
      integer, parameter :: json=5000

      double precision, intent(in) :: composition(15,json)
      double precision, intent(in) :: log_density(json)
      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: num_zones, use_current_model

! common/const1/: only ln10 and cc13 are used here. Naming matches
! dburn.f90.
      double precision :: ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi
      common/const1/ ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi

! common/newrat/: lithium/beryllium burning rates at the end of the
! timestep. Naming matches liburn.f90.
      double precision :: rate_li6(json), rate_li7(json), rate_be9(json)
      common/newrat/ rate_li6, rate_li7, rate_be9

! common/oldmod/: previous-timestep model snapshot; only
! old_composition, old_density, and old_temperature are used here.
! Naming matches dburn.f90.
      double precision :: old_pressure(json), old_temperature(json), &
           old_radius(json), old_luminosity(json), old_density(json), &
           old_composition(15,json), old_shell_mass(json), old_teff
      logical :: old_convective_flag(json), old_cz_flag(json)
      integer :: old_num_zones
      common/oldmod/ old_pressure, old_temperature, old_radius, &
           old_luminosity, old_density, old_composition, old_shell_mass, &
           old_convective_flag, old_cz_flag, old_teff, old_num_zones

! common/oldrat/: lithium/beryllium burning rates at the start of the
! timestep. Naming matches liburn.f90.
      double precision :: rate_li6_start(json), rate_li7_start(json), &
           rate_be9_start(json)
      common/oldrat/ rate_li6_start, rate_li7_start, rate_be9_start

! common/burnscs/: G Somers 6/14, light-element burning rate scale
! factors (used to scale the cross sections relative to the CF88
! rates below).
      double precision :: li6_rate_scale, li7_rate_scale, be9_pg_rate_scale, &
           be9_pd_rate_scale, be9_palpha_rate_scale
      common/burnscs/ li6_rate_scale, li7_rate_scale, be9_pg_rate_scale, &
           be9_pd_rate_scale, be9_palpha_rate_scale
! G Somers END

      double precision :: tlim
      data tlim/6.0d0/
      save

      integer :: zone_idx, tail_idx
      double precision :: rhox, t9, t913, t923, t943, t953, t932, t934
      double precision :: fli6, fli7, fbe91, fbe92, fbe93, fx, ex, fsbe9
      double precision :: t9a, c56

      do 100 zone_idx = 1,num_zones
         if(log_temperature(zone_idx).lt.tlim.and.old_temperature(zone_idx).lt.tlim)goto 110
         if(use_current_model.eq.1)then
            rhox = exp(ln10*log_density(zone_idx))*composition(1,zone_idx)
            t9=exp(ln10*(log_temperature(zone_idx)-9.0d0))
         else
            rhox = exp(ln10*old_density(zone_idx))*old_composition(1,zone_idx)
            t9=exp(ln10*(old_temperature(zone_idx)-9.0d0))
         endif
         t913=t9**cc13
         t923=t913*t913
         t943=t923*t923
! MHP 10/91 ADDED DEFINITION
         t953 = t943*t913
         t932=t9**1.5d0
         t934=t9**7.5d-1
!
! LI6(P,HE3)ALPHA
!
!           FLI6=3.73D10/T923*EXP(-8.413/T913-(T9/5.50)**2)*
!    1      (1.+%TYPO%0.50*T913-0.061*T923-0.021*T9+0.006*T943+0.006*T953)
!    2      +1.33D10/T932*EXP(-17.763/T9)+1.29D9/T9*EXP(-21.820/T9)
! ***TYPO IN OLD RATE***
         fli6=3.73d10/t923*exp(-8.413d0/t913-(t9/5.50d0)**2)* &
         (1.d0+0.050d0*t913-0.061d0*t923-0.021d0*t9+0.006d0*t943 &
          +0.005d0*t953)+1.33d10/t932*exp(-17.763d0/t9) &
          +1.29d9/t9*exp(-21.820d0/t9)
!
! LI7(P,ALFA)HE4
!
!           FLI7=8.04D08/T923*EXP(-8.471/T913-(T9/30.068)**2)*
!    1      (1.+0.049*T913+0.230*T923+0.079*T9-0.027*T943-0.023*T953)
!    2      +1.54D06/T932*EXP(-4.479/T9)+1.07D10/T932*EXP(-30.443/T9)
         t9a = t9/(1.0d0+0.759d0*t9)
         c56 = 1.25d0*cc23
         fli7=1.096d9/t923*exp(-8.472d0/t913)-4.830d8*t9a**c56/ &
         t932*exp(-8.472d0/t9a**cc13)+1.06d10/t932*exp(-30.442d0/t9)
!
! BE9(P,GAMMA)B10
!
         fbe91=1.33d+07/t923*exp(-10.359d0/t913-(t9/0.846d0)**2) &
         *(1.0d0+4.0d-2*t913+1.52d0*t923+4.28d-1*t9+2.15d0*t943+ &
         1.54d0*t953) + 9.64d+4/t932*exp(-3.445d0/t9) + &
         2.72d+6/t932*exp(-10.62d0/t9)
!    3      2.76D+6/T932*EXP(-10.62D0/T9) %OLD LAST LINE%
!
! BE9(P,D)2HE4 - UNCHANGED.
!
         fx=2.11d+11/t923*exp(-10.359d0/t913-(t9/5.2d-1)**2) &
         *(1.0d0+4.0d-2*t913+1.09d0*t923+3.07d-1*t9+3.21d0*t943 &
         +2.3d0*t953)
         ex = exp(-3.046d0/t9)/t9
         fbe92=fx+5.79d+8*ex+8.5d+8/t934*exp(-5.8d0/t9)
!
! BE9(P,ALPHA)LI6
!
         fbe93=fx+4.51d+8*ex+6.7d+8/t934*exp(-5.16d0/t9)
! G Somers 6/14, SCALE BY THE NEW CROSS SECTIONS
         fli6=li6_rate_scale*fli6
         fli7=li7_rate_scale*fli7
         fbe91=be9_pg_rate_scale*fbe91
         fbe92=be9_pd_rate_scale*fbe92
         fbe93=be9_palpha_rate_scale*fbe93
! G Somers END
! SUM RATES
         fsbe9 = fbe91 + fbe92 + fbe93
         if(use_current_model.eq.1)then
            rate_li6(zone_idx) = fli6*rhox
            rate_li7(zone_idx) = fli7*rhox
            rate_be9(zone_idx) = fsbe9*rhox
         else
            rate_li6_start(zone_idx) = fli6*rhox
            rate_li7_start(zone_idx) = fli7*rhox
            rate_be9_start(zone_idx) = fsbe9*rhox
         endif
  100 continue
  110 continue
      do 120 tail_idx = zone_idx,num_zones
         rate_li6(tail_idx)=0.0d0
         rate_li7(tail_idx)=0.0d0
         rate_be9(tail_idx)=0.0d0
         rate_li6_start(tail_idx)=0.0d0
         rate_li7_start(tail_idx)=0.0d0
         rate_be9_start(tail_idx)=0.0d0
  120 continue
      return
end subroutine lirate88
