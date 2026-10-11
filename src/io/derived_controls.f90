!----------------------------------------------------------------------
! derived_controls
!----------------------------------------------------------------------
! The controls derived from the ones read (formerly part of parmin):
! tenv, the rotation flags, the wind law's coefficients and exponents,
! the composition mixture, and the other settings computed from the
! inlist rather than set by it.
module derived_controls
      use star_info_lib, only: star, control_nml_override, physics_nml_override
      use physical_constants
      use io_units
      use legacy_controls
      implicit none
      private
      public :: echo_settings

contains

! ---------------------------------------------------------------
! Write the full settings echo to the short/log stream (legacy
! layout, format statements local to this block).
subroutine echo_settings
      use math_lib
! 2026 log redesign: the PT TOL / O.S. / wind-index echo is deleted
! (settings provenance = the inlist_used copy in the output dir).
      star%tenv = 0.5d0*(tenv0 + tenv1)
      if(lrot) then
         lnew0 = .true.
         if(walpcz.lt.-2.0d0) walpcz = -2.0d0
         if(walpcz.gt.0.0d0) walpcz = 0.0d0
! JNT 09/2025 LSOLID OVERWRITES IMPJMOD
         if(lsolid) impjmod = 1
!CCCCC OLD OR NEW WINDLAW.
         if(.not.lmwind)then
!CCCCC INSTRUCTIONS FOR THE OLD WINDLAW
!CCCCC
!CCCCC SET UP COEFFICIENTS FOR ANGULAR MOMENTUM LOSS VIA WINDS.
!CCCCC GIVEN THE INDEX ALFA, THE FORMULA FOR JDOT IS
!CCCCC JDOT = FK*2.036D33*1.452D9**ALFA*(MDOT/1.0D-14)**(1-2ALFA/3)
!CCCCC *OMEGA**(1+4ALFA/3)*(R/RSUN)**(2-ALFA)*(M/MSUN)**-ALFA/3
!CCCCC EXMD = EXPONENT OF MDOT TERM; EXW SAME FOR OMEGA;EXR FOR R;EXM FOR M.
!CCCCC FK IS A FUDGE FACTOR,SET TO 1 TO REPRODUCE THE OBSERVED SOLAR ANGULAR
!CCCCC MOMENTUM LOSS.
            one_third = 1.0d0/3.0d0
            two_thirds = 2.0d0/3.0d0
            constfactor = fk*2.036d33*pow(1.452d9, alfa)
            exmd = 1.0d0 - two_thirds*alfa
            wind_law_omega_exponent = 1.0d0 + 2.0d0*two_thirds*alfa
            exr = 2.0d0 - alfa
            exm = -one_third*alfa
!
! INSTRUCTIONS FOR THE NEW WINDLAW
!
! SET UP INDICES FOR LOSS LAW IN TERMS OF PMM A, B, C, M.
! INCLUDE A PMMD AS WELL, FOR TURNING TAUCZ and CENTRIFUG DEPENDENCE ON/OFF.
! EXCEN = EXPONENT FOR CENTRIFUGAL TERM, {K2^2/(K2^2+0.5*W^2 R^3/GM)}^PMMM
!
! AWIND = 'K97' ENFORCES KRISHNAMURTHI (1997) WIND LAW.
! AWIND = 'V13' ENFORCES VAN SADERS + PINSONNEAULT (2013) WIND LAW, AS ADAPTED
!         FROM MATT ET AL. (2008,2012).
! AWIND = 'CUS' ADOPTS THE M-D VALUES GIVEN IN THE NAMELIST.
!
!
! THEN CONVERT THE PMM VALUES INTO EXPONENTS FOR THE TORQUE CALCULATION.
! BASIC EQUATIONS, WHERE EVERY TERM IS SCALED TO SOLAR, ARE:
!
! Jdot = w * Bmag^4m * Mdot^1-2m * R^5m+2 * M^-m * Fcen^md
! Bmag = Pphot^0.5 * (w * R^c)^b * Tcz^d
! GM/R * Mdot = (Lx/Lbol * Lbol)^a -> Mdot = R * M^-1 * Lbol * w^a
! ...where Lx/Lbol is assumed to scale with w to the power of a.
!
! THIS GIVES:
! JDOT ~ w ^ 1+a-2ma+4mb
!      ~ Tcz ^ 4md
!      ~ R ^ 3+3m+4mbc
!      ~ M ^ m-1
!      ~ Lbol ^ 1-2m
!      ~ Pphot ^ 2m
!      ~ Fcen ^ md
!
! For moment, ignore GM/R term in Mdot. This makes:
! JDOT ~ R ^ 2+5m+4mbc
! JDOT ~ M ^ -m
!
         else
            if(awind.eq.'K97')then
               pmmm = alfa/3.0
! MHP 8/17 CORRECT DEFAULT FOR A = 2, NOT 0
               pmma = 2.0
               pmmb = 1.0
               pmmc = 2.0
               pmmd = 0.0
               lbscale = .false.
            elseif(awind.eq.'V13')then
               pmmm = 0.22
               pmma = 2.0
               pmmb = 1.0
               pmmc = 0.0
               pmmd = 1.0
               lbscale = .true.
            endif
            wind_law_omega_exponent   = 1.0d0 + pmma - 2.0d0*pmma*pmmm + 4.0d0*pmmm*pmmb
! G Somers 8/17 ZERO'D OUT EXTAU. TAUCZ TERM NOW COMPUTED IN
! MWIND/MCOWIND, NOT IN AMCALC.
!            EXTAU = 4.0D0*PMMM
            extau = 0.0
! JvS 09/25 REMOVED TYPO IN EXR = 2.0D0+5.0D0*PMMM-4.0D0*PMMM*PMMB*PMMC
            exr   = 2.0d0+5.0d0*pmmm-4.0d0*pmmm*pmmc
            exm   = -pmmm
            exl   = 1.0d0 - 2.0d0*pmmm
! G Somers 11/17, ADDED PMMD TO SWTICH OFF IN K97 FORM.
            expr  = 2.0d0*pmmm*pmmd
            excen = pmmm*pmmd
! INITIALIZE CONSTANT FACTOR FOR CENTRIFUGAL TERM
            c_2 = 0.0506
! SET THE CONSTANT FACTOR
            constfactor = fk*pmmjd/pow(pmmsolw, wind_law_omega_exponent)
! IF RELEVANT RESET THE SATURATION THRESHOLD IN
! TERMS OF THE SOLAR ROTATION RATE.  WMAX_SUN<1000
! INDICATES SATURATION (AT THE SUN), SO
! WSAT = WMAX_SUN*PMMSOLW
            if(wmax_sun.lt.1.0e3)then
               wmax = wmax_sun*pmmsolw
            endif
         endif
      endif
!     WINDLAW END
! impjmod/pmma/pmmb/pmmc/pmmd/pmmm/lbscale/wmax can all be overridden
! above (the LSOLID-overwrites-IMPJMOD line and the K97/V13 PMM-windlaw
! branch), after they were already copy-assigned into their const_lib
! canonical names (solid_body_mode_flag/pmm_exponent_a/pmm_exponent_b/
! pmm_exponent_c/pmm_exponent_d/pmm_exponent_m/scale_by_b_field/
! wind_saturation_omega) earlier in this subroutine -- re-sync now so
! the const_lib values reflect any override from this block, not just
! the raw namelist read.
      solid_body_mode_flag = impjmod
      pmm_exponent_a = pmma
      pmm_exponent_b = pmmb
      pmm_exponent_c = pmmc
      pmm_exponent_d = pmmd
      pmm_exponent_m = pmmm
      scale_by_b_field = lbscale
      wind_saturation_omega = wmax
!
      if(lnewcp) then
       value_relative_to_h = .true.
       if(atmp.eq.'ABS') value_relative_to_h = .false.
! DECIDE WHICH ELEMENT IN ARRAY HCOMP TO BE RESCALED
! USING CHARACTER ARRAY AID AND INPUT CHARACTER VARIABLE ANEWCP
       do i = 1,12
          if(anewcp.eq.element_id(i)) then
! INEWCP IS THE INDEX OF THE ELEMENT BEING ALTERED
            new_species_index = i + 3
            exit
          endif
       end do
       if (i > 12) then
! ANEWCP NOT A RECOGNIZED ELEMENT
       lnewcp = .false.
       write(run_log_unit,20) anewcp
      20    format(1x,'VARIABLE',a4,1x,'NOT A RECOGNIZED ELEMENT'/1x, &
           &    'RESCALING NOT PERFORMED')
       end if
      endif
      change_cno_mixture_active = .false.
      change_isotope_ratios_active = .false.
      if(isetmix.eq.1)then
! IF DEFAULT MIX (GS98) CNOFRACS ARE ALREADY SET.
! FOR A CUSTOM MIX,DISABLE IF THE SUM OF CNO MASS FRACTIONS EXCEEDS ONE
! OR IF ANY MASS FRACTION IS NEGATIVE
! (Restructured 2026: the goto 602/606 aborts left
! change_cno_mixture_active false, which already suppresses the
! label-606 print, so structured fall-through is equivalent.)
         if(amix.eq.'CUS')then
            if(frac_c.lt.0.0d0.or.frac_n.lt.0.0d0.or.frac_o.lt.0.0d0)then
               write(*,591)frac_c,frac_n,frac_o
               write(run_log_unit,591)frac_c,frac_n,frac_o
      591          format('NEGATIVE INPUT CNO FRACTION ',3e12.4, &
           &              ' MIX NOT MODIFIED')
            else
            sum_frac=frac_c+frac_n+frac_o
            if(sum_frac.ge.1.0d0)then
               write(*,598)frac_c,frac_n,frac_o
               write(run_log_unit,598)frac_c,frac_n,frac_o
      598          format('INPUT CNO FRACTION ',3e12.4, &
           &              ' EXCEEDS 1. MIX NOT MODIFIED')
            else
! VALID MIXTURE, USE CUSTOM ENTRIES FROM .NML1
            change_cno_mixture_active = .true.
            endif
            endif
         else
! SEARCH THROUGH OTHER VALID MIXTURE ENTRIES;IF FOUND,ASSIGN
         do i = 1,4
            if(amix.eq.mixture_id_table(i))then
               zxmix = zx_mix_table(i)
               frac_c = frac_c_table(i)
               frac_n = frac_n_table(i)
               frac_o = frac_o_table(i)
               change_cno_mixture_active = .true.
               exit
            endif
         end do
         if (.not. change_cno_mixture_active) then
!     NO VALID MIX SPECIFIED
         write(*,589)amix
         write(run_log_unit,589)amix
      589    format(1x,'warning: CNO mixture ',a8,' not recognized; mixture not altered')
         end if
      endif
      endif
      if(change_cno_mixture_active)then
         write(*,604)amix,frac_c,frac_n,frac_o
         write(run_log_unit,604)amix,frac_c,frac_n,frac_o
      604    format(1x,'CNO mixture ',a8,': C =',es12.4,'  N =',es12.4, &
           &         '  O =',es12.4,' (fractions of Z); applied to the starting model')
      endif
!     CHECK IF ISOTOPE RATIOS NEED TO BE ALTERED
      if(isetiso.eq.1)then
! FOR A CUSTOM MIX,DISABLE IF THE SUM OF CNO MASS FRACTIONS EXCEEDS ONE
! OR IF ANY MASS FRACTION IS NEGATIVE
         if(aiso.eq.'CUS')then
            if(r12_13.lt.0.0d0 .or. r14_15.lt.0.0d0 .or. r16_17.lt.0.0d0 &
           &  .or. r16_18.lt.0.0d0 .or. xh2_ini.lt.0.0d0 .or. xhe3_ini.lt.0.0d0 &
           &  .or. xli6_ini.lt.0.0d0 .or. xli7_ini.lt.0.0d0 .or. &
           &  xbe9_ini.lt.0.0d0.or.xb10_ini.lt.0.0d0 .or.xb11_ini.lt.0.0d0)then
               write(*,596)r12_13,r14_15,r16_17,r16_18, &
           &  xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini,xb10_ini,xb11_ini
               write(run_log_unit,596)r12_13,r14_15,r16_17,r16_18, &
           &  xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini,xb10_ini,xb11_ini
      596          format(1x,'warning: negative isotope ratio or light-element', &
           &   ' mass fraction (',11es12.4,'); mixture not modified')
            else
            sum_frac= xh2_ini+xhe3_ini+xli6_ini+xli7_ini+xbe9_ini+ &
           &                 xb10_ini+xb11_ini
            if(sum_frac.ge.1.0d0)then
               write(*,595)xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini, &
           &  xb10_ini,xb11_ini
               write(run_log_unit,595)xh2_ini,xhe3_ini,xli6_ini,xli7_ini, &
           &  xbe9_ini,xb10_ini,xb11_ini
      595          format(1x,'warning: light-element mass fractions sum above 1', &
           &  ' (',11es12.4,'); mixture not modified')
            else
!     PASSED ALL CHECKS - THE CUSTOM SETTINGS WILL BE APPLIED
            change_isotope_ratios_active = .true.
            endif
            endif
!     CURRENTLY THERE ARE ONLY 2 VALID OPTIONS - THE DEFAULT (L21) OR
! A CUSTOM MIXTURE (CUS) - IF NEITHER IS TRUE, EXIT
         else if(aiso.eq.'L21')then
!     THE DEFAULT SETTINGS WILL BE APPLIED
            change_isotope_ratios_active = .true.
         endif
      endif
      if(change_isotope_ratios_active)then
         write(*,605)aiso,r12_13,r16_18, &
           &  xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini
         write(run_log_unit,605)aiso,r12_13,r16_18, &
           &  xh2_ini,xhe3_ini,xli6_ini,xli7_ini,xbe9_ini
      605    format(1x,'isotope/light-element mixture ',a8,': C12/C13 =', &
           &    es12.4,'  O16/O18 =',es12.4,'  H2 =',es12.4,'  He3 =',es12.4, &
           &    '  Li6 =',es12.4,'  Li7 =',es12.4,'  Be9 =',es12.4, &
           &    '; applied to the starting model')
      endif
! DBG 12/95 ENSURE CORRECT PARAMETERS FOR Z DIFFUSION
      if (ldifz) then
           ldify=.true.
         lthoul=.true.
      end if

! G SOMERS 04/15; ENSURE CORRECT PARAMETERS FOR LIGHT ELEMENT DIFFUSION.
      if(ldifli)then
          lnewdif=.true.
          ldify=.true.
          ldifz=.true.
          lthoul=.true.
          lthoulfit=.false.
          ilambda=4
      endif
! ldify/ldifz/lnewdif/lthoulfit/ilambda can all be overridden above,
! after they were already copy-assigned into their const_lib canonical
! names (diffuse_helium_active/use_diffusion_z/
! use_new_diffusion_routines/use_thoul_fit/coulomb_log_choice) earlier
! in this subroutine -- re-sync now so the const_lib values reflect
! any override from this block, not just the raw namelist read.
      diffuse_helium_active = ldify
      use_diffusion_z = ldifz
      use_new_diffusion_routines = lnewdif
      use_thoul_fit = lthoulfit
      coulomb_log_choice = ilambda
      use_thoul_diffusion = lthoul
! lnewcp/frac_c/frac_n/frac_o/zxmix can likewise be overridden above
! (the ANEWCP-rescaling and CNO-mixture validation blocks), after
! already being copy-assigned into their const_lib canonical names
! earlier in this subroutine -- re-sync for the same reason as above.
      rescale_species_active = lnewcp
      target_carbon_cno_fraction = frac_c
      target_nitrogen_cno_fraction = frac_n
      target_oxygen_cno_fraction = frac_o
      target_metal_fraction = zxmix

! 2026 log redesign: the STANDARD/CURRENT settings tables (LINE 1-11)
! are deleted -- the verbatim inlist copy in the output directory
! (inlist_used) is the settings provenance now.
      if(npenv.le.0) npenv = 9999
      if(pulse_gyre_interval.lt.0) pulse_gyre_interval = 0
      if(kttau .eq. 0) then
           write(run_log_unit, 197)
      else if (kttau .eq. 1) then
           write(run_log_unit, 198)
      else if (kttau .eq. 2) then
           write(run_log_unit, 1999)
      else if (kttau .eq. 3) then
           write(run_log_unit, 1888)
      else if (kttau .eq. 4) then
           write(run_log_unit, 1889)
! JNT 6/14 ADD FOR NEW KURUCZ/CASTELLI ATMOSPHERE TABLES
      else if (kttau .eq. 5) then
           write(run_log_unit, 1887)
      end if
      197 format(1x,'surface boundary: Eddington gray T(tau) relation')
      198 format(1x,'surface boundary: Krishna-Swamy T(tau) relation')
      1999 format(1x,'surface boundary: Harvard-Smithsonian reference atmosphere')
      1888 format(1x,'surface boundary: Kurucz atmosphere tables')
      1889 format(1x,'surface boundary: Allard atmosphere tables')
      1887 format(1x,'surface boundary: Kurucz/Castelli atmosphere tables')

! DBG PULSE
      if (lpurez) then
          write(run_log_unit,'(1x,a)') 'opacity: pure C and N tables enabled'
      end if

      write(run_log_unit,314)
      314    format(/,1x,100('='))
      write(run_log_unit,version_fmt) yrec_version_string, git_hash_string
      write(run_log_unit,'(1x,2a)') 'description: ', trim(descrip(1))
      if (len_trim(descrip(2)) > 0) &
           write(run_log_unit,'(14x,a)') trim(descrip(2))
      write(run_log_unit,'(1x,100(''=''))')


end subroutine echo_settings

end module derived_controls
