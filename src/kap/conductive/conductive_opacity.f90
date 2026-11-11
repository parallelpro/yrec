!----------------------------------------------------------------------
! conductive_opacity
!----------------------------------------------------------------------
! Electron-conduction opacities from the table of Potekhin et al. (1999,
! A&A 346, 345; condall06.d, http://www.ioffe.rssi.ru/astro/conduct).
!     read_conductive_table        read the table, once per run (kap_init)
!     compute_conductive_opacity   the mixture's conductive opacity: the
!                                  table's conductivities for H, He and
!                                  an O proxy for the metals, weighted by
!                                  each species' ionized fraction
!     condopacp, cinterp3          Potekhin's interpolation routines,
!                                  under their own names
! The table is a conductive_table_set (opacity_tables) passed in. The
! formulas the F77 code kept as comments are in
! kap/retired/conductive_opacity.f90.txt; the sources are at the end.
module conductive_opacity
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use math_lib
      use numerics_lib, only: hunt
      use opacity_tables, only: cond_n_rho, cond_n_temp, cond_n_z, conductive_table_set, &
           kap_settings_t, kap_table_paths_t, opacity_table_state
      use opacity_source, only: kap_point, kap_result, kap_source
      use kap_messages, only: kap_terminal, kap_terminal_values, message_length, &
           unwritten_line
      use conductive_reader, only: read_conductive_table
      implicit none
      private
      public :: conductive_source
      public :: compute_conductive_opacity

! kappa_c = 16 sigma T^3/(3 rho K) for a thermal conductivity K: log10 of
! 16 sigma/3 (CGS), as the F77 code rounded it (log10(16 sigma/3) from
! physical_constants' sigma differs at 1e-5: a measured change)
      real(dp), parameter :: log10_16_sigma_over_3 = -3.5194d0
! the table's reach: log T in [conductive_logt_min, conductive_logt_max],
! log rho up to conductive_logrho_max; below conductive_logrho_min (the
! table's lowest density) the opacity is extrapolated where log R >= 0
      real(dp), parameter :: conductive_logt_min = 3.0d0, conductive_logt_max = 9.0d0, &
           conductive_logrho_min = -6.0d0, conductive_logrho_max = 9.75d0
! log10 kappa_c where the table gives none (no conduction, kappa ~ 1e10)
      real(dp), parameter :: no_conduction_log10_opacity = 9.9998d0
! the ion charges of the table's species: H, He and O, the metals' proxy
      real(dp), parameter :: charge_h1 = 1.0d0, charge_he4 = 2.0d0, charge_ox = 8.0d0

contains


!----------------------------------------------------------------------
! compute_conductive_opacity (formerly condopacpint)
!----------------------------------------------------------------------
! The mixture's conductive opacity and its log derivatives (kap_eval's
! conduction step, through the family's lookup), from the table's
! conductivities K_i (condopacp) for fully ionized H, He and O, the
! metals' proxy:
!   kappa_c = 16 sigma T^3/(3 rho K),  1/K = sum_i w_i/K_i,
! with weights w_i = n_i f_i: n_i the species' number fraction, f_i its
! ionized fraction (ion_fraction: H+, and (He+ + 2 He++)/2 for He; the
! metals fully ionized). got_conductive_opacity is .false. outside the
! table's reach (log T 3-9, log rho up to 9.75, and below log rho = -6
! only where log R >= 0), where kappa_c = 10**no_conduction_log10_opacity
! is returned unused.
subroutine compute_conductive_opacity(tbl, log10_density, log10_temperature, &
     hydrogen_fraction, metal_fraction, conductive_opacity, &
     conductive_log10_opacity, conductive_dlnkap_dlnrho, &
     conductive_dlnkap_dlnt, ion_fraction, got_conductive_opacity, ierr)


      type(conductive_table_set), intent(in) :: tbl
      real(dp), intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      real(dp), intent(out) :: conductive_opacity, &
           conductive_log10_opacity, conductive_dlnkap_dlnrho, &
           conductive_dlnkap_dlnt
      real(dp), intent(in) :: ion_fraction(:)
      logical, intent(out) :: got_conductive_opacity
      integer, intent(out) :: ierr
! atomic weights of H, He4 and the oxygen proxy for the metals
      real(dp), parameter :: atomic_weight_h1 = 1.008D0, &
           atomic_weight_he4 = 4.004D0, atomic_weight_ox = 16D0

! --- locals ---
      real(dp) :: helium_fraction
      real(dp) :: number_fraction_h1, number_fraction_he4, &
           number_fraction_ox
      real(dp) :: weight_h1, weight_he4, weight_ox
      real(dp) :: extrapolation_indicator
      real(dp) :: log10_cond_h1, dlnkap_dlnrho_h1, dlnkap_dlnt_h1
      real(dp) :: log10_cond_he4, dlnkap_dlnrho_he4, dlnkap_dlnt_he4
      real(dp) :: log10_cond_ox, dlnkap_dlnrho_ox, dlnkap_dlnt_ox
      real(dp) :: cond_h1, cond_he4, cond_ox
      real(dp) :: mix_log10_cond
      real(dp) :: unused_deriv1, unused_deriv2

! the number fractions of H1, He4 and the metals (as O)
      ierr = 0
      helium_fraction=1.0D0-hydrogen_fraction-metal_fraction
      number_fraction_h1  = (hydrogen_fraction/atomic_weight_h1)  / &
           (hydrogen_fraction/atomic_weight_h1 + &
           helium_fraction/atomic_weight_he4 + metal_fraction/atomic_weight_ox)
      number_fraction_he4 = (helium_fraction/atomic_weight_he4) / &
           (hydrogen_fraction/atomic_weight_h1 + &
           helium_fraction/atomic_weight_he4 + metal_fraction/atomic_weight_ox)
      number_fraction_ox  = (metal_fraction/atomic_weight_ox)  / &
           (hydrogen_fraction/atomic_weight_h1 + &
           helium_fraction/atomic_weight_he4 + metal_fraction/atomic_weight_ox)

! each weighted by its ionized fraction: ion_fraction(1) of H is ionized,
! ion_fraction(2) of He singly and ion_fraction(3) doubly; the metals (few
! by number) are taken as fully ionized
       weight_h1 = number_fraction_h1 * ion_fraction(1) ! (FXION(1)/1D0 is the fraction
       weight_he4 = number_fraction_he4 * (ion_fraction(2) + 2D0*ion_fraction(3))/2D0
       weight_ox = number_fraction_ox

! no conduction until the table gives one
       got_conductive_opacity=.FALSE.
       conductive_log10_opacity = no_conduction_log10_opacity
       conductive_opacity = pow(10D0, conductive_log10_opacity)
       conductive_dlnkap_dlnrho = conductive_opacity
       conductive_dlnkap_dlnt = conductive_opacity

       extrapolation_indicator = log10_density - 3.0D0*(log10_temperature-6.0D0)
! within the table's reach, each species fully ionized
      if(log10_temperature >= conductive_logt_min.and.log10_temperature <= conductive_logt_max &
           .and.log10_density <= conductive_logrho_max) then
         if (log10_density >= conductive_logrho_min) then
           call condopacp(tbl, charge_h1,log10_temperature,log10_density, &
                log10_cond_h1,dlnkap_dlnrho_h1,dlnkap_dlnt_h1,ierr)
           if (ierr /= 0) return
           call condopacp(tbl, charge_he4,log10_temperature,log10_density, &
                log10_cond_he4,dlnkap_dlnrho_he4,dlnkap_dlnt_he4,ierr)
           if (ierr /= 0) return
           call condopacp(tbl, charge_ox,log10_temperature,log10_density, &
                log10_cond_ox,dlnkap_dlnrho_ox,dlnkap_dlnt_ox,ierr)
           if (ierr /= 0) return

           cond_h1 = exp10((-log10_cond_h1))
           cond_he4 = exp10((-log10_cond_he4))
           cond_ox = exp10((-log10_cond_ox))

! mix = log10 K of the mixture, -log10(sum_i w_i 10^-K_i) for the
! species' log10 K_i; its derivative is the 1/K-weighted mean of theirs,
! d mix/dx = sum_i w_i 10^-K_i dK_i/dx / sum_i w_i 10^-K_i
           mix_log10_cond = -log10(weight_h1*cond_h1 + weight_he4*cond_he4 + &
                weight_ox*cond_ox)
           conductive_dlnkap_dlnrho = (weight_h1*cond_h1*dlnkap_dlnrho_h1 &
                + weight_he4*cond_he4*dlnkap_dlnrho_he4 + weight_ox*cond_ox*dlnkap_dlnrho_ox) &
                / (weight_h1*cond_h1 + weight_he4*cond_he4 + weight_ox*cond_ox)

           conductive_dlnkap_dlnt = (weight_h1*cond_h1*dlnkap_dlnt_h1 &
                + weight_he4*cond_he4*dlnkap_dlnt_he4 + weight_ox*cond_ox*dlnkap_dlnt_ox) &
                / (weight_h1*cond_h1 + weight_he4*cond_he4 + weight_ox*cond_ox)

         conductive_log10_opacity = log10_16_sigma_over_3+3.0D0*log10_temperature-log10_density-mix_log10_cond
           conductive_dlnkap_dlnrho = -1.0D0-conductive_dlnkap_dlnrho
           conductive_dlnkap_dlnt = 3.0D0-conductive_dlnkap_dlnt
           got_conductive_opacity=.TRUE.
           conductive_opacity = exp10(conductive_log10_opacity)
         end if
! below the table's lowest density, where log R >= 0: the main branch's
! formula at conductive_logrho_min, extrapolated linearly in log rho to
! the point with its slope there (the H table's slope, not the mixture's)
         if (log10_density < conductive_logrho_min.and.extrapolation_indicator >= 0.0D0) then
           call condopacp(tbl, charge_h1,log10_temperature,conductive_logrho_min,log10_cond_h1, &
                conductive_dlnkap_dlnrho,conductive_dlnkap_dlnt,ierr)
           if (ierr /= 0) return
           call condopacp(tbl, charge_he4,log10_temperature,conductive_logrho_min,log10_cond_he4, &
                unused_deriv1,unused_deriv2,ierr)
           if (ierr /= 0) return
           call condopacp(tbl, charge_ox,log10_temperature,conductive_logrho_min,log10_cond_ox, &
                unused_deriv1,unused_deriv2,ierr)
           if (ierr /= 0) return

           cond_h1 = exp10((-log10_cond_h1))
           cond_he4 = exp10((-log10_cond_he4))
           cond_ox = exp10((-log10_cond_ox))

           mix_log10_cond = -log10(weight_h1*cond_h1 + weight_he4*cond_he4 + &
                weight_ox*cond_ox)

           conductive_log10_opacity = log10_16_sigma_over_3+3.0D0*log10_temperature-conductive_logrho_min-mix_log10_cond
           conductive_dlnkap_dlnrho = -1.0D0-conductive_dlnkap_dlnrho
           conductive_dlnkap_dlnt = 3.0D0-conductive_dlnkap_dlnt
           conductive_log10_opacity = conductive_log10_opacity + &
                conductive_dlnkap_dlnrho*(log10_density-conductive_logrho_min)
           conductive_opacity = exp10(conductive_log10_opacity)
           got_conductive_opacity=.TRUE.
         end if
      end if

end subroutine compute_conductive_opacity

!----------------------------------------------------------------------
! condopacp
!----------------------------------------------------------------------
!  This subroutine interpolates the electron thermal conductivity
!               from the data file "condall.d" (YREC ships condall06.d)
!  ---------------------------------------------------  Version 23.05.99
! Input: ion_charge - ion charge, log10_temperature - lg(T[K]),
!        log10_density - lg(rho[g/cc])
! Output: log10_conductivity - Log_{10} thermal conductivity (kappa) [CGS units]
!         dlnkappa_dlnrho - d log kappa / d log rho
!         dlnkappa_dlnt - d log kappa / d log T
!!!      (it is also possible to obtain all second derivatives)      !!!
!
!     This Subroutine interpolates thermal conductivity from
!     a precalculated table condall.d available at
!     http://www.ioffe.rssi.ru/astro/conduct
!     For theoretical base of this calculation, see
!      A.Y.Potekhin, D.A.Baiko, P.Haensel, D.G.Yakovlev, 1999
!                  Astron. Astrophys. 346, 345.
!     Extension from strongly- to weakly-degenerate regime
!     has been done using the thermal averaging - for example
!         A.Y.Potekhin, 1999, Astron. Astrophys. 351, 787.
!     Please quote these publication when using this program
!     (Here the table is tbl, read once per run by read_conductive_table.)
!     Please address any questions/comments to Alexander Potekhin:
!     e-mail: palex@astro.ioffe.rssi.ru
subroutine condopacp(tbl, ion_charge, log10_temperature, log10_density, &
     log10_conductivity, dlnkappa_dlnrho, dlnkappa_dlnt, ierr)
      type(conductive_table_set), intent(in) :: tbl

      real(dp), intent(in) :: ion_charge, log10_temperature, &
           log10_density
      real(dp), intent(out) :: log10_conductivity, &
           dlnkappa_dlnrho, dlnkappa_dlnt
      integer, intent(out) :: ierr
! (the grid sizes cond_n_* are opacity_tables', and must be the table
! file's)
      integer :: z_index, t_index, r_index
      real(dp) :: log10_ion_charge
      integer :: t_index_lo, t_index_hi, r_index_lo, r_index_hi
      real(dp) :: cktmz0, drktmz0, dr2ktmz0, xr
      real(dp) :: ckt0z0, drkt0z0, dr2kt0z0
      real(dp) :: ckt1z0, drkt1z0, dr2kt1z0
      real(dp) :: cktpz0, drktpz0, dr2ktpz0
      real(dp) :: cktmz1, drktmz1, dr2ktmz1
      real(dp) :: ckt0z1, drkt0z1, dr2kt0z1
      real(dp) :: ckt1z1, drkt1z1, dr2kt1z1
      real(dp) :: cktpz1, drktpz1, dr2ktpz1
      real(dp) :: xz1, xz0
      real(dp) :: cktm, drktm, dr2ktm
      real(dp) :: ckt0, drkt0, dr2kt0
      real(dp) :: ckt1, drkt1, dr2kt1
      real(dp) :: cktp, drktp, dr2ktp
      real(dp) :: dt2k, drtk, drt2k, xt
      character(len=message_length) :: lines(1)

      if (.not. tbl%loaded) then
         lines = unwritten_line
         write(lines,*) 'condopacp: the conductive-opacity table is not loaded (kap_init)'
         call kap_terminal(lines)
         ierr = 1
         return
      end if
      ierr = 0
      log10_ion_charge=log10(ion_charge)
! hunt's guesses: none (each call bisects; the indices are locals, so
! there is no previous call's index to start from)
      z_index = 0
      t_index = 0
      r_index = 0
      call hunt(tbl%z_grid,cond_n_z,log10_ion_charge,z_index)
      if (z_index == 0.or.z_index == cond_n_z) then
         lines = unwritten_line
         write(lines,*) 'CONINTER: Z out of range'
         call kap_terminal(lines)
         ierr = 1
         return
      end if
      call hunt(tbl%temp_grid,cond_n_temp,log10_temperature,t_index)

      if (t_index == 0.or.t_index == cond_n_temp) then
                  call kap_terminal_values(tbl%z_grid)
                  lines = unwritten_line
                  write(lines,*) cond_n_z, log10_temperature, z_index
                  call kap_terminal(lines)
                  lines = unwritten_line
                  write(lines,*) 'CONINTER: T out of range'
                  call kap_terminal(lines)
                  ierr = 1
                  return
      end if

      call hunt(tbl%rho_grid,cond_n_rho,log10_density,r_index)
      if (r_index == 0.or.r_index == cond_n_rho) then
         lines = unwritten_line
         write(lines,*) 'CONINTER: rho out of range'
         call kap_terminal(lines)
         ierr = 1
         return
      end if
      t_index_lo=max(1,t_index-1)
      t_index_hi=min(cond_n_temp,t_index+2)
      r_index_lo=max(1,r_index-1)
      r_index_hi=min(cond_n_rho,r_index+2)
! Cubic interpolation in RLG:
! Z0:
      call cinterp3(tbl%rho_grid(r_index_lo),tbl%rho_grid(r_index),tbl%rho_grid(r_index+1), &
           tbl%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           tbl%log10_kappa_table(t_index_lo,r_index_lo,z_index), &
           tbl%log10_kappa_table(t_index_lo,r_index,z_index), &
           tbl%log10_kappa_table(t_index_lo,r_index+1,z_index), &
           tbl%log10_kappa_table(t_index_lo,r_index_hi,z_index), &
           cktmz0,drktmz0,dr2ktmz0,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(tbl%rho_grid(r_index_lo),tbl%rho_grid(r_index),tbl%rho_grid(r_index+1), &
           tbl%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           tbl%log10_kappa_table(t_index,r_index_lo,z_index), &
           tbl%log10_kappa_table(t_index,r_index,z_index), &
           tbl%log10_kappa_table(t_index,r_index+1,z_index), &
           tbl%log10_kappa_table(t_index,r_index_hi,z_index), &
           ckt0z0,drkt0z0,dr2kt0z0,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(tbl%rho_grid(r_index_lo),tbl%rho_grid(r_index),tbl%rho_grid(r_index+1), &
           tbl%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           tbl%log10_kappa_table(t_index+1,r_index_lo,z_index), &
           tbl%log10_kappa_table(t_index+1,r_index,z_index), &
           tbl%log10_kappa_table(t_index+1,r_index+1,z_index), &
           tbl%log10_kappa_table(t_index+1,r_index_hi,z_index), &
           ckt1z0,drkt1z0,dr2kt1z0,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(tbl%rho_grid(r_index_lo),tbl%rho_grid(r_index),tbl%rho_grid(r_index+1), &
           tbl%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           tbl%log10_kappa_table(t_index_hi,r_index_lo,z_index), &
           tbl%log10_kappa_table(t_index_hi,r_index,z_index), &
           tbl%log10_kappa_table(t_index_hi,r_index+1,z_index), &
           tbl%log10_kappa_table(t_index_hi,r_index_hi,z_index), &
           cktpz0,drktpz0,dr2ktpz0,xr,ierr)
      if (ierr /= 0) return
! Z1:
      call cinterp3(tbl%rho_grid(r_index_lo),tbl%rho_grid(r_index),tbl%rho_grid(r_index+1), &
           tbl%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           tbl%log10_kappa_table(t_index_lo,r_index_lo,z_index+1), &
           tbl%log10_kappa_table(t_index_lo,r_index,z_index+1), &
           tbl%log10_kappa_table(t_index_lo,r_index+1,z_index+1), &
           tbl%log10_kappa_table(t_index_lo,r_index_hi,z_index+1), &
           cktmz1,drktmz1,dr2ktmz1,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(tbl%rho_grid(r_index_lo),tbl%rho_grid(r_index),tbl%rho_grid(r_index+1), &
           tbl%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           tbl%log10_kappa_table(t_index,r_index_lo,z_index+1), &
           tbl%log10_kappa_table(t_index,r_index,z_index+1), &
           tbl%log10_kappa_table(t_index,r_index+1,z_index+1), &
           tbl%log10_kappa_table(t_index,r_index_hi,z_index+1), &
           ckt0z1,drkt0z1,dr2kt0z1,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(tbl%rho_grid(r_index_lo),tbl%rho_grid(r_index),tbl%rho_grid(r_index+1), &
           tbl%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           tbl%log10_kappa_table(t_index+1,r_index_lo,z_index+1), &
           tbl%log10_kappa_table(t_index+1,r_index,z_index+1), &
           tbl%log10_kappa_table(t_index+1,r_index+1,z_index+1), &
           tbl%log10_kappa_table(t_index+1,r_index_hi,z_index+1), &
           ckt1z1,drkt1z1,dr2kt1z1,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(tbl%rho_grid(r_index_lo),tbl%rho_grid(r_index),tbl%rho_grid(r_index+1), &
           tbl%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           tbl%log10_kappa_table(t_index_hi,r_index_lo,z_index+1), &
           tbl%log10_kappa_table(t_index_hi,r_index,z_index+1), &
           tbl%log10_kappa_table(t_index_hi,r_index+1,z_index+1), &
           tbl%log10_kappa_table(t_index_hi,r_index_hi,z_index+1), &
           cktpz1,drktpz1,dr2ktpz1,xr,ierr)
      if (ierr /= 0) return
! Linear interpolation in ZLG:
      xz1=(log10_ion_charge-tbl%z_grid(z_index))/(tbl%z_grid(z_index+1)-tbl%z_grid(z_index))
      xz0=1d0-xz1
      cktm=xz0*cktmz0+xz1*cktmz1
      drktm=xz0*drktmz0+xz1*drktmz1
      dr2ktm=xz0*dr2ktmz0+xz1*dr2ktmz1
      ckt0=xz0*ckt0z0+xz1*ckt0z1
      drkt0=xz0*drkt0z0+xz1*drkt0z1
      dr2kt0=xz0*dr2kt0z0+xz1*dr2kt0z1
      ckt1=xz0*ckt1z0+xz1*ckt1z1
      drkt1=xz0*drkt1z0+xz1*drkt1z1
      dr2kt1=xz0*dr2kt1z0+xz1*dr2kt1z1
      cktp=xz0*cktpz0+xz1*cktpz1
      drktp=xz0*drktpz0+xz1*drktpz1
      dr2ktp=xz0*dr2ktpz0+xz1*dr2ktpz1
! Cubic interpolation in TLG:
      call cinterp3(tbl%temp_grid(t_index_lo),tbl%temp_grid(t_index),tbl%temp_grid(t_index+1), &
           tbl%temp_grid(t_index_hi),log10_temperature,t_index,cond_n_temp, &
           cktm,ckt0,ckt1,cktp, & ! input: values of lg kappa
           log10_conductivity,dlnkappa_dlnt,dt2k,xt,ierr) ! lg kappa, d lg k / d lg T, d2 lg k / d2 lg T
      if (ierr /= 0) return
      call cinterp3(tbl%temp_grid(t_index_lo),tbl%temp_grid(t_index),tbl%temp_grid(t_index+1), &
           tbl%temp_grid(t_index_hi),log10_temperature,t_index,cond_n_temp, &
           drktm,drkt0,drkt1,drktp, & ! input: values of d lg k / d lg rho
           dlnkappa_dlnrho,drtk,drt2k,xt,ierr) ! d lg k / d lg rho, d2 lgk/(d lgT d lg rho)
      if (ierr /= 0) return
end subroutine condopacp

! Given 4 values of Z and 4 values of V, find VF corresponding to 5th Z
!                                                       Version 23.05.99
!   Output: interp_value - interpolated value of function
!           interp_deriv - interpolated derivative
!           interp_deriv2 - interpolated second derivative
!           interp_frac - fraction of the path from grid_index to
!                         grid_index+1
subroutine cinterp3(grid_lo2,grid_lo1,grid_hi1,grid_hi2,grid_target, &
     grid_index,grid_size,val_lo2,val_lo1,val_hi1,val_hi2, &
     interp_value,interp_deriv,interp_deriv2,interp_frac, ierr)


      real(dp), intent(in) :: grid_lo2, grid_lo1, grid_hi1, &
           grid_hi2, grid_target
      integer, intent(in) :: grid_index, grid_size
      real(dp), intent(in) :: val_lo2, val_lo1, val_hi1, val_hi2
      real(dp), intent(out) :: interp_value, interp_deriv, &
           interp_deriv2, interp_frac
      integer, intent(out) :: ierr

      real(dp) :: x, h, hm, v01, hp, v11, c2, c3
      character(len=message_length) :: lines(1)

      ierr = 0
      if (grid_index <= 0.or.grid_index >= grid_size) then
         lines = unwritten_line
         write(lines,*) 'CINTERP: N0 out of range'
         call kap_terminal(lines)
         ierr = 1
         return
      end if
      x=grid_target-grid_lo1
      h=grid_hi1-grid_lo1   ! basic interval
      interp_frac=x/h
      if (grid_index > 1) then
         hm=grid_lo1-grid_lo2  ! left adjoint interval
         v01=((val_hi1-val_lo1)/h**2+(val_lo1-val_lo2)/hm**2)/ &
         (1d0/h+1d0/hm) ! left derivative
      end if
      if (grid_index < grid_size-1) then
         hp=grid_hi2-grid_hi1 ! right adjoint interval
         v11=((val_hi1-val_lo1)/h**2+(val_hi2-val_hi1)/hp**2)/ &
         (1d0/h+1d0/hp) ! right derivative
      end if
      if (grid_index > 1.and.grid_index < grid_size-1) then   ! Cubic interpolation
         c2=3d0*(val_hi1-val_lo1)-h*(v11+2d0*v01)
         c3=h*(v01+v11)-2d0*(val_hi1-val_lo1)
         interp_value=val_lo1+v01*x+c2*interp_frac**2+c3*interp_frac**3
         interp_deriv=v01+(2d0*c2*interp_frac+3d0*c3*interp_frac**2)/h
         interp_deriv2=(2d0*c2+6d0*c3*interp_frac)/h**2
         return
      end if
      if (grid_index == 1) then   ! Quadratic interpolation
         c2=val_lo1-val_hi1+v11*h
         interp_value=val_hi1-v11*(h-x)+c2*(1d0-interp_frac)**2
         interp_deriv=v11-2d0*c2*(1d0-interp_frac)/h
         interp_deriv2=2d0*c2/h**2
      else  ! N0=MXNV-1
         c2=val_hi1-val_lo1-v01*h
         interp_value=val_lo1+v01*x+c2*interp_frac**2
         interp_deriv=v01+2d0*c2*interp_frac/h
         interp_deriv2=2d0*c2/h**2
      end if
end subroutine cinterp3


!----------------------------------------------------------------------
! The family's descriptor (opacity_source) and its procedures
!----------------------------------------------------------------------
! The Potekhin conductive opacity: in_table is .false. where the table
! gives none (kap_eval then tries Hubbard & Lampe's fit).
function conductive_source(settings) result(source)
      type(kap_settings_t), intent(in) :: settings
      type(kap_source) :: source

      source%name = 'Potekhin conductive'
      source%selected = settings%use_conductive_opacity
      source%load => load_conductive
      source%lookup => lookup_conductive
end function conductive_source

subroutine load_conductive(settings, paths, envelope_hydrogen_fraction, tables, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(kap_table_paths_t), intent(in) :: paths
      real(dp), intent(in) :: envelope_hydrogen_fraction
      type(opacity_table_state), intent(inout) :: tables
      integer, intent(out) :: ierr

      call read_conductive_table(tables%conductive, paths%conductive, ierr)
end subroutine load_conductive

subroutine lookup_conductive(settings, tables, point, res, in_table, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(opacity_table_state), intent(inout) :: tables
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      logical, intent(out) :: in_table
      integer, intent(out) :: ierr

      call compute_conductive_opacity(tables%conductive, point%log10_density, &
           point%log10_temperature, point%hydrogen_fraction, point%metal_fraction, &
           res%opacity, res%log10_opacity, res%dlnkap_dlnrho, res%dlnkap_dlnt, &
           point%ion_fraction, in_table, ierr)
end subroutine lookup_conductive

!----------------------------------------------------------------------
! Sources
!----------------------------------------------------------------------
! The table and condopacp, cinterp3: A. Y. Potekhin's (condall06.d and
! its interpolation code, version 23.05.99; Potekhin et al. 1999, A&A
! 346, 345; Potekhin 1999, A&A 351, 787), kept under their names.
! compute_conductive_opacity: LLP's CONDOPACPINT (5/31/04),
! after Chris Burke's conductive-opacity code. Since its conversion to
! Fortran 90 (2026): the derivatives are the mixture's (the 2026-08 bug
! sweep: they multiplied by mix and took the rho derivatives for T), the
! table is read once per run at kap_init and passed in as a table set,
! and the range stops return ierr.

end module conductive_opacity
