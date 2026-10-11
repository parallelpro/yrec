!----------------------------------------------------------------------
! conductive_opacity
!----------------------------------------------------------------------
! Electron-conduction opacities from the table of Potekhin et al. (1999,
! A&A 346, 345; condall06.d, http://www.ioffe.rssi.ru/astro/conduct).
!     read_conductive_table        read the table, once per run (kap_init)
!     compute_conductive_opacity   (formerly condopacpint, LLP 5/04) the
!                                  mixture's conductive opacity: the
!                                  table's values for H, He and an O
!                                  proxy for the metals, weighted by
!                                  each species' ionization fraction
!     condopacp, cinterp3, hunt    Potekhin's interpolation routines,
!                                  under their own names
module conductive_opacity
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use math_lib
      implicit none
      private
      public :: read_conductive_table, compute_conductive_opacity
      integer, parameter :: cond_n_temp = 19, cond_n_rho = 64, &
           cond_n_z = 15
! table_loaded_flag holds this value once condopacp has read the table.
      integer, parameter :: cond_table_loaded_marker = 12345

      type :: conductive_table_state
            real(dp) :: temp_grid(cond_n_temp)
            real(dp) :: rho_grid(cond_n_rho)
            real(dp) :: z_grid(cond_n_z)
            real(dp) :: log10_kappa_table(cond_n_temp, cond_n_rho, &
                 cond_n_z)
            integer :: table_loaded_flag = -1
      end type conductive_table_state

      type(conductive_table_state), save :: cond_table

contains

!----------------------------------------------------------------------
! read_conductive_table
!----------------------------------------------------------------------
! Read the Potekhin table (condall06.d) at path, once per run from
! kap_init. (2026-10: read_controls opened it on a
! unit and condopacp read it on its first lookup.)
subroutine read_conductive_table(path, ierr)
      character(len=*), intent(in) :: path
      integer, intent(out) :: ierr
      integer :: file_unit, z_index, t_index, r_index
      real(dp) :: z_grid_value

      open(newunit=file_unit, file=path, status='OLD', iostat=ierr)
      if (ierr /= 0) then
         write(*,*) 'read_conductive_table: cannot open ', trim(path)
         ierr = 1
         return
      end if
      read(file_unit,'(A)') ! skip the first line
      do z_index=1,cond_n_z
         read(file_unit,*) z_grid_value,(cond_table%temp_grid(t_index),t_index=1,cond_n_temp)
         cond_table%z_grid(z_index)=log10(z_grid_value)
         do r_index=1,cond_n_rho
            read(file_unit,*) cond_table%rho_grid(r_index), &
                 (cond_table%log10_kappa_table(t_index,r_index,z_index),t_index=1,cond_n_temp)
         enddo
      enddo
      close(file_unit)
      cond_table%table_loaded_flag=cond_table_loaded_marker
end subroutine read_conductive_table

!----------------------------------------------------------------------
! compute_conductive_opacity (formerly condopacpint)
!----------------------------------------------------------------------
!     compute_conductive_opacity -  Calculate  Conductive Opacities based on Chirs
!                     Burke's conductive opacity code and callin
!                     on Potehkin code and table.
!
!                      LLP    5/31/04
!
! Combines the fully-ionized-species Potekhin conductive opacity
! (condopacp.f90, called for pure H, He, and a Z=8 metal proxy) into a
! single mixture value, weighted by each species' actual ionization
! fraction (ion_fraction, from the compute_fully_ionized_eos solve
! upstream). Called from kap_lib.f90's kap_get whenever conductive
! opacity is enabled.
subroutine compute_conductive_opacity(log10_density, log10_temperature, &
     hydrogen_fraction, metal_fraction, conductive_opacity, &
     conductive_log10_opacity, conductive_dlnkap_dlnrho, &
     conductive_dlnkap_dlnt, ion_fraction, got_conductive_opacity, ierr)


      real(dp), intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      real(dp), intent(out) :: conductive_opacity, &
           conductive_log10_opacity, conductive_dlnkap_dlnrho, &
           conductive_dlnkap_dlnt
      real(dp), intent(in) :: ion_fraction(3)
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

! Get fractional abundances in nubers of atoms of H1, He4 and metals (Ox)
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

! In addition we now weight the H1, He and Ox by oncly
! considering those atoms H1 and He4 atomos which are
! charged. FXION(1) is the fraction of H1 atoms that are
! (singly) charged. FXION(2) is the fraction of He4 atoms
! that are singly charged and FXION(3) the fractoin that
! are doubly charged.  It is presumed that the metals (Ox)
! are fully ionized in first approximation. Because of
! their low abundance (tyically less than 2% divided by
! the atomic mass
       weight_h1 = number_fraction_h1 * ion_fraction(1) ! (FXION(1)/1D0 is the fraction
       weight_he4 = number_fraction_he4 * (ion_fraction(2) + 2D0*ion_fraction(3))/2D0
       weight_ox = number_fraction_ox

! We initially presume no significant conductive opacity
       got_conductive_opacity=.FALSE.
!       OC = 1D50
!       OCL = LOG10(OC)
       conductive_log10_opacity = 9.9998D0
       conductive_opacity = pow(10D0, conductive_log10_opacity)
       conductive_dlnkap_dlnrho = conductive_opacity
       conductive_dlnkap_dlnt = conductive_opacity

       extrapolation_indicator = log10_density - 3.0D0*(log10_temperature-6.0D0)
!     DO CONDUCTIVE OPACITY CORRECTION in fully ionized approximation
      if(log10_temperature >= 3.0D0.and.log10_temperature <= 9.0D0 &
           .and.log10_density <= 9.75D0) then
         if (log10_density >= -6.0D0) then
           call condopacp(1.0D0,log10_temperature,log10_density, &
                log10_cond_h1,dlnkap_dlnrho_h1,dlnkap_dlnt_h1,ierr)
           if (ierr /= 0) return
           call condopacp(2.0D0,log10_temperature,log10_density, &
                log10_cond_he4,dlnkap_dlnrho_he4,dlnkap_dlnt_he4,ierr)
           if (ierr /= 0) return
           call condopacp(8.0D0,log10_temperature,log10_density, &
                log10_cond_ox,dlnkap_dlnrho_ox,dlnkap_dlnt_ox,ierr)
           if (ierr /= 0) return

           cond_h1 = exp10((-log10_cond_h1))
           cond_he4 = exp10((-log10_cond_he4))
           cond_ox = exp10((-log10_cond_ox))

           mix_log10_cond = -log10(weight_h1*cond_h1 + weight_he4*cond_he4 + &
                weight_ox*cond_ox)
!           CONDL = LOG10( (X + 1D0) / 2D0 ) -
!     *             LOG10( X*condx + .5d0*Y*condy + .5d0*z*condz)

! 2026 (bugsweep Batch 2): with mix = -log10(sum_i w_i 10^-K_i) the
! derivative is the conductivity-weighted mean of the per-species
! dK_i/dx: d mix/dx = sum_i w_i cond_i dK_i/dx / sum_i w_i cond_i.
! The inherited code (a) multiplied by mix_log10_cond (~12-17)
! instead of dividing by the summed conductivity and (b) fed the
! rho-derivatives into the T-derivative as well (dlnkap_dlnt_* were
! computed but never read). Only the derivatives change; the
! conductive opacity itself is untouched.
           conductive_dlnkap_dlnrho = (weight_h1*cond_h1*dlnkap_dlnrho_h1 &
                + weight_he4*cond_he4*dlnkap_dlnrho_he4 + weight_ox*cond_ox*dlnkap_dlnrho_ox) &
                / (weight_h1*cond_h1 + weight_he4*cond_he4 + weight_ox*cond_ox)
!           QODC = 2D0 * CONDL / ( 1D0 + X ) *
!     *          (X*condx*qodcx +.5d0*Y*condy*qodcy +.5d0*Z*cobdz*qodcz)

           conductive_dlnkap_dlnt = (weight_h1*cond_h1*dlnkap_dlnt_h1 &
                + weight_he4*cond_he4*dlnkap_dlnt_he4 + weight_ox*cond_ox*dlnkap_dlnt_ox) &
                / (weight_h1*cond_h1 + weight_he4*cond_he4 + weight_ox*cond_ox)
!           QOTC = 2D0 * CONDL / ( 1D0 + X ) *
!     *          (X*condx*qotcx +.5d0*Y*condy*qotcy +.5d0*Z*condz*qotcz)

         conductive_log10_opacity = -3.5194D0+3.0D0*log10_temperature-log10_density-mix_log10_cond
           conductive_dlnkap_dlnrho = -1.0D0-conductive_dlnkap_dlnrho
           conductive_dlnkap_dlnt = 3.0D0-conductive_dlnkap_dlnt
           got_conductive_opacity=.TRUE.
           conductive_opacity = exp10(conductive_log10_opacity)
         end if
         if (log10_density < -6.0D0.and.extrapolation_indicator >= 0.0D0) then
! Extrapolate Conductive opacity
           call condopacp(1.0D0,log10_temperature,-6.0D0,log10_cond_h1, &
                conductive_dlnkap_dlnrho,conductive_dlnkap_dlnt,ierr)
           if (ierr /= 0) return
           call condopacp(2.0D0,log10_temperature,-6.0D0,log10_cond_he4, &
                unused_deriv1,unused_deriv2,ierr)
           if (ierr /= 0) return
           call condopacp(8.0D0,log10_temperature,-6.0D0,log10_cond_ox, &
                unused_deriv1,unused_deriv2,ierr)
           if (ierr /= 0) return

           cond_h1 = exp10((-log10_cond_h1))
           cond_he4 = exp10((-log10_cond_he4))
           cond_ox = exp10((-log10_cond_ox))

           mix_log10_cond = log10(weight_h1*cond_h1 + weight_he4*cond_he4 + &
                weight_ox*cond_ox)
!           CONDL = LOG10( (X + 1D0) / 2D0 ) -
!     *            LOG10( X*condx + .5d0*Y*condy + .5d0*Z*condz)

           conductive_log10_opacity = -3.5194D0+3.0D0*log10_temperature+6.0D0-mix_log10_cond
           conductive_dlnkap_dlnrho = -1.0D0-conductive_dlnkap_dlnrho
           conductive_dlnkap_dlnt = 3.0D0-conductive_dlnkap_dlnt
           conductive_log10_opacity = conductive_log10_opacity + &
                conductive_dlnkap_dlnrho*(-6.0D0-log10_density)
           conductive_opacity = exp10(conductive_log10_opacity)
           got_conductive_opacity=.TRUE.
         end if
      end if

end subroutine compute_conductive_opacity

!----------------------------------------------------------------------
! condopacp
!----------------------------------------------------------------------
!  This subroutine interpolates the electron thermal conductivity
!               from the data file "condall.d"
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
!     (The table is read once per run by kap_init:
!      conductive_opacity_table's read_conductive_table.)
!     Please address any questions/comments to Alexander Potekhin:
!     e-mail: palex@astro.ioffe.rssi.ru
subroutine condopacp(ion_charge, log10_temperature, log10_density, &
     log10_conductivity, dlnkappa_dlnrho, dlnkappa_dlnt, ierr)

      real(dp), intent(in) :: ion_charge, log10_temperature, &
           log10_density
      real(dp), intent(out) :: log10_conductivity, &
           dlnkappa_dlnrho, dlnkappa_dlnt
      integer, intent(out) :: ierr
! (the grid sizes cond_n_* and the table itself, cond_table, live in
! kap/conductive/conductive_opacity.f90 and must be consistent with "condall.d")
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

      if (cond_table%table_loaded_flag /= cond_table_loaded_marker) then
         write(*,*) 'condopacp: the conductive-opacity table is not loaded (kap_init)'
         ierr = 1
         return
      end if
      ierr = 0
      log10_ion_charge=log10(ion_charge)
      call hunt(cond_table%z_grid,cond_n_z,log10_ion_charge,z_index)
      if (z_index == 0.or.z_index == cond_n_z) then
! 2026 audit (section 7): stop -> ierr (runtime table-range guard;
! see kap_lib's kap_get).
         write(*,*) 'CONINTER: Z out of range'
         ierr = 1
         return
      end if
      call hunt(cond_table%temp_grid,cond_n_temp,log10_temperature,t_index)

!      if (IT.eq.0.or.IT.eq.MAXT) stop 'CONINTER: T out of range'
      if (t_index == 0.or.t_index == cond_n_temp) then
                  print*, cond_table%z_grid
                  print*, cond_n_z, log10_temperature, z_index
! 2026 (docs/ROADMAP.md stage 3): stop 'CONINTER: T out of range'
! converted to ierr (see kap_lib's kap_get); message preserved.
                  print*, 'CONINTER: T out of range'
                  ierr = 1
                  return
      endif

      call hunt(cond_table%rho_grid,cond_n_rho,log10_density,r_index)
      if (r_index == 0.or.r_index == cond_n_rho) then
! 2026 audit (section 7): stop -> ierr (runtime table-range guard;
! see kap_lib's kap_get).
         write(*,*) 'CONINTER: rho out of range'
         ierr = 1
         return
      end if
      t_index_lo=max0(1,t_index-1)
      t_index_hi=min0(cond_n_temp,t_index+2)
      r_index_lo=max0(1,r_index-1)
      r_index_hi=min0(cond_n_rho,r_index+2)
! Cubic interpolation in RLG:
! Z0:
      call cinterp3(cond_table%rho_grid(r_index_lo),cond_table%rho_grid(r_index),cond_table%rho_grid(r_index+1), &
           cond_table%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           cond_table%log10_kappa_table(t_index_lo,r_index_lo,z_index), &
           cond_table%log10_kappa_table(t_index_lo,r_index,z_index), &
           cond_table%log10_kappa_table(t_index_lo,r_index+1,z_index), &
           cond_table%log10_kappa_table(t_index_lo,r_index_hi,z_index), &
           cktmz0,drktmz0,dr2ktmz0,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(cond_table%rho_grid(r_index_lo),cond_table%rho_grid(r_index),cond_table%rho_grid(r_index+1), &
           cond_table%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           cond_table%log10_kappa_table(t_index,r_index_lo,z_index), &
           cond_table%log10_kappa_table(t_index,r_index,z_index), &
           cond_table%log10_kappa_table(t_index,r_index+1,z_index), &
           cond_table%log10_kappa_table(t_index,r_index_hi,z_index), &
           ckt0z0,drkt0z0,dr2kt0z0,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(cond_table%rho_grid(r_index_lo),cond_table%rho_grid(r_index),cond_table%rho_grid(r_index+1), &
           cond_table%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           cond_table%log10_kappa_table(t_index+1,r_index_lo,z_index), &
           cond_table%log10_kappa_table(t_index+1,r_index,z_index), &
           cond_table%log10_kappa_table(t_index+1,r_index+1,z_index), &
           cond_table%log10_kappa_table(t_index+1,r_index_hi,z_index), &
           ckt1z0,drkt1z0,dr2kt1z0,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(cond_table%rho_grid(r_index_lo),cond_table%rho_grid(r_index),cond_table%rho_grid(r_index+1), &
           cond_table%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           cond_table%log10_kappa_table(t_index_hi,r_index_lo,z_index), &
           cond_table%log10_kappa_table(t_index_hi,r_index,z_index), &
           cond_table%log10_kappa_table(t_index_hi,r_index+1,z_index), &
           cond_table%log10_kappa_table(t_index_hi,r_index_hi,z_index), &
           cktpz0,drktpz0,dr2ktpz0,xr,ierr)
      if (ierr /= 0) return
! Z1:
      call cinterp3(cond_table%rho_grid(r_index_lo),cond_table%rho_grid(r_index),cond_table%rho_grid(r_index+1), &
           cond_table%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           cond_table%log10_kappa_table(t_index_lo,r_index_lo,z_index+1), &
           cond_table%log10_kappa_table(t_index_lo,r_index,z_index+1), &
           cond_table%log10_kappa_table(t_index_lo,r_index+1,z_index+1), &
           cond_table%log10_kappa_table(t_index_lo,r_index_hi,z_index+1), &
           cktmz1,drktmz1,dr2ktmz1,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(cond_table%rho_grid(r_index_lo),cond_table%rho_grid(r_index),cond_table%rho_grid(r_index+1), &
           cond_table%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           cond_table%log10_kappa_table(t_index,r_index_lo,z_index+1), &
           cond_table%log10_kappa_table(t_index,r_index,z_index+1), &
           cond_table%log10_kappa_table(t_index,r_index+1,z_index+1), &
           cond_table%log10_kappa_table(t_index,r_index_hi,z_index+1), &
           ckt0z1,drkt0z1,dr2kt0z1,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(cond_table%rho_grid(r_index_lo),cond_table%rho_grid(r_index),cond_table%rho_grid(r_index+1), &
           cond_table%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           cond_table%log10_kappa_table(t_index+1,r_index_lo,z_index+1), &
           cond_table%log10_kappa_table(t_index+1,r_index,z_index+1), &
           cond_table%log10_kappa_table(t_index+1,r_index+1,z_index+1), &
           cond_table%log10_kappa_table(t_index+1,r_index_hi,z_index+1), &
           ckt1z1,drkt1z1,dr2kt1z1,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(cond_table%rho_grid(r_index_lo),cond_table%rho_grid(r_index),cond_table%rho_grid(r_index+1), &
           cond_table%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           cond_table%log10_kappa_table(t_index_hi,r_index_lo,z_index+1), &
           cond_table%log10_kappa_table(t_index_hi,r_index,z_index+1), &
           cond_table%log10_kappa_table(t_index_hi,r_index+1,z_index+1), &
           cond_table%log10_kappa_table(t_index_hi,r_index_hi,z_index+1), &
           cktpz1,drktpz1,dr2ktpz1,xr,ierr)
      if (ierr /= 0) return
! Linear interpolation in ZLG:
      xz1=(log10_ion_charge-cond_table%z_grid(z_index))/(cond_table%z_grid(z_index+1)-cond_table%z_grid(z_index))
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
      call cinterp3(cond_table%temp_grid(t_index_lo),cond_table%temp_grid(t_index),cond_table%temp_grid(t_index+1), &
           cond_table%temp_grid(t_index_hi),log10_temperature,t_index,cond_n_temp, &
           cktm,ckt0,ckt1,cktp, & ! input: values of lg kappa
           log10_conductivity,dlnkappa_dlnt,dt2k,xt,ierr) ! lg kappa, d lg k / d lg T, d2 lg k / d2 lg T
      if (ierr /= 0) return
      call cinterp3(cond_table%temp_grid(t_index_lo),cond_table%temp_grid(t_index),cond_table%temp_grid(t_index+1), &
           cond_table%temp_grid(t_index_hi),log10_temperature,t_index,cond_n_temp, &
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

      ierr = 0
      if (grid_index <= 0.or.grid_index >= grid_size) then
! 2026 audit (section 7): stop -> ierr (runtime table-range guard;
! see kap_lib's kap_get).
         write(*,*) 'CINTERP: N0 out of range'
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
      endif
      if (grid_index < grid_size-1) then
         hp=grid_hi2-grid_hi1 ! right adjoint interval
         v11=((val_hi1-val_lo1)/h**2+(val_hi2-val_hi1)/hp**2)/ &
         (1d0/h+1d0/hp) ! right derivative
      endif
      if (grid_index > 1.and.grid_index < grid_size-1) then   ! Cubic interpolation
         c2=3d0*(val_hi1-val_lo1)-h*(v11+2d0*v01)
         c3=h*(v01+v11)-2d0*(val_hi1-val_lo1)
         interp_value=val_lo1+v01*x+c2*interp_frac**2+c3*interp_frac**3
         interp_deriv=v01+(2d0*c2*interp_frac+3d0*c3*interp_frac**2)/h
         interp_deriv2=(2d0*c2+6d0*c3*interp_frac)/h**2
         return
      endif
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
      endif
end subroutine cinterp3

!   W.H.Press, B.P.Flannery, S.A.Teukolsky, W.T.Vetterling
!   Numerical Receipes(Cambridge Univ., 1986)
!     Given an array XX of length N, and given a value X,
!     returns a value JLO such that X is between XX(JLO) and XX(JLO+1).
!     XX must be monotonic, either increasing or decreasing.
!     JLO=0 or JLO=N is returned to indicate that X is out of range.
!     JLO on input is taken as the initial guess for JLO on output.
pure subroutine hunt(table,n,target_value,index)


      integer, intent(in) :: n
      real(dp), intent(in) :: table(*)
      real(dp), intent(in) :: target_value
      integer, intent(inout) :: index

      logical :: ascending
      integer :: index_hi, increment, index_mid

      ascending=table(n) > table(1) ! true if ascending order, false otherwise
      if (index <= 0.or.index > n) then ! Input guess not useful.
         index=0
         index_hi=n+1  ! go immediately to bisection
      else
      increment=1 ! set the hunting increment
      if (target_value >= table(index).eqv.ascending) then ! Hunt up:
        do
          index_hi=index+increment
        if (index_hi > n) then ! Done hunting, since off end of table
           index_hi=n+1
           exit
        elseif (target_value >= table(index_hi).eqv.ascending) then ! Not done hunting
           index=index_hi
           increment=increment+increment
        else
           exit
        endif
        end do
      else ! Hunt down:
         index_hi=index
         do
          index=index_hi-increment
        if (index < 1) then ! Done hunting, since off end of table
           index=0
           exit
        elseif (target_value < table(index).eqv.ascending) then ! Not done hunting
           index_hi=index
           increment=increment+increment ! so double the increment
        else
           exit
        endif ! Done hunting, value bracketed
         end do
      endif
      endif
!   Hunt is done, so begin the final bisection phase:
      do
      if (index_hi-index == 1) return
      index_mid=(index_hi+index)/2
      if (target_value >= table(index_mid).eqv.ascending) then
         index=index_mid
      else
         index_hi=index_mid
      endif
      end do
end subroutine hunt

end module conductive_opacity
