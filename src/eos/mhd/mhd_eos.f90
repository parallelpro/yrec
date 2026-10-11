!----------------------------------------------------------------------
! mhd_eos
!----------------------------------------------------------------------
! The Mihalas-Hummer-Daeppen equation of state, from tables in
! (log P, log T) for several compositions.
!     read_mhd_tables  open the 3 ZAMS-type and up to 5 centre-type
!                      table files and read them (at eos_init)
!     compute_mhd_eos  density and the thermodynamic derivatives at
!                      (log T, log P, X, Z)
! Under these, Daeppen's routines under their original names: mhdst,
! mhdst1, rabu and rtab read the tables; mhdpx, mhdpx1 and mhdpx2
! interpolate them. The tables live in mhd_eos_tables (mhd).
module mhd_eos
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use mhd_eos_tables
      use eos_mixture, only: compute_mean_weights
      implicit none
      private
      public :: compute_mhd_eos, read_mhd_tables

contains

!----------------------------------------------------------------------
! compute_mhd_eos (formerly meqos)
!----------------------------------------------------------------------
! Mihalas, Hummer, Dappen equation of state. Implemented by YC,
! March 7, 1990. Given log10(T), log10(P), X, and Z, calls the MHD
! table interpolator (mhdpx) and converts its output
! (mhd%mhd_output) into the same thermodynamic derivative set
! used by the rest of the EOS machinery.
subroutine compute_mhd_eos(log10_temperature, temperature, log10_pressure, &
     pressure, log10_density, density, hydrogen_fraction, metal_fraction, &
     beta, beta_inverse, beta14, ion_fraction, specific_gas_constant, &
     ion_mean_weight_inverse, electron_mean_weight_inverse, &
     electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
     specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
     dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
     specific_heat_cp_dt, specific_heat_cp_dp, ierr)

      use io_units
      use math_lib

      real(dp), intent(in) :: log10_temperature, log10_pressure
      real(dp), intent(out) :: temperature, pressure, &
           log10_density, density
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      real(dp), intent(out) :: beta, beta_inverse, beta14
      real(dp), intent(out) :: ion_fraction(3)
      real(dp), intent(out) :: specific_gas_constant, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, &
           electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
           dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp

      real(dp), parameter :: cnvs = 0.434294481d0
      real(dp) :: mhdpx_r10
      real(dp) :: chi_rho, chi_t, log10_specific_heat_cp
      real(dp) :: specific_gas_constant_check
      integer :: ion_idx

      integer, intent(out) :: ierr

      ierr = 0

      temperature = exp10(log10_temperature)
      pressure = exp10(log10_pressure)
      call mhdpx(log10_pressure, log10_temperature, hydrogen_fraction, &
           mhdpx_r10, ierr)
      if (ierr /= 0) return
      log10_density = mhd%mhd_output(i_mhd_log10_rho)
      density = exp10(log10_density)
      dlnrho_dlnp = 1.0d0/mhd%mhd_output(i_mhd_chi_rho)
      chi_rho = mhd%mhd_output(i_mhd_chi_rho)
      chi_t = mhd%mhd_output(i_mhd_chi_t)
      dlnrho_dlnt = -chi_t/chi_rho
      log10_specific_heat_cp = mhd%mhd_output(i_mhd_log10_cp)
      specific_heat_cp = exp10(log10_specific_heat_cp)
      specific_heat_cp_dp = dlnrho_dlnp*mhd%mhd_output(i_mhd_dlogcp_dlogrho)
      specific_heat_cp_dt = mhd%mhd_output(i_mhd_dlogcp_dlogt) + mhd%mhd_output(i_mhd_dlogcp_dlogrho)*dlnrho_dlnt
      adiabatic_gradient = mhd%mhd_output(i_mhd_grad_ad)
      adiabatic_gradient_dp = dlnrho_dlnp*mhd%mhd_output(i_mhd_dgrad_ad_dlog10_rho)*cnvs/adiabatic_gradient
      adiabatic_gradient_dt = (mhd%mhd_output(i_mhd_dgrad_ad_dlog10_t) + mhd%mhd_output(i_mhd_dgrad_ad_dlog10_rho)*dlnrho_dlnt)* &
           cnvs/adiabatic_gradient
      dlnrho_dlnp_dt = dlnrho_dlnt*(adiabatic_gradient_dp - 1.0d0 + &
           dlnrho_dlnp + specific_heat_cp_dp)
      dlnrho_dlnt_dt = dlnrho_dlnt*(adiabatic_gradient_dt + 1.0d0 + &
           dlnrho_dlnt + specific_heat_cp_dt)
      beta = exp10((mhd%mhd_output(i_mhd_log10_pgas) - mhd%mhd_output(i_mhd_log10_p)))
      beta_inverse = 1.0d0/beta
      beta14 = 1.0d0 - beta
      do ion_idx = 1, 3
      ion_fraction(ion_idx) = mhd%mhd_output(ion_idx+i_mhd_ion_frac_1-1)
      end do
      electron_degeneracy_parameter = mhd%mhd_output(i_mhd_eta)
      specific_gas_constant = exp10(mhd%mhd_output(i_mhd_log10_pgas))/(density*temperature)
      call compute_mean_weights(temperature, pressure, density, hydrogen_fraction, &
           metal_fraction, specific_gas_constant_check, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, beta)
      if (abs((specific_gas_constant-specific_gas_constant_check)/ &
           specific_gas_constant).gt.5.0d-7) then
          write(terminal_unit,*)' ERROR(MHD) IN MEAN WEIGHTS ... '
          write(terminal_unit,*) specific_gas_constant, &
               specific_gas_constant_check
        write(terminal_unit,*) 'ERROR (MHD): CHECK MU'
          write(run_log_unit,*)' ERROR(MHD): IN MEAN WEIGHTS ... '
          write(run_log_unit,*) specific_gas_constant, &
               specific_gas_constant_check
        write(run_log_unit,*) 'ERROR (MHD): CHECK MU'
          ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
          ! facades stop when their caller passes no ierr.
          ierr = 1
          return
      end if
      return
end subroutine compute_mhd_eos

!----------------------------------------------------------------------
! read_mhd_tables (formerly mhdtbl)
!----------------------------------------------------------------------
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! Opens the 8 MHD equation-of-state table files (3 ZAMS-type tables
! and optionally 5 centre-type tables) on the unit numbers held in
! star%ctrl, calls mhdst to read them, and closes the files again.
subroutine read_mhd_tables(zams_a_table_path, zams_b_table_path, zams_c_table_path, &
     centre1_table_path, centre2_table_path, centre3_table_path, &
     centre4_table_path, centre5_table_path, ierr)

      character(len=256), intent(in) :: zams_a_table_path, zams_b_table_path, &
           zams_c_table_path, centre1_table_path, centre2_table_path, &
           centre3_table_path, centre4_table_path, centre5_table_path
      integer, intent(out) :: ierr
      integer :: unit_zams_a, unit_zams_b, unit_zams_c, unit_centre1, unit_centre2, unit_centre3, unit_centre4, unit_centre5
      ierr = 0

      open(newunit=unit_zams_a, file=zams_a_table_path, status='OLD', &
          form='UNFORMATTED')
      open(newunit=unit_zams_b, file=zams_b_table_path, status='OLD', &
          form='UNFORMATTED')
      open(newunit=unit_zams_c, file=zams_c_table_path, status='OLD', &
          form='UNFORMATTED')
      open(newunit=unit_centre1, file=centre1_table_path, status='OLD', &
          form='UNFORMATTED')
      open(newunit=unit_centre2, file=centre2_table_path, status='OLD', &
          form='UNFORMATTED')
      open(newunit=unit_centre3, file=centre3_table_path, status='OLD', &
          form='UNFORMATTED')
      open(newunit=unit_centre4, file=centre4_table_path, status='OLD', &
          form='UNFORMATTED')
      open(newunit=unit_centre5, file=centre5_table_path, status='OLD', &
          form='UNFORMATTED')

      call mhdst(unit_zams_a, unit_zams_b, unit_zams_c, unit_centre1, &
                 unit_centre2, unit_centre3, unit_centre4, unit_centre5, ierr)
      if (ierr /= 0) return
      close(unit_zams_a)
      close(unit_zams_b)
      close(unit_zams_c)
      close(unit_centre1)
      close(unit_centre2)
      close(unit_centre3)
      close(unit_centre4)
      close(unit_centre5)
      return
end subroutine read_mhd_tables

!----------------------------------------------------------------------
! mhdst
!----------------------------------------------------------------------
!
! Reads the 8 MHD equation-of-state table files (3 ZAMS-type unit
! numbers covering compositions A, B, C, plus up to 5 optional
! centre-type unit numbers) into mhd_eos_tables's mhd state and sets
! the table temperature limits used by mhdpx1.
subroutine mhdst(unit_zams_a, unit_zams_b, unit_zams_c, unit_centre1, &
     unit_centre2, unit_centre3, unit_centre4, unit_centre5, ierr)
! `use const_lib` removed (2026): unused (nothing in this file's body
! reaches a const_lib symbol), and its dummy arguments above now share
! names with the const_lib lunum members, which would otherwise be an
! ambiguous reference.

      integer, parameter :: ivarc = mhd_ivarc
      integer, parameter :: ivarx = mhd_ivarx
      integer, parameter :: nchem0 = mhd_nchem0
      integer, parameter :: nt1m = mhd_nt1m
      integer, parameter :: nt2m = mhd_nt2m
      integer, parameter :: ntxm = mhd_ntxm
      integer, parameter :: nr1m = mhd_nr1m
      integer, parameter :: nr2m = mhd_nr2m
      integer, parameter :: nrxm = mhd_nrxm

      integer, intent(in) :: unit_zams_a, unit_zams_b, unit_zams_c, &
           unit_centre1, unit_centre2, unit_centre3, unit_centre4, &
           unit_centre5

      real(dp) :: unused_zams_table(nt1m,nr1m,ivarc), unused_zams_log10t(nt1m)
!     WORKING STORAGE FOR NUMERICAL X-DERIVATIVES. COMMON TO ALL
!     X-TABLES
      real(dp) :: log10t_down(nt2m), log10t_up(nt2m)
      real(dp) :: table_centre_vars(nt2m,nr2m,ivarc)
      real(dp) :: table_down_vars(nt2m,nr2m,ivarc)
      real(dp) :: table_up_vars(nt2m,nr2m,ivarc)
      real(dp) :: atomic_weight_down(nchem0), number_abundance_down(nchem0), &
           mass_fraction_down(nchem0)
      real(dp) :: atomic_weight_up(nchem0), number_abundance_up(nchem0), &
           mass_fraction_up(nchem0)
!     Placeholder actuals for the mhdst1 arguments that the centre-table
!     calls do not need (the ZAMS lower-table slots).
      integer :: num_chem_species, unused_num_t, unused_num_r, table_index
      real(dp) :: unused_drho
!     Readability W3 (2026): the eight unit numbers gathered per table
!     family so the eight mhdst1 calls (formerly written out one per
!     composition, in this same order) are two loops over the
!     composition index of mhd's arrays.
      integer :: unit_zams(mhd_n_zams), unit_centre(mhd_n_centre), ic

      integer, intent(out) :: ierr

      ierr = 0

      num_chem_species = 0
      unused_drho = 0.d0
      unused_num_t = 0
      unused_num_r = 0
      unit_zams(imhd_zams_a) = unit_zams_a
      unit_zams(imhd_zams_b) = unit_zams_b
      unit_zams(imhd_zams_c) = unit_zams_c
      unit_centre(imhd_centre_1) = unit_centre1
      unit_centre(imhd_centre_2) = unit_centre2
      unit_centre(imhd_centre_3) = unit_centre3
      unit_centre(imhd_centre_4) = unit_centre4
      unit_centre(imhd_centre_5) = unit_centre5
!     READ ZAMS TABLES
      do ic = 1, mhd_n_zams
         if (unit_zams(ic).gt.0) then
            table_index = 0
            call mhdst1(unit_zams(ic),table_index,nt1m,nr1m,ivarc,nt2m,nr2m,ivarc,nchem0, &
                        mhd%zams_lower_num_t,mhd%zams_lower_num_r,mhd%zams_upper_num_t,mhd%zams_upper_num_r, &
                        mhd%zams_lower_log10t,mhd%zams_upper_log10t, &
                        mhd%zams_lower_table(:,:,:,ic),mhd%zams_upper_table(:,:,:,ic), &
                        mhd%zams_lower_drho,mhd%zams_upper_drho,num_chem_species, &
                        mhd%zams_atomic_weight(:,ic),mhd%zams_number_abundance(:,ic), &
                        mhd%zams_mass_fraction(:,ic),mhd%zams_mean_molecular_weight(ic), &
                        log10t_down,log10t_up,table_centre_vars,table_down_vars,table_up_vars, &
                        atomic_weight_down,atomic_weight_up, &
                        number_abundance_down,number_abundance_up,mass_fraction_down,mass_fraction_up, ierr)
            if (ierr /= 0) return
         end if
      end do
!     READ CENTRE TABLES
      do ic = 1, mhd_n_centre
         if (unit_centre(ic).gt.0) then
            table_index = 1
            call mhdst1(unit_centre(ic),table_index,nt1m,nr1m,ivarc,ntxm,nrxm,ivarx,nchem0, &
                        unused_num_t,unused_num_r,mhd%centre_num_t,mhd%centre_num_r,unused_zams_log10t, &
                        mhd%centre_log10t,unused_zams_table,mhd%centre_table(:,:,:,ic), &
                        unused_drho,mhd%centre_drho,num_chem_species, &
                        mhd%centre_atomic_weight(:,ic),mhd%centre_number_abundance(:,ic), &
                        mhd%centre_mass_fraction(:,ic),mhd%centre_mean_molecular_weight(ic), &
                        log10t_down,log10t_up,table_centre_vars,table_down_vars,table_up_vars, &
                        atomic_weight_down,atomic_weight_up, &
                        number_abundance_down,number_abundance_up,mass_fraction_down,mass_fraction_up, ierr)
            if (ierr /= 0) return
         end if
      end do
!     TEMPERATURE LIMITS
      mhd%table_log10t_min = mhd%zams_lower_log10t(  1)
      mhd%zams_lower_upper_boundary_log10t = mhd%zams_lower_log10t(mhd%zams_lower_num_t)
      if (unit_centre1.le.0) then
         mhd%zams_centre_boundary_log10t = mhd%zams_upper_log10t(mhd%zams_upper_num_t)
         mhd%table_log10t_max = mhd%zams_centre_boundary_log10t
      else
         mhd%zams_centre_boundary_log10t = mhd%centre_log10t(  1)
         mhd%table_log10t_max = mhd%centre_log10t(mhd%centre_num_t)
      end if
      return
end subroutine mhdst

!----------------------------------------------------------------------
! mhdst1
!----------------------------------------------------------------------
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! Reads one MHD equation-of-state table unit. table_kind=0 reads a
! single ZAMS-type table pair (lower/upper T region). table_kind=1
! reads a centre-type table plus the down/up-composition tables
! needed to build numerical X-derivatives, and stores the extended
! variable set (indices ivar1+1..ivar2) used in T-rho regions with
! inhomogeneous composition.
subroutine mhdst1(table_unit,table_kind,nt1m,nr1m,ivar1,nt2m,nr2m,ivar2,nchem0, &
                  num_t1,num_r1,num_t2,num_r2,log10t1,log10t2,table_vars1,table_vars2, &
                  drho1,drho2,num_chem,atomic_weight,number_abundance,mass_fraction,mean_molecular_weight, &
                  log10t_down,log10t_up,table_vars_centroid,table_vars_down,table_vars_up, &
                  atomic_weight_down,atomic_weight_up, &
                  number_abundance_down,number_abundance_up,mass_fraction_down,mass_fraction_up, ierr)
      use mhd_eos_tables, only: i_mhd_log10_p, i_mhd_log10_u, i_mhd_grad_ad, &
           i_mhd_log10_cp, i_mhd_dlog10p_dx, i_mhd_dlog10u_dx, &
           i_mhd_dgrad_ad_dx, i_mhd_dlog10cp_dx, i_mhd_placeholder
      integer, intent(in) :: table_unit, table_kind, nt1m, nr1m, ivar1, &
           nt2m, nr2m, ivar2, nchem0
      integer, intent(inout) :: num_t1, num_r1, num_t2, num_r2, num_chem
      real(dp), intent(inout) :: log10t1(nt1m), &
           table_vars1(nt1m,nr1m,ivar1), &
           log10t2(nt2m), table_vars2(nt2m,nr2m,ivar2)
      real(dp), intent(inout) :: drho1, drho2
      real(dp), intent(inout) :: atomic_weight(nchem0), &
           number_abundance(nchem0), mass_fraction(nchem0), mean_molecular_weight

      real(dp) :: log10t_down(nt2m), log10t_up(nt2m)
      real(dp) :: table_vars_centroid(nt2m,nr2m,ivar1)
      real(dp) :: table_vars_down(nt2m,nr2m,ivar1)
      real(dp) :: table_vars_up(nt2m,nr2m,ivar1)
      real(dp) :: atomic_weight_down(nchem0), number_abundance_down(nchem0), &
           mass_fraction_down(nchem0)
      real(dp) :: atomic_weight_up(nchem0), number_abundance_up(nchem0), &
           mass_fraction_up(nchem0)
      integer :: num_composition_reads, composition_pass, num_vars_read, &
           table_kind_read, composition_flag_read, composition_flag_expected
      real(dp) :: delta_x, composition_tolerance
      real(dp) :: unused_mean_molecular_weight_down, &
           unused_mean_molecular_weight_up
      integer :: species_index, temp_check_index, temp_deriv_index, &
           density_deriv_index, var_index

      integer, intent(out) :: ierr

      ierr = 0

      if (table_kind.eq.0) then
          num_composition_reads = 1
      else
          num_composition_reads = 3
      end if
      do composition_pass=1,num_composition_reads
      read(table_unit,   end=1000) num_vars_read,table_kind_read,composition_flag_read,delta_x
      if (ivar1.lt.num_vars_read) then
         ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      if (table_kind.ne.table_kind_read) then
         ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      if (composition_pass.eq.1) then
          composition_flag_expected= 0
      else if (composition_pass.eq.2) then
          composition_flag_expected=-1
      else if (composition_pass.eq.3) then
          composition_flag_expected= 1
      end if
      if (composition_flag_expected.ne.composition_flag_read) then
         ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      if (composition_pass.eq.1) call rabu(table_unit,nchem0,num_chem,atomic_weight,number_abundance,mass_fraction,mean_molecular_weight,ierr)
      if (composition_pass.eq.2) call rabu(table_unit,nchem0,num_chem,atomic_weight_down,number_abundance_down,mass_fraction_down,unused_mean_molecular_weight_down,ierr)
      if (composition_pass.eq.3) call rabu(table_unit,nchem0,num_chem,atomic_weight_up,number_abundance_up,mass_fraction_up,unused_mean_molecular_weight_up,ierr)
      if (ierr /= 0) return
      read(table_unit     ) num_t1,num_t2,drho1,drho2
      if (table_kind.eq.1 .and. num_t1.ne.0) then
          ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
          ! facades stop when their caller passes no ierr.
          ierr = 1
          return
      end if
      if (num_t1.gt.0) then
         call rtab(table_unit,nt1m,nr1m,ivar1,num_t1,num_r1,log10t1,table_vars1, ierr)
         if (ierr /= 0) return
      end if
      if (table_kind.eq.1) then
       if (composition_pass.eq.1) call rtab(table_unit,nt2m,nr2m,ivar1,num_t2,num_r2,log10t2,table_vars_centroid, ierr)
       if (ierr /= 0) return
       if (composition_pass.eq.2) call rtab(table_unit,nt2m,nr2m,ivar1,num_t2,num_r2,log10t_down,table_vars_down, ierr)
       if (ierr /= 0) return
       if (composition_pass.eq.3) call rtab(table_unit,nt2m,nr2m,ivar1,num_t2,num_r2,log10t_up,table_vars_up, ierr)
       if (ierr /= 0) return
      else if (table_kind.eq.0) then
       call rtab(table_unit,nt2m,nr2m,ivar2,num_t2,num_r2,log10t2,table_vars2, ierr)
       if (ierr /= 0) return
      end if
      end do
      if (table_kind.ne.0) then
!     IF IDX=1: CHECK TABLES FOR CORRECT COMPOSITION CONSTRUCTION
!     AND PERFORM NUMERICAL DERIVATIVES W.R.T. X
      composition_tolerance = 0.05d0*abs(delta_x)
      if ( abs(mass_fraction(1)-mass_fraction_down(1)-delta_x).gt.composition_tolerance  .or. abs(mass_fraction(1)-mass_fraction_up(1)+delta_x).gt.composition_tolerance  .or. abs(mass_fraction(2)-mass_fraction_down(2)+delta_x).gt.composition_tolerance  .or. abs(mass_fraction(2)-mass_fraction_up(2)-delta_x).gt.composition_tolerance ) then
         ierr = 1
         return
      end if
      do species_index=3,num_chem
      if ( abs(mass_fraction(species_index)-mass_fraction_up(species_index)).gt.composition_tolerance ) then
         ierr = 1
         return
      end if
      if ( abs(mass_fraction(species_index)-mass_fraction_down(species_index)).gt.composition_tolerance ) then
         ierr = 1
         return
      end if
      end do
      do temp_check_index=1,num_t2
      if (log10t2(temp_check_index).ne.log10t_down(temp_check_index)) then
         ierr = 1
         return
      end if
      if (log10t2(temp_check_index).ne.log10t_up(temp_check_index)) then
         ierr = 1
         return
      end if
      end do
!     NUMERICAL DERIVATIVES W.R.T. X
      do temp_deriv_index =1,num_t2
      do density_deriv_index =1,num_r2
      do var_index=1,ivar1
        table_vars2(temp_deriv_index,density_deriv_index,var_index)=table_vars_centroid(temp_deriv_index,density_deriv_index,var_index)
      end do
!
!     EXTENDED SET OF VARIABLES (TDVAR2(N,M,IVAR1+1...IVAR2))
!     FOR T-RHO REGIONS WITH INHOMOGENEOUS COMPOSITION.
!     IN THE COMMENTS,R AND T DENOTE LOG10(RHO) AND LOG10(T).
!     log10(P)/DX,log10(U)/DX,DDELAD/DX,log10(CP)/DX
      table_vars2(temp_deriv_index,density_deriv_index,i_mhd_dlog10p_dx)=(table_vars_up(temp_deriv_index,density_deriv_index,i_mhd_log10_p)-table_vars_down(temp_deriv_index,density_deriv_index,i_mhd_log10_p))/(2.d0*delta_x)
      table_vars2(temp_deriv_index,density_deriv_index,i_mhd_dlog10u_dx)=(table_vars_up(temp_deriv_index,density_deriv_index,i_mhd_log10_u)-table_vars_down(temp_deriv_index,density_deriv_index,i_mhd_log10_u))/(2.d0*delta_x)
      table_vars2(temp_deriv_index,density_deriv_index,i_mhd_dgrad_ad_dx)=(table_vars_up(temp_deriv_index,density_deriv_index,i_mhd_grad_ad)-table_vars_down(temp_deriv_index,density_deriv_index,i_mhd_grad_ad))/(2.d0*delta_x)
      table_vars2(temp_deriv_index,density_deriv_index,i_mhd_dlog10cp_dx)=(table_vars_up(temp_deriv_index,density_deriv_index,i_mhd_log10_cp)-table_vars_down(temp_deriv_index,density_deriv_index,i_mhd_log10_cp))/(2.d0*delta_x)
!     SPACE-HOLDER VARIABLE (LIKE VAR(20))
      table_vars2(temp_deriv_index,density_deriv_index,i_mhd_placeholder)=8888844444.d0
      end do
      end do
      end if
      return
!     EOF while reading a table header
 1000 continue
      ierr = 1
      return
end subroutine mhdst1

!----------------------------------------------------------------------
! rabu
!----------------------------------------------------------------------
! Reads one composition-table record (element count, atomic weights,
! number abundances, mass fractions, and mean molecular weight) from
! an unformatted file unit.
subroutine rabu(table_unit, nchem0, num_chem, atomic_weight, &
     number_abundance, mass_fraction, mean_molecular_weight, ierr)


      integer, intent(in) :: table_unit, nchem0
      integer, intent(out) :: num_chem
      real(dp), intent(out) :: atomic_weight(nchem0), &
           number_abundance(nchem0), mass_fraction(nchem0)
      real(dp), intent(out) :: mean_molecular_weight
      integer, intent(out) :: ierr

      integer :: ic

      ierr = 0

      read(table_unit   ) num_chem,(atomic_weight(ic),number_abundance(ic), &
           mass_fraction(ic), ic=1,num_chem),mean_molecular_weight
      if (nchem0.lt.num_chem) then
! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      return
end subroutine rabu

!----------------------------------------------------------------------
! rtab
!----------------------------------------------------------------------
!
! Reads an unformatted opacity/EOS-style table of shape
! (num_t_points, num_rho_points, num_vars) from file_unit: for each of
! num_t_points temperature records, reads the record's density-point
! count (asserted equal for every record) and log10(T), then reads
! num_vars variables at each of its density points.
subroutine rtab(file_unit,max_t_points,max_rho_points,num_vars, &
     num_t_points,num_rho_points,log_t,table_data, ierr)

      integer, intent(in) :: file_unit, max_t_points, max_rho_points, &
           num_vars, num_t_points
      integer, intent(out) :: num_rho_points

!     NT IS INPUT; NR,TL,TDVAR ARE OUTPUT
      real(dp), intent(out) :: log_t(max_t_points), &
           table_data(max_t_points,max_rho_points,num_vars)
! --- locals ---
      integer :: t_idx, rho_idx, var_idx, rho_count_read

      integer, intent(out) :: ierr

      ierr = 0

      do t_idx = 1, num_t_points
      read(file_unit     ) rho_count_read,log_t(t_idx)
      if(t_idx.eq.1) num_rho_points=rho_count_read
!     CHECK LIMITS AND NUMBER OF DENSITY POINTS OF TABLE
      if ( num_t_points.gt.max_t_points .or.  num_rho_points.gt.max_rho_points ) then
         ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      if ( t_idx.gt.1 .and. rho_count_read.ne.num_rho_points ) then
         ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      do rho_idx = 1, num_rho_points
      read(file_unit     ) (table_data(t_idx,rho_idx,var_idx),var_idx=1,num_vars)
      end do
      end do
      return
end subroutine rtab

!----------------------------------------------------------------------
! mhdpx
!----------------------------------------------------------------------
!
!     MHD INTERPOLATION PACKAGE
!     PRESSURE AS ARGUMENT INTERPOLATION IN X, Z FIXED
!
!     INPUT   PGL = LOG10 (GAS PRESSURE [DYN/CM2])
!             TL  = LOG10 (TEMPERATURE [K])
!             XC  = HYDROGEN ABUNDANCE (BY MASS)
!
!     OUTPUT  RL  = LOG10 (DENSITY [G/CM3])        : ARGUMENT
!             OTHER (SEE SEPARATE INSTRUCTIONS)    : mhd%mhd_output
!
!     ERROR   IERR = 1 SIGNALS PGL,TL OUTSIDE THE DOMAIN OF TABLES
!             (OTHERWISE IERR = 0).
!
! Top-level entry point to the MHD (Mihalas, Hummer, Dappen)
! equation-of-state table interpolator. S/R MHDST (via read_mhd_tables) must
! be called in main to load the tables first. Dispatches to the
! variable-X interpolation in mhdpx1, then returns log10(density).
subroutine mhdpx(log10_pressure, log10_temperature, hydrogen_fraction, &
     log10_density, ierr)

      real(dp), intent(in) :: log10_pressure, log10_temperature, &
           hydrogen_fraction
      real(dp), intent(out) :: log10_density
      integer, intent(out) :: ierr

      ierr = 0

      call mhdpx1(log10_pressure, log10_temperature, hydrogen_fraction, ierr)
      if (ierr /= 0) return
      log10_density = mhd%mhd_output(i_mhd_log10_rho)
      return
end subroutine mhdpx

!----------------------------------------------------------------------
! mhdpx1
!----------------------------------------------------------------------
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     MHDST MUST BE CALLED IN MAIN.
!     INTERPOLATION IN TABLES WITH DIFFERENT X AND FIXED Z
!     TL < TLIM1:         LOWER PART OF ZAMS TABLES
!     TLIM1 < TL < TLIM2: UPPER PART OF ZAMS TABLES
!     TL > TLIM2:         CENTRE TABLES
!     TMINI,TMAXI:        TEMPERATURE INTERVAL COVERED
!                         BY THE TABLES
!     (TLIM1/TLIM2/TMINI/TMAXI are mhd%zams_lower_upper_boundary_log10t,
!     zams_centre_boundary_log10t, table_log10t_min, table_log10t_max.)
!     An out-of-range TL returns silently with ierr = 0 and
!     mhd%mhd_output untouched (historical behaviour, preserved).
subroutine mhdpx1(log10_pressure, log10_temperature, hydrogen_fraction, ierr)
      use io_units
      use numerics_lib
      integer, parameter :: ivarx = mhd_ivarx
      integer, parameter :: ndimt = 8

      real(dp), intent(in) :: log10_pressure, log10_temperature, &
           hydrogen_fraction

      real(dp) :: table_vars(ndimt,ivarx), table_hfrac(ndimt)

!     QUANTITIES FOR INTERPOLATION IN X
      real(dp) :: cubic_vars(ivarx,4), cubic_x_nodes(4)
      integer :: cubic_table_index(4)
! LIR's "type" flag (an integer, despite the original's L-prefixed name).
      integer :: lir_type_flag
!     READ FROM APPROPRIATE TABLES
!     AND FILL ARRAYS VAROUT AND XC
!     CAN WE DO IT?
!
      integer :: i, itbl, iv, ixmin, num_vars, var_leading_dim, &
           num_points, interp_mode
      real(dp) :: x_grid_origin, x_grid_spacing, &
           var_at_x0, var_at_x1, var_at_x2

      integer, intent(out) :: ierr

      ierr = 0

      if (log10_temperature.lt.mhd%table_log10t_min .or. &
           log10_temperature.gt.mhd%table_log10t_max) then
          return
      end if
!     LOWER ZAMS TABLES
      if (log10_temperature.lt.mhd%zams_lower_upper_boundary_log10t) then
         do i=1,3
         itbl = -i
         call mhdpx2(log10_pressure, log10_temperature, itbl, table_vars, table_hfrac, ndimt)
         end do
      else if (log10_temperature.lt.mhd%zams_centre_boundary_log10t) then
!     UPPER ZAMS TABLES
         do i=1,3
         itbl = i
         call mhdpx2(log10_pressure, log10_temperature, itbl, table_vars, table_hfrac, ndimt)
         end do
      else
!     CENTRE TABLES
         do i=1,5
!        OFFSET (+3) TO ACCESS CENTER TABLES
         itbl = i + 3
         call mhdpx2(log10_pressure, log10_temperature, itbl, table_vars, table_hfrac, ndimt)
         end do
      end if
!     INTERPOLATION IN X
      if (log10_temperature.le.mhd%zams_centre_boundary_log10t) then
!        QUADRATIC NEWTON (EQUIDISTANT XC'S)
         x_grid_origin = table_hfrac(1)
         x_grid_spacing  = table_hfrac(2) - table_hfrac(1)
         if (abs(table_hfrac(3)-table_hfrac(2)-x_grid_spacing).gt.1.d-4) then
            write(terminal_unit,*) 'ERROR (MHD): NON-EQUIDISTANT ZAMS TABLES.'
            write(terminal_unit,*) 'XC(1-3)= ',table_hfrac(1),table_hfrac(2),table_hfrac(3)
            write(run_log_unit,*) 'ERROR (MHD): NON-EQUIDISTANT ZAMS TABLES.'
            write(run_log_unit,*) 'XC(1-3)= ',table_hfrac(1),table_hfrac(2),table_hfrac(3)
            ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
            ! facades stop when their caller passes no ierr.
            ierr = 1
            return
         end if
         do iv=1,ivarx
         var_at_x0 = table_vars(1,iv)
         var_at_x1 = table_vars(2,iv)
         var_at_x2 = table_vars(3,iv)
         call quint(hydrogen_fraction, x_grid_origin, x_grid_spacing, &
              var_at_x0, var_at_x1, var_at_x2, mhd%mhd_output(iv))
         end do
      else
!     CUBIC LAGRANGIAN (ARBITRARILY SPACED XC'S)
!     USE 4 HIGHEST XC(I)'S, IF POSSIBLE.
!     DISTINGUISH CASE WHERE THE XC(I), I=4..8,
!     INCREASE OR DECREASE.
         if (table_hfrac(4).lt.table_hfrac(8)) then
            ixmin = 5
            if (hydrogen_fraction.lt.table_hfrac(5)) ixmin=4
         else
            ixmin = 4
            if (hydrogen_fraction.lt.table_hfrac(7)) ixmin=5
         end if
         do i=1,4
            cubic_table_index(i)=ixmin+i-1
         end do
         do i =1,4
            cubic_x_nodes(i) = table_hfrac(cubic_table_index(i))
            do iv=1,ivarx
               cubic_vars(iv,i) = table_vars(cubic_table_index(i),iv)
            end do
         end do
         num_vars=ivarx
         var_leading_dim=ivarx
         num_points=4
         lir_type_flag=1
         interp_mode=1
         call lir(hydrogen_fraction, cubic_x_nodes, mhd%mhd_output, cubic_vars, &
              num_vars, var_leading_dim, num_points, lir_type_flag, interp_mode)
      end if
      return
end subroutine mhdpx1

!----------------------------------------------------------------------
! mhdpx2
!----------------------------------------------------------------------
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! Interpolates one table (selected by table_selector: -1/-2/-3 lower
! ZAMS A/B/C, 1/2/3 upper ZAMS A/B/C, 4..8 centre 1..5) at
! (log10 P, log10 T) and stores the result as row |table_selector|
! (ZAMS) or table_selector (centre) of table_vars, with that table's
! hydrogen mass fraction in table_hfrac. The tables are the
! composition-indexed arrays of mhd (imhd_zams_a..c,
! imhd_centre_1..5 in mhd_eos_tables).
subroutine mhdpx2(log10_pressure, log10_temperature, table_selector, &
     table_vars, table_hfrac, ndimt)
      use numerics_lib
      integer, parameter :: ivarc = mhd_ivarc
      integer, parameter :: ivarx = mhd_ivarx
      integer, parameter :: nt1m = mhd_nt1m
      integer, parameter :: nt2m = mhd_nt2m
      integer, parameter :: ntxm = mhd_ntxm
      integer, parameter :: nr1m = mhd_nr1m
      integer, parameter :: nr2m = mhd_nr2m
      integer, parameter :: nrxm = mhd_nrxm
!     ZAMS TABLES (LABELLED BY A,B,C)
      integer, intent(in) :: table_selector, ndimt
      real(dp), intent(in) :: log10_pressure, log10_temperature
      real(dp), intent(inout) :: table_vars(ndimt,ivarx), &
           table_hfrac(ndimt)
!     NOMENCLATURE FOR ACCESSING THE TABLES
!     ZAMS TABLES
!     CONTRAL (VARIABLE X) TABLES
!     NOTE: THE SAME OUTPUT VARIABLES VAROUT(ITBL,.) IS USED
!           FOR ITBL=1 AND -1, 2 AND -2, 3 AND -3, RESPECTIVELY.
      real(dp) :: intpt_work1(ivarx,4), intpt_work2(ivarx,4), &
           intpt_y(ivarx), interpolated_vars(ivarx)

      integer :: i, ic

      call zero(interpolated_vars, ivarx)
!     MAIN SELECTION OF TABLES
!     Readability W3 (2026): the former 11-way ladder over the named
!     per-composition members is now index arithmetic into mhd's
!     composition-indexed arrays: selector -ic is lower ZAMS composition
!     ic (row ic), +ic is upper ZAMS composition ic (row ic), and
!     mhd_n_zams+ic is centre composition ic (row mhd_n_zams+ic). Any
!     other selector still does nothing, as before.
!     ZAMS TABLES
      if (table_selector.le.-1 .and. table_selector.ge.-mhd_n_zams) then
         ic = -table_selector
         call interpolate_and_store(mhd%zams_lower_table(:,:,:,ic), &
              nt1m, nr1m, ivarc, mhd%zams_lower_log10t, mhd%zams_lower_num_t, &
              mhd%zams_lower_num_r, ic, mhd%zams_mass_fraction(1,ic))
      else if (table_selector.ge.1 .and. table_selector.le.mhd_n_zams) then
         ic = table_selector
         call interpolate_and_store(mhd%zams_upper_table(:,:,:,ic), &
              nt2m, nr2m, ivarc, mhd%zams_upper_log10t, mhd%zams_upper_num_t, &
              mhd%zams_upper_num_r, ic, mhd%zams_mass_fraction(1,ic))
!     CENTER TABLES
      else if (table_selector.ge.mhd_n_zams+1 .and. &
           table_selector.le.mhd_n_zams+mhd_n_centre) then
         ic = table_selector - mhd_n_zams
         call interpolate_and_store(mhd%centre_table(:,:,:,ic), &
              ntxm, nrxm, ivarx, mhd%centre_log10t, mhd%centre_num_t, &
              mhd%centre_num_r, table_selector, mhd%centre_mass_fraction(1,ic))
      end if
!     END SELECTION OF TABLES
      return

contains

! Interpolate one table at (log10 P, log10 T) and store its num_vars
! results as row `row` of table_vars, with the table's hydrogen mass
! fraction in table_hfrac(row). The work arrays and the (P,T) point are
! host-associated from mhdpx2.
subroutine interpolate_and_store(table_data, table_dim_t, table_dim_r, &
     num_vars, table_log10t, num_t, num_r, row, hydrogen_mass_fraction)
      integer, intent(in) :: table_dim_t, table_dim_r, num_vars, num_t, num_r, row
      real(dp), intent(in) :: table_data(table_dim_t,table_dim_r,num_vars)
      real(dp), intent(in) :: table_log10t(table_dim_t)
      real(dp), intent(in) :: hydrogen_mass_fraction

      call intpt (log10_pressure, log10_temperature, table_data, &
           table_dim_t, table_dim_r, num_vars, table_log10t, num_t, &
           num_r, intpt_work1, intpt_work2, intpt_y, interpolated_vars)
      do i=1,num_vars
         table_vars(row,i) = interpolated_vars(i)
      end do
      table_hfrac(row)    = hydrogen_mass_fraction
end subroutine interpolate_and_store

end subroutine mhdpx2

end module mhd_eos
